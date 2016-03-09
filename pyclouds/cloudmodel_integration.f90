module cloudmodel_integration
   use integrators, only: integrate_with_message

   use cloudmodel_common, only: nvars

   use mo_ccfm_clouddata,  only: t_clinfo, t_cldata
   use mo_ccfm_parameter,  only: dp, clno, r_min, g

   implicit none
   public ccfm_interface

   logical, parameter :: debug = .true.

   real(dp), allocatable, dimension(:,:) :: F_ambient
   real(dp), allocatable, dimension(:) :: z_profile

   ! for holding the last height-index used for integrating
   integer :: m__integrate_last_idx = 0

contains
   !> Initiate the cloud-model integration. Primarily serves to init the
   !> microphysics routine and set all the state-mapping up correctly
   subroutine init()
      use microphysics_register, only: idx_pressure, idx_temp
      use microphysics_register, only: idx_water_vapour, idx_cwater, idx_rain
      use microphysics_register, only: idx_cice, idx_graupel
      use microphysics_register, only: mphys_nvars => n_variables
      use microphysics_initialisation, only: unified_microphysics_init => init

      use cloudmodel_common, only: i_p, i_T, i_qv, i_ql, i_qr, i_qi, i_r, i_w

      integer :: n_assigned_vars

      call unified_microphysics_init('no_ice', 'isobaric')
      n_assigned_vars = mphys_nvars

      i_p = idx_pressure
      i_T = idx_temp
      i_qv = idx_water_vapour
      i_ql = idx_cwater
      i_qr = idx_rain
      i_qi = idx_cice

      if (idx_graupel > 0) then
         print *, "Error: the CCFM cloud-model does not currently represent graupel"
         call exit(0);
      endif

      if (i_p == 0) then
         print *, "Error: the pressure variable has not been assigned correctly"
         call exit(0);
      endif

      if (i_T == 0) then
         print *, "Error: the temperature variable has not been assigned correctly"
         call exit(0);
      endif

      if (i_qv == 0) then
         print *, "Error: the water vapour variable has not been assigned correctly"
         call exit(0);
      endif

      if (i_ql == 0) then
         print *, "Error: the cloud water variable has not been assigned correctly"
         call exit(0);
      endif

      if (i_qr == 0) then
         print *, "Error: the rain variable has not been assigned correctly"
         call exit(0);
      endif

      if (i_qi == 0) then
         print *, "Warning: The selected microphysics implementation does not represent ice"
         n_assigned_vars = 1 + n_assigned_vars
         i_qi = n_assigned_vars
      endif

      ! the microphysics code doesn't use the radius and vertical velocity
      n_assigned_vars = 1 + n_assigned_vars
      i_r = n_assigned_vars

      n_assigned_vars = 1 + n_assigned_vars
      i_w = n_assigned_vars


      if (debug) then
         print *, "cloud-model internal state mapping"
         print *, 'i_r', i_r
         print *, 'i_w', i_w
         print *, 'i_p', i_p
         print *, 'i_T', i_T
         print *, 'i_qv', i_qv
         print *, 'i_ql', i_ql
         print *, 'i_qr', i_qr
         print *, 'i_qi', i_qi
      endif
   end subroutine

   !> subroutine which has the same call interface as the old `cloudmodel`
   !facilitating an identical call from CCFM, this will replace the
   !`cloud_driver` sub-routine that integrates all cloud profiles in a given
   !environment
   !>
   !! @todo what are the units of the water mixing rate argument?
   !! @todo what are the units of the turbulent kinetic energy argument?
   !! @param level number of vertical levels [1]
   !! @param pres pressure profile [Pa]
   !! @param geog geopotential height profile (sorted top down) [gpm]
   !! @param T_ls temperature profile [K]
   !! @param q_ls humidity profile [kg/kg]
   !! @param T_base temperature at cloud base [K]
   !! @param q_base water mix rate at cloud base [?]
   !! @param cl_base cloud base height index [1] 
   !! @param tke turbulent kinetic energy profile [?]
   !! @param r_max maximum cloud radius [m]  
   !! @param h_wind horizon wind speed profile [m/s]
   !! @param cdnc cloud droplet number concentration profile [1/m^3]
   !! @param land boolean representing calculation over land or not, if `true` then over land
   !! @param area grid box area at surface [m^2]
   !! @return cl_data vertical profile of each discrete cloud
   !! @return cl_info info about each discrete cloud
   !! @return dr  cloud radius increment used (calculate from input variables) [m]
   subroutine ccfm_interface(level,         & ! in
         pres,          & ! in
         geog,          & ! in
         T_ls,          & ! in 
         q_ls,          & ! in
         cl_base,       & ! in     
         T_base,        & ! in
         q_base,        & ! in
         tke,           & ! in  
         r_max,         & ! in
         h_wind,        & ! in  
         cdnc,          & ! in
         land,          & ! in  
         area,          & ! in  
         cl_info,       & ! out
         cl_data,       & ! out 
         dr         ) ! out

      use cloudmodel_common, only: i_r, i_w, i_T, i_p, i_qv

      integer,       intent(in)  :: level
      real(dp),      intent(in)  :: pres    (level)        ! pressure [Pa]
      real(dp),      intent(in)  :: geog    (level)        ! height - sorted top down [gpm]      
      real(dp),      intent(in)  :: T_ls    (level)        ! (Initial)temperatur [K]
      real(dp),      intent(in)  :: q_ls    (level)        ! (Initial)humidity [kg/kg]

      real(dp),      intent(in)  :: T_base                 ! temperature at cloud base
      real(dp),      intent(in)  :: q_base                 ! water mix rat at cloud base

      integer,       intent(in)  :: cl_base                ! cloudbase [level] 
      real(dp),      intent(in)  :: tke     (level)       
      real(dp),      intent(in)  :: r_max                  ! max cloud radius [m]  
      real(dp),      intent(in)  :: h_wind  (level)        ! horizantal wind speed [m/s]
      real(dp),      intent(in)  :: cdnc    (level)
      logical,       intent(in)  :: land                   ! land-sea mask
      real(dp),      intent(in)  :: area                   ! grid box area at surface [m^2]

      type(t_cldata),intent(out) :: cl_data (clno, level)
      type(t_clinfo),intent(out) :: cl_info (clno)
      real(dp),      intent(out) :: dr                 ! cloud radius increment [m]


      real(dp), dimension(nvars) :: F0 ! Initial condition, will be set for each cloud type
      integer :: i
      real(dp) :: z = 0.0


      ! Init ambient
      call init_ambient_state(n_levels=level, T_e=T_ls, qv_e=q_ls, p_e=pres, phi_geog=geog)



      if ( clno == 1 )  then 
         dr = 0.
      else
         dr = (r_max - r_min) / (clno - 1)
      end if

      ! for all diff. cloud radii
      do i = 1, clno
         ! the profiles are provided top-down so we'll need be interpolating from
         ! the "end"
         m__integrate_last_idx = level

         ! Create initial condition
         z = geog(cl_base)/g

         F0 = 0.0
         F0(i_w) = 1.0
         F0(i_r) = r_min + dr*(i-1)
         F0(i_T) = T_base
         F0(i_qv) = q_base

         call integrate_profile(F0=F0, cloud_number=clno)

      enddo


      ! Integrate
   end subroutine

   subroutine integrate_profile(F0, cloud_number)
      real(dp), intent(in), dimension(nvars) :: F0
      integer, intent(in) :: cloud_number
   end subroutine

   function dFdz(z, F)
      use cloudmodel_common, only: i_T, i_p, i_qv
      use spec_conc_eqns, only: model_dFdz => dFdz

      real(dp), intent(in), dimension(nvars) :: F
      real(dp), intent(in) :: z

      real(dp), dimension(nvars) :: dFdz, F_e

      real(dp) :: T_e, p_e, qv_e

      F_e = interpolate_ambient_state(z)

      dFdz = model_dFdz(F=F, F_e=F_e)
   end function

   subroutine init_ambient_state(n_levels, T_e, qv_e, p_e, phi_geog)
      use cloudmodel_common, only: i_T, i_p, i_qv

      integer, intent(in) :: n_levels
      real(dp), intent(in), dimension(n_levels) :: T_e, p_e, qv_e, phi_geog

      allocate(F_ambient(n_levels, nvars))
      F_ambient = 0.0_dp
      F_ambient(:, i_T) = T_e
      F_ambient(:, i_qv) = qv_e
      F_ambient(:, i_p) = p_e

      allocate(z_profile(n_levels))
      z_profile = phi_geog/g
      
   end subroutine

   recursive function interpolate_ambient_state(z) result(F_e)
      real(dp), intent(in) :: z
      real(dp), dimension(nvars) :: F_e, dF_e__dz

      real(dp) :: z_bottom, z_top
      integer :: k_b, k_t

      k_b = m__integrate_last_idx
      k_t = m__integrate_last_idx-1

      z_top = z_profile(k_t)

      if (z > z_top) then
         m__integrate_last_idx = m__integrate_last_idx - 1
         F_e = interpolate_ambient_state(z)
      endif

      z_bottom = z_profile(k_b)
      if (z < z_bottom) then
         m__integrate_last_idx = m__integrate_last_idx + 1
         F_e = interpolate_ambient_state(z)
      endif

      dF_e__dz = (F_ambient(k_t,:) - F_ambient(k_b,:))/(z_top - z_bottom)
      F_e = F_ambient(k_b,:) + (z - z_bottom)*dF_e__dz

   end function
end module

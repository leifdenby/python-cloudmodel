module ccfm
   use mo_ccfm_parameter

   implicit none 

contains

!> Calculate the vertical derivative of vertical velocity. Calculation is done
!> using the Lagrangian equation steady-state equation in Till's thesis (p 70).
!!
!! @returns vertical derivative of vertical velocity [m/s/m]
real(dp) function ddv_z(v_z, q_l, q_r, q_i, q_s, T_vc, T_ve, mu)
 
   !use mo_ccfm_parameter,      only: g
   !use mo_ccfm_clouddata, only: t_clinfo, t_cldata, t_cldata_size

      
   !type(t_cldata), intent(in) :: y     ! cloud data      
   
   real(dp),       intent(in) :: v_z   ! vertical velocity [m/s]
   real(dp),       intent(in) :: q_l, q_r, q_i, q_s


   real(dp),       intent(in) :: T_vc  ! virtual temperature in the cloud [K]
   real(dp),       intent(in) :: T_ve  ! virtual temperature environment  [K]
   real(dp),       intent(in) :: mu    ! entrainment factor               [1/m]
   real(dp)                   :: buoy  ! buoyancy acceleration            [1/s]
   real(dp)                   :: grav  ! gravity drag                     [1/s]
   real(dp)                   :: entr  ! entrainmenment factor            [1/s] 
 
   buoy = 0._dp
   grav = 0._dp
   entr = 0._dp
   
   if ( v_z > 0._dp ) then
       
      ! equation (3) [Weinstein... 1969]
      ! virtual mass coefficient: 
      buoy = g * ( T_vc - T_ve ) / T_ve / v_z / 1.4_dp   
      grav = - g * ( q_l + q_r + q_i + q_s ) / v_z
      entr = - mu / 2._dp * v_z

   end if
   
   ddv_z = buoy + grav + entr
     
   return

end function ddv_z


!> Calculate vertical derivative of cloud radius. Calculation is done using
!> Lagrangian form of 1D continuity equation, describing a steady-state cloud,
!> as in Till's thesis (p 73).
!! @warning there may be a typo in the "proportionality factor for the entrainment rate"
!! @param y struct of cloud data already known (vertical velocity and specific concentration of water phases)
!! @param ddT_z vertical derivative of temperature [K/m]
!! @param ddv_z vertical derivative of vertical velocity [1/s]
!! @param T_ve virtual temperature environment [K]
!! @param T_vc virtual temperature in the cloud [K]
!! @returns vertical derivative of cloud radius [m/m]
real(dp) function ddr_z(ddT_z, ddv_z, T_ve, T_vc, rad, v_z)

   !use mo_ccfm_parameter,      only: C_mu, g, Rd
   !use mo_ccfm_clouddata, only: t_clinfo, t_cldata, t_cldata_size

   implicit none 
      
   real(dp),       intent(in) :: rad, v_z

   !type(t_cldata), intent(in) :: y      ! cloud data   
   real(dp),       intent(in) :: ddT_z  ! temperature derivation           [K/m]
   real(dp),       intent(in) :: ddv_z  ! vertical velocity derivation     [1/s]
   real(dp),       intent(in) :: T_ve   ! virtual temperature environment  [K]
   real(dp),       intent(in) :: T_vc   ! virtual temperature in the cloud [K]
   real(dp)                   :: dr_1, dr_2, dr_3                    !     [-]
 
   ! radius - equation [Kreitzberg and Perkey 1976]

   dr_1 =  rad / 2._dp * (g / Rd / T_ve  +  ddT_z / T_vc)

   if ( v_z > 0._dp ) then
      dr_2 = - rad / (2._dp * v_z) * ddv_z
   else
      dr_2 = 0._dp
   end if

   dr_3 = C_mu / 2._dp

   ddr_z = dr_1 + dr_2 + dr_3

   return

end function ddr_z

 real(dp) function ddT_z(q_satw, T_c, mu, T_e, q_ve, p_e, dpdz)
      
   implicit none 

   real(dp), intent(in) :: q_satw  ! saturation mixing ratio  [kg/kg]
   real(dp), intent(in) :: T_c     ! cloud temperature        [K]
   real(dp), intent(in) :: mu      ! entrainment rate         [1/m]
   real(dp), intent(in) :: T_e     ! environment temperature  [K]
   real(dp), intent(in) :: q_ve    ! environment mixing ratio [kg/kg]
   real(dp), intent(in) :: p_e     ! pressure                 [Pa] 
   real(dp), intent(in) :: dpdz    ! pressure derivative      [Pa/m]
      
   real(dp)             :: dT_1, dT_2, dT_3 !                 [K/m]
   real(dp)             :: c
   
   ! c    : divisor of equation (1) in [Weinstein, MacCready 1969]
   ! dT_1 : first  summand of equation (1) in [Weinstein, MacCready 1969] 
   !            -> heating due to latent heat release due to condensation and
   !               cooling due to  lifting  (dry adiabatical lapse rate)
   ! dT_2 : second  -"-
   !            -> cooling due to entrainment of cooler environmental air
   ! dT_3 : third   -"-
   !            -> cooling due to entrainment of dry air and therefore 
   !               release  of latent heat due to evaporation of cloud water 

   c     =  1._dp + epsi * aLv*aLv * q_satw / (cpd * Rd * T_c*T_c)
   
   ! hydrostatic assumption  
   dT_1 = -g/cpd * ( 1._dp + aLv * q_satw / ( Rd * T_c ) )
   ! pressure gradient known 
   !dT_1  = -g/cpd  + ( aLv * q_satw * dpdz ) / ( cpd * p_e ) 
   
   dT_2  = -mu * ( T_c - T_e )
      
   dT_3  = -mu * ( q_satw - q_ve ) * aLv / cpd
      
   ddT_z = ( dT_1 + dT_2 + dT_3 ) / c
   
 end function ddT_z

end module ccfm

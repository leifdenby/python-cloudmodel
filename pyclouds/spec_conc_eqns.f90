module spec_conc_eqns
   use cloudmodel_common, only: i_r, i_w, i_T, i_p, i_qv, i_ql, i_qr, i_qi
   use cloudmodel_common, only: dp, nvars

   implicit none


   real(dp), parameter :: R_d = 287.05_dp
   real(dp), parameter :: R_v = 461.51_dp
   real(dp), parameter :: g = 9.80665_dp
   real(dp), parameter :: rho_i = 1000._dp
   real(dp), parameter :: rho_l = 500._dp
   real(dp), parameter :: cp_d = 1005.46_dp
   real(dp), parameter :: cp_v = 1859.0_dp
   real(dp), parameter :: L_v = 2.5008e6_dp
   real(dp), parameter :: L_s = 2.8345e6_dp

   real(dp), parameter :: m__beta = 0.2_dp
   real(dp), parameter :: m__gamma = 1.0_dp
   real(dp), parameter :: m__D = 1.0_dp



contains
   !> Compute vertical derivate of full state-vector describing the cloud-parcel state `F` 
   !> given the environmental conditions
   function dFdz(F, F_e)
      use microphysics_register, only: i_T_mphys => idx_temp
      use microphysics_register, only: i_p_mphys => idx_pressure
      use microphysics_register, only: i_qv_mphys => idx_water_vapour
      use microphysics_register, only: i_ql_mphys => idx_cwater
      use microphysics_register, only: i_qr_mphys => idx_rain
      use microphysics_register, only: i_qi_mphys => idx_cice

      ! TODO: get this from some pointer instead
      use isobaric_integration_helpers, only: dydt_mphys => dydt_isobaric


      real(dp), dimension(nvars), intent(in) :: F, F_e
      real(dp) :: p_e, rho_e, T_e

      real(dp), dimension(nvars) :: dFdz

      real(dp) :: r, w, T, p
      real(dp) :: q_v, q_l, q_i, q_r, q_d

      ! for holding the two contributions to temperature changes:
      ! adiabatic expansion + entrainment and phase changes (microphysics)
      real(dp) :: dTdz_ae, dTdz_mphys

      ! For holding the state given to and return from the microphysics
      ! TODO: the explicit number here shouldn't be necessary
      real(dp) :: y(8), dydt(8)

      real(dp) :: cm_m_p

      ! TODO: For now the moisture in the enviroment is not considered when
      ! calculating the changes wrt entrainment
      !qv_e = F_e(i_qv)
      T_e = F_e(i_T)
      p_e = F_e(i_p)
      rho_e = mod__mixture_density_from_eos(p=p_e, T=T_e, qd=1.0_dp, qv=0.0_dp, ql=0.0_dp, qi=0.0_dp)


      r = F(i_r)
      w = F(i_w)
      T = F(i_T)
      q_v = F(i_qv)
      q_l = F(i_ql)
      q_i = F(i_qi)
      q_r = F(i_qr)
      q_d = 1.0_dp - q_v - q_l - q_i - q_r

      ! cloud is assumed to be at same pressure as in environment
      p = p_e

      ! initiate to zero
      dFdz = 0.0_dp

      ! 1. Estimate change in vertical velocity with initial state
      dFdz(i_w) = dw_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, rho_e=rho_e)

      ! 2. Estimate changes from microphysics
      ! (need to create state vector to supply to microphysics)
      y = 0.0
      y(i_p_mphys) = p
      y(i_T_mphys) = T
      y(i_qv_mphys) = q_v
      y(i_ql_mphys) = q_l
      y(i_qr_mphys) = q_r
      y(i_qi_mphys) = q_i

      dydt = dydt_mphys(t=0.0_dp, y=y)

      ! w = dz/dt => dy/dt * 1/w = dy/dt * dt/dz = dy/dz
      dFdz(i_qv) = dydt(i_qv_mphys)/w
      dFdz(i_ql) = dydt(i_ql_mphys)/w
      dFdz(i_qr) = dydt(i_qr_mphys)/w
      dFdz(i_qi) = dydt(i_qi_mphys)/w


      ! 3. Estimate temperature change forgetting about phase-changes for now (i.e. considering only adiabatic adjustment and entrainment)
      dTdz_ae = dT_dz(r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, &
                      qi_c=q_i, dql_c__dz=0.0_dp, dqi_c__dz=0.0_dp, T_e=T_e)


      ! 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
      dFdz(i_r) = dr_dz(p=p, w_c=w, r_c=r, T_c=T, qd_c=q_d, qv_c=q_v, ql_c=q_l, qi_c=q_i, &
                        dql_c__dz=dFdz(i_ql), dqi_c__dz=dFdz(i_qi), dTc_dz=dFdz(i_T), dw_dz=dFdz(i_w))


   end function


   !> Compute the mixture density from the full equation of state
   pure function mod__mixture_density_from_eos(p, T, qd, qv, ql, qi) result(rho)
      real(dp), intent(in) :: p, T, qd, qv, ql, qi
      real(dp) :: rho
      real(dp) :: rho_inv

      rho_inv = (qd*R_d + qv*R_v)*T/p + ql/rho_l + qi/rho_i

      rho = 1.0_dp/rho_inv
   end function

   !> Compute the gas density from equation of state. This is the same as the
   !> full equation of state, but rearranged to have the form of a
   !> traditional ideal gas equation of state.
   pure function mod__cloud_gas_density_from_eos(p, T_c, qd_c, qv_c) result(rho_g)
      real(dp), intent(in) :: p, T_c, qd_c, qv_c
      real(dp) :: rho_g

      real(dp) :: R_s

      R_s = (R_v*qv_c + R_d*qd_c)/(qd_c + qv_c)

      rho_g = p/(T_c*R_s)
   end function

   !> Momentum equation
   !>
   !>     State variables:
   !>         w: verticaly velocity
   !>         r: cloud radius
   !>         rho_c: cloud density
   !>
   !> rho_e: environment density
   pure function dw_dz(p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, rho_e)
      ! TODO: remove this when I've found out how to get f2py to use module variables
      integer, parameter :: dp = selected_real_kind(12)

      real(dp), intent(in) :: p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, rho_e

      real(dp) :: dw_dz
      real(dp) :: rho_c, B, mu

      rho_c = mod__mixture_density_from_eos(p=p, T=T_c, qd=qd_c, qv=qv_c, ql=ql_c, qi=qi_c)

      B = (rho_e - rho_c)/rho_c

      mu = m__beta/r_c

      dw_dz = 1.0_dp/w_c * (g/(1.0_dp+m__gamma)*B - mu*w_c**2.0_dp - m__D*w_c**2.0_dp/r_c)

   end function


   !> Temperature (energy) equation
   !>
   !>      State variables:
   !>          qd_c, qv_c, ql_c, qi_c: in-cloud dry air, water vapour, liquid water, ice
   !>          T_c: in-cloud absolute temp
   !>
   !>          dql_c__dz, qdi_c__dz: vertical in
   !>
   !>          qd_c, qv_c, ql_c, qi_c: environment (constant) dry air, water vapour, liquid water, ice
   !>          T_e: environment absolute temp
   pure function dT_dz(r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz, T_e)
      real(dp), intent(in) :: r_c, T_c, qd_c, qv_c, ql_c, qi_c, dql_c__dz, dqi_c__dz, T_e

      real(dp) :: c_em_p, c_cm_p
      real(dp) :: qd_e, qv_e, ql_e, qi_e
      real(dp) :: Ds, mu

      real(dp) :: dT_dz

      ! heat capacity of cloud mixture
      c_cm_p = cp_d*qd_c + cp_v*(qv_c + ql_c + qi_c)

      ! heat capacity of environment mixture
      ! qd_e, qv_e, ql_e, qi_e = self.qd_e, self.qv_e, self.ql_e, self.qi_e
      ! TODO: Allow definition of moist environment
      qd_e = 1.0_dp
      qv_e = 0.0_dp
      ql_e = 0.0_dp
      qi_e = 0.0_dp

      c_em_p = cp_d*qd_e + cp_v*(qv_e + ql_e + qi_e)

      ! difference in environment and cloud moist static energy
      ! XXX: There appears to be something wrong with the formulation that
      ! includes liquid and ice, use just the liquid formulation for now
      ! Ds = (c_em_p*T_e - ql_e*L_v - qi_e*L_s)\
      ! -(c_cm_p*T_c - ql_c*L_v - qi_c*L_s)

      Ds = (c_em_p*T_e + ql_e*L_v) &
          -(c_cm_p*T_c + ql_c*L_v)

      mu = m__beta/r_c

      dT_dz = -g/c_cm_p + mu*Ds/c_cm_p + L_v/c_cm_p*dql_c__dz + L_s/c_cm_p*dqi_c__dz
   end function


   !> Mass conservation equation
   pure function dr_dz (p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c,  dql_c__dz, dqi_c__dz, dTc_dz, dw_dz)
      real(dp), intent(in) :: p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c,  dql_c__dz, dqi_c__dz, dTc_dz, dw_dz
      real(dp) :: dr_dz

      real(dp) :: dqd_c__dz
      real(dp) :: dqv_c__dz
      real(dp) :: rho_c
      real(dp) :: rho_cg
      real(dp) :: Rs_c
      real(dp) :: qg_c
      real(dp) :: mu

      dqd_c__dz = 0.0_dp
      dqv_c__dz = 0.0_dp
      rho_c = 0.0_dp
      rho_cg = 0.0_dp
      Rs_c = 0.0_dp
      qg_c = 0.0_dp
      mu = 0.0_dp

      ! XXX: for now assume that specific concentration of in-cloud dry does not change
      dqd_c__dz = 0.0_dp

      ! XXX: assume that changes in water vapour are exactly opposite to the increase in liquid water and ice
      dqv_c__dz = - dql_c__dz - dqi_c__dz

      ! in-cloud mixture density
      rho_c = mod__mixture_density_from_eos(p=p, T=T_c, qd=qd_c, qv=qv_c, ql=ql_c, qi=qi_c)

      ! in-cloud gas density
      rho_cg = mod__cloud_gas_density_from_eos(p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c)

      ! effective gas constant
      Rs_c = (R_v*qv_c + R_d*qd_c)/(qv_c + qd_c)

      ! in-cloud specific constant of gas constituents
      qg_c = qd_c + qv_c

      mu = m__beta/r_c

      dr_dz = r_c/2.0_dp*( &
         qg_c*rho_c/rho_cg * rho_c/rho_cg * g/(Rs_c*T_c) &
         + qg_c*rho_c/rho_cg * 1.0_dp/T_c*dTc_dz &
         + rho_c/(rho_cg*Rs_c) * (dqv_c__dz*R_v + dqd_c__dz*R_d) &
         + rho_c/rho_i*dqi_c__dz &
         + rho_c/rho_l*dql_c__dz &
         + mu & !  1./M*dM_dz
         - 1.0_dp/w_c*dw_dz &
         )

   end function

end module spec_conc_eqns


module mo_ccfm_parameter

  implicit none
  !> Definition double precision used with ECHAM
  integer,  parameter :: dp = selected_real_kind(12,307)  

  ! model parameters
  !> maxium cloud radius [m]
  real(dp), parameter :: r_min = 200._dp
  !> number of different cloud types
  integer,  parameter :: clno  = 10
 
  ! switches
  logical,  parameter :: emic  = .true.           ! ECHAM5 microphysics parametrisation (instead of Kessler Scheme)
  logical,  parameter :: shear = .false.          ! correction of v_z due to wind shear

  ! empirical parameters
  real(dp), parameter :: C_mu  =  0.2_dp          ! entrainment rate factor [-]

  ! physical and mathematical constants
  real(dp), parameter :: pi    = 3.14159265358_dp !                                             [-]
  real(dp), parameter :: g     = 9.80665_dp       ! acceleration of gravity                     [m/s^2]
  real(dp), parameter :: Tnull = 273.16_dp        ! absolute, zero not wodka                    [K]
  real(dp), parameter :: Rd    = 287.05_dp        ! gas constant for dry air                    [J/(kg*K)]
  real(dp), parameter :: Rv    = 461.51_dp        ! gas constant for water vapour               [J/(kg*K)]
  real(dp), parameter :: Cpd   = 1005.46_dp       ! dry air specific heat capacity at const. p. [J/(kg*K)]
  !> latent heat of evaporation for water [J/kg]
  real(dp), parameter :: aLv   = 2.5008e6_dp      ! latent heat of vaporisation (at 0 C)        [J/kg]
  !> latent heat of sublimation for water [J/kg]
  real(dp), parameter :: aLs   = 2.8345e6_dp      ! latent heat of sublimation  (at 0 C)        [J/kg]
  real(dp), parameter :: aLf   = aLs - aLv        ! latent heat of fusion       (at 0 C)        [J/kg]  
  real(dp), parameter :: Rho_w = 1000._dp         ! density of water                            [kg/m^3]
  real(dp), parameter :: Rho_0 = 1.3_dp           ! air density at surface                      [kg/m^3]
  real(dp), parameter :: Rho_i = 500._dp          ! density of cloud ice                        [kg/m^3]
  real(dp), parameter :: Rho_s = 100._dp          ! density of snow                             [kg/m^3]
  !> Ratio of gas constants for dry air and water vapor
  real(dp), parameter :: Epsi  = Rd / Rv          ! (molecular weight of water vapour)/(molecular weight of dry air) [-]

end module mo_ccfm_parameter

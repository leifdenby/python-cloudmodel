
module mo_ccfm_cloudbase

  use mo_ccfm_parameter, only: dp, Tnull, Rv, Rd, cpd, aLv, aLs, epsi
  
  private
  
  real(dp), parameter :: BASEBUO = 0.5_dp

  real(dp), parameter :: C1ES  = 610.78_dp
  real(dp), parameter :: C2ES  = C1ES*RD/RV 
!> "a" constant for saturation vapor pressure over liquid water in Tetens equation
  real(dp), parameter :: C3LES = 17.269_dp 
!> "b" constant for saturation vapor pressure over liquid in Tetens equation
  real(dp), parameter :: C4LES = 35.86_dp
!> "a" constant for saturation vapor pressure over ice in Tetens equation
  real(dp), parameter :: C3IES = 21.875_dp
!> "b" constant for saturation vapor pressure over ice in Tetens equation
  real(dp), parameter :: C4IES = 7.66_dp 

  real(dp), parameter :: C5LES = C3LES*(TNULL-C4LES)
  real(dp), parameter :: C5IES = C3IES*(TNULL-C4IES)
  
  public cloudbase

  ! required for f2py:
  public moist_adjust, lua, lub, luc


contains

!> Given an input profile describing the atmospheric state compute the cloud
!> base height, temperature and mixing ratio.
!>
!> The cloud base is defined as cumulus condensation level. It is determined by
!> lifting an air parcel of environmental air from the surface until condensation
!> occurs and the air parcel reaches a level where the virtual temperature of the
!> lifted air parcel is maximally 0.5 K lower than the environmental temperature,
!> which means that an initial buoyancy perturbation for the lifted air parcel is
!> assumed.
!!
!! @param level number of levels in host model
!! @param T_e temperature profile [K]
!! @param q_e relative humidity profile [kg/kg]
!! @param geo geopotential height * g (unclear) in profile [m]
!! @param prs pressure profile [Pa]
!! @return cl_base height of cloudbase [profile height index]
!! @return T_c temperature at cloud base [K]
!! @return q_c water mixing ratio at cloud base [?]
subroutine cloudbase ( level,   &
                       T_e,     &
                       q_e,     &  
                       geo,     & 
                       prs,     &
                       cl_base, & 
                       T_c,     & 
                       q_c)
 
  implicit none      

  integer,                    intent(in)  :: level     ! no vert levels [-]
  real(dp), dimension(level), intent(in)  :: prs       ! pressure
  real(dp), dimension(level), intent(in)  :: geo       ! geo pot height 
  real(dp), dimension(level), intent(in)  :: T_e       ! environmental temperature
  real(dp), dimension(level), intent(in)  :: q_e       ! environmental mix ratio 

  integer,                    intent(out) :: cl_base   ! cloud base
  real(dp),                   intent(out) :: T_c       ! cloud temperature
  real(dp),                   intent(out) :: q_c       ! cloud mix ratio
  
  real(dp)                                :: buo
  real(dp)                                :: q_old
  integer                                 :: l
  logical                                 :: cond
  
  ! init
  cl_base = 1
  cond    = .false.

  T_c = T_e(level)
  q_c = q_e(level)

  base: do l = level-1, 2 ,-1
        
     ! dry adiabatic ascent
     T_c =  ( cpd * T_c + geo(l+1) - geo(l) ) / cpd
     q_old = q_c

     ! adjustmend for moist ascent
     call moist_adjust( T_c, q_c, prs(l) )
 
     ! if condensation occurred
     if ( q_c < q_old ) cond = .true.

     buo =   T_c    * ( 1._dp + epsi * q_c )     &
           - T_e(l) * ( 1._dp + epsi * q_e(l) )
    
     if ( cond .and. (buo >= - basebuo) ) then

         cl_base = l
         return
            
     end if
        
     if ( cond .and. (buo <= -(basebuo + 0.2_dp)) ) return
        
  end do base

end subroutine cloudbase



!> Assuming contant pressure adjust temperature and saturation specific
!> concentration of water vapor.
!>
!> Saturation water vapor pressure is calculated using Tetens (1930) formula.
!!
!! @param tem initial temperature [K]
!! @param q_v initial specific concentration of water vapor [kg/kg]
!! @param prs pressure [Pa]
!! @returns tem temperature after adjustment [K]
!! @return q_v saturation specification concentration of water vapor [kg/kg]
subroutine moist_adjust( tem, q_v, prs )

  implicit none

  real(dp), intent(inout) :: tem, q_v
  real(dp), intent(in)    :: prs

  real(dp) :: q_sat ! saturation water vapour over water or ice
  real(dp) :: cor
  real(dp) :: cond  ! condensate
  
 
  q_sat = lua( max( 180._dp, min( tem, 370._dp))) / prs
    
  q_sat = min( 0.5_dp, q_sat )
  cor   = 1._dp / ( 1._dp - ( Rv/Rd - 1._dp ) * q_sat )
  q_sat = q_sat * cor
  
  cond  = (q_v - q_sat) /                                                &
          (1._dp + q_sat * cor * lub( max(180._dp,min( tem ,370._dp))))
  
  cond  =  max( cond, 0._dp )

  ! increase temperature through latent heat of condensation
  tem = tem + luc(tem) * cond
      
  ! remove condensed water from water vapor
  q_v  = q_v - cond

  if ( cond > 0._dp ) then
  
     ! calculate specific humidity from saturation water vapor pressure using Tetens formula
     q_sat = lua( max( 180._dp, min( tem, 370._dp))) / prs
    
     q_sat = min( 0.5_dp, q_sat )
     cor   = 1._dp / ( 1._dp - ( Rv/Rd - 1._dp ) * q_sat )
     q_sat = q_sat * cor
  
     cond  = (q_v - q_sat) /                                             &
             (1._dp + q_sat * cor * lub( max(180._dp,min(tem,370._dp))) )
  
     cond  =  max( cond, 0._dp )

     tem = tem + luc(tem) * cond
      
     q_v  = q_v - cond
     
  end if

end subroutine moist_adjust

real(dp) function lua(T)
! lua(T) = es(T) * Rd/Rv
! in echam: tlucua

  implicit none

  real(dp) :: T
  real(dp) :: CVM3, CVM4 

  if ( T > Tnull ) then
     CVM3 = C3LES
     CVM4 = C4LES
  else
     CVM3 = C3IES
     CVM4 = C4IES
  end if

  lua = c2es * exp( CVM3 * (T - TNULL)/(T - CVM4)) 
  
end function lua


real(dp) function lub(T)

  implicit none
    
  real(dp)       :: T
  real(dp)       :: CVM4, CVM5
            
  if ( T > Tnull ) then
     CVM4 = C4LES
     CVM5 = C5LES * ALV / CPD
  else
     CVM4 = C4IES
     CVM5 = C5IES * ALS / CPD
  end if

  lub = CVM5 / ( T - CVM4 )**2

end function lub


!> Calculate temperature dependent specific temperature change due to phase-transition.
!> 
!> If temperature is below freezing latent heat of freezing is used, otherwise
!> latent heat of condensation is used.
!! @param T temperature at which to calculate latent heat release [K]
!! @return specific temperature change [K*kg/kg]
real(dp) function luc(T)
! latent heat

  implicit none

  real(dp) :: T
      
  if ( T > Tnull ) then
     ! latent heat of evaporation / dry air spec. heat cap. at contant pressure
     ! J/kg / ( J/kg/K ) = [K*kg/kg]
     luc = ALV / CPD
  else
     luc = ALS / CPD
  end if

end function luc


end module mo_ccfm_cloudbase

module ccfm_microphysics

contains

!> ECHAM microphysics scheme
!> cp. Zhang and Lohmann, 2004
! @todo Find equations and verify this
 !subroutine ecmic(derive, x, y, dx, pe, q_ve, rho,  mu, q_satw, q_sati, dpdz, cdnc, land)                           

   !use mo_ccfm_parameter, only: pi, g, epsi, alv, als, tnull, rd,  &
                                !rho_w, rho_0, rho_i, rho_s
   !use mo_ccfm_clouddata
   
   !implicit none

   !!type(t_cldata), intent(inout) :: derive  ! derivation of cloud data    
   !!type(t_cldata), intent(in)    :: y       ! cloud data      
   !real(dp),       intent(in)    :: x       ! height [m] - only for debugging
      
   !real(dp),       intent(in)    :: dx      ! stepsize (or stephight)
   !real(dp),       intent(in)    :: pe      ! pressure                             [Pa]
   !real(dp),       intent(in)    :: q_ve    ! mixing ratio environment             [kg/kg]      
   !real(dp),       intent(in)    :: q_satw  ! saturation mixing ratio (over water) [kg/kg] 
   !real(dp),       intent(in)    :: q_sati  ! saturation mixing ratio (over ice)   [kg/kg]
   !real(dp),       intent(in)    :: dpdz    ! pressure derivative                  [Pa/m]
   !real(dp),       intent(in)    :: rho     ! dry air density                      [kg/m^3]
   !real(dp),       intent(in)    :: mu      ! entrainment rate

   !real(dp),       intent(in)    :: cdnc    ! cloud droplet number concentration   [1/m^3]
   !logical,        intent(in)    :: land    ! land sea mask                        [-]
   
   !real(dp) :: T_c    ! cloud Temp                           [K]      
   !real(dp) :: Q_cond ! condensation rate                    [kg/kg/m]
   !real(dp) :: Q_dep  ! deposition rate                      [kg/kg/m]
   !real(dp) :: Q_ifrl ! freezing rate (cloud water to ice)   [kg/kg/m]
   !real(dp) :: Q_mlt  ! melting rate                         [kg/kg/m]
   !real(dp) :: Q_raul ! autoconversion rate (liquid -> rain) [kg/kg/m]
   !real(dp) :: Q_sagi ! aggregation (ice crystals to snow)   [kg/kg/m]
   !real(dp) :: Q_racl ! accretion (cloud water to raindrops) [kg/kg/m]
   !real(dp) :: Q_sacl ! accretion (cloud water to snow)      [kg/kg/m]
   !real(dp) :: Q_saci ! accretion (ice crystals to snow)     [kg/kg/m]
   !real(dp) :: Q_falr ! fallout of rain                      [kg/kg/m]
   !real(dp) :: Q_ifrr ! freezing of rain to cloud ice        [kg/kg/m]     
   !real(dp) :: Q_iacr ! accretion (rain to cloud ice)        [kg/kg/m]

   !real(dp) :: Q_detv ! detrainment rate water vapour        [kg/kg/m]
   !real(dp) :: Q_detl ! detrainment rate cloud water         [kg/kg/m]
   !real(dp) :: Q_deti ! detrainment rate cloud ice           [kg/kg/m]
   !real(dp) :: Q_detr ! detrainment rate rain                [kg/kg/m]
   !real(dp) :: Q_dets ! detrainment rate snow                [kg/kg/m]

   !real(dp) ::  Q_dilvl = 0._dp
   !real(dp) ::  Q_dilvi = 0._dp

   !real(dp) :: Q_ifrl_hetstat ! heterogeneous / statistical freezing rate   [kg/kg/m] 
   !real(dp) :: Q_ifrl_contact ! contact freezing rate        [kg/kg/m]


   !real(dp) :: a_1, a_2, a_4, a_5, a_6   ! empirical constants
   !real(dp) :: b_1, b_3, b_4, b_5        ! empirical constants
   !real(dp) :: c_1                       ! empirical constants
   !real(dp) :: dt_1
   !real(dp) :: E_ii, E_ll, E_ri, E_rl, E_sl, E_si  ! collection efficiencies
   !real(dp) :: f_1
   !real(dp) :: g_1, g_2
   !real(dp) :: gamma_2, gamma_3            ! tunable constants
   !real(dp) :: D_b
   !real(dp) :: M_i
   !real(dp) :: N_lx, n_0s, n_0r, n_a       ! concentration number [1/m^3]
   !real(dp) :: mue 
   !real(dp) :: r_ei, r_iv, r_so, r_lv, r_m ! hydrometeor radii [m]
   !real(dp) :: lambda_s, lambda_r 
   !real(dp) :: p_satw                      ! saturation water vapour pressure [kg/kg]
   !real(dp) :: q_l_thres                   ! autoconversion thershold [kg/kg/]
   !real(dp) :: XX
   
   !! temp
   !integer  :: time
   !type(t_cldata) :: y_tmp
   !real(dp) :: sq_xl, fq_l, &
               !sq_xi, fq_i


!!------------------------------- initialisation ------------------------
   !T_c = y%tem

   !Q_cond  = 0._dp
   !Q_dep   = 0._dp 
   !Q_mlt   = 0._dp
   !Q_raul  = 0._dp
   !Q_sagi  = 0._dp
   !Q_racl  = 0._dp
   !Q_sacl  = 0._dp
   !Q_saci  = 0._dp
   !Q_ifrr  = 0._dp
   !Q_iacr  = 0._dp
   !Q_dilvl = 0._dp
   !Q_dilvi = 0._dp
   !Q_detl  = 0._dp
   !Q_deti  = 0._dp
   !Q_falr  = 0._dp
   !Q_ifrl  = 0._dp
   !Q_ifrl_hetstat = 0._dp
   !Q_ifrl_contact = 0._dp
   !Q_detv  = 0._dp
 
   !a_1    = 1000._dp         ! empirical constant [1/(s*m^3)] 
   !a_2    =  700._dp         ! empirical constant [1/s]  
   !a_4    =    4.83_dp       ! empirical constant
   !a_5    =   83.81_dp       ! empirical constant [/]
   !a_6    =    2.115_dp      ! 2115._dp*.01_dp**.2_dp ! empirical constant - 2115 cm^(1-b_3)/s = 841.996666
      
   !b_1    = 0.66_dp          ! empirical constant [1/K]  
   !b_3    = 0.8_dp           ! empirical constant [-]
   !b_4    = 0.25_dp          ! empirical constant
   !b_5    = 0.216_dp         ! empirical constant [-]      
      
   !E_ii   = 0.1_dp           ! collection efficiency between ice crystals
   !E_ri   = 1._dp            ! collection efficiency of rain for cloud ice
   !E_sl   = 1._dp            ! collection efficiency of snow for cloud droplets
  
   !g_1    = 2.54926_dp       ! Gamma_function(3 + b_4)
   !g_2    = 496.606078_dp    ! Gamma_function(6 + b_3)
  
   !M_i    = 4.19e-13_dp      ! mass of a single ice crystal [kg] 
    
   !n_0s   = 3.e6_dp          ! intercept parameter of the snow size distribution [1/m^4]
   !n_0r   = 8.e6_dp          ! intercept parameter of the snow size distribution [1/m^4]

   !r_so   = 1.e-4_dp         ! smallest radius of particle in snow class [m]

   !XX      = 0.25_dp         ! dispersion of fall velocity 

   !gamma_2 = 95._dp          ! efficiency of snow formation - tunable constant   
   !gamma_3 = 0.1_dp          ! tunable constant

!!-------------------------- end: initialisation -------------------------------

!!==== warm cloud  processes =========================================================
   !! condensation = change of sat vap rat due to temp and pressure change:
   
   !! hydrostatic assumption 
   !!Q_cond  = - q_satw / Rd / T_c  * (epsi * ALV / T_c * derive%tem + g)

   !! if pressure change is not taken into account:
   !!Q_cond  = - q_satw * ( epsi * aLv / Rd / y%tem/y%tem * derive%tem)
   
   !!if dpdz is given ( instead of hydrostatic assumption: g/(Rd*T) = - dp/dz * 1/p ):
   !Q_cond  = - q_satw * ( epsi * ALV / Rd / T_c / T_c * derive%tem - dpdz / pe )
 
   !!q_swnext = q_satw(y%tem + derive%tem * dx, pe + dpdz * dx)
   !!Q_cond  = ( q_sw - q_swnext ) / dx
 
   
   !if ( y%q_l > 1.e-15_dp ) then
         
      !!--- autoconversion (collision and coalesence) - of cloud water to rain --!
      
      !! Beheng and Berry relatively similar - Berry a bit better
      
      !!--- Beheng1994
      !! 6e28 * n**(-1.7) = 1.2e27   , where n=10
      !N_lx   = (1.e-6_dp*cdnc)**(-3.3_dp)
      !Q_raul = 15._dp * 1.2e27_dp * N_lx * (1.e-3_dp*rho*y%q_l)**(4.7_dp) / rho
          
      !!--- Lin1983 (a modified Berry1968 version)
      !!q_l_thres = 0.002_dp ! autoconversion threshold [kg/kg]
      !!if ( y%q_l > q_l_thres ) then                           
      !!    Q_raul = 1.e-3_dp * (rho * (y%q_l - q_l_thres)**2)              &
      !!             / ( 1.2e-4_dp + ( 1.569e-18_dp * cdnc                  &
      !!                              / (0.15_dp * rho * (y%q_l-q_l_thres))))
      !!end if

      !!--- Berry1986
      !D_b = 0.366_dp              ! maritime    cloud
      !if ( land )  D_b = 0.146_dp ! continental cloud
          
      !Q_raul = y%q_l * y%q_l * rho * 1.e3_dp                        &
               !/ (60._dp * (5._dp + .0366_dp * cdnc * 1.e-9_dp      &
                                    !/ (y%q_l * rho * D_b) ) )

      !!--- Ferrier1994 (also modified Berry1968 with cdnc dependend threshold)
      
      
      !!--- Tripoli/Cotton 1980
      !!q_l_thres = 0.0005_dp  ! autoconversion threshold [kg/kg]
      !!if ( y%q_l > q_l_thres ) then
      !!   E_ll   = 0.55_dp       ! collection efficiency
      !!   ! dynamic viscosity [kg/m/s]
      !!   mue    = 1.72e-5_dp * (393._dp / (y%tem + 120._dp))        &
      !!            * ( y%tem / 273._dp )**1.5_dp 
      !! 
      !!   Q_raul = 0.104 * g * E_ll * rho**1.33_dp                   &
      !!            / ( mue * (cdnc * rho_w)**0.33_dp)                &
      !!            * y%q_l**2.33_dp * ( y%q_l - q_l_thres )
      !!end if
           
      !! kg/kg/s -> kg/kg/m
      !Q_raul = Q_raul / y%v_z

      !!----- accretion - of cloud water by rain ------------------------------!
      !if ( y%q_r > 1e-15 ) then
  
       !! Beheng1994: (works ok)
       !Q_racl = 6._dp * rho * y%q_l * y%q_r 
       
       !! Tripoli1980: (works good for 1995 - works with Lin for 1995/1997 not for darwin2005)
       !!r_m  = 2.7e-4_dp ! characteristic drop radius [m] 
       !!E_rl = 1.0_dp    ! collection efficiency of rain for cloud water
       !!Q_racl = 0.884_dp * E_rl * sqrt( g * rho / (rho_w * r_m) ) * y%q_l * y%q_r 

       !! [kg/kg/s] -> [kg/kg/m]
       !Q_racl = Q_racl / y%v_z
   
      !end if ! q_r > 0 - collection
      
   !end if  ! q_l >  0  - autoconversion + collection
    
!!= mixed phase cloud and ice cloud processes ================================
   !mixed_cl: if ( Tnull > T_c  ) then
            
      !!-------------------  deposition  ------------------------------------!
                     
      !if ( y%q_i > 5.e-7_dp .or. T_c < 261._dp ) then    ! Bergeron-Findeisen - not really
        
          !!Q_dep  = - q_sati / RD / T_c * (epsi * ALS / T_c * derive%tem + g)

          !Q_dep  = - q_sati * ( epsi * ALV / Rd / T_c / T_c * derive%tem - dpdz / pe )
          
          !Q_cond = 0._dp
          
      !end if
    
            
      !!----- freezing of cloud droplets ---------------------------------!
         
      !if ( y%q_l > 1.e-15_dp ) then
         
           !!-- heterogeneous freezing rate (cloud water to ice) [kg/kg/s] ----!
           !! Murakami1990, Levkov1992, Ferrier1994, Lohmann1996
           !Q_ifrl_hetstat =  a_1 * (exp( b_1 * (Tnull - T_c) ) - 1._dp)        &
                            !* rho / (Rho_w * cdnc) * y%q_l*y%q_l
            
           !!-- contact       freezing rate (cloud water to ice) [kg/kg/s] ----!
           !! Levkov1992, Cotton1986, Lohmann1996
            
           !! mean volume cloud droplet radius                  
           !r_lv           = ( 0.75_dp * y%q_l * rho                           &
                              !/ ( pi * Rho_w * cdnc) )**(1._dp/3._dp)
              
           !! active contact nuclei [1/m^3] - Cotton1986
           !! tmw30: that's rather a lot than rather not
           !n_a            = (Tnull - 3._dp - T_c ) * 2.e5_dp 
            
           !! ice nuclei -  Meyers1992, PruppacherKlett1997
           !!n_a = exp( 0.262_dp * (Tnull - T_c) - 2.8_dp)  
             
           !! ice nuclei -  Fletcher1962, PruppacherKlett1997 
           !! ic = 1.0e-2_dp * ( 0.6_dp * (Tnull - T_c) )

           !f_1             = 4._dp * pi * r_lv * cdnc * n_a / rho

           !f_1             = max( 0._dp, f_1 )
         
           !Q_ifrl_contact = 1.4e-20_dp * f_1
            
           !Q_ifrl = Q_ifrl_hetstat + Q_ifrl_contact

           !! [kg/kg/s] -> [kg/kg/m]
           !Q_ifrl = Q_ifrl / y%v_z
          
      !end if 
            
            
      !!--- aggregation  - cloud ice to snow [kg/kg/m] ----------------------!
      !if ( y%q_i > 1.e-15_dp ) then
   
            !c_1   = y%q_i * rho * a_2 * E_ii * XX / Rho_i * (Rho_0/rho)**(0.33_dp)
            
            !! effective radius of ice crystals [10e-6*m]
            !r_ei  = a_5 * (1.0e3_dp * rho * y%q_i)**b_5
            !r_ei  = min( max(r_ei, 10._dp), 150._dp)
            
            !! mean volume ice crystal radius [m]
            !r_iv  = 1.0e-6_dp * &
                   !(sqrt(2809._dp * r_ei**3 +  5113188._dp) - 2261._dp)**(1/3._dp)

            !dt_1  = - 6._dp / c_1 * log10(r_iv/r_so) ! log_10((x)^3) = 3 * log_10(x)
           
            !Q_sagi = gamma_2 * y%q_i / dt_1
            
            !Q_sagi = Q_sagi/ y%v_z
            
      !end if


      !!-- accretion - of cloud ice and cloud water by snow  [kg/kg/m] ------!
      !if ( y%q_s > 1.0e-15_dp ) then   ! if snow exist (to collect ice and water)
            
            !! slope of particle distribution
            !lambda_s = ( PI * Rho_s * n_0s / rho / y%q_s )**0.25_dp
            !lambda_s = lambda_s**(3._dp+b_4)

            !if ( y%q_l > 1.e-15_dp) then

             !Q_sacl   = gamma_3 * PI * E_sl * n_0s * a_4 * y%q_l * g_1    &
                        !/ (4.0_dp * lambda_s) * ( Rho_0 / rho )**0.5_dp

             !Q_sacl = Q_sacl / y%v_z

            !end if
            
            !if ( y%q_i > 1.e-15_dp) then 
             !! collection efficiency of snow for ice crystals
             !E_si   = exp( 0.025_dp * (T_c - TNULL) ) 
            
             !Q_saci =           PI * E_si * n_0s * a_4 * y%q_i * g_1    &
                        !/ (4._dp * lambda_s) * ( Rho_0 / rho )**0.5_dp
            
             !Q_saci = Q_saci / y%v_z
            !end if
            
      !end if
         
         
      !!-- freezing of rain ---!
      !if ( y%q_r > 1e-15_dp ) then
         
            !! slope of rain particle distribution
            !lambda_r = ( PI * Rho_w * n_0r / rho / y%q_r )**0.25_dp

            !!if ( lambda_r < 1.e5_dp ) then      

            !!   lambda_r = lambda_r**7._dp

            !!   freezing of raindrops (should be hail but is here considered as ice) [kg/kg/]
            !!
            !!   (cp. LinFarleyOrville83)   - tmw210706: hier geht's los   
            !!    Q_ifrr   = a_1 * (exp( b_1 * (TNULL - T_c) ) - 1) / lambda_r &
            !!                * 20 * PI*PI * n_0r * Rho_w / rho 
               
            !!   lambda_r = lambda_r**(6._dp + b_3)
               
            !!end if
            
      !end if
     
   !end if mixed_cl
   

   !!------------------------------ dilution ----------------------------------!

   !Q_detl = mu / (1._dp + mu) * y%q_l
   !Q_deti = mu / (1._dp + mu) * y%q_i
   !Q_detr = mu / (1._dp + mu) * y%q_r
   !Q_dets = mu / (1._dp + mu) * y%q_s
 
   !! dilution effect (entrainment of dry air) on cloud liquid water or ice
   !Q_dilvl = mu * ( q_satw - q_ve ) 
   !Q_dilvi = mu * ( q_sati - q_ve ) 
 
   !if ( y%q_i > 1.e-5_dp  ) then
      !Q_dilvl = 0._dp              ! if deposition - substract vapour dilution from ice
   !else
      !Q_dilvi = 0._dp              ! warm phase cloud - substract vapour dilution from liquid water
   !end if
   
   !sq_xl =          (Q_raul + Q_racl + Q_sacl + Q_ifrl)*dx
   !! future q_l, considering only cond and dilution
   !fq_l  = y%q_l  + (Q_cond - Q_detl - Q_dilvl)*dx

   !if ( fq_l < 0._dp ) then        
      !derive% v_z = -100._dp     ! stop cloudmodel calculation
   !else if ( sq_xl > fq_l ) then
      !fq_l = fq_l * 0.7_dp
      !Q_raul = Q_raul / sq_xl * fq_l
      !Q_racl = Q_racl / sq_xl * fq_l
      !Q_sacl = Q_sacl / sq_xl * fq_l
      !Q_ifrl = Q_ifrl / sq_xl * fq_l
   !end if
    
   !sq_xi =         (Q_sagi +Q_saci)*dx
   !! future q_l, considering only cond and dilution
   !fq_i  = y%q_i + (Q_dep -Q_deti -Q_dilvi +Q_ifrl)*dx

   !if ( fq_i < 0._dp ) then        
      !derive% v_z = -100._dp     ! stop cloudmodel calculation
   !else if ( sq_xi > fq_i ) then
      !fq_i = fq_i * 0.7_dp
      !Q_sagi = Q_sagi / sq_xi * fq_i
      !Q_saci = Q_saci / sq_xi * fq_i
   !end if
   
   !derive% q_l = -Q_raul -Q_racl -Q_sacl                  -Q_ifrl                 +Q_cond        -Q_detl -Q_dilvl 
   !derive% q_i =                           Q_ifrr +Q_iacr +Q_ifrl -Q_sagi -Q_saci         +Q_dep -Q_deti -Q_dilvi      
   !derive% q_r =  Q_raul +Q_racl          -Q_ifrr -Q_iacr                                        -Q_detr
   !derive% q_s =                  Q_sacl                          +Q_sagi +Q_saci                -Q_dets
   !derive% q_v =                                                                  -Q_cond -Q_dep -Q_dilvl-Q_dilvi
   
   
   !!--------  data check  -----------------------------------------------------!
   !y_tmp = y + derive * dx

   !if ( y%q_i       <  1.e-15_dp  .and. &      ! i.e. no depostion
        !y_tmp% q_l  <  0._dp      .and. &
        !derive% v_z > -100._dp          ) then
!#ifndef AUTONOM         
      !derive% v_z = -100._dp     ! stop cloudmodel calculation
      !time = get_time_step()
      !print '("exit: q_l < 0 (v -> -100), height: ",f10.1," temp ",f0.1," time ",i0," rad ",f0.0," v ",f0.2)', &
              !x, y%tem, time, y%rad, y%v_z
      !print *, 'q_l:', y%q_l,'cond:',Q_cond*dx, 'dilv:',-Q_dilvl*dx, 'dill:', -Q_detl*dx, '->q_r:',(-Q_raul-Q_racl)*dx
      !print *
!#endif
   !end if

   !if ( y_tmp%q_i < 0._dp     )  then
      !derive%v_z = - 10000._dp
!#ifndef AUTONOM
      !time = get_time_step()
      !print *, 'q_i:', y%q_i, 'dep:',Q_dep*dx,'dilv:',-Q_dilvi*dx, 'dili:', -Q_deti*dx, '->q_p:',(-Q_sagi-Q_saci)*dx, &
             !"fq_i",fq_i
      !print '("exit: q_i < 0 (v -> -10000), height: ",f10.1," temp ",f0.1," time ",i0," rad ",f0.0," v ",f0.2)', &
              !x, y%tem, time, y%rad, y%v_z
      !print *
!#endif
   !end if
   !!---- end: data check  -----------------------------------------------------!

!!-----------  test output -----------------------------------------------------!
!#ifdef AUTONOM
  !write(40,'(2f8.1, f8.2, 13es15.7)') x,dx,y%tem,Q_cond,Q_dep,Q_dilvl+Q_dilvi, &
            !Q_detl,Q_ifrl,Q_mlt,Q_raul,Q_sagi,Q_racl,Q_sacl,Q_saci,Q_ifrr,Q_iacr
   !write(41,'(2f8.1, f8.2, 8es11.3)') x,dx,y%tem, y%q_l, Q_cond*dx, Q_detv*dx, & 
            !Q_detl*dx,  Q_ifrl*dx, Q_raul*dx, Q_racl*dx, Q_sacl*dx
   !write(42,'(2f8.1, f8.2, 8es11.3)') x,dx,y%tem, y%q_i, Q_dep*dx,  Q_detv*dx, &
            !Q_ifrl*dx, Q_sagi*dx, Q_saci*dx, Q_ifrr*dx, Q_iacr*dx
   !write(43,'(2f8.1, f8.2, 6es11.3)') x,dx,y%tem, y%q_r, Q_raul*dx, Q_racl*dx, & 
            !Q_ifrr*dx, Q_iacr*dx, Q_falr*dx
   !write(44,'(2f8.1, f8.2, 4es11.3)') x, dx, y%tem, y%q_s, Q_sagi*dx,          & 
            !Q_sacl*dx, Q_saci*dx
!#endif
!!------ end: test output ------------------------------------------------------!

!end subroutine ecmic

subroutine kessler(tem, q_r, q_l, v_z, rho, mu, pe, q_ve, q_satw, dpdz, cdnc, land, dTdz,   &
                   Q_conv, Q_coll, Q_entr, Q_cond)

   use mo_ccfm_parameter
             
   implicit none 
      
   real(dp),       intent(in)  :: rho    ! dry air density [kg/m^3]
   real(dp),       intent(in)  :: mu     ! entrainment rate
   real(dp),       intent(in)  :: pe     ! pressure [Pa]
   real(dp),       intent(in)  :: q_ve   ! environment water vapour mixing ratio [kg/kg] 
   real(dp),       intent(in)  :: q_satw ! saturation  water vapour mixing ratio [kg/kg] 
   real(dp),       intent(in)  :: dpdz   ! pressure derivative                   [Pa/m]
   
   real(dp),       intent(in)  :: cdnc    ! cloud droplet n c        [1/m^3]
   logical,        intent(in)  :: land   ! land sea mask            [-]  

   real(dp),       intent(in)  :: dTdz   ! change of temp [K/m]
    
   real(dp),       intent(out) :: Q_conv ! conversion rate 
   real(dp),       intent(out) :: Q_coll ! collision rate
   real(dp),       intent(out) :: Q_cond ! condensation rate
   real(dp),       intent(out) :: Q_entr ! entrainment
   
   real(dp) :: K_1    ! conversion const [1/s]
   real(dp) :: K_2    ! collection const [1/s]
   real(dp) :: K_3    !       
   real(dp) :: K_a    ! autoconversion threshold [kg/kg]     
   real(dp) :: D_0    ! sqrt(D_0), where D_0 = median volume drop diameter 
   real(dp) :: v_t    ! terminal vertical velocity (of falling rain drops) [m/s] 
   real(dp) :: D_b    ! droplet relative dispersion [/]


   ! previously in `y`:
   real(dp), intent(in) :: tem ! temperature
   real(dp), intent(in) :: q_r
   real(dp), intent(in) :: q_l
   real(dp), intent(in) :: v_z

   ! autoconversion threshold [kg/kg]  -  cp. (WeinsteinMacCready 1969):
   K_a = 0.0005

   ! setting Kessler constants - cp. WeinsteinMacCready1969 p.943
   if ( tem > 258.15_dp ) then   ! if temperature > -15 C
                                   !  i.e. before ice nucleation occurs
      K_1 = 0.00075_dp             ! autoconversion const [1/s]
      K_2 = 0.0052_dp              ! collection const [1/s]
      K_3 = 15.39_dp               ! terminal velocity constant
   else                            ! if ice nucleation occurs, temp < -15 C
      K_1 = 0.0015_dp
      K_2 = 0.00696_dp
      K_3 = 11.58_dp
   end if
        
   ! D_0 = sqrt( median volume drop diameter ) (cp.  WeinsteinMcCready69 p. 940 + 943)
   !    q_r/rho = rain water [kg/m^3] 
   D_0 = K_3 * (1.e3_dp * max(0._dp,q_r))**0.125_dp
   ! terminal vertical velocity of falling rain relative to updraft		
   v_t = - 130._dp * D_0
   v_t = - K_3 * (max(0._dp,q_r))**0.125_dp  

   ! (auto)conversion rate: (dQ_h/dz)_conv = - dQ_c/dz [kg/kg/m]
   !  change of cloud liquid water (Q_c) -> rain water (hydrometeor water) (Q_h)
   Q_conv = 0._dp

   ! collision rate: (dQ_h/dz)_coll [kg/kg/m]
   Q_coll = 0._dp
      
   !---- autoconversion and collision
   if ( q_l > 1.e-20_dp) then
   
      ! Berry scheme [Simpson, Wiggert (1969)]
    
      if ( land ) then  ! i.e. continental cloud
            D_b = 0.146_dp
      else                 ! i.e. maritime    cloud
         D_b = 0.366_dp
      end if

      Q_conv = q_l * q_l * rho * 1.e3_dp / v_z                 &
                / (60._dp * (5._dp + .0366_dp * cdnc * 1.e-9_dp      &
                                 / (q_l * rho * D_b) ) )

     ! Kessler scheme
     !if ( y%q_l * rho > K_a ) then 	! Kessler scheme
         
         ! third summand of equation (4a) in [Weinstein, MacCready (1969)]
         ! threshold value K_a as in [Weinstein, MacCready (1969)], 
         ! i.e. K_a = 0.0005 [kg/kg] or [kg/m^3] (as in [Simpson, Wiggert (1969)])
     !    Q_conv = K_1 * ( y%q_l - K_a ) / y%v_z		! [kg/kg/m]
          
     !end if
      
      ! collection rate - accretion of cloud droplets by rain
      ! fourth summand of equation (4a) in [Weinstein, MacCready (1969)]
      Q_coll = K_2/(v_z - v_t) *  rho**(-0.875_dp)              &
               * max(0._dp, q_l) * max(0._dp, q_r)**(0.875_dp)

   end if

   ! entrainment rate: (dQ_h/dz)_entr [kg/kg/m]    
   !   second summand of equation (4a) in [Weinstein, MacCready (1969)]
   Q_entr = mu * ( q_satw - q_ve + q_l )       ! [kg/kg/m]
      
   ! change due to condesation (lifting decreases temp + press): 
   !   (dQ_h/dz)_lift [kg/kg/m]
   !  first summand of equation (4a) in [Weinstein, MacCready (1969)]
   
   !Q_cond  = - q_satw / Rd / T_c  * (epsi * ALV / T_c * derive%tem + g)
   
   ! often, pressure change is not taken into account, in this case:
   Q_cond  = - q_satw * epsi * aLv / Rd / tem/tem * dTdz
   
   ! without hydrostatic assumption, i.e known pressure change 
   !Q_cond  = - q_satw * ( epsi * aLv / Rd / y%tem/y%tem * dTdz - dpdz / pe )
   
end subroutine kessler

end module ccfm_microphysics


module radiation

  !--------------------------------------------------------
  !  
  ! Set up and run all of the radiative calculations.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use fund_consts,        only : stefan, pi
  use dimensions,         only : nGas, nLay, nLev, nAng, nS, nTem
  use cross_section,      only : read_line_data

  use composition,        only : composition_init 
  use composition,        only : iGas_CO2, iGas_H2O, iGas_H2, iGas_N2 
  use composition,        only : iGas_CH4 

  use cia,                only : setup_cia

  use temperature_k_grid, only : setup_T_k_grid, calculate_T_k_grid, save_T_k_grid
  use longwave,           only : setup_longwave, calculate_longwave, save_longwave
  use shortwave,          only : setup_shortwave, calculate_shortwave, save_shortwave, Omega_stel, ISR

  implicit none
  private

  public :: setup_radiation, update_radiation, save_radiation

  !------- general variables ----------
  integer iLay, iLev, iAng                    ! do loop variables
  integer iS

  real(8) cosa(nAng)                         ! propagation angle cosine []
  real(8) ang_wt(nAng)                       ! propagation angle cosine weight []
  real(8) T_k_grid(nLay,nTem)                ! temperature array for cross-sections [K]

  ! for saving spectral results to file
  real(8) nu_lw(nS)                          ! longwave wavenumber [cm^-1]
  real(8) nu_sw(nS)                          ! shortwave wavenumber [cm^-1]
  real(8) OLRnu(nS)                          ! spectral outgoing longwave flux [W/m2/cm^-1] 
  real(8) OLRnu_temp(nS)                     ! spectral outgoing longwave irradiance [W/m2/cm^-1/sr] 
  real(8) OLR_temp                           ! outgoing longwave irradiance [W/m2/sr]
  real(8) OSRnu(nS)                          ! spectral outgoing shortwave flux [W/m2/cm^-1] 
  real(8) OSRnu_temp(nS)                     ! spectral outgoing shortwave irradiance [W/m2/cm^-1/sr] 
  real(8) OSR_temp                           ! outgoing shortwave irradiance [W/m2/sr]
  real(8) GSRnu(nS)                          ! spectral surface shortwave flux [W/m2/cm^-1] 
  real(8) GSRnu_temp(nS)                     ! spectral surface shortwave irradiance [W/m2/cm^-1/sr] 
  real(8) GSR_temp                           ! surface shortwave irradiance [W/m2/sr] 
  real(8) OSR                                ! outgoing shortwave flux [W/m2] 
  real(8) GSR                                ! surface shortwave flux [W/m2] 

  real(8) Ilev_up_lw(nLev)                   ! upward longwave irradiance in each layer [W/m2/sr]
  real(8) Ilev_dn_lw(nLev)                   ! downward longwave irradiance in each layer [W/m2/sr]
  real(8) Flev_up_lw(nLev)                   ! upward longwave flux in each layer [W/m2]
  real(8) Flev_dn_lw(nLev)                   ! downward longwave flux in each layer [W/m2]
  real(8) Ilev_up_sw(nLev)                   ! upward shortwave irradiance in each layer [W/m2/sr]
  !real(8) Ilev_dn_sw(nLev)                   ! downward shortwave irradiance in each layer [W/m2/sr]
  real(8) Flev_up_sw(nLev)                   ! upward shortwave flux in each layer [W/m2]
  real(8) Flev_dn_sw(nLev)                   ! downward shortwave flux in each layer [W/m2]

  real(8) dF_lw(nLay)                        ! longwave flux difference across each layer [W/m2]
  real(8) dF_sw(nLay)                        ! shortwave flux difference across each layer [W/m2]
  real(8) dTdt_lw(nLay)                      ! longwave heating rate in each layer [K/s]
  real(8) dTdt_sw(nLay)                      ! shortwave heating rate in each layer [K/s]

  !------- program options --------------  
  logical, parameter :: verbose = .false.

  save

contains
  ! ==================================================================================

  subroutine setup_radiation(calc_sigma,mu_i,f_i,play,Tlay,grav,ps)

    implicit none

    logical, intent(in) :: calc_sigma(nGas) ! calculate sigma for gas in question
    real(8), intent(in) :: mu_i(nGas)       ! molar mass of species [g/mol]
    real(8), intent(in) :: f_i(nLay,nGas)   ! species molar concentration [mol/mol]
    real(8), intent(in) :: play(nLay)       ! pressure layers [Pa]
    real(8), intent(in) :: Tlay(nLay)       ! temperature in layer  [K]
    real(8), intent(in) :: grav             ! gravity [m/s/s]
    real(8), intent(in) :: ps               ! surface pressure [Pa]

    ! mu_i is used here only to calculate Doppler broadening
    
    !-------- set up collision-induced absorption --------
    if(iGas_CO2.ne.-1 .and. iGas_H2.ne.-1)then
       call setup_cia('CO2H2_')
       !call setup_cia('N2_H2_')
       !print*,'Using N2-H2 CIA in place of CO2-H2!'
    endif
    if(iGas_CO2.ne.-1 .and. iGas_CH4.ne.-1)then
       call setup_cia('CO2CH4')
       !call setup_cia('N2_CH4')
       !print*,'Using N2-CH4 CIA in place of CO2-CH4!'
    endif
    if(iGas_CO2.ne.-1)then
       call setup_cia('CO2CO2')
    endif
    if(iGas_N2.ne.-1)then
       call setup_cia('N2_N2_')
    endif
    if(iGas_H2O.ne.-1)then
       call setup_cia('H2OH2O')
    endif

    !-------- read line data --------
    ! line data is only read if calc_sigma(iGas) = .true.
    call read_line_data(calc_sigma)

    call system('mkdir saved_sigma_data/')

    !-------- set up temperature grid on which cross-sections are defined --------
    call setup_T_k_grid(Tlay,T_k_grid)

    !-------- set up radiative transfer calculations --------
    call setup_shortwave(T_k_grid,play,ps,f_i,mu_i,grav,calc_sigma,nu_sw)
    call setup_longwave(T_k_grid,play,f_i,mu_i,calc_sigma,nu_lw)
    call define_cosa_quad(cosa,ang_wt)
    
  end subroutine setup_radiation

  ! ==================================================================================

  subroutine update_radiation(report_results,grav,cp_heat,kappa_gray,Ts,Tlay,Tlev,ps,dp,play,plev,f_i,mu_avg,dTdt,ASR,OLR)

    implicit none

    logical, intent(in)  :: report_results ! summarize the results ?
    real(8), intent(in)  :: grav           ! gravity [m/s/s]
    real(8), intent(in)  :: cp_heat        ! specific heat capacity for heating rate [J/K/kg]
    real(8), intent(in)  :: kappa_gray     ! gray gas mass absorption cross-section (if used) [m2/kg]
    real(8), intent(in)  :: Ts             ! surface temperature [K]
    real(8), intent(in)  :: Tlay(nLay)     ! temperature in layer  [K]
    real(8), intent(in)  :: Tlev(nLev)     ! temperature at level  [K]
    real(8), intent(in)  :: ps             ! surface pressure [Pa]
    real(8), intent(in)  :: dp(nLay)       ! pressure difference [Pa]
    real(8), intent(in)  :: play(nLay)     ! pressure layers [Pa]
    real(8), intent(in)  :: plev(nLev)     ! pressure at level [Pa]
    real(8), intent(in)  :: f_i(nLay,nGas) ! species molar concentration [mol/mol]
    real(8), intent(in)  :: mu_avg(nLay)   ! average molar mass in layer [g/mol]
    real(8), intent(out) :: dTdt(nLay)     ! net heating rate in each layer [K/s]
    real(8), intent(out) :: ASR            ! absorbed shortwave flux [W/m2]
    real(8), intent(out) :: OLR            ! outgoing longwave flux [W/m2]

    integer iT1(nLay)                      ! T-grid array points for linear interpolation [] 
    real(8) T_lin_weight(nS,nLay)          ! temperature weigting for linear interpolation [] 

    !-------- temperature grid interpolation calculation --------
    T_lin_weight = -1.0e8
    if(nTem>1) call calculate_T_k_grid(Tlay,iT1,T_lin_weight)

    !-------- visible radiative transfer calculation --------
    
    OSRnu(:)      = 0.0d0
    OSR           = 0.0d0
    GSR           = 0.0d0
    Flev_up_sw(:) = 0.0d0
    Flev_dn_sw(:) = 0.0d0

    if(verbose)then
       write(*,*) ' Calculating shortwave direct beam attenuation.'
    endif
    
    sw_angle_loop: do iAng = 1, nAng ! integration over propagation angle

       if(verbose)then
          write(*,'(a44,f8.2,a9)') ' Calculating shortwave with emission angle = ', &
               acos(cosa(iAng))*180.0d0/pi, ' degrees.'
       endif
       call calculate_shortwave(iT1,cosa(iAng),T_lin_weight,ps,dp,f_i,mu_avg,grav,kappa_gray, &
            OSRnu_temp,OSR_temp,GSRnu_temp,GSR_temp,Ilev_up_sw,Flev_dn_sw)

       ! c.f. e.g. Stamnes et al. (2000), DISORT tech. report, eqn. (9a)
       OSR        = OSR        + 2*pi*OSR_temp  *cosa(iAng)*ang_wt(iAng)
       OSRnu      = OSRnu      + 2*pi*OSRnu_temp*cosa(iAng)*ang_wt(iAng)
       Flev_up_sw = Flev_up_sw + 2*pi*Ilev_up_sw*cosa(iAng)*ang_wt(iAng)

       ! at the moment GSR does not depend on propagation angle as we are in the
       ! no atmospheric scattering regime. So no need to integrate it over cosa.
       !GSR        = GSR        + 2*pi*GSR_temp  *cosa(iAng)*ang_wt(iAng)
       !GSRnu      = GSRnu      + 2*pi*GSRnu_temp*cosa(iAng)*ang_wt(iAng)

    end do sw_angle_loop
    ! these two are fluxes in units of W/m2 and W/m2/cm^-1 already
    GSR        = GSR_temp
    GSRnu      = GSRnu_temp

    ! absorbed stellar radiation = ISR - OSR [W/m2]
    ASR = ISR - OSR

    !-------- IR radiative transfer calculation --------
    
    OLRnu(:)      = 0.0d0
    OLR           = 0.0d0
    Flev_up_lw(:) = 0.0d0
    Flev_dn_lw(:) = 0.0d0
    lw_angle_loop: do iAng = 1, nAng ! integration over propagation angle

       if(verbose)then
          write(*,'(a44,f8.2,a9)') ' Calculating longwave with emission angle = ', &
               acos(cosa(iAng))*180.0d0/pi, ' degrees.'
       endif

       call calculate_longwave(iT1,cosa(iAng),Tlay,Tlev,plev,Ts,ps,T_lin_weight,play,dp,f_i,mu_avg, &
            grav,kappa_gray,OLRnu_temp,OLR_temp,Ilev_up_lw,Ilev_dn_lw)

       OLR        = OLR        + 2*pi*OLR_temp  *cosa(iAng)*ang_wt(iAng)
       OLRnu      = OLRnu      + 2*pi*OLRnu_temp*cosa(iAng)*ang_wt(iAng)
       Flev_up_lw = Flev_up_lw + 2*pi*Ilev_up_lw*cosa(iAng)*ang_wt(iAng)
       Flev_dn_lw = Flev_dn_lw + 2*pi*Ilev_dn_lw*cosa(iAng)*ang_wt(iAng)

    end do lw_angle_loop

    !-------- heating rate calculation --------
    do ilay=1,nlay
       iLev = iLay + 1
       dF_lw(iLay) = (Flev_up_lw(iLev) - Flev_up_lw(iLev-1)) - (Flev_dn_lw(iLev) - Flev_dn_lw(iLev-1))
       dF_sw(iLay) = (Flev_up_sw(iLev) - Flev_up_sw(iLev-1)) - (Flev_dn_sw(iLev) - Flev_dn_sw(iLev-1))
    end do
    dTdt_lw(:) = -(grav/cp_heat)*dF_lw(:)/dp(:)
    dTdt_sw(:) = -(grav/cp_heat)*dF_sw(:)/dp(:)
    dTdt(:)    = dTdt_lw(:) + dTdt_sw(:)

    !-------- report results --------
    if(report_results)then
       write(*,*) '---------------------------------'
       write(*,*) 'Results:'
       write(*,'(a9,f14.7,a5)') '  ISR  = ',ISR,' W/m2'
       write(*,'(a9,f14.7,a5)') '  OLR  = ',OLR,' W/m2'
       write(*,'(a9,f14.7,a5)') '  ASR  = ',ASR,' W/m2'
       write(*,'(a9,f14.7,a2)') '  Tem  = ',(OLR/stefan)**0.25d0,' K'
       write(*,'(a9,f14.7,a2)') '  Tsrf = ',Ts,' K'
       write(*,'(a9,f14.7)')    '  OLRe = ',OLR/(stefan*Ts**4)
       write(*,'(a9,f14.7)')    '  Apla = ',OSR/ISR
       write(*,'(a9,f14.7)')    '  Aeqi = ',1.0d0 - OLR/ISR
       ! Aeqi is the albedo required for thermal equilibrium
    endif

  end subroutine update_radiation
  
  ! ==================================================================================

  subroutine save_radiation(verbose_output,play,f_i,OLR)

    implicit none

    logical, intent(in) :: verbose_output ! output files contain headers with variable info and units 
    real(8), intent(in) :: play(nLay)     ! pressure layers [Pa]
    real(8), intent(in) :: f_i(nLay,nGas) ! species molar concentration [mol/mol]
    real(8), intent(in) :: OLR            ! outgoing longwave flux [W/m2]

    !-------- save results to output file --------

    ! call save T_k_grid routine
    if(nTem>1) call save_T_k_grid(verbose_output,play)

    ! call save shortwave and longwave routines
    ! these save arrays internal to shorwave and longwave
    call save_shortwave(verbose_output,f_i)
    call save_longwave(verbose_output)
    
    ! save all other results
    
    open(unit=3,file='results/OLR.out')
    if(verbose_output)then
       write(3,*) 'Outgoing Longwave Radiation (OLR) [W/m2]'
    end if
    write(3,'(e12.5)') OLR
    close(3)

    open(unit=10,file='results/Flev_up.out')
    open(unit=11,file='results/Flev_dn.out')
    if(verbose_output)then
       write(10,*) 'Upwards broadband flux at each level [W/m2]'
       write(11,*) 'Downwards broadband flux at each level [W/m2]'
    end if
    do iLev = 1, nLev
       write(10,'(f14.4)') Flev_up_lw(iLev)
       write(11,'(f14.4)') Flev_dn_lw(iLev)
    end do
    close(10)
    close(11)

    open(unit=14,file='results/dTdt_lw.out')
    open(unit=15,file='results/dTdt_sw.out')
    if(verbose_output)then
       write(14,*) 'Longwave heating rate [K/s]'
       write(15,*) 'Shortwave heating rate [K/s]'
    end if
    do iLay = 1, nLay
       write(14,'(e12.4)') dTdt_lw(iLay)
       write(15,'(e12.4)') dTdt_sw(iLay)
    end do
    close(14)
    close(15)

    open(unit=2,file='results/nu_lw.out')
    open(unit=5,file='results/OLRnu.out')
    if(verbose_output)then
       write(2,*) 'Longwave spectral wavenumber [cm^-1]'
       write(5,*) 'Spectral outgoing longwave flux [W/m2/cm^-1]'
    end if
    do iS = 1, nS
       write(2,'(f12.4)') nu_lw(iS)
       write(5,'(e12.4)') OLRnu(iS)
    end do
    close(2)
    close(5)

    open(unit=2,file='results/nu_sw.out')
    open(unit=4,file='results/OSRnu.out')
    open(unit=5,file='results/GSRnu.out')
    if(verbose_output)then
       write(2,*) 'Shortwave spectral wavenumber [cm^-1]'
       write(4,*) 'Spectral outgoing shortwave flux [W/m2/cm^-1]'
       write(5,*) 'Spectral ground shortwave flux [W/m2/cm^-1]'
    end if
    do iS = 1, nS
       write(2,'(f12.4)') nu_sw(iS)
       if(OSRnu(iS)>1.0e-30)then
          write(4,'(e12.4)') OSRnu(iS)
       else
          write(4,'(e12.4)') 0.0d0
       end if
       if(GSRnu(iS)>1.0e-30)then
          write(5,'(e12.4)') GSRnu(iS)
       else
          write(5,'(e12.4)') 0.0d0
       end if
    end do
    close(2)
    close(5)


  end subroutine save_radiation

  ! ==================================================================================

  subroutine define_cosa_quad(cosa,ang_wt)

    ! calculate the gaussian nodes and weights for
    ! the mean emission angle cosine quadrature

    implicit none

    real(8), intent(inout) :: cosa(nAng)
    real(8), intent(inout) :: ang_wt(nAng)

    ! gaussian quadrature nodes and weights
    real(8) xi(nAng), ci(nAng)

    ! gauss interval (cosa = 0 to 1 here, for a stream in one hemisphere)
    real(8), parameter :: a = 0.0d0
    real(8), parameter :: b = 1.0d0

    ! move to dynamic allocation eventually

    if(nAng == 1)then
       xi(1) = 0.0d0
       ci(1) = 2.0d0
    elseif(nAng == 2)then
       xi(1) = +sqrt(1.0d0/3.0d0)
       xi(2) = -sqrt(1.0d0/3.0d0)
       ci(:) = 1.0d0
    elseif(nAng == 4)then
       xi(1) = -sqrt(3./7. + (2./7.)*sqrt(6./5.))
       xi(2) = -sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(3) = +sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(4) = +sqrt(3./7. + (2./7.)*sqrt(6./5.))
       ci(1) = (18.0d0-sqrt(30.0d0))/36.0d0
       ci(2) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(3) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(4) = (18.0d0-sqrt(30.0d0))/36.0d0
    else
       write(*,*) 'Error, emission angle weighting not defined for nAng = ', nAng
    endif

    cosa   = ((b-a)*xi + (b+a))/2
    ang_wt = ci*(b-a)/2
    
    return
  end subroutine define_cosa_quad

  ! ==================================================================================

end module radiation

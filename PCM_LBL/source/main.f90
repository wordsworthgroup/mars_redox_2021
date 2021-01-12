program PCM_LBL

  !--------------------------------------------------------
  !
  ! Planetary Climate Model (PCM-LBL)
  !
  ! A 1D line-by-line radiative convective model designed
  ! for diverse applications.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use dimensions,      only : nGas, nLay, nLev
  use composition,     only : composition_init
  use atmos_structure, only : atmos_structure_init, set_Tadia_fvar
  use radiation,       only : setup_radiation, update_radiation, save_radiation
 
  ! for RH diagnostic
  use composition,     only : igas_H2O
  use thermodynamics,  only : get_psat_H2O

  implicit none

  !------- namelist parameters --------
  logical calc_sigma(nGas)                   ! calculate sigma for gas in question
  integer nt                                 ! number of timesteps. Set to 1 for oneshot calculation.
  real(8) cp_heat                            ! specific heat capacity for heating rate [J/K/kg]
  real(8) rho_H_cp                           ! surface layer rho*height*cp [J/m2/K]
  real(8) delta_t                            ! timestep interval [s]
  real(8) kappa_gray                         ! gray gas mass absorption cross-section (if used) [m2/kg]

  !------- general variables ----------
  integer iLay, iLev                         ! do loop variables
  integer ierr, it, iGas

  real(8) f_i(nLay,nGas)                     ! species molar concentration [mol/mol]
  real(8) plev(nLev)                         ! pressure levels [Pa]
  real(8) Tlev(nLev)                         ! temperature at level  [K]
  real(8) Tlev_adia(nLev)                    ! adiabatic profile temperature at level  [K]
  real(8) play(nLay)                         ! pressure layers [Pa]
  real(8) dp(nLay)                           ! pressure difference [Pa]
  real(8) Tlay(nLay)                         ! temperature in layer  [K]
  real(8) Tlay_adia(nLay)                    ! adiabatic profile temperature in layer  [K]
  real(8) dTdt(nLay)                         ! net heating rate in each layer [K/s]
  real(8) ASR                                ! absorbed shortwave flux [W/m2]
  real(8) OLR                                ! outgoing longwave flux [W/m2]
  real(8) mu_i(nGas)                         ! molar mass of species [g/mol]
  real(8) mu_avg(nLay)                       ! average molar mass in layer [g/mol]
  real(8) ps                                 ! surface pressure [Pa]
  real(8) grav                               ! gravity [m/s/s]
  real(8) Ts                                 ! surface temperature [K]
  real(8), allocatable :: Ts_ar(:)           ! surface temperature array [K]
  real(8), allocatable :: TTOA_ar(:)         ! TOA temperature array [K]
  real(8), allocatable :: OLR_ar(:)          ! outgoing longwave radiation array [W/m2]
  real(8), allocatable :: ASR_ar(:)          ! absorbed shortwave radiation array [W/m2]
 
  ! for RH diagnostic
  real(8) psat_tmp                           ! saturation vapour pressure [K]

  !------- program options --------------  
  ! some of these should be in namelist file
  logical, parameter :: save_results   = .true.   ! save results to file?
  logical, parameter :: verbose_output = .true.  ! output files contain headers with variable info and units 
  logical add_update_Tsurf                        ! update surface temperature?
  logical relax_adiabat                           ! do dry or moist convection?

  namelist/gen_params_nml/nt,cp_heat,rho_H_cp,delta_t,calc_sigma,kappa_gray,add_update_Tsurf,relax_adiabat

  write(*,*) ""
  write(*,*) ""
  write(*,*) "-------------------------------"
  write(*,*) "PCM LBL version 0.7"
  write(*,*) "-------------------------------"
  write(*,*) ""
  write(*,*) "Robin Wordsworth (2016)"
  write(*,*) ""
  write(*,*) ""
  write(*,*) ""
  write(*,*) ""

  write(*,*) "loading parameters..."
  open(90,file='input.nml',status='old',form='formatted',iostat=ierr)
  if (ierr/=0) then
     print*, 'Cannot find required input.nml file, aborting.'
     call abort
  else
     write(*,*) 'ngas = ', ngas
     read(90,gen_params_nml)
     close(90)
  endif

  if(save_results)then   
     call system('mkdir results/')
     allocate(Ts_ar(nT))
     allocate(TTOA_ar(nT))
     allocate(OLR_ar(nT))
     allocate(ASR_ar(nT))
  end if
  
  write(*,*) "setting up atmospheric profile and spectral grid..."

  !-------- set up atmospheric composition --------
  call composition_init

  !-------- set up pressure, temperature and condensate profiles --------
  call atmos_structure_init(verbose_output,kappa_gray,Ts,Tlay,Tlev,play,plev,dp,ps,f_i,mu_i,mu_avg,grav)
      
  Tlay_adia = Tlay
  Tlev_adia = Tlev
   
  !-------- set up all radiation modules --------
  call setup_radiation(calc_sigma,mu_i,f_i,play,Tlay,grav,ps)
 
  time_loop: do it = 1,nt

     ! update fluxes and heating rates
     call update_radiation(.true.,grav,cp_heat,kappa_gray,Ts,Tlay,Tlev,ps,dp,play,plev,f_i,mu_avg,dTdt,ASR,OLR)

     ! update surface temperature and other variables
     if(save_results) Ts_ar(it)   = Ts
     if(save_results) TTOA_ar(it) = Tlay(nLay)
     if(save_results) OLR_ar(it)  = OLR
     if(save_results) ASR_ar(it)  = ASR

     ! exit loop now if in one-shot mode
     if (nt==1) exit

    ! update surface temperature
     if(add_update_Tsurf)then
        call update_T_surf(rho_H_cp,delta_t,ASR,OLR,Ts)
     end if

     ! recalculate adiabatic profile given new surface temperature
     if(relax_adiabat)then
        call set_Tadia_fvar(ps,Ts,Tlay,play,plev,dp,Tlay_adia,Tlev_adia,f_i)
     end if
     
     ! update atmospheric temperature profile
     call update_T_profile(it,dTdt,Tlay_adia,Tlev_adia,Ts,relax_adiabat,Tlay,Tlev)

     ! recalculate volatile molar concentration, if desired
     !if(vgas>0) call set_fvar(Tlay,play,Tlay_adia,f_i)

     if(relax_adiabat)then
        call set_Tadia_fvar(ps,Ts,Tlay,play,plev,dp,Tlay_adia,Tlev_adia,f_i)
     end if
     
     write(*,'(a15,i3,a1)') ' Finished step ',it,'.'

  end do time_loop

  ! call one more time to report final results
  call update_radiation(.true.,grav,cp_heat,kappa_gray,Ts,Tlay,Tlev,ps,dp,play,plev,f_i,mu_avg,dTdt,ASR,OLR)
  
  !-------- save results to output file --------

  if(save_results)then

     call save_radiation(verbose_output,play,f_i,OLR)

     open(unit=9, file='results/Tsurf.out')
     open(unit=10,file='results/TTOA.out')
     open(unit=11,file='results/OLR.out')
     open(unit=12,file='results/ASR.out')
     if(verbose_output)then
        write(9,*)  'Surface temperature vs. time [K]'
        write(10,*) 'TOA temperature vs. time [K]'
        write(11,*) 'Outgoing longwave flux vs. time [W/m2]'
        write(12,*) 'Absorbed shortwave flux vs. time [W/m2]'
     end if
     do it = 1, nt
        write(9,'(f12.4)' ) Ts_ar(it)
        write(10,'(f12.4)') TTOA_ar(it)
        write(11,'(e12.4)') OLR_ar(it)
        write(12,'(f12.4)') ASR_ar(it)
     end do
     close(9)
     close(10)
     close(11)
     close(12)

     open(unit=111,file='results/Tsurf_final.out')
     write(111,'(f12.4)') Ts_ar(nt)
     close(111)

     deallocate(Ts_ar)
     deallocate(TTOA_ar)
     deallocate(OLR_ar)
     deallocate(ASR_ar)
     
     open(unit=12,file='results/plev.out')
     open(unit=13,file='results/Tlev.out')
     if(verbose_output)then
        write(12,*) 'Pressure at each level [Pa]'
        write(13,*) 'Temperature at each level [K]'
     end if
     do iLev = 1, nLev
        write(12,'(e12.4)') plev(iLev)
        write(13,'(f12.4)') Tlev(iLev)
     end do
     close(12)
     close(13)
     
     open(unit=11,file='results/play.out')
     open(unit=12,file='results/Tlay.out')
     open(unit=13,file='results/Tlay_adia.out')
     open(unit=14,file='results/flay.out')
     if(verbose_output)then
        write(11,*) 'Pressure at each layer [Pa]'
        write(12,*) 'Temperature at each layer [K]'
        write(13,*) 'Adiabat temperature at each layer [K]'
        write(14,*) 'Species molar concentration at each layer [mol/mol]'
     end if
     do iLay = 1, nLay
        write(11,'(e12.4)') play(iLay)
        write(12,'(f12.4)') Tlay(iLay)
        write(13,'(f12.4)') Tlay_adia(iLay)
        do iGas = 1, nGas
           write(14,'(e12.4)') f_i(iLay,iGas)
        end do
     end do
     close(11)
     close(12)
     close(13)
     close(14)
    
     if(igas_H2O>0)then
        open(unit=15,file='results/psat.out')
        do iLay=1,nLay
           call get_psat_H2O(Tlay(iLay),psat_tmp)
           write(15,'(e12.4)') psat_tmp
        end do
        close(15)
     end if
 
  end if

end program PCM_LBL

  ! ==================================================================================

  subroutine update_T_profile(it,dTdt,Tlay_adia,Tlev_adia,Ts,relax_adiabat,Tlay,Tlev)

    use dimensions,      only : nLay, nLev
    
    implicit none

    integer, intent(in)    :: it                                 ! timestep []
    real(8), intent(in)    :: dTdt(nLay)                         ! net heating rate in each layer [K/s]
    real(8), intent(in)    :: Tlev_adia(nLev)                    ! adiabatic profile temperature at level  [K]
    real(8), intent(in)    :: Tlay_adia(nLay)                    ! adiabatic profile temperature in layer  [K]
    real(8), intent(in)    :: Ts                                 ! surface temperature [K]
    logical, intent(in)    :: relax_adiabat                      ! do dry or moist convection?
    real(8), intent(inout) :: Tlev(nLev)                         ! temperature at level  [K]
    real(8), intent(inout) :: Tlay(nLay)                         ! temperature in layer  [K]

    !------- general variables ----------
    integer iLay, iLev                                           ! do loop variables

    !--------  update atmospheric temperature profile --------

    Tlay    = Tlay + max(10.0d0/it,0.1d0)*dTdt/abs(dTdt + 1.0d-16)
    Tlev(1) = Ts
    do iLev=2,nLev-1
       iLay       = iLev
       Tlev(iLev) = (Tlay(iLay-1) + Tlay(iLay))/2.0d0
       ! could also weight this by path
    end do
    Tlev(nLev) = Tlay(nLay)
    !Tlev(nLev) = 10.0d0
    ! BCs for top Tlev don't seem to matter much in practice
    
    ! now relax to adiabatic profile where unstable
    if(relax_adiabat)then
       do iLay=1,nLay
          if(Tlay(iLay)<Tlay_adia(iLay)) Tlay(iLay) = Tlay_adia(iLay)
       end do
       do iLev=1,nLev
          if(Tlev(iLev)<Tlev_adia(iLev)) Tlev(iLev) = Tlev_adia(iLev) 
       end do
    end if

  end subroutine update_T_profile
  
  ! ==================================================================================
      
  subroutine update_T_surf(rho_H_cp,delta_t,ASR,OLR,Ts)

    implicit none

    real(8), intent(in)    :: rho_H_cp                           ! surface layer rho*height*cp [J/m2/K]
    real(8), intent(in)    :: delta_t                            ! timestep interval [s]
    real(8), intent(in)    :: ASR                                ! absorbed shortwave flux [W/m2]
    real(8), intent(in)    :: OLR                                ! outgoing longwave flux [W/m2]
    real(8), intent(inout) :: Ts                                 ! surface temperature [K]

      !------- subroutine options --------------  
    logical, parameter :: verbose      = .true.                  ! print more info for diagnostic

    !--------  update surface temperature --------

    ! ASR, OLR = W/m2
    ! rho*H*cp = kg/m3 m J/kg/K = J/m2/K
    ! W/m2 / J/m2/K = K/s
    
    Ts = Ts + (ASR - OLR)*delta_t/rho_H_cp

    if(verbose) write(*,'(a9,f7.2)') ' Tsurf = ',Ts
    
  end subroutine update_T_surf
  
  ! ==================================================================================
  

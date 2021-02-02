module shortwave

  !--------------------------------------------------------
  !  
  ! Solve the RTE in the shortwave portion of the spectrum.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use dimensions, only : nGas, nLay, nLev, nS, nTem

  implicit none
  private

  public :: setup_shortwave, calculate_shortwave, save_shortwave, Omega_stel, ISR

  !------- general variables ----------

  real(8) dnu_sw                             ! shortwave wavenumber interval [cm^-1]
  real(8) nu_sw(nS)                          ! shortwave wavenumber [cm^-1]
  real(8) sigma_sw(nS)                       ! absorption cross-section [cm2 / molecules of species] 
  real(8) sigma_UV(nS)                       ! absorption cross-section due to UV processes [cm2 / molecules of species] 
  real(8) sigma_sw_ar(nS,nLay,nGas,nTem)     ! absorption cross-section array by species and layer [cm2 / molecules of species] 
  real(8) dtau_sw(nLay,nS,nTem)              ! vertical path optical depth of each layer []
  real(8) log_dtau_sw(nLay,nS,2)             ! log10 of vertical path optical depth of each layer []
  real(8) dtau_sw_i(nLay,nS)                 ! optical depth of each layer at local temperature []
  real(8) dTran0(nLay,nS)                    ! transmission across each layer (direct beam) []
  real(8) dTran(nLay,nS)                     ! transmission across each layer (diffuse beam) []
  real(8) I_lev_up(nLev,nS)                  ! upward spectral irradiance at each level [W/m2/cm^-1/sr]
  real(8) I_lev_dn(nLev,nS)                  ! downward spectral irradiance at each level [W/m2/cm^-1/sr]
  real(8) I_stel_dn(nS)                      ! stellar spectral irradiance at TOA [W/m2/cm^-1/sr]
  real(8) F_stel_dn(nS)                      ! stellar spectral flux at TOA [W/m2/cm^-1]
  real(8) tau_sw_inf(nS)                     ! vertical path optical depth at top of atmosphere [] 
  real(8) Rs(nS)                             ! surface spectral reflectivity (albedo) [] 
  real(8) Rs_out(nS)                         ! surface spectral reflectivity (albedo) [] 
  real(8) Rpla(nS)                           ! planetary spectral reflectivity (albedo) [] 
  real(8) Rpla_out(nS)                       ! planetary spectral reflectivity (albedo) [] 
  real(8) ISR                                ! stellar flux (global average) [W/m2]

  real(8), parameter :: Omega_stel = 6.8d-5  ! stellar solid angle [sr]
  ! Sun at Earth by default, e.g.: http://www.geo.mtu.edu/~scarn/teaching/GE4250/radiation_lecture_slides.pdf
  
  integer iS, iLay, iLev, iTem               ! do loop variables
  integer iRev, iGas, ierr

  !------- subroutine options --------------  
  logical, parameter :: include_UV_abs  = .false. ! include UV absorption cross-sections?
  logical, parameter :: gray_debug      = .false. ! uniform opacity vs. wavenumber

  !------- namelist parameters --------
  real(8) nu_sw1                               ! shortwave starting wavenumber [cm^-1]
  real(8) nu_sw2                               ! shortwave finishing wavenumber [cm^-1]
  real(8) Asurf                                ! surface albedo []
  real(8) Fstel0                               ! stellar flux [W/m2]
  real(8) cosa0                                ! stellar zenith angle cosine []
  logical rayleigh_top                         ! apply Rayleigh scattering at TOA (increases albedo)?

  namelist/shortwave_nml/nu_sw1,nu_sw2,Asurf,Fstel0,cosa0,rayleigh_top

  save
  
contains

  ! ==================================================================================

  subroutine setup_shortwave(T_k_grid,play,ps,f_i,mu_i,grav,calc_sigma,nu_sw_out)

    use cross_section,  only : read_line_data, line_strength_width, get_line_abs_cross_section
    use composition,    only : gas_name
    !use ultraviolet,    only : setup_UV, calculate_UV

    implicit none

    logical, intent(in)  :: calc_sigma(nGas)    ! calculate sigma for gas in question
    real(8), intent(in)  :: T_k_grid(nLay,nTem) ! temperature array for cross-sections [K]
    real(8), intent(in)  :: play(nLay)          ! pressure layers [Pa]
    real(8), intent(in)  :: f_i(nLay,nGas)      ! species molar concentration [mol/mol]
    real(8), intent(in)  :: ps                  ! surface pressure [Pa]
    real(8), intent(in)  :: mu_i(nGas)          ! [g/mol]
    real(8), intent(in)  :: grav                ! gravity [m/s/s]
    real(8), intent(out) :: nu_sw_out(nS)       ! shortwave wavenumber [cm^-1]
    
    real(8) I_temp(nS)                          ! temporary spectral irradiance array [W/m2/cm^-1/sr]
    integer nTem_check                          ! number of cross-section temperatures at each layer
    real(8) play_check                          ! pressure layer [Pa]
    real(8) nu_sw_check                         ! shortwave wavenumber [cm^-1]
    
    logical, parameter   :: verbose = .true.

    ! read namelist 
    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
       print*, 'Cannot find required input.nml file, aborting.'
       call abort
    else
       read(10,shortwave_nml)
       close(10)
    endif

    ! divide by four to get flux averaged over the sphere
    ISR = Fstel0/4.0d0
        
    ! create fine spectral grid and initialize abs array
    nu_sw(1)  = nu_sw1
    dnu_sw = (nu_sw2-nu_sw1)/dble(nS-1)  ! Wavenumber interval [cm^-1]
    do iS = 2, nS
       nu_sw(iS) = nu_sw(iS-1) + dnu_sw
    end do
    nu_sw_out = nu_sw

    ! calculate surface reflectance (albedo)
    Rs(:) = Asurf

    ! calculate stellar Planck spectral irradiance vs. wavenumber
    call init_stellar_spectrum(I_temp)
    
    F_stel_dn(:) = I_temp(:)*ISR/sum(I_temp(:)*dnu_sw) ! stellar spectral flux [W/m2/cm^-1]
    I_stel_dn(:) = F_stel_dn(:)/(Omega_stel*cosa0)     ! stellar spectral irradiance [W/m2/cm^-1/sr]
    ! from classic flux defn. (e.g. Goody & Yung)

    ! paint the ground blue with Rayleigh scattering
    ! we can do this when vib-rot bands extend to the near-IR only
    if(.not.rayleigh_top) then 
       call calculate_rayleigh(ps,grav,Rs,Rs_out)
       Rs = Rs_out
    end if

    if(verbose)then
       open(unit=4,file='results/ISRnu.out')
       do iS = 1, nS
          write(4,*) F_stel_dn(iS)
       end do
       close(4)
    end if
    
    !-------- shortwave radiative transfer calculation --------

    ! calculate absorption cross section at each level
    sigma_sw(:)          = 0.0d0
    tau_sw_inf(:)        = 0.0d0
    sigma_sw_ar(:,:,:,:) = 0.0d0
    dtau_sw(:,:,:)       = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low
    
    gray_debug_if : if(gray_debug)then

       ! don't calculate anything in this case

    else

       big_species_loop : do iGas = 1, nGas

          ! read gas absorption cross section data from file if requested
          if(.not.calc_sigma(iGas))then

             ! check the number of temperatures for each level match
             open(unit=1,file='saved_sigma_data/nTem.dat')
             read(1,*) nTem_check
             if(nTem_check.ne.nTem)then
                write(*,*) 'Error: nTem in saved_sigma_data does not match!'
                stop
             end if
             close(1)

             ! check the pressure arrays match
             ! inefficient to do this for each gas but time cost is negligible; whatever
             ! we only do this in shortwave, not in longwave
             open(unit=2,file='saved_sigma_data/play.dat')
             do iLay = 1, nLay
                read(2,*) play_check
                if(abs(play_check-play(iLay))>1.0d-8)then
                   write(*,*) 'Error: play data in saved_sigma_data does not match!'
                   print*,play_check,' vs. ',play(iLay)
                   stop
                end if
             end do
             close(2)             
             
             ! check the spectral arrays match
             open(unit=2,file='saved_sigma_data/nu_sw.dat')
             do iS = 1, nS
                read(2,*) nu_sw_check
                if(abs(nu_sw_check - nu_sw(iS))>1.0d-8)then
                   write(*,*) 'Error: nu_sw data in saved_sigma_data does not match!'
                   print*,nu_sw_check,' vs. ',nu_sw(iS)
                   stop
                end if
             end do
             close(2)
             
             open(unit=3,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_sw.dat')
             do iS = 1, nS
                do iLay = 1, nLay
                   do iTem = 1,nTem
                      read(3,'(e12.4)') sigma_sw_ar(iS,iLay,iGas,iTem)
                   end do
                end do
             end do
             close(3)
             
          else

             do iLay = 1,nLay
                if(verbose)then
                   print*,'For shortwave at layer ',iLay,': '
                endif

                ! note no UV-cross-section temperature dependence included for now
                sigma_UV = 0.0d0
                if(include_UV_abs)then
                   if(gas_name(iGas)=='CO2' .or. gas_name(iGas)=='O3_') then
                      do iS = 1,nS
                         call calculate_UV(gas_name(iGas),nu_sw(iS),sigma_UV(iS),.false.)
                      end do
                   end if
                end if

                ! calculate cross-sections from scratch
                ! note these line widths will become inaccurate for a gas if
                ! a) it is a major consituent (f_i>0.1 or so) and
                ! b) its abundance varies with time
                do iTem = 1,nTem
                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),T_k_grid(iLay,iTem),iGas)
                   call get_line_abs_cross_section(iGas,nu_sw,sigma_sw,T_k_grid(iLay,iTem))

                   sigma_sw_ar(:,iLay,iGas,iTem) = sigma_sw + sigma_UV
                   ! cross-section per molecule of species
                end do

             end do

             ! save sigma data for future use
             ! save nu_sw, play and nTem to allow consistency checks
             open(unit=1,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_sw.dat')
             open(unit=2,file='saved_sigma_data/nu_sw.dat')
             do iS = 1, nS
                do iLay = 1, nLay
                   do iTem = 1,nTem
                      write(1,'(e12.4)') sigma_sw_ar(iS,iLay,iGas,iTem)
                   end do
                end do
                write(2,*) nu_sw(iS)
             end do
             close(1)
             close(2)
             open(unit=3,file='saved_sigma_data/nTem.dat')
             write(3,*) nTem
             close(3)
             open(unit=4,file='saved_sigma_data/play.dat')
             do iLay = 1, nLay
                write(4,*) play(iLay)
             end do
             close(4)

          end if

       end do big_species_loop

    end if gray_debug_if


  end subroutine setup_shortwave

  ! ==================================================================================

  subroutine calculate_shortwave(iT1,cosa_i,T_lin_weight,ps,dp,f_i,mu_avg,grav,kappa_gray, &
       OSRnu,OSR,GSRnu,GSR,I_lev_up_int,F_lev_dn_dir)

    use fund_consts, only : pi, mpr
    
    implicit none

    integer iT2(nLay)

    integer, intent(in)  :: iT1(nLay)             ! T-grid array points for linear interpolation [] 
    real(8), intent(in)  :: cosa_i                ! emission angle cosine []
    real(8), intent(in)  :: T_lin_weight(nS,nLay) ! temperature weigting for linear interpolation [] 
    real(8), intent(in)  :: ps                    ! surface pressure [Pa]
    real(8), intent(in)  :: dp(nLay)              ! pressure difference [Pa]
    real(8), intent(in)  :: f_i(nLay,nGas)        ! species molar concentration [mol/mol]
    real(8), intent(in)  :: mu_avg(nLay)          ! average molar mass in layer [g/mol]
    real(8), intent(in)  :: grav                  ! gravity [m/s/s]
    real(8), intent(in)  :: kappa_gray            ! gray gas mass absorption cross-section for debug [m2/kg]
    
    real(8), intent(out) :: GSRnu(nS)             ! ground shortwave spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: GSR                   ! ground shortwave irradiance [W/m2/sr] 
    real(8), intent(out) :: OSRnu(nS)             ! outgoing shortwave spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: OSR                   ! outgoing shortwave irradiance [W/m2/sr] 
    real(8), intent(out) :: I_lev_up_int(nLev)    ! upward irradiance in each layer [W/m2/sr]
    real(8), intent(out) :: F_lev_dn_dir(nLev)    ! downward direct flux in each layer [W/m2]

    real(8) I_lev_dn_int(nLev)                    ! downward (direct only) irradiance in each layer [W/m2/sr]

    ! note: we are currently neglecting H2O continuum effects in the shortwave
    ! fine for cold climates, will cause problems for RG calculations etc.

    ! as a result dtau_sw calculation is done a little differently from that for dtau_lw
    ! it is not at all optimized for speed yet.

    ! ----- beginning of part that requires nTem --------------

    dtau_sw(:,:,:) = 0.0d0

    if(gray_debug)then
       
       if(nTem>1)then
          write(*,*) 'nTem>1 and gray_debug incompatible.'
          stop
       end if
       do iLay = 1,nLay
          dtau_sw(iLay,:,1) = +kappa_gray*dp(ilay)/grav
          !dtau_sw(iLay,:,1) = 0.0d0
       end do

    else
    
       ! compute total layer optical depth (changes every timestep)
       do iGas = 1,nGas
          do iLay = 1, nLay
             ! sigma*f_i gives us cross-section per unit molecule of air.
             ! Hence to convert to dtau we need to use mu_avg(iLay), not mu_i.
             dtau_sw(iLay,:,:) = dtau_sw(iLay,:,:) + &
                  1.0d-4*sigma_sw_ar(:,iLay,iGas,:)*f_i(iLay,iGas)*dp(ilay)/(grav*mu_avg(iLay)*mpr)
          end do
       end do
       
    end if
    
    ! either just select 1st values in array, or
    ! interpolate over T grid to get cross-section
    ! at actual atmospheric temperature
    if(nTem==1)then   
       dtau_sw_i(:,:) = dtau_sw(:,:,1)
    else

       iT2(:) = iT1(:) + 1

       ! do log-space linear interpolation
       ! calculate logarithm of cross-section arrays
       ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
       ! note array flip [nS,nLay] ==> [nLay,nS]
       do iLay = 1,nLay
          log_dtau_sw(iLay,:,1:2) = log10(dtau_sw(iLay,:,iT1(iLay):iT2(iLay)) + 1.0d-50)
          dtau_sw_i(iLay,:)       = 10.0d0**(log_dtau_sw(iLay,:,1) + &
               (log_dtau_sw(iLay,:,2) - log_dtau_sw(iLay,:,1))*T_lin_weight(:,iLay))
       end do
       
    end if
    ! ----- end of part that requires nTem --------------

    ! calculate vertical path total optical depth at TOA
    ! somewhat inefficient to calculate it for every propagation angle
    tau_sw_inf = sum(dtau_sw_i(:,:),1) 

    ! calculate irradiance boundary conditions
    I_lev_dn(nLev,:) = I_stel_dn(:)
    
    ! calculate transmission
    ! note this is vertical path optical depth here
    dTran0 = exp(-dtau_sw_i/cosa0)  ! layer transmission for incoming stellar radiation (direct beam)

    ! Calculate downwards (direct beam) spectral irradiance 
    ! starting from the TOA down
    do iRev = 1,nLay
       iLev             = nLev - iRev ! start at nLev-1, finish at 1
       iLay             = iLev
       I_lev_dn(iLev,:) = I_lev_dn(iLev+1,:)*dTran0(iLay,:) 
    end do

    ! --- up to this point, nothing depends on emission angle ---
    ! --- make this into a separate subroutine eventually ---
    
    dTran  = exp(-dtau_sw_i/cosa_i) ! layer transmission for radiation scattered from surface

    ! make sure transmission does not cause floating point exception
    do iLay = 1,nLay
       do iS = 1,nS
          if(dTran(iLay,iS) < 1.0d-16)then           
             dTran(iLay,iS) = 1.0d-16
          end if
       end do
    end do

    ! Calculate upwards spectral irradiance 
    ! starting from the BOA up
    ! assuming a Lambertian surface

    ! F_dn = Omega_stel*I_dn*cosa0
    ! upwards is isotropic so by defn. I not a fn. of angle
    ! then F_up = 2*pi*I_up*sum(cosa_i*ang_wt) = pi*I_up
    ! conservation of energy: F_up = F_dn*Rs
    ! so  pi*I_up = Omega_stel*I_dn*cosa0*Rs

    I_lev_up(1,:) = I_lev_dn(1,:)*Rs(:)*(Omega_stel*cosa0/pi)
    
    do iLev = 2, nLev
       iLay             = iLev - 1
       I_lev_up(iLev,:) = I_lev_up(iLev-1,:)*dTran(iLay,:)
    end do
    
    ! calculate spectral fluxes and averaged irradiances for output
    ! trapezoidal quadrature for the irradiances
    GSRnu(:)        = I_lev_dn(1,:)*Omega_stel*cosa0 ! [W/m2/cm^-1]
    OSRnu(:)        = I_lev_up(nLev,:)               ! [W/m2/cm^-1/sr] (notation not ideal! this is an irradiance)

    ! calculate Rayleigh scattering at TOA
    if(rayleigh_top)then
       print*,'WARNING: rayleigh_top currently only works for nAng = 1.'
       Rpla(:)  = pi*OSRnu(:)/(F_stel_dn(:)+1.0d-16)
       call calculate_rayleigh(ps,grav,Rpla,Rpla_out)
       Rpla(:)  = Rpla_out(:)
       OSRnu(:) = Rpla(:)*I_lev_dn(nLev,:)*(Omega_stel*cosa0/pi)
    end if

    I_lev_up_int(:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_up(:,1) + I_lev_up(:,nS) + 2.0d0*sum(I_lev_up(:,2:nS-1),2) )
    I_lev_dn_int(:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_dn(:,1) + I_lev_dn(:,nS) + 2.0d0*sum(I_lev_dn(:,2:nS-1),2) )

    F_lev_dn_dir(:) = I_lev_dn_int(:)*Omega_stel*cosa0
    ! for now, this is correct as entire downwards irradiance is direct beam only
    ! downward direct flux in each layer [W/m2]

    ! note we could have done whole direct beam calculation with fluxes (Omega_stel value is irrelevant).
    
    GSR = 0.0d0
    OSR = 0.0d0
    do iS = 1, nS
       GSR = GSR + GSRnu(iS)
       OSR = OSR + OSRnu(iS)
    end do
    GSR = GSR*dnu_sw
    OSR = OSR*dnu_sw

    !-------- end shortwave radiative transfer calculation --------

  end subroutine calculate_shortwave

  ! ==================================================================================

  subroutine save_shortwave(verbose_output,f_i)
    
    implicit none

    logical, intent(in) :: verbose_output ! output files contain headers with variable info and units 
    real(8), intent(in) :: f_i(nLay,nGas) ! species molar concentration [mol/mol]

    real(8) sigma_total(nS,nLay,nTem)     ! total absorption cross-section [cm2 / molecules of air] 
    real(8) log_sigma_total(nS,nLay,nTem) ! log10 of total absorption cross-section [cm2 / molecules of air] 
    real(8) sigma_tot_save(nS,nLay)       ! total absorption cross-section [cm2 / molecules of air] 

    ! to be re-done at some point...
    
    !-------- save results to output file --------
    sigma_total(:,:,:) = 0.d0
    do iGas=1,nGas
       do iLay=1,nLay
          sigma_total(:,iLay,:) = sigma_total(:,iLay,:) + sigma_sw_ar(:,iLay,iGas,:)*f_i(iLay,iGas)
       end do
    end do
    log_sigma_total = log10(sigma_total + 1.0d-50)
    
    ! linear interpolation to scale for temperature (if required).
    if(nTem==1)then   
       sigma_tot_save = sigma_total(:,:,1)
    else
       do iLay = 1, nLay
          ! log-space linear interpolation
          sigma_tot_save(:,iLay) = 0.0d0
          !10.0d0**(log_sigma_total(iS,:,1) + & (log_sigma_total(iS,:,nTem)-log_sigma_total(iS,:,1))*(Tlay(:)-Tcold(:))/(2.0d0*dTlay))
          !sigma_tot_save(:,iLay) = 10.0d0**(sigma_total(iLay,:,1) + &
          !     sigma_total(iLay,:,2) - sigma_total(iLay,:,1))*T_lin_weight(:,iLay))
       end do
       ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
    end if

    open(unit=2,file='results/tau_sw_inf.out')
    open(unit=3,file='results/albedo.out')
    open(unit=4,file='results/sigma_total_sw.out')
    if(verbose_output)then
       write(2,*) 'Longwave optical depth [dimless]'
       write(3,*) 'Surface albedo [dimless]'
       write(4,*) 'Total absorption cross-section [cm2 / molecules of air]'
    end if
    do iS = 1, nS
       write(2,'(e12.4)') tau_sw_inf(iS)
       write(3,'(e12.4)') Rs(iS)
       do iLay = 1, nLay
          write(4,'(e12.4)') sigma_tot_save(iS,iLay)
       end do
    end do
    close(2)
    close(3)
    close(4)

    open(unit=5,file='results/Ilev_dn_sw.out')
    open(unit=6,file='results/Ilev_up_sw.out')
    if(verbose_output)then
       write(5,*) 'downward shortwave spectral irradiance at each level [W/m2/cm^-1/sr]'
       write(6,*) 'upward shortwave spectral irradiance at each level [W/m2/cm^-1/sr]'
    end if
    do iS = 1, nS
       do iLev = 1, nLev
          if(I_lev_dn(iLev,iS)<1.0d-30) I_lev_dn(iLev,iS) = 0.0d0
          if(I_lev_up(iLev,iS)<1.0d-30) I_lev_up(iLev,iS) = 0.0d0
          ! set very small values to zero so that saved quantities are not messed up
          write(5,'(e12.4)') I_lev_dn(iLev,iS)
          write(6,'(e12.4)') I_lev_up(iLev,iS)
       end do
    end do
    close(5)
    close(6)
    
  end subroutine save_shortwave

  ! ==================================================================================

  subroutine init_stellar_spectrum(I_temp)

    use dimensions, only : datadir
    use planck,     only : B_nu
    
    implicit none

    integer iS_stel, icount

    logical, parameter :: stellar_blackbody = .false.  ! use a blackbody curve for the stellar spectrum?
    integer, parameter :: nS_stel = 1000000             ! number of points in stellar spectrum

    real(8), intent(out) :: I_temp(nS)                ! spectral irradiance array [W/m2/cm^-1/sr]

    real(8) nu_stel(nS_stel)                          ! stellar spectrum wavenumber [cm^-1]
    real(8) F_stel(nS_stel)                           ! stellar spectrum spectral flux (unnormalized) [W/m2/cm^-1]

    I_temp(:) = 0.0d0
    
    if(stellar_blackbody)then
       do iS = 1, nS
          call B_nu(nu_sw(iS),3200.0d0,I_temp(iS))
          !call B_nu(nu_sw(iS),5800.0d0,I_temp(iS))
       end do
    else
       open(unit=3,file=trim(datadir)//'stellar_data/nu.dat')
       open(unit=4,file=trim(datadir)//'stellar_data/Fsol_3p8_Ga.dat')
       ! choice of spectra to be added later (just have filename in input.nml)
       do iS_stel = 1, nS_stel
          read(3,*) nu_stel(iS_stel)
          read(4,*) F_stel(iS_stel)
       end do
       close(3)       
       close(4)

       ! bin the stellar data by band
       ! we assume stellar resolution is higher, and both
       ! wavenumber grids are even
       ! techincally in this code we are using wavenumber bins
       ! cross-sections are evaluated at bin centers
       ! stellar spectrum is added to entire bins
       ! everything is normalized later so units are irrelevant here

       iS     = 1
       icount = 0
       do iS_stel = 1,nS_stel

          if(nu_stel(iS_stel)<nu_sw(iS)-dnu_sw/2.0d0)then
             ! if stellar spectrum point is below starting bin lower boundary, do nothing
       
          elseif(nu_stel(iS_stel)>nu_sw(iS)+dnu_sw/2.0d0)then
             ! if stellar spectrum point is above current bin upper boundary,
             ! move iS up one bin and add to that, or exit if we were at last bin.
             if(iS<nS)then
                I_temp(iS) = I_temp(iS)/dble(icount)      ! normalise amount in previous bin
                iS         = iS + 1                       ! move up one bin
                I_temp(iS) = I_temp(iS) + F_stel(iS_stel) ! add to total in new bin
                icount     = 1                            ! icount is reset
             else
                exit
             end if
          else
             ! if stellar spectrum point falls in current bin, just add to total.
             I_temp(iS) = I_temp(iS) + F_stel(iS_stel)
             icount     = icount + 1 ! one point added to current bin
          end if

       end do
       
    end if

    return

  end subroutine init_stellar_spectrum
  
  ! ==================================================================================

  subroutine calculate_rayleigh(ps,grav,Rbelow,Rabove)

    implicit none

    ! no temperature dependence so we only need to do this once
    ! unless we have really variable major constituents in the atmosphere

    ! general variables
    !real(8), parameter :: gamma = 0.5d0*sqrt(3.0d0) ! quadrature approximation
    real(8), parameter   :: gamma = 3.0d0/4.0d0       ! Eddington approximation

    real(8), intent(in)  :: ps         ! surface pressure [Pa]
    real(8), intent(in)  :: grav       ! gravity [m/s/s]
    real(8), intent(in)  :: Rbelow(nS) ! reflectivity/albedo below Rayleigh layer []
    real(8), intent(out) :: Rabove(nS) ! reflectivity/albedo above Rayleigh layer []

    real(8) beta(nS)
    real(8) f(nS)
    real(8) R_dn(nS)    ! reflectivity for downward beam
    real(8) R_up(nS)    ! reflectivity for upward beam
    !real(8) Rpla(nS)    ! planetary (net) reflectivity
    real(8) tau_Ray(nS) ! total vertical path Rayleigh optical depth []
    real(8) lam_sw(nS)  ! wavelength [um]
    
    ! asymmetry parameter is zero for Rayleigh scattering

    lam_sw  = 1.0d4/nu_sw 
    ! kg/m2 x m2/kg = []
    !tau_Ray = (ps/grav)*1.2d-6*lam_sw**(-4) ! PPC Table 5.2
    tau_Ray = 1.4098d-6*lam_sw**(-4)*(1.0d0 + 0.013d0*lam_sw**(-2))*(ps/grav) ! Hansen & Travis eqn. (2.32)
    ! (1.527/(93*101325))*8.7 = 1.4098e-6
    ! At 1.5 bar CO2 on Mars albedo is 0.02 greater for 2nd formula

    ! R_dn and R_up were previously incorrectly reversed here.
    ! Note however when gamma = 3/4 and cosa0 = 0.666, R_dn = R_up.

    beta    = 1.0d0 - exp(-tau_Ray/cosa0)
    f       = gamma*tau_Ray
    R_dn    = ((0.5d0 - gamma*cosa0)*beta + f)/(1.0d0 + f)                             ! downwards beam albedo
    R_up    = f/(1.0d0 + f)                                                            ! upwards beam albedo
    Rabove  = 1.0d0 - (1 - R_up)*(1 - Rbelow)/((1.0d0 - Rbelow)*R_dn + (1.0d0 - R_dn)) ! planetary albedo

    ! should work for clear-sky atmospheres because Rayleigh scattering and
    ! molecular absorption occur in separate regions.
    ! literally we are painting the ground (or TOA) blue.
   
  end subroutine calculate_rayleigh

  ! ==================================================================================

  real(8) function m_real(lam)

    ! just make dominant gas be the main rayleigh scatterer
    ! not used for now
    
    implicit none
    
    real(8), intent(in) :: lam ! wavelength [um]

    !if(iGas_CO2==1)
    ! pure CO2 real refractive index
    ! http://refractiveindex.info/?shelf=main&book=CO2&page=Bideau-Mehu
    ! first two terms from Bidea-Mehu et al., Opt. Commun., (1973)
    m_real = 1.0d0 + 6.99100d-2/(166.175d0 - 1.0d0/lam**2) + 1.44720d-3/(79.609 - 1.0d0/lam**2)
    !else
    !endif
    
  end function m_real
   
  ! ==================================================================================
  
end module shortwave

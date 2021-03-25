module longwave
  
  !--------------------------------------------------------
  !  
  ! Solve the RTE in the longwave portion of the spectrum.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use dimensions,  only : nGas, nLay, nLev, nAng, nS, nTem
  use composition, only : vgas

  implicit none
  private

  public :: setup_longwave, calculate_longwave, save_longwave

  !------- general variables ----------

  real(8) dnu_lw                             ! longwave wavenumber interval [cm^-1]  
  real(8) nu_lw(nS)                          ! longwave wavenumber [cm^-1]
  real(8) sigma_lw(nS)                       ! absorption cross-section [cm2 / molecules of species] 
  real(8) sigma_CIA(nS)                      ! absorption cross-section due to CIA [cm2 / molecules of air] 

  real(8) sigma_lw_ar(nS,nLay,nGas,nTem)     ! absorption cross-section array by species and layer [cm2 / molecules of species] 
  real(8) sigma_total_dry(nS,nLay,nTem)      ! dry total absorption cross-section (no variable gas) [cm2 / molecules of air] 
  real(8) sigma_t_i(nS,nLay)                 ! interpolated total absorption cross-section [cm2 / molecules of air] 

  real(8) dtau_lw(nLay,nS)                   ! vertical path optical depth of each layer []

  real(8) dtau_lw_a(nLay,nS)                 ! emission angle weigted optical depth of each layer at local temperature []
  real(8) dTran(nLay,nS)                     ! transmission across each layer []
  real(8) Bnu_lay(nLay,nS)                   ! Planck spectral irradiance in each layer [W/m2/cm^-1/sr]
  real(8) Bnu_lev(nLev,nS)                   ! Planck spectral irradiance in each level [W/m2/cm^-1/sr]
  real(8) Bnu_s(nS)                          ! Planck spectral irradiance at the surface [W/m2/cm^-1/sr]
  real(8) tau_lw_inf(nS)                     ! vertical path optical depth at top of atmosphere [] 
  real(8) I_lev_up(nLev,nS)                  ! upward spectral irradiance in each layer [W/m2/cm^-1/sr]
  real(8) I_lev_dn(nLev,nS)                  ! downward spectral irradiance in each layer [W/m2/cm^-1/sr]
  
  ! we also need to define Bnu_s because in some test cases such as
  ! gray radiative equilibrium, there is a thermal discontinuity
  ! at the surface

  integer iS, iLay, iLev, iTem               ! do loop variables
  integer iRev, iGas, ierr

  !------- program options --------------  
  logical, parameter :: gray_debug   = .false. ! uniform opacity vs. wavenumber
  logical, parameter :: save_results = .true.  ! save results to file?

  !------- namelist parameters --------
  real(8) nu_lw1                               ! longwave starting wavenumber [cm^-1]
  real(8) nu_lw2                               ! longwave finishing wavenumber [cm^-1]

  namelist/longwave_nml/nu_lw1,nu_lw2

  save
  
contains

  ! ==================================================================================

  subroutine setup_longwave(T_k_grid,play,f_i,mu_i,calc_sigma,nu_lw_out)

    use planck,         only : B_nu
    use cross_section,  only : line_strength_width
    use cross_section,  only : get_line_abs_cross_section
    use cross_section,  only : get_CIA_cross_section
    use composition,    only : gas_name
    
    implicit none

    logical, intent(in)  :: calc_sigma(nGas)    ! calculate sigma for gas in question
    real(8), intent(in)  :: T_k_grid(nLay,nTem) ! temperature array for cross-sections [K]

    real(8), intent(in)  :: play(nLay)          ! pressure layers [Pa]
    real(8), intent(in)  :: f_i(nLay,nGas)      ! species molar concentration [mol/mol]
    real(8), intent(in)  :: mu_i(nGas)          ! molar mass of species [g/mol]
    real(8), intent(out) :: nu_lw_out(nS)       ! longwave wavenumber [cm^-1]

    real(8) nu_lw_check                         ! longwave wavenumber [cm^-1]

    !------- subroutine options --------------  
    logical, parameter   :: verbose = .true.

    ! read namelist 
    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
       print*, 'Cannot find required input.nml file, aborting.'
       call abort
    else
       read(10,longwave_nml)
       close(10)
    endif
    
    ! create fine spectral grid and initialize abs array
    nu_lw(1)  = nu_lw1
    dnu_lw = (nu_lw2-nu_lw1)/dble(nS-1)  ! Wavenumber interval [cm^-1]
    do iS = 2, nS
       nu_lw(iS) = nu_lw(iS-1) + dnu_lw
    end do
    nu_lw_out = nu_lw
    
    ! calculate absorption cross section at each layer

    sigma_lw(:)            = 0.0d0
    tau_lw_inf(:)          = 0.0d0
    sigma_lw_ar(:,:,:,:)   = 0.0d0
    sigma_total_dry(:,:,:) = 0.0d0
    sigma_t_i(:,:)         = 0.0d0
    dtau_lw(:,:)           = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low
    dtau_lw(:,:)           = +1.0d-10 ! 

    gray_debug_if : if(gray_debug)then

       ! don't calculate anything in this case
       write(*,*) 'Assuming gray gas in longwave.'

    else

       big_species_loop : do iGas = 1, nGas
          
          ! read gas absorption cross section data from file if requested
          if(.not.calc_sigma(iGas))then

             ! no need to check that the number of temperatures for each layer match
             ! this was done in setup_shortwave
             
             ! check the spectral arrays match
             open(unit=2,file='saved_sigma_data/nu_lw.dat')
             do iS = 1, nS
                read(2,*) nu_lw_check
                if(abs(nu_lw_check-nu_lw(iS))>1.0d-8)then
                   write(*,*) 'Error: nu_lw data in saved_sigma_data does not match!'
                   print*,nu_lw_check,' vs. ',nu_lw(iS)
                   stop
                end if
             end do
             close(2)

             ! no need to check that the pressure layers match
             ! this was done in setup_shortwave
             
             open(unit=3,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
             do iS = 1, nS
                do iLay = 1, nLay
                   do iTem = 1,nTem
                      read(3,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
                   end do
                end do
             end do
             close(3)
             
          else

             ! calculate cross-sections from scratch
             do iLay = 1,nLay
                if(verbose)then
                   print*,'For longwave at layer ',iLay,': '
                endif

                do iTem = 1,nTem
                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),T_k_grid(iLay,iTem),iGas)
                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,T_k_grid(iLay,iTem))
                   sigma_lw_ar(:,iLay,iGas,iTem) = sigma_lw
                   ! cross-section per molecule of species
                end do
                
             end do

             ! save sigma data for future use
             ! save nu_lw to allow consistency checks
             open(unit=1,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
             open(unit=2,file='saved_sigma_data/nu_lw.dat')
             do iS = 1, nS
                do iLay = 1, nLay
                   do iTem = 1,nTem
                      write(1,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
                   end do
                end do
                write(2,*) nu_lw(iS)
             end do
             close(1)
             close(2)

          end if

          ! multiply the species cross-section by molar concentration to get the
          ! contribution to the total cross-section, in cm2/molecule of air.
          ! at this point we also add the CIA and compute the total 'dry' cross-section.
          ! do not add variable gas cross-section to the total yet.
          not_vgas : if(.not. iGas == vgas)then
             if (any(f_i(:,iGas)*play(:)>1.0d1)) then ! include CIA if partial pressure > 10 Pa
                do iLay = 1, nLay

                   sigma_CIA(:) = 0.0d0

                   do iTem = 1,nTem
                      call get_CIA_cross_section(iGas,nu_lw,sigma_CIA,T_k_grid(iLay,iTem),play(iLay),f_i(iLay,:))
                      ! note that f_i is sent as a nGas array to get_CIA_cross_section even when we are inside
                      ! an iGas loop. this is because we need to know the abundance of the partner gas.
                      sigma_total_dry(:,iLay,iTem) = sigma_total_dry(:,iLay,iTem) + &
                           sigma_lw_ar(:,iLay,iGas,iTem)*f_i(iLay,iGas) + sigma_CIA
                   end do
                
                end do
             else
                do iLay = 1, nLay
                   do iTem = 1,nTem
                      sigma_total_dry(:,iLay,iTem) = sigma_total_dry(:,iLay,iTem) + sigma_lw_ar(:,iLay,iGas,iTem)*f_i(iLay,iGas)
                   end do
                end do
             end if
          end if not_vgas
          
       end do big_species_loop
    end if gray_debug_if

  end subroutine setup_longwave

  ! ==================================================================================

  subroutine calculate_longwave(iT1,cosa_i,Tlay,Tlev,plev,Ts,ps,T_lin_weight,play,dp,f_i,mu_avg, &
       grav,kappa_gray,OLRnu,OLR,I_lev_up_int,I_lev_dn_int)

    use fund_consts,    only : mpr
    use planck,         only : B_nu
    use cross_section,  only : get_CIA_cross_section

    implicit none

    integer iT2(nLay)

    integer, intent(in)  :: iT1(nLay)             ! T-grid array points for linear interpolation [] 
    real(8), intent(in)  :: cosa_i                ! emission angle cosine []
    real(8), intent(in)  :: Tlay(nLay)            ! temperature in layer [K]
    real(8), intent(in)  :: Tlev(nLev)            ! temperature at level [K]
    real(8), intent(in)  :: plev(nLev)            ! pressure at level [Pa]
    real(8), intent(in)  :: Ts                    ! surface temperature [K]
    real(8), intent(in)  :: ps                    ! surface pressure [Pa]
    real(8), intent(in)  :: T_lin_weight(nS,nLay) ! temperature weigting for linear interpolation [] 
    real(8), intent(in)  :: play(nLay)            ! pressure layers [Pa]
    real(8), intent(in)  :: dp(nLay)              ! pressure difference [Pa]
    real(8), intent(in)  :: f_i(nLay,nGas)        ! species molar concentration [mol/mol]
    real(8), intent(in)  :: mu_avg(nLay)          ! average molar mass in layer [g/mol]
    real(8), intent(in)  :: grav                  ! gravity [m/s/s]
    real(8), intent(in)  :: kappa_gray            ! gray gas mass absorption cross-section for debug [m2/kg]

    real(8), intent(out) :: OLRnu(nS)             ! outgoing longwave spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: OLR                   ! total outgoing longwave irradiance [W/m2/sr] 
    real(8), intent(out) :: I_lev_up_int(nLev)    ! upward irradiance in each layer [W/m2/sr]
    real(8), intent(out) :: I_lev_dn_int(nLev)    ! downward irradiance in each layer [W/m2/sr]

    real(8) log_sig_d(nS,nLay,2)                  ! logarithm of dry absorption cross-section at grid temperatures surrounding layer T [] 
    real(8) log_sig_v(nS,nLay,2)                  ! logarithm of variable gas absorption cross-section at grid temperatures surrounding layer T [] 
    real(8) sigma_d_i(nS,nLay)                    ! interpolated dry gas absorption cross-section [cm2 / molecules of air] 
    real(8) sigma_v_i(nS,nLay)                    ! interpolated variable gas absorption cross-section [cm2 / molecules of variable gas] 
  
    real(8) Beff(nS)                              ! Effective Planck function across each level [W/m2/cm^-1/sr]
    real(8) Cterm(nLay,nS)                        ! dimensionless fn. of tau, varies between 0.5 and zero

    sigma_d_i(:,:) = 0.0d0
    sigma_v_i(:,:) = 0.0d0
    dtau_lw(:,:)   = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low
    dtau_lw(:,:)   = +1.0d-10 ! this value gives better results for very high atm 
    
    ! ----- beginning of part that requires nTem --------------

    gray_debug_loop : if(gray_debug)then

       ! nTem = 1 for this case verified in setup_shortwave
       do iLay = 1,nLay
          do iS = 1,nS
             ! adjust if loop to give semigray behaviour
             if(nu_lw(iS)>600.0d0 .and. nu_lw(iS)<650.0d0)then
                dtau_lw(iLay,iS) = +kappa_gray*dp(ilay)/grav
                !kappa_semi       = kappa_gray*10.0d0**(12.0d0*(nu_lw(iS)-500.0d0)/100.0d0)
                !dtau_lw(iLay,iS) = kappa_semi*dp(ilay)/grav
                !dtau_lw(iLay,iS) = +kappa_gray*((nu_lw(iS)-500.0d0)/100.0d0)*dp(ilay)/grav
             end if
          end do
          dtau_lw(iLay,:) = +kappa_gray*dp(ilay)/grav
       end do
       !dtau_lw(:,:) = +1.0d-6

    else

       ! either just select 1st values in array, or
       ! interpolate over T grid to get cross-section
       ! at actual atmospheric temperature
       if(nTem==1)then   
          sigma_d_i = sigma_total_dry(:,:,1)
          if(vgas.ne.0) sigma_v_i = sigma_lw_ar(:,:,vgas,1)
       else
          iT2(:) = iT1(:) + 1

          ! do log-space linear interpolation
          ! calculate logarithm of cross-section arrays
          ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
          do iLay = 1,nLay
             log_sig_d(:,iLay,:) = log10(sigma_total_dry(:,iLay,iT1(iLay):iT2(iLay)) + 1.0d-50)
             if(vgas.ne.0) log_sig_v(:,iLay,:) = log10(sigma_lw_ar(:,iLay,vgas,iT1(iLay):iT2(iLay)) + 1.0d-50)
          end do
          sigma_d_i(:,:) = 10.0d0**(log_sig_d(:,:,1) + (log_sig_d(:,:,2) - log_sig_d(:,:,1))*T_lin_weight(:,:))
          if(vgas.ne.0) sigma_v_i(:,:) = 10.0d0**(log_sig_v(:,:,1) + (log_sig_v(:,:,2) - log_sig_v(:,:,1))*T_lin_weight(:,:))
       end if
    
       ! get sigma_CIA for variable gas (continuum absorption if H2O)
       ! and add variable gas cross-section (single-molecule + CIA) to the total
       ! for speed, we assume that _only_ the mixing ratio of the volatile gas
       ! 'iGas_var' may change at each timestep.
       if(vgas==0)then
          sigma_t_i = sigma_d_i
       else
          do iLay = 1, nLay

             sigma_CIA(:) = 0.0d0

             ! calculate variable gas continuum at exact instantaneous atmospheric temperature, so no T interpolation
             call get_CIA_cross_section(vgas,nu_lw,sigma_CIA,Tlay(iLay),play(iLay),f_i(iLay,:))
             sigma_t_i(:,iLay) = sigma_d_i(:,iLay) + sigma_v_i(:,iLay)*f_i(iLay,vgas) + sigma_CIA

          end do
       endif

       ! compute total layer optical depth
       ! note array flip [nS,nLay] ==> [nLay,nS]
       ! also note g/mol == molecular mass in Daltons in denominator
       do iLay = 1, nLay
          dtau_lw(iLay,:) = dtau_lw(iLay,:) + 1.0d-4*sigma_t_i(:,iLay)*dp(ilay)/(grav*mu_avg(iLay)*mpr)
       end do

    end if gray_debug_loop
    
    ! calculate Planck function vs. height and wavenumber
    do iS = 1, nS
       do iLev = 1, nLev
          call B_nu(nu_lw(iS),Tlev(iLev),Bnu_lev(iLev,iS))
       end do
       do iLay = 1, nLay
          call B_nu(nu_lw(iS),Tlay(iLay),Bnu_lay(iLay,iS))
       end do
       call B_nu(nu_lw(iS),Ts,Bnu_s(iS))
    end do

    ! calculate irradiance boundary conditions
    I_lev_up(1,:)    = Bnu_s(:)
    I_lev_dn(nLev,:) = 0.0d0

    ! scale dtau_lw by emission angle cosine
    dtau_lw_a(:,:) = dtau_lw(:,:)/cosa_i

    ! ----- end of part that requires nTem --------------

    ! vertical path optical depth at TOA
    ! somewhat inefficient to calculate it for every propagation angle
    tau_lw_inf(:) = sum(dtau_lw(:,:),1)
    
    ! calculate dTran
    ! the transmission through each layer,
    ! individually, for a given propagation angle.
    dTran = exp(-dtau_lw_a)

    ! make sure transmission does not cause floating point exception
    do iLay = 1,nLay
       do iS = 1,nS
          if(dTran(iLay,iS) < 1.0d-16)then           
             dTran(iLay,iS) = 1.0d-16
          end if
       end do
    end do

    ! calculate Cterm (not dependent on flux direction)
    Cterm(:,:) = 1.0d0/dtau_lw_a(:,:) - dTran(:,:)/(1.0d0 - dTran(:,:))

    ! Calculate upwards spectral irradiance using Clough et al. (1992) method
    ! starting from the BOA up
    do iLev = 2, nLev
       iLay             = iLev - 1
       Beff             = Bnu_lev(iLev,:) + 2.0d0*(Bnu_lay(iLay,:) - Bnu_lev(iLev,:))*Cterm(iLay,:)
       I_lev_up(iLev,:) = I_lev_up(iLev-1,:)*dTran(iLay,:) + Beff*(1.0d0 - dTran(iLay,:))
    end do
    ! Calculate downwards spectral irradiance using Clough et al. (1992) method
    ! starting from the TOA down
    do iRev = 1,nLay
       iLev             = nLev - iRev ! start at nLev-1, finish at 1
       iLay             = iLev
       Beff             = Bnu_lev(iLev,:) + 2.0d0*(Bnu_lay(iLay,:) - Bnu_lev(iLev,:))*Cterm(iLay,:)
       I_lev_dn(iLev,:) = I_lev_dn(iLev+1,:)*dTran(iLay,:) + Beff*(1.0d0 - dTran(iLay,:))
    end do

    ! calculate OLR and irradiances for output
    ! trapezoidal quadrature for the irradiances
    OLRnu(:)        = I_lev_up(nLev,:) 
    I_lev_up_int(:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * (I_lev_up(:,1) + I_lev_up(:,nS) + 2.0d0*sum(I_lev_up(:,2:nS-1),2) )
    I_lev_dn_int(:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * (I_lev_dn(:,1) + I_lev_dn(:,nS) + 2.0d0*sum(I_lev_dn(:,2:nS-1),2) )

    OLR = 0.0d0
    do iS = 1, nS
       OLR = OLR + OLRnu(iS)
    end do
    OLR = OLR*dnu_lw

    !-------- end longwave radiative transfer calculation --------

  end subroutine calculate_longwave

  ! ==================================================================================

  subroutine save_longwave(verbose_output)
    
    implicit none

    logical, intent(in) :: verbose_output ! output files contain headers with variable info and units 

    !-------- save results to output file --------

    open(unit=2,file='results/tau_lw_inf.out')
    open(unit=3,file='results/sigma_total_lw.out')
    if(verbose_output)then
       write(2,*) 'Longwave optical depth [dimless]'
       write(3,*) 'Total absorption cross-section [cm2 / molecules of air]'
    end if
    do iS = 1, nS
       write(2,'(e12.4)') tau_lw_inf(iS)
       do iLay = 1, nLay
          write(3,'(e12.4)') sigma_t_i(iS,iLay)
       end do
    end do
    close(2)
    close(3)

    ! for now, we do not save the irradiances for every single propagation angle
    ! so right now the final values saved corresponed to the final cosai value
    ! in the integration.
    
    open(unit=11,file='results/Ilev_up.out')
    if(verbose_output)then
       write(11,*) 'Upwards spectral irradiance at each level [W/m2]'
    end if
    do iLev = 1, nLev
       do iS = 1, nS
          write(11,'(f14.4)') I_lev_up(iLev,iS)
       end do
    end do
    close(11)

    open(unit=5,file='results/Ilev_dn_lw.out')
    open(unit=6,file='results/Ilev_up_lw.out')
    if(verbose_output)then
       write(5,*) 'downward longwave spectral irradiance at each level [W/m2/cm^-1/sr]'
       write(6,*) 'upward longwave spectral irradiance at each level [W/m2/cm^-1/sr]'
    end if
    do iS = 1, nS
       do iLev = 1, nLev
          write(5,'(e12.4)') I_lev_dn(iLev,iS)
          write(6,'(e12.4)') I_lev_up(iLev,iS)
       end do
    end do
    close(5)
    close(6)
    
  end subroutine save_longwave

  ! ==================================================================================

end module longwave

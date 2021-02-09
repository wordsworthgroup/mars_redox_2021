module cross_section

  !--------------------------------------------------------
  !  
  ! Calculate molecular absorption cross-sections.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use dimensions, only : nGas, nlinesMAX
  use composition, only : gas_name,iGas_CO2,iGas_H2O,iGas_CH4,iGas_O3,iGas_SO2,iGas_H2S,iGas_H2,iGas_NH3,iGas_N2

  implicit none
  private

  public :: read_line_data, line_strength_width, get_line_abs_cross_section, get_CIA_cross_section

  integer nlines(nGas)             ! actual number of lines
  integer, parameter :: nlines_beguier = 7786 ! number of lines in Beguier+ JQRST (2015) CH4 data

  !------- inputs --------
  logical HITEMP                   ! use HITEMP line lists?
  real(8) deltanu_trunc(nGas)      ! truncation wavenumber [cm^-1]
  real(8) Sref_cutoff              ! cutoff line intensity [cm^-1 / cm^2 molecule^-1]

  ! constants once loaded
  integer mol(nGas,nlinesMAX)      ! molecule number
  integer iso(nGas,nlinesMAX)      ! isotope number
  real(8) nu0(nGas,nlinesMAX)      ! wavenumber [cm^-1]
  real(8) Sref(nGas,nlinesMAX)     ! reference line intensity [cm^-1 / cm^2 molecule^-1]
  real(8) einstein(nGas,nlinesMAX) ! Einstein A coefficient [s^-1]
  real(8) gam_air(nGas,nlinesMAX)  ! Air-broadened half width [cm^-1 atm^-1]
  real(8) gam_self(nGas,nlinesMAX) ! Self-broadened half width [cm^-1 atm^-1]
  real(8) Egnd(nGas,nlinesMAX)     ! Ground state energy [cm^-1]
  real(8) ncoeff(nGas,nlinesMAX)   ! Exponent for pressure broadening T-dependence
  real(8) delta(nGas,nlinesMAX)    ! Air-broadened pressure shift [cm^-1/atm]

  ! functions of temperature
  real(8) gam_lore(nGas,nlinesMAX) ! Lorentz half width [cm^-1 atm^-1]
  real(8) gam_dopp(nGas,nlinesMAX) ! Doppler half width [cm^-1 atm^-1]
  real(8) Strue(nGas,nlinesMAX)    ! true line intensity [cm^-1 / cm^2 molecule^-1]

  save

  !------- program options --------------  
  logical, parameter :: use_beguier_2015 = .false. ! include Beguier+ JQRST (2015) CH4 data in near-IR?

contains

  ! ==================================================================================

  subroutine read_line_data(calc_sigma)

    ! read all line data

    use dimensions, only : datadir, nGas, nlinesMAX

    implicit none

    integer kk, il, ierr, iGas

    logical, intent(in) :: calc_sigma(nGas)                   ! calculate sigma for gas in question

    integer mol_temp      ! molecule number
    integer iso_temp      ! isotope number
    real(8) nu0_temp      ! wavenumber [cm^-1]
    real(8) Sref_temp     ! reference line intensity [cm^-1 / cm^2 molecule^-1]
    real(8) einstein_temp ! Einstein A coefficient [s^-1]
    real(8) gamair_temp   ! Air-broadened half width [cm^-1 atm^-1]
    real(8) gamself_temp  ! Self-broadened half width [cm^-1 atm^-1]
    real(8) Egnd_temp     ! Ground state energy [cm^-1]
    real(8) ncoeff_temp   ! Exponent for pressure broadening T-dependence
    real(8) delta_temp    ! Air-broadened pressure shift [cm^-1/atm]

    
    namelist/crosssec_nml/deltanu_trunc, HITEMP, Sref_cutoff

    ! read lines
    write(*,*) "creating line data..."

    open(90,file='input.nml',status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
       print*, 'Cannot find required input.nml file, aborting.'
       call abort
    else
       read(90,crosssec_nml)
       close(90)
    endif

    molec_loop : do iGas = 1, nGas
       if(calc_sigma(iGas))then

          if(iGas==iGas_H2O)then
             if(HITEMP)then
                open(unit=111,file=trim(datadir)//'HITEMP_H2O.par')
             else
                open(unit=111,file=trim(datadir)//'HITRAN_H2O.par')
             endif
          elseif(iGas==iGas_CO2)then
             if(HITEMP)then
                open(unit=111,file=trim(datadir)//'HITEMP_CO2.par')
             else
                open(unit=111,file=trim(datadir)//'HITRAN_CO2.par')
             endif
          elseif(iGas==iGas_O3)then
             open(unit=111,file=trim(datadir)//'HITRAN_O3.par')
          elseif(iGas==iGas_CH4)then
             open(unit=111,file=trim(datadir)//'HITRAN_CH4.par')
          elseif(iGas==iGas_SO2)then
             open(unit=111,file=trim(datadir)//'HITRAN_SO2.par')
          elseif(iGas==iGas_H2S)then
             open(unit=111,file=trim(datadir)//'HITRAN_H2S.par')
          elseif(iGas==iGas_NH3)then
             open(unit=111,file=trim(datadir)//'HITRAN_NH3.par')             
          else
             write(*,*) 'Note: spectral dataset not found for iGas = ',iGas,' .'
          end if

          il = 1
          read_loop : do 
             if(HITEMP)then
                read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2)',iostat=kk) mol_temp, iso_temp, &
                     nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp!, delta_temp
             else
                read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2, f8.2)',iostat=kk) mol_temp, iso_temp, &
                     nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp, delta_temp
             endif

             if(Sref_temp>Sref_cutoff)then ! don't allow weak lines
                mol(iGas,il)      = mol_temp
                iso(iGas,il)      = iso_temp
                nu0(iGas,il)      = nu0_temp
                Sref(iGas,il)     = Sref_temp
                einstein(iGas,il) = einstein_temp
                gam_air(iGas,il)  = gamair_temp
                gam_self(iGas,il) = gamself_temp
                Egnd(iGas,il)     = Egnd_temp
                ncoeff(iGas,il)   = ncoeff_temp
                delta(iGas,il)    = delta_temp
                il = il + 1
             endif

             if(kk == -1)then
                exit read_loop
             endif
             if(il>nlinesMAX)then
                exit read_loop
             endif

          end do read_loop
          close(111)

          nlines(iGas) = il - 1

          ! special extra to add unassigned near-IR methane lines from
          ! Beguier+, JQRST (2015)
          if(use_beguier_2015 .and. iGas==iGas_CH4)then
             write(*,*) 'Adding Beguier+ (2015) unassigned CH4 lines in near-IR...'
             open(unit=112,file=trim(datadir)//'beguier_2015.par')
             do il=nlines(iGas)+1,nlines(iGas)+nlines_beguier
                read(112,*) nu0_temp, Sref_temp
                mol(iGas,il)      = 6
                iso(iGas,il)      = 1
                nu0(iGas,il)      = nu0_temp
                Sref(iGas,il)     = Sref_temp
                einstein(iGas,il) = -1.0e10
                gam_air(iGas,il)  = 1.8d-2 ! from Beguier+ (2015) Section 3.1.
                gam_self(iGas,il) = 0.0d0
                Egnd(iGas,il)     = 1.0d0
                ncoeff(iGas,il)   = 0.7d0 ! based on eyeball comparison with ncoeff values in adjacent bands
                delta(iGas,il)    = 0.0d0
             end do
             close(112)

             nlines(iGas) = il - 1
          end if


       end if
    end do molec_loop

  end subroutine read_line_data

  ! ==================================================================================

  subroutine line_strength_width(mu_i,p,p_i,T,iGas)

    ! calculate line strengths, widths and frequency pressure shifts
    ! for a given species at a given temperature

    use composition, only : iGas_CO2, iGas_H2O, iGas_H2, iGas_N2
    use fund_consts, only : c, c2, atm, Rstar, Tref

    implicit none

    real(8), intent(in) :: T    ! temperature [K]
    real(8), intent(in) :: p    ! total pressure [Pa]
    real(8), intent(in) :: p_i  ! partial pressure [Pa]
    real(8), intent(in) :: mu_i ! molar mass [g/mol]
    integer, intent(in) :: iGas

    real(8) Tfact, Tfact1       ! temperature-dependent factors []
    real(8) QrefQ               ! total internal partition sum
    !real(8) gi                  ! state independent degeneracy factor

    integer il

    logical, parameter :: do_air_shift = .false. ! calculate pressure-shifted line center
    logical, parameter :: verbose      = .true. ! print more info for diagnostic

    if(T<20.0d0)then
       write(*,*) 'In line_strength_width T = ',T,', exiting.'
       stop
    end if

    if(iGas==iGas_H2 .or. iGas==iGas_N2)then ! todo: generalize this
       if(verbose) write(*,*) 'Treating ', gas_name(iGas), ' as vib-rot inactive.'
    else

       if(verbose) write(*,*) "calculating line strengths for ",gas_name(iGas),"..."

       read_loop : do il = 1, nlines(iGas)

          ! calculate (air) pressure-shifted nu
          if(do_air_shift)then
             !nu_shift(iGas,il) = nu0(iGas,il) + delta(il)*p_i(iGas)/atm
          end if

          ! calculate Lorentz and Doppler linewidths [cm^-1]
          ! Lorentz linewidth is the Lorentz HWHM
          ! Doppler linewidth is the Doppler HWHM / sqrt(ln(2))
          gam_lore(iGas,il) = (Tref/T)**ncoeff(iGas,il) &
               *(gam_air(iGas,il)*(p - p_i)/atm + gam_self(iGas,il)*p_i/atm)
          gam_dopp(iGas,il) = nu0(iGas,il)*sqrt(2*Rstar*T/(mu_i/1.0d3))/c
          ! note mu_i not mu_avg here because it is the mean speed of the molecule doing the
          ! absorbing that counts.

          ! get Total Internal Partition Sum (TIPS) ratio Qref / Q
          ! simple approach based on BD_TIPS data
          ! robust up to about 1000 K for H2O and CO2
          ! up to 500-600 K for CH4
          ! up to about 350 K for O3
          if(iGas==iGas_H2O .or. iGas==iGas_CO2 .or. iGas==iGas_O3 .or. iGas==iGas_CH4)then
             QrefQ = (Tref/T)**1.5d0
          else
             write(*,*) 'Molecule Q(T) needs to be assessed still!'
             stop
          end if

          ! get actual line intensity
          Tfact1          = (1.0d0 - exp(-c2*nu0(iGas,il)/T)) &
               / (1.0d0 - exp(-c2*nu0(iGas,il)/Tref))
          Tfact           = exp( -c2*Egnd(iGas,il)*(1.0d0/T - 1.0d0/Tref) )
          Strue(iGas,il)  = Sref(iGas,il)*Tfact*Tfact1*QrefQ

          ! no temperature scaling for Beguier+ (2015) lines
          if(use_beguier_2015 .and. iGas.eq.iGas_CH4 .and. il>(nlines(iGas)-nlines_beguier)) Strue(iGas,il) = Sref(iGas,il)

       end do read_loop

       nlines(iGas) = il - 1

    end if

  end subroutine line_strength_width

  ! ==================================================================================

  subroutine get_line_abs_cross_section(iGas,nu,sigma,T)

    ! calculate line absorption cross-section vs. nu for a given species
    ! in a given wavenumber range

    use line_profs, only : lorentz, doppler, voigt
    use dimensions, only : nS 

    implicit none

    integer, intent(in)  :: iGas
    real(8), intent(in)  :: nu(nS)     ! wavenumber [cm^-1]
    real(8), intent(in)  :: T          ! temperature [K]: for CO2 chi-factor calculation
    real(8), intent(out) :: sigma(nS)  ! total absorption cross-section [cm2/molec]

    logical mask(nS)
    integer iS, il
    real(8) f_temp, gam_temp

    real(8) nu1  ! start wavenumber [cm^-1]
    real(8) nu2  ! finish wavenumber [cm^-1]

    logical, parameter :: use_voigt   = .true.  ! use voigt function for line profiles
    logical, parameter :: debug       = .false. ! print additional text for diagnostic

    if(debug) write(*,*) "calculating absorption cross sections for ",gas_name(iGas),"..."

    sigma(:) = 0.0d0

    nu1 = nu(1)
    nu2 = nu(nS)

    if(.not.iGas==iGas_H2)then ! todo: transfer update from generate_kmatrix

       ! add the lines, one by one
       ! for now we use all isotopes (with Earth-like abundances!)
       line_loop : do il = 1, nlines(iGas)

          ! this is the very slow part
          trunc : if(nu0(iGas,il)>nu1-deltanu_trunc(iGas) .and. nu0(iGas,il)<nu2+deltanu_trunc(iGas))then

             mask = abs(nu - nu0(iGas,il)) < deltanu_trunc(iGas)

             write_sig_loop : do iS = 1, nS

                ! update sigma only when mask = T
                if(mask(iS))then

                   if(iGas==iGas_CO2)then
                      gam_temp = gam_lore(iGas,il)*chi_factor(nu(iS),nu0(iGas,il),T)
                   else
                      gam_temp = gam_lore(iGas,il) ! no sublorentzian lineshift
                   end if

                   if(use_voigt)then
                      call voigt(gam_temp, gam_dopp(iGas,il), nu(iS)-nu0(iGas,il), f_temp)
                   else
                      call lorentz(nu(iS)-nu0(iGas,il),gam_temp,f_temp)
                      !call doppler(nu(iS)-nu0(iGas,il),gam_dopp(iGas,il),f_temp)
                   endif

                   if(debug)then
                      call voigt(0.0d0,1.0d0,0.0d0,f_temp)
                      ! f_V(0) when Lorentz HWHM = 0: should yield f_D(0) = 1/sqrt(pi)
                      print*,f_temp
                      call voigt(1.0d0,1.0d-8,0.0d0,f_temp)
                      ! f_V(0) when Doppler HWHM -> 0: should yield f_L(0) = 1/pi
                      print*,f_temp
                      stop
                   endif

                   sigma(iS) = sigma(iS) + Strue(iGas,il)*f_temp

                   ! included as a (non-ideal) fix to -ve numbers from the Beguier+ 2015 data.
                   if(sigma(iS)<0.0d0)then
                      sigma(iS) = 0.0d0 
                   endif
                   
                endif
             end do write_sig_loop
          end if trunc
       end do line_loop

    end if

  end subroutine get_line_abs_cross_section

  ! ==================================================================================

  subroutine get_CIA_cross_section(iGas,nu,sigma_CIA,T,p,f_i)

    ! calculate cross-section vs. nu for a given species
    ! in a given wavenumber range

    use dimensions, only  : nGas, nS
    use cia, only         : calculate_CIA
    use cia, only         : interpolateH2Ocont_PPC
    use cia, only         : interpolateH2Ocont_CKD

    implicit none

    integer, intent(in)  :: iGas
    real(8), intent(in)  :: nu(nS)        ! wavenumber [cm^-1]
    real(8), intent(in)  :: T             ! temperature [K]
    real(8), intent(in)  :: p             ! pressure [Pa]
    real(8), intent(in)  :: f_i(nGas)     ! molar concentration [mol/mol]
    real(8), intent(out) :: sigma_CIA(nS) ! CIA absorption cross-section [cm2/molecule of air]

    integer iS
    real(8) sigma_temp                  ! absorption cross-section [cm2/molecule of air]

    sigma_CIA(:) = 0.0d0

    !--------------------------------------------------------------
    ! add any relevant collision-induced absorption (longwave only)

    if(iGas==iGas_H2O)then
       do iS = 1, nS
          sigma_temp = 0.0d0
          call interpolateH2Ocont_CKD(nu(iS),T,p*f_i(iGas_H2O),p*(1.0d0-f_i(iGas_H2O)),sigma_temp,.false.)
          sigma_CIA(iS) = sigma_temp
       end do
    end if

    if(iGas==iGas_CO2)then
       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>0.0d0 .and. nu(iS)<3000.0d0)then
             call calculate_cia('CO2CO2',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_CO2),sigma_temp,.false.)
          end if
          sigma_CIA(iS) = sigma_temp
       end do
    end if

    if(iGas==iGas_N2)then
       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>0.0d0 .and. nu(iS)<500.0d0)then
             call calculate_cia('N2_N2_',nu(iS),T,p,p*f_i(iGas_N2),p*f_i(iGas_N2),sigma_temp,.false.)
          end if
          sigma_CIA(iS) = sigma_temp
       end do
    end if
    
    if(iGas==iGas_H2 .and. iGas_CO2.ne.-1)then ! H2-CO2 CIA gets called once, for H2
       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>10.0d0 .and. nu(iS)<1880.0d0)then
             !call calculate_cia('N2_H2_',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_H2),sigma_temp,.false.)
             !print*,'using N2-H2 for CO2-H2!'
             call calculate_cia('CO2H2_',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_H2),sigma_temp,.false.)
          end if
          sigma_CIA(iS) = sigma_temp
       end do
    end if

    if(iGas==iGas_CH4 .and. iGas_CO2.ne.-1)then ! CO2-CH4 gets called once, for CH4

       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>10.0d0 .and. nu(iS)<1370.0d0)then
             !call calculate_cia('N2_CH4',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_CH4),sigma_temp,.false.)
             !print*,'using N2-CH4 for CO2-CH4!'
             call calculate_cia('CO2CH4',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_CH4),sigma_temp,.false.)
          end if
          sigma_CIA(iS) = sigma_temp
       end do
    end if

  end subroutine get_CIA_cross_section

  ! ==================================================================================

  real(8) function chi_factor(nu_temp,nu_L,T)

    ! calculate sublorentzian line profile chi factor using Perrin & Hartman (1989) assumptions
    ! currently for CO2-dominated atmospheres only
    ! (note they have separate data for CO2-N2 if needed)

    implicit none

    real(8), intent(in) :: T       ! temperature [K]
    real(8), intent(in) :: nu_L    ! line center wavenumber [cm^-1]
    real(8), intent(in) :: nu_temp ! wavenumber [cm^-1]

    ! parameters for the empirical scheme
    ! see Table 3, p. 314
    real(8), parameter, dimension(3) :: alpha = [0.0888d0,  0.0d0,     0.0232d0]
    real(8), parameter, dimension(3) :: beta  = [-0.160d0,  0.0526d0,  0.0d0]
    real(8), parameter, dimension(3) :: epsi  = [0.00410d0, 0.00152d0, 0.0d0]
    real(8), parameter, dimension(3) :: sigPH = [3.0,       30.0,      120.0] ! cm^-1

    real(8), dimension(3) :: B
    real(8) deltanu ! wavenumber separation [cm^-1]

    ! equation (3), p. 314
    B          = alpha + beta*exp(-epsi*T)
    ! no change by default
    chi_factor = 1.0d0

    deltanu   = abs(nu_temp - nu_L)
    ! see Table 2, p. 314
    if(deltanu<sigPH(1))then
       chi_factor = 1.0d0
    elseif(deltanu>sigPH(1) .and. deltanu<sigPH(2))then
       chi_factor = exp(-B(1)*(deltanu - sigPH(1)))
    elseif(deltanu>sigPH(2) .and. deltanu<sigPH(3))then
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1))-B(2)*(deltanu - sigPH(2)))
    else
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1)) - B(2)*(sigPH(3) - sigPH(2)) &
            - B(3)*(deltanu - sigPH(3)))
    end if

    return
  end function chi_factor

  ! ==================================================================================

end module cross_section

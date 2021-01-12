module composition

  !-----------------------------------------------------------------------
  !  
  ! Defines atmospheric composition and calculates relevant thermodynamic 
  ! and radiative quantities.
  !
  ! Robin Wordsworth (2016)
  !
  !-----------------------------------------------------------------------

  ! if cp_air etc. are defined in namelist, composition calculation is overridden

  ! this should be renamed eventually - after all, condensable component is a climate
  ! variable - not just atmosphere
  
  use fund_consts, only : Rstar
  use dimensions, only : nGas
  
  implicit none
  private

  integer, public, save :: ngas_check
  integer, public, save :: vgas

  character(len=3), save, allocatable, dimension(:) :: gas_name
  real(8), save, allocatable, dimension(:)          :: gas_molarconc
  character(len=3), dimension(10)                   :: gas_name_MAX
  real(8),          dimension(10)                   :: gas_molarconc_MAX

  logical :: initialized = .false.

  ! initialize all these variables to -1
  integer, public, save :: iGas_H2   = -1
  integer, public, save :: iGas_He   = -1
  integer, public, save :: iGas_H2O  = -1
  integer, public, save :: iGas_CO2  = -1
  integer, public, save :: iGas_CO   = -1
  integer, public, save :: iGas_N2   = -1
  integer, public, save :: iGas_O2   = -1
  integer, public, save :: iGas_O3   = -1
  integer, public, save :: iGas_SO2  = -1
  integer, public, save :: iGas_H2S  = -1
  integer, public, save :: iGas_CH4  = -1
  integer, public, save :: iGas_NH3  = -1

  !------------ compound properties and thermodynamic constants ---------------
  ! data from http://pubchem.ncbi.nlm.nih.gov/ and Principles of Planetary Climate p. 92, Table 2.1.
  ! where not specified values are defined at the triple point

  real(8), parameter :: mu_SO2      = 64.06        ! molar mass of SO2 [g/mol]
  real(8), parameter :: mu_H2S      = 34.08        ! molar mass of H2S [g/mol]
  real(8), parameter :: mu_CH4      = 16.04        ! molar mass of CH4 [g/mol]
  real(8), parameter :: mu_O3       = 48.00        ! molar mass of O3 [g/mol]
  real(8), parameter :: mu_NH3      = 17.031       ! molar mass of NH3 [g/mol]

  real(8), parameter :: mu_N2       = 28.014       ! molar mass of N2 [g/mol]
  real(8), parameter :: cp_N2       = 1.040        ! constant pressure specific heat of N2 at 300 K [kJ/kg/K]

  real(8), parameter :: mu_H2       = 2.0159       ! molar mass of H2 [g/mol]
  real(8), parameter :: cp_H2       = 14.31        ! constant pressure specific heat of H2 at 300 K [kJ/kg/K]

  real(8), parameter :: mu_H2O      = 18.01528     ! molar mass of H2O [g/mol]
  real(8), parameter :: T0_H2O      = 273.16       ! triple point temperature [K]
  real(8), parameter :: p0_H2O      = 610.78       ! triple point pressure [Pa]
  real(8), parameter :: hlv_H2O     = 2.255e6      ! latent heat of evaporation [J/kg]
  real(8), parameter :: hlf_H2O     = 3.34e5       ! latent heat of fusion {J/kg]
  real(8), parameter :: rho_H2O_ice = 920.0e0      ! density of solid [kg/m3]
  
  real(8), parameter :: R_H2O       = Rstar/mu_H2O ! H2O specific gas constant for thermodynamics.f90 (non-SI units) [J/g/K]

  real(8), parameter :: mu_CO2      = 44.0095      ! molar mass of CO2 [g/mol]
  real(8), parameter :: T0_CO2      = 216.54       ! triple point temperature [K]
  real(8), parameter :: p0_CO2      = 5.185e5      ! triple point pressure [Pa]
  real(8), parameter :: hlv_CO2     = 3.97e5       ! latent heat of evaporation [J/kg]
  real(8), parameter :: hlf_CO2     = 1.96e5       ! latent heat of fusion [J/kg]
  real(8), parameter :: rho_CO2_ice = 1620.0e0     ! density of solid [kg/m3]

  real(8), parameter :: rho_dust    = 2500.0d0     ! density of solid [kg/m3]
  
  public cp_N2, cp_H2, R_H2O
  public rho_H2O_ice, rho_CO2_ice, rho_dust
  public mu_CO2, mu_H2O, mu_N2, mu_H2, mu_SO2, mu_H2S, mu_CH4, mu_O3, mu_NH3

  real(8), public, save :: T0 ! Triple point temperature of variable gas [K]
  real(8), public, save :: p0 ! Triple point pressure of variable gas [Pa]

  real(8), public, save :: rvgas ! specific gas constant [J/kg/K]
  real(8), public, save :: hlv
  real(8), public, save :: hlf
  real(8), public, save :: hls   ! latent heat of sublimation {J/kg]

  logical :: verbose = .true.

  ! all this MAX stuff could be removed eventually
  namelist /composition_nml/ ngas_check, gas_name_MAX, gas_molarconc_MAX

  !-----------------------------------------------------------------------
  public :: composition_init, composition_end, get_vgas_params
  public :: gas_molarconc, gas_name

contains

  subroutine composition_init 

    !-----------------------------------------------------------------------
    !
    ! routine for initializing the planetary parameters
    !
    !-----------------------------------------------------------------------

    integer ierr, iGas, count

    ! read namelist to check number of gases is correct
    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
       print*, 'Cannot find required input.nml file, aborting.'
       call abort
    else
       read(10,composition_nml)
       close(10)
    endif
    if(ngas_check/=nGas)then
       write(*,*) "ngas_check =",ngas_check
       write(*,*) "nGas       =",nGas
       write(*,*) "ngas_check does not match nGas."
       write(*,*) "Use a different composition or modify nGas in params_consts.f90."
       stop
    end if

    ! create gas name and molar concentration arrays
    if (.not.allocated(gas_name))      allocate(gas_name(nGas))
    if (.not.allocated(gas_molarconc)) allocate(gas_molarconc(nGas))

    gas_name(1:nGas)      = gas_name_MAX(1:nGas)
    gas_molarconc(1:nGas) = gas_molarconc_MAX(1:nGas)

    !-----------------------------------------------------------------------

    if(verbose)then
       write(*,*) ''
       write(*,*) 'Atmospheric composition:'
       write(*,*) 'nGas          = ',nGas
       write(*,*) 'gas_name 1    = ',gas_name(1)
       if(nGas>1) write(*,*) 'gas_name 2    = ',gas_name(2)
       if(nGas>2) write(*,*) 'gas_name 3    = ',gas_name(3)
       if(nGas>3) write(*,*) 'gas_name 4    = ',gas_name(4)
       write(*,*) 'gas_molarconc = ',gas_molarconc
       write(*,*) ''
    endif

    vgas = 0
    ! identify the variable gas
    do iGas = 1, nGas
       if(gas_molarconc(iGas)<0.0d0)then
          print*,'Variable gas is ',gas_name(iGas)
          vgas = iGas
          call get_vgas_params
          exit
       end if
    end do
    if(vgas==0)then
       print*,'No variable gas found!'
    end if

    ! assign the 'iGas_X' labels
    count=0
    do iGas=1,nGas
       if (gas_name(iGas)=="H2_") then
          iGas_H2=iGas
          count=count+1
       elseif (gas_name(iGas)=="He_") then
          iGas_He=iGas
          count=count+1
       elseif (gas_name(iGas)=="H2O") then
          iGas_H2O=iGas
          count=count+1
       elseif (gas_name(iGas)=="CO2") then
          iGas_CO2=iGas
          count=count+1
       elseif (gas_name(iGas)=="CO_") then
          iGas_CO=iGas 
          count=count+1
       elseif (gas_name(iGas)=="N2_") then
          iGas_N2=iGas
          count=count+1
       elseif (gas_name(iGas)=="O2_") then
          iGas_O2=iGas
          count=count+1
       elseif (gas_name(iGas)=="O3_") then
          iGas_O3=iGas
          count=count+1
       elseif (gas_name(iGas)=="SO2") then
          iGas_SO2=iGas
          count=count+1
       elseif (gas_name(iGas)=="H2S") then
          iGas_H2S=iGas
          count=count+1
       elseif (gas_name(iGas)=="CH4") then
          iGas_CH4=iGas
          count=count+1
       elseif (gas_name(iGas)=="NH3") then
          iGas_NH3=iGas
          count=count+1
       endif
    enddo

    if(count/=nGas)then
       print*,'Mismatch between nGas and number of recognised gases in composition_init.'
       print*,'Either we haven`t managed to assign all the gases, or there are duplicates.'
       print*,'Please try again.'
       call abort
    endif

    initialized = .true.

  end subroutine composition_init


  ! ==================================================================================
  subroutine get_vgas_params

    use fund_consts, only: Rstar

    ! sets the thermodynamic properties of the variable gas

    if(gas_name(vgas)=='H2O')then
       p0    = p0_H2O
       T0    = T0_H2O
       hlv   = hlv_H2O
       hlf   = hlf_H2O
       rvgas = Rstar/(mu_H2O/1.0e3)
    elseif(gas_name(vgas)=='CO2')then
       p0    = p0_CO2
       T0    = T0_CO2
       hlv   = hlv_CO2
       hlf   = hlf_CO2
       rvgas = Rstar/(mu_CO2/1.0e3)
    else
       print*,'We do not have thermodynamic data for variable gas ',gas_name(vgas),', exiting.'
       call abort
    endif
    
    hls = hlv + hlf

  end subroutine get_vgas_params
  ! ==================================================================================

  subroutine composition_end

    deallocate (gas_name, gas_molarconc)

  end subroutine composition_end

  ! ==================================================================================

end module composition



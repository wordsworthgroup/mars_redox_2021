module fund_consts
  
  !------------------------------------------------------------------
  !  
  ! Defines fundamental constants and parameters needed by the model.
  !
  ! Robin Wordsworth (2016)
  !
  !------------------------------------------------------------------
  
  implicit none

  real(8), parameter :: pi     = 3.14159265d0  
  complex(8), parameter :: i   = dcmplx(0.0d0,+1.0d0)
  real(8), parameter :: Rstar  = 8.3144621       ! ideal gas constant [J/mol/K]
  real(8), parameter :: c      = 2.99792458d8    ! speed of light [m/s]
  real(8), parameter :: p0     = 1013.25d0       ! atm to mbar conversion
  real(8), parameter :: Tref   = 296.d0          ! hitran standard temperature [K]
  real(8), parameter :: mpr    = 1.672621e-27    ! mass of proton [kg]
  real(8), parameter :: h      = 6.6260693d-34   ! Planck constant [J s]
  real(8), parameter :: NA     = 6.0221415d23    ! Avogadro's number [molecules/mol]
  real(8), parameter :: kB     = 1.3806488d-23   ! Boltzmann constant [m2/kg/s2/K]
  real(8), parameter :: atm    = 101325.0d0      ! Pa to atm conversion
  real(8), parameter :: c2     = (h*c/kB)*1.0d2  ! a constant [cm/K]
  real(8), parameter :: losch  = 2.6867774d19    ! Loschmidt's number [molecules cm^-3]
  real(8), parameter :: stefan = 5.670367d-8     ! Stefan-Boltzmann constant [W/m2/K4]
  
end module fund_consts

module dimensions
  
  implicit none

  ! Moderate resolution calculation with three radiatively active species  
  integer, parameter :: nGas      = 3            ! number of species
  integer, parameter :: nLay      = 50           ! number of layers
  integer, parameter :: nLev      = nLay+1       ! number of levels
  integer, parameter :: nAng      = 4            ! number of angles to integrate over
  integer, parameter :: nlinesMAX = 3000000      ! maximum number of lines to read
  integer, parameter :: nS        = 2000         ! spectral resolution (same in sw and lw for now)
  integer, parameter :: nTem      = 1            ! number of cross-section temperatures at each layer
  
  ! Set nTem > 1 to calculate results more accurately (runs more slowly)

  ! directory containing data needed by the model
  character(len=200), save :: datadir='/Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/data/'
  
end module dimensions

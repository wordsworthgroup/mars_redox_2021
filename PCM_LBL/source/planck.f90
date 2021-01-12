module planck
  
  !------------------------------------------
  !  
  ! Calculates the Planck spectral irradiance
  ! for a given wavenumber and temperature
  !
  ! Robin Wordsworth (2016)
  !
  !------------------------------------------

  implicit none
  private

  public :: B_nu
  
contains

  subroutine B_nu(nu_cm,T,Bnu)

    use fund_consts, only : h, c, kB
    !use fund_consts, only : c2

    implicit none

    real(8), intent(in) :: nu_cm ! wavenumber [cm^-1]
    real(8), intent(in) :: T     ! temperature [K]
    real(8), intent(out) :: Bnu  ! Planck spectral irradiance [W/m2/cm^-1/sr]

    !real(8), parameter :: c1 = 2.0d8*h*c**2 ! W/m2/cm^-4/sr
    real(8) nu, Btemp

    nu = 1.0d2*c*nu_cm ! cm^-1 --> Hz

    !Bnu = c1*nu_cm**3/(exp(c2*nu_cm/T) - 1.0d0)
    Btemp = (2*h*nu**3/c**2)/(exp(h*nu/(kB*T))  -  1.0d0) ! W/m2/sr/Hz
    Bnu   = 1.0d2*c*Btemp  ! W/m2/sr/cm^-1
    
  end subroutine B_nu

end module planck

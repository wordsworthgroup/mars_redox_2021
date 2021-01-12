module line_profs
  
  !--------------------------------------------------------
  !  
  ! Calculate Doppler, Lorentz and Voigt lineshapes.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  implicit none

  private
  
  public :: lorentz, doppler, voigt
  
contains

  ! ==================================================================================

  subroutine doppler(delta_nu,gammD,fD)

    use fund_consts, only : pi

    implicit none

    real(8), parameter   :: sqrtpi    = sqrt(pi)
    real(8), intent(in)  :: delta_nu ! (nu-nu0) [cm^-1]
    real(8), intent(in)  :: gammD    ! Doppler linewidth (not HWHM!) [cm^-1]
    real(8), intent(out) :: fD       ! Doppler line profile [1/cm^-1]
    
    fD = exp(-delta_nu**2.0d0 / gammD**2.0d0)/(gammD*sqrtpi)

    return
  end subroutine doppler
  
  ! ==================================================================================

  subroutine lorentz(delta_nu,gammL,fL)

    use fund_consts, only : pi

    implicit none

    real(8), intent(in)  :: delta_nu ! (nu-nu0) [cm^-1]
    real(8), intent(in)  :: gammL    ! Lorentz linewidth [cm^-1]
    real(8), intent(out) :: fL       ! Lorentz line profile [1/cm^-1]

    fL = gammL/(pi*(gammL**2.0d0 + delta_nu**2.0d0))

    return
  end subroutine lorentz

  ! ==================================================================================

  subroutine voigt(gammL, gammD, delta_nu, fV)

    use fund_consts, only : i, pi

    implicit none

    ! inputs / outputs
    real(8), intent(in)  :: gammL    ! Lorentz linewidth [cm^-1]
    real(8), intent(in)  :: gammD    ! Doppler linewidth (not HWHM!) [cm^-1]
    real(8), intent(in)  :: delta_nu ! (nu-nu0) [cm^-1]
    real(8), intent(out) :: fV       ! Voigt profile [1/cm^-1]

    ! general variables
    real(8) x, y
    complex(4) z, w

    x  = delta_nu/gammD
    y  = gammL/gammD
    z  = cmplx(x + i*y)
    w  = humilcek_w4(z)
    fV = real(w)*sqrt(1.0d0/pi)/gammD
    ! see Schreier, JQRST (1992), eqn. (1-4)

    ! note that we define gammD as the Doppler linewidth, not the Doppler HWHM.
    ! this removes the need for a lot of superfluous sqrt(log(2)) terms.

    return
  end subroutine voigt
    
  ! ==================================================================================

  complex(4) function humilcek_w4(z)

    ! Humilcek, JQRST (1982) algorithm for the complex probability function
    ! code copied and pasted (with F77 -> .f90 conversion) from his Appendix.
    ! computes w(z) = exp(-z**2)*erfc(-i*z)
    
    implicit none

    complex(4) z, t, u
    
    real(4) x, y, s
    
    x = real(z)
    y = aimag(z)
    t = cmplx(y,-x)
    s = abs(x) + y

    ! s = |nu-nu0|/gammD + gammL/gammD
    ! large s means Lorentz profile...

    if(s>=15.0e0)then
       ! region I: Lorentz profile basically

       humilcek_w4 = t*0.5641896e0/(0.5e0 + t*t)

    elseif(s<15.0e0 .and. s>=5.5e0)then       
       ! region II
       u = t*t
       
       humilcek_w4 = t*(1.410474 + u*0.5641896e0)/(0.75e0 + u*(3.0e0 + u))

    else!if(s<5.5e0 .and. y>=0.195e0*abs(x) - 0.176e0)then
       ! region III

       humilcek_w4 = (16.4955e0 + t*(20.20933e0 + t*(11.96482e0 + t*(3.778987e0 + t*0.5642236e0)))) / &
            (16.4955e0 + t*(38.82363e0 + t*(39.27121e0 + t*(21.69274e0 + t*(6.699398e0 + t)))))

       ! note we are not using Humilcek's recommended 'Region IV' here - it caused -ve fV values in this
       ! implementation - because of typos in the values in his table perhaps?
       ! plotting vs. matlab results shows close correspondence to exact Voigt results in practice.
       ! should be fine for all but the most sensitive remote sensing applications - not our focus here.

       !u = t*t
       !humilcek_w4 = cexp(u) - &
       !     t*(36183.31e0 - u*(3321.9905e0 - u*(1540.787e0 - u*(219.0313e0 - u*(35.76683e0 - u*(1.320522e0 - u*0.56419e0)))))) / &
       !     (32066.6e0 - u*(24322.84e0 - u*9022.228e0 - u*(2186.181e0 - u*(364.2191e0 - u*(61.57037e0 - u*(1.841439e0 - u))))))

    end if
    
    return
  end function humilcek_w4
  
  ! ==================================================================================

end module line_profs

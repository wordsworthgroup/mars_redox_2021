module atmos_profile
  
  !---------------------------------------------------------------------------------
  !  
  ! Calculates the general moist adiabat.
  !
  ! Robin Wordsworth (2016)
  ! Some of the code in this module developed from routines originally shared by
  ! Emmanuel Marcq (Marcq 2012).
  !
  !---------------------------------------------------------------------------------

  use dimensions,     only  : nLay, nLev
  use thermodynamics, only  : get_psat_H2O, get_Tsat_CO2
  use thermodynamics, only  : therm, tdpsdt, cp_CO2, cp_H2O

  implicit none
  private

  public :: make_Tp_profile

contains  

  ! ==================================================================================

  subroutine make_Tp_profile(condense_CO2,variable_cp,m_n,m_v,Rmn,cp_n_in,Tstra,Ts,ptop,ps,Psurf_v,play,plev,Tlay,Tlev,flay)

    ! this subroutine creates profiles of T, rho_v, rho_n, Pv and Pn following 
    ! a procedure similar to that described in Kasting 1988

    ! it is not very computationally efficient.
    ! but this matters little: even at low spectral resolution the
    ! time-consuming part will be the radiative tranfer calculation.
    
    use fund_consts, only : Rstar
    use composition, only : gas_name, vgas

    implicit none

    integer iLayF, iLay, iLev, iLay_final

    integer, parameter   :: nlayF = 100000 ! number of vertical layers

    logical, intent(in)  :: condense_CO2   ! CO2 bulk condensation
    logical, intent(in)  :: variable_cp    ! allow cp(T)
    real(8), intent(in)  :: m_n            ! molar mass of dry atmospheric component [kg/mol]
    real(8), intent(in)  :: m_v            ! molar mass of moist atmospheric component [kg/mol]
    real(8), intent(in)  :: Rmn            ! specific gas constant of dry atmospheric component [J/kg/K]
    real(8), intent(in)  :: Tstra          ! lowest temperature we can send adiabat to [K]
    real(8), intent(in)  :: Ts             ! surface temperature [K]
    real(8), intent(in)  :: cp_n_in        ! specific heat capacity at constant pressure of non-condensable component [J/kg/K]
    real(8), intent(in)  :: Ptop           ! pressure at TOA [Pa]
    real(8), intent(inout) :: ps           ! total surface pressure [Pa]
    real(8), intent(inout) :: Psurf_v      ! surface partial pressure (condensible species) [Pa]
    real(8), intent(out) :: play(1:nLay)
    real(8), intent(out) :: plev(1:nLev)
    real(8), intent(out) :: Tlay(1:nLay)
    real(8), intent(out) :: Tlev(1:nLev)
    real(8), intent(out) :: flay(1:nLay)

    real(8) lnp1,lnp2,lnpnew
    real(8) dlogp_rcm

    real(8) Dp
    real(8) dlogp, Psat_max

    real(8) T(1:nlayF)                       ! temperature [K]
    real(8) Pv(1:nlayF),Pn(1:nlayF)          ! pressures [Pa]
    real(8) rho_v(1:nlayF), rho_n(1:nlayF)   ! density [kg/m3]
    real(8) f_v(1:nlayF)                     ! = Pv/P [mol/mol]
    real(8) a_v(1:nlayF)                     ! = rho_v/rho_n [kg/kg]
    real(8) mtot(1:nlayF)                    ! = (rho_v+rho_n)/(n_v+n_n) [g/mol]
    real(8) cp_n                             ! specific heat capacity at constant pressure of non-condensable component [J/kg/K]

    integer profil_flag(1:nlayF) ! 0 = dry, 1 = moist, 2 = isothermal

    real(8) Psurf_n     ! surface partial pressure (incondensible species) [Pa]
    real(8) dTdp        ! [K/Pa]
    real(8) dPvdp,dPndp ! [Pa/Pa]
    real(8) psat_v      ! local Psat_H2O value [Pa]
    real(8) Tcond       ! condensation temperature [K]

    logical, parameter :: verbose           = .false.
    logical, parameter :: add_Pvar_to_total = .false.
    ! careful: do not set this to true until we can interpolate cross-sections in pressure space!
    
    ! initialise flags
    profil_flag(:) = 0

    !-----------------------
    ! assign input variables

    Psat_max = 9.0d6 ! maximum vapour pressure [Pa]
    
    cp_n = -1.0d8
    if(.not.variable_cp) cp_n = cp_n_in
    
    if(vgas<1)then
       if(Psat_max>0.0d0)then
          print*,'Must have Psat_max = 0 if no variable species'
          Psat_max = 0.0d0
       endif
    end if

    ! calculate local saturation vapour pressure
    psat_v = Psat_max
    if(vgas>0)then
       if(gas_name(vgas)=='H2O')then
          call get_psat_H2O(Ts,psat_v)
       endif
    endif

    ! Moist adiabat unless greater than or equal to psat_max
    if(Psurf_v<psat_v)then
       profil_flag(1) = 0 ! dry adiabat
    else
       if(verbose)then
          print*,'Surface is super-saturated, decreasing Psurf_v in atmos_profile.'
       end if
       Psurf_v        = psat_v
       profil_flag(1) = 1 ! moist adiabat
    endif

    if(add_Pvar_to_total)then
       Psurf_n = ps
       ps      = Psurf_n + Psurf_v
    else
       Psurf_n = ps - Psurf_v

       if(Psurf_n<0.0d0)then
          write(*,*) 'Error, high atmospheric condensible gas content has'
          write(*,*) 'driven non-condensing partial pressure to negative values!'
          write(*,*) 'Psurf_n = ',Psurf_n,' Pa'
          stop
       endif

    endif

    if(verbose)then
       print*,'Psat_v  =',psat_v,' Pa'
       print*,'Ts      =',Ts,' K'
       print*,'Tstra   =',Tstra,' K'
       print*,'Psurf_v =',Psurf_v,' Pa'
       print*,'Psurf_n =',Psurf_n,' Pa'
       print*,'m_n     =',m_n,' kg/mol'
       print*,'m_v     =',m_v,' kg/mol'
       print*,'Rstar   =',Rstar,' J/mol/K'
    endif

    ! define pressure grid
    dlogp_rcm = -real(log(ps)-log(ptop))/real(nLay)

    play(1)  = ps*exp(dlogp_rcm)
    do iLay=1,nLay-1
       play(iLay+1) = play(iLay)*exp(dlogp_rcm)
    enddo

    plev(1) = ps
    do iLev=2,nLev-1
       ! log-linear interpolation
       iLay = iLev
       plev(iLev) = exp( log( play(iLay)*play(iLay-1) )/2.0d0 )
    enddo

    !-------------------------------
    ! Layer 1
    T(1)     = Ts
    Pv(1)    = Psurf_v
    Pn(1)    = Psurf_n
    rho_n(1) = m_n*Pn(1)/(Rstar*Ts)
    rho_v(1) = m_v*Pv(1)/(Rstar*Ts)
    a_v(1)   = rho_v(1)/rho_n(1)

    ! log pressure grid spacing (constant)
    dlogp = -(log(Pn(1)+Pv(1))-log(ptop))/(nlayF-1)

    ! update cp_n if necessary
    if(variable_cp) then
        write(*,*) 'variable_cp: assuming CO2 is non-condensable constituent.'
        cp_n = cp_CO2(Ts)*1.0d3
    endif

    call gradients_atmos(variable_cp,profil_flag(1),m_n,m_v,Rmn,rho_v(1),rho_n(1),Ts,cp_n,dTdp,dPvdp,dPndp)
    if(verbose)then
       print*, 'dT/dp     ground [K/Pa] =',dTdp
       print*, 'dlnT/dlnp ground [K/Pa] =',dTdp*((Pv(1)+Pn(1))/Ts)
    endif

    ! initial delta p, delta z
    Dp = (Pn(1) + Pv(1))*(exp(dlogp) - 1.0d0)

    !-------------------------------
    ! Layer 2
    T(2)     = Ts + dTdp*Dp 
    Pv(2)    = Pv(1) + dPvdp*Dp
    Pn(2)    = Pn(1) + dPndp*Dp
    rho_n(2) = m_n*Pn(2)/(Rstar*T(2))
    rho_v(2) = m_v*Pv(2)/(Rstar*T(2))
    a_v(2)   = rho_v(2)/rho_n(2)

    !-------------------------------
    ! start vertical ascent
    do iLayF=2,nlayF-1

       ! 1st assume next layer same as last one
       profil_flag(iLayF) = profil_flag(iLayF-1)

       ! update delta p
       Dp = (Pn(iLayF) + Pv(iLayF))*(exp(dlogp) - 1.0d0)

       ! update cp_n if necessary
       if(variable_cp) then
          cp_n = cp_CO2(T(iLayF))*1.0d3
       end if

       ! intial gradients call to calculate temperature at next level
       call gradients_atmos(variable_cp,profil_flag(iLayF),m_n,m_v,Rmn,rho_v(iLayF),rho_n(iLayF),T(iLayF),cp_n,dTdp,dPvdp,dPndp)

       T(iLayF+1) = T(iLayF) + dTdp*Dp

       ! test for moist adiabat at next level
       ! calculate local saturation vapour pressure
       psat_v = psat_max
       if(vgas>0)then
          if(gas_name(vgas)=='H2O')then
             call get_psat_H2O(T(iLayF+1),psat_v)
          endif
       endif

       ! if we are supersaturated, recalculate
       if (psat_v < Pv(iLayF) + dPvdp*Dp) then
          profil_flag(iLayF) = 1
          call gradients_atmos(variable_cp,profil_flag(iLayF),m_n,m_v,Rmn,rho_v(iLayF),rho_n(iLayF),T(iLayF),cp_n,dTdp,dPvdp,dPndp)
       endif

       ! test for stratosphere at next level
       if (T(iLayF+1) <= Tstra) then
          profil_flag(iLayF) = 2
          call gradients_atmos(variable_cp,profil_flag(iLayF),m_n,m_v,Rmn,rho_v(iLayF),rho_n(iLayF),T(iLayF),cp_n,dTdp,dPvdp,dPndp)
       endif

       ! calculate pressures at next level
       Pn(iLayF+1) = Pn(iLayF) + dPndp*Dp
       Pv(iLayF+1) = Pv(iLayF) + dPvdp*Dp

       ! if on the moist adiabat, recalculate Pv to avoid supersaturation
       if(profil_flag(iLayF) == 1)then

          ! calculate local saturation vapour pressure
          psat_v = psat_max 
          if(vgas>0)then
             if(gas_name(vgas)=='H2O')then
                call get_psat_H2O(T(iLayF+1),psat_v)
             endif
          endif
          Pv(iLayF+1) = psat_v

       endif

       ! calculate gas densities at next level (assume ideal)
       rho_n(iLayF+1) = m_n*Pn(iLayF+1)/(Rstar*T(iLayF+1)) 
       select case(profil_flag(iLayF)) 
       case(2) ! isothermal
          rho_v(iLayF+1) = (rho_v(iLayF)/rho_n(iLayF))*rho_n(iLayF+1)
       case(1) ! moist
          rho_v(iLayF+1) = m_v*psat_v/(Rstar*T(iLayF+1))
       case(0) ! dry
          rho_v(iLayF+1) = m_v*Pv(iLayF+1)/(Rstar*T(iLayF+1)) 
       end select

    enddo

    !-------------------------------
    ! save to regular grid 

    ! surface pressure
    ps = Pn(1) + Pv(1)

    ! create f_v, q_v, mtot for saving
    do iLayF=1,nlayF
       f_v(iLayF)  = Pv(iLayF)/(Pv(iLayF) + Pn(iLayF))
       !mtot(iLayF) = kg_2_g*(rho_v(iLayF) + rho_n(iLayF))/(rho_v(iLayF)/m_v + rho_n(iLayF)/m_n)
       mtot(iLayF) = 1.0d3*(rho_v(iLayF) + rho_n(iLayF))/(rho_v(iLayF)/m_v + rho_n(iLayF)/m_n)
    enddo

    ! convert to lower-res grid for radiative code
    Tlay(:) = 0.0d0
    flay(:) = 0.0d0
    iLay = 1
    do iLayF=2,nlayF

       if(iLay<=nLay)then
          ! interpolate rcm variables

          if(Pn(iLayF)+Pv(iLayF) < play(iLay))then

             if(iLayF==1)then
                print*,'Error in atm_profile: Psurf here less than Psurf in main code!'
                call abort
             endif

             lnp1   = log(Pn(iLayF-1)+Pv(iLayF-1))
             lnp2   = log(Pn(iLayF)+Pv(iLayF))
             lnpnew = dble(log(play(iLay)))

             Tlay(iLay) = T(iLayF-1)*(lnp2-lnpnew)/(lnp2-lnp1)    + T(iLayF)*(lnpnew-lnp1)/(lnp2-lnp1)
             flay(iLay) = f_v(iLayF-1)*(lnp2-lnpnew)/(lnp2-lnp1)  + f_v(iLayF)*(lnpnew-lnp1)/(lnp2-lnp1)

             iLay = iLay + 1

          endif

       endif
    enddo

    iLay_final=iLay-1
    if(iLay_final<nLay)then
       if(verbose)then
          print*,'Interpolation in atm_profile stopped at layer',iLay,'!'
       endif

       do iLay=iLay_final+1,nLay
          Tlay(iLay) = Tlay(iLay-1)
          flay(iLay) = flay(iLay-1)
       enddo
    endif

    ! CO2 condensation 'haircut' of temperature profile if necessary
    ! not consistent with moist profile necessarily
    if(condense_CO2)then
       do iLay = 1,nLay
          call get_Tsat_CO2(play(iLay),Tcond)
          if(Tlay(iLay)<Tcond)then
             Tlay(iLay) = Tcond
             if(gas_name(vgas)=='H2O')then
                call get_psat_H2O(Tlay(iLay),psat_v)
                flay(iLay) = psat_v/play(iLay)
             endif
          end if
       end do
    end if

    ! finally, get Tlev via interpolation
    Tlev(:) = 0.0d0
    Tlev(1) = Ts
    do iLev=2,nLev-1
       ! linear interpolation
       iLay = iLev 
       Tlev(iLev) = (Tlay(iLay) + Tlay(iLay-1))/2.0d0 
    enddo
    Tlev(nLev) = Tlay(nLay)

    return
  end subroutine make_Tp_profile
  
  ! ==================================================================================

  subroutine gradients_atmos(variable_cp,profil_flag,m_n,m_v,Rmn,rho_v,rho_n,T,cp_n,dTdp,dPvdp,dPndp)

    use fund_consts, only : Rstar
    use composition, only : gas_name
    use dimensions,  only : nGas
    
    implicit none

    logical, intent(in)  :: variable_cp ! allow cp(T)
    integer, intent(in)  :: profil_flag ! 0 = dry, 1 = moist, 2 = isothermal
    real(8), intent(in)  :: m_n         ! molar mass of non-condensable component [kg/mol]
    real(8), intent(in)  :: m_v         ! molar mass of condensable component [kg/mol]
    real(8), intent(in)  :: Rmn         ! specific gas constant of non-condensable component [J/kg/K]
    real(8), intent(in)  :: rho_v       ! density of condensable component [kg/m3]
    real(8), intent(in)  :: rho_n       ! density of non-condensable component [kg/m3]
    real(8), intent(in)  :: T           ! temperature [K]
    real(8), intent(in)  :: cp_n        ! specific heat capacity at constant pressure of non-condensable component [J/kg/K]
    real(8), intent(out) :: dTdp        ! gradient of temperature wrt pressure [K/Pa]
    real(8), intent(out) :: dPndp       ! gradient of non-cond partial pressure wrt total pressure [Pa/Pa]
    real(8), intent(out) :: dPvdp       ! gradient of cond partial pressure wrt total pressure [Pa/Pa]

    ! internal variables
    real(8) cp_v                        ! specific heat capacity of condensable component [J/kg/K]
    real(8) a_v                         ! ratio of condensable and non-condensable densities []

    real(8) press, rho_plus, rho_minus
    real(8) dlnr, dlna, dpsat, dsv
    real(8) s_minus, s_plus
    real(8) s_v, s_c                   ! specific entropy variables []
    real(8) psat_plus, psat_minus, Pn
    real(8) nul ! dummy variable

    select case(profil_flag)
    case(2) ! isothermal

       dTdp  = 0.0d0
       a_v   = rho_v/rho_n ! constant here
       dPndp = 1.0d0/(1.0d0 + (m_n/m_v)*a_v)
       dPvdp = 1.0d0 - dPndp

    case(1) ! moist

       Pn = rho_n*T*rmn

       if(ngas==1)then
          print*,'Cannot have moist adiabat with one gas...'
          stop
       endif

       if( .not. gas_name(ngas)=='H2O')then
          write(*,*) 'Variable gas type not recognized, exiting'
          stop
       endif
       
       call get_psat_H2O(T - 2.0d-1,psat_minus)
       call get_psat_H2O(T + 2.0d-1,psat_plus)
       call get_psat_H2O(T,press)

       rho_minus = m_v*psat_minus/(Rstar*(T - 2.0d-1))
       rho_plus  = m_v*psat_plus/(Rstar*(T + 2.0d-1))
       
       rho_minus = rho_minus + 1.0e-12
       rho_plus  = rho_plus  + 1.0e-12
       ! therm crashes with low rho_pm values
       
       call therm(T-2.0d-1,rho_minus*1.0d-3,nul,nul,nul,nul,nul,nul,nul,nul,nul,press,s_minus,nul)
       call therm(T+2.0d-1,rho_plus*1.0d-3, nul,nul,nul,nul,nul,nul,nul,nul,nul,press,s_plus,nul)
       s_c = 2.06 * log(T/273.15d0) ! 

       s_plus  = s_plus  * 1.0d3
       s_minus = s_minus * 1.0d3
       s_c     = s_c     * 1.0d3 ! convert to SI

       if(T < 280.0d0)then 
          dpsat = press*( 1730.63d0*log(10.0d0) / (T - 39.714d0)**2.0d0 )
       else
          call tdpsdt(T,dpsat)
          dpsat = dpsat / T
       endif
       
       dsv  = (s_plus - s_minus)/4.0d-1  ! dsv*T = ds / d ln[T]
       s_v  = (s_plus + s_minus)/2.0d0
       dlnr = (T/rho_v) * (rho_plus - rho_minus)/4.0d-1  ! d ln[rho_v] / d ln[T]

       if(rho_n/rho_v<1.0d-5)then
          dlna = -T*dsv/(s_v - s_c)
       else
          a_v  = rho_v/rho_n
          dlna = (rmn*dlnr - cp_n + rmn - a_v*T*dsv)/(a_v*(s_v - s_c) + rmn)  
          ! d ln[alpha_v] / d ln[T]
          ! note cp_n + rmn = cv_n, which is what's required
       endif

       dTdp  = 1.0d0 / (dpsat + rho_n*rmn*(1.0d0 + dlnr - dlna)) ! c.f. Marcq+ (2012) S2.2.2
       dPvdp = dTdp*dpsat
       dPndp = 1.0d0 - dPvdp ! from p = p_v + p_n

    case(0) ! dry

       if(rho_v>1.0d-16)then ! dry runaway case

          if(variable_cp)then
             cp_v  = cp_H2O(T)*1.0d3
          else
             cp_v  = cp_H2O(300.0d0)*1.0d3
          endif

          dTdp  = 1.0d0/(rho_n*cp_n + rho_v*cp_v)
          dPndp = 1.0d0/(1.0d0 + m_n/m_v*rho_v/rho_n)
          dPvdp = 1.0d0/(1.0d0 + m_v/m_n*rho_n/rho_v)
       else                ! entirely dry case
          dTdp  = 1.0d0/(rho_n*cp_n)
          dPndp = 1.0d0
          dPvdp = 0.0d0
       endif  

    end select

  end subroutine gradients_atmos

  ! ==================================================================================

end module atmos_profile





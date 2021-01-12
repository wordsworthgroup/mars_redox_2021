  module atmos_structure

    !---------------------------------------------------------------------------------
    !  
    ! Defines atmospheric temperature, pressure and composition variation with height.
    !
    ! Robin Wordsworth (2016)
    !
    !---------------------------------------------------------------------------------


    use dimensions,      only : nLay, nLev, nGas
    use composition,     only : gas_molarconc,  iGas_H2O
    use thermodynamics,  only : get_Tsat_H2O, get_Tsat_H2O_GB, get_Tsat_CO2, get_psat_H2O

    use atmos_profile,   only : make_Tp_profile

    use fund_consts,     only : Rstar        ! for atmos_profile. a little messy.
    use composition,     only : mu_H2O, vgas ! for atmos_profile. a little messy.
    use thermodynamics,  only : cp_CO2


    implicit none
    private

    logical :: initialized = .false.

    logical, parameter :: strato_isotherm   = .true.   ! isothermal stratosphere
    logical, parameter :: strato_invert     = .false.  ! temperature inversion in upper atmosphere
    logical, parameter :: condense_H2O      = .false.  ! H2O bulk condensation
    logical, parameter :: variable_cp       = .false.  ! allow cp(T)
    logical, parameter :: gray_rad_equi     = .false.  ! gray radiative equilibrium profile for debugging
    logical, parameter :: load_prof         = .false.  ! load initial T and f_H2O profiles from results/
    logical, parameter :: verbose           = .false.

    logical condense_CO2                               ! CO2 bulk condensation
    logical use_atmos_profile                          ! use moist adiabat subroutines in atmos_profile

    real(8) mu_avg_dry                                 ! average dry molar mass [g/mol]
    real(8) R_dry                                      ! dry specific gas constant [J/kg/K]
    real(8) cp_dry                                     ! dry specific heat capacity [J/kg/K]

    !------- general namelist parameters --------
    real(8) ptop                                       ! top-of-atmosphere pressure [Pa]
    real(8) Tstra                                      ! stratosphere temperature (isothermal case) [K]
    real(8) RH                                         ! relative humidity for condensible species []
    real(8) K                                          ! R/cp
    real(8) pv_atmos_profile                           ! surface saturation vapor pressure for atmos_profile [Pa]

    public :: atmos_structure_init, set_Tadia_fvar

    save

  contains

    subroutine atmos_structure_init(verbose_output, kappa_gray, Ts, Tlay, Tlev, play, plev, &
         dp, ps, f_i, mu_i, mu_avg, grav)

      use fund_consts,     only : Rstar, stefan
      use composition,     only : composition_init
      use composition,     only : iGas_CO2, iGas_N2, iGas_H2, iGas_O3, mu_CO2, mu_N2, mu_H2O, mu_H2, mu_O3
      use composition,     only : iGas_CH4, iGas_SO2, iGas_H2S, iGas_NH3, mu_CH4, mu_SO2, mu_H2S, mu_NH3
      use thermodynamics,  only : cp_H2O

      implicit none

      !------- general variables ----------
      integer iLay, iLev, iGas, iFine             ! do loop variables
      integer ierr
      integer, parameter :: nFine = 1000          ! number of fine steps to take per level in vary_cp calculation []
      real(8) psat_tmp                            ! saturation vapour pressure for condensible species [Pa]
      real(8) ps_v                                ! surface pressure of condensible species [Pa]
      real(8) f_trap                              ! cold trap molar concentration for condensible species [mol/mol]
      real(8) Tcond                               ! condensation temperature (to be moved) [K]
      real(8) plev_temp(nLev)                     ! pressure levels [Pa]
      real(8) Ttemp                               ! temporary temperature for variable_cp calculation [K]
      real(8) Ktemp                               ! temporary adiabatic ratio for variable_cp calculation [K]
      real(8) logp1, logp2                        ! assorted parameters for defining T-p profile []
      real(8) log_dT, log_dp                      ! assorted parameters for defining T-p profile []
      real(8) logplev                             ! assorted parameters for defining T-p profile []
      real(8) p0strat                             ! pressure at which stratoinversion begins [Pa]
      real(8) ASR                                 ! absorbed stellar flux [W/m2]
      real(8) tau_i
      real(8) taulay(nLay)
      real(8) taulev(nLev)
      real(8) dtau(nLay)
      real(8) rholay(nLay)                        ! density in layer [kg/m3]
      real(8) f_temp(nLay,nGas)                   ! molar concentration of species [mol/mol]
      real(8) vgas_mask(nGas)                     ! mask to remove vgas from mu dry calculation
      character(128) blah                         ! for reading and discarding output file headers

      logical, intent(in)  :: verbose_output      ! output files (to be read as input here) contain headers
      real(8), intent(in)  :: kappa_gray          ! gray gas mass absorption cross-section for debug [m2/kg]
      real(8), intent(out) :: Tlay(nLay)          ! temperature in layer  [K]
      real(8), intent(out) :: Tlev(nLev)          ! temperature at level  [K]
      real(8), intent(out) :: play(nLay)          ! pressure layers [Pa]
      real(8), intent(out) :: plev(nLev)          ! pressure levels [Pa]
      real(8), intent(out) :: dp(nLay)            ! pressure difference [Pa]
      real(8), intent(out) :: mu_i(nGas)          ! molar mass of species [g/mol]
      real(8), intent(out) :: f_i(nLay,nGas)      ! molar concentration of species [mol/mol]
      real(8), intent(out) :: mu_avg(nLay)        ! average molar mass in layer [g/mol]

      !------- non-general namelist parameters --------
      real(8), intent(out) :: ps                  ! surface pressure [Pa]
      real(8), intent(out) :: Ts                  ! surface temperature [K]
      real(8), intent(out) :: grav                ! gravity [m/s/s]

      namelist/atmos_structure_nml/ps,Ts,ptop,Tstra,grav,RH,condense_CO2,use_atmos_profile,pv_atmos_profile

      ! read namelist
      open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
      if (ierr/=0) then
         print*, 'Cannot find required input.nml file, aborting.'
         call abort
      else
         read(10,atmos_structure_nml)
         close(10)
      endif

      !-------- set up dry molar concentrations, atmospheric mass and heat capacity --------

      if(iGas_CO2/=-1) call set_f_i_mu_i(iGas_CO2,mu_CO2,mu_i(iGas_CO2),f_i(:,iGas_CO2))
      if(iGas_N2/=-1)  call set_f_i_mu_i(iGas_N2, mu_N2, mu_i(iGas_N2), f_i(:,iGas_N2) )
      if(iGas_H2/=-1)  call set_f_i_mu_i(iGas_H2, mu_H2, mu_i(iGas_H2), f_i(:,iGas_H2) )
      if(iGas_O3/=-1)  call set_f_i_mu_i(iGas_O3, mu_O3, mu_i(iGas_O3), f_i(:,iGas_O3) )
      if(iGas_CH4/=-1) call set_f_i_mu_i(iGas_CH4,mu_CH4,mu_i(iGas_CH4),f_i(:,iGas_CH4))
      if(iGas_SO2/=-1) call set_f_i_mu_i(iGas_SO2,mu_SO2,mu_i(iGas_SO2),f_i(:,iGas_SO2))
      if(iGas_H2S/=-1) call set_f_i_mu_i(iGas_H2S,mu_H2S,mu_i(iGas_H2S),f_i(:,iGas_H2S))
      if(iGas_NH3/=-1) call set_f_i_mu_i(iGas_NH3,mu_NH3,mu_i(iGas_NH3),f_i(:,iGas_NH3))

      ! H2O : not constant!
      ! we only set mu_i here: f_i is done later.
      if(iGas_H2O/=-1)then
         mu_i(iGas_H2O)  = mu_H2O
         if(.not.iGas_H2O==vgas)then! .and. use_atmos_profile)then
            write(*,*) 'Please generalize atmos_structure...'
            f_i(:,iGas_H2O) = gas_molarconc(iGas_H2O)
            !stop
         end if
      end if
      
      ! check mu_i and calculate mean molar mass in layer
      if(any(mu_i<1.0d-8))then
         write(*,*) 'Some values of mu_i are zero, check all gases are recognized.'
         stop
      endif

      vgas_mask(:) = 1.0d0
      if(vgas>0) vgas_mask(vgas) = 0.0d0
      do iLay = 1,nlay
         mu_avg(:) = sum(f_i(iLay,:)*mu_i(:)*vgas_mask(:))/sum(f_i(iLay,:)*vgas_mask(:))
      end do

      ! assign mean molar mass values
      ! this section needs to be generalized still.
      mu_avg_dry = mu_avg(1)
      R_dry      = Rstar/(mu_avg_dry/1.0d3)
      cp_dry     = 1040.0d0 ! for N2 by default
      if(variable_cp)then
         cp_dry = -1.0d8
      elseif(iGas_CO2>-1 .and. iGas_N2>-1)then
            cp_dry = cp_CO2(270.0d0)*1.0e3*(mu_CO2/mu_avg_dry)*f_i(nLay,iGas_CO2) + 1040.0d0*(mu_N2/mu_avg_dry)*f_i(nLay,iGas_N2)
      elseif(iGas_CO2>-1)then
         if (gas_molarconc(iGas_CO2)>0.8d0)then
            cp_dry = cp_CO2(270.0d0)*1.0e3
         endif
      elseif(iGas_N2>-1)then
         if(gas_molarconc(iGas_N2)>0.8d0)then
            cp_dry = 1040.0d0
         endif
      elseif(iGas_H2O>-1)then
         if(gas_molarconc(iGas_H2O)>0.8d0)then
            cp_dry = 2080.0d0
            ! TO BE REPLACED
         endif
      else
         print*,'Your cp_dry scheme in atmos_structure is still too primitive...'
         stop
      endif
      
      K = R_dry/cp_dry

      !-------- set up initial pressure and temperature profile --------
      logp1 = log10(ptop)
      logp2 = log10(ps)  
      do iLev = 1,nLev
        logplev    = logp1 + dble(iLev-1)*(logp2 - logp1)/dble(nLev-1)
        plev(iLev) = 1.0d1**logplev
      end do
      ! now reverse the data
      plev_temp = plev
      do iLev = 1,nLev
        plev(iLev) = plev_temp(nLev - iLev + 1)
      end do

      log_dp = log(plev(2)) - log(plev(1))

      do iLay = 1,nLay
        iLev       = iLay
        dp(iLay)   = plev(iLev) - plev(iLev+1)
        !play(iLay) = plev(iLev) - dp(iLay)/2.0d0
        play(iLay) = exp(log(plev(iLev)) + log_dp/2.0d0)
      end do

      if(variable_cp .and. .not.use_atmos_profile) then

        write(*,*) 'WARNING - simple variable_cp currently works for pure composition only...'
        write(*,*) 'working for H2O by default'
        !write(*,*) 'working for CO2 by default'

        Tlev(1) = Ts ! first level is the surface
        do iLev = 1,nLev-1
          Ttemp = Tlev(iLev)
          do iFine = 1,nFine
            Ktemp  = (Rstar/mu_H2O)/cp_H2O(Ttemp)
            !Ktemp  = (Rstar/mu_CO2)/cp_CO2(Ttemp)
            log_dT = (log_dp/dble(nFine))*Ktemp
            Ttemp  = Ttemp*exp(log_dT)
          end do
          ! add to array if above stratosphere temperature, otherwise shut down calculation
          if(Ttemp>Tstra)then
             Tlev(iLev+1) = Ttemp
           else
             Tlev(iLev+1:nLev) = Tstra
             exit
          endif
        end do

        Ttemp   = (Tlev(1)+Tlev(2))/2.0d0
        Ktemp   = (Rstar/mu_H2O)/cp_H2O(Ttemp)
        !Ktemp   = (Rstar/mu_CO2)/cp_CO2(Ttemp)
        Tlay(1) = Ts*(play(1)/ps)**Ktemp ! should be pretty good estimate
        do iLay = 1,nLay-1
          Ttemp = Tlay(iLay)
          do iFine = 1,nFine
            Ktemp  = (Rstar/mu_H2O)/cp_H2O(Ttemp)
            !Ktemp  = (Rstar/mu_CO2)/cp_CO2(Ttemp)
            log_dT = (log_dp/dble(nFine))*Ktemp
            Ttemp  = Ttemp*exp(log_dT)
          end do
          Tlay(iLay+1) = Ttemp
        end do

     else

       Tlay = Ts*(play/ps)**K
       Tlev = Ts*(plev/ps)**K

     end if

      ! set to svp of H2O where necessary, assuming it is the major constituent
      if(condense_H2O)then
        
         if(use_atmos_profile)then
            write(*,*) 'condense_H2O not compatible with use_atmos_profile'
            stop
         end if
         
         call get_Tsat_H2O_GB(ps,Tcond)
         if(Ts<Tcond) Ts = Tcond
         do iLay = 1,nLay
            call get_Tsat_H2O_GB(play(iLay),Tcond)
            if(Tlay(iLay)<Tcond) Tlay(iLay) = Tcond
         end do

         do iLev = 1,nLev
            call get_Tsat_H2O_GB(plev(iLev),Tcond)
            Tlev(iLev) = Tcond
            if(Tlev(iLev)<Tcond) Tlev(iLev) = Tcond
         end do

      end if

      ! set to s.v.p. of CO2 where necessary, assuming it is the major constituent
      if(condense_CO2)then

         if(use_atmos_profile)then
            write(*,*) 'condense_CO2 may not be compatible with use_atmos_profile...'
            !stop
         end if
         
         do iLay = 1,nLay
            call get_Tsat_CO2(play(iLay),Tcond)
            if(Tlay(iLay)<Tcond) Tlay(iLay) = Tcond
         end do

         do iLev = 1,nLev
            call get_Tsat_CO2(plev(iLev),Tcond)
            if(Tlev(iLev)<Tcond) Tlev(iLev) = Tcond
         end do

      end if

      ! create an isothermal stratosphere
      if(strato_isotherm)then
         do iLay = 1,nLay
            if(Tlay(iLay)<Tstra)then
               Tlay(iLay) = Tstra
            end if
         end do
         do iLev = 1,nLev
            if(Tlev(iLev)<Tstra)then
               Tlev(iLev) = Tstra
            end if
         end do
         if(Ts<Tstra)then
            Ts = Tstra
         end if
      end if

      ! create a stratospheric inversion
      if(strato_invert)then
         p0strat = 5.0d2
         do iLay = 1,nLay
            if(play(iLay)<p0strat)then
               Tlay(iLay) = Tstra*(play(iLay)/p0strat)**(-0.3d0)
            end if
            if(Tlay(iLay)>320.0d0)then
               Tlay(iLay) = 320.0d0
            end if
         end do
         do iLev = 1,nLev
            if(plev(iLev)<p0strat)then
               Tlev(iLev) = Tstra*(plev(iLev)/p0strat)**(-0.3d0)
            end if
            if(Tlev(iLev)>320.0d0)then
               Tlev(iLev) = 320.0d0
            end if
         end do
      end if

      ! set up gray radiative equilibrium temperature profile
      if(gray_rad_equi)then
         ASR  = 200.0d0 ! absorbed stellar flux [W/m2]
         dtau = +kappa_gray*dp/(grav*0.5d0)
         taulev(1) = 0.0d0
         do iLev = 2,nLev
            iLay = iLev - 1
            taulev(iLev) = taulev(iLev-1) + dtau(iLay)
         end do
         taulay(:) = taulev(1:nLay) + dtau/2.0d0
         !tau_i  = +kappa_gray*ps/(grav*0.5d0) ! assumes cosa = 0.5
         tau_i  = sum(dtau)

         Ts   = (ASR*(1.0d0 + tau_i/2.0d0)/stefan)**0.25d0
         Tlay = (ASR*(1.0d0 + (tau_i-taulay))/(2.0d0*stefan))**0.25d0
         Tlev = (ASR*(1.0d0 + (tau_i-taulev))/(2.0d0*stefan))**0.25d0
      end if

      ! set up H2O molar concentration when it is condensing
      if(iGas_H2O/=-1)then ! first check if H2O is present
         
         if(gas_molarconc(iGas_H2O)<0.0 .and. .not. use_atmos_profile)then
            
            ! a value of -1 for the molar concentration in input.nml
            ! is used to indicate that the gas should be assigned to
            ! its saturation vapour pressure curve.
            
            f_trap = 1.0d0
            do iLay = 1,nLay

               ! moist adiabat with fixed relative humidity
               call get_psat_H2O(Tlay(iLay),psat_tmp)
               f_i(iLay,iGas_H2O) = RH*(psat_tmp/play(iLay))

               if(iLay==1 .and. f_i(1,iGas_H2O)>1.0d0)then
                  print*,'H2O molar concentration at surface > 1!'
                  stop
               end if

               f_trap = min(f_trap,f_i(iLay,iGas_H2O))        
            end do
            do iLay = 1,nLay ! ensure that f_i does not increase in the stratosphere
               if(Tlay(iLay)<=Tstra)then
                  f_i(iLay,iGas_H2O) = f_trap
               end if
            end do
            mu_i(iGas_H2O) = mu_H2O
         else
            ! well-mixed throughout
            f_i(:,iGas_H2O) = gas_molarconc(iGas_H2O)
            mu_i(iGas_H2O)  = mu_H2O
         end if
      end if

      ! just load T and f_var profiles from a file
      if(load_prof)then
         open(unit=12,file='results/Tlay.out')
         open(unit=14,file='results/flay.out')
         ! note flay is irrelevant to the (single-molecule) cross-section calcs
         ! we load it just to have consistency for the 1st step calculation
         if(verbose_output)then
            read(12,*) blah
            read(14,*) blah
         end if
         do iLay = 1, nLay
            read(12,'(f12.4)') Tlay(iLay)
            if(iGas_H2O/=-1)then
               do iGas = 1,nGas
                  read(14,'(e12.4)') f_temp(iLay,iGas)
               end do
            end if
         end do
         close(12)
         close(14)
         if(iGas_H2O/=-1) f_i(:,iGas_H2O) = f_temp(:,iGas_H2O)
      end if
     
      ! use the general atmos_profile routine to calculate p, T and fH2O
   
      ! once code is generalized:
      ! move this upwards, calculate mu_avg, mu_H2O in a separate routine
      if(use_atmos_profile)then

         ps_v = pv_atmos_profile

         call make_Tp_profile(condense_CO2,variable_cp,mu_avg_dry/1.0d3,mu_H2O/1.0d3,R_dry,cp_dry,&
              Tstra,Ts,ptop,ps,ps_v,play,plev,Tlay,Tlev,f_i(:,vgas))

         ! remember also to recalculate dp afterwards! 
         do iLay = 1,nLay
            iLev       = iLay
            dp(iLay)   = plev(iLev) - plev(iLev+1)
         end do
      end if
      
      initialized = .true.

    end subroutine atmos_structure_init

    ! ==================================================================================

    subroutine set_Tadia_fvar(ps,Ts,Tlay,play,plev,dp,Tlay_adia,Tlev_adia,f_i)

      implicit none

      real(8), intent(inout) :: ps                  ! surface pressure [Pa]
      real(8), intent(in)    :: Ts                  ! surface temperature [K]
      real(8), intent(in)    :: Tlay(nLay)          ! temperature layers [K]
      real(8), intent(inout) :: play(nLay)          ! pressure layers [Pa]
      real(8), intent(inout) :: plev(nLev)          ! pressure levels [Pa]
      real(8), intent(inout) :: dp(nLay)            ! pressure difference [Pa]
      real(8), intent(inout) :: f_i(nLay,nGas)      ! molar concentration of species [mol/mol]
      real(8), intent(out)   :: Tlay_adia(nLay)     ! adiabatic temperature in layer  [K]
      real(8), intent(out)   :: Tlev_adia(nLev)     ! adiabatic temperature at level  [K]


      !------- general variables ----------
      integer iLay, iLev,iStrat                     ! do loop variables
      real(8) Tcond                                 ! condensation temperature (to be moved) [K]
      real(8) psat_tmp                              ! saturation vapour pressure for condensible species [Pa]
      real(8) ps_v                                  ! surface pressure of condensible species [Pa]

      atm_prof_if : if(use_atmos_profile)then

         ps_v = pv_atmos_profile

         ! right now we are non-general in treatment of mu_avg, mu_var etc.
         ! need to create a subroutine to handle all of this
         ! within atmos_structure
         call make_Tp_profile(condense_CO2,variable_cp,mu_avg_dry/1.0d3,mu_H2O/1.0d3,R_dry,cp_dry, &
              100.0d0,Ts,ptop,ps,ps_v,play,plev,Tlay_adia,Tlev_adia,f_i(:,vgas))

         f_i(:,vgas) = RH*f_i(:,vgas)
         
         ! note RH variations apply to the radiative transfer, not the moist adiabat
         ! this is a bit inconsistent, but a reasonable approximation for a 1d model.

         do iLay = 1,nLay
            iLev       = iLay
            dp(iLay)   = plev(iLev) - plev(iLev+1)
         end do
         
      else

         !-------- recalculate adiabatic temperature profile (crude version) --------

         ! use atmos_profile if iterative RG calcs are required
         
         Tlay_adia = Ts*(play/ps)**K
         Tlev_adia = Ts*(plev/ps)**K

         ! note we assume condensing gas is major constituent in each case
         if(condense_H2O)then
            do iLay = 1,nLay
               call get_Tsat_H2O_GB(play(iLay),Tcond)
               if(Tlay_adia(iLay)<Tcond) Tlay_adia(iLay) = Tcond
            end do
            do iLev = 1,nLev
               call get_Tsat_H2O_GB(plev(iLev),Tcond)
               Tlev_adia(iLev) = Tcond
               if(Tlev_adia(iLev)<Tcond) Tlev_adia(iLev) = Tcond
            end do
         end if

         if(condense_CO2)then
            do iLay = 1,nLay
               call get_Tsat_CO2(play(iLay),Tcond)
               if(Tlay_adia(iLay)<Tcond) Tlay_adia(iLay) = Tcond
            end do
            do iLev = 1,nLev
               call get_Tsat_CO2(plev(iLev),Tcond)
               if(Tlev_adia(iLev)<Tcond) Tlev_adia(iLev) = Tcond
            end do
         end if

      end if atm_prof_if

      if(vgas>0)then

         !-------- recalculate stratosphere height --------
         iStrat = nLay
         do iLay = 2,nLay
            if (Tlay_adia(iLay)<Tlay(iLay)) then
               iStrat = iLay
               if(verbose) write(*,*) 'pStrat = ',play(iStrat)/1.0d5,' bars'
               exit
            end if
         end do
         
         !-------- recalculate variable gas molar concentration --------
         if(vgas==iGas_H2O)then

            ! first, set the atmosphere to the same relative humidity everywhere
            do iLay = 1,nLay
               call get_psat_H2O(Tlay(iLay),psat_tmp)
               f_i(iLay,iGas_H2O) = RH*(psat_tmp/play(iLay))
            end do

            ! now initially assume well-mixed in the stratosphere
            f_i(iStrat+1:nLay,iGas_H2O) = f_i(iStrat,iGas_H2O)

            ! now go all the way through the stratosphere, setting profile
            ! to fixed RH whenever necessary
            do iLay = iStrat+1,nLay
               call get_psat_H2O(Tlay(iLay),psat_tmp)
               if(f_i(iLay,iGas_H2O)>RH*(psat_tmp/play(iLay)))then
                  f_i(iLay:nLay,iGas_H2O) = RH*(psat_tmp/play(iLay))
               end if
            end do

            ! it would probably be better to set f_var outside this routine

         else
            write(*,*) 'set_Tadia_fvar not set up for use with given vgas'
            stop
         end if

      end if
         
    end subroutine set_Tadia_fvar

    ! ==================================================================================

    subroutine set_f_i_mu_i(iGas,mu,mu_i,f_i)
        
      implicit none

      integer, intent(in)  :: iGas      ! species number []
      real(8), intent(in)  :: mu        ! molar mass of species (input) [g/mol]
      real(8), intent(out) :: mu_i      ! molar mass of species [g/mol]
      real(8), intent(out) :: f_i(nLay) ! molar concentration of species [mol/mol]

      f_i(:) = gas_molarconc(iGas)
      mu_i   = mu
      
    end subroutine set_f_i_mu_i

    ! ==================================================================================

  end module atmos_structure



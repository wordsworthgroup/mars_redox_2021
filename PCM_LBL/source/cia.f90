module cia

  !--------------------------------------------------------
  !  
  ! Initialise and calculate collision-induced absorption.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  implicit none
  private

  ! public interfaces
  public :: setup_cia, calculate_cia
  public :: interpolateH2Ocont_PPC
  public :: interpolateH2Ocont_CKD

  ! module variables
  character(LEN=*), parameter :: fmat1 = "(A20,F10.3,F10.3,I7,F7.1,E10.3,F8.3)"

  integer, parameter :: nSmax = 3000 ! max. possible spectral intervals
  integer, parameter :: nTmax = 10   ! max. possible temperature intervals
  integer, parameter :: nPmax = 12   ! total number of possible gas pairs

  real(8), save :: wn_arr(nSmax,nPmax)        ! wavenumber [cm^-1]
  real(8), save :: temp_arr(nTmax,nPmax)      ! temperature [K]
  real(8), save :: abs_arr(nSmax,nTmax,nPmax) ! absorption coefficient [cm^-1/amagat^2]

  integer, save :: nS(nPmax)    ! number of spectral intervals for each pair
  integer, save :: nT(nPmax)    ! number of temperature intervals for each pair

  real(8), parameter :: cm2m = 1.0d-2 ! to be moved elsewhere

  logical :: initialized(nPmax) = .false.

contains

  ! ==================================================================================

  subroutine setup_cia(pairtype)

    use dimensions, only : datadir
    use fund_consts, only: losch

    ! Set up collision-induced absorption parameters used by the radiation code.

    implicit none

    ! input & output
    character(6), intent(in) :: pairtype          ! gas pair type
    character(200) dt_file
    integer ios

    character(20) bleh
    real(8) blah, Ttemp
    integer nres

    integer ipair,iS,iT

    logical, parameter :: verbose = .false.

    logical, parameter :: use_turbet_data = .true.

    ipair = 0

    write(*,*) '-----------------------------------------------------'
    write(*,*) 'Initialising ',pairtype,' continuum from HITRAN database...'

    ! match pair number to pair type
    if(pairtype=='H2__H2')then
       ipair     = 1
       nT(ipair) = 10
       nS(ipair) = 2428
    elseif(pairtype=='H2__He')then
       ipair     = 2
       nT(ipair) = 10
       nS(ipair) = 2428
    elseif(pairtype=='N2_N2_')then
       ipair     = 3
       nT(ipair) = 10
       nS(ipair) = 582
    elseif(pairtype=='N2_H2_')then
       ipair     = 4
       nT(ipair) = 10
       nS(ipair) = 1914
    elseif(pairtype=='N2_CH4')then
       ipair     = 5
       nT(ipair) = 10
       nS(ipair) = 1407
    elseif(pairtype=='CO2H2_')then
       ipair     = 6
       if(use_turbet_data.eqv..true.)then
         nT(ipair) = 6
         nS(ipair) = 300
       else
         nT(ipair) = 4
         nS(ipair) = 601
       endif
    elseif(pairtype=='CO2CH4')then
       ipair     = 7
       nT(ipair) = 4
       nS(ipair) = 601
    elseif(pairtype=='CO2CO2')then
       ipair     = 8
       nT(ipair) = 7
       nS(ipair) = 3000
       !call get_baranov(200.0d0,1200.0d0,blah,.true.)
    elseif(pairtype=='H2OH2O')then
       ipair     = 9
       nS(ipair) = 1001
       nT(ipair) = 11
       call interpolateH2Ocont_CKD(1000.0d0,400.0d0,1.0d4,0.0d4,blah,.true.)
    else
       write(*,*) 'Type ',pairtype,' unknown in setup_cia.'
       call abort
    endif
    initialized(ipair) = .true.

    ! Open the ASCII files
    if(ipair==1)then
       dt_file=trim(datadir)//'cia/H2-H2_norm_2011.cia'
    elseif(ipair==2)then
       dt_file=trim(datadir)//'cia/H2-He_norm_2011.cia'
    elseif(ipair==3)then
       dt_file=trim(datadir)//'cia/N2-N2_2011.cia'
    elseif(ipair==4)then
       dt_file=trim(datadir)//'cia/N2-H2_2011.cia'
    elseif(ipair==5)then
       dt_file=trim(datadir)//'cia/N2-CH4_2011.cia'
    elseif(ipair==6)then
      if(use_turbet_data)then
         !dt_file=trim(datadir)//'cia/co2-h2_turbet2019_interp.dat'
         dt_file=trim(datadir)//'cia/CO2-H2_TURBET2020_V2.cia'
      else
         dt_file=trim(datadir)//'cia/CO2-H2_200_250_300_350.cia'
      endif
    elseif(ipair==7)then
         dt_file=trim(datadir)//'cia/CH4-CO2_200_250_300_350_new.cia'
    elseif(ipair==8)then

      ! HITRAN 2012 database only includes induced dipole data up to 250 cm^-1
      ! so here we use an interpolation created from Gruszka and Borysow data from 
      ! 0 to 500 cm^-1 and Baranov et al. (2004) data from 700 to 2000 cm^-1
      ! It should be pretty reliable for temperatures in the 200-350 K range.
      ! Bets are off outside of that regime.

      ! References: 
      !
      ! Wordsworth, R., Forget, F. and Eymet, V., 2010. Infrared collision-induced 
      ! and far-line absorption in dense CO2 atmospheres. Icarus, 210(2), pp. 992-997.
      !
      ! Baranov, Y.I., Lafferty, W.J. and Fraser, G.T., 2004. Infrared spectrum of the 
      ! continuum and dimer absorption in the vicinity of the O2 vibrational fundamental 
      ! in O2/CO2 mixtures. Journal of Molecular Spectroscopy, 228(2), pp. 432-440.
      !
      ! Gruszka, M. and Borysow, A., 1997. Roto-translational collision-induced absorption 
      ! of CO2 for the atmosphere of Venus at frequencies from 0 to 250 cm− 1, at 
      ! temperatures from 200 to 800 K. Icarus, 129(1), pp.172-177.

         dt_file=trim(datadir)//'cia/CO2-CO2_interp.cia'

    endif

    if(ipair<6 .or. (ipair == 6 .and. use_turbet_data.eqv..true.) .or. ipair == 8)then 

       open(33,file=dt_file,form='formatted',status='old',iostat=ios)
       if (ios/=0) then        ! file not found
          write(*,*) 'Error from setup_cia'
          write(*,*) 'Data file ',trim(dt_file),' not found.'
          write(*,*) 'Check that your path is defined correctly.'
          call abort
       else

          do iT=1,nT(ipair)

             read(33,fmat1) bleh,blah,blah,nres,Ttemp
             if(verbose)then
                write(*,*) 'Resolution given in file: ',trim(dt_file),' is ',nres,'.'
             endif
             if(nS(ipair)/=nres)then
                write(*,*) 'nS(ipair) = ',nS(ipair),', we have a problem. Exiting.'
                call abort
             endif
             temp_arr(iT,ipair)=Ttemp

             do iS=1,nres
                read(33,*) wn_arr(iS,ipair),abs_arr(iS,iT,ipair)
             end do

          end do

       endif
       close(33)

       if(ipair == 6 .and. use_turbet_data.eqv..true.)then
         abs_arr(1:nS(ipair),:,ipair) = abs_arr(1:nS(ipair),:,ipair)
       else
         abs_arr(1:nS(ipair),:,ipair) = (abs_arr(1:nS(ipair),:,ipair)*losch)*losch
         ! cm^5/molecule^2 * (molecule/cm3 / amagat)^2 = cm^-1/amagat^2
       endif

       !call calculate_cia(pairtype,200.0d0,250.0d0,1.0d5,0.95d5,5.0d03,blah,.true.)
       call calculate_cia(pairtype,200.0d0,250.0d0,1.0d5,1.0d5,1.0d05,blah,.true.)
      
    elseif(ipair == 6 .or. ipair == 7)then

       ! CO2CH4 and CO2H2_ data from Andrei Vigasin is stored in a slightly different
       ! format. Also it has units of cm^-1/amagat^2 not cm^5/molecule^2.

       temp_arr(1:4,ipair) = [200.0d0, 250.0d0, 300.0d0, 350.0d0]
       wn_arr(:,ipair)     = 0.0d0
       abs_arr(:,:,ipair)  = 0.0d0

       open(33,file=dt_file,form='formatted',status='old',iostat=ios)
       if (ios/=0) then        ! file not found
          write(*,*) 'Error from setup_cia'
          write(*,*) 'Data file ',trim(dt_file),' not found.'
          write(*,*) 'Check that your path is defined correctly.'
          call abort
       else

          if(verbose)then
             write(*,*) 'Resolution set in code for: ',trim(dt_file),' is ',nS(ipair),'.'
          endif

          do iS=1,nS(ipair)
             read(33,*) wn_arr(iS,ipair),abs_arr(iS,1,ipair), &
                 abs_arr(iS,2,ipair),abs_arr(iS,3,ipair),abs_arr(iS,4,ipair)
          end do

       endif
       close(33)

       ! no need for conversion here as the data is already in cm^-1/amagat^2
       call calculate_cia(pairtype,400.0d0,250.0d0,1.0d5,0.95d5,5.0d03,blah,.true.)
       
    end if

    return
  end subroutine setup_cia

  ! ==================================================================================

  subroutine calculate_cia(pairtype,wn,temp,pres,presA,presB,abcoef,verbose)

    ! Calculate collision-induced absorption for a given gas pair, using lookup tables from HITRAN 2012

    use fund_consts, only: kB, atm

    implicit none

    ! input & output
    logical,      intent(in)  :: verbose
    character(6), intent(in)  :: pairtype          ! gas pair type
    real(8),      intent(in)  :: wn                ! wavenumber [cm^-1]
    real(8),      intent(in)  :: temp              ! temperature [K]
    real(8),      intent(in)  :: pres              ! total pressure [Pa]
    real(8),      intent(in)  :: presA             ! partial pressure A [Pa]
    real(8),      intent(in)  :: presB             ! partial pressure B [Pa]
    real(8),      intent(out) :: abcoef            ! absorption coefficient [cm^2/molecule of air]

    real(8) amagatA
    real(8) amagatB
    real(8) n_molec_cm3

    integer ipair

    ! Loschmit's number (molecule cm^-3 at STP) 
    ! converts cm^5 molecule^-2 --> cm^-1 amagat^-2
    ! see Richard et al. 2011, JQSRT for details

    ! match pair number to pair type
    ! note: this section needs optimization
    if(pairtype=='H2_H2_')then
       ipair     = 1
    elseif(pairtype=='H2_He_')then
       ipair     = 2
    elseif(pairtype=='N2_N2_')then
       ipair     = 3
    elseif(pairtype=='N2_H2_')then
       ipair     = 4
    elseif(pairtype=='N2_CH4')then
       ipair     = 5
    elseif(pairtype=='CO2H2_')then
       ipair     = 6
    elseif(pairtype=='CO2CH4')then
       ipair     = 7
    elseif(pairtype=='CO2CO2')then
       ipair     = 8
    else
       write(*,*) 'Type ',pairtype,' unknown in calculate_cia.'
       call abort
    endif

    if(initialized(ipair).eqv..false.)then
       print*,'cia function called before initialization, exiting.'
       stop
    endif
    
    ! derive abcoef in cm^-1/amagat^2
    call bilinear_cia(wn_arr(:,ipair),temp_arr(:,ipair), &
         abs_arr(:,:,ipair),wn,temp,nS(ipair),nT(ipair),abcoef)

    ! convert to cm^-1
    amagatA = (273.15d0/temp)*(presA/atm)
    amagatB = (273.15d0/temp)*(presB/atm)
    abcoef  = (abcoef)*amagatA*amagatB
    
    ! convert to cm^2/molecule (of air)
    n_molec_cm3 = pres/(kB*temp)/1.0d6 ! molecule/cm^3
    abcoef      = abcoef/n_molec_cm3

    ! unlike for Rayleigh scattering, we do not currently weight by the BB function
    ! however our bands are normally thin, so this is no big deal.

    if(verbose)then
       write(*,*) 'At wavenumber ',wn,' cm^-1'
       write(*,*) 'temperature ',temp,' K'
       write(*,*) 'partial pressure A',presA,' Pa'
       write(*,*) 'partial pressure B',presB,' Pa'
       write(*,*) 'for gas pair ',pairtype
       write(*,*) 'We have ',amagatA,' amagats of ',pairtype(1:3)
       write(*,*) 'and ',amagatB,' amagats of ',pairtype(4:6)
       write(*,*) 'So the absorption is ',abcoef,' cm^2/molecule'
       write(*,*) '                  or ',abcoef*n_molec_cm3*1.0d2,' m^-1'
    endif

    return
  end subroutine calculate_cia

  ! ==================================================================================

  subroutine bilinear_cia(x_arr,y_arr,f2d_arr,x_in,y_in,nX,nY,f)

    ! necessary for interpolation of continuum data
    ! does not need to know pairtype

    implicit none

    integer, intent(in) :: nX, nY

    integer i,j,a,b
    integer, parameter :: nXmax = nSmax
    integer, parameter :: nYmax = nTmax

    real(8), intent(in)  :: x_in,y_in
    real(8), intent(in)  :: x_arr(nXmax)
    real(8), intent(in)  :: y_arr(nYmax)
    real(8), intent(in)  :: f2d_arr(nXmax,nYmax)
    real(8), intent(out) :: f

    real(8) x,y,x1,x2,y1,y2
    real(8) f11,f12,f21,f22,fA,fB

    logical, parameter :: verbose = .false.

    x = x_in
    y = y_in
    a = -10
    b = -10

    ! 1st check we're within the wavenumber range, exit if not
    if ((x<x_arr(2)).or.(x>x_arr(nX-2))) then
       if(verbose)then
          write(*,*) 'Warning from bilinear_cia:'
          write(*,*) 'Outside continuum wavenumber range!'
       endif
       f = 0.0D+0
       return
    endif

    ! in the x (wavenumber) direction 1st
    i  = 1
    x1 = 0.0d0
    x2 = 0.0d0
    x_loop : do while (i<(nX+1))
       if (x_arr(i)>x) then
          x1 = x_arr(i-1)
          x2 = x_arr(i)
          a  = i-1
          i  = 9999
       endif
       i = i + 1
    end do x_loop

    ! now check we're within the temperature range, just warn if not
    if ((y<y_arr(1)).or.(y>y_arr(nY))) then
       if(verbose)then
          write(*,*) 'Warning from bilinear_cia:'
          write(*,*) 'Outside continuum temperature range!'
       endif

       if(y<y_arr(1))then
          y = y_arr(1)  + 0.01d0
       endif
       if(y>y_arr(nY))then
          y = y_arr(nY) - 0.01d0
       endif
    endif

    ! in the y (temperature) direction 2nd
    j  = 1
    y1 = 0.0d0
    y2 = 0.0d0
    y_loop : do while (j<(nY+1))
       if (y_arr(j)>y) then
          y1 = y_arr(j-1)
          y2 = y_arr(j)
          b  = j-1
          j  = 9999
       endif
       j = j + 1
    end do y_loop

    f11 = f2d_arr(a,b)
    f21 = f2d_arr(a+1,b)
    f12 = f2d_arr(a,b+1)
    f22 = f2d_arr(a+1,b+1)

    ! first in the x direction
    fA  = f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)
    fB  = f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1)

    ! then in the y direction
    f   = fA*(y2-y)/(y2-y1) + fB*(y-y1)/(y2-y1)

    return
  end subroutine bilinear_cia

  ! ==================================================================================

  subroutine interpolateH2Ocont_PPC(wn,temp,presS,abcoef)

    ! Calculates the H2O continuum opacity, using the formulae
    ! provided in Pierrehumbert, PPC (2010). As this is based on
    ! the CKD continuum, it provides a useful check for the 
    ! implementation of the more general interpolateH2Ocont_CKD.F90.

    use fund_consts, only: mpr
    use composition, only: mu_H2O

    implicit none

    ! input
    real(8) wn                 ! wavenumber             [cm^-1]
    real(8) temp               ! temperature            [K]
    real(8) presS              ! self-pressure          [Pa]
    !real(8) presF              ! foreign (air) pressure [Pa]

    ! parameters
    real(8), parameter :: T0 = 296.0
    real(8), parameter :: p0 = 1.0d+4

    ! variables
    real(8) x

    ! output
    real(8) abcoef             ! absorption coefficient [m^-1]

    logical, parameter :: debug = .false.

    x = wn - 2500.0d0

    if(debug)then
       write(*,*) '----------------------------------------------------'
       write(*,*) 'Testing H2O continuum...'

       write(*,*) 'interpolateH2Ocont: At wavenumber ',wn,' cm^-1'
       write(*,*) '   temperature ',temp,' K'
       write(*,*) '   H2O pressure ',presS,' Pa'

       if(wn>500.0d0 .and. wn<1400.0d0)then
          abcoef = exp(12.167 - 0.050898*wn + 8.3207e-5*wn**2 - 7.0748e-8*wn**3 + 2.3261e-11*wn**4)*(T0/temp)**4.25*(presS/p0)
       elseif(wn>2100.0d0 .and. wn<3000.0d0)then
          abcoef = exp(-6.0055 - 0.0021363*x + 6.4723e-7*x**2 - 1.493e-8*x**3 &
               + 2.5621e-11*x**4 + 7.328e-14*x**5)*(T0/temp)**4.25*(presS/p0)
       else
          abcoef = 0.0
       endif

       write(*,*) 'The H2O self-absorption is ',abcoef,' m2/kg or '
       abcoef = abcoef*mpr*mu_H2O*1.0d4 ! convert m2/kg --> cm2/molecule
       write(*,*) abcoef, ' cm2/molecule'

    else

       if(wn>500.0d0 .and. wn<1400.0d0)then
          abcoef = exp(12.167 - 0.050898*wn + 8.3207e-5*wn**2 - 7.0748e-8*wn**3 + 2.3261e-11*wn**4)*(T0/temp)**4.25*(presS/p0)
       elseif(wn>2100.0d0 .and. wn<3000.0d0)then
          abcoef = exp(-6.0055 - 0.0021363*x + 6.4723e-7*x**2 - 1.493e-8*x**3 &
               + 2.5621e-11*x**4 + 7.328e-14*x**5)*(T0/temp)**4.25*(presS/p0)
       else
          abcoef = 0.0d0
       endif

       abcoef = abcoef*mpr*mu_H2O*1.0d4 ! convert m2/kg --> cm2/molecule

    endif

    return
  end subroutine interpolateH2Ocont_PPC

  ! ==================================================================================

  subroutine interpolateH2Ocont_CKD(wn,temp,presS,presF,abcoef,firstcall)

    ! Calculates the H2O continuum opacity, using a lookup table from Clough (2005).

    ! the self-continuum scales as ( Rself/Ro )
    ! the foreign continuum scales as ( (Rtot-Rself)/Ro ) 
    ! where R is the density rho [ R = (P/Po)*(To/T) ]. 

    use dimensions,  only : datadir
    use fund_consts, only : atm

    implicit none

    ! input
    real(8) wn                 ! wavenumber             [cm^-1]
    real(8) temp               ! temperature            [K]
    real(8) presS              ! self-pressure          [Pa]
    real(8) presF              ! foreign (air) pressure [Pa]

    ! output
    real(8) abcoef             ! absorption coefficient [m^-1]

    integer nS,nT
    parameter(nS=1001)
    parameter(nT=11)

    real(8) amagatS, amagatF, abcoefS, abcoefF
    real(8) wn_arr(nS)
    real(8) temp_arr(nT)
    real(8) abs_arrS(nS,nT)
    real(8) abs_arrF(nS,nT)
    real(8) data_tmp(nT)

    integer k
    logical firstcall

    save wn_arr, temp_arr, abs_arrS, abs_arrF ! read by master

    character(200) dt_file
    integer ios

    logical, parameter :: debug = .false.

    amagatS = (273.15/temp)*(presS/atm)
    amagatF = (273.15/temp)*(presF/atm)

    ! todo: replace 'firstcall' with a modular init function 
    if(firstcall)then
       print*,'----------------------------------------------------'
       print*,'Initialising H2O continuum from MT_CKD data...'

       ! open the ASCII files

       ! nu array
       dt_file=TRIM(datadir)//'cia/h2o/H2O_CONT_NU.dat'
       open(33,file=dt_file,form='formatted',status='old',iostat=ios)
       if (ios.ne.0) then        ! file not found
          write(*,*) 'Error from interpolateH2O_cont'
          write(*,*) 'Data file ',trim(dt_file),' not found.'
          call abort
       else
          do k=1,nS
             read(33,*) wn_arr(k)
          enddo
       endif
       close(33)

       ! self broadening
       dt_file=TRIM(datadir)//'cia/h2o/H2O_CONT_SELF.dat'
       open(34,file=dt_file,form='formatted',status='old',iostat=ios)
       if (ios.ne.0) then        ! file not found
          write(*,*) 'Error from interpolateH2O_cont'
          write(*,*) 'Data file ',trim(dt_file),' not found.'
          call abort
       else
          do k=1,nS
             read(34,*) data_tmp
             abs_arrS(k,1:nT)=data_tmp(1:nT)
          end do
       endif
       close(34)

       ! foreign (N2+O2+Ar) broadening
       dt_file=TRIM(datadir)//'cia/h2o/H2O_CONT_FOREIGN.dat'
       open(35,file=dt_file,form='formatted',status='old',iostat=ios)
       if (ios.ne.0) then        ! file not found
          write(*,*) 'Error from interpolateH2O_cont'
          write(*,*) 'Data file ',trim(dt_file),' not found.'
          call abort
       else
          do k=1,nS
             read(35,*) data_tmp
             abs_arrF(k,1:nT)=data_tmp(1:nT)
          end do
       endif
       close(35)

       temp_arr(1)  = 200.0d0
       temp_arr(2)  = 250.0d0
       temp_arr(3)  = 300.0d0
       temp_arr(4)  = 350.0d0
       temp_arr(5)  = 400.0d0
       temp_arr(6)  = 450.0d0
       temp_arr(7)  = 500.0d0
       temp_arr(8)  = 550.0d0
       temp_arr(9)  = 600.0d0
       temp_arr(10) = 650.0d0
       temp_arr(11) = 700.0d0

       write(*,*) 'interpolateH2Ocont: At wavenumber ',wn,' cm^-1'
       write(*,*) '   temperature  ',temp,' K'
       write(*,*) '   H2O pressure ',presS,' Pa'
       write(*,*) '   air pressure ',presF,' Pa'
    endif

    call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arrS,wn,temp,abcoefS)!,ind)
    if(debug) write(*,*) 'the self absorption is ',abcoefS,' cm^2 molecule^-1'

    call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arrF,wn,temp,abcoefF)!,ind)
    if(debug) write(*,*) 'the foreign absorption is ',abcoefF,' cm^2 molecule^-1'


    abcoef = abcoefS*amagatS + abcoefF*amagatF ! Eq. (15) in Clough (1989)
    abcoef = abcoef*(presS/(presF + presS))    ! take H2O mixing ratio into account
    ! abs coeffs are given per molecule of H2O

    return
  end subroutine interpolateH2Ocont_CKD

  ! ==================================================================================

  subroutine bilinearbig(nX,nY,x_arr,y_arr,f2d_arr,x_in,y_in,f)

    ! Necessary for interpolation of continuum data
    ! based on some code by A. Spiga (01/2013) 

    implicit none

    integer, intent(in)    :: nX, nY   
    real(8), intent(in)    :: x_in, y_in
    real(8), intent(out)   :: f
    real(8), intent(in)    :: x_arr(nX)
    real(8), intent(in)    :: y_arr(nY)
    real(8), intent(in)    :: f2d_arr(nX,nY)
    real(8), save          :: x, y

    integer i, j, b

    real(8) x1, x2, y1, y2
    real(8) f11, f12, f21, f22

    integer ind    

    logical, parameter :: verbose = .false.

    ind = -9999

    x = x_in
    y = y_in

    ! important to optimize here because the array is quite large
    ! ... and actually calculations only need to be done once
    if ( ind == -9999) then
       ! 1st check we're within the wavenumber range
       if ((x.lt.x_arr(2)).or.(x.gt.x_arr(nX-2))) then
          ind=-1
       else
          i=1
          x2=x_arr(i)
          do while ( x2 .le. x )
             x1=x2
             i=i+1
             x2=x_arr(i)
             ind=i-1
          end do
       endif
    endif

    ! either we already saw we are out of wavenumber range
    ! ... and we just have to set f=0 and exit
    if ( ind == -1) then 
       f = 0.0D+0
       return
       ! or we already determined ind -- so we just proceed
    else
       x1 = x_arr(ind)
       x2 = x_arr(ind+1)
    endif

    ! ... and for y within the temperature range
    if(y.lt.y_arr(1))then
       if(verbose)then
          write(*,*) 'Warning from bilinearbig:'
          write(*,*) 'Below continuum temperature range!'
          write(*,*) 'Setting absorption to minimum T value.'
       endif
       y = y_arr(1) + 0.01d0
       b = 1
       y1 = y_arr(1)
       y2 = y_arr(2)
    elseif(y.gt.y_arr(nY))then
       if(verbose)then
          write(*,*) 'Warning from bilinearbig:'
          write(*,*) 'Above continuum temperature range!'
          write(*,*) 'Setting absorption to maximum T value.'
       endif
       y  = y_arr(nY) - 0.01d0
       b  = nY - 1
       y1 = y_arr(nY-1)
       y2 = y_arr(nY)
    else
       j  = 1
       b  = 1
       y1 = 0.0d0
       y2 = y_arr(j)
       do while ( y2 .le. y )
          y1 = y2
          j  = j + 1
          y2 = y_arr(j)
          b  = j - 1
       end do
    endif

    f11 = f2d_arr(ind,b)
    f21 = f2d_arr(ind+1,b)
    f12 = f2d_arr(ind,b+1)
    f22 = f2d_arr(ind+1,b+1)

    call bilinear(f,f11,f21,f12,f22,x,x1,x2,y,y1,y2)

    return
  end subroutine bilinearbig

  ! ==================================================================================

  subroutine bilinear(f,f11,f21,f12,f22,x,x1,x2,y,y1,y2)

    ! Used for interpolation of continuum data
    
    ! 3 x interpolation routines... merging and optimizing in future would be great.
    
    implicit none

    real(8), intent(in)  :: x, y, x1, x2, y1, y2
    real(8), intent(in)  :: f11, f12, f21, f22
    real(8), intent(out) :: f
    real(8) fA, fB
    
    ! 1st in x-direction
    fA = f11*(x2 - x)/(x2 - x1) + f21*(x - x1)/(x2 - x1)
    fB = f12*(x2 - x)/(x2 - x1) + f22*(x - x1)/(x2 - x1)

    ! then in y-direction
    f = fA*(y2 - y)/(y2 - y1) + fB*(y - y1)/(y2 - y1)

    return
  end subroutine bilinear

  ! ==================================================================================

end module cia




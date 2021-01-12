module temperature_k_grid

  !--------------------------------------------------------
  !  
  ! Handles the temperature interpolation of vib-rot
  ! cross-sections for both longwave and shortwave.
  !
  ! Robin Wordsworth (2016)
  !
  !--------------------------------------------------------

  use dimensions,  only : nS, nLay, nTem

  implicit none
  private

  public :: setup_T_k_grid, calculate_T_k_grid, save_T_k_grid 

  !------- general variables ----------
  integer iLay, iTem                           ! do loop variables
  integer iT1(nLay), iT2(nLay)
  real(8) T_k_grid(nLay,nTem)                  ! temperature array for cross-sections [K]

  !------- program options --------------  
  logical, parameter :: verbose  = .false.      ! more output for debug

  !------- namelist parameters --------
  real(8) dTlay                                ! temperature shift for hot/cold cross-sections [K]

  namelist/temperature_k_grid_nml/dTlay

  save

contains

  ! ==================================================================================

  subroutine setup_T_k_grid(Tlay,T_k_grid_out)

    implicit none

    integer ierr

    real(8), intent(in)  :: Tlay(nLay)              ! temperature in layer  [K]
    real(8), intent(out) :: T_k_grid_out(nLay,nTem) ! temperature array for cross-sections [K]


    ! read namelist 
    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
       print*, 'Cannot find required input.nml file, aborting.'
       call abort
    else
       read(10,temperature_k_grid_nml)
       close(10)
    endif
    
    ! calculate temperature grid on which cross-sections are defined
    do iLay = 1,nLay
       do iTem = 1,nTem
          T_k_grid(iLay,iTem) = Tlay(iLay) + dble(iTem-1-(nTem-1.)/2)*dTlay
       end do
    end do

    ! default values for iT1 and iT2
    iT1(:) = 1
    iT2(:) = iT1(:) + 1

    T_k_grid_out = T_k_grid

  end subroutine setup_T_k_grid

  ! ==================================================================================

  subroutine calculate_T_k_grid(Tlay,iT1_out,T_lin_weight)

    implicit none

    real(8), intent(in)  :: Tlay(nLay)            ! temperature in layer  [K]
    integer, intent(out) :: iT1_out(nLay)         ! T-grid array points for linear interpolation [] 
    real(8), intent(out) :: T_lin_weight(nS,nLay) ! temperature weigting for linear interpolation [] 

    iT1_out(:) = -1
    
    ! either just select 1st values in array, or
    ! interpolate over T grid to get cross-section
    ! at actual atmospheric temperature

    ! find lower / upper T values at each layer

    do iLay = 1,nLay

       find_iT : do 

          if(Tlay(iLay)<T_k_grid(iLay,iT1(iLay)))then
             ! move to lower T_k_grid values
             iT1(iLay) = iT1(iLay) - 1
             ! halt program if Tlay is below lowest T_k_grid value
             if(iT1(iLay)==0)then
                write(*,*) 'Temperature at iLay = ',iLay,' too cold.'
                write(*,*) 'Tlay(iLay)          = ',Tlay(iLay),' K.'
                write(*,*) 'min T_k_grid        = ',T_k_grid(iLay,1),' K.'
                stop
             end if
          elseif(Tlay(iLay)>T_k_grid(iLay,iT1(iLay)+1))then
             ! move to higher T_k_grid values
             iT1(iLay) = iT1(iLay) + 1
             ! halt program if Tlay is above highest T_k_grid value
             if(iT1(iLay)==nTem)then
                write(*,*) 'Temperature at iLay = ',iLay,' too hot.'
                write(*,*) 'Tlay(iLay)          = ',Tlay(iLay),' K.'
                write(*,*) 'max T_k_grid        = ',T_k_grid(iLay,nTem),' K.'
                !stop
             end if
          else
             ! T_k_grid values are correct, exit
             exit
          end if

       end do find_iT

    end do
    iT2(:) = iT1(:) + 1

    do iLay = 1,nLay
       T_lin_weight(:,iLay) = (Tlay(iLay) - T_k_grid(iLay,iT1(iLay)))/(T_k_grid(iLay,iT2(iLay)) - T_k_grid(iLay,iT1(iLay))) 
    end do

    if(verbose)then
       do iLay = 1,nLay
          write(*,'(a15,i3)')    'iT1          = ', iT1(iLay)
       end do
       do iLay = 1,nLay
          write(*,'(a15,f8.3)')  'Tlay         = ', Tlay(iLay)
       end do
       do iLay = 1,nLay
          write(*,'(a15,f8.3)')  'T_lin_weight = ', T_lin_weight(1,iLay)
       end do
    end if

    iT1_out(:) = iT1

  end subroutine calculate_T_k_grid

  ! ==================================================================================

  subroutine save_T_k_grid(verbose_output,play)
    
    implicit none

    logical, intent(in) :: verbose_output ! output files contain headers with variable info and units 
    real(8), intent(in) :: play(nLay)     ! pressure layers [Pa]

    !-------- save results to output file --------

    open(unit=2,file='results/p_k_grid.out')
    open(unit=3,file='results/T_k_grid.out')
    if(verbose_output)then
       write(2,*) 'Temperature grid for cross-section array [K]'
       write(3,*) 'Pressure grid for cross-section array [Pa]'
    end if
    do iLay = 1, nLay
       do iTem = 1, nTem
          write(2,'(e12.4)') play(iLay)
          write(3,'(e12.4)') T_k_grid(iLay,iTem)
       end do
    end do
    close(2)
    close(3)
    
  end subroutine save_T_k_grid

  ! ==================================================================================

  
end module temperature_k_grid

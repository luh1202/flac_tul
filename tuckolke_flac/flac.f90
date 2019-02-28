
! --------- Flac ------------------------- 

subroutine flac

include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

! Update Thermal State
! Skip the therm calculations if itherm = 3

if(time-time_t.gt.dtmax_therm/10 .and. iynts.ne.1) call fl_therm

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Update stresses by constitutive law (and mix isotropic stresses)
call fl_rheol

! Heat Dike Injection 
if (ny_inject.gt.0 .and. iynts.ne.1) then
  call fl_injectheat
endif

! update stress boundary conditions
if (ynstressbc.eq.1.) call bc_update

! Calculations in a node: forces, balance, velocities, new coordinates
call fl_node

! New coordinates
call fl_move

! Adjust real masses due to temperature
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses

! Adjust inertial masses or time step due to deformations
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass

! Adjust time Step 
call dt_adjust


!if (mod(nloop,100).eq.0) &
!     write(*,*) nloop, 'Dev XX',0.25*(stress0(1,1,1,98)+stress0(1,2,1,98)+stress0(1,3,1,98)+stress0(1,4,1,98))-stressI(98,1), &
!     0.25*(stress0(1,1,1,100)+stress0(1,2,1,100)+stress0(1,3,1,100)+stress0(1,4,1,100))-stressI(100,1)
!if (mod(nloop,100).eq.0) &
!     write(*,*) nloop, 'Tot XX',0.25*(stress0(1,1,1,98)+stress0(1,2,1,98)+stress0(1,3,1,98)+stress0(1,4,1,98)), &
!     0.25*(stress0(1,1,1,100)+stress0(1,2,1,100)+stress0(1,3,1,100)+stress0(1,4,1,100))

!if (mod(nloop,100).eq.0) &
!     write(*,*) nloop, 'Dev ZZ',0.25*(stress0(2,1,1,98)+stress0(2,2,1,98)+stress0(2,3,1,98)+stress0(2,4,1,98))-stressI(98,1), &
!     0.25*(stress0(2,1,1,100)+stress0(2,2,1,100)+stress0(2,3,1,100)+stress0(2,4,1,100))-stressI(100,1)
!if (mod(nloop,100).eq.0) &
!     write(*,*) nloop, 'Tot ZZ',0.25*(stress0(2,1,1,98)+stress0(2,2,1,98)+stress0(2,3,1,98)+stress0(2,4,1,98)), &
!     0.25*(stress0(2,1,1,100)+stress0(2,2,1,100)+stress0(2,3,1,100)+stress0(2,4,1,100))

!if (mod(nloop,100).eq.0) write(*,*) ' '

500 continue

return
end

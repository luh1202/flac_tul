!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

! Update Thermal State
! Skip the therm calculations if itherm = 3

if((time-time_t.gt.dtmax_therm/10) .and. (iynts.ne.1)) call fl_therm
!call fl_therm
if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Changing marker phases
! XXX: change_phase is slow, don't call it every loop
!Hao commented following line to skip phase change 5.15.2018
if( mod(nloop, 10).eq.0 ) call change_phase_dike
! Update stresses by constitutive law (and mix isotropic stresses)
!depl = 0.
call fl_rheol

if (ny_inject.gt.0 .and. iynts.ne.1) then
  call fl_injectheat
end if

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

500 continue


return
end subroutine flac

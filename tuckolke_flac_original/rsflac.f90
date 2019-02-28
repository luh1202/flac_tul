
! Re-starting FLAC

subroutine rsflac
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


parameter( kindr=4, kindi=2 )

real(kindr), allocatable :: dum2(:,:),dum3(:,:,:),dum4(:,:,:,:)
real*8 rtime, rdt

! Try to open 'restart.rec'. 
! If file does not exist - restart from last record in '_contents.rs'
! If exists - read nrec from it.
open(1,file='restart.rec', status='old', err=10)
read(1,*) nrec
close(1)
goto 20

10 continue
open( 1, file='_contents.rs', status='old' )
do while (.TRUE.)
    read( 1, *, end=30 ) nrec, nloop, time_my
end do
30 close(1)
goto 40

20 continue
open( 1, file='_contents.rs', status='old' )
do while (.TRUE.)
    read( 1, *, end=60 ) nrecf, nloop, time_my
    if( nrecf .eq. nrec ) goto 50
end do
60 call SysMsg('RESTART: could not find record number given in RESTART.REC in _CONTENTS.RS')
stop

50 endfile(1)
goto 30


40 continue

! Read time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
read (1,rec=nrec) rtime, rdt
close (1)
time = rtime
dt = rdt
time_t = time

! Coordinates and velocities
allocate( dum3(nz,nx,2) )

nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum3
close (1)
cord(1:nz,1:nx,1:2) = dum3(1:nz,1:nx,1:2)
open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum3
close (1)
vel(1:nz,1:nx,1:2) = dum3(1:nz,1:nx,1:2)

deallocate( dum3 )

! Strain
allocate( dum3(3,nz-1,nx-1) )

nwords = 3*(nz-1)*(nx-1)

open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum3
close (1)
strain(1:3,1:nz-1,1:nx-1) = dum3(1:3,1:nz-1,1:nx-1)

deallocate( dum3 )


! Stress
allocate( dum4(4,4,nz-1,nx-1) )

nwords = 4*4*(nz-1)*(nx-1)

open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum4
close (1)
stress0(1:4,1:4,1:nz-1,1:nx-1) = dum4(1:4,1:4,1:nz-1,1:nx-1)

deallocate( dum4 )


! 2-D (nx*nz) arrays - nodes defined
allocate( dum2(nz,nx) )

nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
temp(1:nz,1:nx) = dum2(1:nz,1:nx)

deallocate( dum2 )


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
open (1,file='phasez.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
phasez(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)
!  Read distribution of the phases from the dat file
if (irphase .gt. 0) then
open(12,file='phasedat.inp')
read(12,*) nphasl
do 333 k=1,nphasl
read(12,*) lphase(k)
333 continue
do 332 i=1,nx-1
do 332 j=1,nz-1
read(12,*) ii,jj,phasez(j,i)
332  continue
!    call SysMsg('INIT_PHASE: Read phases from file: Option not ready yet!')
!    stop 21
close(12)
endif

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(i,j,phasez(j,i))
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do

! Markers
open (1,file='rmarker.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
rmarker(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
aps(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
source(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

deallocate( dum2 )


! Pressure at the bottom: pisos 
if( nyhydro .eq. 2 ) then
    open(1,file='pisos.rs')
    read(1,*) pisos
    close (1)
endif

! Calculate AREAS (Important: phasez is needed to calculate area!)
call init_areas

! Distribution of REAL masses to nodes
call rmasses

! Boundary conditions
call init_bc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass

return
end

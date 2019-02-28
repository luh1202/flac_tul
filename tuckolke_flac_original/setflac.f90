
! Setup some parameters (rmass,amass,initial stress,vel,viscosity)

subroutine setflac
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


nloop = 0
mloop = 0
time = 0.

! Mesh generator
call init_cord

! Initial accumulated plastic strain
aps = 0

! Initial injection counter
!inj_count = 0

! Phases in the mesh
call init_phase

! Inverse Areas of triangles
call init_areas

! Initiate temperature field
call init_temp

! Calculation of the initial STRESSES (as hydrostatic)
call init_stress

! Setup boundary conditions
call init_bc

!------------Debugging
!do j=1,nz
!  write(*,*) j,cord(j,1,2), bc(j,1,1), ncod(j,1,1)
!enddo
!------------Debugging

! Distribution of REAL masses to nodes
call rmasses

! Initialization of viscosity
if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic and maxwell)
call dt_mass
dt = min( dt_elastic, dt_maxwell )
time_t = 0

return
end

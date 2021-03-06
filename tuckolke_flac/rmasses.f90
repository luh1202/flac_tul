
!  Distribution of real masses in nodes

subroutine rmasses
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


c1d12 = 1./12.


!   Calcualtion of the TRUE GRAVITATIONAL ZONE MASSES
!-----------------------------------
! THE area(n,it) is inverse of "real" DOUBLE area (=1./det) =>
! area (n,it) ( in program) = 1./(2 real_area)
! real_area = 0.5* (1./area(n,t))
!-----------------------------------

rmass = 0
!$DIR LOOP_PARALLEL
!$DIR LOOP_PRIVATE(j,dens)

do 1 i = 1, nx-1
    do 1 j = 1, nz-1

        !  Area and densities of zones
        dens = Eff_dens( i,j )

        ! Distribution 1/3 of the mass of each element to the nodes 
        ! *0.5 - becuase 1/area is double of real area; *0.5 - 2 grids
        ! (1) Element A:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(1,j,i)*dens
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(1,j,i)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(1,j,i)*dens 
        ! (2) Element B:
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(2,j,i)*dens 
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(2,j,i)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(2,j,i)*dens 

        ! (3) Element C:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(3,j,i)*dens
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(3,j,i)*dens 
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(3,j,i)*dens 

        ! (4) Element D:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(4,j,i)*dens
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(4,j,i)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(4,j,i)*dens 
1 continue

return
end


! Move grid and adjust stresses due to rotation

subroutine fl_move
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

      
! Move Grid
if (movegrid .eq. 0) return


! UPDATING COORDINATES
!$DIR PREFER_PARALLEL
cord = cord + vel*dt

! Diffuse topography
if( topo_kappa.gt.0. .OR. bottom_kappa.gt.0. ) call diff_topo    

!--- Adjusting Stresses And Updating Areas Of Elements
!$DIR PREFER_PARALLEL
do 13  i = 1,nx-1
    do 13  j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! Velocities
        vx1 = vel (j  ,i  ,1)
        vy1 = vel (j  ,i  ,2)
        vx2 = vel (j+1,i  ,1)
        vy2 = vel (j+1,i  ,2)
        vx3 = vel (j  ,i+1,1)
        vy3 = vel (j  ,i+1,2)
        vx4 = vel (j+1,i+1,1)
        vy4 = vel (j+1,i+1,2)

        ! (1) Element A:
        oldvol = 1./2/area(1,j,i)
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        area(1,j,i) = 1./det
        dvol(1,j,i) = det/2/oldvol - 1

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0 (1,1,j,i)
        s22 = stress0 (2,1,j,i)
        s12 = stress0 (3,1,j,i)
        stress0 (1,1,j,i) = s11 + s12*2.*dw12
        stress0 (2,1,j,i) = s22 - s12*2.*dw12
        stress0 (3,1,j,i) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain (1,j,i)
        s22 = strain (2,j,i)
        s12 = strain (3,j,i)
        strain (1,j,i) = s11 + s12*2.*dw12
        strain (2,j,i) = s22 - s12*2.*dw12
        strain (3,j,i) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        oldvol = 1./2/area(2,j,i)
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        area(2,j,i) = 1./det
        dvol(2,j,i) = det/2/oldvol - 1

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0 (1,2,j,i)
        s22 = stress0 (2,2,j,i)
        s12 = stress0 (3,2,j,i)
        stress0 (1,2,j,i) = s11 + s12*2.*dw12
        stress0 (2,2,j,i) = s22 - s12*2.*dw12
        stress0 (3,2,j,i) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        oldvol = 1./2/area(3,j,i)
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        area(3,j,i) = 1./det
        dvol(3,j,i) = det/2/oldvol - 1

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0 (1,3,j,i)
        s22 = stress0 (2,3,j,i)
        s12 = stress0 (3,3,j,i)
        stress0 (1,3,j,i) = s11 + s12*2.*dw12
        stress0 (2,3,j,i) = s22 - s12*2.*dw12
        stress0 (3,3,j,i) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        oldvol = 1./2/area(4,j,i)
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        area(4,j,i) = 1./det
        dvol(4,j,i) = det/2/oldvol - 1

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0 (1,4,j,i)
        s22 = stress0 (2,4,j,i)
        s12 = stress0 (3,4,j,i)
        stress0 (1,4,j,i) = s11 + s12*2.*dw12
        stress0 (2,4,j,i) = s22 - s12*2.*dw12
        stress0 (3,4,j,i) = s12 + dw12*(s22-s11)

13 continue
 
return
end


!============================================================
! Diffuse topography
!============================================================
subroutine diff_topo
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension dh(mnx+1)

if( topo_kappa .gt. 0. ) then             
    xtgtmax = 0.
    dzmax = 0.
    do i = 2, nx-1
        xtgt = abs((cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1))) 
        xtgtmax = max(xtgt,xtgtmax)
        snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dh(i) = topo_kappa * dt * snder
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
    end do

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
        call SysMsg('DIFF_TOPO: divergence!')
        stop
    endif
    if(xtgtmax.gt.0.4) then
    do i = 2, nx-1
        cord(1,i,2) = cord(1,i,2) + dh(i)
    end do
    endif
    cord(1,1 ,2) = cord(1,2   ,2)
    cord(1,nx,2) = cord(1,nx-1,2)

endif


! diffuse also bottom boundary
if( bottom_kappa .gt. 0. ) then
    
    dzmax = 0.
    dhsum = 0.
    do i = 2, nx-1
        snder = ( (cord(nz,i+1,2)-cord(nz,i  ,2))/(cord(nz,i+1,1)-cord(nz,i  ,1)) - &
            (cord(nz,i  ,2)-cord(nz,i-1,2))/(cord(nz,i  ,1)-cord(nz,i-1,1)) ) / &
            (cord(nz,i+1,1)-cord(nz,i-1,1))
        dh(i) = bottom_kappa * dt * snder
        dhsum = dhsum + dh(i)
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
    end do

    new_cord = cord(nz,2,2) + dh(2)
    dh(1) = new_cord - cord(nz,1,2)
    dhsum = dhsum + dh(1)
    if( dabs(dh(1)).gt.dzmax) dzmax = dabs(dh(1))

    new_cord = cord(nz,nx-1,2) + dh(nx-1)
    dh(nx) = new_cord - cord(nz,nx,2)
    dhsum = dhsum + dh(nx)
    if( dabs(dh(nx)).gt.dzmax) dzmax = dabs(dh(nx))

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
        call SysMsg('DIFF_TOPO: divergency at the bottom!')
        stop
    endif

    deltah = dhsum/nx
    do i = 1, nx
        cord(nz,i,2) = cord(nz,i,2) + dh(i) - deltah
    end do

endif
       
return
end


subroutine diff_topo_old
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension dh(mnx+1)

             
dzmax = 0.

do i = 2, nx-1
    dx = cord(1,i,1)-cord(1,i-1,1)
    snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
              (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
    dh(i) = topo_kappa * dt * snder
    if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
end do

if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
    call SysMsg('DIFF_TOPO: divergency!')
    stop
endif

do i = 2, nx-1
    cord(1,i,2) = cord(1,i,2) + dh(i)
end do


! diffuse also bottom boundary
if( bottom_kappa .ne. 0. ) then
    do i = 2, nx-1
        dx = cord(nz,i,1)-cord(nz,i-1,1)
        snder = (cord(nz,i+1,2)-2*cord(nz,i,2)+cord(nz,i-1,2))/4/dx/dx
        dh(i) = bottom_kappa * dt * snder
        if( dabs(dh(i)).gt.dzmax) dzmax = dabs(dh(i))
    end do

    if( dzmax*(nz-1) .gt. dabs(rzbo) ) then
        call SysMsg('DIFF_TOPO: divergency!')
        stop
    endif

    do i = 2, nx-1
        cord(nz,i,2) = cord(nz,i,2) + dh(i)
    end do

endif
       
return
end

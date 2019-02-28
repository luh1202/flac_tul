!======================================================================
  subroutine particle_move
! Move particles every time step
! G. Ito 3/07
!======================================================================

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

!integer iorder(np)
!real*8 xpord(np)
!--------------------------------------------------------------------------
!Interpolate velocities to particles and move them
!--------------------------------------------------------------------------
dxp=rxbo/dble(np)
xpmax=x0
xpmin=x0+rxbo


do 100 i=1,np
  vxp=999.
  vzp=999.
  if (xp(i).lt.cord(1,1,1)) then		!particle outside left side of nodes
    vxp=vel(1,1,1)
    vzp=vel(1,1,2)
  elseif (xp(i).gt.cord(1,nx,1)) then		!particle outside right side of nodes
    vxp=vel(1,nx,1)
    vzp=vel(1,nx,2)
  else						!particle within bounds of nodes
    do 110 ii=1,nx-1
      xl=cord(1,ii,1)
      xr=cord(1,ii+1,1)
      if (xp(i).ge.xl.and.xp(i).le.xr) then
        vxl=vel(1,ii,1)
        vxr=vel(1,ii+1,1)
        vzl=vel(1,ii,2)
        vzr=vel(1,ii+1,2)
        vxp=vxl + ((vxr-vxl)/(xr-xl)) * (xp(i)-xl)
        vzp=vzl + ((vzr-vzl)/(xr-xl)) * (xp(i)-xl)
        exit
      endif
    110 continue
  endif
  
  xp(i)=xp(i) + vxp*dt
  zp(i)=zp(i) + vzp*dt
  
  if (xp(i).gt.xpmax) then
    ipmax=i
    xpmax=xp(i);
  endif
  if (xp(i).lt.xpmin) then
    ipmin=i
    xpmin=xp(i)
  endif
  if (vxp.eq.999.or.vzp.eq.999.0) then
    write(*,*) 'Problems in particle_move, setting particle velocity'
    stop
  endif
!  write(*,*) xp(i), zp(i)
100  continue

!--------------------------------------------------------------------------
!Re-seed near the axis if particle out of original model domain
!--------------------------------------------------------------------------
!write(*,'(2i4,5E15.3)') ipmin, ipmax, xp(ipmin), xp(ipmax), x0, rxbo, dxp

if (xp(ipmax).gt.(x0+rxbo+dxp)) then
  zp(ipmax)=999.0
  xp(ipmax)=cord(1,((nx-1)/2)+1,1)+0.5*dxp  ! M.Behn 9/12/07
  do 200 ii=1,nx-1
    xl=cord(1,ii,1)
    xr=cord(1,ii+1,1)
    if (xp(ipmax).ge.xl.and.xp(ipmax).le.xr) then
      zl=cord(1,ii,2)
      zr=cord(1,ii+1,2)
      zp(ipmax)=zl + ((zr-zl)/(xr-xl)) * (xp(ipmax)-xl)
      exit
    endif
200 continue
  if (zp(ipmax).eq.999.) then
    write(*,*) 'Problems in particle_move, re-seeding particle off the right side'
    stop
  endif
endif

if (xp(ipmin).lt.(x0-dxp).and.x0.lt.0.0) then
  zp(ipmin)=999.0
  xp(ipmin)=cord(1,((nx-1)/2)+1,1)-0.5*dxp  ! M.Behn 9/12/07
  do 300 ii=1,nx-1
    xl=cord(1,ii,1)
    xr=cord(1,ii+1,1)
    if (xp(ipmin).ge.xl.and.xp(ipmin).le.xr) then
      zl=cord(1,ii,2)
      zr=cord(1,ii+1,2)
      zp(ipmin)=zl + ((zr-zl)/(xr-xl)) * (xp(ipmin)-xl)
      exit
    endif
300 continue
  if (zp(ipmin).eq.999.) then
    write(*,*) 'Problems in particle_move, re-seeding particle off the right side'
    stop
  endif
  
endif

! Testing sortrx
!do i=1,np/2
!  xpord(i)=xp(np-i+1)
!  xpord(np/2+i)=xp(i)
!enddo
!  xpord(10)=xpord(11)
!
!call sortrx(np,xpord,iorder)
!
!do i=1,np
!  ii=iorder(i)
!  write(*,*) xpord(ii)/1000, xp(i)/1000, iorder(i)
!enddo
!stop

return
end


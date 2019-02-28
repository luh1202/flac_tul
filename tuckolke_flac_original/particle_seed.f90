!======================================================================
  subroutine particle_seed
! Initial locations of surface particles
! G. Ito 3/07
!======================================================================

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dxp=rxbo
do i = 1,nzonx
  dxzone=rxbo*sizez_x(i)/nelz_x(i)/npelem
  if (dxzone.lt.dxp) dxp=dxzone
enddo
np=rxbo/dxp
  
if (mod(np,2).ne.0) then
  np=np+1
endif

if (np.gt.mnp) then
  write(*,*) 'Too many particles, np > mnp :', np, mnp
endif

write(*,*) '****Surface Particles for RIDGE models: np=',np

dxp=rxbo/dble(np)

do 100 i=1,np/2
   xp(i)=x0+dxp*(i-1)
   xp(np-i+1)=x0+rxbo-dxp*(i-1)
100 continue

do 210 i=1,np
  do 200 ii=1,nx-1
    xl=cord(1,ii,1)
    xr=cord(1,ii+1,1)
    if (xp(i).ge.xl.and.xp(i).le.xr) then
      zl=cord(1,ii,2)
      zr=cord(1,ii+1,2)
      zp(i)=zl + ((zr-zl)/(xr-xl)) * (xp(i)-xl)
      exit
    endif
  200 continue
210  continue


return
end

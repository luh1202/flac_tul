!=======================================================
! HYDROTHERMAL ALTERATION OF THERMAL DIFFUSIVITY
! Modified from Luc's code, G. Ito 8/7/06
!=======================================================
!-------------------------------------------------------------
subroutine ReadHydro()
!-------------------------------------------------------------
include 'precision.inc'
include 'params.inc'
!common /hydroth/ xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2

open( 9, file='hydrother.inp',status='old',err=2001 )
if_hydro = 1
call AdvanceToNextInputLine( 9 )
read (9,*) xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2
write(*,*) '>>Hydrothermal effects on Diffusivity<<<'
if (xenhc2.lt.xenhc1) then
  write(*,*) 'WARNING:  xenhc2 corrected to be >/= xenhc1'
endif
write(*,*) ' xmaxdepth,     xmaxt,   xmaxstr,    xenhc1,    xenhc2'
write(*,'(5f10.2)') xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2
write(*,*) '>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<'
close( 9 )


return

2001 if_hydro = 0
write(*,*) '>>NO Hydrothermal effects on Diffusivity<<<'
return

end

!-------------------------------------------------------------
function HydroCond(j,i)
!!-------------------------------------------------------------
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'
!  common /hydroth/ xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
yc = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))

if( tmpr.le.xmaxt .and. yc.ge.xmaxdepth) then
cdum=dmin1(aps(j,i)/xmaxstr, 1.0)
HydroCond=(xenhc1 + cdum*(xenhc2-xenhc1))*conduct(iph)
else
HydroCond=conduct(iph)
endif

return


end function HydroCond


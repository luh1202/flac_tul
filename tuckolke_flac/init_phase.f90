!--------------------------------------------------------------
! Initialization of phases in each zone 
!--------------------------------------------------------------

subroutine init_phase
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


!character  phasedat*20

pi = 3.14159265358979323846
degrad = pi/180.




! Main phase
phasez = mphase

! Other phases in horizontal layers 
do k = 1,nphasl
    do i = 1,nx-1
        do j = ltop(k),lbottom(k)
            phasez(j,i) = lphase(k)
        end do
    end do
end do

!  Read distribution of the phases from the dat file
if (irphase .gt. 0) then
open(12,file='phasedat.inp')
read(12,*) nphasl
do 333 k=1,nphasl
read(12,*) lphase(k)
333 continue
do 332 i=1,nx-1
do 332 j=1,nz-1
!write(*,*) nx,nz
read(12,*) ii,jj,phasez(j,i)
!if(i.lt.201.and.j.le.30) phasez(j,i)=1.
!if(i.lt.120.and.j.lt.40) phasez(j,i)=2.
!if(i.ge.120.and.j.lt.40.and.j.gt.30) phasez(j,i)=1.
332  continue
!    call SysMsg('INIT_PHASE: Read phases from file: Option not ready yet!')
!    stop 21
close(12)
endif                        

!   Put different rheologies for inclusions 
do i = 1,inhom
    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        do j = ix1(i),ix2(i)
            do k = iy1(i),iy2(i)
                phasez(k,j) = inphase(i)
            end do
        end do
    endif

    ! Gauss shape:
    if (igeom(i).eq.1.or.igeom(i).eq.2) then
        ! symmetric case:
        if (igeom(i).eq.1) then
            ixc  = (ix1(i)+ix2(i))/2  
            iwidth = (ix2(i)-ix1(i))
        else
            ixc    = ix1(i)
            iwidth = ix2(i)-ix1(i)  
        endif
 
        iamp = iy2(i)-iy1(i)
  
        do j = ix1(i),ix2(i)
            itop = itop_geom(j,ixc,iwidth,iamp) 
            do k = iy1(i),iy2(i)
                if (k .ge. (iy2(i)-itop)) then 
                    phasez(k,j) = inphase(i)
                endif
            end do
        end do
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = int(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            phasez(k,j) = inphase(i)
        end do
    endif
    
    ! Weak zone in accumulated plastic strain at 45 degree        
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k = int(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            aps(k,j)=xinitaps
            phasez(k,j) = inphase(i)
        end do
    endif
end do

10 continue

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(i,j,phasez(j,i))
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do
    
! Put markers
do i = 1, nx-1
    do j = 1, nz-1
        rmarker(j,i) = i + (nx-1)*(j-1)
    end do
end do

return
end


!==========================================================
! Gauss perturbation at the top of heterogenity 
function itop_geom(j,ixc,iwidth,iamp) 
    
    itop_geom = iamp*exp(-(float(j-ixc)/(0.5*float(iwidth)))**2.)

return
end 


!==========================================================
! convert REAL to INTEGER phases for use as array indexes
!==========================================================
function iphase( ki,kj,phas )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

!iphase = nint( phas )
!if( iphase .le. 0 ) iphase = 1
!if( iphase .gt. nphase ) iphase = nphase

iphase = mphase
pdist_min = abs( phas - mphase )
!if(ki.eq.1.or.ki.eq.nx-1.or.kj.le.2.or.kj.ge.nz-2) then 
do iii = 1, nphasl
    pdist = abs( phas - lphase(iii) )
    if( pdist .lt. pdist_min ) then
        pdist_min = pdist
        iphase = lphase(iii)
    end if
end do
!goto 112
!endif
!phmax=0.
!phmin=100.
!iny=3
!if(kj.gt.43) iny = 1
!do jj = kj-iny,kj+iny
!   phmax = max(phasez(jj,ki),phmax)
!   phmin = min(phasez(jj,ki),phmin)
!   xdiff = abs(phmax-phmin)
!enddo
!if(abs(phmax-phmin).gt.1.) then
!   xdist1=abs(phas-phmin)/abs((phmax-phmin))
!   xdist2=abs(phas-phmax)/abs((phmax-phmin))
!   if(xdist1.le.xdist2) then
!   iphase = int(phmin)
!   else
!   iphase = int(phmax)
!   endif
!  if(iphase.eq.0) iphase = 1 
!   else
!iphase = mphase
!pdist_min = abs( phas - mphase )
!do iii = 1, nphasl
!    pdist = abs( phas - lphase(iii) )
!    if( pdist .lt. pdist_min ) then
!        pdist_min = pdist
!        iphase = lphase(iii)
!    end if
!end do
!endif
!if(iphase.eq.0) write(*,*) ki,kj,iphase
!112  continue
return
end


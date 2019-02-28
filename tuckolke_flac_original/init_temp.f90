!=========================================================
! Initiate temperature profile
!=========================================================
subroutine init_temp
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

allocatable :: a(:,:),b(:)

allocate( a(nz,nz),b(nz) )
!  Read distribution of temperatures from the dat file
if (mod(nloop,100).eq.0.and.nloop<1000) write(*,*) '>>>>>>>> INIT_TEMP <<<<<<<<<<', nloop  !G.Ito
if (igeotherm.eq.10) then
   if (nloop<2) write(*,*) '>> Initial Error Function Temperature Profile << '  !G.Ito
   uspread=0.0;
   do ii=1,nofbc
     if ((nofside(ii).eq.1).or.(nofside(ii).eq.3)) then
	uspread=dmax1(dabs(bca(ii)),uspread)
     endif
     if (nofside(ii).eq.4) then
       uspread=dmax1(dabs(bca(ii)),uspread)
     endif
   enddo	
  if (uspread.eq.0.0) then
    write(*,*) 'igeotherm=10 for error function temperature profile requires a '
    write(*,*) 'kinematically driven-spreading rate.  Stopping in INIT_TEMP'
    stop
  endif
endif
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
!     if(temp(j,i).ge.1000.) temp(j,i) = 1000.
    enddo
    enddo
    close(1)

    do i = 1, nx
     if(temp(16,i).gt.450.) then
         do 7 k= 1,16
           temp(k,i)= k*28.125
 7    continue
     endif
    enddo
 
   goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

! Temperature structure for ridges
! uses setup for viscosity from Alexei
!if (iynts.eq.1) then		G.Ito (commented out)
do i = 1,nx
    do j = 1,nz
       xc = cord(j,i,1) 
       yc = cord(j,i,2)
       yc0 = cord(1,i,2)
!       Line
       if (igeotherm .eq.0) then 
         geoth = g_y0c
!       Gauss perturbation
       elseif (igeotherm .eq.1 ) then
         geoth = g_y0c + g_amplitude*exp(-((xc-g_x0)/g_width)**2.)
!        Linear perturbation
       elseif (igeotherm .eq.2) then
         if ( abs(g_x0-xc).lt.1500) geoth = g_y0c+ g_amplitude
         if ( abs(g_x0-xc).lt.g_width) geoth = g_y0c+ &
            g_amplitude*(1.-abs(g_x0-xc)/g_width)
         if ( abs(g_x0-xc).ge.g_width) geoth = g_y0c
         
         !if(j.eq.1.and.i.eq.1) write(*,*) xc,yc0,geoth
       endif
! Temperatures
       te0 = tbos
       efold = efoldc
! E-fold depends on x (correction due to lateral change in geotherm)

       if((yc-yc0).ge.geoth) temp(j,i)=t_top+((te0-t_top)/geoth)*(yc-yc0)
       if((yc-yc0).lt.geoth) temp(j,i)=te0 + ((te0-t_top)/(0.05*geoth))*((yc-yc0)-geoth)
       
!----------------------------------------------------------------------------------------
! Error function temperature profile. G. Ito 7/04
!----------------------------------------------------------------------------------------
       if (igeotherm.eq.10) then 
	 diff = 1e-6;
	 diff = conduct(1)/cp(1)/den(1)
	 difft=dabs(xc)/uspread
	 difft=dmax1(difft,1e+11) !limit temperatures near ridge axis
         zz=dabs(yc)/(2*dsqrt(diff*difft))
	 temp(j,i)=t_bot-t_bot*erfc(zz)
!	 if (i.eq.nx/2.or.i.eq.2) write(*,'(I6,5e12.3)') j,yc,xc,difft,uspread,temp(j,i)
       endif
	 	
       if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
enddo
enddo
goto 10
!endif	G.Ito

! estimate initial temperature as linear (for first approx. of conductivities)
do j = 1,nz
    temp(j,1:nx) = (t_bot-t_top)/abs(rzbo)*abs(cord(j,1,2)-z0) + t_top
end do

dz = abs(cord(2,1,2)-cord(1,1,2))/1000

irep = 0
do while (.true.)
    a = 0; b = 0;
    a(1,1) = 1; b(1) = t_top;
    do j = 2, nz-1
        a(j,j-1) = cnd(j-1)+4*cnd(j)-cnd(j+1)
        a(j,j  ) = -8*cnd(j)
        a(j,j+1) = cnd(j+1)+4*cnd(j)-cnd(j-1)
        b(j) = -4*dz*dz*htgen(j)
    enddo
    a(nz,nz) = 1; b(nz) = t_bot;

    call Gauss(a,nz,b)

    tdiff = 0
    do j = 1,nz
        tdiff = max( abs( b(j)-temp(j,1) ), tdiff )
        temp(j,1:nx) = b(j)
    end do

    if( tdiff .lt. 0.1 ) exit

    irep = irep+1
    if( irep .gt. 1000 ) then
        call SysMsg('INIT_TEMP: No convergence !')
        stop
    endif

end do

deallocate( a,b )

10 continue

open( 1, file='temp0.dat' )
do j = 1,nz
    write(1,'(f5.1,1x,f6.1,1x,f6.1,1x,f6.1)') -cord (j,1,2)*1.e-3, temp(j,1)
end do
close(1)


! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

! Initial quadrilateral temperature perturbation
if( temp_per.gt.0 ) then
    temp(iy1t:iy2t,ix1t:ix2t) = temp(iy1t:iy2t,ix1t:ix2t) + temp_per
endif              


!call RedefineTemp

return
end


!=========================================================
function cnd( j )
include 'precision.inc'

cnd = Eff_conduct(1,j)

return
end


!=========================================================
function htgen( j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

y = - cord(j,1,2)*1.e-3

iph = iphase(1,j,phasez(j,1))

htgen = den(iph)*hs*exp(-y/hr) * 1.e+6

return
end


!=========================================================
subroutine RedefineTemp
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

write(*,*) 'ATTENTION! Special form of initial temperature distribution !'

return
end


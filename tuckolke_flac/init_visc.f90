
subroutine init_visc
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
irh = irheol(mphase) 
do i = 1,nx-1  
    do j = 1,nz-1 

        visn(j,i) = Eff_visc(i,j)

       if (irh.eq.11) then
       xc = 0.25*(cord (j,i  ,1) + cord(j+1,i  ,1) + &
                cord (j,i+1,1) + cord(j+1,i+1,1))

       yc = 0.25*(cord (j,i  ,2) + cord(j+1,i  ,2) + &
                cord (j,i+1,2) + cord(j+1,i+1,2))

       yc0 = 0.25*(cord (1,i  ,2) + cord(2,i  ,2) + &
                cord (1,i+1,2) + cord(2,i+1,2))


! Crust 
!       Line

       if (igeotherm .eq.0) geoth = g_y0c 

!       Gauss perturbation

       if (igeotherm .eq.1 ) then
       geoth = g_y0c + g_amplitude*exp(-((xc-g_x0)/g_width)**2.)
       endif


!       Linear perturbation

        if (igeotherm .eq.2) then
if ( abs(g_x0-xc).lt.g_width) geoth = g_y0c+ &
g_amplitude*(1.-abs(g_x0-xc)/g_width)- yc0
if ( abs(g_x0-xc).ge.g_width) geoth = g_y0c - yc0
        endif                 

!       geoth_c = geoth

! Viscosities
        vis0 = visc(mphase)
        efold = efoldc

! E-fold depends on x (correction due to lateral change in geotherm)

       if(ivis_shape.eq.0) visn(j,i)=vis0
       if(ivis_shape.eq.1) visn(j,i)=vis0+(yc-geoth)*efold
       if(ivis_shape.eq.2) visn(j,i)=vis0*exp((yc-geoth)/efold)

! min and max cutoffs

       if (visn(j,i).lt.v_min) visn(j,i)=v_min
       if (visn(j,i).gt.v_max) visn(j,i)=v_max               
       endif

!if ( abs(g_x0-xc).lt.g_width) write(*,*) v_max,v_min,visn(j,i)
    end do
end do
return
end

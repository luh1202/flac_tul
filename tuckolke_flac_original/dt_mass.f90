!----------------------------------------------------------------
! Inertial Masses and time steps - elastic and maxwell
!---------------------------------------------------------------

!    j,i           j,i+1
!     1-------------3
!     ! \  1     4/ !
!     !  1\     / 4 !
!     !     \ /     !
!     !     / \  2  !
!     !  3/    2\   !
!     ! / 3       \ !
!     2-------------4
!   j+1,i        j+1,i+1

subroutine dt_mass

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


c1d3 = 1./3.
c4d3 = 4./3.
c1d12 = 1./12.

if (nloop.le.1000) then
   strninert = strain_inert0
elseif (nloop.gt.1000.and.nloop.le.2000) then
   strninert = strain_inert0 + (strain_inert - strain_inert0)/(2001-nloop)
elseif (nloop.gt.2000.and.mloop.le.1000) then
   strninert = strain_inert0 + (strain_inert - strain_inert0)/(1001-mloop)
else
   strninert = strain_inert
end if


! minimal propagation distance
dlmin = dlmin_prop()

if (idt_scale .eq. 0) then
    ! find dt below
elseif (idt_scale.eq.1) then 
    ! choosing dt_elastic from sup.dat (non-automatic) 
    dt_elastic = dt_scale
elseif (idt_scale.eq.2) then
    ! choosing dt_elastic from tolerance conditions (automatic)
   if( vbc .eq. 0 ) then  !G.Ito: I dont think this option will happen because of "vbcal"
      ! let it be 1 cm/year
      dt_elastic = dlmin*frac*strninert/(1./sec_year/100)
!      if(nloop.le.1000) then
!         dt_elastic = dlmin*frac*strain_inert0/(1./sec_year/100)
!      else
!         dt_elastic = dlmin*frac*strain_inert/(1./sec_year/100)
!      end if
   else
      if(ny_inject.gt.0.and.ny_inject.le.2) then
         vbc = vbc - (rate_inject/2.)
      endif
      ! make sure dt_elastic does not go to infinity when vbc = rate inject
      ! added by M. Behn 3/01/03
      if(vbc.gt.0) then
         dt_elastic = dlmin*frac*strninert/vbc
!         if(nloop.le.1000) then
!            dt_elastic = dlmin*frac*strain_inert0/vbc
!         else
!            dt_elastic = dlmin*frac*strain_inert/vbc
!         end if
!         dt_elastic = dlmin*frac*strain_inert/(.05/sec_year/100)
!         write(*,*) 'dt_elastic (ext) = ', dt_elastic/sec_year
!         write(*,*) 'vbc = ', vbc
!         write(*,*) 'inj_count = ', inj_count
      else
         dt_elastic = dlmin*frac*strninert/(.05/sec_year/100)
!         if(nloop.le.1000) then
!            dt_elastic = dlmin*frac*strain_inert0/(.05/sec_year/100)
!         else
!            dt_elastic = dlmin*frac*strain_inert/(.05/sec_year/100)
!         end if
!         write(*,*) 'dt_elastic (inj) = ', dt_elastic
!         write(*,*) 'vbc = ', vbc
!         write(*,*) 'inj_count = ', inj_count
      endif
      if(ny_inject.gt.0.and.ny_inject.le.2) then
         vbc = vbc + (rate_inject/2.)
      endif
   endif
endif

amass = 0
dtmax_therm = 1.e+28
dt_maxwell = 1.e+28
!$DIR LOOP_PARALLEL
!$DIR PREFER_PARALLEL
!!$DIR LOOP_PRIVATE(j,iph,pwave,dens,vel_sound,rho_inert,am3,dte,diff,dtt,rmu,dt_maxwell,dt_m)
do 1 i = 1,nx-1
    do 1 j = 1,nz-1

        iph     = iphase(i,j,phasez(j,i))
        pwave   = rl(iph) + 0.6666*rm(iph)  
        dens    = den(iph)
        vel_sound = dlmin*frac/dt_elastic  
        rho_inert = pwave/(vel_sound*vel_sound)  

        ! Find the inert. density for given geometry, elas_mod and dt_scale
        ! idt_scale = 0 (dt = frac*dx_min * sqrt(dens/pwave) )
        ! idt_scale = 1 (dt is taken from sup.dat: dt = dt_scale)
        ! idt_scale = 2 (dt = frac*dx_min * tolerance/Vbc_max)
        if (idt_scale.gt.0) then 

            ! Distribution 1/3 of the inertial mass of each element to the nodes 
            ! am3=c1d12*area(ii,j,i)*pwave*(dlmax*dt_scale/frac)**2
            ! 1/12 = 1/3 * 1/2 (2 meshes) * 1/2 (1/area_num = 2 area_real) 
            am3=c1d12*rho_inert/area(1,j,i)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

            am3=c1d12*rho_inert/area(2,j,i)
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

            am3=c1d12*rho_inert/area(3,j,i)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3

            am3=c1d12*rho_inert/area(4,j,i)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

        else  ! idt_scale=0
            amass(j,i) = rmass(j,i)
            ! Find the dtime for given geometry, density and elas_mod
            dte = frac*dlmin*sqrt(dens/pwave)
            dt_elastic = min(dt_elastic,dte)
        endif 

        ! Find the maximum THERMAL time step from Stability Criterion
        ! dtmax = dxmin^2/diffusivity = dx^2/(lyamda/cp*dens)  
        diff = Eff_conduct(i,j)/den(iph)/Eff_cp(i,j)
        dtt = dlmin*dlmin/diff

        dtmax_therm =min (dtmax_therm,dtt)

        ! Calculate maxwell time step
        if (ivis_present .eq. 1) then
            !dt_m =visn(j,i)/rm(iph)*fracm
            if( (irheol(iph).eq.3 .OR. irheol(iph).eq.12) .AND. rm(iph).lt.1.e+11 ) then
                visc_cut = 1.e+19
                if( v_min .lt. visc_cut ) then
                    rmu = rm(iph) * v_min/visc_cut
                else
                    rmu = rm(iph)
                endif
                dt_m =v_min/rmu * fracm
                dt_maxwell = min (dt_m,dt_maxwell)
            endif
        endif

1 continue 

!write(*,*) 'dtmax_therm =', dtmax_therm, 'dt_maxwell =', dt_maxwell


return
end


!==============================================================
!   Global (in the whole domain) calculation of the shortest side - 
!   i.e. minimal propagation distance
!   dlmin = Area/Dmax for each triangle

function dlmin_prop()
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dlmin = 1.e+28

do 1 i = 1,nx-1
    do 1 j = 1,nz-1

        ! side 1-2 (triangles 1 and 3)
        dl = sqrt( (cord(j+1,i  ,1)-cord(j  ,i  ,1))**2 + (cord(j+1,i  ,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(1,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(3,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 2-4 (triangles 2 and 3)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j+1,i  ,1))**2 + (cord(j+1,i+1,2)-cord(j+1,i  ,2))**2 )
        dlm = 1./(area(2,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(3,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 4-3 (triangles 2 and 4)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j  ,i+1,1))**2 + (cord(j+1,i+1,2)-cord(j  ,i+1,2))**2 )
        dlm = 1./(area(2,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(4,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 3-1 (triangles 1 and 4)
        dl = sqrt( (cord(j  ,i+1,1)-cord(j  ,i  ,1))**2 + (cord(j  ,i+1,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(1,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(4,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! diagonal 1-4 (triangles 3 and 4)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j  ,i  ,1))**2 + (cord(j+1,i+1,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(3,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(4,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! diagonal 2-3 (triangles 1 and 2)
        dl = sqrt( (cord(j+1,i  ,1)-cord(j  ,i+1,1))**2 + (cord(j+1,i  ,2)-cord(j  ,i+1,2))**2 )
        dlm = 1./(area(1,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(2,j,i)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

1 continue

dlmin_prop = dlmin

return
end

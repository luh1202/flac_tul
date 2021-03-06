
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  
    
subroutine fl_rheol
include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
dimension njTinj(nz)
logical rh_sel
real*8 sarc11,sarc22

c1d3 = 1./3.
c1d6 = 1./6.

! initialization of random generators
jran = 13

if( mod(nloop,10).eq.0 .OR. ireset.eq.1 ) then
    rh_sel = .true.
else
    rh_sel = .false.
endif
rh_sel = .true.
irh=irheol(mphase)
if(irh.eq.11) call init_visc
if(iynts.eq.1) call init_temp 


! Initial stress boundary condition
! Accretional Stresses

inj_count = inj_count + 1
!if (ny_inject.eq.1) iinj = 1
!if (ny_inject.ge.2) iinj = (nx/2) + 1 

ninjtop=jinj1
ninjbot=jinj2
ninj=iinj2-iinj1+1
nelem_inject=jinj2-jinj1+1


irh_mark = 0
!!!DIR PREFER_PARALLEL
!$DIR LOOP_PARALLEL
!$DIR LOOP_PRIVATE(j,k,irh,iph,stherm,tmpr,bulkm,rmu,vis,coh,phi,psi,ipls,diss,hardn,sII_plas,sII_visc,quad_area,s0,s0a,s0b)
!$DIR LOOP_PRIVATE(de11,de22,de12,de33,dv,s11p,s22p,s12p,s33p,s11v,s22v,s12v,s33v,depl)
do 3 i = 1,nx-1
   rogh = 0. !Reset for pressure calculation [M. Behn April 4, 2006]
   ! Calculate ninjbot to be the minimum of ninjbot and the depth to the maximum cutoff temp (Tinj)
   jcnt=1
   njTinj(1) = ninjbot
    do j = 1,nz-1
        dcord = cord(1,i,2) - cord(j+1,i,2)
         if(dcord.gt.Tcinj) then
           njTinj(jcnt) = j
           jcnt = jcnt+1
         endif
    end do
    ninjbot = min(ninjbot,njTinj(1))

    do 3 j = 1,nz-1
        ! phasez (j,i) is number of a phase NOT a rheology 
        iph = iphase(i,j,phasez(j,i))
        irh = irheol(iph)

        ! Pressure - for calculating lithostatic stress [M. Behn April 4, 2006]
        if (ny_inject.eq.3) then
           tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
           densT = den(iph) * ( 1 - alfa(iph)*tmpr )
           dh1 = cord (j,i  ,2) - cord (j+1,i  ,2)
           dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
           dh  = 0.5 * (dh1+dh2)
           dPT = densT * g * dh
           dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )
           press = rogh + 0.5*dP
           rogh = rogh + dP
        endif
        ! End Pressure - for calculating lithostatic stress [M. Behn April 4, 2006]

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.*rm(iph)/3.
        rmu   = rm(iph)
        coh   = coha(iph)
        phi   = phimean(iph)
        psi   = psia(iph)
        poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
        young = rm(iph)*2.*(1.+poiss) 

        if (ny_inject.gt.0.and.ny_inject.le.2.and.inj_count.gt.ext_period) then      
	if (i.ge.iinj1.and.i.le.iinj2.and.j.le.ninjbot.and.j.ge.ninjtop) then
           dxinj=0.5*(cord(j,i+1,1)-cord(j,i,1)+cord(j+1,i+1,1)-cord(j+1,i,1))
	   fdum=ratfac*rate_inject*dt/(dxinj*dble(ninj))
           sarc1v = -young/(1.-poiss*poiss)*fdum                        ! Plane Stress 
           sarc2v = poiss*sarc1v
           sarc3v = 0.

           sarc1p = -(1.-poiss)*young/(1.+poiss)/(1.-2.*poiss)*fdum     ! Plane Strain Dikes - M. Behn 04/04/05
           sarc2p = sarc1p*poiss/(1.-poiss)
           sarc3p = poiss * (sarc1p + sarc2p)                                           ! Correct s33 - M. Behn 04/19/05     
        endif
        endif

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))
         ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,jran,hardn)
        ! Re-evaluate viscosity
        if (irh.eq.3 .or. irh.eq.12) then 
           if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(i,j)
        endif
        vis = visn(j,i)

        ! Cycle by triangles
        do 4 k = 1,4   

            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.
            dv = dvol(k,j,i)
            s11p(k) = stress0(1,k,j,i) + stherm 
            s22p(k) = stress0(2,k,j,i) + stherm 
            s12p(k) = stress0(3,k,j,i) 
            s33p(k) = stress0(4,k,j,i) + stherm 
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)
            if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                irheol_fl(j,i) = 0  
                stress0(1,k,j,i) = s11p(k)
                stress0(2,k,j,i) = s22p(k)
                stress0(3,k,j,i) = s12p(k)
                stress0(4,k,j,i) = s33p(k)

            elseif (irh.eq.3) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv)
                irheol_fl(j,i) = -1  
                stress0(1,k,j,i) = s11v(k)
                stress0(2,k,j,i) = s22v(k)
                stress0(3,k,j,i) = s12v(k)
                stress0(4,k,j,i) = s33v(k)

            elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                     s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12)
                irheol_fl(j,i) = 1
                stress0(1,k,j,i) = s11p(k)
                stress0(2,k,j,i) = s22p(k)
                stress0(3,k,j,i) = s12p(k)
                stress0(4,k,j,i) = s33p(k)

            elseif (irh.ge.11) then 
                ! Mixed rheology (Maxwell or plastic)
                if( rh_sel ) then
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                        s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12)
                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12)
                        stress0(1,k,j,i) = s11p(k)
                        stress0(2,k,j,i) = s22p(k)
                        stress0(3,k,j,i) = s12p(k)
                        stress0(4,k,j,i) = s33p(k)
                    else  ! irheol_fl(j,i) = -1
                        call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                            de11,de22,de33,de12,dv)
                        stress0(1,k,j,i) = s11v(k)
                        stress0(2,k,j,i) = s22v(k)
                        stress0(3,k,j,i) = s12v(k)
                        stress0(4,k,j,i) = s33v(k)
                    endif
                endif
             endif

4       continue

             if( irh.ge.11 .AND. rh_sel ) then
                ! deside - elasto-plastic or viscous deformation
                sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2
                sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

                if (sII_plas .lt. sII_visc) then
                   do k = 1, 4
                      stress0(1,k,j,i) = s11p(k)
                      stress0(2,k,j,i) = s22p(k)
                      stress0(3,k,j,i) = s12p(k)
                      stress0(4,k,j,i) = s33p(k)

                      if(ny_inject.gt.0.and.ny_inject.le.2.and.inj_count.gt.ext_period.and.j.le.ninjbot.and.j.ge.ninjtop) then

                         if(i.ge.iinj1.and.i.le.iinj2) then
                            !A0 = cord(j,i+1,1) - cord(j,i,1)
                            !dacr = rate_inject * dt
                            !if(Tinj.ge.0.and.k.eq.1) then
                              !temp(j,i) = temp(j,i) + ((Tinj - temp(j,i))*dacr) / A0
			      !temp(j,i) = temp(j,i) + Tinj*dacr/A0 
			      !temp(j,i) = (temp(j,i)*A0 + dacr*Tinj) / (A0 + dacr)
                            !endif
                            if(ntap.eq.0) then
                               dtpr=1.
                            elseif(ntap.eq.1) then 
                               dtpr = (j + 0.5 - ninjtop) * (2./nelem_inject)
                            elseif(ntap.eq.2) then 
                               dtpr = sqrt( (nelem_inject/2.)**2 - ((j-0.5) - ((ninjtop+ninjbot-1.)/2.))**2 ) / (nelem_inject/2.)
                            endif
                            stress0(1,k,j,i) = s11p(k) +sarc1p*dtpr
                            stress0(2,k,j,i) = s22p(k) +sarc2p*dtpr
                            stress0(4,k,j,i) = s33p(k) +sarc3p*dtpr   ! Reset s33 for accretion M. Behn 04/19/05
                         endif

                      elseif(ny_inject.eq.3.and.inj_count.gt.ext_period.and.j.le.ninjbot.and.j.ge.ninjtop) then

                         if(i.eq.iinj1) then
                            !A0 = cord(j,i+1,1) - cord(j,i,1)
                            !dacr = 2. * rate_inject * dt !Heat input doubled because only one column of elements
                            !if(Tinj.ge.0.and.k.eq.1) then
                               !temp(j,i) = temp(j,i) + ((Tinj - temp(j,i))*dacr) / A0
			       !temp(j,i) = temp(j,i) + Tinj*dacr/A0
			       !temp(j,i) = (temp(j,i)*A0 + dacr*Tinj) / (A0 + dacr)
                            !endif
                            ! Set stress to lithostatic in dike elements [M. Behn, April 4, 2006]
                            stress0 (1,k,j,i) = -press
                            stress0 (2,k,j,i) = -press
                            stress0 (3,k,j,i) = 0.
                            stress0 (4,k,j,i) = -press
                            ! Set stress to lithostatic in dike elements [M. Behn, April 4, 2006]
                         endif
                          
                       endif
                    end do
                   irheol_fl (j,i) = 1

                else 

                   do k = 1, 4
                      stress0(1,k,j,i) = s11v(k)
                      stress0(2,k,j,i) = s22v(k)
                      stress0(3,k,j,i) = s12v(k)
                      stress0(4,k,j,i) = s33v(k)

                      if(ny_inject.gt.0.and.ny_inject.le.2.and.inj_count.gt.ext_period.and.j.le.ninjbot.and.j.ge.ninjtop) then
                                               
                         if(i.ge.iinj1.and.i.le.iinj2) then
                            !A0 = cord(j,i+1,1) - cord(j,i,1)
                            !dacr = rate_inject * dt
                            !if(Tinj.ge.0.and.k.eq.1) then
                               !temp(j,i) = temp(j,i) + ((Tinj - temp(j,i))*dacr) / A0
			       !temp(j,i) = temp(j,i) + Tinj*dacr/A0 
			       !temp(j,i) = (temp(j,i)*A0 + dacr*Tinj) / (A0 + dacr)
                            !endif
                            if(ntap.eq.0) then
                               dtpr=1.
                            elseif(ntap.eq.1) then 
                               dtpr = (j + 0.5 - ninjtop) * (2./nelem_inject)
                            elseif(ntap.eq.2) then
                               dtpr = sqrt( (nelem_inject/2.)**2 - ((j-0.5) - ((ninjtop+ninjbot-1.)/2.))**2 ) / (nelem_inject/2.)
                            endif
                            stress0(1,k,j,i) = s11v(k) +sarc1v*dtpr
                            stress0(2,k,j,i) = s22v(k) +sarc2v*dtpr
                            stress0(4,k,j,i) = s33v(k) +sarc3v*dtpr ! Reset s33 for accretion M. Behn 04/19/05
                          endif

                      elseif(ny_inject.eq.3.and.inj_count.gt.ext_period.and.j.le.ninjbot.and.j.ge.ninjtop) then

                         if(i.eq.iinj1) then
                            !A0 = cord(j,i+1,1) - cord(j,i,1)
                            !dacr = 2. * rate_inject * dt !Heat input doubled because only one column of elements
                             !if(Tinj.ge.0.and.k.eq.1) then
                                !temp(j,i) = temp(j,i) + ((Tinj - temp(j,i))*dacr) / A0
			       !temp(j,i) = temp(j,i) + Tinj*dacr/A0
			       !temp(j,i) = (temp(j,i)*A0 + dacr*Tinj) / (A0 + dacr)
                             !endif
                            ! Set stress to lithostatic in dike elements [M. Behn, April 4, 2006]
                            stress0 (1,k,j,i) = -press
                            stress0 (2,k,j,i) = -press
                            stress0 (3,k,j,i) = 0.
                            stress0 (4,k,j,i) = -press
                            ! Set stress to lithostatic in dike elements [M. Behn, April 4, 2006]
                          endif

                      endif
                   end do
                   irheol_fl (j,i) = -1
                endif
             endif

        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then
        
            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1./(area(1,j,i)+area(2,j,i))
            s0a=0.5*(stress0(1,1,j,i)+stress0(2,1,j,i))
            s0b=0.5*(stress0(1,2,j,i)+stress0(2,2,j,i))
            s0=(s0a*area(2,j,i)+s0b*area(1,j,i))*quad_area
            stress0(1,1,j,i) = stress0(1,1,j,i) - s0a + s0
            stress0(2,1,j,i) = stress0(2,1,j,i) - s0a + s0
            stress0(1,2,j,i) = stress0(1,2,j,i) - s0b + s0
            stress0(2,2,j,i) = stress0(2,2,j,i) - s0b + s0

            ! For C and D couple:
            quad_area = 1./(area(3,j,i)+area(4,j,i))
            s0a=0.5*(stress0(1,3,j,i)+stress0(2,3,j,i))
            s0b=0.5*(stress0(1,4,j,i)+stress0(2,4,j,i))
            s0=(s0a*area(4,j,i)+s0b*area(3,j,i))*quad_area
            stress0(1,3,j,i) = stress0(1,3,j,i) - s0a + s0
            stress0(2,3,j,i) = stress0(2,3,j,i) - s0a + s0
            stress0(1,4,j,i) = stress0(1,4,j,i) - s0b + s0
            stress0(2,4,j,i) = stress0(2,4,j,i) - s0b + s0

        endif

        !  ACCUMULATED PLASTIC STRAIN
        ! Average the strain for pair of the triangles
        ! Note that area (n,it) is inverse of double area !!!!!
        if (sII_plas .lt. sII_visc) then
           aps(j,i) = aps(j,i) &
                + 0.5*( depl(1)*area(2,j,i)+depl(2)*area(1,j,i) ) / (area(1,j,i)+area(2,j,i)) &
                + 0.5*( depl(3)*area(4,j,i)+depl(4)*area(3,j,i) ) / (area(3,j,i)+area(4,j,i))
           if( aps(j,i) .lt. 0. ) aps(j,i) = 0.
        endif
        ! LINEAR HEALING OF THE PLASTIC STRAIN
        if (tau_heal .ne. 0.) &
            aps (j,i) = aps (j,i)/(1.+dt/tau_heal) 
        ! TOTAL FINITE STRAIN  
        strain (1,j,i) = strain(1,j,i) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain (2,j,i) = strain(2,j,i) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain (3,j,i) = strain(3,j,i) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
3 continue

if (inj_count.eq.inj_period+ext_period) inj_count = 0

return
end


!------ Elasto-Plastic

subroutine plastic(bulkm,rmu,coh,phi,psi,depls,ipls,diss,hardn,s11,s22,s33,s12,de11,de22,de33,de12)
      
include 'precision.inc'
include 'params.inc'
       
logical fsflg

pi = 3.14159265358979323846
degrad = pi/180.
c1d3 = 1./3.
c4d3 = 4./3.
c2d3 = 2./3.

! press_add formaely was passed by a parameter. in my case it is always zero.
press_add = 0.

! ------------------------------
! Initialization section
! ------------------------------
depls = 0. 
diss = 0. 
ipls = 0 
      
sphi  = dsin(phi * degrad)
spsi  = dsin(psi * degrad)
anphi = (1.+ sphi) / (1.- sphi)
anpsi = (1.+ spsi) / (1.- spsi)
amc   = 2.0 * coh * sqrt (anphi)
e1    = bulkm + c4d3 * rmu
e2    = bulkm - c2d3 * rmu
sh2   = 2.0 * rmu
x1    = (e1 - e2*anpsi + e1*anphi*anpsi - e2*anphi)

if (phi.eq. 0.) then
    ten_max=ten_off
else
    ten_max=min(ten_off,coh/(tan(phi*degrad)))
end if

! ---------------
! Running section
! ---------------

!---- get new trial stresses from old, assuming elastic increment
!---- add press (which is positive press = - (sxx+syy)*0.5, 
!---- which has 2 components: add pressure due to application of forces from the top 
!---- and subtract pressure of the fluid
s11i = s11 + (de22 + de33) *e2  + de11 *e1 - press_add
s22i = s22 + (de11 + de33) *e2  + de22 *e1 - press_add
s12i = s12 + de12 * sh2
s33i = s33 + (de11 + de22) *e2  + de33 *e1 - press_add
sdif = s11i - s22i
s0   = 0.5 * (s11i + s22i)
rad  = 0.5 * sqrt(sdif*sdif + 4.0 *s12i*s12i)
! principal stresses
si  = s0 - rad
sii = s0 + rad
psdif = si - sii
if (irh_mark.eq.1) then 
    s11 = s11i + press_add
    s22 = s22i + press_add
    s33 = s33i + press_add
    s12 = s12i
 return 
endif

!--------------------------------------------------------- 
!                         3D version 
!--------------------------------------------------------- 
if (ndim.eq.3) then
    !-- determine case ---
    if (s33i .gt. sii) then
        !- s33 is minor p.s. --
        icase = 3
        s1 = si
        s2 = sii
        s3 = s33i
    elseif (s33i .lt. si) then
        !- s33 is major p.s. --
        icase = 2
        s1 = s33i
        s2 = si
        s3 = sii
    else
        !- s33 is intermediate --
        icase = 1
        s1 = si
        s2 = s33i
        s3 = sii
    endif
endif

!------------------------------------------------------- 
!         2D version 
!------------------------------------------------------- 
if (ndim.eq.2) then 
    icase = 1
    s1 = si 
    s2 = s33i
    s3 = sii
endif 

!--------------------------------------------------------
! Check for tensional failure before the shear failure
!-------------------------------------------------------
 
!----- general tension failure
if (s1 .ge. ten_max) then
    ipls = -5
    goto 800
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and. ndim .eq.3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
endif

!- partial failure (only if s3 is greater than ten_max) 
if (s3 .ge. ten_max) then
    s3 = ten_max 
    ipls = -7
endif

!- check for shear yield (if fs<0 -> plastic flow)
fs = s1 - s3 * anphi + amc
fsflg = fs .lt. 0.0
if (fsflg) then
    !-- yielding in shear ----
    if (icase .eq. 1) ipls = -2
    if (icase .eq. 2) ipls = -3
    if (icase .eq. 3) ipls = -4
    alams = fs/(x1+hardn)
    s1 = s1 - alams * (e1 - e2 * anpsi )
    s2 = s2 - alams * e2 * (1.0 - anpsi )
    s3 = s3 - alams * (e2 - e1 * anpsi )

    ! Increment of the plastic strain (2nd Invariant)
    dep1 = alams
    dep3 = -alams*anpsi

    ! FOR 2D caculations
    depm = 0.5*(dep1+dep3)
    depls = 0.5*abs(dep1-dep3)

    ! Dissipation rate
    diss = s1*dep1+s3*dep3
else
    !-- no failure at all (elastic behaviour)
    s11 = s11i + press_add
    s22 = s22i + press_add
    s33 = s33i + press_add
    s12 = s12i
    return
endif

!- general tension failure?
if (s1 .ge. ten_max) then
    ipls = -5
    goto 800
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and.ndim.eq.3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
    goto 205
endif

!- uniaxial tension ... minor p.s. ---
if (s3 .ge. ten_max) then
    ipls = -7
    s3 = ten_max 
endif

!- direction cosines
205 continue
if ( psdif .eq. 0. ) then
    cs2 = 1.
    si2 = 0.
else
    cs2 = sdif / psdif
    si2 = 2.0 * s12i / psdif
endif

!- resolve back to global axes
goto (210,220,230), icase

210 continue
dc2 = (s1-s3) * cs2
dss = s1 + s3
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s1 - s3) * si2
s33 = s2
goto 240

220 continue
dc2 = (s2-s3) * cs2
dss = s2 + s3
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s2 - s3) * si2
s33 = s1
goto 240

230 continue
dc2 = (s1-s2) * cs2
dss = s1 + s2
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s1 - s2) * si2
s33 = s3

240 continue

s11 = s11 + press_add
s22 = s22 + press_add
s33 = s33 + press_add

return

!-- set stresses to plastic apex ---
800   continue
s11        = ten_max
s22        = ten_max
s12        = 0.0
s33        = ten_max

s11 = s11 + press_add
s22 = s22 + press_add
s33 = s33 + press_add

return
end



!==================================================================
! Prepare plastic properties depending on softening or randomness

subroutine pre_plast (i,j,coh,phi,psi,jran,hardn)
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
character*100 msg

data im /259200/, ia /7141/, ic /54773/
iph = iphase(i,j,phasez(j,i))
pls_curr = aps(j,i)
! Hardening modulus
hardn = 0.

! Piece-wise linear softening
isoft = 0

! Find current properties from linear interpolation
!if (i.le.10.or.i.ge.290) then
!     plstrain(2) = 0.5
!endif
!   cohesion(1)=34.e+6
!if (i.eq.iinj.or.i.eq.iinj-1) then
!   cohesion(1) = 24.e+6
!   cohesion(2) = 24.e+6
!   cohesion(3) = 24.e+6
!   cohesion(4) = 24.e+6
!endif
!endif
do is = 1,nsegments-1
    pl1 = plstrain(is)
    pl2 = plstrain(is+1)
 
    if (pls_curr .ge. pl1 .and. pls_curr .le. pl2) then
        isoft = 1

        ! Friction
        tgf = (fric(is+1)-fric(is))/(pl2-pl1)
        phi =  fric(is) + tgf*(pls_curr-pl1)

        ! Cohesion and Dilatation
        tgd = (dilat(is+1)-dilat(is))/(pl2-pl1)
        psi =  dilat(is) + tgd*(pls_curr-pl1)
 
        tgc = (cohesion(is+1)-cohesion(is))/(pl2-pl1)
        coh =  cohesion(is) + tgc*(pls_curr-pl1)

        ! Hardening Modulus (for COhesion ONLY)
        hardn = (cohesion(is+1)-cohesion(is))/(pl2-pl1) 
    endif
end do
 
if (isoft.eq.0 ) then
    write( msg, * ) 'Pre_plast: No segment for current plastic strain (j,i,aps): ', iph,j,i,pls_curr
    call SysMsg( msg )
    stop 25
endif

 
! Random plastic properties
if (phidisp(iph).gt.0.0001) then
    jran = mod (jran*ia+ic, im)
    ran  = float (jran)/float(im)
    randphi = phidisp(iph)*( ran -0.5)
    phi = phimean(iph) + randphi
endif

if (cohdisp(iph).gt.0.0001) then
    jran = mod (jran*ia+ic, im)
    ran  = float (jran)/float(im)
    randcoh = cohdisp(iph)*( ran -0.5)
    coh = coha(iph) + randcoh
endif

return
end


!==============================================
! Density
function Eff_dens( i, j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


iph = iphase(i,j,phasez(j,i))
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
press = 0
do ii = 1, 4
    press = press - (stress0(1,ii,j,i)+stress0(2,ii,j,i)+stress0(4,ii,j,i))/3
enddo
press = press / 4
dens = den(iph) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )

! Effect of melt
fmelt(j,i) = Eff_melt( i,j )
dens = dens * ( 1.-0.1*fmelt(j,i) )
Eff_dens = dens

return
end


!==============================================
! Melt fraction
!==============================================
function Eff_melt( i, j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


iph = iphase(i,j,phasez(j,i))
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

! Effect of melting on density (here - max 10%)
if( tmpr .lt. ts(iph) ) then
    fm = 0.
elseif( tmpr .lt. tk(iph) ) then
    fm = fk(iph)/(tk(iph)-ts(iph)) * (tmpr-ts(iph))
elseif( tmpr .lt. tl(iph) ) then
    fm = (1.-fk(iph))/(tl(iph)-tk(iph))*(tmpr-tk(iph)) + fk(iph)
else
    fm = 1.
endif

if( fm .lt. 0 ) fm = 0.
if( fm .gt. 1 ) fm = 1.

Eff_melt = fm

return
end


!=================================================
! Effective Heat Capacity incorporating latent heat]
! G. Ito 8/2006, Modified to use Tliq, Tsol, and cp(1) 
!=================================================
function Eff_cp( i,j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


!iph = iphase(i,j,phasez(j,i))
!Eff_cp = cp(iph)


tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
  yc = 0.5*(cord(1,i,2)+cord(1,i+1,2)) -  &
       0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  

if(tmpr.ge.Tsol.and.tmpr.le.Tliq.and.yc.le.Tcinj) then			!G.Ito 8/2/06
    Eff_cp = cp(1) + xlatheat/(Tliq-Tsol)			!G.Ito 8/2/06
else								!G.Ito 8/2/06
    Eff_cp = cp(1)						!G.Ito 8/2/06
endif	

!if( tmpr .lt. ts(iph) ) then
!    Eff_cp = cp(iph)
!elseif( tmpr .lt. tk(iph) ) then
!    Eff_cp = cp(iph) + xlatheat * fk(iph)/(tk(iph)-ts(iph))
!elseif( tmpr .lt. tl(iph) ) then				!G.Ito 8/2/06
!    Eff_cp = cp(iph) + xlatheat * (1.-fk(iph))/(tl(iph)-tk(iph))
!else
!    Eff_cp = cp(iph)
!endif


! HOOK
! Intrusions - melting effect - see user_ab.f90
!if( if_intrus .eq. 1 ) then
!
!    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!    if( tmpr .lt. ts(iph) ) then
!        Eff_cp = cp(iph)
!    elseif( tmpr .lt. tl(iph)+1 ) then
!        Eff_cp = cp(iph) + xlatheat/(tl(iph)-ts(iph))
!    else
!        Eff_cp = cp(iph)
!    endif
!endif

return
end


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( i,j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


iph = iphase(i,j,phasez(j,i))
Eff_conduct = conduct(iph)

!if( den(iph) .lt. 3000. ) then  ! for crustal material
!    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!    if( tmpr.lt.25 ) tmpr = 25.
!    Eff_conduct = -0.38*dlog(tmpr) + 4.06
!endif

! HOOK
! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
if( if_hydro .eq. 1 ) then
    Eff_conduct = HydroCond(i,j)				!G.Ito
!   if (aps(j,i).gt.xmaxstr) write(*,*) Eff_conduct, aps(j,i)	
endif

return
end



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( i,j )
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

r=8.31448e0

iph = iphase(i,j,phasez(j,i))
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

pow  =  1./pln(iph) - 1.
pow1 = -1./pln(iph) 

vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
        exp(eactiv(iph)/(pln(iph)*r*(tmpr+273.)))*1.e+6

! Effect of melt
fmelt_crit = 0.05
if( fmelt(j,i) .gt. 0. ) then
    if( fmelt(j,i) .lt. fmelt_crit ) then
        vislog = fmelt(j,i)/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
        vis = 10.**vislog
    else
        vis = v_min
    endif
endif


! limiting from above (quasi-Peierls)
!sIImax = 5.e+8
!vis_peierls = sIImax / srat / 2
!if( vis .gt. vis_peierls ) vis = vis_peierls


! Final cut-off
if (vis .lt. v_min) vis = v_min
if (vis .gt. v_max) vis = v_max

Eff_visc = vis

return
end

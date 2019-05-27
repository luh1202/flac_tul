!==============================================
! Density
function Eff_dens( j, i)
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  include 'phases.inc'

  zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))


  press = 0
  do ii = 1, 4
      press = press - (stress0(j,i,1,ii)+stress0(j,i,2,ii)+stress0(j,i,4,ii))/3
  enddo
  press = press / 4

  if (iint_marker.ne.1) then
      iph = iphase(j,i)
      dens = den(iph) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )

  else
      Eff_dens = 0.
      do k = 1, nphase
          ratio = phase_ratio(k,j,i)
          ! when ratio is small, it won't affect the density
          if(ratio .lt. 0.01) cycle

          dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )


          press = dens*g*zcord

          dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

          Eff_dens = Eff_dens + ratio*dens

      enddo
  endif
  return
end function Eff_dens

!=================================================
! Effective Heat Capacity incorporating latent heat
!=================================================
function Eff_cp( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
yc = 0.5*(cord(1,i,2)+cord(1,i+1,2)) -  &
0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))


if(tmpr.ge.Tsol.and.tmpr.le.Tliq.and.yc.le.Tcinj) then            !G.Ito 8/2/06
Eff_cp = cp(iph) + xlatheat/(Tliq-Tsol)            !G.Ito 8/2/06
else                                !G.Ito 8/2/06
Eff_cp = cp(iph)                        !G.Ito 8/2/06
endif

iph = iphase(j,i)
Eff_cp = cp(iph)

return
end function Eff_cp


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'
!if_hydro = 1
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
yc = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))

if (iint_marker.ne.1) then
iph = iphase(j,i)
cond = conduct(iph)

! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
!if( if_hydro .eq. 1 ) then
if( tmpr.le.xmaxt .and. yc.ge.xmaxdepth) then
cdum=dmin1(aps(j,i)/xmaxstr, 1.0)
HydroCond=xenhc1*cond
else
HydroCond=cond
endif
cond = HydroCond
Eff_conduct = cond

else
Eff_conduct = 0.
do k = 1 , nphase

! when ratio is small, it won't affect the density
if(phase_ratio(k,j,i) .lt. 0.01) cycle

cond = conduct(k)

! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
!if( if_hydro .eq. 1 ) then
!print *,'conduct(k) =',conduct(k)
!print *,'xmaxt =',xmaxt
!print *,'yc =',yc
!print *,'xmaxdepth =',xmaxdepth
if( tmpr.le.xmaxt .and. yc.ge.xmaxdepth) then
!print *, 'hi'
!cdum=dmin1(aps(j,i)/xmaxstr, 1.0)
HydroCond=xenhc1*cond
!print *,'Hydrocond1 =',HydroCond
else
HydroCond=cond
!print *,'Hydrocond2 =',HydroCond
endif
cond = HydroCond
Eff_conduct = phase_ratio(k,j,i)*cond
!print*, 'Eff_conduct', Eff_conduct
enddo
endif
!write(*,*) Eff_conduct, cond
return
end function Eff_conduct



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

Eff_visc = 0.
r=8.31448e0
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))

srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

if (iint_marker.ne.1) then
    iph = iphase(j,i)

    pow  =  1./pln(iph) - 1.
    pow1 = -1./pln(iph)

    vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
         exp(eactiv(iph)/(pln(iph)*r*(tmpr+273.)))*1.e+6



    ! Final cut-off
    if (vis .lt. v_min) vis = v_min
    if (vis .gt. v_max) vis = v_max

    Eff_visc = vis

else

    do k = 1, nphase
        if(phase_ratio(k,j,i) .lt. 0.01) cycle

        pow  =  1./pln(k) - 1.
        pow1 = -1./pln(k)

        vis = 0.25 * srat**pow*(0.75*acoef(k))**pow1* &
             exp(eactiv(k)/(pln(k)*r*(tmpr+273.)))*1.e+6

        ! Final cut-off
        if (vis .lt. v_min) vis = v_min
        if (vis .gt. v_max) vis = v_max

        ! harmonic mean
        Eff_visc = Eff_visc + phase_ratio(k,j,i) / vis
        !write(*,*) i,j, Eff_visc, vis, tmpr,phase_ratio(k,j,i)
    enddo

    Eff_visc = 1 / Eff_visc
endif

return
end function Eff_visc

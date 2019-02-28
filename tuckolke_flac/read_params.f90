! Reads problem parameters from input file

subroutine read_params
include 'precision.inc'
include 'params.inc'


iu =4 
open( iu, file='drezina.inp' )

! write platform file
open (1,file='platform.s')
write(1, '(a)') 'ieee-le'
close(1)

! MESH
call AdvanceToNextInputLine( 4 )
read(4,*) nx, nz
open(11,file='nxnz.0')
write(11,*) nx, nz
close(11)
nq = nx*nz
nx = nx+1
nz = nz+1
call AdvanceToNextInputLine( 4 )
read(4,*) x0, z0
call AdvanceToNextInputLine( 4 )
read(4,*) rxbo, rzbo
call AdvanceToNextInputLine( 4 )
read(4,*) nzonx 
if (nzonx .eq.0) then
    nzonx = 1
    nelz_x(1) = nx - 1
    sizez_x(1) = 1.
    go to 166
endif
do i = 1, nzonx
    call AdvanceToNextInputLine( 4 )
    read(4,*) nelz_x(i), sizez_x(i)
end do
166 continue

call AdvanceToNextInputLine( 4 )
read(4,*) nzony
if (nzony .eq.0) then
    nzony = 1
    nelz_y(1)    = nz - 1
    sizez_y(1)   = 1.
    go to 177
endif
do i = 1,nzony
    sizez_y(i) = 1.
    call AdvanceToNextInputLine( 4 )
    read(4,*) nelz_y(i), sizez_y(i)
end do
177 continue
call AdvanceToNextInputLine( 4 )
read(4,*) npelem			!Surface particles

! MECHANICAL CONDITIONS
call AdvanceToNextInputLine( 4 )
read(4,*) ynstressbc,ydrsides
call AdvanceToNextInputLine( 4 )
read(4,*) nofbc
call AdvanceToNextInputLine( 4 )
        do 21 i = 1,nofbc
        read(iu,*) nofside(i),nbc1(i),nbc2(i),nbc(i),  &
        bca(i),bcb(i),bcc(i),  &
        bcd(i),bce(i),bcf(i),bcg(i),bch(i),bci(i)
 21     continue
call AdvanceToNextInputLine( 4 )
! hydrostatic pressure applied at the bottom
call AdvanceToNextInputLine( 4 )
read(4,*) nyhydro,pisos,iphsub,drosub,damp_vis
! gravity
call AdvanceToNextInputLine( 4 )
read(4,*) g



! THERMAL CONDITIONS
call AdvanceToNextInputLine( 4 )
read(4,*) itherm 
call AdvanceToNextInputLine( 4 )
read(4,*) istress_therm
call AdvanceToNextInputLine( 4 )             ! thermal stresses
read (4,*) ishearh                           ! shear heating
call AdvanceToNextInputLine( 4 )
read (4,*) t_top
call AdvanceToNextInputLine( 4 )
read (4,*) t_bot  
call AdvanceToNextInputLine( 4 )
read (4,*) hs, hr 
! boundary conditions at the bottom (1-T,2-Flux) 
call AdvanceToNextInputLine( 4 )
read (4,*) itemp_bc, bot_bc
if( itemp_bc.eq.2 ) bot_bc = bot_bc/1000  ! convert in W/m3
! temperature pertrubation (rectangular)
call AdvanceToNextInputLine( 4 )
read (4,*) temp_per, ix1t, ix2t, iy1t, iy2t
! Predefined distributions
call AdvanceToNextInputLine( 4 )
read(4,*) irtemp
if ( irtemp .gt. 0 ) then
    call AdvanceToNextInputLine( 4 )
    read(4,*) tempfile
else
    call AdvanceToNextInputLine( 4 )
    read(4,*)
endif
! time scale
call AdvanceToNextInputLine( 4 )
read (4,*) time_scale
! temp structure
call AdvanceToNextInputLine( 4 )
read (4,*) iynts,tbos
call AdvanceToNextInputLine( 4 )
read (4,*) iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1 


! RHEOLOGY
call AdvanceToNextInputLine( 4 )
read(4,*) nphase
do i = 1, nphase 
    call AdvanceToNextInputLine( 4 )
    read(4,*) irheol(i),visc(i),den(i),alfa(i),beta(i),pln(i),acoef(i),eactiv(i),rl(i),rm(i), &
    coha(i),cohdisp(i),phimean(i),phidisp(i),psia(i),conduct(i),cp(i),ts(i),tl(i),tk(i),fk(i)
end do
! Flag to take initial phase distribution from a file
call AdvanceToNextInputLine( 4 )
read(4,*) irphase
write(*,*) irphase

if ( irphase .gt. 0 ) then
    call AdvanceToNextInputLine( 4 )
    read(4,*) phasefile
else
    call AdvanceToNextInputLine( 4 )
    read(4,*)
endif
! main phase
call AdvanceToNextInputLine( 4 )
read(4,*) mphase
! number of horizontal layers
call AdvanceToNextInputLine( 4 )
read(4,*) nphasl
! layers
do i = 1, nphasl
    call AdvanceToNextInputLine( 4 )
    read(4,*) ltop(i), lbottom(i), lphase(i)
end do
! inclusions
call AdvanceToNextInputLine( 4 )
read(4,*) inhom
if( inhom .gt. 9 ) then
    call SysMsg('Read_params: Increase arrays for inhomogenities')
    stop 26
endif
do i = 1, inhom
    call AdvanceToNextInputLine( 4 )
    read(4,*) ix1(i), ix2(i), iy1(i), iy2(i), inphase(i), igeom(i), xinitaps
end do
! Tension cut-off
call AdvanceToNextInputLine( 4 )
read(4,*) ten_off
! softening for plasticity
call AdvanceToNextInputLine( 4 )
read(4,*) nysoft
call AdvanceToNextInputLine( 4 )
read(4,*) nsegments
if( nsegments .gt. 5 ) then
    call SysMsg('Read_params: Increase arrays for softening model')
    stop 27
endif
do i = 1, nsegments
    call AdvanceToNextInputLine( 4 )
    read(4,*) plstrain(i), fric(i), dilat(i), cohesion(i)
end do
!linear healing parameter
call AdvanceToNextInputLine( 4 )
read(4,*) tau_heal
! viscosity limits
call AdvanceToNextInputLine( 4 )
read(4,*) v_min, v_max, ivis_shape,efoldc
call AdvanceToNextInputLine( 4 )
read(4,*)igeotherm,g_x0,g_y0c, g_amplitude,g_width 
call AdvanceToNextInputLine( 4 )
read(4,*) ny_inject,ntap,iinj1,iinj2,jinj1,jinj2,xlatheat,Tcinj  !G.Ito 8/06
call AdvanceToNextInputLine( 4 )
read(4,*) Tsol,Tliq, rate_inject, inj_period, ext_period, inj_count !G.Ito 8/06
!call AdvanceToNextInputLine( 4 )
!read(4,*) ny_inject, ntap, ntop_inject, nelem_inject, Tinj, Tcinj, rate_inject, inj_period, ext_period, inj_count
if (ny_inject.eq.1) then
  ratfac=1.0
else 
  ratfac=2.0
endif
! REMESHING
call AdvanceToNextInputLine( 4 )
read(4,*)  ny_rem, mode_rem, ntest_rem, angle_rem
if ( mode_rem.ne.1 .and. mode_rem.ne.11 .and. mode_rem.ne.3 ) then
    call SysMsg('Illegal remeshing mode! Allowable - 1, 3 or 11')
    stop
endif
! dx_rem - remeshing criteria for mode_rem=11
call AdvanceToNextInputLine( 4 )
read(4,*)  dx_rem
! diffusion of topography
call AdvanceToNextInputLine( 4 )
read(4,*)  topo_kappa, bottom_kappa


! PROCESS CONTROL
! inertial mass scaling
call AdvanceToNextInputLine( 4 )
read(4,*)  idt_scale
call AdvanceToNextInputLine( 4 )
read(4,*)  dt_scale, strain_inert, strain_inert0
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_rmasses
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_imasses
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_visc
call AdvanceToNextInputLine( 4 )
read(4,*)  ifreq_avgsr
! acceleration parameters
call AdvanceToNextInputLine( 4 )
read(4,*) amul, ratl, ratu
call AdvanceToNextInputLine( 4 )
read(4,*) frac, fracm
call AdvanceToNextInputLine( 4 )
read(4,*) n_boff_cutoff
call AdvanceToNextInputLine( 4 )
read(4,*) movegrid,ndim
call AdvanceToNextInputLine( 4 )
read(4,*) demf, mix_strain, mix_stress

! OUTPUT PARAMETERS
call AdvanceToNextInputLine( 4 )
read(4,*) time_max
time_max = time_max * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read (4,*) dtout_screen
dtout_screen = dtout_screen * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read(4,*) dtout_file
dtout_file = dtout_file * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
!read(4,*) nprofil
!do 135 i=1,nprofil
!read(4,*) hv_out(i:i), iprof_out(i) 
!135 continue
!call AdvanceToNextInputLine( 4 )
read(4,*) io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz,io_sxz,io_pres, &
    io_temp,io_melt,io_visc,io_phas,io_mark,io_src,io_cond,io_diss,io_forc, &
    io_hfl,io_topo
call AdvanceToNextInputLine( 4 )
read(4,*) lastout
call AdvanceToNextInputLine( 4 )
read(4,*) dtsave_file
dtsave_file = dtsave_file * 1000 * sec_year
call AdvanceToNextInputLine( 4 )
read(4,*) lastsave

close (iu)


! ADDITIONAL PARTICULAR INPUT
call ReadMoreParams()

call Verify_Input

return
end


!========================================================
subroutine AdvanceToNextInputLine( iu )
character*1 buf

10    read(iu, '(A1)') buf
      if( buf(1:1).eq.';' ) then
       goto 10 
       else
        backspace( iu )
       return
       endif

20 continue
print *, 'AdvanceToNextInputLine: EOF reached!'
stop

30 continue
print *, 'AdvanceToNextInputLine: Error reading file!'
stop

return
end

!========================================================
subroutine Verify_Input
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

if (nx.gt.mnx+1.or.nz.gt.mnz+1) then
  write (*,*) 'Increase maximum array size in "arrays.inc"!!!'
  write (*,*) 'nx, nz, mnx+1, mnz+1 =', nx, nz, mnx+1, mnz+1
  stop
endif

if (ny_inject.gt.0) then
  write(*,*) '------------------------------------------------------------------------'
  write(*,'(a26,4i5)') 'Dike injects in elems iinj1,iinj2,jinj1,jinj2:', &
                            iinj1,iinj2,jinj1,jinj2
  write(*,'(a28,4ES12.4)') 'Tcinj, Tsol, Tliq, xlatheat=', Tcinj, Tsol, Tliq, xlatheat
  if (Tsol.gt.t_bot) then
    write(*,*) 'Tsol>t_bot so NO LATENT HEAT added'
    if (xlatheat.ne.0) then
      write(*,*) 'Verify_Input: xlatheat not zero!! xlatheat=', xlatheat
      stop
    endif
  endif
  if (xlatheat.eq.0.and.Tsol.lt.t_bot) then
    write(*,*) 'Verify_Input: xlatheat=0 but Tsol < t_bot!!'
    stop
  endif
  fdum=0.
  do 30 i=1,nofbc
    if ((nofside(i).eq.3.or.nofside(i).eq.1).and.nbc(i).eq.10) then
        fdum=fdum+abs(bca(i))
    endif
30 continue

  write(*,'(a26,f5.3)') 'Input rate_inject gives M=', ratfac*rate_inject/fdum
  write(*,*) '------------------------------------------------------------------------'
endif

if (ny_inject.eq.3) then
  write(*,*) 'ny_inject=3 not tested yet'
  stop
endif

if (iinj2.eq.nx.or.jinj2.eq.nz) then
  write(*,*) 'Verify_Input: iinj2=nx or jinj2=nz not allowed (yet)'
  stop
endif

return
end


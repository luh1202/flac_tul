common /mesh/ x0,z0,rxbo,rzbo,sizez_x(10),sizez_y(10), &
    dx_rem,angle_rem,anglemin1,anglemint,topo_kappa,bottom_kappa, &
    velbc_l,velbc_r,nx,nz,nq,nzonx,nzony,nelz_x(10),nelz_y(10), &
    ny_rem,mode_rem,ntest_rem,v_min,v_max,ivis_shape,efoldc, &
    igeotherm,g_x0,g_y0c,g_amplitude,g_width,ydrsides,ynstressbc, &
    ny_inject,rate_inject,jinj1, jinj2, Tsol, Tliq, xlatheat,inj_period,ext_period, &
    ntap, nmass_update,nopbmax,iinj, iinj1,iinj2, nxinj, inj_count, Tcinj, ratfac, &
    nofbc,nofside(20),nbc1(20),nbc2(20),nbc(20),iac_rem, &
    bca(20),bcb(20),bcc(20),dx_init,iacret, &
    bcd(20),bce(20),bcf(20),bcg(20),bch(20),bci(20),iph_int, & 
    iynts,tbos, iax1,iay1,ibx1,iby1,icx1,icy1,idx1,idy1

common /particles/ np, npelem  ! G.Ito particles for tracing top 8/07

common /flac_ctrl/ dt_scale,strain_inert,strain_inert0,vbc,amul, &
    ratl,ratu,frac,dtmax_therm,dt_maxwell,fracm,srate0, &
    dt_elastic,demf,boff,mix_strain,mix_stress, &
    ivis_present,n_boff_cutoff,idt_scale,ifreq_imasses,ifreq_rmasses, &
	movegrid,ndim,ifreq_visc
    
character phasefile*20,tempfile*20,hv_out*20
common /inout/ dtout_screen,dtout_file,dtsave_file,lastsave,lastout, &
        io_vel,io_srII,io_eII,io_aps,io_sII,io_sxx,io_szz, &
		io_sxz,io_pres,io_temp,io_melt,io_visc,io_phas,io_mark, &
		io_src,io_cond,io_diss,io_forc,io_hfl,io_topo, &
		irphase,irtemp,phasefile,tempfile,hv_out,nprofil

common /rheol/ visc(20),den(20),alfa(20),beta(20),pln(20),acoef(20),eactiv(20),rl(20),rm(20), &
    coha(20),cohdisp(20),phimean(20),phidisp(20),psia(20),conduct(20),cp(20), &
	ts(20),tl(20),tk(20),fk(20),ten_off, &
	plstrain(5),fric(5),dilat(5),cohesion(5),xinitaps,tau_heal, &
	nphase,mphase,irheol(20),nysoft,nsegments, &
	inhom,ix1(20),ix2(20),iy1(20),iy2(20),inphase(20),igeom(20),nphasl, &
    ltop(20),lbottom(20),lphase(20),iprof_out(20),irh_mark

common /therm/ t_top,t_bot,hs,hr,temp_per,bot_bc,time_t,dt_therm, &
	itherm,istress_therm,initial_geoth,itemp_bc,ix1t,ix2t,iy1t,iy2t,ishearh, &
	ntherm

common /stress/ g,pisos,drosub,damp_vis,nyhydro,iphsub

common /process/ sec_year,time,dt,nloop,mloop, irestart,time_max, &
	ifreq_avgsr,ireset

common /hooks/ if_intrus,if_hydro

common /temporary/ devmax, dvmax

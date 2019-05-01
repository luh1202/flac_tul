program DREZINA

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

!Note: All references to command 'secnds' have been removed
!      for Linux compilation.  Mark Behn, 10.17.03.
!real*4 secnds,time0
real*4 time0

! Area
open( 33, file='area.dat' )

!time0 = secnds(0.0)

! seconds in a year
sec_year = 3.1558e+7

! Read task parameters
call read_params()

! Try to read save-file contents. If file exist - restart, othewise - new start
open(1,file='_contents.rs',status='old',err=10)

irestart = 1
close(1)
goto 20

10 irestart = 0

20 continue

if ( irestart .eq. 1 ) then  !file exists - restart
    call rsflac
    
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you CONTINUE from  ', nloop+1, ' step'
    else
        call SysMsg('you CONTINUE the execution')
    endif
else ! file does not exist - new start
    if( dtout_screen .ne. 0 ) then
        write(6,*) 'you have NEW start conditions'
    else
        call SysMsg('you have NEW start conditions')
    endif

    call setflac
    ! Output of initial configuration
    call outflac
end if


!      ****************** running ********************************
ireset = 1
dtacc_screen = 0
dtacc_file = 0
dtacc_save = 0
dtacc_remesh = 0
nrec = 1 

!if (npelem.gt.0) call particle_seed		!G.Ito 3/07

do while( time .lt. time_max )
    nloop = nloop + 1
    mloop = mloop + 1
 
    if( mod(nloop,1000) .eq. 0 ) then					!G.Ito
        write(*,'(a19,i12,f10.3)') 'nloop, time (kyr)=',nloop, time/sec_year/1.e+3 	!G. Ito
    endif								!G. Ito

    if( dtout_screen .ne. 0 ) then
    if( dtacc_screen .gt. dtout_screen ) then
       write(*,'(I7,A,I7,A,F9.6,A,F9.6,A,I7,A)') nloop,'''s step.',mloop,'''s mstep Time[My]=', &
                 time/sec_year/1.e+6,', dt=', dt/sec_year, ' inj_cnt=', inj_count
!      write(*,*) 'inj_count = ', inj_count, 'ext_period = ', ext_period, 'inj_period = ', inj_period
!      ', dt=', dt/sec_year, ',  elapsed-', secnds(time0)/60, ' min'

! Forces at the boundaries

      if( io_forc.eq.1 ) then
        force_l=0.
        force_r=0.
        do j = 1,nz-1
          sxx = 0.25 * (stress0(1,1,j,1)+stress0(1,2,j,1)+stress0(1,3,j,1)+stress0(1,4,j,1) )
          sxxd = sxx-stressI(1,j)
          dl = cord(j+1,1,2)-cord(j,1,2)
          force_l = force_l+abs(sxxd)*abs(dl)

          sxx = 0.25 * (stress0(1,1,j,nx-1)+stress0(1,2,j,nx-1)+stress0(1,3,j,nx-1)+stress0(1,4,j,nx-1) )
          sxxd = sxx-stressI(nx-1,j)
          dl = cord(j+1,nx-1,2)-cord(j,nx-1,2)
          force_r = force_r+abs(sxxd)*abs(dl)
        end do
        open (1,file='forc.0',access='direct',form='formatted',recl=28)
        write (1,'(f6.2,1x,e10.2,1x,e10.2)',rec=nrec) time/sec_year/1.e6, force_l, force_r
        nrec = nrec + 1
        close (1)
      endif

      dtacc_screen = 0
    endif
    endif

    ! FLAC
    call flac
 
!    if (npelem.gt.0) call particle_move		!G.Ito 3/07

    if( ireset.eq.1 ) ireset = 0

    ! Remeshing

    if( ny_rem.eq.1 .and. itherm.ne.2 ) then
        if( itest_mesh() .eq. 1 ) then
           if(iynts.ne. 1) call fl_therm
             call re_mesh
            ireset = 1
        endif
    endif

    ! REMESHING WITH KINEMATIC ACCRETION
    if(ny_inject.eq.3) then
    dt_remesh = (2.* dx_init)/(rate_inject)
    if (dtacc_remesh.ge.dt_remesh) then
       iacret = 1
        call re_mesh
        dtacc_remesh = 0   
    endif
    endif

    ! OUTPUT  
    if( dtout_file .ne. 0 ) then 
        if( dtacc_file .gt. dtout_file ) then
           write(*,*) nloop, 'out', 0.25*(stress0(1,1,1,100)+stress0(1,2,1,100)+stress0(1,3,1,100)+stress0(1,4,1,100))
           call outflac
           dtacc_file = 0
        endif
    endif

    ! SAVING
    if( dtsave_file .ne. 0 ) then 
        if( dtacc_save .gt. dtsave_file ) then
            call saveflac
            dtacc_save = 0
        endif
    endif

! Area
    if( mod(nloop,1000) .eq. 0 ) then
	    area_diff = total_area(0)/abs(rzbo*rxbo) - 1
!        write( *,'(i6,1x,e9.2,1x,e9.2,1x,e9.2)' ) nloop, area_diff, devmax, dvmax
        write(33,'(i6,1x,e9.2,1x,e9.2,1x,e9.2)' ) nloop, area_diff, devmax, dvmax
        devmax = 0; dvmax = 0;
        !call flush(33)
    endif

    time = time + dt

    dtacc_screen = dtacc_screen + dt
    dtacc_file   = dtacc_file   + dt
    dtacc_save   = dtacc_save   + dt
    dtacc_remesh = dtacc_remesh + dt   
 
end do

! Area
close(33)

call SysMsg('Congratulations !')
stop 999
end

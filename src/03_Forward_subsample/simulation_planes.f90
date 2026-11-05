!=======================================================================
! Highest algorithm level; controls overall simulation and runs the time loop
!
! simulation_init()
! simulation_read()
! simulation_write()
! simulation_run()
! simulation_info(iunit)
!

module simulation

	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	use data_export
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	implicit none
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	INTERFACE
		real*8 function realClock()
		end function realClock
	END INTERFACE

	logical iexist
	CHARACTER(LEN=47) ::   FN2
	CHARACTER(LEN=41) ::   NAME

	real cpu_cflcal , cpu_bc   , cpu_advx , cpu_invx    , cpu_advy
	real cpu_invy   , cpu_advz , cpu_invz , cpu_poisson , cpu_update
	real cpu_divchk , cpu_toxyz, cpu_init , cpu_output  , cpu_total_i , cpu_total
	real Re         , restart

	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	integer nsteps       ! number of steps to run
	real    timeLimit    ! stopping time
	
	integer ianim        ! whether to save data for animation
	real    saveInterval ! time interval to save data
	integer nsave        ! interval number for saving data
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	real	:: TimeLoop_time

        integer nsaveInterval_Rstrt !JL

        ! MZW for reading coarse-grained data
        integer :: istart, iend, iskip, jstart, jend, jskip, kstart, kend, kskip
        
end module

!=======================================================================
subroutine simulation_init()
	use simulation
        ! for computing flags
        integer :: flag_proc, my_istart, my_iend, my_jstart, my_jend
        integer :: loc_istart, loc_iend, loc_jstart, loc_jend        
        integer :: i, j, k
        integer :: ii, jj, kk

	cpu_cflcal =0.0
	cpu_bc     =0.0
	cpu_advx   =0.0
	cpu_invx   =0.0
	cpu_advy   =0.0
	cpu_invy   =0.0
	cpu_advz   =0.0
	cpu_invz   =0.0
	cpu_poisson=0.0
	cpu_update =0.0
	cpu_divchk =0.0
	cpu_toxyz  =0.0
	cpu_init   =0.0
	cpu_stat   =0.0 !JL
	cpu_write  =0.0 !JL
	cpu_output =0.0
	cpu_total_i=0.0
	cpu_total  =0.0
	t_init     = realClock()
	t_over     = realClock() - t_init
	before     = realClock()

	call arrayzero

	!--------------------------------------------------
	!MZW parameters for reading observation files
        call readInt("iassim",iassim)
        call readInt("ifric",ifric)
        !call readInt("ipre",ipre)
        call readInt("isubdomain",isubdomain)
	call readInt("ObsStart",ObsStart)
	call readInt("ObsEnd",ObsEnd)
	call readInt("ObsNt",ObsNt)
        if (isubdomain == 1) then
                call readInt("istart",istart)
                call readInt("iend",iend)
                call readInt("iskip",iskip)
                call readInt("jstart",jstart)
                call readInt("jend",jend)
                call readInt("jskip",jskip)
                call readInt("kstart",kstart)
                call readInt("kend",kend)
                call readInt("kskip",kskip)
        end if
	call readFloat("sigma",sigma)
	!--------------------------------------------------
        ! uflag, vflag, wflag marks which points have observation data
        uflag = 0;  vflag = 0;  wflag = 0

        if(isubdomain == 1) then
        my_istart = iTmin_1(iblock) 
        my_iend   = iTmax_1(iblock)
        my_jstart = jTmin_1(jblock)
        my_jend   = jTmax_1(jblock)

        if (iend.lt.my_istart .or. istart.gt.my_iend .or. jend.lt.my_jstart .or. jstart .gt. my_jend) then 
                flag_proc = 0
        else
                flag_proc = 1
        end if

        if (flag_proc == 1) then
                !-------------x direction-------------
                if (my_istart .le. istart) then
                        loc_istart = imap_1(istart)
                else
                        loc_istart = imap_1(istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
                end if

                if (my_iend .ge. iend) then
                        loc_iend = imap_1(iend)
                else
                        loc_iend = imap_1(my_iend - mod(my_iend-istart,iskip))
                end if
                !-------------y direction-------------
                if (my_jstart .le. jstart) then
                        loc_jstart = jmap_1(jstart)
                else
                        loc_jstart = jmap_1(jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
                end if

                if (my_jend .ge. jend) then
                        loc_jend   = jmap_1(jend)
                else
                        loc_jend   = jmap_1(my_jend - mod(my_jend-jstart,jskip))
                end if

        
        	do j = loc_jstart, loc_jend, jskip
        	do i = loc_istart, loc_iend, iskip
        	do k = kstart, kend, kskip
        	        uflag(k,i,j) = 1
        	        vflag(k,i,j) = 1
        	        wflag(k,i,j) = 1
        	end do
        	end do
        	end do
        end if

        if (istart == 1 .and. iblock == npx) then
                i = imap_1(ngx)
                do j = loc_jstart, loc_jend, jskip
                do k = kstart, kend, kskip
                        uflag(k,i,j) = 1
                end do
                end do
        end if

        if(ifric == 1 .and. jblock == 1) then
                j = jmap_1(1)
                do i = loc_istart, loc_iend, iskip
                do k = kstart, kend, kskip
                        uflag(k,i,j) = 1
                        wflag(k,i,j) = 1
                end do
                end do
        end if
        if(ifric == 1 .and. jblock == npy) then
                j = jmap_1(ngy-1)
                do i = loc_istart, loc_iend, iskip
                do k = kstart, kend, kskip
                        uflag(k,i,j) = 1
                        wflag(k,i,j) = 1
                end do
                end do
        end if
        end if  ! if isubdomain == 1
	!--------------------------------------------------
	call readFloat("Re",Re)         ! Reynolds number   
	!call readFloat("Fr",Fr)        ! Froude number
	!call readFloat("We",We)        ! Weber number

	visc = 1.0/Re                   !  inverse Reynolds number

	call readInt("nsteps", nsteps)
	call readInt("ipout", ipout)
	call readFloat("timeLimit", timeLimit)

	dt    = 0.
	stime = 0.
	ntime = 0

	call readInt("ianim", ianim)
	if (ianim > 0) call readFloat("save_dt", saveInterval)
	call readInt("nsave_Rstrt", nsaveInterval_Rstrt)  !JL
	nsave = 0

	call poisson_init_fft()  ! Fast Fourier transform
	write(6,*), 'FFT in Z and X directions!!!'

	call CFLcontrol_init()
	call advance_init()
        dt_old = dt  ! JL

	return
end

subroutine simulation_read()
	use simulation

	call readFloat("time", stime)
	call readInt("ntime", ntime)

	!T  Gets (dt)
	call CFLcontrol_read()

	return
end

subroutine simulation_write()
	use simulation

	call writeInt("nsteps", nsteps)
	call writeFloat("time", stime)
	call writeInt("ntime", ntime)
	call writeFloat("timeLimit", timeLimit)
	call writeInt("ianim", ianim)

	call writeFloat("Re",Re)         ! Reynolds number   
	!call writeFloat("Fr",Fr)         ! Froude number
	!call writeFloat("We",We)         ! Weber number

	call advance_write()
	call CFLcontrol_write()

	return
end

!=======================================================================
subroutine simulation_run()
	use simulation
        use CFLcontrol, only: iCFL       !JL
        character(len=7) :: fn_archive_backup !JL
        logical :: dir_e

	stime_ = stime
	ntime_ = ntime

	TimeLoop_time = realClock()

	! write flux at initial time
        if(iassim == 2) then
                call FluxBC()
                call WriteFluxes()
        end if

        open(101,file='finish_for')
        finish = 0
        write(101,*) finish
        close(101)

        open(33,file="cost_terms.dat")
        open(111,file='TKE_check.txt')
        call Flux_to_Cart()
        call CartesianBC()
        call energy_xz()

        cost = 0.0
	ntime_obs = ObsStart
	if(iassim == 1 .and. ntime_obs == ntime) then
                !if(ipre == 1) then
		!call WriteSubDomain_pcc('./obs/PreObs_bott.', 1,ngx-1,1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_pcc('./obs/PreObs_uppr.', 1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
                !end if

                !if(ifric == 1) then
		!call friction()
		!! friction at ghost point is zero!!!
		!call WriteSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
		!call WriteSubDomain_dudy('./obs/duyObs_bott.',1,ngx,1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_dudy('./obs/duyObs_uppr.',1,ngx,1, ngy-1,ngy-1,1, 1,ngz-1,1)
                !if(irank==1) write(*,*) 'dudy=',maxval(dudy_tmp)
                !end if

                if(isubdomain == 1) then
                call WriteSubDomain_ffc(istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip)
                end if
		ntime_obs = ntime_obs + ObsNt
	end if 

	if(iassim == 2 .and. ntime_obs == ntime) then
                !if(ipre == 1) then
		!call WriteSubDomain_pcc('/for/PreTmp_bott.', 1,ngx-1,1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_pcc('/for/PreTmp_uppr.', 1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
		!! ReadSubDomain_pcc ???
		!call ReadSubDomain_pcc('./obs/PreObs_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1)
		!p_obs(:,:,1) = p_obs(:,:,2)
		!call ReadSubDomain_pcc('./obs/PreObs_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1)
                !end if

                !if(ifric == 1) then
		!dwdy_obs = 0.0
		!dudy_obs = 0.0
		!call ReadSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
		!call ReadSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
		!call ReadSubDomain_dudy('./obs/duyObs_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1)
		!call ReadSubDomain_dudy('./obs/duyObs_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1)
                !if(irank==1) write(*,*) 'dudy=',maxval(dudy_obs)

		!call friction()
		!! friction at ghost point is zero!!!
		!call WriteSubDomain_dwdy('./for/dwyTmp_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_dwdy('./for/dwyTmp_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
		!call WriteSubDomain_dudy('./for/duyTmp_bott.',1,ngx,  1, 1,1,1, 1,ngz-1,1)
		!call WriteSubDomain_dudy('./for/duyTmp_uppr.',1,ngx,  1, ngy-1,ngy-1,1, 1,ngz-1,1)
                !end if

                if(isubdomain == 1) then
                uob = 0.0; vob = 0.0; wob = 0.0
                call ReadSubDomain_ffc(istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip)
                call FluxBC_uob()
                if(irank==1) write(*,*) 'finish reading subdomain files at ntime=',ntime
                end if
		call cost_function 
                

		ntime_obs = ntime_obs + ObsNt
	end if 

	do while (ntime < ntime_ + nsteps)
		!--------------------
		!    Compute CFL
		!--------------------
		before01 = realClock()
		if(mod(ntime-ntime_,ipout).eq.0 .and. iCFL/=1) & !JL
			call CFLcontrol_maxCFL()  ! fixed dt
                if(iCFL==1) &   !JL
			call CFLcontrol_adjust()  ! fixed CFL
		after01 = realClock()
		cpu_cflcal= (after01 - before01) - t_over  

		!----------------------------------------
		!    Set Boundary Conditions
		!----------------------------------------
		before02 = realClock()
		!	call timeDependent_Inlet_BC(dt)	        !--- Inlet FST + Convective outflow
		!	call timeDependent_Wall_BC()		!--- Wall oscillation
			call CartesianBC(dt)                    !--- No Slip   + Blasius upper surface
			call FluxBC()                           !--- From cartesian
			call CopyBC_out()                       !--- Copy BC to free (uc,vc,wc)
		after02 = realClock()
		cpu_bc = (after02 - before02) - t_over

		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! I. Advance the convection, explicit diffusion and
		!    implicit diffusion terms to get intermediate velocity
		! Note: 'invertu' saves the results of (ust) in (uc)
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!------------- u-equation -------------
		before04 = realClock()
			call advancex
		after04 = realClock()
		cpu_advx = (after04 - before04) - t_over

		before05 = realClock()  
			call invertu
		after05 = realClock()
		cpu_invx = (after05 - before05) - t_over

		!------------- v-equation -------------
		before06 = realClock()
			call advancey
		after06 = realClock()
		cpu_advy = (after06 - before06) - t_over
		
		before07 = realClock()
			call invertv
		after07 = realClock()
		cpu_invy = (after07 - before07) - t_over

		!------------- w-equation -------------
		before08 = realClock()
			call advancez
		after08 = realClock()
		cpu_advz = (after08 - before08) - t_over

		before09 = realClock()
			call invertw
		after09 = realClock() 
		cpu_invz = (after09 - before09) - t_over

		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	II.	POISSON	 SOLVER
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!------------ B.C. on (ust)  -------------
		call poissonRhsBC()

		!--------------- Solver ------------------
		before10 = realClock()
		!-------------------------------------------
		! FFT in z-dir. + FFT in x-dir.
		!-------------------------------------------
		call poisson_fftz_fftx(cpu_poisson)

		after10 = realClock()
		cpu_poisson = (after10 - before10) - t_over

		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	III.	Projection Step
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	Update (u,v,w) from (ust,vst,wst) & P
		!-------------------------------------------------
		before11 = realClock()
			call Update_Flux()
		after11 = realClock()
		cpu_update = (after11 - before11) - t_over

		before12 = realClock()
			call Divergence_Check()
		after12 = realClock()
		cpu_divchk = (after12 - before12) - t_over

		!------- Compute Cartesian Vel. --------
		before13 = realClock()   
			call CopyBC_in()        !--- Copy BC to free (uc,vc,wc)
			call Flux_to_Cart()
		after13 = realClock()
		cpu_toxyz = (after13 - before13) - t_over

		!------- Compute flowAnalysis -----
		! Write flow analysis file
		before14 = realClock()   
		call flowAnalysis_blayerinfo()

		!SJL ------- Vortex visualization (2nd invariance of velocity gradients) -------
		!call vortex_calc(vort_C,dudy_C)

                !------- Compute Statistics --------
                call statistics_sample
		after14 = realClock()
		cpu_stat = (after14 - before14) - t_over


		!---------------------------
		!  INCREMENT TIME COUNTERS
		!---------------------------
		stime = stime + dt
		ntime = ntime + 1
                dt_old = dt   !JL for constant CFL + Adams-Bashforth scheme
		!--------Compute cost function------
                call CartesianBC()
                call energy_xz()

		before15 = realClock()

		! write observations
		if(iassim == 1 .and. ntime_obs == ntime) then
                        !if(ipre==1) then
			!call WriteSubDomain_pcc('./obs/PreObs_bott.', 1,ngx-1,1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_pcc('./obs/PreObs_uppr.', 1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
                        !end if

                        !if(ifric == 1) then
			!call friction()
			!! friction at ghost point is zero!!!
			!call WriteSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
			!call WriteSubDomain_dudy('./obs/duyObs_bott.',1,ngx,1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_dudy('./obs/duyObs_uppr.',1,ngx,1, ngy-1,ngy-1,1, 1,ngz-1,1)
                        !if(irank==1) write(*,*) 'dudy=',maxval(dudy_tmp)
                        !end if

                        if(isubdomain == 1) then
                        call WriteSubDomain_ffc(istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip)
                        end if
			ntime_obs = ntime_obs + ObsNt
		end if

		! read observations and compute cost function
		if(iassim == 2 .and. ntime_obs == ntime) then
                        if(irank==1) write(*,*) 'write observation at time=',ntime_obs
                        !if(ipre == 1) then
			!call WriteSubDomain_pcc('./for/PreTmp_bott.', 1,ngx-1,1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_pcc('./for/PreTmp_uppr.', 1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
			!call ReadSubDomain_pcc('./obs/PreObs_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1)
			!p_obs(:,:,1) = p_obs(:,:,2)
			!call ReadSubDomain_pcc('./obs/PreObs_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1)
                        !if(irank == 1) then
                        !        write(*,*) 'p_tmp = ', maxval(p(:,:,1))
                        !        write(*,*) 'p_obs = ', maxval(p_obs)
                        !end if
                        !end if

                        !if(ifric == 1) then
			!dwdy_obs = 0.0
			!dudy_obs = 0.0
			!call ReadSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
			!call ReadSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
			!call ReadSubDomain_dudy('./obs/duyObs_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1)
			!call ReadSubDomain_dudy('./obs/duyObs_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1)

			!call friction()
			!! friction at ghost point is zero!!!
			!call WriteSubDomain_dwdy('./for/dwyTmp_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_dwdy('./for/dwyTmp_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1)
			!call WriteSubDomain_dudy('./for/duyTmp_bott.',1,ngx,  1, 1,1,1, 1,ngz-1,1)
			!call WriteSubDomain_dudy('./for/duyTmp_uppr.',1,ngx,  1, ngy-1,ngy-1,1, 1,ngz-1,1)
                        !if(irank==1) write(*,*) 'dudy=',maxval(dudy_obs)
                        !end if

                        if(isubdomain == 1) then
                        uob = 0.0; vob = 0.0; wob = 0.0
                        call ReadSubDomain_ffc(istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip)
                        call FluxBC_uob()
                        if(irank==1) write(*,*) 'finish reading subdomain files at ntime=',ntime
                        end if

			call cost_function()
			ntime_obs = ntime_obs + ObsNt
		end if

                if(iassim == 2) call WriteFluxes()

		!cccccccccccccccccccccccccccccccccccccccccc
		!  WRITE RESTART FILE for SCAVENGER QUEUE
		!cccccccccccccccccccccccccccccccccccccccccc
		!------ Write out data for restart file-------
                if (mod(ntime,nsaveInterval_Rstrt) == 0) then 
 			call WriteRstrt()
                        call statistics_write()
                        ! Write archive file
                        call writeFloat("time", stime)
                        call writeInt("ntime", ntime)
                        irestart = 1
                        call writeInt("irestart", irestart)
		
                	! Write archive file
                	call data_isRoot(iflag)
                	if (iflag == 1) then
                                iunit = iopen()
                               !JL-- write archive
                                open (iunit, file="archive", form="unformatted")
                                call archive_write(iunit, 1)
                                close (iclose(iunit))

                		iunit = iopen()
                               !JL-- write archive backup
                                write(fn_archive_backup,'(i7.7)') ntime
                		open (iunit, file="archive."//fn_archive_backup, form="unformatted")
                		call archive_write(iunit, 1)
                		close (iclose(iunit))
                	end if
		end if

		!cccccccccccccccccccccccccccccccccccccccccc
		!		ANIMATION
		!cccccccccccccccccccccccccccccccccccccccccc
		!------ Write out data for animation-------
		if (ianim == 1) then
			if (nint(stime/saveInterval) >= nsave) then
				nsave = nint(stime/saveInterval) + 1
				!JACOBS	call data_writeZdata()
				!call WriteZplane()
			end if
		end if
		if (ianim == 2) then
			if (nint(stime/saveInterval) >= nsave) then
				nsave = nint(stime/saveInterval) + 1
				!JACOBS	call data_writeYdata()
				!	call WriteYplane()
			end if
		end if
		!------------------------------------------
		if (ianim == 3) then
			if (nint(stime/saveInterval) >= nsave) then
                            if( nsave == 0 ) then
				nsave = nint(stime/saveInterval) + 1
                            else
				nsave = nint(stime/saveInterval) + 1
				call WriteSubDomain    ('SubDom_dble.', 1,ngx,2, 1,ngy,2, 1,ngz,2)
                            end if
			end if
		end if
		!------------------------------------------
		after15 = realClock()
		cpu_write = (after15 - before15) - t_over

		!cccccccccccccccccccccccccccccccccccccccccc
		!	Checks for Stopping the code
		!cccccccccccccccccccccccccccccccccccccccccc
		!--------- Check for stopping time --------
		if (timeLimit .ne. 0. .and. stime >= timeLimit) then
			call signalShell("done")
			exit
		end if

		!------- Check for user intervention ------
		if(mod(ntime-ntime_,20).eq.0) then
			inquire(file="resubmit", exist=iexist)
			if (iexist) then
				call signalShell("resubmit")
				! Remove file (only one process actually does this)
				! call system('rm resubmit')
				exit
			end if
			inquire(file="done", exist=iexist)
			if (iexist) then
				call signalShell("done")
				! Remove file (only one process actually does this)
				! call system('rm done')
				exit
			end if
		end if

		!--------- Check for CPU time --------
		tremain = timeRemaining()
		call globalMax(tremain, 1)
		if (tremain < 0.) then
			call signalShell("resubmit")
			exit
		end if

		cpu_total_i = cpu_cflcal+cpu_bc+cpu_advx+cpu_invx+cpu_advy+cpu_invy+ &
		cpu_advz+cpu_invz+cpu_update+cpu_divchk+cpu_toxyz+cpu_poisson+cpu_stat+cpu_write

		if (irank.eq.1) then
			write(6,'(A80)') & 
				'CPU cflcal    bc      advx      invx      advy      invy      advz             '
			write(6,'(7(1PE10.3),A15)') &
				cpu_cflcal,cpu_bc,cpu_advx,cpu_invx,cpu_advy,cpu_invy,cpu_advz,'             #-'
			write(6,'(A80)') & 
				'   invz     update    divchk     toxyz   poisson   Statistics  Write   total_i '
			write(6,'(8(1PE10.3),A5)') &
				cpu_invz,cpu_update,cpu_divchk,cpu_toxyz,cpu_poisson,cpu_stat,cpu_write,cpu_total_i,'   #_'
		end if

		cpu_total = cpu_total + cpu_total_i

	end do

	! write cost function into file
	open(100,file='cost_function')
	write(100,*) cost
	close(100)

        open(101,file='finish_for')
        finish = 1
        write(101,*) finish
        close(101)
        close(33)
        close(111)
        !call CartesianBC(dt)
        !call WriteRstrt()
	!print*,'Time loop timing = ', TimeLoop_time - realClock()

	if (ntime >= ntime_ + nsteps) then
		if (timeLimit .ne. 0. .and. stime + 0.5*dt < timeLimit) then
			call signalShell("continue")
		else
			call signalShell("done")
		end if
	end if

	if (irank.eq.iroot) then 
		write(6,'(A13,1PE10.3)') 'CPU total  = ',cpu_total
	end if

	return
end

!=======================================================================
subroutine simulation_info(iunit)
	use simulation

	call sectionTitle(iunit, "Simulation info")
	write (iunit, 11) ntime-ntime_, stime_, stime, stime-stime_
	
 11  format (4x, "Number of time steps this run : ", I13   / &
	     4x, "Time at beginning of run      : ", f13.6 / &
	     4x, "Time at end of run            : ", f13.6 / &
	     4x, "Time advancement              : ", f13.6)

	!T	call equations_info(iunit)

	return
end


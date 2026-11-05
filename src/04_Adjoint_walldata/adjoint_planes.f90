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
	integer ntime_obs    ! observation time
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
	! MZW copied from QW adjoint source code-----------
	call readInt("veloflag",veloflag)
        !call readInt("AdjointOperator",AdjointOperator)
	if((veloflag .eq. 2) .or. (veloflag .eq. 3))then
		call readInt("RestoreStart",RestoreStart)
		call readInt("RestoreEnd",RestoreEnd)
		call readInt("Restoregap",Restoregap)
		call readInt("AdjointNt",AdjointNt)
	endif
        if(veloflag .ne. 0)then
                call readInt("Adjointflag",Adjointflag)
	!	Adjointflag = 2
        endif
	!--------------------------------------------------

	!MZW parameters for reading observation files
	call readInt("ObsStart",ObsStart)
	call readInt("ObsEnd",ObsEnd)
	call readInt("ObsNt",ObsNt)
	call readFloat("sigma",sigma)
        call readInt("ifric",ifric)
        call readInt("ipre",ipre)
	!--------------------------------------------------
        ! uflag, vflag, wflag marks which points have observation data
        uflag = 0;  vflag = 0;  wflag = 0


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
!	write(6,*), 'FFT in Z and X directions!!!'

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
subroutine adjoint_run()
	use simulation
        use CFLcontrol, only: iCFL       !JL
        character(len=7) :: fn_archive_backup !JL

	stime_ = stime
	ntime_ = ntime

	TimeLoop_time = realClock()
	
	! for convenience, set initial adjoint variables the same as forward variables
	! MZW --------------------------------------------------------
	!u_adj = u; v_adj = v; w_adj = w
	!uc_adj= uc; vc_adj = vc; wc_adj = wc
	!ust_adj = ust; vst_adj = vst; wst_adj = wst
	!u = 0.0; v=0.0; w=0.0           ! for debug
	!uc = 0.0; vc = 0.0; wc=0.0      ! for debug
	! ------------------------------------------------------------
	ntime_obs = ObsEnd
	iobs = 0
	if(ntime_obs == ntime_+nsteps-ntime) then
		iobs = 1
                if(ipre == 1) then
			call ReadSubDomain_pcc('./for/PreTmp_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1,ntime_obs)
			p_tmp(:,:,1) = p_obs(:,:,2)
			p_obs(:,:,2) = 0.0
			call ReadSubDomain_pcc('./for/PreTmp_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1,ntime_obs)
			p_tmp(:,:,2) = p_obs(:,:,2)
			p_obs(:,:,2) = 0.0
			call ReadSubDomain_pcc('./obs/PreObs_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1,ntime_obs)
			p_obs(:,:,1) = p_obs(:,:,2)
			p_obs(:,:,2) = 0.0
			call ReadSubDomain_pcc('./obs/PreObs_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1,ntime_obs)

			call advancep_Adjoint()		! since u_adj v_adj w_adj = 0.0
                end if
		!-------------read observation and forward skin friction----------
                if(ifric == 1) then
			dwdy_obs = 0.0
			dudy_obs = 0.0
			call ReadSubDomain_dwdy('./for/dwyTmp_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dwdy('./for/dwyTmp_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./for/duyTmp_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./for/duyTmp_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			dwdy_tmp = dwdy_obs
			dudy_tmp = dudy_obs

			dwdy_obs = 0.0
			dudy_obs = 0.0
			call ReadSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./obs/duyObs_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./obs/duyObs_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)

			!call advancexyz_Adjoint(cpu_advx,cpu_advy,cpu_advz,iobs)
			!if(irank==1) then
                	!        write(*,*) 'p_obs', maxval(p_obs)
                	!        write(*,*) 'p_tmp', maxval(p_tmp)
                	!        write(*,*) 'dwdy_obs', maxval(dwdy_obs)
                	!        write(*,*) 'dwdy_tmp', maxval(dwdy_tmp)
			!	write(*,*) 'w_adj',maxval(w_adj)
			!end if
			if(irank==1) then
				write(*,*) 'reading forward pressure and observation at',ntime_obs
			end if
                end if

                call advancexyz_Adjoint(cpu_advx,cpu_advy,cpu_advz,iobs)
		ntime_obs = ntime_obs - ObsNt
	end if

	do while (ntime < ntime_ + nsteps)
		!--------------------
		!    Compute CFL
		!--------------------
		!!before01 = realClock()
		!!if(mod(ntime-ntime_,ipout).eq.0 .and. iCFL/=1) & !JL
		!!	call CFLcontrol_maxCFL()  ! fixed dt
                !!if(iCFL==1) &   !JL
		!!	call CFLcontrol_adjust()  ! fixed CFL
		!!after01 = realClock()
		!!cpu_cflcal= (after01 - before01) - t_over  
		!write(*,*) 'ntime=, ntime_+nsteps',ntime, ntime_ + nsteps
!		write(*,*) 'stime=',stime
		stime = stime + dt
		ntime = ntime + 1
                dt_old = dt   !JL for constant CFL + Adams-Bashforth scheme
		!----------------------------------------
		!    Set Boundary Conditions
		!----------------------------------------
		before02 = realClock()
		!	call timeDependent_Inlet_BC(dt)	        !--- Inlet FST + Convective outflow
		!	call timeDependent_Wall_BC()		!--- Wall oscillation
			!call CartesianBC_Adjoint(dt)           ! MZW
			call FluxBC_Adjoint()                  ! MZW actually we should deal with B.C. for adjoint later (after invert)
			! don't need CopyBC currently since ust_adj and uc_adj are not equavalent
			! might need to add something in the future
		after02 = realClock()
		cpu_bc_adj = (after02 - before02) - t_over

		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	I.	POISSON	 SOLVER
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!------------no need to deal with B.C. on (u_adj)  -------------
		before03 = realClock()
		!-------------------------------------------
		! FFT in z-dir. + FFT in x-dir.
		!-------------------------------------------
		call poisson_Adjoint_fftz_fftx()
		after03 = realClock()
		cpu_poisson = (after03 - before03) - t_over

		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	II.	Read forward flux and observation 
		!ccccccccccccccccccccccccccccccccccccccccccccccccci
		!-------------Read forward flux field------------
		before04 = realClock()
			call FluxTmpRstrt()
		after04 = realClock()
		cpu_readflux = (after04 - before04) - t_over
		!-------------Read observation data------------
		before05 = realClock()
			!call ObsTmpRstrt(1)	! to be developed in the future
		after05 = realClock()
		cpu_readobs = (after05 - before05) - t_over
		before06 = realClock()
			call FluxBC()
		after06 = realClock()
		cpu_bc = (after06 - before06) - t_over
!!!		u=u_save_adj; v=v_save_adj; w=w_save_adj
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	III.	Projection (Inverse) Step 
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!------------- u-equation -------------
		before07 = realClock() 
			call invertu_Adjoint()
		after07 = realClock()
		cpu_invx = (after07 - before07) - t_over

		!------------- v-equation -------------		
		before08 = realClock()
			call invertv_Adjoint()
		after08 = realClock()
		cpu_invy = (after08 - before08) - t_over

		!------------- w-equation -------------
		before09 = realClock()
			call invertw_Adjoint()
		after09 = realClock() 
		cpu_invz = (after09 - before09) - t_over

		!-------------B.C. for ust_adj, vst_adj, wst_adj---------
		before10 = realClock()
			call FluxBC_Adjoint_st()
		after10 = realClock()
		cpu_bc_st = after10 - before10
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	IV.	Advance Step 
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	solve (ust_adj,vst_adj,wst_adj) from (u_adj,v_adj,w_adj) & P
		!-------------------------------------------------
		! select the scheme (Adams-Bashforth or Euler) for advance_Adjoint
		before11 = realClock() 
			call advance_pre_Adjoint(nsteps)
		!------------- u-equation -------------


		!-------------read observation and forward skin friction----------
		iobs = 0
		if(ntime_obs == ntime_+nsteps-ntime .and. ntime_obs .ge. ObsStart .and. ifric==1) then
			iobs = 1
			dwdy_obs = 0.0
			dudy_obs = 0.0
			call ReadSubDomain_dwdy('./for/dwyTmp_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dwdy('./for/dwyTmp_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./for/duyTmp_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./for/duyTmp_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			dwdy_tmp = dwdy_obs
			dudy_tmp = dudy_obs
	
			dwdy_obs = 0.0
			dudy_obs = 0.0
			call ReadSubDomain_dwdy('./obs/dwyObs_bott.',1,ngx-1,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dwdy('./obs/dwyObs_uppr.',1,ngx-1,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./obs/duyObs_bott.',1,ngx  ,1, 1,1,1, 1,ngz-1,1,ntime_obs)
			call ReadSubDomain_dudy('./obs/duyObs_uppr.',1,ngx  ,1, ngy-1,ngy-1,1, 1,ngz-1,1,ntime_obs)
		else
			dwdy_tmp = 0.0;  dwdy_obs = 0.0
			dudy_tmp = 0.0;  dudy_obs = 0.0				
		end if


		call advancexyz_Adjoint(cpu_advx,cpu_advy,cpu_advz,iobs)
		after11 = realClock()
		cpu_advxyz = (after11 - before11) - t_over
		
		!------------- v-equation -------------		
		!before11 = realClock()
		!	call advancey_Adjoint()
		!after11 = realClock()
		!cpu_advy = (after11 - before11) - t_over

		!------------- w-equation -------------
		!before12 = realClock()
		!	call advancez_Adjoint()
		!after12 = realClock() 
		!cpu_advz = (after12 - before12) - t_over

		!--------------p_adj(only for dp form)-----------
		before13 = realClock()

		if(ntime_obs == ntime_+nsteps-ntime .and. ntime_obs .ge. ObsStart) then
                        if(ipre == 1) then
				call ReadSubDomain_pcc('./for/PreTmp_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1,ntime_obs)
				p_tmp(:,:,1) = p_obs(:,:,2)
				call ReadSubDomain_pcc('./for/PreTmp_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1,ntime_obs)
				p_tmp(:,:,2) = p_obs(:,:,2)
				call ReadSubDomain_pcc('./obs/PreObs_bott.',1,ngx-1,1,1,1,1,1,ngz-1,1,ntime_obs)
				p_obs(:,:,1) = p_obs(:,:,2)
				call ReadSubDomain_pcc('./obs/PreObs_uppr.',1,ngx-1,1,ngy-1,ngy-1,1,1,ngz-1,1,ntime_obs)
				if(irank==1) then
					write(*,*) 'reading forward pressure and observation at', ntime_obs
				end if
                        end if
			ntime_obs = ntime_obs - ObsNt
		else
			p_tmp = 0.0
			p_obs = 0.0				
		end if
			! advancep_Adjoint() does not include information about ntime_obs
			call advancep_Adjoint()
		after13 = realClock()
		cpu_advp = (after13 - before13)


		cpu_total_i = cpu_bc_adj+cpu_poisson+cpu_readflux+cpu_readobs+cpu_bc+cpu_invx+ &
		cpu_invy+cpu_invz+cpu_bc_st+cpu_advxyz+cpu_advp+cpu_write

		if (irank.eq.1) then
			write(6,'(A75)') & 
				'CPU bc_adj poisson   readflux  readobs   bc        invx      invy      invz   '
			write(6,'(8(1PE10.3),A15)') &
				cpu_bc_adj,cpu_poisson,cpu_readflux,cpu_readobs,cpu_bc,cpu_invx,cpu_invy,cpu_invz,'      #-'
			write(6,'(A70)') & 
				' bc_st     advx      advy      advz      advp      Write     total_i   '
			write(6,'(7(1PE10.3),A5)') &
				cpu_bc_st,cpu_advx,cpu_advy,cpu_advz,cpu_advp,cpu_write,cpu_total_i,'   #_'
		end if

		cpu_total = cpu_total + cpu_total_i

	end do

        umax = maxval(abs(u_adj))
        vmax = maxval(abs(v_adj))
        wmax = maxval(abs(w_adj))
        call globalMax(umax,1)
        call globalMax(vmax,1)
        call globalMax(wmax,1)
        if(irank==1) write(*,*) 'uadj,vadj,wadj=',umax,vmax,wmax

	call WriteAdjFlux()
	call Flux_to_Cart_adj()
	call CartesianBC_Adjoint(dt)
!!!
        umax = maxval(abs(uc_adj))
        vmax = maxval(abs(vc_adj))
        wmax = maxval(abs(wc_adj))
        call globalMax(umax,1)
        call globalMax(vmax,1)
        call globalMax(wmax,1)
        if(irank==1) write(*,*) 'ucadj,vcadj,wcadj=',umax,vmax,wmax

	call readInt("isymproj",isymproj)
	if(isymproj == 1) then
                if(irank==1) write(*,*) 'use symmetric projector after adjoint!'
                call sym_proj_poisson_init_vel()
                call Cart_to_flux_adj ! now u_adj becomes S * uc_adj, not the u_adj used in the main loop
                call FluxBC_Adjoint()
                call Divergence_check_adj()
                ust_adj = uc_adj
                vst_adj = vc_adj
                wst_adj = wc_adj
                call sym_proj_poisson_vel_adj(cpu_poisson_total)
                call sym_proj_update_vel_adj()

                call CartesianBC_Adjoint(dt)
                call Cart_to_flux_adj()
                call FluxBC_Adjoint()
                call Divergence_check_adj()

        end if

	!call WriteAdjFlux()
        call WriteAdjCart()
	!-----------------------------------------------
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

subroutine Divergence_Check_adj()
!ccccccccccccccccccccccccccccc
!  Compute the divergence
!ccccccccccccccccccccccccccccc
	use data_export

        integer :: idivmax, jdivmax, kdivmax

        real    :: divmax2 !JL
        integer :: idivmax2, jdivmax2, kdivmax2 !JL

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the maximum divergence
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	divmax=-1.e6
	divmax2=-1.e6 !JL
	do j=j1_T,j2_T
	do i=i1_T,i2_T
	do k=1   ,ngzm
		ii = ibmap_1(i)
		jj = jbmap_1(j)

		! divmax=max(divmax, &
		! 	abs(u(k,i+1,j)-u(k,i,j)+v(k,i,j+1)-v(k,i,j)+w(k+1,i,j)-w(k,i,j))  /  &
		! 	( abs(u(k,i+1,j))+abs(u(k,i,j))+abs(v(k,i,j+1))+abs(v(k,i,j))+abs(w(k+1,i,j))+abs(w(k,i,j)) )    )
		!divmax=max(divmax, &
		!	abs(u(k,i+1,j)-u(k,i,j)+v(k,i,j+1)-v(k,i,j)+w(k+1,i,j)-w(k,i,j))/vp(ii,jj)  )

		divmaxtmp  = abs( u_adj(k,i+1,j)-u_adj(k,i,j)+v_adj(k,i,j+1)-v_adj(k,i,j)+w_adj(k+1,i,j)-w_adj(k,i,j))/vp(ii,jj)

		if (divmaxtmp.gt.divmax) then
			divmax = divmaxtmp;
			idivmax = ii ;  jdivmax = jj ;  kdivmax = k ;
		end if
		! SJL: commented out for fixing spurious oscillations in spanwise variation on 10/26/2014
		!if(abs(w(k,i,j)).le.1.e-12) w(k,i,j)=0.0
	end do
	end do
	end do

	!call globalMax(divmax, 1)
	divmaxtmp = divmax
	call globalMax(divmaxtmp, 1)

	if(mod(ntime,ipout).eq.0) then

		!if(abs(divmax-divmaxtmp)/divmax .le. 1e-10) then

		!SJL: prevented a potential bug in debugging mode on 10/26/2014
		if(abs(divmax-divmaxtmp) .le. 1e-10*abs(divmax)) then
			write(*,'(a,es12.5,3(a,i5))') 'Max. divergence is' , divmax, ' @', idivmax,',', jdivmax,',', kdivmax
		end if

	!if (irank.eq.1)  write(6,3000) divmax
 !3000     format(1x,' Maximum divergence is ', 1pe15.8)
	end if

	return
end

! only correct for channel flow
subroutine Cart_to_Flux_adj
        use data_export

        do jj=j1_u,j2_u
        do ii=i1_u,i2_u
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                u_adj(kk,ii,jj) = uc_adj(kk,ii,jj) * sqrt(suxix(i,j)**2+suxiy(i,j)**2)
        end do
        end do
        end do

        do jj=j1_v,j2_v
        do ii=i1_v,i2_v
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                v_adj(kk,ii,jj) = vc_adj(kk,ii,jj) * sqrt(svetax(i,j)**2+svetay(i,j)**2)
        end do
        end do
        end do

        do jj=j1_w,j2_w
        do ii=i1_w,i2_w
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                w_adj(kk,ii,jj) = wc_adj(kk,ii,jj) * swz(i,j)
        end do
        end do
        end do

end

subroutine Flux_to_Cart_adj
        use data_export

        do jj=j1_u,j2_u
        do ii=i1_u,i2_u
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                uc_adj(kk,ii,jj) = u_adj(kk,ii,jj) * sqrt(suxix(i,j)**2+suxiy(i,j)**2)
        end do
        end do
        end do

        do jj=j1_v,j2_v
        do ii=i1_v,i2_v
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                vc_adj(kk,ii,jj) = v_adj(kk,ii,jj) * sqrt(svetax(i,j)**2+svetay(i,j)**2)
        end do
        end do
        end do

        do jj=j1_w,j2_w
        do ii=i1_w,i2_w
        do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
                wc_adj(kk,ii,jj) = w_adj(kk,ii,jj) * swz(i,j)
        end do
        end do
        end do

end

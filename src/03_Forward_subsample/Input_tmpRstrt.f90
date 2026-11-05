!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   SAVE THE RESTART IMAGE TO ONE SINGLE FILE USING              *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================


module Input_mod_tmp
	use messenger
	use data_export
	implicit none
	include "mpif.h"
	
	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)      :: global_indices, local_indices
	INTEGER                    :: filetype, memtype
	!CHARACTER(LEN=18)          :: FN2
	CHARACTER*19     FN2
	CHARACTER(LEN=12)          :: NAME
	
	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0
        real*8  :: Interpratio
        integer :: Nbefore, Nafter, Nobs
end module

! ! MZW: only used when observing u, v, w
!subroutine ObsTmpRstrt(flag)
!        use Input_mod_tmp
!
!	!----------------------------------------------------
!	! (InBuffer) is allocated of exact size as input data.
!	! It circumvents a bug in ALC MPI-IO
!	!---------------------------------------------------
!	real , allocatable :: InBuffer(:,:,:)
!	integer :: flag
!
!        ! here we assume 
!	! ObsENd = RestoreEnd
!	! ObsStart = RestoreStart
!	if (flag == 0) then		!read during forward simulation
!		Nobs = ntime/ObsNt
!	else if(flag == 1) then
!		Nobs = ObsEnd - ntime/ObsNt	!read during adjoint simulation 
!	end if
!
!	call MPI_TYPE_SIZE(MPI_REAL8, FloatSize, ierr)
!	
!	!********************************************************
!	! Note (ntime, ntime_, stime, stime_) are global variables
!	! ntime  :  time step at end of last simulation ==> begining of this simulation
!	! ntime_ :  time step at begining of this simulation
!	CALL MPI_BARRIER(icomm_grid,IERR)
!
!       !=========================================================
!       !                 Observation
!       !=========================================================
!        NAME='04uobs.dble.'
!	call InFileName_tmp(Nobs)
!        CALL MPI_FILE_OPEN(icomm_grid, FN2, &
!                           MPI_MODE_RDONLY, &
!                           MPI_INFO_NULL, fh, ierr)
!        
!        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
!        global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
!        if (iblock.eq. 1 ) global_indices(2) = 0
!        if (jblock.eq. 1 ) global_indices(3) = 0
!        
!        !ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
!        !ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
!        !ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
!        		local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
!			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
!			if (iblock.eq. 1 ) is = 0					!ALC_FIX
!			if (jblock.eq. 1 ) js = 0					!ALC_FIX
!        !------------------------ u  ---------------------------
!        global_cnt = 0
!        offset     = 0 
!        
!        gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
!        lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
!        if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
!        if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
!	if ( npx  .eq. 1 )   lsizes(2)=ngx+2
!        if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
!        if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
!	if ( npy  .eq. 1 )   lsizes(3)=ngy+1
!       
!	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)	
!			allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
!			memsizes = lsizes							!ALC_FIX
!
!        CALL Create_commit_inputview_tmp()
!        CALL Create_commit_in_subarray_tmp()
!        CALL MPI_FILE_READ_ALL(fh, InBuffer , 1, memtype, status, ierr)
!        CALL MPI_FILE_CLOSE(fh, ierr)
!		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
!		u_obs(ks:ke, is:ie, js:je) = InBuffer !dble(InBuffer)					!ALC_FIX
!		deallocate(InBuffer)							!ALC_FIX
!
!        NAME='05vobs.dble.'
!	call InFileName_tmp(Nobs)
!!        if (irank.eq.iroot) write(*,*) 'reading  ', FN2
!
!        CALL MPI_FILE_OPEN(icomm_grid, FN2, &
!                           MPI_MODE_RDONLY, &
!                           MPI_INFO_NULL, fh, ierr)	
!        
!        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
!        global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
!        if (iblock.eq. 1 ) global_indices(2) = 0
!        if (jblock.eq. 1 ) global_indices(3) = 0
!        
!        !ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
!        !ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
!        !ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
!        		local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
!			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
!			if (iblock.eq. 1 ) is = 0					!ALC_FIX
!			if (jblock.eq. 1 ) js = 0					!ALC_FIX
!
!       !---------------- v -----------------
!        global_cnt = 0
!        offset     = 0 
!
!        gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
!        lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
!        if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
!        if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
!	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
!        if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
!        if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
!	if ( npy  .eq. 1 )   lsizes(3)=ngy+2
!
!	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
!			allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))		!ALC_FIX
!			memsizes = lsizes						!ALC_FIX
!        CALL Create_commit_inputview_tmp()
!        CALL Create_commit_in_subarray_tmp()
!        CALL MPI_FILE_READ_ALL(fh, InBuffer , 1, memtype, status, ierr)
!        CALL MPI_FILE_CLOSE(fh, ierr)
!		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
!		v_obs(ks:ke, is:ie, js:je) = InBuffer !dble(InBuffer)					!ALC_FIX	
!		deallocate(InBuffer)							!ALC_FIX	
!
!        NAME='06wobs.dble.'
!	call InFileName_tmp(Nobs)
!!        if (irank.eq.iroot) write(*,*) 'reading  ', FN2
!        CALL MPI_FILE_OPEN(icomm_grid, FN2, &
!                           MPI_MODE_RDONLY, &
!                           MPI_INFO_NULL, fh, ierr)
!        
!        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
!        global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
!        if (iblock.eq. 1 ) global_indices(2) = 0
!        if (jblock.eq. 1 ) global_indices(3) = 0
!        
!        !ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
!        !ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
!        !ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
!        		local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
!			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
!			if (iblock.eq. 1 ) is = 0					!ALC_FIX
!			if (jblock.eq. 1 ) js = 0					!ALC_FIX
!
!       !--------------- w  ------------------
!        global_cnt = 0
!        offset     = 0 
!
!        gsizes = (/ ngz+2 , ngx+1 , ngy+1 /) 
!        lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
!        if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
!        if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
!	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
!        if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
!        if (jblock.eq.npy)   lsizes(3)=nlyb-noy1-0+1
!	if ( npy  .eq. 1 )   lsizes(3)=ngy+1
!
!	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
!			allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))		!ALC_FIX
!			memsizes = lsizes						!ALC_FIX
!        CALL Create_commit_inputview_tmp()
!        CALL Create_commit_in_subarray_tmp()
!        CALL MPI_FILE_READ_ALL(fh, InBuffer , 1, memtype, status, ierr)
!        CALL MPI_FILE_CLOSE(fh, ierr)
!		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
!		w_obs(ks:ke, is:ie, js:je) = InBuffer !dble(InBuffer)					!ALC_FIX
!		deallocate(InBuffer)							!ALC_FIX
!
!	!CALL MPI_FILE_CLOSE(fh, ierr)
!        CALL MPI_BARRIER(icomm_grid,IERR)
!
!return
!
!
!end



subroutine  Create_commit_inputview_tmp()
        use Input_mod_tmp

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, global_indices, &
                        MPI_ORDER_FORTRAN, MPI_REAL8,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_REAL8, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	FILE_FLAG = 1
        return
end

subroutine Create_commit_in_subarray_tmp()
        use Input_mod_tmp

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes, lsizes, local_indices, &
                        MPI_ORDER_FORTRAN, MPI_REAL8,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1
        return
end

subroutine InFileName_tmp(ntmp)
	use Input_mod_tmp
        integer ntmp
        IF(ntmp.LE.9                            ) &
        WRITE(FN2,'(A12,A6,I1)') NAME,'000000',ntmp
        IF(ntmp.GE.10     .AND. ntmp.LE.99     ) &
        WRITE(FN2,'(A12,A5,I2)') NAME,'00000' ,ntmp
        IF(ntmp.GE.100    .AND. ntmp.LE.999    ) &
        WRITE(FN2,'(A12,A4,I3)') NAME,'0000'  ,ntmp
        IF(ntmp.GE.1000   .AND. ntmp.LE.9999   ) &
        WRITE(FN2,'(A12,A3,I4)') NAME,'000'   ,ntmp
        IF(ntmp.GE.10000  .AND. ntmp.LE.99999  ) &
        WRITE(FN2,'(A12,A2,I5)') NAME,'00'    ,ntmp
        IF(ntmp.GE.100000 .AND. ntmp.LE.999999 ) &
        WRITE(FN2,'(A12,A1,I6)') NAME,'0'     ,ntmp
        IF(ntmp.GE.1000000.AND. ntmp.LE.9999999) &
        WRITE(FN2,'(A12,   I7)') NAME,        ntmp
	!MZW
	if(irank == 1) then
		write(*,*) 'reading:', FN2 
	end if
	return
end

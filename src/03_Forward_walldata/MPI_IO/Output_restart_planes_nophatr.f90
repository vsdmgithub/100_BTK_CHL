!=======================================================================
!     ******************************************************************
!     *                                                                *
!     *   SAVE THE RESTART IMAGE TO ONE SINGLE FILE USING              *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================


module Output_mod
	use messenger
	use data_export
	implicit none
	include "mpif.h"
	! include "mpiof.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER			   :: filetype, memtype
	CHARACTER(LEN=25)          :: FN2
	CHARACTER(LEN=29)	   :: FN3
	CHARACTER(LEN=18)          :: NAME
	CHARACTER(LEN=12)	   :: NAME_obs

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan
end module



subroutine WriteRstrt()
use Output_mod

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	double precision, allocatable :: OutBuffer(:,:,:)

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!************************************************
	! Note (ntime, ntime_, stime, stime_) are global variables
	! ntime  :  time step at end of this simulation
	! ntime_ :  time step at begining of this simulation

	aan=1.0*ntime

	!===================================
	!   Output (ntime, stime, aan)
	!===================================
	IF (irank.EQ.iroot) THEN 
		NAME='ucvcwc.data.'
		call OutFileName()
		write(*,*) FN2
		open(29,file=FN2,form='formatted')
		write(29,'(i10,2e16.8)') ntime,stime,aan
		close(29)
	END IF

	! synchronize all processors
	CALL MPI_BARRIER(icomm_grid,IERR)

	!=======================================================
	!         CARTESIAN VELOCITIES
	!=======================================================
	NAME='01uc00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	call MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ uc ---------------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +2
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = uc(ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX




	NAME='02vc00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	call MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0


	!------------------------ vc ---------------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +2

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = vc(ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX




	NAME='03wc00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	call MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0


	!--------------- wc ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1
     
	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = wc(ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



	!=========================================================
	!                 FLUXES
	!=========================================================
	NAME='04uu00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ u  ---------------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +2
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = u (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX




	NAME='05vv00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!---------------- v -----------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +2

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = v (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



	NAME='06ww00.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!--------------- w  ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = w (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



	!=========================================================
	!                  CONVECTIVE TERMS 
	!=========================================================
	NAME='conold.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	call MPI_FILE_OPEN(icomm_grid, FN2,                 & 
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	!  Note:  MPI assumes here that numbering starts from zero
	!  Since my numbering starts from (one), I have to subtract
	!  (one) from every indeex
	!-------------------------------------------------------
	global_indices = (/ 0 , iTmin_1(iblock)-1 , jTmin_1(jblock)-1 /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1-1 , noy1-1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 1;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 1					!ALC_FIX
			if (jblock.eq. 1 ) js = 1					!ALC_FIX

       !--------------- conx ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz-1+1 , ngx-1+1 , ngy-1+1 /)					!TJ: xperiodic fix
	lsizes = (/ ngz-1+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)				!TJ: xperiodic fix
	if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(1   )+1
	if (iblock.eq.npx)   lsizes(2)=(nlxb-1   )-(nox1)+1+1				!TJ: xperiodic fix
	if ( npx  .eq. 1 )   lsizes(2)=(ngx -1+1 )					!TJ: xperiodic fix
	if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(1   )+1
	if (jblock.eq.npy)   lsizes(3)=(nlyb-1   )-(noy1)+1+1				!TJ: xperiodic fix
	if ( npy  .eq. 1 )   lsizes(3)=(ngy -1+1 )					!TJ: xperiodic fix

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = conx(ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz-1 , nlx-1 , nly-1 /)
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

       !---------------- cony -----------------
	global_cnt = global_cnt + (ngx-1+1)*(ngy-1+1)*(ngz-1+1)				!TJ: xperiodic fix
	offset = global_cnt * FloatSize

	OutBuffer = cony(ks:ke, is:ie, js:je)					!ALC_FIX
	CALL Create_commit_fileview()
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

       !--------------- conz ------------------
	global_cnt = global_cnt + (ngx-1+1)*(ngy-1+1)*(ngz-1+1)				!TJ: xperiodic fix
	offset = global_cnt*FloatSize

	!OR  global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
	!OR  global_cnt = global_cnt + product(gsizes)

	OutBuffer = conz(ks:ke, is:ie, js:je)					!ALC_FIX
	CALL Create_commit_fileview()
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	deallocate(OutBuffer)							!ALC_FIX

	!----------- codetime ------------
	!---- Should all be in Archive ---
	!  global_cnt = global_cnt + (ngx-1)*(ngy-1)*(ngz-1)
	!  offset = global_cnt*FloatSize
	!  IF (irank.EQ.iroot) THEN
	!     call MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, &
	!                            MPI_DOUBLE_PRECISION, 'native',   &
	!                            MPI_INFO_NULL, ierr)
	!     print *, 'stime = ', stime
	!     call MPI_FILE_WRITE(fh, stime, 1, MPI_DOUBLE_PRECISION, &
	!                             status, ierr)
	!     print *, 'aan      = ', aan
	!     call MPI_FILE_WRITE(fh, aan, 1, MPI_DOUBLE_PRECISION, &
	!                             status, ierr)
	!  END IF

	CALL MPI_FILE_CLOSE(fh, ierr)



	!=========================================================
	!                 PRESSURE
	!=========================================================
	NAME='pres_p.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	call MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!--------------- p ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(0)   +1
	if (iblock.eq.npx)   lsizes(2)=(nlxb+0)   -(nox1)+1
	if ( npx  .eq. 1 )   lsizes(2)=(ngx+1)
	if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(0)   +1
	if (jblock.eq.npy)   lsizes(3)=(nlyb+0)   -(noy1)+1
	if ( npy  .eq. 1 )   lsizes(3)=(ngy+1)

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = p(ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
       ! CALL MPI_BARRIER(icomm_grid,IERR)
       ! if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from processor=',irank

       !==========================================================
       !     Dump Out Pressure (phatr) to check I/O working properly
       !----------------------------------------------------------
       !
       !  NAME='pressure_ph.'
       !  WRITE(FN2,'(A12,A5,I1)') NAME,'00000',irank
       !  write(*,*) FN2
       !  open(29,file=FN2,form='unformatted')
       !  write(29) phatr
       !  close(29)
  
        call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
        call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	!=========================================================
	!                 DATA
	!=========================================================
	IF (irank.EQ.iroot) THEN
		open(28,file='data',form='formatted')
		write(28,*) ntime
		close(28)
	END IF

	aan=1.0*ntime
	!========================================================


       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
        CALL MPI_BARRIER(icomm_grid,IERR)
        if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from  processor=',irank

return
end

subroutine WriteFluxes()
use Output_mod

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	double precision, allocatable :: OutBuffer(:,:,:)

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!************************************************
	! Note (ntime, ntime_, stime, stime_) are global variables
	! ntime  :  time step at end of this simulation
	! ntime_ :  time step at begining of this simulation

	aan=1.0*ntime

	!===================================
	!   Output (ntime, stime, aan)
	!===================================
	IF (irank.EQ.iroot) THEN 
		NAME='./for/ucvcwc.data.'
		call OutFileName_long()
		write(*,*) FN2
		open(29,file=FN2,form='formatted')
		write(29,'(i10,2e16.8)') ntime,stime,aan
		close(29)
	END IF

	! synchronize all processors
	CALL MPI_BARRIER(icomm_grid,IERR)


	!=========================================================
	!                 FLUXES
	!=========================================================
	NAME='./for/04uu00.dble.'
	call OutFileName_long()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ u  ---------------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +2
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = u (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX




	NAME='./for/05vv00.dble.'
	call OutFileName_long()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!---------------- v -----------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +2

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = v (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



	NAME='./for/06ww00.dble.'
	call OutFileName_long()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!--------------- w  ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = w (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



					!ALC_FIX



       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
       ! CALL MPI_BARRIER(icomm_grid,IERR)
       ! if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from processor=',irank

       !==========================================================
       !     Dump Out Pressure (phatr) to check I/O working properly
       !----------------------------------------------------------
       !
       !  NAME='pressure_ph.'
       !  WRITE(FN2,'(A12,A5,I1)') NAME,'00000',irank
       !  write(*,*) FN2
       !  open(29,file=FN2,form='unformatted')
       !  write(29) phatr
       !  close(29)
  
        call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
        call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	!=========================================================
	!                 DATA
	!=========================================================
	IF (irank.EQ.iroot) THEN
		open(28,file='data',form='formatted')
		write(28,*) ntime
		close(28)
	END IF

	aan=1.0*ntime
	!========================================================


       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
        CALL MPI_BARRIER(icomm_grid,IERR)
        if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from  processor=',irank

return
end

subroutine WriteObsFluxes()
! write observation
use Output_mod

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	double precision, allocatable :: OutBuffer(:,:,:)
	double precision mean, noise
	integer seed
	!double precision sigma
	!----------------------------------------------------
	! Set noise parameter
	! Gaussian distribution
	!----------------------------------------------------
	seed = 12345*(irank+1)		! different random series in different processors
	mean = 0.0
	!sigma = 1.e-3			in this version, sigma comes from input file



	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!************************************************
	! Note (ntime, ntime_, stime, stime_) are global variables
	! ntime  :  time step at end of this simulation
	! ntime_ :  time step at begining of this simulation

	aan=1.0*ntime

	!===================================
	!   Output (ntime, stime, aan)
	!===================================
	IF (irank.EQ.iroot) THEN 
		NAME='ucvcwc.data.'
		call OutFileName()
		write(*,*) FN2
		open(29,file=FN2,form='formatted')
		write(29,'(i10,2e16.8)') ntime,stime,aan
		close(29)
	END IF

	! synchronize all processors
	CALL MPI_BARRIER(icomm_grid,IERR)


	!=========================================================
	!                 FLUXES
	!=========================================================
	NAME='04uobs.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ u  ---------------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +2
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = u (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX

	!----------------------------------------------------
	! Add noise
	! Gaussian distribution
	!----------------------------------------------------
	do jj=1,je-js+1
	do ii=1,ie-is+1
	do kk=1,ke-ks+1
		noise = r8_normal_ab(mean,sigma,seed)
		OutBuffer(kk,ii,jj) = OutBuffer(kk,ii,jj) + noise
	end do
	end do
	end do

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX




	NAME='05vobs.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!---------------- v -----------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +2

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = v (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	!----------------------------------------------------
	! Add noise
	! Gaussian distribution
	!----------------------------------------------------
	do jj=1,je-js+1
	do ii=1,ie-is+1
	do kk=1,ke-ks+1
		noise = r8_normal_ab(mean,sigma,seed)
		OutBuffer(kk,ii,jj) = OutBuffer(kk,ii,jj) + noise
	end do
	end do
	end do

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



	NAME='06wobs.dble.'
	call OutFileName()
	if (irank.eq.iroot) write(*,*) FN2

	CALL MPI_FILE_OPEN(icomm_grid, FN2,                 &
			   MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
	local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
	ks = 0;		is = nox1;		js = noy1		!ALC_FIX
	if (iblock.eq. 1 ) is = 0					!ALC_FIX
	if (jblock.eq. 1 ) js = 0					!ALC_FIX


	!--------------- w  ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb  -nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx +1
	if (jblock.eq. 1 )   lsizes(3)=nlyb  -noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy +1

	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
	ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
	OutBuffer = w (ks:ke, is:ie, js:je)					!ALC_FIX
	memsizes = lsizes							!ALC_FIX
	!----------------------------------------------------
	! Add noise
	! Gaussian distribution
	!----------------------------------------------------
	do jj=1,je-js+1
	do ii=1,ie-is+1
	do kk=1,ke-ks+1
		noise = r8_normal_ab(mean,sigma,seed)
		OutBuffer(kk,ii,jj) = OutBuffer(kk,ii,jj) + noise
	end do
	end do
	end do

	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
	deallocate(OutBuffer)							!ALC_FIX



					!ALC_FIX



       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
       ! CALL MPI_BARRIER(icomm_grid,IERR)
       ! if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from processor=',irank

       !==========================================================
       !     Dump Out Pressure (phatr) to check I/O working properly
       !----------------------------------------------------------
       !
       !  NAME='pressure_ph.'
       !  WRITE(FN2,'(A12,A5,I1)') NAME,'00000',irank
       !  write(*,*) FN2
       !  open(29,file=FN2,form='unformatted')
       !  write(29) phatr
       !  close(29)
  
        call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
        call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	!=========================================================
	!                 DATA
	!=========================================================
	IF (irank.EQ.iroot) THEN
		open(28,file='data',form='formatted')
		write(28,*) ntime
		close(28)
	END IF

	aan=1.0*ntime
	!========================================================


       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
        CALL MPI_BARRIER(icomm_grid,IERR)
        if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from  processor=',irank

return
end


!===========================================================
!    Gaussian random number generator (start)
!===========================================================
function r8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
  r8_normal_ab = a + b * x
  return
end

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
!===========================================================
!    Gaussian random number generator (end)
!===========================================================

subroutine Create_commit_fileview()
	use Output_mod

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, global_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	FILE_FLAG = 1
	return
end

subroutine Create_commit_subarray()
	use Output_mod

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes, lsizes, local_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1
	return
end

subroutine OutFileName()
        use Output_mod

        WRITE(FN2,'(A12, I7.7)') NAME,        ntime !JL
        return
end

subroutine OutFileName_long()
        use Output_mod

        WRITE(FN2,'(A18, I7.7)') NAME,        ntime !JL
        return
end
subroutine OutFileName_obs()
        use Output_mod

        WRITE(FN3,'(A16, I7.7)') NAME_obs,        ntime !JL
        return
end

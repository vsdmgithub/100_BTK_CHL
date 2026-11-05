!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   Output X-planes of the instantaneous field                   *
!     *   Each plane is (ngy,ngz):                                     * 
!     *                   No halos in order to ease post-processing.   *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!     *   This code only works when each processor makes the same      *
!     *   number of write, or file access, calls -- otherwise          *
!     *   the code doesn't work                                        *
!     ******************************************************************
!=======================================================================


module Output_OneXplane_mod
      use messenger
      use data_export
      implicit none
	include "mpif.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER			   :: filetype, memtype
	CHARACTER(LEN=19)          :: FN2
	CHARACTER(LEN=12)          :: NAME

	INTEGER, PARAMETER            :: n_Iplanes = 1
	INTEGER, DIMENSION(n_Iplanes) :: Iplane

	INTEGER	:: WRITER_COMM

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan
end module



subroutine WriteXplane(ii)
use Output_OneXplane_mod

	implicit none
	REAL	, DIMENSION( 1:ngz, 1, 1:nly ) :: u_Xplane
	REAL	, DIMENSION( 1:ngz, 1, 1:nly ) :: v_Xplane
	REAL	, DIMENSION( 1:ngz, 1, 1:nly ) :: w_Xplane

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	integer	:: is, js, ks, ii
	integer	:: ie, je, ke
	double precision, allocatable	:: OutBuffer(:,:,:)

	INTEGER	:: my_istart  , my_iend
	INTEGER	:: global_indx, loc_indx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	:: i,j,k

        character(len=4) :: char_ii

        Iplane(1) = ii
        write(char_ii,'(i4.4)') ii

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output X-planes on processes) -------
	my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx

	!------ Create WRITER_COMM -------
	do i = 1,n_Iplanes
		global_indx = Iplane(i)
		if ( global_indx.gt.my_istart  .and.  global_indx.lt.my_iend ) then
			FLAG_RPOC_WRITES = 1
		end if
	end do
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	NAME='Iplane_'//char_ii//'.'
        	call OneXplane_FileName()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	        global_indices = (/  0  , 0  , jTmin_1(jblock)-1 /)
	        if (jblock.eq. 1 ) global_indices(3) = 0

		!ALC_BUG	local_indices = (/  0  , nox1-1 , 0 /)
		!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
				local_indices = (/  0  , 0 , 0 /)					!ALC_FIX
				ks = 1;		is = nox1;	js = 1					!ALC_FIX
			        if (jblock.eq. 1 ) js = 1					!ALC_FIX

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/ ngz ,  1  , ngy /)
		lsizes = (/ ngz ,  1  , nlyb-2*noy1+1 /)
	        if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(1)   +1
	        if (jblock.eq.npy)   lsizes(3)=(nlyb+0)   -(noy1)+1
	        if ( npy  .eq. 1 )   lsizes(3)=ngy


		!ALC_BUG	!--------- MEMORY SIZE OF (u_Xplane, v_Xplane, w_Xplane) -------
		!ALC_BUG	memsizes = (/ ngz , 1 , nly /)

        	global_cnt = 0
        	offset     = 0

			allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
			ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
			memsizes = lsizes							!ALC_FIX

		!------------------------ uc ---------------------------
		do i = 1,n_Iplanes
			global_indx = Iplane(i)
			loc_indx    = imap_1(global_indx)
				if ( global_indx.gt.my_istart  .and.  global_indx.lt.my_iend ) then
				do j=1,nly
				do k=1,ngz
					u_Xplane(k,1,j) = 0.25*( uc( k-1, loc_indx, j-1)+uc( k  , loc_indx, j-1) &
								+uc( k  , loc_indx, j  )+uc( k-1, loc_indx, j  )  )
				end do
				end do

				OutBuffer = u_Xplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneXplane_fileview()
				CALL Create_commit_OneXplane_subarray() 
				CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
		
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
		end do
		!------------------------ vc ---------------------------
		do i = 1,n_Iplanes
			global_indx = Iplane(i)
			loc_indx    = imap_1(global_indx)
	
			if ( global_indx.gt.my_istart  .and.  global_indx.lt.my_iend ) then
				do j=1,nly
				do k=1,ngz
					v_Xplane(k,1,j) = 0.25*( vc( k-1, loc_indx-1, j)+vc( k  , loc_indx-1, j) &
								+vc( k  , loc_indx  , j)+vc( k-1, loc_indx  , j)  )
				end do
				end do

				OutBuffer = v_Xplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneXplane_fileview()
				CALL Create_commit_OneXplane_subarray() 
				CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
	
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
		end do
	
       		!--------------- wc ------------------
		do i = 1,n_Iplanes
			global_indx = Iplane(i)
			loc_indx    = imap_1(global_indx)
	
			if ( global_indx.gt.my_istart  .and.  global_indx.lt.my_iend ) then
				do k=1,ngz
				do j=1,nly
					w_Xplane(k,1,j)= 0.25*(  wc(k, loc_indx-1, j-1)+wc(k, loc_indx-1, j  ) &
								+wc(k, loc_indx  , j  )+wc(k, loc_indx  , j-1)  )
				end do
				end do

				OutBuffer = w_Xplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneXplane_fileview()
				CALL Create_commit_OneXplane_subarray() 
        			CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
	
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
		end do
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)								!ALC_FIX

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	if (irank.eq.iroot) PRINT *,'Snapshot (I-plane) dumped by processes'

        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end


subroutine Create_commit_OneXplane_fileview()
	use Output_OneXplane_mod

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

subroutine Create_commit_OneXplane_subarray()
	use Output_OneXplane_mod

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

subroutine OneXplane_FileName()
        use Output_OneXplane_mod

        WRITE(FN2,'(A12, I7.7)') NAME,        ntime

        return
end

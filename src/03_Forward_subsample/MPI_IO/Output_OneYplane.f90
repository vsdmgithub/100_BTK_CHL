!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   Output Y-planes of the instantaneous field                   *
!     *   Each plane is (ngx,ngz):                                     * 
!     *                   No halos in order to ease post-processing.   *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!     *   This code only works when each processor makes the same      *
!     *   number of write, or file access, calls -- otherwise          *
!     *   the code doesn't work                                        *
!     ******************************************************************
!=======================================================================


module Output_OneYplane_mod
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

!JL	INTEGER, PARAMETER            :: n_Jplanes = 10
	! INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 5 , 10 , 15 , ngy-4 , ngy-9 , ngy-14 /)
!JL	  INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 15 , 25 , 35 , 45 , 55 , 60 , 65 , 70 , 75 , 80 /)
	INTEGER, PARAMETER            :: n_Jplanes = 5
	INTEGER, DIMENSION(n_Jplanes) :: Jplane

	INTEGER	:: WRITER_COMM

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan
end module

subroutine WriteYplane(jjout1)
use Output_OneYplane_mod

	implicit none
	REAL	, DIMENSION( 1:ngz, 1:nlx, 1 ) :: u_Yplane
	REAL	, DIMENSION( 1:ngz, 1:nlx, 1 ) :: v_Yplane
	REAL	, DIMENSION( 1:ngz, 1:nlx, 1 ) :: w_Yplane

        integer :: jjout1

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	integer	:: is, js, ks
	integer	:: ie, je, ke
	double precision, allocatable	:: OutBuffer(:,:,:)

	INTEGER	:: my_jstart  , my_jend
	INTEGER	:: global_jndx, loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	:: i,j,k

        character(len=4) :: char_jj !JL

        write(char_jj,'(i4.4)') jjout1

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output Y-planes on processes) -------
	my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
		global_jndx = jjout1
		if ( global_jndx.gt.my_jstart  .and.  global_jndx.lt.my_jend ) then
			FLAG_RPOC_WRITES = 1
		end if
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	NAME='Jplane_'//char_jj//'.' !JL
        	call OneYplane_FileName()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , iTmin_1(iblock)-1 , 0 /)
		if (iblock.eq. 1 ) global_indices(2) = 0

		!ALC_BUG	local_indices = (/  0  , nox1-1 , 0 /)
		!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
				local_indices = (/  0  , 0 , 0 /)					!ALC_FIX
				ks = 1;		is = nox1;	js = 1					!ALC_FIX
				if (iblock.eq. 1 ) is = 1						!ALC_FIX

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/ ngz , ngx ,  1  /)
		lsizes = (/ ngz , nlxb-2*nox1+1 ,  1  /)
		if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(1)   +1
		if (iblock.eq.npx)   lsizes(2)=(nlxb+0)   -(nox1)+1
		if ( npx  .eq. 1 )   lsizes(2)=ngx

		!ALC_BUG	!--------- MEMORY SIZE OF (u_Yplane, v_Yplane, w_Yplane) -------
		!ALC_BUG	memsizes = (/ ngz , nlx ,  1  /)

        	global_cnt = 0
        	offset     = 0

				allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
				ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
				memsizes = lsizes							!ALC_FIX

		!------------------------ uc ---------------------------
			global_jndx = jjout1
			loc_jndx    = jmap_1(global_jndx)
			if ( global_jndx.gt.my_jstart  .and.  global_jndx.lt.my_jend ) then
				do i=1,nlx
				do k=1,ngz
					u_Yplane(k,i,1) = 0.25*( uc( k-1, i, loc_jndx-1)+uc( k  , i, loc_jndx-1) &
								+uc( k  , i, loc_jndx  )+uc( k-1, i, loc_jndx  )  )
				end do
				end do

				OutBuffer = u_Yplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneYplane_fileview()
				CALL Create_commit_OneYplane_subarray() 
				CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
		
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
		!------------------------ vc ---------------------------
			global_jndx = jjout1
			loc_jndx    = jmap_1(global_jndx)
	
			if ( global_jndx.gt.my_jstart  .and.  global_jndx.lt.my_jend ) then
				do i=1,nlx
				do k=1,ngz
					v_Yplane(k,i,1) = 0.25*( vc( k-1, i-1, loc_jndx)+vc( k  , i-1, loc_jndx) &
								+vc( k  , i  , loc_jndx)+vc( k-1, i  , loc_jndx)  )
				end do
				end do

				OutBuffer = v_Yplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneYplane_fileview()
				CALL Create_commit_OneYplane_subarray() 
				CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
	
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
	
       		!--------------- wc ------------------
			global_jndx = jjout1
			loc_jndx    = jmap_1(global_jndx)
	
			if ( global_jndx.gt.my_jstart  .and.  global_jndx.lt.my_jend ) then
				do i=1,nlx
				do k=1,ngz
					w_Yplane(k,i,1)= 0.25*(  wc(k, i-1, loc_jndx-1)+wc(k, i-1, loc_jndx  ) &
								+wc(k, i  , loc_jndx  )+wc(k, i  , loc_jndx-1)  )
				end do
				end do

				OutBuffer = w_Yplane(ks:ke, is:ie, js:je)				!ALC_FIX
	
				CALL Create_commit_OneYplane_fileview()
				CALL Create_commit_OneYplane_subarray() 
        			CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
	
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)								!ALC_FIX

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	if (irank.eq.iroot) PRINT *,'Snapshot (J-plane) dumped by processes'
	FLAG_RPOC_WRITES = 0
        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end

subroutine WriteYplane_pcc(jjout1)
use Output_OneYplane_mod
! write pressure at cell center 
	implicit none
	!REAL	, DIMENSION( 1:ngz-1, 0:nlx, 1 ) :: p_Yplane

        integer :: jjout1
	integer :: ii,jj,kk

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	integer	:: is, js, ks
	integer	:: ie, je, ke
	double precision, allocatable	:: OutBuffer(:,:,:)

	INTEGER	:: my_jstart  , my_jend
	INTEGER	:: global_jndx, loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	:: i,j,k

        character(len=4) :: char_jj !JL

        write(char_jj,'(i4.4)') jjout1

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output Y-planes on processes) -------
	! Modified by MZW
	my_jstart = jTmin_1(jblock)
	my_jend   = jTmax_1(jblock)

	!------ Create WRITER_COMM -------
		global_jndx = jjout1
		loc_jndx    = jmap_1(global_jndx)
		if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
			FLAG_RPOC_WRITES = 1
		end if
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then

		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	NAME='PreTmp_'//char_jj//'.' !MZW
        	call OneYplane_FileName()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , iTmin_1(iblock), 0 /) !iTmin_1(iblock)-1 ???
		if (iblock.eq. 1 ) global_indices(2) = 0

		!ALC_BUG	local_indices = (/  0  , nox1-1 , 0 /)
		!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
				local_indices = (/  0  , 0 , 0 /)					!ALC_FIX
				ks = 0;		is = nox1;	js = loc_jndx					!ALC_FIX
				if (iblock.eq. 1 ) is = 0						!ALC_FIX

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/ ngz+1 , ngx+1 ,  1  /)
		lsizes = (/ ngz+1 , nlxb-2*nox1+1 ,  1  /)  ! nlxb-2*nox1+1 ??
		if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(0)   +1
		if (iblock.eq.npx)   lsizes(2)=(nlxb+0)   -(nox1)+1
		if ( npx  .eq. 1 )   lsizes(2)=ngx+1

		!ALC_BUG	!--------- MEMORY SIZE OF (u_Yplane, v_Yplane, w_Yplane) -------
		!ALC_BUG	memsizes = (/ ngz , nlx ,  1  /)

        	global_cnt = 0
        	offset     = 0

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
		memsizes = lsizes							!ALC_FIX

		!------------------------ p ---------------------------

			if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
				OutBuffer = reshape(p(ks:ke,is:ie,js:je),(/lsizes(1),lsizes(2),lsizes(3)/))	!ALC_FIX

				CALL Create_commit_OneYplane_fileview()
				CALL Create_commit_OneYplane_subarray() 
				CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
			end if
		
        		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        		offset = global_cnt * FloatSize
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)								!ALC_FIX

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
	FLAG_RPOC_WRITES = 0
       	if (irank.eq.iroot) PRINT *,'Snapshot (J-plane) dumped by processes'
        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end

subroutine Create_commit_OneYplane_fileview()
	use Output_OneYplane_mod

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

subroutine Create_commit_OneYplane_subarray()
	use Output_OneYplane_mod

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

subroutine OneYplane_FileName()
        use Output_OneYplane_mod

        WRITE(FN2,'(A12, I7.7)') NAME,        ntime

        return
end

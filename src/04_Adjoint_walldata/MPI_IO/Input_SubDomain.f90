! Input Subdomain
! MZW 07/23/17
! only for reading cell-centered value




module Input_SubDomain_mod
      use messenger
      use data_export
      implicit none
	include "mpif.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER*8                  :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER			   :: filetype, memtype
	CHARACTER(LEN=26)          :: FN2
	CHARACTER(LEN=18)          :: NAME

	INTEGER, PARAMETER            :: n_Jplanes = 10
	! INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 5 , 10 , 15 , ngy-4 , ngy-9 , ngy-14 /)
	  INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 15 , 25 , 35 , 45 , 55 , 60 , 65 , 70 , 75 , 80 /)

	INTEGER	:: WRITER_COMM

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0
	integer :: ntime_name

	REAL aan
end module

! MZW: read subdomain of fluxes at face center
subroutine ReadSubDomain_ffc(istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip,ntime_obs)
use Input_SubDomain_mod
!use simulation, only: Pr !JL

	implicit none
	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=18)          :: FNAME
	integer	:: istart, iend, iskip
	integer	:: jstart, jend, jskip
	integer	:: kstart, kend, kskip
	double precision, allocatable	:: InBuffer(:,:,:)

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k
	INTEGER	:: ii,jj,kk,ntime_obs

	!NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

        ntime_name = ntime_obs
	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	!my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	!my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx
	!my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	!my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

        ! MZW: this part is different from Write_Subdomain()
	my_istart = iTmin_1(iblock)	!;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	!;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	!;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jTmax_1(jblock)	!;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	!call SubDomain_FileName()
                !
        	!call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                !           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                !           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
				!loc_istart	  = imap_1( my_istart + iskip - mod(my_istart-istart,iskip) )
				 loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
				!loc_jstart	  = jmap_1( my_jstart + jskip - mod(my_jstart-jstart,jskip) )
				 loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)

                ! if skin friction is observed, it is equivalant to observe the
                ! first grid point from the wall
                if(ifric == 1) then
                        if(jblock == 1 .or. jblock == npy) lsizes(3) = lsizes(3) + 1 
                        gsizes(3) = gsizes(3) + 2
                        if(my_jstart > jstart) then
			        global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
                                global_indices(3) = global_indices(3) + 1
                        end if
                end if
		memsizes = lsizes
		!------------------------ u ---------------------------
                NAME='./obs/04uobs.dble.'
        	call SubDomain_inFileName_long()
                if(irank==1) write(*,*) 'reading:',FN2
        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))
        	global_cnt = 0
        	offset     = 0

		CALL Create_commit_SubDomain_inputview()
		CALL Create_commit_SubDomain_in_subarray() 
		CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
                CALL MPI_FILE_CLOSE(fh,ierr)

                if(ifric == 1 .and. jblock == 1) then
                        j = jmap_1(1)
                        jj = 1
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                uob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do

		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                uob(k,i,j) = InBuffer(kk,ii,jj+1)

		        end do
		        end do
		        end do
                end if
                if(ifric == 1 .and. jblock == npy) then
		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                uob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
		        end do
                        
                        j = jmap_1(ngy-1)
                        jj = lsizes(3)
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                uob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
                end if
                if(ifric == 0 .or. (ifric == 1 .and. jblock /= 1 .and. jblock /= npy)) then
		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                uob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
		        end do
                end if
                deallocate(InBuffer) 
        
                !global_cnt = int8(global_cnt) + int8(gsizes(1))*int8(gsizes(2))*int8(gsizes(3))
                !offset = int8(global_cnt) * int8(FloatSize)

		!------------------------ v ---------------------------
                if(my_jstart > jstart) then
		        global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
                end if
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes

                NAME='./obs/05vobs.dble.'
        	call SubDomain_inFileName_long()
                if(irank==1) write(*,*) 'reading:',FN2
                
        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))

		CALL Create_commit_SubDomain_inputview()
		CALL Create_commit_SubDomain_in_subarray() 
		CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		CALL MPI_FILE_CLOSE(fh,ierr)

		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

                        vob(k,i,j) = InBuffer(kk,ii,jj)

		end do
		end do
		end do
                deallocate(InBuffer)

                !global_cnt = int8(global_cnt) + int8(gsizes(1))*int8(gsizes(2))*int8(gsizes(3))
                !offset = int8(global_cnt) * int8(FloatSize)

		!------------------------ wc ---------------------------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
                ! if skin friction is observed, it is equivalant to observe the
                ! first grid point from the wall
                if(ifric == 1) then
                        if(jblock == 1 .or. jblock == npy) lsizes(3) = lsizes(3) + 1 
                        gsizes(3) = gsizes(3) + 2
                        if(my_jstart > jstart) then
			        global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
                                global_indices(3) = global_indices(3) + 1
                        end if
                end if
		memsizes = lsizes

                NAME='./obs/06wobs.dble.'
                call SubDomain_inFileName_long()
                if(irank==1) write(*,*) 'reading:',FN2

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))

		CALL Create_commit_SubDomain_inputview()
		CALL Create_commit_SubDomain_in_subarray() 
		CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		CALL MPI_FILE_CLOSE(fh,ierr)

                if(ifric == 1 .and. jblock == 1) then
                        j = jmap_1(1)
                        jj = 1
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                wob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do

		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                wob(k,i,j) = InBuffer(kk,ii,jj+1)

		        end do
		        end do
		        end do
                end if
                if(ifric == 1 .and. jblock == npy) then
		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                wob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
		        end do
                        
                        j = jmap_1(ngy-1)
                        jj = lsizes(3)
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                wob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
                end if
                if(ifric == 0 .or. (ifric == 1 .and. jblock /= 1 .and. jblock /= npy)) then
		        do j = loc_jstart, loc_jend, jskip
		        do i = loc_istart, loc_iend, iskip
		        do k =     kstart,     kend, kskip
		        	jj = (j-loc_jstart)/jskip + 1
		        	ii = (i-loc_istart)/iskip + 1
		        	kk = (k-    kstart)/kskip + 1

                                wob(k,i,j) = InBuffer(kk,ii,jj)

		        end do
		        end do
		        end do
                end if

		deallocate(InBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	if (irank.eq.iroot) PRINT *,'Snapshot (3D) dumped by processes'

        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end

subroutine ReadSubDomain_pcc(FNAME, istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip,ntime_obs)
use Input_SubDomain_mod
!use simulation, only: Pr !JL

	implicit none
	!----------------------------------------------------
	! (InBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=18)          :: FNAME
	integer	:: istart, iend, iskip
	integer	:: jstart, jend, jskip
	integer	:: kstart, kend, kskip
	double precision, allocatable	:: InBuffer(:,:,:)

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k
	INTEGER	:: ii,jj,kk
	integer :: ntime_obs

	ntime_name = ntime_obs

	NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	!my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	!my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx
	!my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	!my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

	my_istart = iTmin_1(iblock)	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jTmax_1(jblock)	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	call SubDomain_InFileName_long()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)
                write(*,*) 'reading ',FN2
		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
				!loc_istart	  = imap_1( my_istart + iskip - mod(my_istart-istart,iskip) )
				 loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
				!loc_jstart	  = jmap_1( my_jstart + jskip - mod(my_jstart-jstart,jskip) )
				 loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))

        	global_cnt = 0
        	offset     = 0

                CALL Create_commit_SubDomain_inputview()
                CALL Create_commit_SubDomain_in_subarray()
                CALL MPI_FILE_READ(fh, InBuffer, 1, memtype, status, ierr)

                !------------------------ P ---------------------------
                do j = loc_jstart, loc_jend, jskip
                do i = loc_istart, loc_iend, iskip
                do k =     kstart,     kend, kskip
                        jj = (j-loc_jstart)/jskip + 1
                        ii = (i-loc_istart)/iskip + 1
                        kk = (k-    kstart)/kskip + 1

                        p_obs( k  , i  , 2  ) = InBuffer(kk,ii,jj) * dt		! only for jend = jstart
                end do
                end do
                end do
                if(irank == 1) then
                        write(*,*) 'In readsubdomain, InBuffer=',maxval(InBuffer),dt
                        write(*,*) 'In readsubdomain,p_obs=',maxval(p_obs(1:ngz-1,1:16,2))
                end if
                global_cnt = int8(global_cnt) + int8(gsizes(1))*int8(gsizes(2))*int8(gsizes(3))
                offset = int8(global_cnt) * int8(FloatSize)

		!------ close file and clean up -----
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(InBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	!if (irank.eq.iroot) PRINT *,'Snapshot (3D) dumped by processes'
        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end

subroutine ReadSubDomain_dudy(FNAME, istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip,ntime_obs)
use Input_SubDomain_mod
!use simulation, only: Pr !JL

	implicit none
	!----------------------------------------------------
	! (InBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=18)          :: FNAME
	integer	:: istart, iend, iskip
	integer	:: jstart, jend, jskip
	integer	:: kstart, kend, kskip
	double precision, allocatable	:: InBuffer(:,:,:)

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k, jout
	INTEGER	:: ii,jj,kk, ntime_obs

	ntime_name = ntime_obs

	NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	!my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	!my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx
	!my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	!my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

	my_istart = iTmin_1(iblock)	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jTmax_1(jblock)	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	call SubDomain_InFileName_long()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
				!loc_istart	  = imap_1( my_istart + iskip - mod(my_istart-istart,iskip) )
				 loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
				!loc_jstart	  = jmap_1( my_jstart + jskip - mod(my_jstart-jstart,jskip) )
				 loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))

        	global_cnt = 0
        	offset     = 0

                CALL Create_commit_SubDomain_inputview()
                CALL Create_commit_SubDomain_in_subarray()
                CALL MPI_FILE_READ(fh, InBuffer, 1, memtype, status, ierr)

                !------------------------ dudy ---------------------------
		! only for output wall information 
		if(jstart .eq. 1) then
			jout = 1
		elseif(jstart .eq. ngy-1) then
			jout = 2
		else
			write(*,*) 'Read_dudy: wrong input of jstart, jend!!!'
			stop
		end if

                j = loc_jstart
                do i = loc_istart, loc_iend, iskip
                do k =     kstart,     kend, kskip
                        jj = (j-loc_jstart)/jskip + 1
                        ii = (i-loc_istart)/iskip + 1
                        kk = (k-    kstart)/kskip + 1

                        dudy_obs( k  , i  , jout  ) = InBuffer(kk,ii,jj) 		! only for jend = jstart
                end do
                end do

                global_cnt = int8(global_cnt) + int8(gsizes(1))*int8(gsizes(2))*int8(gsizes(3))
                offset = int8(global_cnt) * int8(FloatSize)

		!------ close file and clean up -----
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(InBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	!if (irank.eq.iroot) PRINT *,'Snapshot (3D) dumped by processes'

        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end


subroutine ReadSubDomain_dwdy(FNAME, istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip, ntime_obs)
use Input_SubDomain_mod
!use simulation, only: Pr !JL

	implicit none
	!----------------------------------------------------
	! (InBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=18)          :: FNAME
	integer	:: istart, iend, iskip
	integer	:: jstart, jend, jskip
	integer	:: kstart, kend, kskip
	double precision, allocatable	:: InBuffer(:,:,:)

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k, jout
	INTEGER	:: ii,jj,kk,ntime_obs

	ntime_name = ntime_obs
	NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	!my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	!my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx
	!my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	!my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

	my_istart = iTmin_1(iblock)	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jTmax_1(jblock)	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
        	call SubDomain_InFileName_long()

        	call MPI_FILE_OPEN(WRITER_COMM, FN2, &
                           	MPI_MODE_RDONLY, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
				!loc_istart	  = imap_1( my_istart + iskip - mod(my_istart-istart,iskip) )
				 loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
				!loc_jstart	  = jmap_1( my_jstart + jskip - mod(my_jstart-jstart,jskip) )
				 loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))

        	global_cnt = 0
        	offset     = 0

                CALL Create_commit_SubDomain_inputview()
                CALL Create_commit_SubDomain_in_subarray()
                CALL MPI_FILE_READ(fh, InBuffer, 1, memtype, status, ierr)

                !------------------------ dwdy ---------------------------
		! only for output wall information 
		if(jstart .eq. 1) then
			jout = 1
		elseif(jstart .eq. ngy-1) then
			jout = 2
		else
			write(*,*) 'Read_dwdy: wrong input of jstart, jend!!!'
			stop
		end if

                j = loc_jstart
                do i = loc_istart, loc_iend, iskip
                do k =     kstart,     kend, kskip
                        jj = (j-loc_jstart)/jskip + 1
                        ii = (i-loc_istart)/iskip + 1
                        kk = (k-    kstart)/kskip + 1

                        dwdy_obs( k  , i  , jout  ) = InBuffer(kk,ii,jj) 	! only for jend = jstart
                end do
                end do


                global_cnt = int8(global_cnt) + int8(gsizes(1))*int8(gsizes(2))*int8(gsizes(3))
                offset = int8(global_cnt) * int8(FloatSize)

		!------ close file and clean up -----
	
        	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(InBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	!if (irank.eq.iroot) PRINT *,'Snapshot (3D) dumped by processes'

        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end


subroutine Create_commit_SubDomain_inputview()
	use Input_SubDomain_mod

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

subroutine Create_commit_SubDomain_in_subarray()
	use Input_SubDomain_mod

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

subroutine SubDomain_InFileName_long()
        use Input_SubDomain_mod

        WRITE(FN2,'(A18, I7.7)') NAME,        ntime_name
	!write(*,*) 'reading ',FN2
        return
end
subroutine SubDomain_InFileName()
        use Input_SubDomain_mod

        WRITE(FN2,'(A12, I7.7)') NAME,        ntime_name
	!write(*,*) 'reading ',FN2
        return
end

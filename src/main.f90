!=======================================================================
! Simulation setup
!
!=======================================================================
program setup

	! Get parameters for setup
	call archive_init()
	iunit = iopen()
	open (iunit, file="input.setup", status="old", iostat=ierr)
	if (ierr .ne. 0) stop "An input file is required: input.setup"
	call archive_read(iunit, 0)
	close (iclose(iunit))

	!-------------------------------------------------------------------
	! Define the mesh	(Now includes reading in "grid.data")
	!-------------------------------------------------------------------
	call mesh_define()

	!-------------------------------------------------------------------
	! Initialize flow field
	!-------------------------------------------------------------------
	call initialField_define()

	!-------------------------------------------------------------------
	! Write data file
	!-------------------------------------------------------------------
	iunit = iopen()
	open (iunit, file="ucvcwc.dble.0000000", form="unformatted")
	call initialField_write(iunit)
	close (iclose(iunit))

	! Write data file
	iunit = iopen()
	open (iunit, file="uuvvww.dble.0000000", form="unformatted")
	call initial_Flux_write(iunit)
	close (iclose(iunit))

	! Write data file
	iunit = iopen()
	open (iunit, file="conold.dble.0000000", form="unformatted")
	call initial_Con_write(iunit)
	close (iclose(iunit))

	! Write data file
	iunit = iopen()
	open (iunit, file="pres_p.dble.0000000", form="unformatted")
	call initial_phatr_write(iunit)
	close (iclose(iunit))

	!-------------------------------------------------------------------
	! Clear archive and write mesh to archive
	!-------------------------------------------------------------------
	call archive_clear()
	call writeInt("irestart", 0)
	!---------------
	!open(100,file="mesh_output.dat")
	call mesh_write(0,0) !call mesh_write(100,1) call mesh_write(0, 0)
	!---------------
	!call txtr_write(0, 0)

	! Write archive file
	iunit = iopen()
	open (iunit, file="archive.0000000", form="unformatted")
	call archive_write(iunit, 1)
	close (iclose(iunit))

	! Output report
	ireport = iopen()
	open (ireport, file="report.0000000")
	call archive_info(ireport)
	call mesh_info(ireport)
	call initialField_info(ireport)
	close (iclose(ireport))

	write(*,*) 'MZW: check iMeshY!!!!!'
	stop "Exited normally"
end


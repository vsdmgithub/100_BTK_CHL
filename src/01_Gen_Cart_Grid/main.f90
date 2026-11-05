!
!  Tamer Zaki
!  Mechanical Engineering Department
!  Stanford University
!  tzaki@stanford.edu
!
!	MAC	:	setenv ABSOFT_RT_FLAGS -littleendian	(changes it to little_endian)
!			f90 -N113 -YEXT_NAMES=LCS -o Gen_grid.data.exe -s main.f90 mesh_tanh_stretch.f90
!
!	COMPAQ	:	f95    -r8       -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Alpha)
!	LAMBDA	:	lf95 --dbl       -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Lahey)
!	ALC	:	ifc -w95 -cm -r8 -o Gen_grid.data.exe main.f90 mesh_tanh_stretch.f90	(Intel)
!	CX1	:	ifc -w95 -cm -r8 -o Gen_grid.data.exe main.f90 mesh_tanh_stretch.f90	(Intel)
!


program main

	implicit none

	integer	:: i,j
	integer	:: iunit, ilength, ierr

	!--------INPUT PARAMETERS------
	real	:: Lx, Ly
	integer	:: ngx, ngy, imeshY
	integer	:: ngx2,ngy2
	real	:: ratio

	!------ GEOMETRY VARIABLES -------
	real :: y_lower, y_upper
	real :: dx, dy
	real, allocatable, dimension( : ) ::  x ,  y
	real, allocatable, dimension(:,:) :: xpg, ypg
	real, allocatable, dimension(:,:) :: xpg2, ypg2

	open(17, file="input.file")
	read(17, *) Lx
	read(17, *) Ly
	read(17, *) ngx
	read(17, *) ngy
	read(17, *) imeshY
	read(17, *) ratio
	close(17)

        ngy = (ngy-1)/2 + 1

  	allocate(  x (  ngx  ) ,  y (  ngy  ) )
  	allocate( xpg ( ngx  , ngy  ) , ypg ( ngx  , ngy  ) )

	!----- X grid ------
	dx = Lx  / (ngx-1.)
	do i=1,ngx
		x(i) = (i-1.) * dx
		xpg(i,:) = x(i)
	end do

	!----- Y grid ------
	y_lower = 0.;
	y_upper = Ly;
	if (imeshY.eq.1) then
        	call tanh_stretch(y_lower,y_upper,ngy,ratio, y)
		do i=1,ngx
			ypg(i,:) = y(:)
		end do
	else
		dy = Ly  / (ngy-1.)
		do j=1,ngy
			y(j) = (j-1.) * dy
			ypg(:,j) = y(j)
		end do
	end if

	ngx2 = ngx
	ngy2 = 2*(ngy-1)+1
  	allocate( xpg2( ngx2 , ngy2 ) , ypg2( ngx2 , ngy2 ) )

	do j=1,ngy
		ypg2(:,j) = ypg(:,j)
	end do
	do j=ngy,ngy2
		ypg2(:, j ) = Ly + ypg(:,ngy) - ypg(:,ngy2-j+1)
	end do
	do j=1,ngy2
		xpg2(:, j ) = xpg(:,1)
	end do

	!----print out y-grid to screen----
!	print*,'Print out ypg2(j=1,ngy2)'
!	do j=1,ngy2
!		print*,'j=',j, ypg2(1,j)
!	end do
!------------generate a tecplot-formatted file to visualize grid-------
!------------Mengze Sep 19 2016 -----------------------------
	open(100,file="grid_tec.dat")
	write(100,*) 'VARIABLES = "X", "Y"'
	write(100,*) 'ZONE I=', ngx2, ', J=', ngy2, ', F=POINT'
	do j=1,ngy2
		do i=1, ngx2
			write(100,*) xpg2(i,j), ypg2(i,j)
		end do
	end do
	close(100)
	open(101,file="grid_stretch_tec.dat")
	write(101,*) 'VARIABLES = "n", "Y", "dY", "dY_n+1/dY_n"'
	do j=1,ngy2-2
		write(101,*) j, ypg2(1,j), ypg2(1,j+1)-ypg2(1,j), (ypg2(1,j+2)-ypg2(1,j+1))/(ypg2(1,j+1)-ypg2(1,j))
	end do
	close(101)

!---------------------------------------------------------------------
	inquire(iolength=ilength) xpg2
	iunit = 19
	open(iunit, file="grid.data", form="unformatted", access="direct", &
		recl=ilength, iostat=ierr)
	if (ierr .ne. 0) stop "grid.data file can not be created"
	write(iunit, rec=1) xpg2
	write(iunit, rec=2) ypg2
	close(iunit)

end program

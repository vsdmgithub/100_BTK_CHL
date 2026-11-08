! MZW 07/13/17
! IC-only version for channel flow: periodic in x and z, no-slip at top/bottom walls.
! Control vector holds the independent interior velocities only:
!   u(1:ngz-1, 2:ngx  , 1:ngy-1)   (u(:,1,:) is the periodic copy of u(:,ngx,:))
!   v(1:ngz-1, 1:ngx-1, 2:ngy-1)   (v = 0 on the walls j=1, ngy)
!   w(1:ngz-1, 1:ngx-1, 1:ngy-1)   (w(ngz,:,:) is the periodic copy of w(1,:,:))


subroutine Initial_guess(dimv,x)
! read initial flux field
	include "param.inc"
	integer,intent(in) :: dimv
	real(kind=8),intent(out) :: x(dimv)
	real(kind=8) :: u(0:ngz,0:ngx+1,0:ngy), v(0:ngz,0:ngx,0:ngy+1), w(0:ngz+1,0:ngx,0:ngy)
	integer :: Instart, Inend
	integer :: ilength

	inquire(iolength=ilength) u
	open(10,file="01uc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "01uc00 file can not be read"
	read(10,rec=1) u
	write(*,*) 'reading 01uc00.dble.0000000'
	close(10)

	inquire(iolength=ilength) v
	open(11,file="02vc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "02vc00 file can not be read"
	read(11,rec=1) v
	write(*,*) 'reading 02vc00.dble.0000000'
	close(11)

	inquire(iolength=ilength) w
	open(12,file="03wc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "03wc00 file can not be read"
	read(12,rec=1) w
	write(*,*) 'reading 03wc00.dble.0000000'
	close(12)

	Instart          = 1
	Inend            = (ngx-1)*(ngy-1)*(ngz-1)
	x(Instart:Inend) = reshape(u(1:ngz-1,2:ngx,1:ngy-1),(/(ngx-1)*(ngy-1)*(ngz-1)/))
	Instart          = Inend+1
	Inend            = Inend+(ngx-1)*(ngy-2)*(ngz-1)
	x(Instart:Inend) = reshape(v(1:ngz-1,1:ngx-1,2:ngy-1),(/(ngx-1)*(ngy-2)*(ngz-1)/))
	Instart          = Inend+1
	Inend            = Inend+(ngx-1)*(ngy-1)*(ngz-1)
	x(Instart:Inend) = reshape(w(1:ngz-1,1:ngx-1,1:ngy-1),(/(ngx-1)*(ngy-1)*(ngz-1)/))

end

subroutine Evaluate_f(dimv,x,f)
! run forward simulation and evalute cost function
	include "param.inc"
	integer,intent(in) :: dimv
	real(kind=8),intent(in) :: x(dimv)
	real(kind=8),intent(out) :: f
	real(kind=8) :: u(0:ngz,0:ngx+1,0:ngy), v(0:ngz,0:ngx,0:ngy+1), w(0:ngz+1,0:ngx,0:ngy)
	integer :: Instart, Inend
	integer :: ilength
	character(LEN=96) :: command_line

	u=0.0; v=0.0; w=0.0
	Instart                    = 1
	Inend                      = (ngx-1)*(ngy-1)*(ngz-1)
	u(1:ngz-1,2:ngx,1:ngy-1)   = reshape(x(Instart:Inend),(/ngz-1,ngx-1,ngy-1/))
	Instart                    = Inend+1
	Inend                      = Inend+(ngx-1)*(ngy-2)*(ngz-1)
	v(1:ngz-1,1:ngx-1,2:ngy-1) = reshape(x(Instart:Inend),(/ngz-1,ngx-1,ngy-2/))
	Instart                    = Inend+1
	Inend                      = Inend+(ngx-1)*(ngy-1)*(ngz-1)
	w(1:ngz-1,1:ngx-1,1:ngy-1) = reshape(x(Instart:Inend),(/ngz-1,ngx-1,ngy-1/))

	! fill periodic copies (same as FluxBC in the forward code), so the
	! written IC file is self-consistent; the forward run resets them anyway
	! x-direction
	u(:,1    ,:) = u(:,ngx  ,:)
	u(:,0    ,:) = u(:,ngx-1,:)
	u(:,ngx+1,:) = u(:,2    ,:)
	v(:,0    ,:) = v(:,ngx-1,:)
	v(:,ngx  ,:) = v(:,1    ,:)
	w(:,0    ,:) = w(:,ngx-1,:)
	w(:,ngx  ,:) = w(:,1    ,:)
	! z-direction
	u(0    ,:,:) = u(ngz-1,:,:)
	u(ngz  ,:,:) = u(1    ,:,:)
	v(0    ,:,:) = v(ngz-1,:,:)
	v(ngz  ,:,:) = v(1    ,:,:)
	w(0    ,:,:) = w(ngz-1,:,:)
	w(ngz  ,:,:) = w(1    ,:,:)
	w(ngz+1,:,:) = w(2    ,:,:)

	inquire(iolength=ilength) u
	open(10,file="01uc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength,iostat=ierr)
	if (ierr .ne. 0) stop "01uc00 file can not be created"
	write(10,rec=1) u
	write(*,*) '01uc00.dble.0000000'
	close(10)

	inquire(iolength=ilength) v
	open(11,file="02vc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength,iostat=ierr)
	if (ierr .ne. 0) stop "02vc00 file can not be created"
	write(11,rec=1) v
	write(*,*) '02vc00.dble.0000000'
	close(11)

	inquire(iolength=ilength) w
	open(12,file="03wc00.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, iostat=ierr)
	if (ierr .ne. 0) stop "03wc00 file can not be created"
	write(12,rec=1) w
	write(*,*) '03wc00.dble.0000000'
	close(12)

	! run forward simulation
	call execute_command_line("cp data_zero data")
	call execute_command_line("cp archive.0000000 archive")
	! call execute_command_line("rm -f for/*")

	write(*,*) '---------------start forward simulation------------'
	! nproc is in param.inc
	write(command_line,'(A,I3,A)') 'mpirun -np',nproc,' ./forward.out > history_forward'
	CALL execute_command_line(trim(command_line)) 
	write(*,*) '---------------end forward simulation------------'

	open(100,file="cost_function")
	read(100,*) f
	write(*,*) 'cost function=',f
	close(100)		! very important!!
end

subroutine Evaluate_g(dimv,g)
! run adjoint simulation and read gradient w.r.t. initial condition
	include "param.inc"
	integer,intent(in) :: dimv
	real(kind=8),intent(out) :: g(dimv)
	real(kind=8) :: u_adj(0:ngz,0:ngx+1,0:ngy), v_adj(0:ngz,0:ngx,0:ngy+1), w_adj(0:ngz+1,0:ngx,0:ngy)
	integer :: Instart, Inend
	integer :: ilength
	character(LEN=96) :: command_line


	! run adjoint
	call execute_command_line("cp data_zero data")
	call execute_command_line("cp archive.0000000 archive")

	write(*,*) '---------------start adjoint simulation------------'
	write(command_line,'(A,I3,A)') 'mpirun -np',nproc,' ./adjoint.out > history_adjoint'
	CALL execute_command_line(trim(command_line)) 
	write(*,*) '---------------end adjoint simulation------------'

	inquire(iolength=ilength) u_adj
	open(10,file="01ucad.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "01ucad file can not be read"
	read(10,rec=1) u_adj
	close(10)
	inquire(iolength=ilength) v_adj
	open(11,file="02vcad.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "02vcad file can not be read"
	read(11,rec=1) v_adj
	close(11)
	inquire(iolength=ilength) w_adj
	open(12,file="03wcad.dble.0000000", form="unformatted", access="direct", &
		recl=ilength, status="old", iostat=ierr)
	if (ierr .ne. 0) stop "03wcad file can not be read"
	read(12,rec=1) w_adj
	close(12)

	Instart          = 1
	Inend            = (ngx-1)*(ngy-1)*(ngz-1)
	g(Instart:Inend) = reshape(u_adj(1:ngz-1,2:ngx,1:ngy-1),(/(ngx-1)*(ngy-1)*(ngz-1)/))
	Instart          = Inend+1
	Inend            = Inend+(ngx-1)*(ngy-2)*(ngz-1)
	g(Instart:Inend) = reshape(v_adj(1:ngz-1,1:ngx-1,2:ngy-1),(/(ngx-1)*(ngy-2)*(ngz-1)/))
	Instart          = Inend+1
	Inend            = Inend+(ngx-1)*(ngy-1)*(ngz-1)
	g(Instart:Inend) = reshape(w_adj(1:ngz-1,1:ngx-1,1:ngy-1),(/(ngx-1)*(ngy-1)*(ngz-1)/))

	write(*,*) 'gradient',norm2(g)
end

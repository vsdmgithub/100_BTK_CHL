!234567
subroutine setup
  use msetup
  use mvariable
  use mgeomdata, only: Lz
  implicit none
  integer(int_p):: i
  character(70):: dummy

  !==== READ INPUT FILE ========
  open(2,file='input_postproc',status='old',action='read')

  do i = 1 , 4
    read (2,300) dummy
    write(*,300) dummy
  end do

  read (2,304) dummy, cfn      ; write(*,304) dummy, cfn        ! data directory
  read (2,301) dummy, inssn    ; write(*,301) dummy, inssn      ! first snapshot
  read (2,301) dummy, insen    ; write(*,301) dummy, insen      ! last snapshot
  read (2,301) dummy, insitv   ; write(*,301) dummy, insitv     ! snapshot interval
  read (2,302) dummy, Lz       ; write(*,302) dummy, Lz         ! spanwise length
  read (2,302) dummy, Re       ; write(*,302) dummy, Re         ! Reynolds number
  read (2,301) dummy, wgrid    ; write(*,301) dummy, wgrid      ! grid: 0 none, 1 compact, 2 full 3D
  read (2,301) dummy, iwrite3d ; write(*,301) dummy, iwrite3d   ! 3D velocity toggle
  read (2,301) dummy, iwrite2d ; write(*,301) dummy, iwrite2d   ! 2D slice toggle
  read (2,301) dummy, k2d      ; write(*,301) dummy, k2d        ! z-plane of 2D slice
  read (2,301) dummy, isingle  ; write(*,301) dummy, isingle    ! all snapshots in one file
  read (2,304) dummy, outdir   ; write(*,304) dummy, outdir     ! output directory

  read (2,300) dummy
  write(*,300) dummy
  close(2)

300 format(a70)
301 format(a50,i15)
302 format(a50,es15.8)
304 format(a15,a50)

  if( k2d == 0 ) k2d = (ngz+1)/2
  if( k2d < 1 .or. k2d > ngz ) stop 'setup: 2D slice index k outside 1..ngz'

  if( isingle == 1 .and. iwrite3d == 1 .and. wgrid /= 0 ) then
    write(*,*) 'NOTE: single-file mode stores the grid inside vel.plt; wgrid -> 0'
    wgrid = 0
  end if

  if( iwrite3d == 1 .and. wgrid == 1 ) then
    write(*,*) 'WARNING: 3D velocity files need the full 3D grid; wgrid 1 -> 2'
    wgrid = 2
  end if

  call mesh

  allocate( u(ngx,ngy,ngz), v(ngx,ngy,ngz), w(ngx,ngy,ngz) )
  allocate( ut(0:ngz  ,0:ngx+1,0:ngy  ) )
  allocate( vt(0:ngz  ,0:ngx  ,0:ngy+1) )
  allocate( wt(0:ngz+1,0:ngx  ,0:ngy  ) )
  if( isingle == 1 .and. iwrite2d == 1 ) then
    allocate( s2d1(ngx,ngy,(insen-inssn)/insitv+1), s2d2(ngx,ngy,(insen-inssn)/insitv+1) )
  end if

  return
end subroutine setup

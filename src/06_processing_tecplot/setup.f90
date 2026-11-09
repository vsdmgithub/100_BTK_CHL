!234567
subroutine setup
  use msetup
  use mvariable
  use mgeomdata, only: Lz
  use mplanes, only: init_planes
  implicit none
  integer(int_p):: i, d
  character(70):: dummy
  character(200):: pline(3)
  integer(int_p):: plist(64,3), pcount(3), ngd(3)

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
  do d = 1, 3                                                   ! i / j / k plane lists
    read (2,305) dummy, pline(d) ; write(*,305) dummy, trim(adjustl(pline(d)))
  end do
  read (2,301) dummy, isingle  ; write(*,301) dummy, isingle    ! all snapshots in one file
  read (2,304) dummy, outdir   ; write(*,304) dummy, outdir     ! output directory

  read (2,300) dummy
  write(*,300) dummy
  close(2)

300 format(a70)
301 format(a50,i15)
302 format(a50,es15.8)
304 format(a15,a50)
305 format(a50,a)

  ngd = [ngx, ngy, ngz]
  pcount = 0
  if( iwrite2d == 1 ) then
    do d = 1, 3
      call parse_planes(pline(d), ngd(d), plist(:,d), pcount(d))
    end do
    if( sum(pcount) == 0 ) write(*,*) 'WARNING: 2D output on, but no planes selected'
  end if

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
  call init_planes(plist, pcount)

  return

contains

  ! Plane list such as "1,64" or "1 33 64": 0 = mid-plane, negative = no plane
  subroutine parse_planes(str, ng, list, n)
    character(*), intent(in):: str
    integer(int_p), intent(in):: ng
    integer(int_p), intent(out):: list(:), n
    character(len(str)):: s
    integer(int_p):: val, p, ios

    s = str
    do p = 1, len(s)
      if( s(p:p) == ',' .or. s(p:p) == ';' ) s(p:p) = ' '
    end do
    n = 0
    s = adjustl(s)
    do while( len_trim(s) > 0 )
      p = index(s,' ')
      read(s(1:p-1),*,iostat=ios) val
      if( ios /= 0 ) stop 'setup: cannot read the 2D plane list in input_postproc'
      if( val == 0 ) val = (ng+1)/2
      if( val > ng ) stop 'setup: 2D plane index larger than the grid size'
      if( val > 0 ) then
        if( n == size(list) ) stop 'setup: too many 2D planes in one direction'
        n = n + 1
        list(n) = val
      end if
      s = adjustl(s(p:))
    end do
  end subroutine parse_planes

end subroutine setup

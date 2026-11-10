!234567
! 2D planes of the velocity field (x y z u v w), one Tecplot file per plane:
!   i-plane (YZ), j-plane (XZ), k-plane (XY), named vel_<i|j|k>_NNNN[_TTTTTTT].plt
! Per-snapshot mode: vel_k_0033_0000200.plt for every snapshot.
! Single-file mode:  slices are buffered and written after the snapshot loop as
!                    vel_k_0033.plt (one zone per snapshot, x y z shared from zone 1),
!                    because this TECIO cannot keep two .plt files open at once.
module mplanes
  use mparameter
  use msetup
  use mgeomdata
  use mtecio
  implicit none
  type tplane
    integer(int_p):: dir, idx          ! dir: 1 = i (YZ), 2 = j (XZ), 3 = k (XY)
    integer(int_p):: ni, nj, nk        ! Tecplot zone size (one of them is 1)
    integer(int_p):: n1, n2            ! in-plane size
    real(real_p), allocatable:: buf(:,:,:,:)   ! (n1,n2,u/v/w,snapshot), single-file mode
  end type tplane
  type(tplane), allocatable:: planes(:)
  character(1), parameter:: cdir(3) = ['i','j','k']

contains

  subroutine init_planes(lists, counts)
    integer(int_p), intent(in):: lists(:,:), counts(3)
    integer(int_p):: d, m, ip, nsnap
    nsnap = (insen-inssn)/insitv + 1
    allocate( planes(sum(counts)) )
    ip = 0
    do d = 1, 3
      do m = 1, counts(d)
        ip = ip + 1
        planes(ip)%dir = d
        planes(ip)%idx = lists(m,d)
        planes(ip)%ni = ngx ; planes(ip)%nj = ngy ; planes(ip)%nk = ngz
        select case(d)
        case(1); planes(ip)%ni = 1; planes(ip)%n1 = ngy; planes(ip)%n2 = ngz
        case(2); planes(ip)%nj = 1; planes(ip)%n1 = ngx; planes(ip)%n2 = ngz
        case(3); planes(ip)%nk = 1; planes(ip)%n1 = ngx; planes(ip)%n2 = ngy
        end select
        if( isingle == 1 ) allocate( planes(ip)%buf(planes(ip)%n1,planes(ip)%n2,3,nsnap) )
        write(*,'(a,a,a,i5)') '    2D plane ', cdir(d), ' =', lists(m,d)
      end do
    end do
  end subroutine init_planes

  function plane_name(p) result(s)
    type(tplane), intent(in):: p
    character(10):: s
    write(s,'(a,a1,a,i4.4)') 'vel_', cdir(p%dir), '_', p%idx
  end function plane_name

  subroutine slice(p, f, s)
    type(tplane), intent(in):: p
    real(real_p), intent(in):: f(ngx,ngy,ngz)
    real(real_p), intent(out):: s(:,:)
    select case(p%dir)
    case(1); s = f(p%idx,:,:)
    case(2); s = f(:,p%idx,:)
    case(3); s = f(:,:,p%idx)
    end select
  end subroutine slice

  subroutine coords(p, x, y, z)
    type(tplane), intent(in):: p
    real(real_p), intent(out):: x(:,:), y(:,:), z(:,:)
    integer(int_p):: k
    select case(p%dir)
    case(1)
      do k = 1, ngz
        x(:,k) = xpg(p%idx,:); y(:,k) = ypg(p%idx,:); z(:,k) = dz*real(k-1,real_p)
      end do
    case(2)
      do k = 1, ngz
        x(:,k) = xpg(:,p%idx); y(:,k) = ypg(:,p%idx); z(:,k) = dz*real(k-1,real_p)
      end do
    case(3)
      x = xpg; y = ypg; z = dz*real(p%idx-1,real_p)
    end select
  end subroutine coords

  subroutine write_xyz(p)
    type(tplane), intent(in):: p
    real(real_p), allocatable:: x(:,:), y(:,:), z(:,:)
    allocate( x(p%n1,p%n2), y(p%n1,p%n2), z(p%n1,p%n2) )
    call coords(p, x, y, z)
    call tec_data(p%n1*p%n2, x)
    call tec_data(p%n1*p%n2, y)
    call tec_data(p%n1*p%n2, z)
  end subroutine write_xyz

  ! Called once per snapshot
  subroutine write_planes(imore)
    use mvariable, only: u, v, w
    integer(int_p), intent(in):: imore
    integer(int_p):: ip, n
    character(7):: cidx
    real(real_p), allocatable:: s(:,:)

    if( size(planes) == 0 ) return
    write(cidx,'(i7.7)') imore
    n = (imore-inssn)/insitv + 1

    do ip = 1, size(planes)
      associate( p => planes(ip) )
      if( isingle == 1 ) then
        call slice(p, u, p%buf(:,:,1,n))
        call slice(p, v, p%buf(:,:,2,n))
        call slice(p, w, p%buf(:,:,3,n))
      else
        write(*,*) '    Writing ', plane_name(p)
        call tec_open('InsPlane','x y z u v w', &
             trim(outdir)//plane_name(p)//'_'//cidx//'.plt', 0)
        call tec_zone('T='//cidx, p%ni, p%nj, p%nk, 6, real(imore,8), 1)
        call write_xyz(p)
        allocate( s(p%n1,p%n2) )
        call slice(p, u, s); call tec_data(p%n1*p%n2, s)
        call slice(p, v, s); call tec_data(p%n1*p%n2, s)
        call slice(p, w, s); call tec_data(p%n1*p%n2, s)
        deallocate( s )
        call tec_close
      end if
      end associate
    end do
  end subroutine write_planes

  ! Single-file mode: called once after the snapshot loop
  subroutine write_planes_all
    integer(int_p):: ip, n, c, imore
    character(7):: cidx

    do ip = 1, size(planes)
      associate( p => planes(ip) )
      write(*,*) '    Writing ', plane_name(p), '.plt (all snapshots)'
      call tec_open('InsPlane','x y z u v w', trim(outdir)//plane_name(p)//'.plt', 0)
      do n = 1, size(p%buf,4)
        imore = inssn + (n-1)*insitv
        write(cidx,'(i7.7)') imore
        if( n == 1 ) then
          call tec_zone('T='//cidx, p%ni, p%nj, p%nk, 6, real(imore,8), 1)
          call write_xyz(p)
        else
          call tec_zone('T='//cidx, p%ni, p%nj, p%nk, 6, real(imore,8), 1, nshare=3)
        end if
        do c = 1, 3
          call tec_data(p%n1*p%n2, p%buf(:,:,c,n))
        end do
      end do
      call tec_close
      end associate
    end do
  end subroutine write_planes_all

end module mplanes

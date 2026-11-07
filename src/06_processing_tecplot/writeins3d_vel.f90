!234567
! Per-snapshot mode: solution-only vel_NNNNNNN.plt (u v w), pairs with the full 3D grid.plt.
! Single-file mode:  all snapshots as zones of one self-contained vel.plt (x y z u v w);
!                    x y z are written in zone 1 only and shared by the later zones.
!                    The file stays open across snapshots and is closed in main.
subroutine writeins3d_vel(imore)
  use msetup
  use mvariable
  use mgeomdata
  use mtecio
  implicit none
  integer(int_p),intent(in):: imore
  integer(int_p):: k
  character(7):: cidx
  character(10):: cre
  real(real_p), allocatable:: Grid(:,:,:)

  write(cidx,fmt='(i7.7)') imore
  write(cre,fmt='(f10.1)') Re
  write(*,*) '    Writing the instantaneous 3D velocity field'

  if( isingle /= 1 ) then
    call tec_open('InsField','u v w',trim(outdir)//'vel_'//cidx//'.plt',2)
    call tec_aux('Re',adjustl(cre))
    call tec_zone('FIELD',ngx,ngy,ngz,3,real(imore,8),1)
  else if( imore == inssn ) then
    call tec_open('InsField','x y z u v w',trim(outdir)//'vel.plt',0)
    call tec_aux('Re',adjustl(cre))
    call tec_zone('T='//cidx,ngx,ngy,ngz,6,real(imore,8),1)
    allocate( Grid(ngx,ngy,ngz) )
    do k = 1, ngz
      Grid(:,:,k) = xpg
    end do
    call tec_data(ngx*ngy*ngz, Grid)
    do k = 1, ngz
      Grid(:,:,k) = ypg
    end do
    call tec_data(ngx*ngy*ngz, Grid)
    do k = 1, ngz
      Grid(:,:,k) = dz * real( k - 1, real_p )
    end do
    call tec_data(ngx*ngy*ngz, Grid)
    deallocate( Grid )
  else
    call tec_zone('T='//cidx,ngx,ngy,ngz,6,real(imore,8),1,nshare=3)
  end if

  call tec_data(ngx*ngy*ngz, u)
  call tec_data(ngx*ngy*ngz, v)
  call tec_data(ngx*ngy*ngz, w)
  if( isingle /= 1 ) call tec_close

  return
end subroutine writeins3d_vel

!234567
! Interpolate the staggered fields (solver order k,i,j) to the grid nodes and
! store them in Tecplot order (i,j,k).
subroutine collocate_on_grid
  use mparameter, only: ngx, ngy, ngz, real_p
  use mvariable
  implicit none
  integer :: i, j, k

!$omp parallel do private(i,k) schedule(static)
  do j = 1, ngy
    do i = 1, ngx
      do k = 1, ngz
        u(i,j,k) = 0.25d0 * ( ut(k,i,j) + ut(k,i,j-1) + ut(k-1,i,j) + ut(k-1,i,j-1) )
        v(i,j,k) = 0.25d0 * ( vt(k,i,j) + vt(k,i-1,j) + vt(k-1,i,j) + vt(k-1,i-1,j) )
        w(i,j,k) = 0.25d0 * ( wt(k,i,j) + wt(k,i-1,j) + wt(k,i,j-1) + wt(k,i-1,j-1) )
      end do
    end do
  end do
!$omp end parallel do

  return
end subroutine collocate_on_grid

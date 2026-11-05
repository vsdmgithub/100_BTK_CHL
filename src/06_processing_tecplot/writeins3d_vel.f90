!234567
! Solution-only file (u v w); pairs with the full 3D grid.plt.
subroutine writeins3d_vel(imore)
  use msetup
  use mvariable
  use mtecio
  implicit none
  integer(int_p),intent(in):: imore
  character(7):: cidx
  character(10):: cre

  write(cidx,fmt='(i7.7)') imore
  write(cre,fmt='(f10.1)') Re
  write(*,*) '    Writing the instantaneous 3D velocity field'

  call tec_open('InsField','u v w',trim(outdir)//'vel_'//cidx//'.plt',2)
  call tec_aux('Re',adjustl(cre))
  call tec_zone('FIELD',ngx,ngy,ngz,3,real(imore,8),1)
  call tec_data(ngx*ngy*ngz, u)
  call tec_data(ngx*ngy*ngz, v)
  call tec_data(ngx*ngy*ngz, w)
  call tec_close

  return
end subroutine writeins3d_vel

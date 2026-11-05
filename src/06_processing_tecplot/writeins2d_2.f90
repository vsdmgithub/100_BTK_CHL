!234567
! Self-contained 2D x-y slice (x y + two fields) of the plane k = k2d.
subroutine writeins2d_2(dat1,dat2,varname1,varname2,imore)
  use mgeomdata
  use msetup
  use mtecio
  implicit none
  real(real_p),intent(in):: dat1(ngx,ngy),dat2(ngx,ngy)
  character(4),intent(in):: varname1,varname2
  integer(int_p),intent(in):: imore
  character(6):: cidx

  write(cidx,fmt='(i6.6)') imore

  call tec_open('InsField','x y '//varname1//' '//varname2, &
       trim(outdir)//'ins_'//varname1//'_'//varname2//'_XY_'//cidx//'.plt',0)
  call tec_zone('FIELD',ngx,ngy,1,4,real(imore,8),1)
  call tec_data(ngx*ngy, xpg)
  call tec_data(ngx*ngy, ypg)
  call tec_data(ngx*ngy, dat1)
  call tec_data(ngx*ngy, dat2)
  call tec_close

  return
end subroutine writeins2d_2

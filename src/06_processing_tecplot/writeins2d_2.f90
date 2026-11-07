!234567
! 2D x-y slice (x y + two fields) of the plane k = k2d.
! Per-snapshot mode: one ins_<v1>_<v2>_XY_NNNNNN.plt per snapshot.
! Single-file mode:  the slice is stored in s2d1/s2d2; writeins2d_2_all writes them after the
!                    loop (TECIO here cannot keep two .plt files open at the same time).
subroutine writeins2d_2(dat1,dat2,varname1,varname2,imore)
  use mgeomdata
  use msetup
  use mvariable, only: s2d1, s2d2
  use mtecio
  implicit none
  real(real_p),intent(in):: dat1(ngx,ngy),dat2(ngx,ngy)
  character(4),intent(in):: varname1,varname2
  integer(int_p),intent(in):: imore
  integer(int_p):: n
  character(6):: cidx

  if( isingle == 1 ) then
    n = (imore-inssn)/insitv + 1
    s2d1(:,:,n) = dat1
    s2d2(:,:,n) = dat2
    return
  end if

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

!===============================================================
! Single-file mode: all stored slices as zones of ins_<v1>_<v2>_XY.plt; x y only in zone 1.
!===============================================================
subroutine writeins2d_2_all(varname1,varname2)
  use mgeomdata
  use msetup
  use mvariable, only: s2d1, s2d2
  use mtecio
  implicit none
  character(4),intent(in):: varname1,varname2
  integer(int_p):: n, imore
  character(6):: cidx

  write(*,*) '    Writing all 2D slices to one file'
  call tec_open('InsField','x y '//varname1//' '//varname2, &
       trim(outdir)//'ins_'//varname1//'_'//varname2//'_XY.plt',0)

  do n = 1, size(s2d1,3)
    imore = inssn + (n-1)*insitv
    write(cidx,fmt='(i6.6)') imore
    if( n == 1 ) then
      call tec_zone('T='//cidx,ngx,ngy,1,4,real(imore,8),1)
      call tec_data(ngx*ngy, xpg)
      call tec_data(ngx*ngy, ypg)
    else
      call tec_zone('T='//cidx,ngx,ngy,1,4,real(imore,8),1,nshare=2)
    end if
    call tec_data(ngx*ngy, s2d1(:,:,n))
    call tec_data(ngx*ngy, s2d2(:,:,n))
  end do

  call tec_close
  return
end subroutine writeins2d_2_all

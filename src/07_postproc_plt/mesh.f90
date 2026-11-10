!234567
subroutine mesh
  use mgeomdata
  implicit none
  integer(int_p):: ierr,ilength

  write(*,*) 'Reading... grid.data'
  inquire(iolength=ilength) xpg
  open(19, file='./grid.data', form='unformatted', access='direct', &
       recl=ilength, status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) stop 'grid.data file is required'
  read(19, rec=1) xpg
  read(19, rec=2) ypg
  close(19)

  dz = Lz / real(ngz-1, real_p)

  return
end subroutine mesh

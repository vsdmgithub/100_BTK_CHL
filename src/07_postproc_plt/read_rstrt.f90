!234567
subroutine read_restart(imore)
  use mvariable
  use msetup
  implicit none
  integer(int_p),intent(in):: imore
  integer(int_p):: length,lengthu,lengthv,lengthw
  character(7):: cidx
  character(100):: fname

  write(cidx,fmt='(i7.7)') imore

  inquire(iolength=length) ut(1,1,1)
  lengthu = length * (ngz+1)*(ngx+2)*(ngy+1)
  lengthv = length * (ngz+1)*(ngx+1)*(ngy+2)
  lengthw = length * (ngz+2)*(ngx+1)*(ngy+1)

  fname = trim(cfn)//'01uc00.dble.'//cidx
  write(*,*) 'Read ', trim(fname)
  open (18,file=trim(fname),form="unformatted",access="direct",action='read',status='old',recl=lengthu)
  read (18,rec=1) ut
  close(18)

  fname = trim(cfn)//'02vc00.dble.'//cidx
  write(*,*) 'Read ', trim(fname)
  open (19,file=trim(fname),form="unformatted",access="direct",action='read',status='old',recl=lengthv)
  read (19,rec=1) vt
  close(19)

  fname = trim(cfn)//'03wc00.dble.'//cidx
  write(*,*) 'Read ', trim(fname)
  open (20,file=trim(fname),form="unformatted",access="direct",action='read',status='old',recl=lengthw)
  read (20,rec=1) wt
  close(20)

  call collocate_on_grid

  write(*,*) ' read done'

  return
end subroutine read_restart

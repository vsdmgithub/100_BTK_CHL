!234567
program main
  use msetup
  use mvariable
  use mtecio
  implicit none
  integer(int_p):: imore

  call setup

  if( wgrid == 1 ) call write_grid_compact
  if( wgrid == 2 ) call write_grid

  do imore = inssn, insen, insitv

    write(*,101)
    write(*,102) imore
    write(*,101)

    call read_restart(imore)
    if( iwrite2d == 1 ) call writeins2d_2(u(:,:,k2d),v(:,:,k2d),'uraw','vraw',imore)
    if( iwrite3d == 1 ) call writeins3d_vel(imore)

  end do

  if( isingle == 1 .and. iwrite3d == 1 ) call tec_close
  if( isingle == 1 .and. iwrite2d == 1 ) call writeins2d_2_all('uraw','vraw')

  write(*,*) 'END SIMULATION! GOOD LUCK!'

101 format(30('*='))
102 format('     ##### Post process for the ',i6,' snapshot #####')
end program main

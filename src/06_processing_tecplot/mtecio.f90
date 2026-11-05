!234567
! Thin wrapper around the TECIO "142" API: return codes are checked, and
! every zone is tagged with a strand id and solution time so Tecplot sees
! the snapshots as one transient data set.
module mtecio
  implicit none
  ! 0 = store single precision on disk, 1 = double precision on disk
  integer, parameter :: plt_precision = 0
  integer, external  :: tecini142, teczne142, tecdat142, tecend142
  integer, external  :: tecauxstr142, teczauxstr142
contains

  subroutine tec_check(ret, what)
    integer, intent(in)      :: ret
    character(*), intent(in) :: what
    if( ret /= 0 ) then
      write(*,*) 'TECIO error in ', what, ' (code ', ret, ')'
      stop 1
    end if
  end subroutine tec_check

  ! filetype: 0 = full, 1 = grid only, 2 = solution only
  subroutine tec_open(title, variables, fname, filetype)
    character(*), intent(in) :: title, variables, fname
    integer, intent(in)      :: filetype
    integer :: ret
    ret = tecini142( trim(title)//char(0), trim(variables)//char(0), &
                     trim(fname)//char(0), '.'//char(0),             &
                     0, filetype, 0, plt_precision )
    call tec_check(ret, 'tecini142 '//trim(fname))
  end subroutine tec_open

  ! strand = 0 for static (grid) zones, >= 1 for time-dependent zones
  subroutine tec_zone(ztitle, ni, nj, nk, nvar, time, strand)
    character(*), intent(in) :: ztitle
    integer, intent(in)      :: ni, nj, nk, nvar, strand
    real(8), intent(in)      :: time
    integer :: ret, passive(nvar), valloc(nvar), share(nvar)
    passive = 0     ! no passive variables
    valloc  = 1     ! all variables nodal
    share   = 0     ! no variable sharing
    ret = teczne142( trim(ztitle)//char(0), 0, ni, nj, nk, 0, 0, 0, time, &
                     strand, 0, 1, 0, 0, 0, 0, 0, passive, valloc, share, 0 )
    call tec_check(ret, 'teczne142 '//trim(ztitle))
  end subroutine tec_zone

  subroutine tec_data(n, arr)
    integer, intent(in) :: n
    real(8), intent(in) :: arr(*)
    integer :: ret
    ret = tecdat142( n, arr, 1 )
    call tec_check(ret, 'tecdat142')
  end subroutine tec_data

  subroutine tec_aux(name, value)
    character(*), intent(in) :: name, value
    integer :: ret
    ret = tecauxstr142( trim(name)//char(0), trim(value)//char(0) )
    call tec_check(ret, 'tecauxstr142 '//trim(name))
  end subroutine tec_aux

  subroutine tec_close
    integer :: ret
    ret = tecend142()
    call tec_check(ret, 'tecend142')
  end subroutine tec_close

end module mtecio

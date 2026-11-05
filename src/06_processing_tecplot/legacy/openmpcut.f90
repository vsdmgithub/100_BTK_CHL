! *******************************************************
! START MODULE OPENMPCUT ********************************
! *******************************************************
  module mopenmpcut
  use mparameter
  implicit none
  integer(int_p)::nthds 
  integer(int_p)::kstart(ngz),kend(ngz)
  integer(int_p)::jstart(ngy),jend(ngy)
  end module mopenmpcut 
! =======================================================
! END MODULE OPENMPCUT ==================================
! =======================================================


! *******************************************************
! START SUBROUTINE OPENMPCUT ****************************
! *******************************************************
  subroutine openmpcut(nn)
  use mparameter
  use mopenmpcut
  use msetup
  implicit none
  integer(int_p)::ndelta,nblock,itmp,thd,nn
   
  nthds = nn
  
  nblock = nzbox / nthds
  ndelta = mod( nzbox , nthds )

  if ( ndelta == 0 ) then
   kstart(1) = 1 
   kend(1)   = kstart(1) + nblock - 1
   do thd = 2 , nthds
    kstart(thd) = kend(thd-1)   + 1
    kend(thd)   = kstart(thd)   + nblock - 1
   end do
  else
   kstart(1) = 1
   kend(1)   = kstart(1) + nblock - 1 + 1
   do thd = 2 , nthds
    itmp = 0
    if ( ndelta >= thd ) itmp = 1
    kstart(thd) = kend(thd-1)   + 1
    kend(thd)   = kstart(thd)   + nblock - 1 + itmp
   end do
  end if

  nblock = nybox / nthds
  ndelta = mod( nybox , nthds )

  if ( ndelta == 0 ) then
   jstart(1) = 1
   jend(1)   = jstart(1) + nblock - 1
   do thd = 2 , nthds
    jstart(thd) = jend(thd-1)   + 1
    jend(thd)   = jstart(thd)   + nblock - 1
   end do
  else
   jstart(1) = 1
   jend(1)   = jstart(1) + nblock - 1 + 1
   do thd = 2 , nthds
    itmp = 0
    if ( ndelta >= thd ) itmp = 1
    jstart(thd) = jend(thd-1)   + 1
    jend(thd)   = jstart(thd)   + nblock - 1 + itmp
   end do
  end if

  return
  end subroutine openmpcut
! =======================================================
! END SUBROUTINE OPENMPCUT ==============================
! =======================================================

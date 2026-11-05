!234567
! Full 3D grid (x,y extruded over k, z uniform): needed so Tecplot can pair
! the solution-only vel_*.plt files with a grid of identical size.
subroutine write_grid
  use mgeomdata
  use msetup
  use mtecio
  implicit none
  integer(int_p):: k
  real(real_p), allocatable:: OutArray(:,:,:)

  write(*,*) 'Start writing grid...'

  call tec_open('Grid','x y z',trim(outdir)//'grid.plt',1)
  call tec_zone('GRID',ngx,ngy,ngz,3,0.d0,0)

  allocate( OutArray(ngx,ngy,ngz) )

  do k = 1, ngz
    OutArray(:,:,k) = xpg
  end do
  call tec_data(ngx*ngy*ngz, OutArray)

  do k = 1, ngz
    OutArray(:,:,k) = ypg
  end do
  call tec_data(ngx*ngy*ngz, OutArray)

  do k = 1, ngz
    OutArray(:,:,k) = dz * real( k - 1, real_p )
  end do
  call tec_data(ngx*ngy*ngz, OutArray)

  deallocate( OutArray )
  call tec_close

  write(*,*) 'End writing grid...'
  return
end subroutine write_grid

!===============================================================
! Compact grid: 2D x-y plane zone + 1D z-axis zone (~1.5 MB).
! It cannot be paired with the 3D vel_*.plt files (different zone size);
! use it for 2D work only.
!===============================================================
subroutine write_grid_compact
  use mgeomdata
  use msetup
  use mtecio
  implicit none
  integer(int_p):: k
  real(real_p), allocatable:: Zero(:), Zaxis(:)

  write(*,*) 'Start writing compact grid...'

  call tec_open('GridCompact','x y z',trim(outdir)//'grid.plt',1)

  ! Zone 1: x-y plane (z placeholder = 0)
  allocate( Zero(max(ngx*ngy,ngz)) )
  Zero = 0.d0
  call tec_zone('XY_PLANE',ngx,ngy,1,3,0.d0,0)
  call tec_data(ngx*ngy, xpg)
  call tec_data(ngx*ngy, ypg)
  call tec_data(ngx*ngy, Zero)

  ! Zone 2: z axis (x,y placeholders = 0)
  allocate( Zaxis(ngz) )
  do k = 1, ngz
    Zaxis(k) = dz * real( k - 1, real_p )
  end do
  call tec_zone('Z_AXIS',1,1,ngz,3,0.d0,0)
  call tec_data(ngz, Zero)
  call tec_data(ngz, Zero)
  call tec_data(ngz, Zaxis)

  deallocate( Zero, Zaxis )
  call tec_close

  write(*,*) 'End writing compact grid.'
  return
end subroutine write_grid_compact

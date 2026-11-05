!=========================================================================
! Analytic profiles for velocity as well as polymeric conformation tensor
! are derived. Newton's method is employed to find correct pressure gradient
! in order to maintain correct mass flow rate.
!
! General derivations for the laminar base state for viscoelastic flow are
! given in Cruz et. al. (JNNFM, 2005) and Richter et. al. (JFM, 2010)
!=========================================================================

!=========================================================================
module analyticField
	integer :: ngx_u, ngx_vw, ngy_v, ngy_uw, ngz_w, ngz_uv
	real    :: Re, pi
	real, dimension(:), allocatable :: x_u, x_vw, y_v, y_uw, z_w, z_uv
	real, dimension(:), allocatable :: Ub
end module
!=========================================================================

subroutine initialField_analytic()
	use initialField
	use analyticField

	integer :: FENE 

	allocate(x_u (0:nixm))
	allocate(x_vw(0:nix ))
	allocate(y_v (0:niym))
	allocate(y_uw(0:niy ))
	allocate(z_w (0:nizm))
	allocate(z_uv(0:niz ))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Generate mean flow
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Read flow parameters
	call readFloat("Re"    ,Re  )

	pi   = acos(-1.0)		! definition of pi

	z_w (:)     =z (1:niz )	!z_w  has ( niz ) elements [ 0  Lz]
	z_uv(1:nizm)=zm(1:nizm)	!z_uv has (niz+1) elements [ 0  Lz]
	z_uv(0)  =0
	z_uv(niz)=zL

	x_u (:)     =x (1:nix )	!x_u  has ( nix ) elements [ 0  Lx]
	x_vw(1:nixm)=xm(1:nixm)	!x_vw has (nix+1) elements [ 0  Lx]
	x_vw(0)  =0
	x_vw(nix)=xL

	y_v (:)     =y (1:niy )-(yL/2.)	!y_v  has ( niy ) elements [-Ly/2  Ly/2]
	y_uw(1:niym)=ym(1:niym)-(yL/2.)	!y_uw has (niy+1) elements [-Ly/2  Ly/2]
	y_uw(0)  =-yL/2.
	y_uw(niy)= yL/2.

	allocate( Ub (0:niy) )

	! initialise the fields
	Ub =0.0

	!end if

	! initialise the field
	U=0.0; V=0.0; W=0.0

	call ComputeInitialFieldStaggeredUVW_mean()

	deallocate(Ub)

	return
end
!=======================================================================

subroutine initialField_info(iunit)
        use analyticField

                call sectionTitle(iunit, "Flow")
                write(iunit, 13) Re

 13     format ( 4x, "Flow parameters :" / &
                12x, "Re   =", f10.3     )

        return
end
!=======================================================================



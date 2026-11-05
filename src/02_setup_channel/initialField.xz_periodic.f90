!=========================================================================
! Initial fields are created here for channel flows
! You can impose one of initial flow conditions:
!
! (1) No flow
! (2) Uniform profile for velocity
! (3) Poiseuille profile for velocity
! (4) Dummy fields
!
! Sang Jin Lee
! Imperial College, 2015
!=========================================================================

module initialField
	include "mesh.inc"
	real, dimension(:,:,:), allocatable :: U, V, W

	integer			:: BC_bottom, BC_top
	integer			:: ifield, iOSS, iPeriodic, ierr
end module

!=========================================================================
subroutine initialField_define()
	use initialField

	allocate( U (0:niz  ,0:nix+1,0:niy  ))
	allocate( V (0:niz  ,0:nix  ,0:niy+1))
	allocate( W (0:niz+1,0:nix  ,0:niy  ))

	! Initialise velocity field to zero
	U = 0.0
	V = 0.0
	W = 0.0

	! Get type of initial field from archive
	call readInt("ifield", ifield)

	select case (ifield)
	case (1)
		! zero initial field: do nothing
	case (2)
		call initialField_uniform()
	case (3)
		! STOP "Option for PARABOLIC initial field is not available"
		call initialField_parabolic()
	case(4) 
		call initialField_dummy()
	end select

	!---- Fluctuations ----
	call readInt("iturb", iturb)
	if (iturb .ne. 0) then
		! STOP "Option for initial fluctuations is not available"
		  call initialField_fluct_v3()
	end if

	!---- bottom & top BCs ----
	call readInt("BC_bottom", BC_bottom)
	call readInt("BC_top",    BC_top)

	!---- perform a check of the BCs ----
	if (BC_bottom.ne.BC_top) then
		stop 'BC_top must equal BC_bottom'
	end if

	!---- impose boundary conditions ----
	call initialField_boundaries_UVW

	return
end

!=======================================================================
subroutine initialField_uniform()
	use initialField

	call readFloat("ubulk", ubulk)
	call readFloat("wbulk", wbulk)

	! Streamwise flow
	do j=jmin,jmax
		U(:,:,j) = ubulk
	end do

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	do j=jmin,jmax
		W(:,:,j) = wbulk
	end do

	return
end

!=======================================================================
subroutine initialField_parabolic()
	use initialField

	real, dimension(0:niy+1) :: U_v, W_v
	double precision :: deltay

	U_v = 0.0; W_v = 0.0

	call readFloat("ubulk", ubulk)
	call readFloat("wbulk", wbulk)

	! Streamwise flow
	y1 = y(jmin)		!lower wall
	y2 = y(jmax)		!upper wall
	yc = 0.5*(y2+y1)+y1	!centre-line position
	yd = 0.5*(y2-y1)	!domain half-height

	! Create streamwise baseflow on collocated grid    
	do j=jmin,jmax
			U_v(j) = 1.5*ubulk*(1.0 - (ypg(1,j)-yc)**2)
	enddo

	! Relocate baseflow on U-staggered grid
	do j=jmin,jmax-1
			U(:,:,j) = 0.5*(U_v(j)+U_v(j+1))
	enddo

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	y1 = y(jmin)		!lower wall
	y2 = y(jmax)		!upper wall
	yc = 0.5*(y2+y1)+y1	!centre-line position
	yd = 0.5*(y2-y1)	!domain half-height

	! Create spanwise baseflow on collocated grid
	do j=jmin,jmax
			W_v(j) = 1.5*wbulk*(1.0 - (ypg(1,j)-yc)**2)
	enddo

	! Relocate baseflow on W-staggered grid
	do j=jmin,jmax-1
			W(:,:,j) = 0.5*(W_v(j)+W_v(j+1))
	enddo

	return
end

!==========================================================================
subroutine initialField_dummy()
	use initialField

	!--- define pi ---
	twopi = 2*2.*acos(0.) 

	!--- compute dz ---
	dz = zL/dble(niz-1)

	!--- create a dummy baseflow ---
	!do k=1,niz-1
	do k=1,niz
		U(k,:,:) =    SIN(16.*twopi*(k-1)*dz/zL)
		V(k,:,:) = 2.*SIN(16.*twopi*(k-1)*dz/zL)
		W(k,:,:) =    COS(16.*twopi*(k-1)*dz/zL)
	end do

	print*,W(1:niz-1,1,1)

	!print*,MAXVAL(U),MAXVAL(V),MAXVAL(W)
	
	print*,W(:,64,64)

	return
end

!=======================================================================
subroutine initialField_fluct_v3()
	use initialField

	        ! Add random fluctuations

        ! Initialize random generator
        iseed = 7
        call randomSeed(iseed)

        ! Fluctuation amplitude is typically 30%
        amplitude = 0.3
        write(*,*) 'amplitude =', amplitude

        do j=jmin,jmax
        do k=kmin,kmax
        do i=imin,imax
                umean = max(U(k,i,j), V(k,i,j), W(k,i,j))
                fluct1 = umean*amplitude*2.*(random()-0.5)
                fluct2 = umean*amplitude*2.*(random()-0.5)
                fluct3 = umean*amplitude*2.*(random()-0.5)
                U(k,i,j) = U(k,i,j) + fluct1
                V(k,i,j) = V(k,i,j) + fluct2
                W(k,i,j) = W(k,i,j) + fluct3
        end do
        end do
        end do

	return
end

!=======================================================================
subroutine initialField_write(iunit)  !JL160206
	use initialField
        integer(8) ilength

	inquire(iolength=ilength) U
	! Reopen file as direct access
	open (iunit, file="01uc00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) U
	close(iunit)

	inquire(iolength=ilength) V
	! Reopen file as direct access
	open (iunit, file="02vc00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) V
	close(iunit)
	inquire(iolength=ilength) W
	! Reopen file as direct access
	open (iunit, file="03wc00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) W
	close(iunit)

	return
end

!=======================================================================
subroutine initial_Flux_write(iunit)  !JL160206
        use initialField
        integer(8) ilength

        U = 0.0
        V = 0.0
        W = 0.0

	inquire(iolength=ilength) U
	! Reopen file as direct access
	open (iunit, file="04uu00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) U
	close(iunit)


	inquire(iolength=ilength) V
	! Reopen file as direct access
	open (iunit, file="05vv00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) V
	close(iunit)

	inquire(iolength=ilength) W
	! Reopen file as direct access
	open (iunit, file="06ww00.dble.0000000", form="unformatted", &
	      access="direct",recl=ilength)
	write (iunit,rec = 1) W
	close(iunit)

        return
end

!=======================================================================
subroutine initial_Con_write(iunit)
        use initialField
        real, dimension(niz-1+1,nix-1+1,niy-1+1) :: con

        con = 0.0
        !TJ:  Use a record length of (niz-1+1)
        inquire(iolength=ilength) con(1,1,1)
        ilength = ilength * ( (nix-1+1)*(niy-1+1)*(niz-1+1) )

        ! Reopen file as direct access
        close(iunit)
        open (iunit, file="conold.dble.0000000", form="unformatted", &
                access="direct",recl=ilength)

        write (iunit,rec = 1) con
        write (iunit,rec = 2) con
        write (iunit,rec = 3) con

        return
end

!=======================================================================
subroutine initial_phatr_write(iunit)
        use initialField
        real, dimension(0:nix, 0:niy , niz-1) :: phatr

        phatr = 0.0
        !T  Use a record length of (nix+1)*(niy+1)*(niz-1)
        inquire(iolength=ilength) phatr(:,:,:)

        ! Reopen file as direct access
        close(iunit)
        open (iunit, file="pressure_ph.0000000", form="unformatted", &
                access="direct",recl=ilength)

        write (iunit,rec = 1) phatr

        return
end

!=======================================================================
subroutine initialField_boundaries_UVW
	use initialField

	!---- Wall BCs ----
	V(:,:, jmin ) =  0.0
	V(:,:, jmax ) =  0.0

	!---- Wall BCs ----
	SELECT CASE(BC_bottom)
	CASE(0)
		U(:,:,jminm1) = -U(:,:, jmin )
		V(:,:,jminm1) =  V(:,:,jminp1)
		W(:,:,jminm1) = -W(:,:, jmin )
	CASE(1)
		U(:,:,jminm1) =  U(:,:, jmin )
		V(:,:,jminm1) =  V(:,:,jminp1)
		W(:,:,jminm1) =  W(:,:, jmin )
	END SELECT

	!---- Wall BCs ----
	SELECT CASE(BC_top)
	CASE(0)
		U(:,:, jmax ) = -U(:,:,jmaxm1)
		V(:,:,jmaxp1) =  V(:,:,jmaxm1)
		W(:,:, jmax ) = -W(:,:,jmaxm1)
	CASE(1)
		U(:,:, jmax ) =  U(:,:,jmaxm1)
		V(:,:,jmaxp1) =  V(:,:,jmaxm1)
		W(:,:, jmax ) =  W(:,:,jmaxm1)
	END SELECT
        
	!---- Periodicity in x ----
	U(:,iminm1,:) = U(:,imaxm1,:)
	U(:,imaxp1,:) = U(:,iminp1,:)
	V(:,iminm1,:) = V(:,imaxm1,:)
	V(:, imax, :) = V(:, imin, :)
	W(:,iminm1,:) = W(:,imaxm1,:)
	W(:, imax, :) = W(:, imin, :)

	!---- Watch out: The following is important ----
	!U(:, imax, :) = 0.5*(U(:, imax, :) + U(:, imin, :))
	U(:, imin, :) = U(:, imax, :)

	!---- Periodicity in z ----
	U(kminm1,:,:) = U(kmaxm1,:,:)
	U( kmax ,:,:) = U( kmin ,:,:)
	V(kminm1,:,:) = V(kmaxm1,:,:)
	V( kmax ,:,:) = V( kmin ,:,:)
	W(kminm1,:,:) = W(kmaxm1,:,:)
	W(kmaxp1,:,:) = W(kminp1,:,:)

	!---- Watch out: The following is important ----
	!W( kmax ,:,:) = 0.5 * (W(kmax,:,:) + W(kmin,:,:))
	W( kmin ,:,:) = W(kmax,:,:)

	return
end

!=======================================================================
subroutine Inflow_BC_write(iunit)
	use initialField

	return
end

!=======================================================================


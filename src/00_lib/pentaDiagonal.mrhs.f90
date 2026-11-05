!=======================================================================
! Penta-diagonal solvers
!
! Charles Pierce, June 1995
!
! pentaDiagonal_1_Mrhs
! pentaDiagonalM_1_Mrhs
! pentaDiagonalC_mrhs (Complex)
!
!=======================================================================
! pentaDiagonal_1_Mrhs
!
! Combined matrix-rhs into single array
!
subroutine pentaDiagonal_1_Mrhs(A,r,n,nrhs)
        integer n, nrhs, i
	real*8 A(5,n), r(nrhs,n)
        real*8 const

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)/A(3,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const = A(2,2)/A(3,1)
		A(3,2) = A(3,2) - A(4,1)*const
                r(:,2) = r(:,2) - r(:,1)*const
		r(:,2) = r(:,2)/A(3,2)
		r(:,1) = (r(:,1) - A(4,1)*r(:,2))/A(3,1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate A(2,i+1)
		const = A(2,i+1)/A(3,i)
		A(3,i+1) = A(3,i+1) - A(4,i)*const
		A(4,i+1) = A(4,i+1) - A(5,i)*const
		r(:,i+1) = r(:,i+1) - r(:,i)*const

		! Eliminate A(1,i+2)
		const = A(1,i+2)/A(3,i)
		A(2,i+2) = A(2,i+2) - A(4,i)*const
		A(3,i+2) = A(3,i+2) - A(5,i)*const
		r(:,i+2) = r(:,i+2) - r(:,i)*const
	end do
	! Eliminate A(2,n)
	const = A(2,n)/A(3,n-1)
	A(3,n) = A(3,n) - A(4,n-1)*const
	r(:,n) = r(:,n) - r(:,n-1)*const

	! Back-substitution
	r(:,n) = r(:,n)/A(3,n)
	r(:,n-1) = (r(:,n-1) - A(4,n-1)*r(:,n))/A(3,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - A(4,i)*r(:,i+1) - A(5,i)*r(:,i+2))/A(3,i)
	end do
    
	return
end

!=======================================================================
! pentaDiagonalM_1_Mrhs
!
! Combined matrix-rhs into single array
!
subroutine pentaDiagonalM_1_Mrhs(A, r, n, lot, nrhs)
        integer n, lot, nrhs, i, j
	real*8 A(lot,5,n), r(nrhs,lot,n)
	real*8 const(lot)

	if (n == 1) then
		! Solve 1x1 system
                do j=1,lot
                    r(:,j,1) = r(:,j,1)/A(j,3,1)
                end do
		return
	else if (n == 2) then
		! Solve 2x2 system
		const(:) = A(:,2,2)/A(:,3,1)
		A(:,3,2) = A(:,3,2) - A(:,4,1)*const(:)
                do j=1,lot
                    r(:,j,2) = r(:,j,2) - r(:,j,1)*const(j)
                    r(:,j,2) = r(:,j,2)/A(j,3,2)
                    r(:,j,1) = (r(:,j,1) - A(j,4,1)*r(:,j,2))/A(j,3,1)
                end do
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate A(2,i+1)
		const(:) = A(:,2,i+1)/A(:,3,i)
		A(:,3,i+1) = A(:,3,i+1) - A(:,4,i)*const(:)
		A(:,4,i+1) = A(:,4,i+1) - A(:,5,i)*const(:)
                do j=1,lot
                    r(:,j,i+1) = r(:,j,i+1) - r(:,j,i)*const(j)
                end do

		! Eliminate A(1,i+2)
		const(:) = A(:,1,i+2)/A(:,3,i)
		A(:,2,i+2) = A(:,2,i+2) - A(:,4,i)*const(:)
		A(:,3,i+2) = A(:,3,i+2) - A(:,5,i)*const(:)
                do j=1,lot
                    r(:,j,i+2) = r(:,j,i+2) - r(:,j,i)*const(j)
                end do
	end do
	! Eliminate A(2,n)
	const(:) = A(:,2,n)/A(:,3,n-1)
	A(:,3,n) = A(:,3,n) - A(:,4,n-1)*const(:)
        do j=1,lot
            r(:,j,n) = r(:,j,n) - r(:,j,n-1)*const(j)
        end do

	! Back-substitution
        do j=1,lot
            r(:,j,n) = r(:,j,n)/A(j,3,n)
        end do

        do j=1,lot
            r(:,j,n-1) = (r(:,j,n-1) - A(j,4,n-1)*r(:,j,n))/A(j,3,n-1)
        end do

	do i=n-2,1,-1
            do j=1,lot
                r(:,j,i) = (r(:,j,i) - A(j,4,i)*r(:,j,i+1) - A(j,5,i)*r(:,j,i+2))/A(j,3,i)
            end do
	end do

	return
end


!=======================================================================
! pentaDiagonalC_mrhs
! Solves a penta-diagonal system of algebraic equations.  The banded
! matrix given by (a,b,c,d,e) is destroyed, and the rhs vector r is
! replaced by the solution vector.
! All variables are complex.
!
! Multiple systems, different (a,b,c,d,e) for each system.
!
subroutine pentaDiagonalC_mrhs(a, b, c, d, e, r, nrhs, n)
	complex	a(n), b(n), c(n), d(n), e(n), r(nrhs,n)
	complex	const

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)/c(1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const  = b(2)/c(1)
		c(2)   = c(2) - d(1)*const
		r(:,2) = r(:,2) - r(:,1)*const
		r(:,2) = r(:,2)/c(2)
		r(:,1) = (r(:,1) - d(1)*r(:,2))/c(1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		const    = b(i+1)/c(i)
		c(i+1)   = c(i+1) - d(i)*const
		d(i+1)   = d(i+1) - e(i)*const
		r(:,i+1) = r(:,i+1) - r(:,i)*const

		! Eliminate a(i+2)
		const    = a(i+2)/c(i)
		b(i+2)   = b(i+2) - d(i)*const
		c(i+2)   = c(i+2) - e(i)*const
		r(:,i+2) = r(:,i+2) - r(:,i)*const
	end do
	! Eliminate b(n)
	const  = b(n)/c(n-1)
	c(n)   = c(n) - d(n-1)*const
	r(:,n) = r(:,n) - r(:,n-1)*const

	! Back-substitution
	r(:,n) = r(:,n)/c(n)
	r(:,n-1) = (r(:,n-1) - d(n-1)*r(:,n))/c(n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - d(i)*r(:,i+1) - e(i)*r(:,i+2))/c(i)
	end do

	return

end

subroutine pentaPeriodic_mrhs_1(A, r, nrhs, n)
        real A(5,n), r(nrhs,n)
        real s1(n), s2(n)
        real const

        if (n == 1) then
                ! Solve 1x1 system
                r = r/A(3,1)
                return
        else if (n == 2) then
                ! Solve 2x2 system
!               stop "error, end up with 2x2 pentadiagonal system!!!"
                A(3,:) = A(3,:) + A(1,:) + A(5,:)
                A(4,1) = A(4,1) + A(2,1)
                A(2,2) = A(2,2) + A(4,2)
                const = A(2,2)/A(3,1)
                A(3,2) = A(3,2) - A(4,1)*const
                r(:,2) = r(:,2) - r(:,1)*const
                r(:,2) = r(:,2)/A(3,2)
                r(:,1) = (r(:,1) - A(4,1)*r(:,2))/A(3,1)
                return
        else if (n == 3) then
!               stop "error, end up with 3x3 pentadiagonal system!!!"
                A(2,:) = A(2,:) + A(5,:)
                A(4,:) = A(4,:) + A(1,:)
                call triDiagPeriodic_mrhs(A(2,:), A(3,:), A(4,:), r, nrhs, n)
                return
        else if (n == 4) then
                A(1,:) = A(1,:) + A(5,:)
                A(5,:) = 0.
        end if

        ! Initialize boundary data
        s1 = 0.
        s1(1) = A(1,1)
        s1(n-3) = s1(n-3) + A(5,n-3)
        s1(n-2) = A(4,n-2)
        s1(n-1) = A(3,n-1)
        s1(n) = A(2,n)
        s2 = 0.
        s2(1) = A(2,1)
        s2(2) = A(1,2)
        s2(n-2) = s2(n-2) + A(5,n-2)
        s2(n-1) = A(4,n-1)
        s2(n) = A(3,n)

        ! Forward elimination
        do i=1,n-2
                ! Eliminate A(2,i+1)
                const = A(2,i+1)/A(3,i)
                A(3,i+1) = A(3,i+1) - A(4,i)*const
                A(4,i+1) = A(4,i+1) - A(5,i)*const
                r(:,i+1) = r(:,i+1) - r(:,i)*const
                s1(i+1) = s1(i+1) - s1(i)*const
                s2(i+1) = s2(i+1) - s2(i)*const

                ! Eliminate A(1,i+2)
                const = A(1,i+2)/A(3,i)
                A(2,i+2) = A(2,i+2) - A(4,i)*const
                A(3,i+2) = A(3,i+2) - A(5,i)*const
                r(:,i+2) = r(:,i+2) - r(:,i)*const
                s1(i+2) = s1(i+2) - s1(i)*const
                s2(i+2) = s2(i+2) - s2(i)*const
        end do

        ! Backward elimination
        do i=n-2,3,-1
                ! Eliminate A(4,i-1)
                const = A(4,i-1)/A(3,i)
                r(:,i-1) = r(:,i-1) - r(:,i)*const
                s1(i-1) = s1(i-1) - s1(i)*const
                s2(i-1) = s2(i-1) - s2(i)*const

                ! Eliminate A(5,i-2)
                const = A(5,i-2)/A(3,i)
                r(:,i-2) = r(:,i-2) - r(:,i)*const
                s1(i-2) = s1(i-2) - s1(i)*const
                s2(i-2) = s2(i-2) - s2(i)*const
        end do
        i=2
                ! Eliminate A(4,i-1)
                const = A(4,i-1)/A(3,i)
                r(:,i-1) = r(:,i-1) - r(:,i)*const
                s1(i-1) = s1(i-1) - s1(i)*const
                s2(i-1) = s2(i-1) - s2(i)*const

        ! Eliminate oddball region
        const = A(5,n-1)/A(3,1)
        r(:,n-1) = r(:,n-1) - r(:,1)*const
        s1(n-1) = s1(n-1) - s1(1)*const
        s2(n-1) = s2(n-1) - s2(1)*const

        const = A(4,n)/A(3,1)
        r(:,n) = r(:,n) - r(:,1)*const
        s1(n) = s1(n) - s1(1)*const
        s2(n) = s2(n) - s2(1)*const

        const = A(5,n)/A(3,2)
        r(:,n) = r(:,n) - r(:,2)*const
        s1(n) = s1(n) - s1(2)*const
        s2(n) = s2(n) - s2(2)*const

        ! Eliminate corner region
        const = s1(n)/s1(n-1)
        r(:,n) = r(:,n) - r(:,n-1)*const
        s2(n) = s2(n) - s2(n-1)*const

        r(:,n) = r(:,n)/s2(n)
        r(:,n-1) = (r(:,n-1) - s2(n-1)*r(:,n))/s1(n-1)
        do i=n-2,1,-1
                r(:,i) = (r(:,i) - s1(i)*r(:,n-1) - s2(i)*r(:,n))/A(3,i)
        end do

        return
end


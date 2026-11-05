!=======================================================================
! Tri-diagonal solvers
!
! Charles Pierce, June 1995
!
! triDiagonal_mrhs
! triDiagonalC_mrhs
! triDiagonalM_mrhs
!

!=======================================================================
! triDiagonal_mrhs
! Solves a tri-diagonal system of algebraic equations using gauss 
! elimination.  The banded matrix given by (a,b,c) is destroyed, 
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
!
subroutine triDiag_mrhs(a, b, c, r, nrhs, n)
	real a(n), b(n), c(n), r(nrhs,n)

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(:,i) = r(:,i) - r(:,i-1)*const
	end do

	! Back substitution
	r(:,n) = r(:,n)/b(n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(i)*r(:,i+1))/b(i)
	end do

	return
end

!=======================================================================
! triDiagonalM
! Solves an array of tri-diagonal systems of algebraic equations using 
! gauss elimination.  The banded matrices given by (a,b,c) are 
! destroyed, and the rhs vectors r are replaced by the solution vectors.
!
! Multiple systems, different (a,b,c) for each system.
!
subroutine triDiagM_mrhs(a, b, c, r, n, lot, nrhs)
	real a(lot,n), b(lot,n), c(lot,n), r(nrhs,lot,n)
	real const(lot)

	! Forward elimination
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
            do j=1,lot
		r(:,j,i) = r(:,j,i) - r(:,j,i-1)*const(j)
            end do
	end do

	! Back-substitution
        do j=1,lot
	   r(:,j,n) = r(:,j,n)/b(j,n)
        end do

	do i=n-1,1,-1
           do j=1,lot
		r(:,j,i) = (r(:,j,i) - c(j,i)*r(:,j,i+1))/b(j,i)
           end do
	end do

	return
end


!=======================================================================
! triDiagonalC_mrhs
! Solves a tri-diagonal system of algebraic equations using gauss
! elimination.  The banded matrix given by (a,b,c) is destroyed,
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
! All Variables are Complex
!
subroutine triDiagC_mrhs(a, b, c, r, nrhs, n)
	complex a(n), b(n), c(n), r(nrhs,n)
	complex const

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(:,i) = r(:,i) - r(:,i-1)*const
	end do

	! Back substitution
	r(:,n) = r(:,n)/b(n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(i)*r(:,i+1))/b(i)
	end do

return
end

subroutine triDiagPeriodic_mrhs(a, b, c, r, nrhs, n)
        real a(n), b(n), c(n), r(nrhs,n)

        if (n == 1) then
                r = r/(a(1) + b(1) + c(1))
                return
        else if (n == 2) then
                ! Solve 2x2 system
                c(1) = c(1) + a(1)
                a(2) = a(2) + c(2)
                const = a(2)/b(1)
                b(2) = b(2) - c(1)*const
                r(:,2) = r(:,2) - r(:,1)*const
                r(:,2) = r(:,2)/b(2)
                r(:,1) = (r(:,1) - c(1)*r(:,2))/b(1)
                return
        end if

        ! Forward elimination
        do i=2,n-2
                const = a(i)/b(i-1)
                b(i) = b(i) - c(i-1)*const
                r(:,i) = r(:,i) - r(:,i-1)*const
                ! Boundary is stored in a(i)
                a(i) = -a(i-1)*const
        end do
        i=n-1
                const = a(i)/b(i-1)
                b(i) = b(i) - c(i-1)*const
                r(:,i) = r(:,i) - r(:,i-1)*const
                a(i) = c(i) - a(i-1)*const
        i=n
                const = a(i)/b(i-1)
                r(:,i) = r(:,i) - r(:,i-1)*const
                a(i) = b(i) - a(i-1)*const

        ! Backward elimination
        do i=n-2,1,-1
                const = c(i)/b(i+1)
                r(:,i) = r(:,i) - r(:,i+1)*const
                a(i) = a(i) - a(i+1)*const
        end do

        ! Eliminate oddball
        const = c(n)/b(1)
        r(:,n) = r(:,n) - r(:,1)*const
        a(n) = a(n) - a(1)*const

        ! Backward substitution
        r(:,n) = r(:,n)/a(n)
        do i=n-1,1,-1
                r(:,i) = (r(:,i) - a(i)*r(:,n))/b(i)
        end do

        return
end

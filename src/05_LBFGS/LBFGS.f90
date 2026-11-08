! MZW 01/14/23

program main

	call LBFGS


end program



subroutine LBFGS
!==============================================================
! Get user-specified parameters
!
      include "param.inc"
!==============================================================
	integer :: i, i1, j , l , m , n, i_start
        integer ::  nbc_opti, m_store, dimv, nf, nfl, line_search, restart
	DOUBLE PRECISION :: J0,J1,J0_back
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: desdir,perICvecstock,gradientICvecstock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: svec,yvec
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: perICvec,gradientICvec	
	DOUBLE PRECISION :: J0_esti,J1_esti,jcostfunction,reduc_gra
        character(len=50) :: folder, fname1, fname2, command_line

	double precision :: alpha,tol, ftol, gtol, xtol, stpmin, stpmax		! for line search 
	integer :: maxfev, nfev, info, nfev_sum				! for line search
	LOGICAL :: exitboucle
	
	! Initialize parameters
	nbc_opti = 10           ! maximum # of iterations
   restart = 0             ! restart from a certain iteration.  0: no restart
	reduc_gra = 1.e-10
	m_store = 5 
        line_search = 2     ! 1 for Armijo-backtracking line search; 
                            ! 2 for line search that satisfy both Armijo and strong Wolfe condition

	! initial condition only (periodic in x, z; no-slip walls in y), see FUNCT.f90
	dimv = (ngx-1)*(ngy-1)*(ngz-1) + (ngx-1)*(ngy-2)*(ngz-1) + (ngx-1)*(ngy-1)*(ngz-1)
	ALLOCATE(desdir(dimv),gradientICvecstock(dimv),perICvecstock(dimv))
	ALLOCATE(svec(m_store,dimv),yvec(m_store,dimv),perICvec(dimv),gradientICvec(dimv))
	

        if(restart == 0) then	
        	!--------------------------------------------------------
        	!        read initial guess and evaluate f, g
        	!--------------------------------------------------------
        	call Initial_guess(dimv,perICvec)
        	
        	call Evaluate_f(dimv,perICvec,jcostfunction)
        	call Evaluate_g(dimv,gradientICvec)
                
        	nf = 1
        	nfev_sum = 0 	! total number of function-gradient evaluations
                i_start = 1
        	J0=jcostfunction
        	J0_back=jcostfunction
        	
        	! Optimization loop
        	open(300,file="LBFGS_result.plt")
        	open(301,file="LBFGS_alpha.plt")
        	open(302,file="line_search_check")
	        write(300,*) 0, jcostfunction,maxval(abs(gradientICvec))

        else
                !-------------------------------------------------------
                ! open files to read variables for restart
                !-------------------------------------------------------
                folder = './LBFGS_backup/'
                fname1 = 'LBFGS_iter.' 
                write(fname2,'(A15,A11,I3.3)') folder, fname1, restart
                open(10,file=fname2,form="formatted")
                read(10,*) i_start
                i_start = i_start + 1
                read(10,*) nfev_sum
                read(10,*) nf
                read(10,*) jcostfunction
                read(10,*) J0
                close(10)
        
	   if ( restart <= m_store) then
                i1 = 0
           else 
                i1 = restart-m_store
           end if
           j = 1
	   do i=i1,restart
                   inquire(iolength=ilength) perICvec
                   fname1 = 'ICvec.'
                   write(fname2,'(A15,A6,I3.3)') folder, fname1, i
                   write(*,*) 'reading:', fname2
                   open(11,file=fname2,form="unformatted",access="direct", &
                           recl=ilength,iostat=ierr)
                   if(ierr .ne. 0) stop "perICvec.dat cannot be created!"
                   read(11,rec=1) perICvec
                   close(11)

                   inquire(iolength=ilength) gradientICvec
                   fname1 = 'gradvec.'
                   write(fname2,'(A15,A8,I3.3)') folder, fname1, i
                   open(12,file=fname2,form="unformatted",access="direct", &
                           recl=ilength,iostat=ierr)
                   if(ierr .ne. 0) stop "grad.dat cannot be created!"
                   read(12,rec=1) gradientICvec
                   close(12)
                   if(i > i1) then
                           svec(j,:) = perICvec(:) - perICvecstock(:)
                           yvec(j,:) = gradientICvec(:) - gradientICvecstock(:)
                           j = j + 1
                   end if
                   perICvecstock = perICvec
                   gradientICvecstock = gradientICvec
           end do

        	open(300,file="LBFGS_result.plt",status="old",position="append")
        	open(301,file="LBFGS_alpha.plt",status="old",position="append")
        	open(302,file="line_search_check",status="old",position="append")
        end if


        ! Optimization loop
	DO i=i_start,nbc_opti	
	!--------------------------------------------------------
	!        check resubmit file to stop
	!--------------------------------------------------------	
	inquire(file="resubmit", exist=iexist)
	if (iexist) then
		write(*,*) 'file exist!'
		exit
	end if
	
	!--------------------------------------------------------
	!        Compute the descent direction
	!--------------------------------------------------------	
	   IF (i == 1) THEN
	       desdir(1:dimv)=-gradientICvec(1:dimv) / norm2(gradientICvec)
	   ELSE
	       CALL Compute_desdir_new ( desdir , gradientICvec , svec , yvec , dimv , i , m_store)           
	   END IF
	
	   DO j=1,dimv
	      gradientICvecstock(j)=gradientICvec(j) ! save gradient
	      perICvecstock(j)=perICvec(j)	! save location
	   END DO

	!--------------------------------------------------------
	!        Line search
	!--------------------------------------------------------	
	if(line_search == 1) then
	   ! Armijo line search
	   !alpha = 1.0
	   !write(*,*) 'norm2',1./norm2(desdir)
           exitboucle = .FALSE.
   	   alpha = 1.0
   	   !if(i==1) alpha =  1.d0/norm2(gradientICvec)
	   !alpha = min(jcostfunction/norm2(desdir)**2, 1.0)
	   tol = 1.e-4 !0.0001
	   call linesearch_Armijo(dimv, perICvec, jcostfunction, gradientICvec, desdir, alpha, tol, nfl, exitboucle)
	   nf = nf + nfl
           nfev_sum = nfev_sum + 1
	   J0_back = jcostfunction
	   write(*,*) 'step size:',alpha
           write(*,*) 'function evaluations:',nfl
	   if(exitboucle .eq. .FALSE.) then
		write(302,*) i,'line search fail!'
	    end if
	   ! Armijo line search only return new perICvec and jcost function, gradient is not updated
	   call Evaluate_g(dimv,gradientICvec)

	elseif(line_search == 2) then
	   ! More & Thuente line search

   	   alpha = 1.d0
           if(i==1) alpha = min(1.0,jcostfunction/norm2(gradientICvec))
   	   !if(i==1) alpha = 1.d0/norm2(gradientICvec)
	   !if(i==1) alpha = min(jcostfunction/norm2(desdir)**2, 1.0)
           write(*,*) 'initial step size alpha = ',alpha
	   !alpha = min(1./norm2(desdir), 1.0e-5)
	   !alpha = min(jcostfunction/norm2(desdir)**2, 1.0)
	   ftol = 1.e-4
	   gtol = 0.9
	   xtol = 1.e-16
	   stpmin = 1.e-16
	   stpmax = 1.e16
	   maxfev = 20
	   ! last input is perICvecstock, will not be changed after line search. It's just a working array. 
	   call cvsrch(dimv,perICvec,jcostfunction,gradientICvec,desdir,alpha,ftol,gtol,xtol, &
			stpmin,stpmax,maxfev,info,nfev,perICvecstock)
           write(*,*) '--------------exit cvsrch---------------'
           !write(301,'(I3,6E11.4)') i,perICvec(1),perICvec(2),perICvec(3),perICvec(4),perICvec(5),perICvec(6)
           write(301,*) i, alpha
	   write(*,*) 'f-g evaluation:',nfev
	   nfev_sum = nfev_sum + nfev
           nf = nfev_sum
	   if (info /= 1) then   ! means something wrong in line search
		   if(info == 0) then
			write(302,*) 'Line search: Improper input parameters!'
		   else if(info == 2) then
			write(302,*) 'Line search: Relative width of the interval of uncertainty is at most xtol!'
		   else if(info == 3) then
			write(302,*) 'Line search: Number of calls to fcn has reached maxfev!'
		   else if(info == 4) then
			write(302,*) 'Line search: The step is at the lower bound stpmin!'
		   else if(info == 5) then
			write(302,*) 'Line search: The step is at the upper bound stpmax!'
		   else if(info == 6) then
			write(302,*) 'Line search: Tolerances may be too small!'
		   end if
                   write(*,*) 'line search failed! recompute the cost function and grad:'
                   do k=1,dimv
                        perICvec(k) = perICvecstock(k) + alpha * desdir(k)
                   end do
                   call Evaluate_f(dimv,perICvec,jcostfunction)
                   call Evaluate_g(dimv,gradientICvec)    
	   end if
	end if
	   ! line search algorithm returns:
	   ! perICvec:  the optimized location
	   ! jcostfunction:  function value at perICvec
	   ! gradientICvec:  gradient at perICvec
	
	!--------------------------------------------------------
	!        Save s and y for updating next Hessian matrix
	!--------------------------------------------------------
	   if ( i <= m_store) then
	   	do j=1,dimv
	      		svec(i,j)=perICvec(j)-perICvecstock(j)
	      		yvec(i,j)=gradientICvec(j)-gradientICvecstock(j)
	   	end do
	   else
	        do m=1,m_store-1
			svec(m,:) = svec(m+1,:)
			yvec(m,:) = yvec(m+1,:)
	        end do
		do j=1,dimv
			svec(m_store,j)=perICvec(j)-perICvecstock(j)
			yvec(m_store,j)=gradientICvec(j)-gradientICvecstock(j)
		end do 
	   END IF	
	
	   WRITE(6,*)"COST FUNCTION"
	   WRITE(6,*) jcostfunction,J0

        !-------------------------------------------------------
        ! open files to save variables for restart
        !-------------------------------------------------------
        folder = './LBFGS_backup/'
        fname1 = 'LBFGS_iter.' 
        write(fname2,'(A15,A11,I3.3)') folder, fname1, i
        open(10,file=fname2,form="formatted")
        write(10,*) i
        write(10,*) nfev_sum
        write(10,*) nf
        write(10,*) jcostfunction
        write(10,*) J0
        close(10)

        inquire(iolength=ilength) perICvec
        fname1 = 'ICvec.'
        write(fname2,'(A15,A6,I3.3)') folder, fname1, i
        open(11,file=fname2,form="unformatted",access="direct", &
                recl=ilength,iostat=ierr)
        if(ierr .ne. 0) stop "perICvec.dat cannot be created!"
        write(11,rec=1) perICvec
        close(11)

        inquire(iolength=ilength) gradientICvec
        fname1 = 'gradvec.'
        write(fname2,'(A15,A8,I3.3)') folder, fname1, i
        open(12,file=fname2,form="unformatted",access="direct", &
                recl=ilength,iostat=ierr)
        if(ierr .ne. 0) stop "grad.dat cannot be created!"
        write(12,rec=1) gradientICvec
        close(12)

        if(i==1) then       
                inquire(iolength=ilength) perICvec
                fname1 = 'ICvec.'
                write(fname2,'(A15,A6,I3.3)') folder, fname1, 0
                open(11,file=fname2,form="unformatted",access="direct", &
                        recl=ilength,iostat=ierr)
                if(ierr .ne. 0) stop "perICvec.dat cannot be created!"
                write(11,rec=1) perICvecstock
                close(11)

                inquire(iolength=ilength) gradientICvec
                fname1 = 'gradvec.'
                write(fname2,'(A15,A8,I3.3)') folder, fname1, 0
                open(12,file=fname2,form="unformatted",access="direct", &
                        recl=ilength,iostat=ierr)
                if(ierr .ne. 0) stop "grad.dat cannot be created!"
                write(12,rec=1) gradientICvecstock
                close(12)
        end if

	write(300,*) i, jcostfunction,maxval(abs(gradientICvec))
	!write(300,'(I3,3E24.17)') i, jcostfunction,gradientICvec(1),gradientICvec(2)
	
	   !WRITE(UNIT=1,FMT=*) jcostfunction-J1_esti,J1_esti,jcostfunction/J0
	   
	!--------------------------------------------------------
	!        stop criteria for LBFGS
	!--------------------------------------------------------	
	   if(maxval(abs(gradientICvec)) < reduc_gra*(1.+abs(jcostfunction))) exit
	
	END DO
        ! call one more forward simulation 
        call execute_command_line("cp data_zero data")
        write(*,*) '---------start final forward simulation----------'
        write(command_line,'(A,I3,A)') 'mpirun -np',nproc,' ./forward.out > history_forward'
        call execute_command_line(trim(command_line))
        write(*,*) '-----------end final forward simulation-----------'
        call execute_command_line("cp archive.0000000 archive")
	close(300)
	close(301)
	close(302)
	   WRITE(6,*) jcostfunction,J0
	! End of the optimization loop
	write(*,*) 'iterations:',i
	write(*,*) 'f evaluations',nf
	write(*,*) 'f-g evaluations:',nfev_sum
	!write(*,*) 'minimum poisition 1-3',perICvec(1),perICvec(2),perICvec(3)
end 



subroutine Compute_desdir_new ( desdir , gradient , svec , yvec , dimv , iter , m)
! this subroutine, the number of rows in svec, yvec is equal to m
! we only use m(=5) previous gradients to update Hessian matrix

!     see Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
implicit none

! dimV is the length of gradient vector
! svec is s in Nocedal's paper
! yvec is y in Nocedal's paper
! nbcda is the maximum number of optimization loops
integer, intent(in) :: dimV,iter,m
double precision, intent(in) :: gradient(dimV)
double precision, intent(in) :: svec(m,dimV),yvec(m,dimV)
double precision, intent(inout) :: desdir(dimV)

integer :: i,j,k,l, bound, incr
double precision, allocatable :: alpha(:), beta(:), rho(:)
double precision :: H0

! choose upper bound for loop
!m = 5

if (iter-1 < m) then
   bound = iter-1
   incr = 0
else
   bound = m
   incr = iter-1-m
end if

allocate(alpha(m),beta(m),rho(m))

desdir(:)= -gradient(:)

do i=bound,1,-1
   j = i !+ incr
   rho(j) = 1.d0/dot_product(yvec(j,:),svec(j,:))
   alpha(i) = rho(j)*dot_product(svec(j,:),desdir(:))
   desdir(:) = desdir(:) - alpha(i) * yvec(j,:)
end do

H0 = dot_product(svec(bound,:),yvec(bound,:))/dot_product(yvec(bound,:),yvec(bound,:))
desdir(:) = H0 * desdir(:)

do i=1,bound
   j = i !+ incr
   beta(j) = rho(j)*dot_product(yvec(j,:),desdir(:))
   desdir(:) = desdir(:) + (alpha(i) - beta(i)) * svec(j,:)
end do



end



!subroutine FUNCT(x, f, g, dimv)
!! simple test case
!implicit none
!
!integer, intent(in) :: dimv
!double precision, intent(in) :: x(dimv)
!double precision, intent(out):: g(dimv), f
!
!f = (x(1) - 1.)**4 + (x(2) - 1.)**4 + (x(3) -2.)**4;
!g(1) = 4 * (x(1) - 1)**3;
!g(2) = 4 * (x(2) - 1)**3;
!g(3) = 4 * (x(3) - 2)**3;
!end

subroutine linesearch_Armijo(n, x0, f0, g0, p, alpha, tol, nfl, info)
! n	the length of gradient vector
! f	value of cost function
! g	gradient vector
! p	new search direction, obtained by BFGS
! alpha input: initial step length
!	output: final step length
! tol	tolerance of termination
implicit none
integer :: n, l, i, j
double precision, intent(in) :: p(n), tol
double precision, intent(inout) :: alpha, f0, x0(n), g0(n) 
logical, intent(inout) :: info
integer, intent(out) :: nfl
double precision, allocatable :: x(:), g(:) 
double precision :: f, scal


allocate(x(n), g(n))

nfl = 0

! Line search
   scal=0.d+00
   DO j=1,n
      scal=scal+p(j)*g0(j)
   END DO
   !if(scal > 0) then
   !   write(*,*) 'not a descent direction!'
   !   stop
   !end if
!       IF (i > 2) THEN
!!          alpha=12.25d+00
!           alpha=1.d+00
!       ELSE
!           alpha=1.d+00
!       END IF

   info = .FALSE.

   DO l=1,20
   write(*,*) 'alpha = ',alpha
	DO j=1,n
           x(j)=x0(j)+alpha*p(j)
        END DO
!        n=1
!        DO m=2,iteration_max
!           DO j=1,dim_gust
!              u_gust(m,j)=perICvec(n)  
!              n=n+1            
!           END DO           
!        END DO
!        DO m=2,iteration_max
!           DO j=1,dim_gust
!              v_gust(m,j)=perICvec(n)  
!              n=n+1            
!           END DO           
!        END DO
!        DO m=1,nts
!           DO j=1,nvs
!              perIC(j,m)=perICvec(n)
!              n=n+1
!           END DO
!        END DO
!        CALL fvmls_s26NS_wg_var_2_obs_stoch(1,.FALSE.,perIC,perFF,HM_obs,.FALSE.,u_gust,v_gust)

        ! calculate cost function

	call Evaluate_f(n,x,f)        ! actually only need forward
	nfl = nfl + 1
        IF (f < f0+tol*alpha*scal) THEN
           info = .TRUE.
           f0 = f
	   x0 = x
	   !g0 = g
           EXIT    ! out the loop for l
        END IF
        alpha=alpha*0.33d+00    
        !alpha=alpha*0.5d+00    
   END DO
! end of line search
	f0 = f
	x0 = x
end


subroutine Test_cvsrch
	double precision x, f, g, s, stp, ftol, gtol, xtol, stpmin, stpmax, wa
	integer n, maxfev, info, nfev

	n = 1
	x = 0.
	stp = 1.e3
	ftol = 1.e-1
	gtol = 0.1
	xtol = 1.e-12
	stpmin = 0.
	stpmax = 1.e4
	maxfev = 100
	!call FUNCT(x,f,g,1)
	s = 1.0
	call cvsrch(n,x,f,g,s,stp,ftol,gtol,xtol,stpmin,stpmax,maxfev,info,nfev,wa)
	write(*,*) 'final step size',stp
	write(*,*) 'info',info
	write(*,*) 'nfev',nfev
	write(*,*) 'g',g


end


!
!     ----------------------------------------------------------------------
!     Line Search Routine
!     ----------------------------------------------------------------------
!
subroutine cvsrch(n,x,f,g,s,stp,ftol,gtol,xtol,		&
                       stpmin,stpmax,maxfev,info,nfev,wa)
      integer n,maxfev,info,nfev
      double precision f,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision x(n),g(n),s(n),wa(n)
!     in Nash & Nocedal's paper, 
!     ftol = 1.e-4
!     gtol = 0.9
!     stp = 1.0

!     **********
!
!     Subroutine cvsrch
!
!     The purpose of cvsrch is to find a step which satisfies 
!     a sufficient decrease condition and a curvature condition.
!     The user must provide a subroutine which calculates the
!     function and the gradient.
!
!     At each stage the subroutine updates an interval of
!     uncertainty with endpoints stx and sty. The interval of
!     uncertainty is initially chosen so that it contains a 
!     minimizer of the modified function
!
!          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
!
!     If a step is obtained for which the modified function 
!     has a nonpositive function value and nonnegative derivative, 
!     then the interval of uncertainty is chosen so that it 
!     contains a minimizer of f(x+stp*s).
!
!     The algorithm is designed to find a step which satisfies 
!     the sufficient decrease condition 
!
!           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
!
!     and the curvature condition
!
!           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
!
!     If ftol is less than gtol and if, for example, the function
!     is bounded below, then there is always a step which satisfies
!     both conditions. If no step can be found which satisfies both
!     conditions, then the algorithm usually stops when rounding
!     errors prevent further progress. In this case stp only 
!     satisfies the sufficient decrease condition.
!
!     The subroutine statement is
!
!        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
!                          stpmin,stpmax,maxfev,info,nfev,wa)
!     where
!
!	fcn is the name of the user-supplied subroutine which
!         calculates the function and the gradient.  fcn must 
!      	  be declared in an external statement in the user 
!         calling program, and should be written as follows.
!
!	  subroutine fcn(n,x,f,g)
!         integer n
!         f
!         x(n),g(n)
!	  ----------
!         Calculate the function at x and
!         return this value in the variable f.
!         Calculate the gradient at x and
!         return this vector in g.
!	  ----------
!	  return
!	  end
!
!       n is a positive integer input variable set to the number
!	  of variables.
!
!	x is an array of length n. On input it must contain the
!	  base point for the line search. On output it contains 
!         x + stp*s.
!
!	f is a variable. On input it must contain the value of f
!         at x. On output it contains the value of f at x + stp*s.
!
!	g is an array of length n. On input it must contain the
!         gradient of f at x. On output it contains the gradient
!         of f at x + stp*s.
!
!	s is an input array of length n which specifies the
!         search direction.
!
!	stp is a nonnegative variable. On input stp contains an
!         initial estimate of a satisfactory step. On output
!         stp contains the final estimate.
!
!       ftol and gtol are nonnegative input variables. Termination
!         occurs when the sufficient decrease condition and the
!         directional derivative condition are satisfied.
!
!	xtol is a nonnegative input variable. Termination occurs
!         when the relative width of the interval of uncertainty 
!	  is at most xtol.
!
!	stpmin and stpmax are nonnegative input variables which 
!	  specify lower and upper bounds for the step.
!
!	maxfev is a positive integer input variable. Termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!	info is an integer output variable set as follows:
!	  
!	  info = 0  Improper input parameters.
!
!	  info = 1  The sufficient decrease condition and the
!                   directional derivative condition hold.
!
!	  info = 2  Relative width of the interval of uncertainty
!		    is at most xtol.
!
!	  info = 3  Number of calls to fcn has reached maxfev.
!
!	  info = 4  The step is at the lower bound stpmin.
!
!	  info = 5  The step is at the upper bound stpmax.
!
!	  info = 6  Rounding errors prevent further progress.
!                   There may not be a step which satisfies the
!                   sufficient decrease and curvature conditions.
!                   Tolerances may be too small.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!	wa is a work array of length n.
!
!     Subprograms called
!
!	user-supplied......fcn
!
!	MINPACK-supplied...cstep
!
!	FORTRAN-supplied...abs,max,min
!	  
!     Argonne National Laboratory. MINPACK Project. June 1983
!     Jorge J. More', David J. Thuente
!
!     **********
      integer infoc,j,finish
      logical brackt,stage1
      double precision dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,	&
            finit,f0,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,	&
            stmin,stmax,width,width1,xtrapf
      data p5,p66,xtrapf/0.5,0.66,4.0/	! this means p5=0.5, p66=-.66, xtrapf=4.0
      info = 0
      infoc = 1

!
!     Check the input parameters for errors.
!


      if (n .le. 0 .or. stp .le. 0.0 .or. ftol .lt. 0.0 .or. 	&
         gtol .lt. 0.0 .or. xtol .lt. 0.0 .or. stpmin .lt. 0.0 	&
         .or. stpmax .lt. stpmin .or. maxfev .le. 0) return

!
!     Compute the initial gradient in the search direction
!     and check that s is a descent direction.
!
      dginit = 0.0
      do j = 1, n
         dginit = dginit + g(j)*s(j)
      end do

      if (dginit .ge. 0.0) then
         write(*,*) 'search direction is not a descent direction!' !MZW
         return
      endif

!
!     Initialize local variables.
!
      brackt = .false.
      stage1 = .true.
      nfev = 0
      finit = f
      dgtest = ftol*dginit
      width = stpmax - stpmin
      width1 = width/0.5d0
      do j = 1, n
         wa(j) = x(j)
      end do
!
!     The variables stx, fx, dgx contain the values of the step, 
!     function, and directional derivative at the best step.
!     The variables sty, fy, dgy contain the value of the step,
!     function, and derivative at the other endpoint of
!     the interval of uncertainty.
!     The variables stp, f, dg contain the values of the step,
!     function, and derivative at the current step.
!
      stx = 0.0
      fx = finit
      dgx = dginit
      sty = 0.0
      fy = finit
      dgy = dginit
!
!     Start of iteration.
!
   30 continue
!
!        Set the minimum and maximum steps to correspond
!        to the present interval of uncertainty.
!
         if (brackt) then
            stmin = min(stx,sty)
            stmax = max(stx,sty)
         else
            stmin = stx
            stmax = stp + xtrapf*(stp - stx)
         endif
!
!        Force the step to be within the bounds stpmax and stpmin.
!
         stp = max(stp,stpmin)
         stp = min(stp,stpmax)
!
!        If an unusual termination is to occur then let 
!        stp be the lowest point obtained so far.
!
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))	&
           .or. nfev .ge. maxfev-1 .or. infoc .eq. 0			&
           .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx
!
!        Evaluate the function and gradient at stp
!        and compute the directional derivative.
!

!        this do while loop is added by MZW
!        when forward simulation blows up, use back-tracking procedure
!        to get a smaller iniital step size         
         write(*,*) '--------------entering cvsrch---------------'
         f = finit
         !f0= f
         stp = 10*stp
         finish = 0
         !do while(abs(f - finit) < 1.e-14 .or. abs(f-f0)<1.e-14 .or. (f.ne.f) ) 
         do while(finish .ne. 1  .or. (f .ne. f))
         ! forward is not finished or cost function is NaN
                 stp = 0.1*stp
                 do j = 1, n
                    x(j) = wa(j) + stp*s(j)
                 end do
                 write(*,*) 'nfev=', nfev 
                 write(*,*) 'stp=', stp
                 !f0 = f
        	 call Evaluate_f(n,x,f)
                 open(101,file='finish_for')
                 read(101,*) finish
                 close(101)
                 nfev = nfev + 1
                 if(nfev .ge. maxfev+1) return
         end do
         call Evaluate_g(n,g)
         write(*,*) 'Adjusted step size =', stp

         dg = 0.0
         do j = 1, n
            dg = dg + g(j)*s(j)
         end do
         ftest1 = finit + stp*dgtest
!
!        Test for convergence.
!
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))	&
           .or. infoc .eq. 0) info = 6
         if (stp .eq. stpmax .and. 					&
            f .le. ftest1 .and. dg .le. dgtest) info = 5
         if (stp .eq. stpmin .and. 					&
            (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4
         if (nfev .ge. maxfev) info = 3
         if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2
         if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit)) info = 1

!
!        check which condition is not satisfied (MZW)
!
         if(f .gt. ftest1) write(*,*) 'Armijo condition is not satisfied!'
         if(abs(dg) .gt. gtol*(-dginit)) write(*,*) 'strong Wolfe condition is not satisfied!'

!
!        Check for termination.
!
         if (info .ne. 0) return
!
!        In the first stage we seek a step for which the modified
!        function has a nonpositive value and nonnegative derivative.
!
         if (stage1 .and. f .le. ftest1 .and. 				&
            dg .ge. min(ftol,gtol)*dginit) stage1 = .false.
!
!        A modified function is used to predict the step only if
!        we have not obtained a step for which the modified
!        function has a nonpositive function value and nonnegative 
!        derivative, and if a lower function value has been  
!        obtained but the decrease is not sufficient.
!
         if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
!
!           Define the modified function and derivative values.
!
            fm = f - stp*dgtest
            fxm = fx - stx*dgtest
            fym = fy - sty*dgtest
            dgm = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest
! 
!           Call cstep to update the interval of uncertainty 
!           and to compute the new step.
!
            call cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,		&
                      brackt,stmin,stmax,infoc)
!
!           Reset the function and gradient values for f.
!
            fx = fxm + stx*dgtest
            fy = fym + sty*dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
         else
! 
!           Call cstep to update the interval of uncertainty 
!           and to compute the new step.
!
            call cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,			&
                      brackt,stmin,stmax,infoc)
            end if
!
!        Force a sufficient decrease in the size of the
!        interval of uncertainty.
!
         if (brackt) then
            if (abs(sty-stx) .ge. p66*width1) 				&
              stp = stx + p5*(sty - stx)
            width1 = width
            width = abs(sty-stx)
            end if
!
!        End of iteration.
!
         go to 30
!
!     Last card of subroutine cvsrch.
!
      end

subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,			&
                      stpmin,stpmax,info)
      integer info
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
      logical brackt,bound
!     **********
!
!     Subroutine cstep
!
!     The purpose of cstep is to compute a safeguarded step for
!     a linesearch and to update an interval of uncertainty for
!     a minimizer of the function.
!
!     The parameter stx contains the step with the least function
!     value. The parameter stp contains the current step. It is
!     assumed that the derivative at stx is negative in the
!     direction of the step. If brackt is set true then a
!     minimizer has been bracketed in an interval of uncertainty
!     with endpoints stx and sty.
!
!     The subroutine statement is
!
!       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
!                        stpmin,stpmax,info)
! 
!     where
!
!       stx, fx, and dx are variables which specify the step,
!         the function, and the derivative at the best step obtained
!         so far. The derivative must be negative in the direction
!         of the step, that is, dx and stp-stx must have opposite 
!         signs. On output these parameters are updated appropriately.
!
!       sty, fy, and dy are variables which specify the step,
!         the function, and the derivative at the other endpoint of
!         the interval of uncertainty. On output these parameters are 
!         updated appropriately.
!
!       stp, fp, and dp are variables which specify the step,
!         the function, and the derivative at the current step.
!         If brackt is set true then on input stp must be
!         between stx and sty. On output stp is set to the new step.
!
!       brackt is a logical variable which specifies if a minimizer
!         has been bracketed. If the minimizer has not been bracketed
!         then on input brackt must be set false. If the minimizer
!         is bracketed then on output brackt is set true.
!
!       stpmin and stpmax are input variables which specify lower 
!         and upper bounds for the step.
!
!       info is an integer output variable set as follows:
!         If info = 1,2,3,4,5, then the step has been computed
!         according to one of the five cases below. Otherwise
!         info = 0, and this indicates improper input parameters.
!
!     Subprograms called
!
!       FORTRAN-supplied ... abs,max,min,sqrt
!                        ... dble
!
!     Argonne National Laboratory. MINPACK Project. June 1983
!     Jorge J. More', David J. Thuente
!
!     **********
      double precision gamma,p,p66,q,r,s,sgnd,stpc,stpf,stpq,theta
      data p66 /0.66/
      info = 0
!
!     Check the input parameters for errors.
!
      if ((brackt .and. (stp .le. min(stx,sty) .or.		&
         stp .ge. max(stx,sty))) .or.				&
         dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return
!
!     Determine if the derivatives have opposite sign.
!
      sgnd = dp*(dx/abs(dx))
!
!     First case. A higher function value.
!     The minimum is bracketed. If the cubic step is closer
!     to stx than the quadratic step, the cubic step is taken,
!     else the average of the cubic and quadratic steps is taken.
!
      if (fp .gt. fx) then
         info = 1
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
           stpf = stpc + (stpq - stpc)/2
           end if
         brackt = .true.
!
!     Second case. A lower function value and derivatives of
!     opposite sign. The minimum is bracketed. If the cubic
!     step is closer to stx than the quadratic (secant) step, 
!     the cubic step is taken, else the quadratic step is taken.
!
      else if (sgnd .lt. 0.0) then
         info = 2
         bound = .false.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
            end if
         brackt = .true.
!
!     Third case. A lower function value, derivatives of the
!     same sign, and the magnitude of the derivative decreases.
!     The cubic step is only used if the cubic tends to infinity 
!     in the direction of the step or if the minimum of the cubic
!     is beyond stp. Otherwise the cubic step is defined to be 
!     either stpmin or stpmax. The quadratic (secant) step is also 
!     computed and if the minimum is bracketed then the the step 
!     closest to stx is taken, else the step farthest away is taken.
!
      else if (abs(dp) .lt. abs(dx)) then
         info = 3
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
!
!        The case gamma = 0 only arises if the cubic does not tend
!        to infinity in the direction of the step.
!
         gamma = s*sqrt(max(0.,(theta/s)**2 - (dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. 0.0 .and. gamma .ne. 0.0) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
            end if
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (brackt) then
            if (abs(stp-stpc) .lt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
         else
            if (abs(stp-stpc) .gt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
            end if
!
!     Fourth case. A lower function value, derivatives of the
!     same sign, and the magnitude of the derivative does
!     not decrease. If the minimum is not bracketed, the step
!     is either stpmin or stpmax, else the cubic step is taken.
!
      else
         info = 4
         bound = .false.
         if (brackt) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
            end if
         end if
!
!     Update the interval of uncertainty. This update does not
!     depend on the new step or the case analysis above.
!
      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. 0.0) then
            sty = stx
            fy = fx
            dy = dx
            end if
         stx = stp
         fx = fp
         dx = dp
         end if
!
!     Compute the new step and safeguard it.
!
      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
      stp = stpf
      if (brackt .and. bound) then
         if (sty .gt. stx) then
            stp = min(stx+p66*(sty-stx),stp)
         else
            stp = max(stx+p66*(sty-stx),stp)
            end if
         end if
      return
!
!     Last card of subroutine cstep.
!
end

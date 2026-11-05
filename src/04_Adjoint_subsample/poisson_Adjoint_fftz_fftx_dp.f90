! MZW 02/21/17
! Test_fftz_fftx()
! poisson_Adjoint_fftz_fftx()
!　diverg_u_Adjoint()
! calculate_GPxi_Adjoint()
! gradx_GPxi_Adjoint()

subroutine Test_fftz_fftx()
	!cccccccccccccccccccccccccc
	! output fft results 
	!cccccccccccccccccccccccccc
	use poisson_fft
        ka = 1
        kb = ngzm
        ia = i1_T
        ib = i2_T
        ja = j1_T
        jb = j2_T

	! here use divergence of u_adj as source term, just for verifying the code
        qr(ka:kb, ia:ib, ja:jb)= u_adj(ka:kb    , ia+1:ib+1, ja:jb    )-u_adj(ka:kb, ia:ib, ja:jb)
        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
                                +v_adj(ka:kb    , ia:ib    , ja+1:jb+1)-v_adj(ka:kb, ia:ib, ja:jb)
        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
                                +w_adj(ka+1:kb+1, ia:ib    , ja:jb    )-w_adj(ka:kb, ia:ib, ja:jb)
        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    ) / dt
	
	write(*,*) '--------------fftz input---------------'
	write(*,*) 'At i=ia=', ia, 'j=ja=', ja
	write(*,*) qr(ka:kb,ia,ja)
	! just calculate once, no loop
	do j = j1_T,j1_T
	do i = i1_T,i1_T
		array_fftz(0:ngzm-1) = qr(1:ngzm,i,j) 

		call realFFTM(array_fftz,ngzm,1,+1)
		ka = ngzm/2 + 1
		kb = 2*(ngzm/2)
		qr(1:ka   , i, j) = array_fftz(0:ngzm:2)
		qr(ka+1:kb, i, j) = array_fftz(3:kb-1:2)
	end do
	end do
	write(*,*) '--------------fftz output---------------'
	write(*,*) 'At i=ia=', ia, 'j=ja=', ja
	write(*,*) array_fftz

	!-----------------------test inverse fft-----------------------
	do j = j1_T,j1_T
	do i = i1_T,i1_T

		call realFFTM(array_fftz,ngzm,1,-1)
		ka = ngzm/2 + 1
		kb = 2*(ngzm/2)
		qr(1:ka   , i, j) = array_fftz(0:ngzm:2)
		qr(ka+1:kb, i, j) = array_fftz(3:kb-1:2)
	end do
	end do
	write(*,*) '--------------inverse fftz output---------------'
	write(*,*) 'At i=ia=', ia, 'j=ja=', ja
	write(*,*) array_fftz

end

!--------------------------------------------------------------------
subroutine poisson_Adjoint_fftz_fftx()
!
! Note: phi and p solved from this subroutine is actually phi_adj * dt and p_adj * dt
!

	!cccccccccccccccccccccccccc
	! solve adjoint poisson equation
	!cccccccccccccccccccccccccc
	use poisson_fft
	real    :: cpu_poisson_total

        cpu_poisson_total=0.
	array_fftz=0.
	coeff_fftz=0.

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!    evaluate the source term using velocity gradient
	!    Source term is nothing but the Divergence of Ust  (D u*)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_source
        ka = 1
        kb = ngzm
        ia = i1_T
        ib = i2_T
        ja = j1_T
        jb = j2_T

!	-------------Commented by MZW, source term of forward Poisson solver-------------
!        qr(ka:kb, ia:ib, ja:jb)= ust(ka:kb    , ia+1:ib+1, ja:jb    )-ust(ka:kb, ia:ib, ja:jb)
!        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
!                                +vst(ka:kb    , ia:ib    , ja+1:jb+1)-vst(ka:kb, ia:ib, ja:jb)
!        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
!                                +wst(ka+1:kb+1, ia:ib    , ja:jb    )-wst(ka:kb, ia:ib, ja:jb)
!        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    ) / dt

	call readInt("iPrj", iPrj)
	if (iPrj.eq.1) then
		GPxi_adj = 0.
		call calculate_GPxi_Adjoint()
		if (irank.eq.1) then
			write(6,*),'MASS FLOW RATE FIXED'
		end if
	else
		if(irank==1) write(*,*) 'fixed pressure gradient !'
		GPxi_adj = 0.
	end if

	!-----------calculate source term------(eq 6.109)-----
	call diverg_u_Adjoint()
	!write(*,*) qr(ka,ia,ja:jb)
	! p_adj is stored as pressure * dt
	qr(ka:kb,  ia:ib, ja:jb)= qr(ka:kb,  ia:ib, ja:jb)*dt - p_adj(ka:kb,  ia:ib, ja:jb)/dt
	! GPxi_adj is stored as GPxi_adj * dx / V (eq 6.113)
	call gradx_GPxi_Adjoint()
!
!
!	after = realClock() ! cpu_source
!	cpu_source= (after - before) - t_over

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Fourier transform the source term q(i,j,k) in z-direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_fftz1
	denom1 = 1./ float(ngzm)
	do j = j1_T,j2_T
	do i = i1_T,i2_T
		array_fftz(0:ngzm-1) = qr(1:ngzm,i,j) 

		call realFFTM(array_fftz,ngzm,1,+1)
		ka = ngzm/2 + 1
		kb = 2*(ngzm/2)
		qr(1:ka   , i, j) = array_fftz(0:ngzm:2)
		!qr(ka+1:kb, i, j) = array_fftz(3:kb-1:2)    !MZW
		!------------MZW adjoint fft start-------------
		qr(ka+1:kb, i, j) = - array_fftz(3:kb-1:2)
		qr(1:kb, i, j) = qr(1:kb, i, j) * ngzm
		!------------MZW adjoint fft end---------------

		!J------------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fftz)
		!J	call dzfft1du(-1,ngzm,array_fftz,1,coeff_fftz)
		!J	!-----------------------------------------------
		!J	!   First 'half' of (qr) has the real values.
		!J	!   Second 'half' has the imaginary component.   
		!J	!-----------------------------------------------
		!J	ka = ngzm/2 + 1
		!J	kb = 2*(ngzm/2)
		!J	qr(1:ka   , i, j) = array_fftz(0:ngzm:2)*denom1
		!J	qr(ka+1:kb, i, j) = array_fftz(3:kb-1:2)*denom1
		!J------------------------------------------------------
	end do
	end do

	after  = realClock() ! cpu_fftz1
	cpu_fftz1= (after - before) - t_over
	before = realClock() ! cpu_tranpose

	!------------------------------------
	!   Time for Transposing qr in qT
	!------------------------------------
	allocate(qT(0:ngx, 0:ngy, nlzm))
	call transpose_qr_qT(qr,qT)

	if (irank.eq.iroot) then
		tmp = sum(qT(1:ngxm,1:ngym,1));
		!print*,'Total sum of RHS for (k=0) = ', iblock,jblock,tmp
	end if
	after  = realClock() ! cpu_transpose
	cpu_transpose= (after - before) - t_over
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Fast Fourier transform of qT in the x direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_fftx1
        ka = kTmin_2(irank)
        kb = kTmax_2(irank)
	do k = ka,kb
	kk = k-ka+1
	do j = 1,ngym
	        array_fftx=0.0
		array_fftx(0:ngxm-1) = qT(1:ngxm,j,kk)	
	
		call realFFTM(array_fftx,ngxm,1,+1) ! See SinCosFFT.f90 in lib directory
		ia = ngxm/2 + 1
		ib = 2*(ngxm/2)
		qT(  1:ia, j,kk) = array_fftx(0:ngxm:2)
		!qT(ia+1:ib,j,kk) = array_fftx(3:ib-1:2)    !MZW

		!------------MZW adjoint fft start-------------
		qT(ia+1:ib,j,kk) = - array_fftx(3:ib-1:2)
		qT(1:ib, j, kk) = qT(1:ib, j, kk) * ngxm
		!------------MZW adjoint fft end---------------
		
	end do
	end do
	after  = realClock() ! cpu_fftx1
	cpu_fftx1= (after - before) - t_over
	before = realClock() ! cpu_helmholtz

!------------------------------------------------------------------------------------------
!			START SOLUTION LOOP OVER WAVE NUMBERS (Kz)
!------------------------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space 
!       real part first (k < nzm/2+1)
!cccccccccccccccccccccccccccccccccccccccccccccc
ka = kTmin_2(irank)
kb = kTmax_2(irank)

do 20 k = ka,kb
	kk = k-ka+1
	allocate(z_fft(memlenap_z))
       
	if(k.le.ngzm/2+1) then
		goto 8000
	else
		ki=k-ngzm/2
		goto 8001
	end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space - real part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8000   nc=0
        !--------------------------
        !  Set up right-hand side
        !--------------------------
        do j=0,ngy
        do i=0,ngx
                ij=i+(ngxm+2)*j
                z_fft(irhs_mt(ngrid)+ij)=qT(i,j,kk) 	! irhs_mt(ngrid)=(ngxm+2)(ngym+2)
                z_fft(ip_mt  (ngrid)+ij)=0.0	! ip_mt(ngrid)=1
        end do
        end do
	nnx=ngxm
	nny=ngym

	!--------------------------------
	!  Solve the equation using tdma
	!--------------------------------
	mm=ngrid
	if(k==1) then
	call triDiagonal_y_k0(z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,             &
			      zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			      zap_mt(ivp_mt (mm)),vkz(k))
	else
	call triDiagonal_y(   z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,k,           &
			      zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			      zap_mt(ivp_mt (mm)),vkz(k))
	end if

	!--------------------------------
	!  Recover "phatr" from "z_fft"
	!--------------------------------
	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_fft(ip_mt(ngrid)+ij)
	end do
	end do
	goto 21

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space - imaginary part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8001   nc=0
	!--------------------------------
	!  Set up right-hand side
	!--------------------------------
	do j=0,ngy 
	do i=0,ngx
		ij=i+(ngxm+2)*j 
                z_fft(irhs_mt(ngrid)+ij)=qT(i,j,kk)
                z_fft(ip_mt  (ngrid)+ij)=0.0
	end do
	end do 
	nnx=ngxm
	nny=ngym
	!--------------------------------
	!  Solve the equation using tdma
	!--------------------------------
	mm=ngrid
	call triDiagonal_y(z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,ki,          &
			   zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			   zap_mt(ivp_mt (mm)),vkz(ki))

	!-------------------------------
	!  Recover "phatr" from "z_fft"
	!-------------------------------
	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_fft(ip_mt(ngrid)+ij)
	end do
	end do

 21   continue
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	deallocate(z_fft)
 20   continue
!------------------------------------------------------------------------------------------
!				FINISHED ALL WAVENUMBERS OF KZ
!------------------------------------------------------------------------------------------
	after  = realClock() ! cpu_helmholtz
	cpu_helmholtz= (after - before) - t_over
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Inverse fast fourier transform of phatr in the x direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_fftx2
        ka = kTmin_2(irank)
        kb = kTmax_2(irank)
	do k = ka,kb
	kk = k-ka+1
	do j = 1,ngym
	        array_fftx=0.0

		ia = ngxm/2 + 1
		ib = 2*(ngxm/2)
		array_fftx(0:ngxm:2) = phatr( 1:ia,  j,kk)
		!array_fftx(3:ib-1:2) = phatr(ia+1:ib,j,kk)     !MZW
		array_fftx(3:ib-1:2) = - phatr(ia+1:ib,j,kk)	!MZW adjoint ifft
		call realFFTM(array_fftx,ngxm,1,-1) ! See SinCosFFT.f90 in lib directory
	        phatr(1:ngxm,j,kk)=array_fftx(0:ngxm-1)
		phatr(1:ngxm,j,kk)= phatr(1:ngxm,j,kk) / real(ngxm) !MZW adjoint ifft

	end do
	end do
	after  = realClock() ! cpu_fftx2
	cpu_fftx2= (after - before) - t_over
	before = realClock() ! cpu_transpose

	!------------------------------------
	!  Time for Transposing phatr in p 
	!------------------------------------
	deallocate(qT)

	call transpose_phat_p(phatr,phi_adj)
	after  = realClock() ! cpu_transpose
	cpu_tranpose= cpu_transpose + (after - before) - t_over
!	if (irank.eq.1) write(6,*) 'Transpose time =', cpu_transpose
	before = realClock() ! cpu_fftz2
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! Fourier transform back to physical space, z-dir first
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do j = j1_T-1, j2_T+1
	do i = i1_T-1, i2_T+1
 		ka = ngzm/2 + 1
 		kb = 2*(ngzm/2)

		array_fftz = 0.
		array_fftz(0:ngzm:2) = phi_adj(1:ka   , i,j)
		!array_fftz(3:kb-1:2) = dp(ka+1:kb, i,j)    !MZW
		array_fftz(3:kb-1:2) = - phi_adj(ka+1:kb, i,j)    !MZW adjoint ifft
		!for2D  array_fftz(1)  = 0.
		!for2D  array_fftz(kb+1) = 0.
	
		call realFFTM(array_fftz,ngzm,1,-1)
		!J---------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fftz)
		!J	call zdfft1du(+1,ngzm,array_fftz,1,coeff_fftz)
		!J---------------------------------------------------
		phi_adj(1:ngzm,i,j) = array_fftz(0:ngzm-1)*dt
		phi_adj(1:ngzm,i,j) = phi_adj(1:ngzm,i,j) / real(ngzm) !MZW adjoint ifft
	end do
	end do

	!-----------------------MZW check adjoint----------------------
!	if (irank == iroot) then
!		write(*,*) '----------print part of ust_adj-----------'
!		write(*,*) ust_adj(:,i1_T,j1_T)
!		write(*,*) '----------print adjoint poisson solver-----------'	
!		write(*,*) ' At i=i1_T, j=j1_T'
!		write(*,*) phi_adj(:,i1_T,j1_T)
!	end if



	!-----------------------MZW check adjoint----------------------

	after = realClock() ! cpu_fftz2
	cpu_fftz2= (after - before) - t_over
	before= realClock() ! cpu_period

        !cccccccccccccccccccccccccccccccccccccccccccc
        ! TJ: make pressure periodic in x directions
        !cccccccccccccccccccccccccccccccccccccccccccc
        p_left_halo  = 0.
        p_right_halo = 0.

        if (iblock.eq.1 .or. iblock.eq.npx) then
	        if (npx.eq.1) then
       	        	p_right_halo(:,:) = phi_adj(:,  1,   :)
       	        	p_left_halo(:,:)  = phi_adj(:,nlxb-1,:)
        	else
                	if (iblock.eq.1) then
                        	p_right_halo(:,:) = phi_adj(:,  1,  :)
                        	call poisson_xperiodic(p_right_halo(0,0),p_left_halo(0,0),((ngzm+2)*(nly+1)))
                	else if (iblock.eq.npx) then
                        	p_left_halo(:,:) = phi_adj(:,nlxb-1,:)
                        	call poisson_xperiodic(p_left_halo(0,0),p_right_halo(0,0),((ngzm+2)*(nly+1)))
                	end if
        	end if
	end if

	if (iblock.eq.1) then
		phi_adj(:,  0, :) = p_left_halo(:,:)
	end if
	if (iblock.eq.npx) then
		phi_adj(:,nlxb,:) = p_right_halo(:,:)
	end if

	!ccccccccccccccccccccccccccccccccccccccccc
	! make pressure periodic in z directions
	!ccccccccccccccccccccccccccccccccccccccccc
	phi_adj(ngz,:,:)=phi_adj(1   ,:,:)
	phi_adj(0  ,:,:)=phi_adj(ngzm,:,:)

    !p_adj(:,:,:) = p_adj(:,:,:) + phi_adj(:,:,:)  no this relation for adjoint

	after = realClock() ! cpu_period
	cpu_period= (after - before) - t_over
	cpu_poisson_total=cpu_source+cpu_fftz1+cpu_transpose+cpu_fftx1+cpu_helmholtz+cpu_fftx2+cpu_fftz2+cpu_period

	! final check for adjoint poisson solver MZW
	!if (norm2(dp - phi_adj)>1e-10) then
	!	write(*,*) 'irank=',irank,': phi_adj is not correct !!!'
	!end if
	!write(*,*) 'complete checking phi_adj'
	!write(*,*) '----------------adjoint--------------'
	!write(*,*) phi_adj
	return
end

subroutine diverg_u_Adjoint()
	use data_export	

	!qr = 0.0

	!------------------------du/dx---------------------------
	do jj=j1_T,j2_T
	do ii=i1_T,i2_T
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!----------------------------------------------------------------
        	! Construct flags to tell begining and end of local domain
        	!----------------------------------------------------------------
        	! iflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )
  	     	! iflag_st3= [ 0,0,1,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )	
        	! iflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= ii <= nlxb-1 )
        	! iflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= ii <= nlxb-1 )
        	! jflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
		! jflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= jj <= nlyb-1 )
        	!-----------------------------------------------------------------------
        	iflag_st2 = 2/ii-2*max(2-ii,0)               ! Only equal to one when (ii=2)	! not used here
		iflag_st3 = 3/ii-max(3-ii,0)-max(2-ii,0)     ! Only equal to one when (ii=3)	! not used here
        	iflag_en1 = ii/(nlxb-1)                      ! Only equal to one when (ii=nlxb-1)	! not used here
		iflag_en2 = ii/(nlxb-2)-max(ii-(nlxb-2),0)   ! Only equal to one when (ii=nlxb-2)	! not used here
        	jflag_st1 = 1/jj                             ! Only equal to one when (jj=1)
		jflag_st2 = 2/jj-2*max(2-jj,0)               ! Only equal to one when (jj=2)
        	jflag_en1 = jj/(nlyb-1)                      ! Only equal to one when (jj=nlyb-1)
		jflag_en2 = jj/(nlyb-2)-max(jj-(nlyb-2),0) ! Only equal to one when (jj=nlyb-2)

		gradx_adj= (1-jflag_st1-jflag_st2-jflag_en1-jflag_en2) * 	&
			  	 ( - (u_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (u_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (u_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  u_adj(k,ii,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (u_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  u_adj(k,ii,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (u_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  	&
			     -  u_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (u_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 		&
			     -  u_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (u_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  u_adj(k,ii+1,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (u_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  u_adj(k,ii+1,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (u_adj(k,ii+1  ,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  	&
			     -  u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (u_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 	&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
		! MZW: have considered that u_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_en1   * &	
				( - (u_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (u_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (u_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *	& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (u_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *	&
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	 				&
			     - (u_adj(k,ii+1,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *	& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (u_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *	&
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )				
		! MZW: have considered that u_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_st1   * &	
				( - (u_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (u_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     + (u_adj(k,ii  ,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *	& 
				(svetax(i-1,j+1  )+svetax(i,j+1  ))/2.*0.25				&
			     + (u_adj(k,ii  ,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *	& 
				(svetay(i-1,j+1  )+svetay(i,j+1  ))/2.*0.25				&
			     + (u_adj(k,ii+1 ,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *	& 
				(svetax(i,j+1  )+svetax(i+1,j+1  ))/2.*0.25				&
			     + (u_adj(k,ii+1 ,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *	& 
				(svetay(i,j+1  )+svetay(i+1,j+1  ))/2.*0.25   )

		gradx_adj= gradx_adj + jflag_en2  * &
		   	    ( - (u_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (u_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (u_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			      )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (u_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			      )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (u_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  	&
			     -  u_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (u_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 		&
			     -  u_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (u_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			      )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (u_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			      )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (u_adj(k,ii+1  ,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  	&
			     -  u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (u_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 	&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )

		gradx_adj= gradx_adj + jflag_st2  * &
		   	    ( - (u_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (u_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (u_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  u_adj(k,ii,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (u_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  u_adj(k,ii,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (  	&
			     -  u_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (	&
			     -  u_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (u_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  u_adj(k,ii+1,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (u_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  u_adj(k,ii+1,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - ( 	&
			     -  u_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (	&
			     -  u_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
				
		!-----------------------------------------------------------------------
		qr(k,ii,jj) = gradx_adj

	end do
	end do
	end do

	!------------------------dv/dy---------------------------
	do jj=j1_T,j2_T
	do ii=i1_T,i2_T
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!----------------------------------------------------------------
        	! Construct flags to tell begining and end of local domain
        	!----------------------------------------------------------------
        	! iflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )
  	     	! iflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )	
        	! iflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= ii <= nlxb-1 )
        	! iflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= ii <= nlxb-1 )
        	! jflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
		! jflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= jj <= nlyb-1 )
        	!-----------------------------------------------------------------------
        	iflag_st1 = 1/ii 	                     ! Only equal to one when (ii=1)
		iflag_st2 = 2/ii-2*max(2-ii,0)	             ! Only equal to one when (ii=2)
        	iflag_en1 = ii/(nlxb-1)                      ! Only equal to one when (ii=nlxb-1)
		iflag_en2 = ii/(nlxb-2)-max(ii-(nlxb-2),0)   ! Only equal to one when (ii=nlxb-2)
        	jflag_st1 = 1/jj                             ! Only equal to one when (jj=1)	! not used here
		jflag_st2 = 2/jj-2*max(2-jj,0)               ! Only equal to one when (jj=2)	! not used here
        	jflag_en1 = jj/(nlyb-1)                      ! Only equal to one when (jj=nlyb-1)	! not used here
		jflag_en2 = jj/(nlyb-2)-max(jj-(nlyb-2),0) ! Only equal to one when (jj=nlyb-2)		! not used here

		grady_adj=	(1-iflag_st1-iflag_st2-iflag_en1-iflag_en2) * 			&
				( - (v_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	v_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (v_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  v_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (v_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  v_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (v_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  v_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (v_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  v_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (v_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  v_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (v_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  v_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (v_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  v_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (v_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  v_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (v_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  v_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )
		!MZW: periodic boundary condition has been considered here, i2_T+1 -> i1_T
		! at inlet and outlet, assume the grid is orthogonal
		grady_adj=    grady_adj + iflag_st1 * 	&
				( - (v_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	v_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (v_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  v_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (				&
			     -  v_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (				&
			     -  v_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25	    				&
			     - (		 		&
			     -  v_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (		 		&
			     -  v_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25		)

		grady_adj=    grady_adj + iflag_en1 * 	&
				( - (v_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	v_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (v_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  v_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (v_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
				     )*				& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (v_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  	&
			     	     )*				& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (v_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     	     )*				& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (v_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     	     )*				& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )

		grady_adj=   grady_adj + iflag_st2 *   &
				( - (v_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	v_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (v_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  v_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (v_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  v_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (v_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  v_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (		  		&
			     -  v_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (			  	&
			     -  v_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (v_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  v_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (v_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  v_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (			 	&
			     -  v_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (			  	&
			     -  v_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )	

		grady_adj=   grady_adj + iflag_en2 *   &
				( - (v_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	v_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (v_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  v_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (v_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     	  )*				& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (v_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     	  )*				& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (v_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  v_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (v_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  v_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (v_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			    	  )*				& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (v_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     	  )*				& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (v_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  v_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (v_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  v_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )
			
		!-----------------------------------------------------------------------
		qr(k,ii,jj) = qr(k,ii,jj) + grady_adj

	end do
	end do
	end do

	!------------------------dw/dz---------------------------
	do jj=j1_T,j2_T
	do ii=i1_T,i2_T
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                !----------------------------------------------------------------
                ! Construct flags to tell begining and end of local domain
                !----------------------------------------------------------------
                ! iflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )
                ! iflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= ii <= nlxb-1 )
                ! jflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
                ! jflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= jj <= nlyb-1 )
                !-----------------------------------------------------------------------
		gradz_adj= (- (w_adj(k,ii,jj)*swz(i,j)/(vp(i,j))					&
			    -w_adj(k+1,ii,jj)*swz(i,j)/(vp(i,j))	)*				&
			    swz(i,j)		)							
		!-----------------------------------------------------------------------
		qr(k,ii,jj) = qr(k,ii,jj) + gradz_adj

	end do
	end do
	end do

end

subroutine calculate_GPxi_Adjoint()  ! (6.170) * dx / V

        use data_export
	use mesh_export

	integer :: flag_cost, my_jstart, my_jend, global_jndx, loc_jndx, jjout1
        real 	:: dxi2, volume  

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! TJ: Initialize GPxi, dVu, Vu
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        GPxi_adj = 0.0
        volume  = 0.0

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! TJ: Start integration here
	!	- integrate w.r.t. xi direction
	!	- summation in eta and z directions
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ka = 1
	kb = ngzm
	do ii=i1_u,i2_u	
	do jj=j1_u,j2_u
	do kk=ka,kb
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
                
		jflag_start = 1/jj              		! Only equal to one when (jj=1)
                jflag_end   = jj/(nlyb-1)       		! Only equal to one when (jj=nlyb-1)

                !-----------------------------------------------------------------------
                ! TJ: Integrate w.r.t. xi across interior points
                !-----------------------------------------------------------------------
  	        GPxi_adj  = GPxi_adj -  (1-jflag_start-jflag_end)*( u_adj(k,ii,jj)* suxix(i,j)   )
                !-----------------------------------------------------------------------
                ! TJ: Integrate w.r.t. xi at jj=nlyb-1
                !-----------------------------------------------------------------------
		GPxi_adj = GPxi_adj -       jflag_end*(   u_adj(k,ii,jj)* suxix(i,j)  )
                !-----------------------------------------------------------------------
                ! TJ: Integrate w.r.t. xi at jj=1
                !-----------------------------------------------------------------------
		GPxi_adj = GPxi_adj -      jflag_start*(  u_adj(k,ii,jj)* suxix(i,j) )
	end do
	end do
	end do
!	!------------------------observation----------------------------------
!	flag_cost = 0
!
!	sigma_inv2 = 1./sigma**2
!
!	my_jstart = jTmin_1(jblock)
!	my_jend   = jTmax_1(jblock)
!
!	jjout1 = 1
!	global_jndx = jjout1
!	loc_jndx    = jmap_1(global_jndx)
!	if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
!		flag_cost = 1
!	end if
!
!	! Here p is p*dt !!!!!!!!!!!!!!!!!!!!	
!	if (flag_cost == 1) then
!		do ii=i1_T,i2_T
!		do kk=1,ngzm
!			i = ibmap_1(ii) 
!			k = kk
!			GPxi_adj = GPxi_adj +  sigma_inv2* (p_tmp(k,ii,1) - p_obs(k,ii,1))/(dt**2)*(xm(i)-xL)
!		end do
!		end do
!	end if
!
!	jjout1 = ngy-1
!	global_jndx = jjout1
!	loc_jndx    = jmap_1(global_jndx)
!	if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
!		flag_cost = 2
!	end if
!	
!	if (flag_cost == 2) then	
!		do ii=i1_T,i2_T
!		do kk=1,ngzm
!			i = ibmap_1(ii) 
!			k = kk
!			GPxi_adj = GPxi_adj + sigma_inv2 * (p_tmp(k,ii,2) - p_obs(k,ii,2))/(dt**2)*(xm(i)-xL)
!		end do
!		end do
!	end if
!
!	flag_cost = 0

	!TJ: !MZW define the uniform xi increment here
	dxi         = (xpg(ngx,1)-xpg(1,1))/DBLE(ngxm)

	!TJ: compute volume here
!	volume = ( xpg(ngx,1)-xpg(1,1) )*( ypg(1,ngy)-ypg(1,1) )*( alz )	
	volume = ( xpg(ngx,1)-xpg(1,1) )*2.0*( alz )	 !JL: The "2.0" indicates Ly=2.0

	!TJ: sum GPxi across the processes
	call globalSum(GPxi_adj,1)

	!TJ: calculate global GPxi
	GPxi_adj = GPxi_adj * dxi /volume
	
!	if (irank.eq.1) then
!		write(6,*),'Adjoint MASS FLOW RATE FIXED: GPxi_adj =', GPxi_adj
!	end if

        return
end

subroutine gradx_GPxi_Adjoint()
	use data_export	

	do jj=j1_T,j2_T
	do ii=i1_T,i2_T
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!----------------------------------------------------------------
        	! Construct flags to tell begining and end of local domain
        	!----------------------------------------------------------------
        	! iflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )
  	     	! iflag_st3= [ 0,0,1,0,......,0,0,0 ]        ( 1<= ii <= nlxb-1 )	
        	! iflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= ii <= nlxb-1 )
        	! iflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= ii <= nlxb-1 )
        	! jflag_st1= [ 1,0,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
		! jflag_st2= [ 0,1,0,0,......,0,0,0 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en1= [ 0,0,0,0,......,0,0,1 ]        ( 1<= jj <= nlyb-1 )
        	! jflag_en2= [ 0,0,0,0,......,0,1,0 ]        ( 1<= jj <= nlyb-1 )
        	!-----------------------------------------------------------------------
        	iflag_st2 = 2/ii-2*max(2-ii,0)               ! Only equal to one when (ii=2)	! not used here
		iflag_st3 = 3/ii-max(3-ii,0)-max(2-ii,0)     ! Only equal to one when (ii=3)	! not used here
        	iflag_en1 = ii/(nlxb-1)                      ! Only equal to one when (ii=nlxb-1)	! not used here
		iflag_en2 = ii/(nlxb-2)-max(ii-(nlxb-2),0)   ! Only equal to one when (ii=nlxb-2)	! not used here
        	jflag_st1 = 1/jj                             ! Only equal to one when (jj=1)
		jflag_st2 = 2/jj-2*max(2-jj,0)               ! Only equal to one when (jj=2)
        	jflag_en1 = jj/(nlyb-1)                      ! Only equal to one when (jj=nlyb-1)
		jflag_en2 = jj/(nlyb-2)-max(jj-(nlyb-2),0) ! Only equal to one when (jj=nlyb-2)

		gradx_adj= (1-jflag_st1-jflag_st2-jflag_en1-jflag_en2) * 	&
			  	 ( - (suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  				&
			     -	suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				&	 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		     		&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			     -  suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*				& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			     -  suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*				& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  				&
			     -  suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 				&
			     -  suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			     -  suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*			& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			     -  suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*			& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  				&
			     -  suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
		! MZW: have considered that u_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_en1   * &	
				( - (suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  				&
			     -	suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *				& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *				&
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	 				&
			     - (suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *			& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *			&
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )				
		! MZW: have considered that u_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_st1   * &	
				( - (suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  				&
			     -	suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     + (suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *				& 
				(svetax(i-1,j+1  )+svetax(i,j+1  ))/2.*0.25				&
			     + (suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *				& 
				(svetay(i-1,j+1  )+svetay(i,j+1  ))/2.*0.25				&
			     + (suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *			& 
				(svetax(i,j+1  )+svetax(i+1,j+1  ))/2.*0.25				&
			     + (suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *			& 
				(svetay(i,j+1  )+svetay(i+1,j+1  ))/2.*0.25   )

		gradx_adj= gradx_adj + jflag_en2  * &
		   	    ( - (suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  				&
			     -	suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			      )*									& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			      )*									& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  				&
			     -  suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 				&
			     -  suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			      )*									&	 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			      )*									& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  				&
			     -  suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )

		gradx_adj= gradx_adj + jflag_st2  * &
		   	    ( - (suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  				&
			     -	suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      				&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*				& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			     -  suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*				& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 				&
			     -  suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*				& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (  									&
			     -  suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (									&
			     -  suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*				& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			     -  suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*			& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 				&
			     -  suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*			& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - ( 									&
			     -  suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (									&
			     -  suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*				& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
				
		!-----------------------------------------------------------------------
		qr(k,ii,jj) = qr(k,ii,jj) + dt * GPxi_adj * gradx_adj

	end do
	end do
	end do
	!write(*,*) 'gradient of GPxi_adj:',gradx_adj
end

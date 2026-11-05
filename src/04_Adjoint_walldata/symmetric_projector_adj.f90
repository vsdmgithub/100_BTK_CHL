! MZW 05/07/18
! symmetric projector, which projects any ust onto divergence-free space
! subroutine sym_proj_poisson_init()
! subroutine sym_proj_poisson()
! subroutine sym_proj_update()
! subroutine sym_proj_preprocess_fft()
! subroutine sym_proj_pcoef_fft()
! subroutine sym_proj_triDiagonal_y
! subroutine sym_proj_triDiagonal_y_k0

!=======================================================================
subroutine sym_proj_poisson_init()
	use poisson_fft

	real    dx2, dz2, pi, twopi, arg
	integer i, k, m

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Compute and store the modified wavenumber (vkz)
	!  for central differencing (used in FFT of Spanwise
	!  direction --> Page 40)
	!  (This replaces "subroutine waveno" in xiahua's code)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	twopi=2.*2.*acos(0.)
	nzh=ngz/2+1
	dz2=(alz/float(ngz-1))**2
	c=twopi/float(ngzm)
	do k=1,ngzm
		if (k .lt. nzh) then
			m = k-1
		else
			m = k-ngzm-1
		end if
		arg = c*float(m)
		!vkz(k)=2.*(1.-cos(arg))/dz2    ! MZW: I modified this to the next line
		vkz(k)=2.*(1-cos(arg))
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Compute and store the modified wavenumber (vkx)
	!  for central differencing (used in Fast Fourier 
	!  Transform(FFT) in the streamwise direction)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	twopi=2.*2.*acos(0.)
	nxh=ngx/2+1
	!write(6,*) 'xpg(ngx,1)-xpg(1,1) = ',xpg(ngx,1)-xpg(1,1)
	dx2=((xpg(ngx,1)-xpg(1,1))/float(ngxm))**2
	c=twopi/float(ngxm)
	do i=1,ngxm
		if (i .lt. nxh) then
			m = i-1
		else
			m = i-ngxm-1
		end if
		arg = c*float(m)
		!vkx(i)=2.*(1.-cos(arg))/dx2     ! MZW: I modified this to the next line
		vkx(i)=2.*(1.-cos(arg))
	end do

        ngrid = 3 !JL (just "3" for 2D FFT pressure solver)

	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Calculate the coefficients of poisson equation
	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	call sym_proj_preprocess_fft()

	return
end

subroutine sym_proj_poisson_adj(cpu_poisson_total)
	!cccccccccccccccccccccccccc
	! solve poisson equation
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
        qr(ka:kb, ia:ib, ja:jb)= ust_adj(ka:kb    , ia+1:ib+1, ja:jb    )-ust_adj(ka:kb, ia:ib, ja:jb)
        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
                                +vst_adj(ka:kb    , ia:ib    , ja+1:jb+1)-vst_adj(ka:kb, ia:ib, ja:jb)
        qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    )  &
                                +wst_adj(ka+1:kb+1, ia:ib    , ja:jb    )-wst_adj(ka:kb, ia:ib, ja:jb)
        !qr(ka:kb, ia:ib, ja:jb)=  qr(ka:kb    , ia:ib    , ja:jb    ) / dt ! MZW: I don't need this
	qr(ka:kb, ia:ib, ja:jb)= - qr(ka:kb, ia:ib, ja:jb)                  ! MZW: I need this

	after = realClock() ! cpu_source
	cpu_source= (after - before) - t_over

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
		qr(ka+1:kb, i, j) = array_fftz(3:kb-1:2)

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
		print*,'Total sum of RHS for (k=0) = ', iblock,jblock,tmp
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
		qT(ia+1:ib,j,kk) = array_fftx(3:ib-1:2)
		
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
	if(k==1) then
	call sym_proj_triDiagonal_y_k0(z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,             &
			      zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			      zap_mt(ivp_mt (mm)),vkz(k))
	else
	call sym_proj_triDiagonal_y(   z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,k,           &
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
	call sym_proj_triDiagonal_y(z_fft (ip_mt  (mm)),z_fft (irhs_mt(mm)),nnx,nny,ki,          &
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
		array_fftx(3:ib-1:2) = phatr(ia+1:ib,j,kk)
		call realFFTM(array_fftx,ngxm,1,-1) ! See SinCosFFT.f90 in lib directory
	        phatr(1:ngxm,j,kk)=array_fftx(0:ngxm-1)

	end do
	end do
	after  = realClock() ! cpu_fftx2
	cpu_fftx2= (after - before) - t_over
	before = realClock() ! cpu_transpose

	!------------------------------------
	!  Time for Transposing phatr in p 
	!------------------------------------
	deallocate(qT)

	call transpose_phat_p(phatr,dp)
	after  = realClock() ! cpu_transpose
	cpu_tranpose= cpu_transpose + (after - before) - t_over
	if (irank.eq.1) write(6,*) 'Transpose time =', cpu_transpose
	before = realClock() ! cpu_fftz2

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! Fourier transform back to physical space, z-dir first
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do j = j1_T-1, j2_T+1
	do i = i1_T-1, i2_T+1
 		ka = ngzm/2 + 1
 		kb = 2*(ngzm/2)

		array_fftz = 0.
		array_fftz(0:ngzm:2) = dp(1:ka   , i,j)
		array_fftz(3:kb-1:2) = dp(ka+1:kb, i,j)
		!for2D  array_fftz(1)  = 0.
		!for2D  array_fftz(kb+1) = 0.
	
		call realFFTM(array_fftz,ngzm,1,-1)
		!J---------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fftz)
		!J	call zdfft1du(+1,ngzm,array_fftz,1,coeff_fftz)
		!J---------------------------------------------------
		!dp(1:ngzm,i,j) = array_fftz(0:ngzm-1)*dt
		! MZW: I don't need dp*dt
		dp(1:ngzm,i,j) = array_fftz(0:ngzm-1)
	end do
	end do

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
       	        	p_right_halo(:,:) = dp(:,  1,   :)
       	        	p_left_halo(:,:)  = dp(:,nlxb-1,:)
        	else
                	if (iblock.eq.1) then
                        	p_right_halo(:,:) = dp(:,  1,  :)
                        	call poisson_xperiodic(p_right_halo(0,0),p_left_halo(0,0),((ngzm+2)*(nly+1)))
                	else if (iblock.eq.npx) then
                        	p_left_halo(:,:) = dp(:,nlxb-1,:)
                        	call poisson_xperiodic(p_left_halo(0,0),p_right_halo(0,0),((ngzm+2)*(nly+1)))
                	end if
        	end if
	end if

	if (iblock.eq.1) then
		dp(:,  0, :) = p_left_halo(:,:)
	end if
	if (iblock.eq.npx) then
		dp(:,nlxb,:) = p_right_halo(:,:)
	end if

	!ccccccccccccccccccccccccccccccccccccccccc
	! make pressure periodic in z directions
	!ccccccccccccccccccccccccccccccccccccccccc
	dp(ngz,:,:)=dp(1   ,:,:)
	dp(0  ,:,:)=dp(ngzm,:,:)

    !p(:,:,:) = p(:,:,:) + dp(:,:,:) ! MZW

	after = realClock() ! cpu_period
	cpu_period= (after - before) - t_over
	cpu_poisson_total=cpu_source+cpu_fftz1+cpu_transpose+cpu_fftx1+cpu_helmholtz+cpu_fftx2+cpu_fftz2+cpu_period

	return
end

subroutine sym_proj_update_adj()
	use data_export

	do jj=j1_u,j2_u
	do ii=i1_u,i2_u
	do k=1,ngz-1
		i  = ibmap_1(ii)
		j  = jbmap_1(jj)

		u_adj(k,ii,jj) = ust_adj(k,ii,jj) - (  dp(k,ii-1,jj)/vp(i-1,j)-dp(k,ii,jj)/vp(i,j)  )

	end do
	end do
	end do

	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
	do k=1,ngz-1
		i  = ibmap_1(ii)
		j  = jbmap_1(jj)

		v_adj(k,ii,jj) = vst_adj(k,ii,jj) - (  dp(k,ii,jj-1)/vp(i,j-1)-dp(k,ii,jj)/vp(i,j)  )

	end do
	end do
	end do

	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do k=1,ngz-1
		i  = ibmap_1(ii)
		j  = jbmap_1(jj)

		w_adj(k,ii,jj) = wst_adj(k,ii,jj) - (  dp(k-1,ii,jj)/vp(i,j)  -dp(k,ii,jj)/vp(i,j)  )

	end do
	end do
	end do


end

subroutine sym_proj_preprocess_fft 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! allocate memory for Poisson coefficients
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export
	dimension zsuv(memlenap_zsuv)
	dimension isux_mt(ngrid),isuy_mt(ngrid),isvx_mt(ngrid),isvy_mt(ngrid)

	if (irank.eq.1) then
	!MZW	write(*,*) 'WARNING: fft int arrays: ngrid  ,nglevel_max  =',ngrid,nglevel_max 
	end if
	memap=0
	m    =ngrid
	nnx  =ngxm 
	nny  =ngym 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for coefficients in poisson equation
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	iaps_mt (m)=malocap_fft(nnx*nny        ,memap,memlenap_zap)
	iapn_mt (m)=malocap_fft(nnx*nny        ,memap,memlenap_zap)
	iapp_mt (m)=malocap_fft(nnx*nny        ,memap,memlenap_zap)
	ivp_mt  (m)=malocap_fft((nnx+2)*(nny+2),memap,memlenap_zap)
	iedge_mt(m)=malocap_fft(2              ,memap,memlenap_zap)
	if (irank.eq.1) then
	!MZW	write(*,*) 'WARNING: zap_fft       : memap+1,memlenap_zap =',memap+1,memlenap_zap 
	end if

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	! temporary pointer for metric
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	isux_mt(m)=malocap_fft((nnx+3)*(nny+2),memap,memlenap_zsuv)
	isuy_mt(m)=malocap_fft((nnx+3)*(nny+2),memap,memlenap_zsuv)
	isvx_mt(m)=malocap_fft((nnx+2)*(nny+3),memap,memlenap_zsuv)
	isvy_mt(m)=malocap_fft((nnx+2)*(nny+3),memap,memlenap_zsuv)
	if (irank.eq.1) then
	!MZW	write(*,*) 'WARNING: zsuv_fft      : memap+1,memlenap_zsuv=',memap+1,memlenap_zsuv
	end if

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for solution, source term
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	ip_mt  (m)=malocap_fft((nnx+2)*(nny+2),memap,memlenap_z)
	irhs_mt(m)=malocap_fft((nnx+2)*(nny+2),memap,memlenap_z)
	if (irank.eq.1) then
	!MZW	write(*,*) 'WARNING: z_fft         : memap+1,memlenap_z   =',memap+1,memlenap_z
	end if

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! assign area vector and volume
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call copy_fft  (zap_mt(ivp_mt (m)),vp(0:ngx,0:ngy)      ,nnx,nny)
	call copysu_fft(zsuv  (isux_mt(m)),suxix ,nnx,nny)
	call copysu_fft(zsuv  (isuy_mt(m)),suxiy ,nnx,nny)
	call copysv_fft(zsuv  (isvx_mt(m)),svetax(0:ngx,0:ngy+1),nnx,nny)
	call copysv_fft(zsuv  (isvy_mt(m)),svetay(0:ngx,0:ngy+1),nnx,nny)

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! evaluate Poisson coefficients
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call sym_proj_pcoef_fft(	zap_mt(iaps_mt (m)),zap_mt(iapn_mt (m)), &
			zap_mt(iapp_mt (m)),zap_mt(ivp_mt  (m)), &
			zsuv  (isux_mt (m)),zsuv  (isuy_mt (m)), &
			zsuv  (isvx_mt (m)),zsuv  (isvy_mt (m)), &
			zap_mt(iedge_mt(m)),nnx,nny) 

      return
end

subroutine sym_proj_pcoef_fft(aps_mt  ,apn_mt  , app_mt   ,vp_mt    , &
                        suxix_mt,suxiy_mt, svetax_mt,svetay_mt, &
                        aedgex  ,nnx,nny)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! compute the Poisson equation coefficients
	! with (nnx, nny) points
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real	aps_mt (nnx,nny), apn_mt(nnx,nny), app_mt(nnx,nny),            &
		suxix_mt (0:nnx+1+1,0:nny+1  ),suxiy_mt (0:nnx+1+1,0:nny+1  ), &
		svetax_mt(0:nnx+1  ,0:nny+1+1),svetay_mt(0:nnx+1  ,0:nny+1+1), &
		aedgex(2),vp_mt(0:nnx+1,0:nny+1)
	integer iedgex(2)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the leading edge and trailing edge i index
	! pressure side and suction side have same index
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	iedgex(1)=1
	iedgex(2)=nnx
	aedgex   =1.0*iedgex

	!=================================================
	!     Initialize all the Poisson Coefficients
	!=================================================
	aps_mt  = 0.0
	app_mt  = 0.0
	apn_mt  = 0.0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,n    (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 33 j=1,nny
      do 33 i=1,nnx
        if(j.eq.nny) then
          aps_mt(i,j)=aps_mt(i,j)
!----------
        else 
            aps_mt(i,j) =aps_mt(i,j)

            !vetan=(vp_mt(i,j)+vp_mt(i,j+1))/2.

            !a1=svetax_mt(i,j+1)
            !a2=svetay_mt(i,j+1)

            !b1=(svetax_mt(i,j+1)+svetax_mt(i,j+2))/2.
            !b2=(svetay_mt(i,j+1)+svetay_mt(i,j+2))/2.
            !apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vetan

            !b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            !b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            !app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vetan
            ! MZW
            apn_mt(i,j)=apn_mt(i,j)+1./vp_mt(i,j+1)
            app_mt(i,j)=app_mt(i,j)-1./vp_mt(i,j)
!--------
        end if         
 33   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,s      (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 34 j=1,nny
      do 34 i=1,nnx
        if(j.eq.1) then
          apn_mt(i,j)=apn_mt(i,j)
!-----
        else 
            apn_mt(i,j)=apn_mt(i,j)

            !vetas=-(vp_mt(i,j)+vp_mt(i,j-1))/2.

            !a1=svetax_mt(i,j)
            !a2=svetay_mt(i,j)

            !b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            !b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            !app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vetas

            !b1=(svetax_mt(i,j)+svetax_mt(i,j-1))/2.
            !b2=(svetay_mt(i,j)+svetay_mt(i,j-1))/2.
            !aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vetas
            ! MZW
            app_mt(i,j)=app_mt(i,j)-1./vp_mt(i,j)
            aps_mt(i,j)=aps_mt(i,j)+1./vp_mt(i,j-1)
!-------------
        end if         
 34   continue

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled equation by TDMA in the y direction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine sym_proj_triDiagonal_y(p_ta,rhs_ta,nnx,nny,k_ta, &
                       aps_ta,apn_ta,app_ta,vp_ta,vkz_ta)
	use data_export
        real, dimension(0:nnx+1,0:nny+1) :: rhs_ta,  p_ta, vp_ta
        real, dimension(  nnx  ,  nny  ) :: aps_ta,apn_ta,app_ta
        real, dimension(  nnx  ,  nny  ) ::   a_ta,  b_ta,  c_ta, d_ta
        integer :: nnx,nny,k_ta

        p_ta=0.0
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        do j=1,nny
        do i=1,nnx

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: distinguish real and imag modified
		!          wavnumbers in x-direction
		!ccccccccccccccccccccccccccccccccccccccccccccc
	
		if(i.le.ngxm/2+1) then
			ii = i
		else
			ii = i - ngxm/2 
		end if	

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: solve using TDMA
		!ccccccccccccccccccccccccccccccccccccccccccccc
		
                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                !b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(ii))
		! MZW
                b_ta(i,j)=-app_ta(i,j)+(vkz_ta+vkx(ii))/vp_ta(i,j)
                a_ta(i,j)=-aps_ta(i,j)
	
        end do
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny,1,nnx)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny
        do i=1,nnx
                p_ta(i,j)=d_ta(i,j)
        end do
        end do

	return
end
!--------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled equation by TDMA in the y direction
! for the zeroth wavenumber of Kz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine sym_proj_triDiagonal_y_k0(p_ta,rhs_ta,nnx,nny, &
                          aps_ta,apn_ta,app_ta,vp_ta,vkz_ta)
	use data_export
        real, dimension(0:nnx+1,0:nny+1) :: rhs_ta,  p_ta, vp_ta
        real, dimension(  nnx  ,  nny  ) :: aps_ta,apn_ta,app_ta
        real, dimension(  nnx  ,  nny  ) ::   a_ta,  b_ta,  c_ta, d_ta
        integer :: nnx,nny

        p_ta=0.0
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! special treatment to remove singularity when i=k=1
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        i=1
        do j=1,nny-1
                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                !b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(i))
		! MZW
                b_ta(i,j)=-app_ta(i,j)+(vkz_ta+vkx(i))/vp_ta(i,j)
                a_ta(i,j)=-aps_ta(i,j)
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny-1,1,1)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny-1
                p_ta(i,j)=d_ta(i,j)
        end do
        !ccccccccccccccccccccccc
        ! solve for i > 1
        !ccccccccccccccccccccccc
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        do j=1,nny
        do i=2,nnx

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: distinguish real and imag modified
		!          wavnumbers in x-direction
		!ccccccccccccccccccccccccccccccccccccccccccccc
	
		if(i.le.ngxm/2+1) then
			ii = i
		else
			ii = i - ngxm/2 
		end if	

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: solve using TDMA
		!ccccccccccccccccccccccccccccccccccccccccccccc

                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                !b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(ii))
		! MZW
                b_ta(i,j)=-app_ta(i,j)+(vkz_ta+vkx(ii))/vp_ta(i,j)
                a_ta(i,j)=-aps_ta(i,j)

        end do
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny,2,nnx)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny
        do i=2,nnx
                p_ta(i,j)=d_ta(i,j)
        end do
        end do

	return
end
!--------------------------------------------------------

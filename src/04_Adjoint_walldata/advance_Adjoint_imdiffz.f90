! MZW 06/13/17
! put all the advance parts together
! advance_pre_Adjoint()
! advancexyz_Adjoint()
! advancex_Adjoint_diff()
! advancey_Adjoint_diff()
! advancez_Adjoint_diff()
! advancex_Adjoint_conv()
! advancey_Adjoint_conv()
! advancez_Adjoint_conv()
! advancex_Adjoint_diffex()
! advancey_Adjoint_diffex()
! advancez_Adjoint_diffex()
! advancep_Adjoint()

module advance_module
	use data_export
	integer	:: flag_newdt
	integer	:: flag_x, flag_y, flag_z
	real	:: x_QUICK
end module

subroutine advance_init()
	use advance_module
	integer	:: iflag
	!----------------------------------------
	! (flag_newdt = 1) ==> changed dt
	! (falg_newdt = 0) ==> same old dt
	! Used in advance.f90 to figure out (r1,r2)
	!----------------------------------------
	call archive_isDefined("flag_newdt",iflag)
	if (iflag == 1) then
		call readInt("flag_newdt", flag_newdt)
	else
		flag_newdt = 1
	endif

	flag_x = flag_newdt
	flag_y = flag_newdt
	flag_z = flag_newdt

	!----------------------------------------
	! (x_QUICK = 1.) ==> Fletcher's QUICK q = 1.*(3/8)/3
	! (x_QUICK = x.) ==> Fletcher's QUICK q = x.*(3/8)/3
	!----------------------------------------
	call archive_isDefined("x_QUICK",iflag)
	if (iflag == 1) then
		call readFloat("x_QUICK", x_QUICK)
	else
		x_QUICK = 1.
	endif
	
	return
end

subroutine advance_write()
	use advance_module
	
	if (flag_x.ne.flag_y  .or.  flag_y.ne.flag_z)  &
		stop "flags in advance_xyz are not consistent"
	flag_newdt = flag_x * flag_y * flag_z
	call writeInt("flag_newdt", flag_newdt)
	
	call writeFloat("x_QUICK", x_QUICK)
	return
end

subroutine advance_pre_Adjoint(nsteps)
	use advance_module
	integer nsteps
	! ntime.eq.0 needs to be changed
	! flag_x should be changed to flag_newdt 
	!if(ntime.eq.0  .or.  flag_x.eq.1) then
	!write(*,*) 'ntime,ntime_+nsteps',ntime,ntime_+nsteps
	if(ntime.eq. ntime_+nsteps  .or.  flag_x.eq.1) then
		r1=1.
		r2=0.5
		s1=1.
		s2=0.
		conx = 0.0
		flag_x = 0
	else
		s1=0.5
		s2=0.5
!JL begin
                r1=1.0 + 0.5*dt/dt_old
                r2=      0.5*dt/dt_old
                if( irank==1 .and. r1/=1.5 ) write(*,*) 'r1 and r2 of Adams-Bashforth =',r1, r2
!JL end
	end if
	!write(*,*) 'ntime,ntime_,steps in adjoint',ntime,ntime_,steps
	ust_adjs = r1*ust_adj - r2*ust_adjs
	vst_adjs = r1*vst_adj - r2*vst_adjs
	wst_adjs = r1*wst_adj - r2*wst_adjs
	!write(*,*) 'r1,r2 in adjoint:',r1,r2
end

subroutine source_cost()
        use advance_module
        real :: weight_f, weight_t
        weight_f = visc ** 2
        if(ifric == 1) then
	        do jj=j1_u,j2_u 
	        do ii=i1_u,i2_u 
	        do kk=1,ngzm
                        i = ibmap_1(ii)
                        j = jbmap_1(jj)
	        	k = kk
	                if(j == 1) then
	                        dJdu = (dudy_tmp(kk,ii,1) - dudy_obs(kk,ii,1)) * surxix(i,j)/ (0.5*(ypg(i,2) - ypg(i,1)))
	                elseif(j == ngy-1) then
	                        dJdu = (dudy_tmp(kk,ii,2) - dudy_obs(kk,ii,2)) * surxix(i,j)/(0.5*(ypg(i,ngy-1) - ypg(i,ngy)))
	                else 
	                        dJdu = 0.0
	                end if

	        	u_adj(k,ii,jj) = u_adj(k,ii,jj) + weight_f * dJdu 
	        end do
	        end do
	        end do
                
	        do jj=j1_w,j2_w 
	        do ii=i1_w,i2_w 
	        do kk=1,ngzm
                        i = ibmap_1(ii)
                        j = jbmap_1(jj)
	        	k = kk
	                if(j == 1) then
	                        dJdw = (dwdy_tmp(kk,ii,1) - dwdy_obs(kk,ii,1)) * swrz(i,j)/ (0.5*(ypg(i,2) - ypg(i,1)))
	                elseif(j == ngy-1) then
	                        dJdw = (dwdy_tmp(kk,ii,2) - dwdy_obs(kk,ii,2)) * swrz(i,j)/(0.5*(ypg(i,ngy-1) - ypg(i,ngy)))
	                else 
	                        dJdw = 0.0
	                end if
                        w_adj(k,ii,jj)= w_adj(k,ii,jj) + weight_f * dJdw
                end do
                end do
                end do

        end if

end

!subroutine calculate_dJdu(dJdu,kk,ii,jj)
!	use advance_module
!	real :: dJdu
!	integer :: kk,ii,jj,jout
!
!	i = ibmap_1(ii)
!	j = jbmap_1(jj)
!	if(j == 1) then
!		dJdu = (dudy_tmp(kk,ii,1) - dudy_obs(kk,ii,1)) * surxix(i,j)/ (0.5*(ypg(i,2) - ypg(i,1)))
!	elseif(j == ngy-1) then
!		dJdu = (dudy_tmp(kk,ii,2) - dudy_obs(kk,ii,2)) * surxix(i,j)/(0.5*(ypg(i,ngy-1) - ypg(i,ngy)))
!	else 
!		dJdu = 0.0
!	end if
!
!
!
!end
!
!subroutine calculate_dJdv(dJdv,kk,ii,jj)
!	use advance_module
!	real :: dJdv
!	integer :: kk,ii,jj
!
!	dJdv = 0.0		! not correct for curvilinear !
!
!end
!
!subroutine calculate_dJdw(dJdw,kk,ii,jj)
!	use advance_module
!	real :: dJdw
!	integer :: kk,ii,jj
!
!	i = ibmap_1(ii)
!	j = jbmap_1(jj)
!	if(j == 1) then
!		dJdw = (dwdy_tmp(kk,ii,1) - dwdy_obs(kk,ii,1)) * swrz(i,j)/ (0.5*(ypg(i,2) - ypg(i,1)))
!	elseif(j == ngy-1) then
!		dJdw = (dwdy_tmp(kk,ii,2) - dwdy_obs(kk,ii,2)) * swrz(i,j)/(0.5*(ypg(i,ngy-1) - ypg(i,ngy)))
!	else 
!		dJdw = 0.0
!	end if
!
!end

subroutine advancexyz_Adjoint(cpu_advx,cpu_advy,cpu_advz,iobs)
	use advance_module
	real, intent(out) :: cpu_advx, cpu_advy, cpu_advz
	integer :: iPrj, iobs
        !real :: weight_f        ! for skin friction 

	u_adj=0.0; v_adj=0.0; w_adj=0.0

	call readInt("iPrj",iPrj)   ! iPrj = 1 for fixed mdot, 0 for fixed dpdx
        !weight_f = visc ** 2
	!weight_f = 1

	beforex = realClock()	
	call advancex_Adjoint_diff()
	call advancex_Adjoint_conv()
	call advancex_Adjoint_diffex()
	!sigma_inv2 = 1./sigma**2
	! update u_adj(k,ii,jj)
	do jj=j1_u,j2_u 
	do ii=i1_u,i2_u 
	do kk=1,ngzm
		k = kk
		!call calculate_dJdu(dJdu,kk,ii,jj)
		!u_adj(k,ii,jj)= u_adj(k,ii,jj) + rhsu(k,ii,jj) + weight_f * dJdu - iPrj*GPxi_adj
		u_adj(k,ii,jj)= u_adj(k,ii,jj) + rhsu(k,ii,jj) - iPrj*GPxi_adj
	end do
	end do
	end do
	afterx = realClock()
	cpu_advx = afterx - beforex - t_over

	beforey = realClock()
	call advancey_Adjoint_diff()
	call advancey_Adjoint_conv()
	call advancey_Adjoint_diffex()
	! update v_adj(k,ii,jj)
	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
	do kk=1,ngzm
		k = kk
		!call calculate_dJdv(dJdv,kk,ii,jj)
		!v_adj(k,ii,jj)= v_adj(k,ii,jj) + rhsv(k,ii,jj) + weight_f * dJdv 
		v_adj(k,ii,jj)= v_adj(k,ii,jj) + rhsv(k,ii,jj) 
	end do
	end do
	end do
	aftery = realClock()
	cpu_advy = aftery - beforey - t_over

	beforez = realClock()
	call advancez_Adjoint_diff()
	call advancez_Adjoint_conv()
	call advancez_Adjoint_diffex()
	! update w_adj(k,ii,jj)
	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do kk=1,ngzm
		k = kk
		!call calculate_dJdw(dJdw,kk,ii,jj)
		!w_adj(k,ii,jj)= w_adj(k,ii,jj) + rhsw(k,ii,jj) + weight_f * dJdw 
		w_adj(k,ii,jj)= w_adj(k,ii,jj) + rhsw(k,ii,jj)
	end do
	end do
	end do
	afterz = realClock()
	cpu_advz = afterz - beforez - t_over

	write(*,*) 'irank, u_adj',irank,maxval(abs(u_adj))
	write(*,*) 'irank, v_adj',irank,maxval(abs(v_adj))
	write(*,*) 'irank, w_adj',irank,maxval(abs(w_adj))
        !-----------source term! the derivative of cost function-----------
        if(iobs == 1)     call source_cost()

	! save ust_adj, vst_adj, wst_adj
	ust_adjs = ust_adj
	vst_adjs = vst_adj
	wst_adjs = wst_adj

end

subroutine advancex_Adjoint_diff()
	use advance_module
	
	!ccccccccccccccccccccc
	!   for u_adj equation
	!ccccccccccccccccccccc
	do jj=j1_u,j2_u
	do ii=i1_u,i2_u
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! implicit diffusion term  D^{\dag}_{\xi}
		!-----------------------------------------------------------------------
		! (eq A-164) first term
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1 = suxix (i,j)
		a2 = suxiy (i,j)
		b1 = surxix(i,j)
		b2 = surxiy(i,j)
		c1 = (a1*b1+b1*a1)
		c2 = (a1*b2+b1*a2)
		c4 = (a2*b1+b2*a1)
		c5 = (a2*b2+b2*a2)

		d1 = (suxix(i,j)+suxix(i+1,j))/2.
		d2 = (suxiy(i,j)+suxiy(i+1,j))/2. 
		e1 = (d1*c1+d2*c4)*   visc/vp(i,j)
		e2 = (d1*c2+d2*c5)*   visc/vp(i,j)
		
		cc = suxix(i+1,j)*e1+suxiy(i+1,j)*e2 
		beta=(vp(i,j)+vp(i+1,j))/2.
		diff_adj= cc*ust_adj(k,ii+1,jj)/beta	
		!--------------------------------------------
		!-----------(eq A-164)--second term----------
		!--------------------------------------------	
		d1 = (suxix(i,j)+suxix(i+1,j))/2.
		d2 = (suxiy(i,j)+suxiy(i+1,j))/2.
		e1 = (d1*c1+d2*c4)*   visc/vp(i,j)
		e2 = (d1*c2+d2*c5)*   visc/vp(i,j)

		d1 = (suxix(i,j)+suxix(i-1,j))/2.
		d2 = (suxiy(i,j)+suxiy(i-1,j))/2.
		e1 = e1 + (d1*c1+d2*c4)*   visc/vp(i-1,j)
		e2 = e2 + (d1*c2+d2*c5)*   visc/vp(i-1,j)
	
		bb = -(suxix(i,j)*e1+suxiy(i,j)*e2)
		beta=(vp(i-1,j)+vp(i,j))/2.
		diff_adj= diff_adj + bb*ust_adj(k,ii,jj)/beta	
		!--------------------------------------------
		!-----------(eq A-164)--third term-----------
		!--------------------------------------------
		d1 = (suxix(i,j)+suxix(i-1,j))/2.
		d2 = (suxiy(i,j)+suxiy(i-1,j))/2. 
		e1 = (d1*c1+d2*c4)*   visc/vp(i-1,j)
		e2 = (d1*c2+d2*c5)*   visc/vp(i-1,j)
	
		aa = suxix(i-1,j)*e1+suxiy(i-1,j)*e2
		beta=(vp(i-2,j)+vp(i-1,j))/2.
		diff_adj= diff_adj + aa*ust_adj(k,ii-1,jj)/beta

		!--------------------------------------------
		!-----------(eq A-164)--fourth term----------
		!--------------------------------------------
	        a1 = (svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
	        a2 = (svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
	        b1 = surxix(i,j)
	        b2 = surxiy(i,j)
	        c1 = (a1*b1+b1*a1)
	        c2 = (a1*b2+b1*a2)
	        c4 = (a2*b1+b2*a1)
	        c5 = (a2*b2+b2*a2)

	        d1 = (svetax(i,j+1)+svetax(i-1,j+1))/2.
	        d2 = (svetay(i,j+1)+svetay(i-1,j+1))/2.
	        tta= (vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
	
	        e1 = (d1*c1+d2*c4)*       visc/tta
	        e2 = (d1*c2+d2*c5)*       visc/tta

	        cc = suxix(i,j+1)*e1+suxiy(i,j+1)*e2
		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		diff_adj= diff_adj + cc*ust_adj(k,ii,jj+1)/beta
		!--------------------------------------------
		!-----------(eq A-164)--fifth term-----------
		!--------------------------------------------	        
		d1 = (svetax(i,j+1)+svetax(i-1,j+1))/2.
	        d2 = (svetay(i,j+1)+svetay(i-1,j+1))/2.
	        tta= (vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
	        e1 = (d1*c1+d2*c4)*       visc/tta
	        e2 = (d1*c2+d2*c5)*       visc/tta
	        
		d1 = (svetax(i,j)+svetax(i-1,j))/2.
	        d2 = (svetay(i,j)+svetay(i-1,j))/2.
	        tta= (vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.	
	        e1 = e1+(d1*c1+d2*c4)*       visc/tta
	        e2 = e2+(d1*c2+d2*c5)*       visc/tta

	        bb = -(suxix(i,j)*e1+suxiy(i,j)*e2)
		beta=(vp(i-1,j)+vp(i,j))/2.
		diff_adj= diff_adj + bb*ust_adj(k,ii,jj)/beta	
		!--------------------------------------------
		!-----------(eq A-164)--sixth term-----------
		!--------------------------------------------
	        d1 = (svetax(i,j)+svetax(i-1,j))/2.
	        d2 = (svetay(i,j)+svetay(i-1,j))/2.
	        tta= (vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
	        e1 = (d1*c1+d2*c4)*       visc/tta
	        e2 = (d1*c2+d2*c5)*       visc/tta

	        aa = suxix(i,j-1)*e1+suxiy(i,j-1)*e2
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		diff_adj= diff_adj + aa*ust_adj(k,ii,jj-1)/beta

		!--------------------------------------------
		!----------(eq A-164)--seventh term----------
		!--------------------------------------------
	        a3 = (swz(i,j)+swz(i-1,j))/2.
	        b1 = surxix(i,j)
	        b2 = surxiy(i,j)
	        c7 = (a3*b1)
	        c8 = (a3*b2)
	
	        d3 = (swz(i,j)+swz(i-1,j))/2.
	        tta= (vp(i,j)+vp(i-1,j))/2.
	        e1 = (d3*c7)*       visc/tta
	        e2 = (d3*c8)*       visc/tta

	        cc = suxix(i,j)*e1+suxiy(i,j)*e2
		beta=(vp(i-1,j)+vp(i,j))/2.
		diff_adj= diff_adj + cc*ust_adj(k+1,ii,jj)/beta
		!--------------------------------------------
		!-----------(eq A-164)--eighth term----------
		!--------------------------------------------	
	        d3 = (swz(i,j)+swz(i-1,j))/2.
	        tta= (vp(i,j)+vp(i-1,j))/2.
	        e1 = (d3*c7)*       visc/tta
	        e2 = (d3*c8)*       visc/tta

	        d3 = (swz(i,j)+swz(i-1,j))/2.
	        tta= (vp(i,j)+vp(i-1,j))/2.
	        e1 = e1+(d3*c7)*       visc/tta
	        e2 = e2+(d3*c8)*       visc/tta

	        bb = -(suxix(i,j)*e1+suxiy(i,j)*e2)
		beta=(vp(i-1,j)+vp(i,j))/2.
		diff_adj= diff_adj + bb*ust_adj(k,ii,jj)/beta		
		!--------------------------------------------
		!-----------(eq A-164)--ninth term----------
		!--------------------------------------------
	        d3 = (swz(i,j)+swz(i-1,j))/2.
	        tta= (vp(i,j)+vp(i-1,j))/2.
	        e1 = (d3*c7)*       visc/tta
	        e2 = (d3*c8)*       visc/tta

	        aa = suxix(i,j)*e1+suxiy(i,j)*e2
		beta=(vp(i-1,j)+vp(i,j))/2.
		diff_adj= diff_adj + aa*ust_adj(k-1,ii,jj)/beta

		u_adj(k,ii,jj) = u_adj(k,ii,jj) + dt*diff_adj

	end do
	end do
	end do

end

subroutine advancey_Adjoint_diff()
	use advance_module
	
	!ccccccccccccccccccccc
	!   for v_adj equation
	!ccccccccccccccccccccc
	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! implicit diffusion term  D^{\dag}_{\eta}
		!-----------------------------------------------------------------------
		! (eq A-166) first term
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	        a1 = (suxix  (i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        	a2 = (suxiy  (i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        	b1 = svretax(i,j)
        	b2 = svretay(i,j)
        	c1 = (a1*b1+b1*a1)
        	c2 = (a1*b2+b1*a2)
        	c4 = (a2*b1+b2*a1)
        	c5 = (a2*b2+b2*a2)
	
        	d1 = (suxix(i+1,j)+suxix(i+1,j-1))/2.
        	d2 = (suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        	tta= (vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        	e1 = (d1*c1+d2*c4)*       visc/tta
        	e2 = (d1*c2+d2*c5)*       visc/tta 

		cc = svetax(i+1,j)*e1+svetay(i+1,j)*e2
		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		diff_adj= cc*vst_adj(k,ii+1,jj)/beta
		!--------------------------------------------
		!-----------(eq A-166)--second term-----------
		!--------------------------------------------
	        d1 = (suxix(i+1,j)+suxix(i+1,j-1))/2.
        	d2 = (suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        	tta= (vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        	e1 = (d1*c1+d2*c4)*       visc/tta 
        	e2 = (d1*c2+d2*c5)*       visc/tta 

        	d1 = (suxix(i,j)+suxix(i,j-1))/2.
        	d2 = (suxiy(i,j)+suxiy(i,j-1))/2.
        	tta= (vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        	e1 = e1+(d1*c1+d2*c4)*       visc/tta 
        	e2 = e2+(d1*c2+d2*c5)*       visc/tta 

		bb = -(svetax(i,j)*e1+svetay(i,j)*e2)
		beta=(vp(i,j-1)+vp(i,j))/2.
		diff_adj= diff_adj + bb*vst_adj(k,ii,jj)/beta
		!--------------------------------------------
		!-----------(eq A-166)--third term-----------
		!--------------------------------------------
    	        d1 = (suxix(i,j)+suxix(i,j-1))/2.
        	d2 = (suxiy(i,j)+suxiy(i,j-1))/2.
        	tta= (vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        	e1 = (d1*c1+d2*c4)*       visc/tta 
        	e2 = (d1*c2+d2*c5)*       visc/tta 

        	aa = svetax(i-1,j)*e1+svetay(i-1,j)*e2
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		diff_adj= diff_adj + aa*vst_adj(k,ii-1,jj)/beta

		!--------------------------------------------
		!-----------(eq A-166)--fourth term----------
		!--------------------------------------------
     		a1 = svetax (i,j)
        	a2 = svetay (i,j)
        	b1 = svretax(i,j)
	        b2 = svretay(i,j)
	        c1 = (a1*b1+b1*a1)
	        c2 = (a1*b2+b1*a2)
	        c4 = (a2*b1+b2*a1)
	        c5 = (a2*b2+b2*a2)

	        d1 = (svetax(i,j)+svetax(i,j+1))/2.
	        d2 = (svetay(i,j)+svetay(i,j+1))/2.
	        tta= vp(i,j)
	        e1 = (d1*c1+d2*c4)*       visc/tta 
	        e2 = (d1*c2+d2*c5)*       visc/tta 

	        cc = svetax(i,j+1)*e1+svetay(i,j+1)*e2
		beta=(vp(i,j)+vp(i,j+1))/2.
		diff_adj= diff_adj + cc*vst_adj(k,ii,jj+1)/beta
		!--------------------------------------------
		!-----------(eq A-166)--fifth term----------
		!--------------------------------------------
        	d1 = (svetax(i,j)+svetax(i,j+1))/2.
       	 	d2 = (svetay(i,j)+svetay(i,j+1))/2.
        	tta= vp(i,j)
        	e1 = (d1*c1+d2*c4)*       visc/tta
        	e2 = (d1*c2+d2*c5)*       visc/tta 

        	d1 = (svetax(i,j)+svetax(i,j-1))/2.
        	d2 = (svetay(i,j)+svetay(i,j-1))/2.
        	tta= vp(i,j-1)
        	e1 = e1+(d1*c1+d2*c4)*       visc/tta
        	e2 = e2+(d1*c2+d2*c5)*       visc/tta 

        	bb = -(svetax(i,j)*e1+svetay(i,j)*e2)
		beta=(vp(i,j-1)+vp(i,j))/2.
		diff_adj= diff_adj + bb*vst_adj(k,ii,jj)/beta	
		!--------------------------------------------
		!-----------(eq A-166)--sixth term----------
		!--------------------------------------------
   	     	d1 = (svetax(i,j)+svetax(i,j-1))/2.
        	d2 = (svetay(i,j)+svetay(i,j-1))/2.
        	tta= vp(i,j-1)
        	e1 = (d1*c1+d2*c4)*       visc/tta 
        	e2 = (d1*c2+d2*c5)*       visc/tta 

        	aa = svetax(i,j-1)*e1+svetay(i,j-1)*e2
		beta=(vp(i,j-2)+vp(i,j-1))/2.
		diff_adj= diff_adj + aa*vst_adj(k,ii,jj-1)/beta

		!--------------------------------------------
		!-----------(eq A-166)--seventh term----------
		!--------------------------------------------
        	a3 = (swz(i,j)+swz(i,j-1))/2.
        	b1 = svretax(i,j)
        	b2 = svretay(i,j)
        	c7 = (a3*b1)
        	c8 = (a3*b2)

        	d3 = (swz(i,j)+swz(i,j-1))/2.
        	tta= (vp(i,j)+vp(i,j-1))/2.
        	e1 = (d3*c7)*       visc/tta
        	e2 = (d3*c8)*       visc/tta

        	cc = svetax(i,j)*e1+svetay(i,j)*e2
		beta=(vp(i,j-1)+vp(i,j))/2.
		diff_adj= diff_adj + cc*vst_adj(k+1,ii,jj)/beta
		!--------------------------------------------
		!-----------(eq A-166)--eighth term----------
		!--------------------------------------------
        	d3 = (swz(i,j)+swz(i,j-1))/2.
        	tta= (vp(i,j)+vp(i,j-1))/2.
        	e1 = (d3*c7)*       visc/tta
	        e2 = (d3*c8)*       visc/tta
	
	        d3 = (swz(i,j)+swz(i,j-1))/2.
	        tta= (vp(i,j)+vp(i,j-1))/2.
	        e1 = e1+(d3*c7)*       visc/tta
	        e2 = e2+(d3*c8)*       visc/tta

	        bb = -(svetax(i,j)*e1+svetay(i,j)*e2)
		beta=(vp(i,j-1)+vp(i,j))/2.
		diff_adj= diff_adj + bb*vst_adj(k,ii,jj)/beta	
		!--------------------------------------------
		!-----------(eq A-166)--ninth term----------
		!--------------------------------------------
	        d3 = (swz(i,j)+swz(i,j-1))/2.
        	tta= (vp(i,j)+vp(i,j-1))/2.
        	e1 = (d3*c7)*       visc/tta
        	e2 = (d3*c8)*       visc/tta

        	aa = svetax(i,j)*e1+svetay(i,j)*e2
		beta=(vp(i,j-1)+vp(i,j))/2.
		diff_adj= diff_adj + aa*vst_adj(k-1,ii,jj)/beta

		v_adj(k,ii,jj) = v_adj(k,ii,jj) + dt*diff_adj

	end do
	end do
	end do

end

subroutine advancez_Adjoint_diff()
	use advance_module
	
	!ccccccccccccccccccccc
	!   for w_adj equation
	!ccccccccccccccccccccc
	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! implicit diffusion term  D^{\dag}_{\xi}
		!-----------------------------------------------------------------------
		! (eq A-168) first term
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	        a1 = (suxix(i,j)+suxix(i+1,j))/2.
        	a2 = (suxiy(i,j)+suxiy(i+1,j))/2.
        	b3 = swrz (i,j)
        	c3 = (a1*b3)
        	c6 = (a2*b3)

	        d1 = suxix(i+1,j)
        	d2 = suxiy(i+1,j)
        	tta= (vp(i,j)+vp(i+1,j))/2. 
        	e3 = (d1*c3+d2*c6)*       visc/tta 

        	cc = swz(i+1,j)*e3
		diff_adj= cc*wst_adj(k,ii+1,jj)/vp(i+1,j)
		!--------------------------------------------
		!-----------(eq A-168)--second term-----------
		!--------------------------------------------
	        d1 = suxix(i+1,j)
        	d2 = suxiy(i+1,j)
        	tta= (vp(i,j)+vp(i+1,j))/2. 
        	e3 = (d1*c3+d2*c6)*       visc/tta 

        	d1 = suxix(i,j)
        	d2 = suxiy(i,j)
        	tta= (vp(i,j)+vp(i-1,j))/2.
        	e3 = e3+(d1*c3+d2*c6)*       visc/tta 

        	bb = -(swz(i,j)*e3)
		diff_adj= diff_adj + bb*wst_adj(k,ii,jj)/vp(i,j)
		!--------------------------------------------
		!-----------(eq A-168)--third term-----------
		!--------------------------------------------
	        d1 = suxix(i,j)
        	d2 = suxiy(i,j)
        	tta= (vp(i,j)+vp(i-1,j))/2.
        	e3 = (d1*c3+d2*c6)*       visc/tta 

        	aa = swz(i-1,j)*e3
		diff_adj= diff_adj + aa*wst_adj(k,ii-1,jj)/vp(i-1,j)

		!--------------------------------------------
		!-----------(eq A-168)--fourth term-----------
		!--------------------------------------------
	        a1 = (svetax(i,j)+svetax(i,j+1))/2.
        	a2 = (svetay(i,j)+svetay(i,j+1))/2.
        	b3 = swrz  (i,j)
        	c3 = (a1*b3)
        	c6 = (a2*b3)
	
        	d1 = svetax(i,j+1)
        	d2 = svetay(i,j+1)
        	tta= (vp(i,j)+vp(i,j+1))/2. 
        	e3 = (d1*c3+d2*c6)*       visc/tta 

        	cc = swz(i,j+1)*e3 
		diff_adj= diff_adj + cc*wst_adj(k,ii,jj+1)/vp(i,j+1)
		!--------------------------------------------
		!-------------(eq A-168)--fifth term-----------
		!--------------------------------------------
	        d1 = svetax(i,j+1)
        	d2 = svetay(i,j+1)
        	tta= (vp(i,j)+vp(i,j+1))/2. 
        	e3 = (d1*c3+d2*c6)*       visc/tta 
	
        	d1 = svetax(i,j)
        	d2 = svetay(i,j)
        	tta= (vp(i,j)+vp(i,j-1))/2.
        	e3 = e3+(d1*c3+d2*c6)*       visc/tta 

        	bb = -(swz(i,j)*e3)
		diff_adj= diff_adj + bb*wst_adj(k,ii,jj)/vp(i,j)
		!--------------------------------------------
		!-------------(eq A-168)--sixth term-----------
		!--------------------------------------------
	        d1 = svetax(i,j)
        	d2 = svetay(i,j)
        	tta= (vp(i,j)+vp(i,j-1))/2.
        	e3 = (d1*c3+d2*c6)*       visc/tta 

        	aa = swz(i,j-1)*e3
		diff_adj= diff_adj + aa*wst_adj(k,ii,jj-1)/vp(i,j-1)

		!--------------------------------------------
		!-------------(eq A-168)--seventh term-----------
		!--------------------------------------------
	        a3 = swz(i,j)
        	b3 = swrz(i,j)
        	c9 = (a3*b3+b3*a3)

	        d3 = swz(i,j)
	        tta= vp(i,j)
	        e3 = (d3*c9)*       visc/tta

	        cc=swz(i,j)*e3
		diff_adj= diff_adj + cc*wst_adj(k+1,ii,jj)/vp(i,j)		
		!--------------------------------------------
		!-------------(eq A-168)--eighth term-----------
		!--------------------------------------------
	        d3 = swz(i,j)
	        tta= vp(i,j)
	        e3 = (d3*c9)*       visc/tta
	
	        d3 = swz(i,j)
	        tta= vp(i,j)
	        e3 = e3+(d3*c9)*       visc/tta

	        bb = -swz(i,j)*e3
		diff_adj= diff_adj + bb*wst_adj(k,ii,jj)/vp(i,j)
		!--------------------------------------------
		!-------------(eq A-168)--ninth term-----------
		!--------------------------------------------
	        d3 = swz(i,j)
	        tta= vp(i,j)
	        e3 = (d3*c9)*       visc/tta

	        aa = swz(i,j)*e3
		diff_adj= diff_adj + aa*wst_adj(k-1,ii,jj)/vp(i,j)

		w_adj(k,ii,jj) = w_adj(k,ii,jj) + dt*diff_adj

	end do
	end do
	end do



end

subroutine advancex_Adjoint_conv()

	use advance_module

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)
       		vst_adjs(:, :, 1) =  0.0			
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)		

       	      	!Extend
		vst_adjs(:, :, 0) = vst_adjs(:, :, 2)		! this halo is not included in following codes
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)
              	vst_adjs(:, :, nlyb  ) =  0.0			
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)	
		
		!Extend
		vst_adjs(:, :, nlyb+1) = vst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
     		
	end if 

	!ccccccccccccccccccccc
	!   for u_adj equation
	!ccccccccccccccccccccc
	do jj=j1_u,j2_u 
	do ii=i1_u,i2_u 
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\xi,1}
		!-----------------------------------------------------------------------
		! first term h1 = a \bu b  in (A-176)
		! a = u^{\xi}_i S_i - u^{\xi}_{i+1} S_{i+1}  adjoint
		! b = u^{\xi}_{i+1/2,j} S_{i+1/2,j} + 1/2 * u^{\eta}_{i,j+1/2} S_{i,j+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h1 = 0.0
		!write(*,*) '---------adjoint convection term-----------'
		beta1=(vp(i-1,j)+vp(i,j))/2.
		beta2=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k,ii,jj)/beta1*suxix(i,j) - ust_adjs(k,ii+1,jj)/beta2*suxix(i+1,j)
		ay = ust_adjs(k,ii,jj)/beta1*suxiy(i,j) - ust_adjs(k,ii+1,jj)/beta2*suxiy(i+1,j)
		bx = (u(k,ii,jj)+u(k,ii+1,jj))/2. * (surxix(i,j)+surxix(i+1,j))/2.
		bx = bx + 0.5 * (v(k,ii,jj)+v(k,ii,jj+1))/2. * (svretax(i,j)+svretax(i,j+1))/2.
		by = (u(k,ii,jj)+u(k,ii+1,jj))/2. * (surxiy(i,j)+surxiy(i+1,j))/2.
		by = by + 0.5 * (v(k,ii,jj)+v(k,ii,jj+1))/2. * (svretay(i,j)+svretay(i,j+1))/2.

		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=' , h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h1 = h1 + a \bu b  in (A-176)
		! a = u^{\xi}_{i-1} S_{i-1} - u^{\xi}_{i} S_{i}   adjoint
		! b = u^{\xi}_{i-1/2,j} S_{i-1/2,j} + 1/2 * u^{\eta}_{i-1,j+1/2} S_{i-1,j+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-2,j)+vp(i-1,j))/2.
		beta2=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k,ii-1,jj)/beta1*suxix(i-1,j) - ust_adjs(k,ii,jj)/beta2*suxix(i,j)
		ay = ust_adjs(k,ii-1,jj)/beta1*suxiy(i-1,j) - ust_adjs(k,ii,jj)/beta2*suxiy(i,j)
		bx = (u(k,ii-1,jj)+u(k,ii,jj))/2. * (surxix(i-1,j)+surxix(i,j))/2.
		bx = bx + 0.5 * (v(k,ii-1,jj)+v(k,ii-1,jj+1))/2. * (svretax(i-1,j)+svretax(i-1,j+1))/2.
		by = (u(k,ii-1,jj)+u(k,ii,jj))/2. * (surxiy(i-1,j)+surxiy(i,j))/2.
		by = by + 0.5 * (v(k,ii-1,jj)+v(k,ii-1,jj+1))/2. * (svretay(i-1,j)+svretay(i-1,j+1))/2.

		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=', h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! third term h1 = h1 + a \bu b  in (A-176)
		! a = u^{\xi}_j S_j - u^{\xi}_{j+1} S_{j+1}  adjoint
		! b = 1/2*u^{\eta}_{i-1/2,j+1,k} S_{i,j+1/2,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j)+vp(i,j))/2.
		beta2=(vp(i-1,j+1)+vp(i,j+1))/2.
		ax = ust_adjs(k,ii,jj)/beta1*suxix(i,j) - ust_adjs(k,ii,jj+1)/beta2*suxix(i,j+1)
		ay = ust_adjs(k,ii,jj)/beta1*suxiy(i,j) - ust_adjs(k,ii,jj+1)/beta2*suxiy(i,j+1)
		bx = 0.5 * (v(k,ii-1,jj+1)+v(k,ii,jj+1) )/2. * (surxix(i,j)+surxix(i,j+1))/2.
		by = 0.5 * (v(k,ii-1,jj+1)+v(k,ii,jj+1) )/2. * (surxiy(i,j)+surxiy(i,j+1))/2.


		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=', h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fourth term h1 = h1 + a \bu b  in (A-176)
		! a = u^{\xi}_{j-1} S_{j-1} - u^{\xi}_{j} S_{j}  adjoint
		! b = 1/2*u^{\eta}_{i-1/2,j,k} S_{i,j-1/2,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j-1)+vp(i,j-1))/2.
		beta2=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k,ii,jj-1)/beta1*suxix(i,j-1) - ust_adjs(k,ii,jj)/beta2*suxix(i,j)
		ay = ust_adjs(k,ii,jj-1)/beta1*suxiy(i,j-1) - ust_adjs(k,ii,jj)/beta2*suxiy(i,j)
		bx = 0.5 * (v(k,ii-1,jj)+v(k,ii,jj) )/2. * (surxix(i,j-1)+surxix(i,j))/2.
		by = 0.5 * (v(k,ii-1,jj)+v(k,ii,jj) )/2. * (surxiy(i,j-1)+surxiy(i,j))/2.

		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=', h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fifth term h1 = h1 + a \bu b  in (A-176)
		! a = u^{\xi}_k S_k - u^{\xi}_{k+1} S_{k+1}  adjoint
		! b = 1/2*u^{z}_{i-1/2,j,k+1} S_{i,j,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k,ii,jj)/beta*suxix(i,j) - ust_adjs(k+1,ii,jj)/beta*suxix(i,j)
		ay = ust_adjs(k,ii,jj)/beta*suxiy(i,j) - ust_adjs(k+1,ii,jj)/beta*suxiy(i,j)
		bx = 0.5 * (w(k+1,ii-1,jj)+w(k+1,ii,jj) )/2. * (surxix(i,j)+surxix(i,j))/2.
		by = 0.5 * (w(k+1,ii-1,jj)+w(k+1,ii,jj) )/2. * (surxiy(i,j)+surxiy(i,j))/2.

		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=', h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sixth term h1 = h1 + a \bu b  in (A-176)
		! a = u^{\xi}_k S_k - u^{\xi}_{k+1} S_{k+1}  adjoint
		! b = 1/2*u^{z}_{i-1/2,j,k-1} S_{i,j,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k-1,ii,jj)/beta*suxix(i,j) - ust_adjs(k,ii,jj)/beta*suxix(i,j)
		ay = ust_adjs(k-1,ii,jj)/beta*suxiy(i,j) - ust_adjs(k,ii,jj)/beta*suxiy(i,j)
		bx = 0.5 * (w(k,ii-1,jj)+w(k,ii,jj) )/2. * (surxix(i,j)+surxix(i,j))/2.
		by = 0.5 * (w(k,ii-1,jj)+w(k,ii,jj) )/2. * (surxiy(i,j)+surxiy(i,j))/2.

		h1 = h1 + ax*bx + ay*by
		!write(*,*) 'h1=', h1
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\xi,2}
		!-----------------------------------------------------------------------
		! first term h1 = a \bu b  in (A-177)
		! a = u^{\eta}_{i-1} S_{i-1} - u^{\eta}_{i} S_{i}  adjoint
		! b = u^{\xi}_{i,j-1/2} S_{i,j-1/2} + 1/2 * u^{\eta}_{i-1/2,j} S_{i-1/2,j}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h2 = 0.0
		
		beta1=(vp(i-1,j-1)+vp(i-1,j))/2.
		beta2=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k,ii-1,jj)/beta1*svetax(i-1,j) - vst_adjs(k,ii,jj)/beta2*svetax(i,j)
		ay = vst_adjs(k,ii-1,jj)/beta1*svetay(i-1,j) - vst_adjs(k,ii,jj)/beta2*svetay(i,j)
		bx = (u(k,ii,jj-1)+u(k,ii,jj))/2. * (surxix(i,j-1)+surxix(i,j))/2.
		bx = bx + 0.5 * (v(k,ii-1,jj)+v(k,ii,jj))/2. * (svretax(i-1,j)+svretax(i,j))/2.
		by = (u(k,ii,jj-1)+u(k,ii,jj))/2. * (surxiy(i,j-1)+surxiy(i,j))/2.
		by = by + 0.5 * (v(k,ii-1,jj)+v(k,ii,jj))/2. * (svretay(i-1,j)+svretay(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h2 = a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,j+1} S_{i-1,j+1} - u^{\eta}_{i,j+1} S_{i,j+1}  adjoint
		! b = u^{\xi}_{i,j+1/2} S_{i,j+1/2} + 1/2 * u^{\eta}_{i-1/2,j+1} S_{i-1/2,j+1}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j)+vp(i-1,j+1))/2.
		beta2=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k,ii-1,jj+1)/beta1*svetax(i-1,j+1) - vst_adjs(k,ii,jj+1)/beta2*svetax(i,j+1)
		ay = vst_adjs(k,ii-1,jj+1)/beta1*svetay(i-1,j+1) - vst_adjs(k,ii,jj+1)/beta2*svetay(i,j+1)
		bx = (u(k,ii,jj)+u(k,ii,jj+1))/2. * (surxix(i,j)+surxix(i,j+1))/2.
		bx = bx + 0.5 * (v(k,ii-1,jj+1)+v(k,ii,jj+1))/2. * (svretax(i-1,j+1)+svretax(i,j+1))/2.
		by = (u(k,ii,jj)+u(k,ii,jj+1))/2. * (surxiy(i,j)+surxiy(i,j+1))/2.
		by = by + 0.5 * (v(k,ii-1,jj+1)+v(k,ii,jj+1))/2. * (svretay(i-1,j+1)+svretay(i,j+1))/2.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! third term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_j S_j - u^{\eta}_{j+1} S_{j+1}  adjoint
		! b = 1/2*u^{\eta}_{i,j+1/2,k} S_{\xi, i+1/2,j,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i,j-1)+vp(i,j))/2.
		beta2=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k,ii,jj)/beta1*svetax(i,j) - vst_adjs(k,ii,jj+1)/beta2*svetax(i,j+1)
		ay = vst_adjs(k,ii,jj)/beta1*svetay(i,j) - vst_adjs(k,ii,jj+1)/beta2*svetay(i,j+1)
		bx = 0.5 * (v(k,ii,jj)+v(k,ii,jj+1) )/2. * (surxix(i,j)+surxix(i+1,j))/2.
		by = 0.5 * (v(k,ii,jj)+v(k,ii,jj+1) )/2. * (surxiy(i,j)+surxiy(i+1,j))/2.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fourth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,j} S_{i-1,j} - u^{\eta}_{i-1,j} S_{i-1,j}  adjoint
		! b = 1/2*u^{\eta}_{i-1,j+1/2,k} S_{\xi, i-1/2,j,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j-1)+vp(i-1,j))/2.
		beta2=(vp(i-1,j)+vp(i-1,j+1))/2.
		ax = vst_adjs(k,ii-1,jj)/beta1*svetax(i-1,j) - vst_adjs(k,ii-1,jj+1)/beta2*svetax(i-1,j+1)
		ay = vst_adjs(k,ii-1,jj)/beta1*svetay(i-1,j) - vst_adjs(k,ii-1,jj+1)/beta2*svetay(i-1,j+1)
		bx = 0.5 * (v(k,ii-1,jj)+v(k,ii-1,jj+1) )/2. * (surxix(i-1,j)+surxix(i,j))/2.
		by = 0.5 * (v(k,ii-1,jj)+v(k,ii-1,jj+1) )/2. * (surxiy(i-1,j)+surxiy(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fifth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_k S_k - u^{\eta}_{k+1} S_{k+1}  adjoint
		! b = 1/8*u^{\eta}_{i,j-1/2,k+1} S_{\xi, i+1/2,j-1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k,ii,jj)/beta*svetax(i,j) - vst_adjs(k+1,ii,jj)/beta*svetax(i,j)
		ay = vst_adjs(k,ii,jj)/beta*svetay(i,j) - vst_adjs(k+1,ii,jj)/beta*svetay(i,j)
		bx = 0.125 * (w(k+1,ii,jj-1)+w(k+1,ii,jj) )/2. * (surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
		by = 0.125 * (w(k+1,ii,jj-1)+w(k+1,ii,jj) )/2. * (surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sixth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{k-1} S_{k-1} - u^{\eta}_{k} S_{k}  adjoint
		! b = 1/8*u^{\eta}_{i,j-1/2,k} S_{\xi, i+1/2,j-1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k-1,ii,jj)/beta*svetax(i,j) - vst_adjs(k,ii,jj)/beta*svetax(i,j)
		ay = vst_adjs(k-1,ii,jj)/beta*svetay(i,j) - vst_adjs(k,ii,jj)/beta*svetay(i,j)
		bx = 0.125 * (w(k,ii,jj-1)+w(k,ii,jj) )/2. * (surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
		by = 0.125 * (w(k,ii,jj-1)+w(k,ii,jj) )/2. * (surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! seventh term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{j+1,k} S_{j+1,k} - u^{\eta}_{j+1,k+1} S_{j+1,k+1}  adjoint
		! b = 1/8*u^{\eta}_{i,j+1/2,k+1} S_{\xi, i+1/2,j+1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k,ii,jj+1)/beta*svetax(i,j+1) - vst_adjs(k+1,ii,jj+1)/beta*svetax(i,j+1)
		ay = vst_adjs(k,ii,jj+1)/beta*svetay(i,j+1) - vst_adjs(k+1,ii,jj+1)/beta*svetay(i,j+1)
		bx = 0.125 * (w(k+1,ii,jj)+w(k+1,ii,jj+1) )/2. * (surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		by = 0.125 * (w(k+1,ii,jj)+w(k+1,ii,jj+1) )/2. * (surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! eighth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{j+1,k-1} S_{j+1,k-1} - u^{\eta}_{j+1,k} S_{j+1,k}  adjoint
		! b = 1/8*u^{\eta}_{i,j+1/2,k} S_{\xi, i+1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k-1,ii,jj+1)/beta*svetax(i,j+1) - vst_adjs(k,ii,jj+1)/beta*svetax(i,j+1)
		ay = vst_adjs(k-1,ii,jj+1)/beta*svetay(i,j+1) - vst_adjs(k,ii,jj+1)/beta*svetay(i,j+1)
		bx = 0.125 * (w(k,ii,jj)+w(k,ii,jj+1) )/2. * (surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		by = 0.125 * (w(k,ii,jj)+w(k,ii,jj+1) )/2. * (surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! ninth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,k} S_{i-1,k} - u^{\eta}_{i-1,k+1} S_{i-1,k+1}  adjoint
		! b = 1/8*u^{\eta}_{i-1,j-1/2,k+1} S_{\xi, i-1/2,j-1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		ax = vst_adjs(k,ii-1,jj)/beta*svetax(i-1,j) - vst_adjs(k+1,ii-1,jj)/beta*svetax(i-1,j)
		ay = vst_adjs(k,ii-1,jj)/beta*svetay(i-1,j) - vst_adjs(k+1,ii-1,jj)/beta*svetay(i-1,j)
		bx = 0.125 * (w(k+1,ii-1,jj-1)+w(k+1,ii-1,jj) )/2. * (surxix(i-1,j)+surxix(i,j)+surxix(i-1,j-1)+surxix(i,j-1))/4.
		by = 0.125 * (w(k+1,ii-1,jj-1)+w(k+1,ii-1,jj) )/2. * (surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j-1)+surxiy(i,j-1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! tenth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,k-1} S_{i-1,k-1} - u^{\eta}_{i-1,k} S_{i-1,k}  adjoint
		! b = 1/8*u^{\eta}_{i-1,j-1/2,k} S_{\xi, i-1/2,j-1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		ax = vst_adjs(k-1,ii-1,jj)/beta*svetax(i-1,j) - vst_adjs(k,ii-1,jj)/beta*svetax(i-1,j)
		ay = vst_adjs(k-1,ii-1,jj)/beta*svetay(i-1,j) - vst_adjs(k,ii-1,jj)/beta*svetay(i-1,j)
		bx = 0.125 * (w(k,ii-1,jj-1)+w(k,ii-1,jj) )/2. * (surxix(i-1,j)+surxix(i,j)+surxix(i-1,j-1)+surxix(i,j-1))/4.
		by = 0.125 * (w(k,ii-1,jj-1)+w(k,ii-1,jj) )/2. * (surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j-1)+surxiy(i,j-1))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! eleventh term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,j+1,k} S_{i-1,j+1,k} - u^{\eta}_{i-1,j+1,k+1} S_{i-1,j+1,k+1}  adjoint
		! b = 1/8*u^{\eta}_{i-1,j+1/2,k+1} S_{\xi, i-1/2,j+1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		ax = vst_adjs(k,ii-1,jj+1)/beta*svetax(i-1,j+1) - vst_adjs(k+1,ii-1,jj+1)/beta*svetax(i-1,j+1)
		ay = vst_adjs(k,ii-1,jj+1)/beta*svetay(i-1,j+1) - vst_adjs(k+1,ii-1,jj+1)/beta*svetay(i-1,j+1)
		bx = 0.125 * (w(k+1,ii-1,jj)+w(k+1,ii-1,jj+1) )/2. * (surxix(i-1,j+1)+surxix(i,j+1)+surxix(i-1,j)+surxix(i,j))/4.
		by = 0.125 * (w(k+1,ii-1,jj)+w(k+1,ii-1,jj+1) )/2. * (surxiy(i-1,j+1)+surxiy(i,j+1)+surxiy(i-1,j)+surxiy(i,j))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! twelfth term h2 = h2 + a \bu b  in (A-177)
		! a = u^{\eta}_{i-1,j+1,k-1} S_{i-1,j+1,k-1} - u^{\eta}_{i-1,j+1,k} S_{i-1,j+1,k}  adjoint
		! b = 1/8*u^{\eta}_{i-1,j+1/2,k} S_{\xi, i-1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		ax = vst_adjs(k-1,ii-1,jj+1)/beta*svetax(i-1,j+1) - vst_adjs(k,ii-1,jj+1)/beta*svetax(i-1,j+1)
		ay = vst_adjs(k-1,ii-1,jj+1)/beta*svetay(i-1,j+1) - vst_adjs(k,ii-1,jj+1)/beta*svetay(i-1,j+1)
		bx = 0.125 * (w(k,ii-1,jj)+w(k,ii-1,jj+1) )/2. * (surxix(i-1,j+1)+surxix(i,j+1)+surxix(i-1,j)+surxix(i,j))/4.
		by = 0.125 * (w(k,ii-1,jj)+w(k,ii-1,jj+1) )/2. * (surxiy(i-1,j+1)+surxiy(i,j+1)+surxiy(i-1,j)+surxiy(i,j))/4.

		h2 = h2 + ax*bx + ay*by
		!write(*,*) 'h2=', h2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\xi,3}
		!-----------------------------------------------------------------------
		! first term h3 = a \bu b  in (A-178)
		! a = u^{z}_{i-1} S_{i-1} - u^{z}_{i} S_{i}  adjoint
		! b = 1/2 * u^{z}_{i-1/2,j} S_{i-1/2,j}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h3 = 0.0
		beta1=vp(i-1,j)
		beta2=vp(i,j)
		az = wst_adjs(k,ii-1,jj)/beta1*swz(i-1,j) - wst_adjs(k,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (w(k,ii-1,jj)+w(k,ii,jj))/2. * (swrz(i-1,j)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!write(*,*) 'h3=', h3
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h3 = a \bu b  in (A-178)
		! a = u^{z}_{i-1,k+1} S_{i-1,k+1} - u^{z}_{i,k+1} S_{i,k+1}  adjoint
		! b = 1/2 * u^{z}_{i-1/2,j,k+1} S_{i-1/2,j}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i-1,j)
		beta2=vp(i,j)
		az = wst_adjs(k+1,ii-1,jj)/beta1*swz(i-1,j) - wst_adjs(k+1,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (w(k+1,ii-1,jj)+w(k+1,ii,jj))/2. * (swrz(i-1,j)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!write(*,*) 'h3=', h3
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sum of all terms
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		con_adj = -h1 - h2 - h3

		u_adj(k,ii,jj) = u_adj(k,ii,jj) + dt*con_adj
	end do
	end do
	end do

end

subroutine advancey_Adjoint_conv()

	use advance_module

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)		! this halo is not included in following codes
       		vst_adjs(:, :, 1) =  0.0			
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)		! this halo is not included in following codes

       	      	!Extend
		vst_adjs(:, :, 0) = vst_adjs(:, :, 2)		! this halo is not included in following codes
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)	! this halo is not included in following codes
              	vst_adjs(:, :, nlyb  ) =  0.0			
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
		
		!Extend
		vst_adjs(:, :, nlyb+1) = vst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
     		
	end if 

	!ccccccccccccccccccccc
	!   for v_adj equation
	!ccccccccccccccccccccc
	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\eta,1}
		!-----------------------------------------------------------------------
		! first term h1 = a \bu b  in (A-180)
		! a = u^{\xi}_{j-1} S_{j-1} - u^{\xi}_{j} S_{j}  adjoint
		! b = 1/2 * u^{\xi}_{i,j-1/2} S_{i,j-1/2} + u^{\eta}_{i-1/2,j} S_{i-1/2,j}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h1 = 0.0

		beta1=(vp(i-1,j-1)+vp(i,j-1))/2.
		beta2=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k,ii,jj-1)/beta1*suxix(i,j-1) - ust_adjs(k,ii,jj)/beta2*suxix(i,j)
		ay = ust_adjs(k,ii,jj-1)/beta1*suxiy(i,j-1) - ust_adjs(k,ii,jj)/beta2*suxiy(i,j) 
		bx = 0.5 * (u(k,ii,jj-1)+u(k,ii,jj))/2. * (surxix(i,j-1)+surxix(i,j))/2.
		bx = bx + (v(k,ii-1,jj)+v(k,ii,jj))/2. * (svretax(i-1,j)+svretax(i,j))/2.
		by = 0.5 * (u(k,ii,jj-1)+u(k,ii,jj))/2. * (surxiy(i,j-1)+surxiy(i,j))/2.
		by = by + (v(k,ii-1,jj)+v(k,ii,jj))/2. * (svretay(i-1,j)+svretay(i,j))/2.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i+1,j-1} S_{i+1,j-1} - u^{\xi}_{i+1,j} S_{i+1,j}  adjoint
		! b = 1/2 * u^{\xi}_{i+1,j-1/2} S_{i+1,j-1/2} + u^{\eta}_{i+1/2,j} S_{i+1/2,j}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i,j-1)+vp(i+1,j-1))/2.
		beta2=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k,ii+1,jj-1)/beta1*suxix(i+1,j-1) - ust_adjs(k,ii+1,jj)/beta2*suxix(i+1,j)
		ay = ust_adjs(k,ii+1,jj-1)/beta1*suxiy(i+1,j-1) - ust_adjs(k,ii+1,jj)/beta2*suxiy(i+1,j) 
		bx = 0.5 * (u(k,ii+1,jj-1)+u(k,ii+1,jj))/2. * (surxix(i+1,j-1)+surxix(i+1,j))/2.
		bx = bx + (v(k,ii,jj)+v(k,ii+1,jj))/2. * (svretax(i,j)+svretax(i+1,j))/2.
		by = 0.5 * (u(k,ii+1,jj-1)+u(k,ii+1,jj))/2. * (surxiy(i+1,j-1)+surxiy(i+1,j))/2.
		by = by + (v(k,ii,jj)+v(k,ii+1,jj))/2. * (svretay(i,j)+svretay(i+1,j))/2.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! third term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j} S_{i,j} - u^{\xi}_{i+1,j} S_{i+1,j}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j} S_{i,j+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j)+vp(i,j))/2.
		beta2=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k,ii,jj)/beta1*suxix(i,j) - ust_adjs(k,ii+1,jj)/beta2*suxix(i+1,j)
		ay = ust_adjs(k,ii,jj)/beta1*suxiy(i,j) - ust_adjs(k,ii+1,jj)/beta2*suxiy(i+1,j)
		bx = 0.5 * (u(k,ii,jj)+u(k,ii+1,jj))/2. * (svretax(i,j)+svretax(i,j+1))/2.
		by = 0.5 * (u(k,ii,jj)+u(k,ii+1,jj))/2. * (svretay(i,j)+svretay(i,j+1))/2.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fourth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j-1} S_{i,j-1} - u^{\xi}_{i+1,j-1} S_{i+1,j-1}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j-1} S_{i,j-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j-1)+vp(i,j-1))/2.
		beta2=(vp(i,j-1)+vp(i+1,j-1))/2.
		ax = ust_adjs(k,ii,jj-1)/beta1*suxix(i,j-1) - ust_adjs(k,ii+1,jj-1)/beta2*suxix(i+1,j-1)
		ay = ust_adjs(k,ii,jj-1)/beta1*suxiy(i,j-1) - ust_adjs(k,ii+1,jj-1)/beta2*suxiy(i+1,j-1)
		bx = 0.5 * (u(k,ii,jj-1)+u(k,ii+1,jj-1))/2. * (svretax(i,j-1)+svretax(i,j))/2.
		by = 0.5 * (u(k,ii,jj-1)+u(k,ii+1,jj-1))/2. * (svretay(i,j-1)+svretay(i,j))/2.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fifth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j,k} S_{i,j,k} - u^{\xi}_{i,j,k+1} S_{i,j,k+1}  adjoint
		! b = 1/8 * u^{z}_{i-1/2,j,k+1} S_{\eta,i-1/2,j+1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k,ii,jj)/beta*suxix(i,j) - ust_adjs(k+1,ii,jj)/beta*suxix(i,j)
		ay = ust_adjs(k,ii,jj)/beta*suxiy(i,j) - ust_adjs(k+1,ii,jj)/beta*suxiy(i,j)
		bx = 0.125 * (w(k+1,ii-1,jj)+w(k+1,ii,jj))/2. * (svretax(i-1,j)+svretax(i,j)+svretax(i-1,j+1)+svretax(i,j+1))/4.
		by = 0.125 * (w(k+1,ii-1,jj)+w(k+1,ii,jj))/2. * (svretay(i-1,j)+svretay(i,j)+svretay(i-1,j+1)+svretay(i,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sixth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j,k-1} S_{i,j,k-1} - u^{\xi}_{i,j,k} S_{i,j,k}  adjoint
		! b = 1/8 * u^{z}_{i-1/2,j,k} S_{\eta,i-1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k-1,ii,jj)/beta*suxix(i,j) - ust_adjs(k,ii,jj)/beta*suxix(i,j)
		ay = ust_adjs(k-1,ii,jj)/beta*suxiy(i,j) - ust_adjs(k,ii,jj)/beta*suxiy(i,j)
		bx = 0.125 * (w(k,ii-1,jj)+w(k,ii,jj))/2. * (svretax(i-1,j)+svretax(i,j)+svretax(i-1,j+1)+svretax(i,j+1))/4.
		by = 0.125 * (w(k,ii-1,jj)+w(k,ii,jj))/2. * (svretay(i-1,j)+svretay(i,j)+svretay(i-1,j+1)+svretay(i,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! seventh term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i+1,j,k} S_{i+1,j,k} - u^{\xi}_{i+1,j,k+1} S_{i+1,j,k+1}  adjoint
		! b = 1/8 * u^{z}_{i+1/2,j,k+1} S_{\eta,i+1/2,j+1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k,ii+1,jj)/beta*suxix(i+1,j) - ust_adjs(k+1,ii+1,jj)/beta*suxix(i+1,j)
		ay = ust_adjs(k,ii+1,jj)/beta*suxiy(i+1,j) - ust_adjs(k+1,ii+1,jj)/beta*suxiy(i+1,j)
		bx = 0.125 * (w(k+1,ii,jj)+w(k+1,ii+1,jj))/2. * (svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
		by = 0.125 * (w(k+1,ii,jj)+w(k+1,ii+1,jj))/2. * (svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! eighth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i+1,j,k-1} S_{i+1,j,k-1} - u^{\xi}_{i+1,j,k} S_{i+1,j,k}  adjoint
		! b = 1/8 * u^{z}_{i+1/2,j,k} S_{\eta,i+1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k-1,ii+1,jj)/beta*suxix(i+1,j) - ust_adjs(k,ii+1,jj)/beta*suxix(i+1,j)
		ay = ust_adjs(k-1,ii+1,jj)/beta*suxiy(i+1,j) - ust_adjs(k,ii+1,jj)/beta*suxiy(i+1,j)
		bx = 0.125 * (w(k,ii,jj)+w(k,ii+1,jj))/2. * (svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
		by = 0.125 * (w(k,ii,jj)+w(k,ii+1,jj))/2. * (svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! ninth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j-1,k} S_{i,j-1,k} - u^{\xi}_{i,j-1,k+1} S_{i,j-1,k+1}  adjoint
		! b = 1/8 * u^{z}_{i-1/2,j-1,k+1} S_{\eta,i-1/2,j-1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		ax = ust_adjs(k,ii,jj-1)/beta*suxix(i,j-1) - ust_adjs(k+1,ii,jj-1)/beta*suxix(i,j-1)
		ay = ust_adjs(k,ii,jj-1)/beta*suxiy(i,j-1) - ust_adjs(k+1,ii,jj-1)/beta*suxiy(i,j-1)
		bx = 0.125 * (w(k+1,ii-1,jj-1)+w(k+1,ii,jj-1))/2. * (svretax(i-1,j-1)+svretax(i,j-1)+svretax(i-1,j)+svretax(i,j))/4.
		by = 0.125 * (w(k+1,ii-1,jj-1)+w(k+1,ii,jj-1))/2. * (svretay(i-1,j-1)+svretay(i,j-1)+svretay(i-1,j)+svretay(i,j))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! tenth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i,j-1,k-1} S_{i,j-1,k-1} - u^{\xi}_{i,j-1,k} S_{i,j-1,k}  adjoint
		! b = 1/8 * u^{z}_{i-1/2,j-1,k} S_{\eta,i-1/2,j-1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		ax = ust_adjs(k-1,ii,jj-1)/beta*suxix(i,j-1) - ust_adjs(k,ii,jj-1)/beta*suxix(i,j-1)
		ay = ust_adjs(k-1,ii,jj-1)/beta*suxiy(i,j-1) - ust_adjs(k,ii,jj-1)/beta*suxiy(i,j-1)
		bx = 0.125 * (w(k,ii-1,jj-1)+w(k,ii,jj-1))/2. * (svretax(i-1,j-1)+svretax(i,j-1)+svretax(i-1,j)+svretax(i,j))/4.
		by = 0.125 * (w(k,ii-1,jj-1)+w(k,ii,jj-1))/2. * (svretay(i-1,j-1)+svretay(i,j-1)+svretay(i-1,j)+svretay(i,j))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! eleventh term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i+1,j-1,k} S_{i,j-1,k} - u^{\xi}_{i+1,j-1,k+1} S_{i+1,j-1,k+1}  adjoint
		! b = 1/8 * u^{z}_{i+1/2,j-1,k+1} S_{\eta,i+1/2,j-1/2,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		ax = ust_adjs(k,ii+1,jj-1)/beta*suxix(i+1,j-1) - ust_adjs(k+1,ii+1,jj-1)/beta*suxix(i+1,j-1)
		ay = ust_adjs(k,ii+1,jj-1)/beta*suxiy(i+1,j-1) - ust_adjs(k+1,ii+1,jj-1)/beta*suxiy(i+1,j-1)
		bx = 0.125 * (w(k+1,ii,jj-1)+w(k+1,ii+1,jj-1))/2. * (svretax(i,j-1)+svretax(i+1,j-1)+svretax(i,j)+svretax(i+1,j))/4.
		by = 0.125 * (w(k+1,ii,jj-1)+w(k+1,ii+1,jj-1))/2. * (svretay(i,j-1)+svretay(i+1,j-1)+svretay(i,j)+svretay(i+1,j))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! twelfth term h1 = h1 + a \bu b  in (A-180)
		! a = u^{\xi}_{i+1,j-1,k-1} S_{i+1,j-1,k-1} - u^{\xi}_{i+1,j-1,k} S_{i+1,j-1,k}  adjoint
		! b = 1/8 * u^{z}_{i+1/2,j-1,k} S_{\eta,i+1/2,j-1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		ax = ust_adjs(k-1,ii+1,jj-1)/beta*suxix(i+1,j-1) - ust_adjs(k,ii+1,jj-1)/beta*suxix(i+1,j-1)
		ay = ust_adjs(k-1,ii+1,jj-1)/beta*suxiy(i+1,j-1) - ust_adjs(k,ii+1,jj-1)/beta*suxiy(i+1,j-1)
		bx = 0.125 * (w(k,ii,jj-1)+w(k,ii+1,jj-1))/2. * (svretax(i,j-1)+svretax(i+1,j-1)+svretax(i,j)+svretax(i+1,j))/4.
		by = 0.125 * (w(k,ii,jj-1)+w(k,ii+1,jj-1))/2. * (svretay(i,j-1)+svretay(i+1,j-1)+svretay(i,j)+svretay(i+1,j))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\eta,2}
		!-----------------------------------------------------------------------
		! first term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{j} S_{j} - u^{\eta}_{j+1} S_{j+1}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j} S_{i+1/2,j} + u^{\eta}_{i,j+1/2} S_{i,j+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h2 = 0.0

		beta1=(vp(i,j-1)+vp(i,j))/2.
		beta2=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k,ii,jj)/beta1*svetax(i,j) - vst_adjs(k,ii,jj+1)/beta2*svetax(i,j+1)
		ay = vst_adjs(k,ii,jj)/beta1*svetay(i,j) - vst_adjs(k,ii,jj+1)/beta2*svetay(i,j+1)
		bx = 0.5 * (u(k,ii,jj)+u(k,ii+1,jj))/2. * (surxix(i,j)+surxix(i+1,j))/2.
		bx = bx + (v(k,ii,jj)+v(k,ii,jj+1))/2. * (svretax(i,j)+svretax(i,j+1))/2.
		by = 0.5 * (u(k,ii,jj)+u(k,ii+1,jj))/2. * (surxiy(i,j)+surxiy(i+1,j))/2.
		by = by + (v(k,ii,jj)+v(k,ii,jj+1))/2. * (svretay(i,j)+svretay(i,j+1))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{j-1} S_{j-1} - u^{\eta}_{j-1} S_{j-1}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j-1} S_{i+1/2,j-1} + u^{\eta}_{i,j-1/2} S_{i,j-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i,j-2)+vp(i,j-1))/2.
		beta2=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k,ii,jj-1)/beta1*svetax(i,j-1) - vst_adjs(k,ii,jj)/beta2*svetax(i,j)
		ay = vst_adjs(k,ii,jj-1)/beta1*svetay(i,j-1) - vst_adjs(k,ii,jj)/beta2*svetay(i,j)
		bx = 0.5 * (u(k,ii,jj-1)+u(k,ii+1,jj-1))/2. * (surxix(i,j-1)+surxix(i+1,j-1))/2.
		bx = bx + (v(k,ii,jj-1)+v(k,ii,jj))/2. * (svretax(i,j-1)+svretax(i,j))/2.
		by = 0.5 * (u(k,ii,jj-1)+u(k,ii+1,jj-1))/2. * (surxiy(i,j-1)+surxiy(i+1,j-1))/2.
		by = by + (v(k,ii,jj-1)+v(k,ii,jj))/2. * (svretay(i,j-1)+svretay(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! third term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{i,j} S_{i,j} - u^{\eta}_{i+1,j} S_{i+1,j}  adjoint
		! b = 1/2 u^{\xi}_{i+1,j-1/2,k}S_{\eta,i+1/2,j,k} forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i,j-1)+vp(i,j))/2.
		beta2=(vp(i+1,j-1)+vp(i+1,j))/2.
		ax = vst_adjs(k,ii,jj)/beta1*svetax(i,j) - vst_adjs(k,ii+1,jj)/beta2*svetax(i+1,j)
		ay = vst_adjs(k,ii,jj)/beta1*svetay(i,j) - vst_adjs(k,ii+1,jj)/beta2*svetay(i+1,j)
		bx = 0.5 * (u(k,ii+1,jj-1)+u(k,ii+1,jj))/2. * (svretax(i,j)+svretax(i+1,j))/2.
		by = 0.5 * (u(k,ii+1,jj-1)+u(k,ii+1,jj))/2. * (svretay(i,j)+svretay(i+1,j))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fourth term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{i-1,j} S_{i-1,j} - u^{\eta}_{i,j} S_{i,j}  adjoint
		! b = 1/2 u^{\xi}_{i,j-1/2,k}S_{\eta,i-1/2,j,k} forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=(vp(i-1,j-1)+vp(i-1,j))/2.
		beta2=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k,ii-1,jj)/beta1*svetax(i-1,j) - vst_adjs(k,ii,jj)/beta2*svetax(i,j)
		ay = vst_adjs(k,ii-1,jj)/beta1*svetay(i-1,j) - vst_adjs(k,ii,jj)/beta2*svetay(i,j)
		bx = 0.5 * (u(k,ii,jj-1)+u(k,ii,jj))/2. * (svretax(i-1,j)+svretax(i,j))/2.
		by = 0.5 * (u(k,ii,jj-1)+u(k,ii,jj))/2. * (svretay(i-1,j)+svretay(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! fifth term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{i,j,k} S_{i,j,k} - u^{\eta}_{i,j,k+1} S_{i,j,k+1}  adjoint
		! b = 1/2 u^{z}_{i,j-1/2,k+1}S_{\eta,i,j,k+1/2} forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k,ii,jj)/beta*svetax(i,j) - vst_adjs(k+1,ii,jj)/beta*svetax(i,j)
		ay = vst_adjs(k,ii,jj)/beta*svetay(i,j) - vst_adjs(k+1,ii,jj)/beta*svetay(i,j)
		bx = 0.5 * (w(k+1,ii,jj-1)+w(k+1,ii,jj))/2. * (svretax(i,j)+svretax(i,j))/2.
		by = 0.5 * (w(k+1,ii,jj-1)+w(k+1,ii,jj))/2. * (svretay(i,j)+svretay(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sixth term h2 = a \bu b  in (A-181)
		! a = u^{\eta}_{i,j,k-1} S_{i,j,k-1} - u^{\eta}_{i,j,k} S_{i,j,k}  adjoint
		! b = 1/2 u^{z}_{i,j-1/2,k}S_{\eta,i,j,k-1/2} forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k-1,ii,jj)/beta*svetax(i,j) - vst_adjs(k,ii,jj)/beta*svetax(i,j)
		ay = vst_adjs(k-1,ii,jj)/beta*svetay(i,j) - vst_adjs(k,ii,jj)/beta*svetay(i,j)
		bx = 0.5 * (w(k,ii,jj-1)+w(k,ii,jj))/2. * (svretax(i,j)+svretax(i,j))/2.
		by = 0.5 * (w(k,ii,jj-1)+w(k,ii,jj))/2. * (svretay(i,j)+svretay(i,j))/2.

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{\eta,3}
		!-----------------------------------------------------------------------
		! first term h3 = a \bu b  in (A-182)
		! a = u^{z}_{j-1} S_{j-1} - u^{z}_{j} S_{j}  adjoint
		! b = 1/2 * u^{z}_{i,j-1/2} S_{i,j-1/2}   forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h3 = 0.0
		
		beta1=vp(i,j-1)
		beta2=vp(i,j)
		az = wst_adjs(k,ii,jj-1)/beta1*swz(i,j-1) - wst_adjs(k,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (w(k,ii,jj-1)+w(k,ii,jj))/2. * (swrz(i,j-1)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h3 = a \bu b  in (A-182)
		! a = u^{z}_{j-1,k+1} S_{j-1,k+1} - u^{z}_{j,k+1} S_{j,k+1}  adjoint
		! b = 1/2 * u^{z}_{i,j-1/2,k+1} S_{i,j-1/2,k+1}   forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i,j-1)
		beta2=vp(i,j)
		az = wst_adjs(k+1,ii,jj-1)/beta1*swz(i,j-1) - wst_adjs(k+1,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (w(k+1,ii,jj-1)+w(k+1,ii,jj))/2. * (swrz(i,j-1)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sum of all terms
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		con_adj = -h1 - h2 - h3

		v_adj(k,ii,jj) = v_adj(k,ii,jj) + dt*con_adj
	end do
	end do
	end do

end

subroutine advancez_Adjoint_conv()

	use advance_module

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)		
       		vst_adjs(:, :, 1) =  0.0			
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)		

       	      	!Extend
		vst_adjs(:, :, 0) = vst_adjs(:, :, 2)		! this halo is not included in following codes
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)	
              	vst_adjs(:, :, nlyb  ) =  0.0			
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)	
		
		!Extend
		vst_adjs(:, :, nlyb+1) = vst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
     		
	end if 

	!ccccccccccccccccccccc
	!   for w_adj equation
	!ccccccccccccccccccccc
	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{z,1}
		!-----------------------------------------------------------------------
		! first term h1 = a \bu b  in (A-184)
		! a = u^{\xi}_{j,k-1} S_{j,k-1} - u^{\xi}_{j,k} S_{j,k}  adjoint
		! b = 1/2 * u^{\xi}_{i,j,k-1/2} S_{i,j,k-1/2} + u^{\eta}_{i-1/2,j+1/2,k-1/2} S_{i-1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h1 = 0.0

		beta=(vp(i-1,j)+vp(i,j))/2.
		ax = ust_adjs(k-1,ii,jj)/beta*suxix(i,j) - ust_adjs(k,ii,jj)/beta*suxix(i,j)
		ay = ust_adjs(k-1,ii,jj)/beta*suxiy(i,j) - ust_adjs(k,ii,jj)/beta*suxiy(i,j)
		bx = 0.5*(u(k-1,ii,jj)+u(k,ii,jj))/2.*surxix(i,j)
		bx = bx + 0.5*(v(k-1,ii-1,jj)+v(k,ii-1,jj)+v(k-1,ii,jj)+v(k,ii,jj)&
			+v(k-1,ii-1,jj+1)+v(k,ii-1,jj+1)+v(k-1,ii,jj+1)+v(k,ii,jj+1))/8. &
			* (svretax(i-1,j)+svretax(i,j)+svretax(i-1,j+1)+svretax(i,j+1))/4.
		by = 0.5*(u(k-1,ii,jj)+u(k,ii,jj))/2.*surxiy(i,j)
		by = by + 0.5*(v(k-1,ii-1,jj)+v(k,ii-1,jj)+v(k-1,ii,jj)+v(k,ii,jj)&
			+v(k-1,ii-1,jj+1)+v(k,ii-1,jj+1)+v(k-1,ii,jj+1)+v(k,ii,jj+1))/8. &
			* (svretay(i-1,j)+svretay(i,j)+svretay(i-1,j+1)+svretay(i,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h1 = a \bu b  in (A-184)
		! a = u^{\xi}_{i+1,j,k-1} S_{i+1,j,k-1} - u^{\xi}_{i+1,j,k} S_{i+1,j,k}  adjoint
		! b = 1/2 * u^{\xi}_{i+1,j,k-1/2} S_{i+1,j,k-1/2} + u^{\eta}_{i+1/2,j+1/2,k-1/2} S_{i+1/2,j+1/2,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i+1,j))/2.
		ax = ust_adjs(k-1,ii+1,jj)/beta*suxix(i+1,j) - ust_adjs(k,ii+1,jj)/beta*suxix(i+1,j)
		ay = ust_adjs(k-1,ii+1,jj)/beta*suxiy(i+1,j) - ust_adjs(k,ii+1,jj)/beta*suxiy(i+1,j)
		bx = 0.5*(u(k-1,ii+1,jj)+u(k,ii+1,jj))/2.*surxix(i+1,j)
		bx = bx + 0.5*(v(k-1,ii,jj)+v(k,ii,jj)+v(k-1,ii+1,jj)+v(k,ii+1,jj)&
			+v(k-1,ii,jj+1)+v(k,ii,jj+1)+v(k-1,ii+1,jj+1)+v(k,ii+1,jj+1))/8. &
			* (svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
		by = 0.5*(u(k-1,ii+1,jj)+u(k,ii+1,jj))/2.*surxiy(i+1,j)
		by = by + 0.5*(v(k-1,ii,jj)+v(k,ii,jj)+v(k-1,ii+1,jj)+v(k,ii+1,jj)&
			+v(k-1,ii,jj+1)+v(k,ii,jj+1)+v(k-1,ii+1,jj+1)+v(k,ii+1,jj+1))/8. &
			* (svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.

		h1 = h1 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{z,2}
		!-----------------------------------------------------------------------
		! first term h2 = a \bu b  in (A-185)
		! a = u^{\eta}_{j,k-1} S_{j,k-1} - u^{\eta}_{j,k} S_{j,k}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j-1/2,k-1/2} S_{i+1/2,j-1/2,k-1/2} + u^{\eta}_{i,j,k-1/2} S_{i,j,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h2 = 0.0

		beta=(vp(i,j-1)+vp(i,j))/2.
		ax = vst_adjs(k-1,ii,jj)/beta*svetax(i,j) - vst_adjs(k,ii,jj)/beta*svetax(i,j)
		ay = vst_adjs(k-1,ii,jj)/beta*svetay(i,j) - vst_adjs(k,ii,jj)/beta*svetay(i,j)
		bx = 0.5*(u(k-1,ii,jj-1)+u(k,ii,jj-1)+u(k-1,ii+1,jj-1)+u(k,ii+1,jj-1)&
			+u(k-1,ii,jj)+u(k,ii,jj)+u(k-1,ii+1,jj)+u(k,ii+1,jj))/8. &
			* (surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j)+surxix(i+1,j))/4.
		bx = bx + 0.5*(v(k-1,ii,jj)+v(k,ii,jj))/2.*svretax(i,j)
		by = 0.5*(u(k-1,ii,jj-1)+u(k,ii,jj-1)+u(k-1,ii+1,jj-1)+u(k,ii+1,jj-1)&
			+u(k-1,ii,jj)+u(k,ii,jj)+u(k-1,ii+1,jj)+u(k,ii+1,jj))/8. &
			* (surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j)+surxiy(i+1,j))/4.
		by = by + 0.5*(v(k-1,ii,jj)+v(k,ii,jj))/2.*svretay(i,j)

		h2 = h2 + ax*bx + ay*by

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! second term h2 = a \bu b  in (A-185)
		! a = u^{\eta}_{j+1,k-1} S_{j+1,k-1} - u^{\eta}_{j+1,k} S_{j+1,k}  adjoint
		! b = 1/2 * u^{\xi}_{i+1/2,j+1/2,k-1/2} S_{i+1/2,j+1/2,k-1/2} + u^{\eta}_{i,j+1,k-1/2} S_{i,j+1,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=(vp(i,j)+vp(i,j+1))/2.
		ax = vst_adjs(k-1,ii,jj+1)/beta*svetax(i,j+1) - vst_adjs(k,ii,jj+1)/beta*svetax(i,j+1)
		ay = vst_adjs(k-1,ii,jj+1)/beta*svetay(i,j+1) - vst_adjs(k,ii,jj+1)/beta*svetay(i,j+1)
		bx = 0.5*(u(k-1,ii,jj)+u(k,ii,jj)+u(k-1,ii+1,jj)+u(k,ii+1,jj)&
			+u(k-1,ii,jj+1)+u(k,ii,jj+1)+u(k-1,ii+1,jj+1)+u(k,ii+1,jj+1))/8. &
			* (surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		bx = bx + 0.5*(v(k-1,ii,jj+1)+v(k,ii,jj+1))/2.*svretax(i,j+1)
		by = 0.5*(u(k-1,ii,jj)+u(k,ii,jj)+u(k-1,ii+1,jj)+u(k,ii+1,jj)&
			+u(k-1,ii,jj+1)+u(k,ii,jj+1)+u(k-1,ii+1,jj+1)+u(k,ii+1,jj+1))/8. &
			* (surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
		by = by + 0.5*(v(k-1,ii,jj+1)+v(k,ii,jj+1))/2.*svretay(i,j+1)

		h2 = h2 + ax*bx + ay*by
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! convection terms  H^{\dag}_{z,3}
		!-----------------------------------------------------------------------
		! first term h3 = a \bu b  in (A-186)
		! a = u^{z}_{j,k} S_{j,k} - u^{z}_{j,k+1} S_{j,k+1}  adjoint
		! b = u^{z}_{i,j,k+1/2}S_{z,i,j,k+1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		h3 = 0.0
		beta=vp(i,j)
		az = wst_adjs(k,ii,jj)/beta*swz(i,j) - wst_adjs(k+1,ii,jj)/beta*swz(i,j)
		bz = (w(k,ii,jj)+w(k+1,ii,jj))/2. * swrz(i,j)

		h3 = h3 + az*bz
		!-----------------------------------------------------------------------
		! second term h3 = a \bu b  in (A-186)
		! a = u^{z}_{j,k-1} S_{j,k-1} - u^{z}_{j,k} S_{j,k}  adjoint
		! b = u^{z}_{i,j,k-1/2}S_{z,i,j,k-1/2}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta=vp(i,j)
		az = wst_adjs(k-1,ii,jj)/beta*swz(i,j) - wst_adjs(k,ii,jj)/beta*swz(i,j)
		bz = (w(k-1,ii,jj)+w(k,ii,jj))/2. * swrz(i,j)

		h3 = h3 + az*bz
		!-----------------------------------------------------------------------
		! third term h3 = a \bu b  in (A-186)
		! a = u^{z}_{j,k} S_{j,k} - u^{z}_{j+1,k} S_{j+1,k}  adjoint
		! b = u^{\eta}_{i,j+1,k-1/2}S_{z,i,j+1/2,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i,j)
		beta2=vp(i,j+1)
		az = wst_adjs(k,ii,jj)/beta1*swz(i,j) - wst_adjs(k,ii,jj+1)/beta2*swz(i,j+1)
		bz = 0.5 * (v(k-1,ii,jj+1)+v(k,ii,jj+1))/2. * (swrz(i,j)+swrz(i,j+1))/2.

		h3 = h3 + az*bz
		!-----------------------------------------------------------------------
		! fourth term h3 = a \bu b  in (A-186)
		! a = u^{z}_{j-1,k} S_{j-1,k} - u^{z}_{j,k} S_{j,k}  adjoint
		! b = u^{\eta}_{i,j,k-1/2}S_{z,i,j-1/2,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i,j-1)
		beta2=vp(i,j)
		az = wst_adjs(k,ii,jj-1)/beta1*swz(i,j-1) - wst_adjs(k,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (v(k-1,ii,jj)+v(k,ii,jj))/2. * (swrz(i,j-1)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!-----------------------------------------------------------------------
		! fifth term h3 = a \bu b  in (A-186)
		! a = u^{z}_{i,j,k} S_{i,j,k} - u^{z}_{i+1,j,k} S_{i+1,j,k}  adjoint
		! b = u^{\xi}_{i+1,j,k-1/2}S_{z,i+1/2,j,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i,j)
		beta2=vp(i+1,j)
		az = wst_adjs(k,ii,jj)/beta1*swz(i,j) - wst_adjs(k,ii+1,jj)/beta2*swz(i+1,j)
		bz = 0.5 * (u(k-1,ii+1,jj)+u(k,ii+1,jj))/2. * (swrz(i,j)+swrz(i+1,j))/2.

		h3 = h3 + az*bz
		!-----------------------------------------------------------------------
		! sixth term h3 = a \bu b  in (A-186)
		! a = u^{z}_{i-1,j,k} S_{i-1,j,k} - u^{z}_{i,j,k} S_{i,j,k}  adjoint
		! b = u^{\xi}_{i,j,k-1/2}S_{z,i-1/2,j,k}  forward
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		beta1=vp(i-1,j)
		beta2=vp(i,j)
		az = wst_adjs(k,ii-1,jj)/beta1*swz(i-1,j) - wst_adjs(k,ii,jj)/beta2*swz(i,j)
		bz = 0.5 * (u(k-1,ii,jj)+u(k,ii,jj))/2. * (swrz(i-1,j)+swrz(i,j))/2.

		h3 = h3 + az*bz
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! sum of all terms
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		con_adj = -h1 - h2 - h3

		w_adj(k,ii,jj) = w_adj(k,ii,jj) + dt*con_adj
	end do
	end do
	end do


end

subroutine advancex_Adjoint_diffex()
	use advance_module

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
        	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)
       		vst_adjs(:, :, 1) =  0.0
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)

       	      	!Extend
		vst_adjs(:, :, 0) = -vst_adjs(:, :, 2)	
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)
              	vst_adjs(:, :, nlyb  ) =  0.0
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)
		
		!Extend
		vst_adjs(:, :, nlyb+1) = -vst_adjs(:, :, nlyb-1)	
     		
	end if 
	
	!ccccccccccccccccccccc
	!   for u_adj equation
	!ccccccccccccccccccccc
	do jj=j1_u,j2_u 
	do ii=i1_u,i2_u
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! explicit diffusion term  D^{\dag}_{\xi,1,ex}
		!-----------------------------------------------------------------------
		! D^{\dag}_{\xi,P,ex}, S^{\xi}_{i,j,k}S_{\xi,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=suxix(i,j)
		a2=suxiy(i,j)
		b1=surxix(i,j)
		b2=surxiy(i,j)
		c1=a1*b1
		c2=a1*b2
		c4=a2*b1
		c5=a2*b2

		!-----------------------------------------------------------------------
		! (u^{\eta}_{i-1,j,k}S^{\eta}_{i-1,j,k}-u^{\eta}_{i-1,j+1,k}S^{\eta}_{i-1,j+1,k})
		! \bu S^{\eta}_{i-1,j+1/2,k}/(Re V_{i-1,j,k})
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i-1,j+1))/2. 		
		a2=(svetay(i-1,j)+svetay(i-1,j+1))/2. 
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.	
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta

		diffex_adj = f1*e1 + f2*e2
		!-----------------------------------------------------------------------
		! (u^{\eta}_{i,j,k}S^{\eta}_{i,j,k}-u^{\eta}_{i,j+1,k}S^{\eta}_{i,j+1,k})
		! \bu S^{\eta}_{i,j+1/2,k}/(Re V_{i,j,k})
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i,j+1))/2. 		
		a2=(svetay(i,j)+svetay(i,j+1))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.		
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! S^{\eta}_{i-1/2,j+1/2,k}S_{\xi,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i-1,j)+svetax(i,j)+svetax(i-1,j+1)+svetax(i,j+1))/4.
		a2=(svetay(i-1,j)+svetay(i,j)+svetay(i-1,j+1)+svetay(i,j+1))/4.
		b1=surxix(i,j)
		b2=surxiy(i,j)
		c1=a1*b1
		c2=a1*b2
		c4=a2*b1
		c5=a2*b2

		!-----------------------------------------------------------------------
		! (u^{\eta}_{i-1,j,k}S^{\eta}_{i-1,j,k}-u^{\eta}_{i,j,k}S^{\eta}_{i,j,k})
		! \bu S^{\xi}_{i,j-1/2,k}/(Re V_{i-1/2,j-1/2,k})
		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i,j))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i,j-1)+vp(i-1,j)+vp(i,j))/4.
		e1=  (a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=  (a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		! (u^{\eta}_{i-1,j+1,k}S^{\eta}_{i-1,j+1,k}-u^{\eta}_{i,j+1,k}S^{\eta}_{i,j+1,k})
		! \bu S^{\xi}_{i,j+1/2,k}/(Re V_{i-1/2,j+1/2,k})
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i,j+1))/2. 		
		a2=(suxiy(i,j)+suxiy(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i,j)+vp(i-1,j+1)+vp(i,j+1))/4.
		e1=  (a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=  (a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! S^{z}_{i-1/2,j,k+1/2}S_{\xi,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=surxix(i,j)
		a2=surxiy(i,j)
		b3=(swz(i-1,j)+swz(i,j))/2.
		c1=a1*b3
		c2=a2*b3

		a1=suxix(i,j)
		a2=suxiy(i,j)
		a3=(vp(i-1,j)+vp(i,j))/2.
		e3= (c1*a1+c2*a2)*visc/a3

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i-1,j)
		f3=f3 - wst_adjs(k+1,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 + wst_adjs(k+1,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,ne,ex}, S^{\xi}_{i+1/2,j+1/2,k}S_{\xi,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j+1)+suxix(i+1,j+1))/4.
		a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j+1)+suxiy(i+1,j+1))/4.
		b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(svetax(i-1,j+1)+svetax(i,j+1))/2. 		
		a2=(svetay(i-1,j+1)+svetay(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i,j)+vp(i-1,j+1)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		f1=f1 - ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=f2 - ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i,j+1))/2. 		
		a2=(suxiy(i,j)+suxiy(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i,j)+vp(i-1,j+1)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(svetax(i,j+1)+svetax(i+1,j+1))/2. 		
		a2=(svetay(i,j+1)+svetay(i+1,j+1))/2. 	
		a3=(vp(i,j)+vp(i+1,j)+vp(i,j+1)+vp(i+1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta
		beta=(vp(i,j+1)+vp(i+1,j+1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj+1)*suxix(i+1,j+1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj+1)*suxiy(i+1,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i+1,j)+suxix(i+1,j+1))/2. 		
		a2=(suxiy(i+1,j)+suxiy(i+1,j+1))/2. 	
		a3=(vp(i,j)+vp(i+1,j)+vp(i,j+1)+vp(i+1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta
		beta=(vp(i+1,j)+vp(i+1,j+1))/2.
		f1=f1 - vst_adjs(k,ii+1,jj+1)*svetax(i+1,j+1)/beta
		f2=f2 - vst_adjs(k,ii+1,jj+1)*svetay(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,ne,ex}, S^{\eta}_{i,j+1,k}S_{\xi,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svetax(i,j+1)
		a2=svetay(i,j+1)
		b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i+1,j))/2. 		
		a2=(suxiy(i,j)+suxiy(i+1,j))/2.
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.	
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i,j+1))/2. 		
		a2=(svetay(i,j)+svetay(i,j+1))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.	
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j+1)+suxix(i+1,j+1))/2. 		
		a2=(suxiy(i,j+1)+suxiy(i+1,j+1))/2.
		a3=vp(i,j+1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		f1=ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta
		beta=(vp(i,j+1)+vp(i+1,j+1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj+1)*suxix(i+1,j+1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj+1)*suxiy(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j+1)+svetax(i,j+2))/2. 		
		a2=(svetay(i,j+1)+svetay(i,j+2))/2. 
		a3=vp(i,j+1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta
		beta=(vp(i,j+1)+vp(i,j+2))/2.	
		f1=f1 - vst_adjs(k,ii,jj+2)*svetax(i,j+2)/beta
		f2=f2 - vst_adjs(k,ii,jj+2)*svetay(i,j+2)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,ne,ex}, S^{z}_{i,j+1/2,k+1/2}S_{\xi,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i,j)+swz(i,j+1))/2.
		b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j)+swz(i,j+1))/2.
		a3=(vp(i,j)+vp(i,j+1))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k-1,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k-1,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 + vst_adjs(k+1,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 + vst_adjs(k+1,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\xi,i+1/2,j+1/2,k}S^{z}_{i,j+1/2,k+1/2}------------
		a1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
		a2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
		b3=(swz(i,j)+swz(i,j+1))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=svetax(i,j+1)
		a2=svetay(i,j+1)
		a3=(vp(i,j)+vp(i,j+1))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i,j+1)
		f3=f3 - wst_adjs(k,ii,jj+1)*swz(i,j+1)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k+1,ii,jj)*swz(i,j)/beta
		beta=vp(i,j+1)
		f3=f3 + wst_adjs(k+1,ii,jj+1)*swz(i,j+1)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,nw,ex}, S^{\xi}_{i-1/2,j+1/2,k}S_{\xi,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i-1,j)+suxix(i,j)+suxix(i-1,j+1)+suxix(i,j+1))/4.
		a2=(suxiy(i-1,j)+suxiy(i,j)+suxiy(i-1,j+1)+suxiy(i,j+1))/4.
		b1=(surxix(i-1,j)+surxix(i,j)+surxix(i-1,j+1)+surxix(i,j+1))/4.
		b2=(surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j+1)+surxiy(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(svetax(i-2,j+1)+svetax(i-1,j+1))/2. 		
		a2=(svetay(i-2,j+1)+svetay(i-1,j+1))/2. 	
		a3=(vp(i-2,j)+vp(i-1,j)+vp(i-2,j+1)+vp(i-1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-2,j)+vp(i-1,j))/2.
		f1=ust_adjs(k,ii-1,jj)*suxix(i-1,j)/beta
		f2=ust_adjs(k,ii-1,jj)*suxiy(i-1,j)/beta
		beta=(vp(i-2,j+1)+vp(i-1,j+1))/2.
		f1=f1 - ust_adjs(k,ii-1,jj+1)*suxix(i-1,j+1)/beta
		f2=f2 - ust_adjs(k,ii-1,jj+1)*suxiy(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j)+suxix(i-1,j+1))/2. 		
		a2=(suxiy(i-1,j)+suxiy(i-1,j+1))/2. 	
		a3=(vp(i-2,j)+vp(i-1,j)+vp(i-2,j+1)+vp(i-1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j)+vp(i-2,j+1))/2.
		f1=vst_adjs(k,ii-2,jj+1)*svetax(i-2,j+1)/beta
		f2=vst_adjs(k,ii-2,jj+1)*svetay(i-2,j+1)/beta
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(svetax(i-1,j+1)+svetax(i,j+1))/2. 		
		a2=(svetay(i-1,j+1)+svetay(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i,j)+vp(i-1,j+1)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		f1=f1 - ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=f2 - ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i,j+1))/2. 		
		a2=(suxiy(i,j)+suxiy(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i,j)+vp(i-1,j+1)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,nw,ex}, S^{\eta}_{i-1,j+1,k}S_{\xi,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svetax(i-1,j+1)
		a2=svetay(i-1,j+1)
		b1=(surxix(i-1,j)+surxix(i,j)+surxix(i-1,j+1)+surxix(i,j+1))/4.
		b2=(surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j+1)+surxiy(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j)+suxix(i,j))/2. 		
		a2=(suxiy(i-1,j)+suxiy(i,j))/2.
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j)+vp(i-1,j))/2.
		f1=ust_adjs(k,ii-1,jj)*suxix(i-1,j)/beta
		f2=ust_adjs(k,ii-1,jj)*suxiy(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i-1,j+1))/2. 		
		a2=(svetay(i-1,j)+svetay(i-1,j+1))/2. 
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j+1)+suxix(i,j+1))/2. 		
		a2=(suxiy(i-1,j+1)+suxiy(i,j+1))/2.
		a3=vp(i-1,j+1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j+1)+vp(i-1,j+1))/2.
		f1=ust_adjs(k,ii-1,jj+1)*suxix(i-1,j+1)/beta
		f2=ust_adjs(k,ii-1,jj+1)*suxiy(i-1,j+1)/beta
		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		f1=f1 - ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=f2 - ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j+1)+svetax(i-1,j+2))/2. 		
		a2=(svetay(i-1,j+1)+svetay(i-1,j+2))/2. 
		a3=vp(i-1,j+1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		beta=(vp(i-1,j+1)+vp(i-1,j+2))/2.
		f1=f1 - vst_adjs(k,ii-1,jj+2)*svetax(i-1,j+2)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+2)*svetay(i-1,j+2)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,nw,ex}, S^{z}_{i-1,j+1/2,k+1/2}S_{\xi,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i-1,j)+swz(i-1,j+1))/2.
		b1=(surxix(i-1,j)+surxix(i,j)+surxix(i-1,j+1)+surxix(i,j+1))/4.
		b2=(surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j+1)+surxiy(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j)+swz(i-1,j+1))/2.
		a3=(vp(i-1,j)+vp(i-1,j+1))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k-1,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k-1,ii-1,jj+1)*svetay(i-1,j+1)/beta
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		f1=f1 + vst_adjs(k+1,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 + vst_adjs(k+1,ii-1,jj+1)*svetay(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\xi,i-1/2,j+1/2,k}S^{z}_{i-1,j+1/2,k+1/2}------------
		a1=(surxix(i-1,j)+surxix(i,j)+surxix(i-1,j+1)+surxix(i,j+1))/4.
		a2=(surxiy(i-1,j)+surxiy(i,j)+surxiy(i-1,j+1)+surxiy(i,j+1))/4.
		b3=(swz(i-1,j)+swz(i-1,j+1))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=svetax(i-1,j+1)
		a2=svetay(i-1,j+1)
		a3=(vp(i-1,j)+vp(i-1,j+1))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i-1,j+1)
		f3=f3 - wst_adjs(k,ii-1,jj+1)*swz(i-1,j+1)/beta
		beta=vp(i-1,j)
		f3=f3 - wst_adjs(k+1,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i-1,j+1)
		f3=f3 + wst_adjs(k+1,ii-1,jj+1)*swz(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,se,ex}, S^{\xi}_{i+1/2,j-1/2,k}S_{\xi,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i,j-1)+suxix(i+1,j-1)+suxix(i,j)+suxix(i+1,j))/4.
		a2=(suxiy(i,j-1)+suxiy(i+1,j-1)+suxiy(i,j)+suxiy(i+1,j))/4.
		b1=(surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j)+surxix(i+1,j))/4.
		b2=(surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j)+surxiy(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(svetax(i-1,j)+svetax(i,j))/2. 		
		a2=(svetay(i-1,j)+svetay(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i,j-1)+vp(i-1,j)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i,j))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i,j-1)+vp(i-1,j)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(svetax(i,j)+svetax(i+1,j))/2. 		
		a2=(svetay(i,j)+svetay(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i+1,j-1)+vp(i,j)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i+1,j-1)+suxix(i+1,j))/2. 		
		a2=(suxiy(i+1,j-1)+suxiy(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i+1,j-1)+vp(i,j)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.		
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i+1,j-1)+vp(i+1,j))/2.	
		f1=f1 - vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=f2 - vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,se,ex}, S^{\eta}_{i,j,k}S_{\xi,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svetax(i,j)
		a2=svetay(i,j)
		b1=(surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j)+surxix(i+1,j))/4.
		b2=(surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j)+surxiy(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i+1,j-1))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i+1,j-1))/2.
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j-1)+svetax(i,j))/2. 		
		a2=(svetay(i,j-1)+svetay(i,j))/2. 
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		if (j-2<0) then
			beta=(vp(i,2)+vp(i,j-1))/2.
		else 
			beta=(vp(i,j-2)+vp(i,j-1))/2.
		end if
		f1=vst_adjs(k,ii,jj-1)*svetax(i,j-1)/beta
		f2=vst_adjs(k,ii,jj-1)*svetay(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i+1,j))/2. 		
		a2=(suxiy(i,j)+suxiy(i+1,j))/2.
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i,j+1))/2. 		
		a2=(svetay(i,j)+svetay(i,j+1))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,se,ex}, S^{z}_{i,j-1/2,k+1/2}S_{\xi,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i,j-1)+swz(i,j))/2.
		b1=(surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j)+surxix(i+1,j))/4.
		b2=(surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j)+surxiy(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j-1)+swz(i,j))/2.
		a3=(vp(i,j-1)+vp(i,j))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k-1,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k-1,ii,jj)*svetay(i,j)/beta
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta
		f1=f1 + vst_adjs(k+1,ii,jj)*svetax(i,j)/beta
		f2=f2 + vst_adjs(k+1,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\xi,i+1/2,j-1/2,k}S^{z}_{i,j-1/2,k+1/2}------------
		a1=(surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j)+surxix(i+1,j))/4.
		a2=(surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j)+surxiy(i+1,j))/4.
		b3=(swz(i,j-1)+swz(i,j))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=svetax(i,j)
		a2=svetay(i,j)
		a3=(vp(i,j-1)+vp(i,j))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i,j-1)
		f3=f3 - wst_adjs(k+1,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 + wst_adjs(k+1,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,sw,ex}, S^{\xi}_{i-1/2,j-1/2,k}S_{\xi,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i-1,j-1)+suxix(i,j-1)+suxix(i-1,j)+suxix(i,j))/4.
		a2=(suxiy(i-1,j-1)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i,j))/4.
		b1=(surxix(i-1,j-1)+surxix(i,j-1)+surxix(i-1,j)+surxix(i,j))/4.
		b2=(surxiy(i-1,j-1)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(svetax(i-2,j)+svetax(i-1,j))/2. 		
		a2=(svetay(i-2,j)+svetay(i-1,j))/2. 	
		a3=(vp(i-2,j-1)+vp(i-1,j-1)+vp(i-2,j)+vp(i-1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-2,j-1)+vp(i-1,j-1))/2.
		f1=ust_adjs(k,ii-1,jj-1)*suxix(i-1,j-1)/beta
		f2=ust_adjs(k,ii-1,jj-1)*suxiy(i-1,j-1)/beta
		beta=(vp(i-2,j)+vp(i-1,j))/2.
		f1=f1 - ust_adjs(k,ii-1,jj)*suxix(i-1,j)/beta
		f2=f2 - ust_adjs(k,ii-1,jj)*suxiy(i-1,j)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j-1)+suxix(i-1,j))/2. 		
		a2=(suxiy(i-1,j-1)+suxiy(i-1,j))/2. 	
		a3=(vp(i-2,j-1)+vp(i-1,j-1)+vp(i-2,j)+vp(i-1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j-1)+vp(i-2,j))/2.
		f1=vst_adjs(k,ii-2,jj)*svetax(i-2,j)/beta
		f2=vst_adjs(k,ii-2,jj)*svetay(i-2,j)/beta
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(svetax(i-1,j)+svetax(i,j))/2. 		
		a2=(svetay(i-1,j)+svetay(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i,j-1)+vp(i-1,j)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i,j))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i,j-1)+vp(i-1,j)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,sw,ex}, S^{\eta}_{i-1,j,k}S_{\xi,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svetax(i-1,j)
		a2=svetay(i-1,j)
		b1=(surxix(i-1,j-1)+surxix(i,j-1)+surxix(i-1,j)+surxix(i,j))/4.
		b2=(surxiy(i-1,j-1)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j-1)+suxix(i,j-1))/2. 		
		a2=(suxiy(i-1,j-1)+suxiy(i,j-1))/2.
		a3=vp(i-1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j-1)+vp(i-1,j-1))/2.
		f1=ust_adjs(k,ii-1,jj-1)*suxix(i-1,j-1)/beta
		f2=ust_adjs(k,ii-1,jj-1)*suxiy(i-1,j-1)/beta
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j-1)+svetax(i-1,j))/2. 		
		a2=(svetay(i-1,j-1)+svetay(i-1,j))/2. 
		a3=vp(i-1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		if (j<2) then
			beta=(vp(i-1,2)+vp(i-1,j-1))/2.
		else
			beta=(vp(i-1,j-2)+vp(i-1,j-1))/2.
		end if
		f1=vst_adjs(k,ii-1,jj-1)*svetax(i-1,j-1)/beta
		f2=vst_adjs(k,ii-1,jj-1)*svetay(i-1,j-1)/beta
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j)+suxix(i,j))/2. 		
		a2=(suxiy(i-1,j)+suxiy(i,j))/2.
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-2,j)+vp(i-1,j))/2.
		f1=ust_adjs(k,ii-1,jj)*suxix(i-1,j)/beta
		f2=ust_adjs(k,ii-1,jj)*suxiy(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i-1,j+1))/2. 		
		a2=(svetay(i-1,j)+svetay(i-1,j+1))/2. 
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\xi,sw,ex}, S^{z}_{i-1,j-1/2,k+1/2}S_{\xi,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i-1,j-1)+swz(i-1,j))/2.
		b1=(surxix(i-1,j-1)+surxix(i,j-1)+surxix(i-1,j)+surxix(i,j))/4.
		b2=(surxiy(i-1,j-1)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j-1)+swz(i-1,j))/2.
		a3=(vp(i-1,j-1)+vp(i-1,j))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k-1,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k-1,ii-1,jj)*svetay(i-1,j)/beta
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		f1=f1 + vst_adjs(k+1,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 + vst_adjs(k+1,ii-1,jj)*svetay(i-1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\xi,i-1/2,j-1/2,k}S^{z}_{i-1,j-1/2,k+1/2}------------
		a1=(surxix(i-1,j-1)+surxix(i,j-1)+surxix(i-1,j)+surxix(i,j))/4.
		a2=(surxiy(i-1,j-1)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i,j))/4.
		b3=(swz(i-1,j-1)+swz(i-1,j))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=svetax(i-1,j)
		a2=svetay(i-1,j)
		a3=(vp(i-1,j-1)+vp(i-1,j))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i-1,j-1)
		f3=wst_adjs(k,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i-1,j)
		f3=f3 - wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i-1,j-1)
		f3=f3 - wst_adjs(k+1,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i-1,j)
		f3=f3 + wst_adjs(k+1,ii-1,jj)*swz(i-1,j)/beta

		diffex_adj = diffex_adj + (f3*e3)

		u_adj(k,ii,jj) = u_adj(k,ii,jj) + dt*diffex_adj
	end do
	end do
	end do

end

subroutine advancey_Adjoint_diffex()
	use advance_module
	
	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)
       		vst_adjs(:, :, 1) =  0.0			
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)		! this halo is not included in following codes

       	      	!Extend
		vst_adjs(:, :, 0) = vst_adjs(:, :, 2)		! this halo is not included in following codes
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)
              	vst_adjs(:, :, nlyb  ) =  0.0			
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
		
		!Extend
		vst_adjs(:, :, nlyb+1) = vst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
     		
	end if 

	!ccccccccccccccccccccc
	!   for v_adj equation
	!ccccccccccccccccccccc
	do jj=j1_v,j2_v 
	do ii=i1_v,i2_v
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! explicit diffusion term  D^{\dag}_{\eta,1,ex}
		!-----------------------------------------------------------------------
		! D^{\dag}_{\eta,P,ex}, S^{\eta}_{i,j,k}S_{\eta,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svetax(i,j)
		a2=svetay(i,j)
		b1=svretax(i,j)
		b2=svretay(i,j)
		c1=a1*b1
		c2=a1*b2
		c4=a2*b1
		c5=a2*b2

		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i+1,j-1))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i+1,j-1))/2. 
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.	
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta

		diffex_adj = f1*e1 + f2*e2
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i+1,j))/2. 		
		a2=(suxiy(i,j)+suxiy(i+1,j))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.		
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! S^{\xi}_{i+1/2,j-1/2,k}S_{\eta,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i,j-1)+suxix(i,j)+suxix(i+1,j-1)+suxix(i+1,j))/4.
		a2=(suxiy(i,j-1)+suxiy(i,j)+suxiy(i+1,j-1)+suxiy(i+1,j))/4.
		b1=svretax(i,j)
		b2=svretay(i,j)
		c1=a1*b1
		c2=a1*b2
		c4=a2*b1
		c5=a2*b2

		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i,j))/2. 		
		a2=(svetay(i-1,j)+svetay(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i-1,j)+vp(i,j-1)+vp(i,j))/4.
		e1=  (a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=  (a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i+1,j))/2. 		
		a2=(svetay(i,j)+svetay(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i,j)+vp(i+1,j-1)+vp(i+1,j))/4.
		e1=  (a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=  (a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! S^{z}_{i,j-1/2,k+1/2}S_{\eta,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=svretax(i,j)
		a2=svretay(i,j)
		b3=(swz(i,j-1)+swz(i,j))/2.
		c1=a1*b3
		c2=a2*b3

		a1=svetax(i,j)
		a2=svetay(i,j)
		a3=(vp(i,j-1)+vp(i,j))/2.
		e3= (c1*a1+c2*a2)*visc/a3

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i,j-1)
		f3=f3 - wst_adjs(k+1,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 + wst_adjs(k+1,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,ne,ex}, S^{\eta}_{i+1/2,j+1/2,k}S_{\eta,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i,j)+svetax(i,j+1)+svetax(i+1,j)+svetax(i+1,j+1))/4.
		a2=(svetay(i,j)+svetay(i,j+1)+svetay(i+1,j)+svetay(i+1,j+1))/4.
		b1=(svretax(i,j)+svretax(i,j+1)+svretax(i+1,j)+svretax(i+1,j+1))/4.
		b2=(svretay(i,j)+svretay(i,j+1)+svretay(i+1,j)+svretay(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(suxix(i+1,j-1)+suxix(i+1,j))/2. 		
		a2=(suxiy(i+1,j-1)+suxiy(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i,j)+vp(i+1,j-1)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		f1=f1 - vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=f2 - vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i+1,j))/2. 		
		a2=(svetay(i,j)+svetay(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i,j)+vp(i+1,j-1)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(suxix(i+1,j)+suxix(i+1,j+1))/2. 		
		a2=(suxiy(i+1,j)+suxiy(i+1,j+1))/2. 	
		a3=(vp(i,j)+vp(i,j+1)+vp(i+1,j)+vp(i+1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta
		beta=(vp(i+1,j)+vp(i+1,j+1))/2.
		f1=f1 - vst_adjs(k,ii+1,jj+1)*svetax(i+1,j+1)/beta
		f2=f2 - vst_adjs(k,ii+1,jj+1)*svetay(i+1,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i,j+1)+svetax(i+1,j+1))/2. 		
		a2=(svetay(i,j+1)+svetay(i+1,j+1))/2. 	
		a3=(vp(i,j)+vp(i,j+1)+vp(i+1,j)+vp(i+1,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta
		beta=(vp(i,j+1)+vp(i+1,j+1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj+1)*suxix(i+1,j+1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj+1)*suxiy(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,ne,ex}, S^{\xi}_{i+1,j,k}S_{\eta,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=suxix(i+1,j)
		a2=suxiy(i+1,j)
		b1=(svretax(i,j)+svretax(i,j+1)+svretax(i+1,j)+svretax(i+1,j+1))/4.
		b2=(svretay(i,j)+svretay(i,j+1)+svretay(i+1,j)+svretay(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i,j+1))/2. 		
		a2=(svetay(i,j)+svetay(i,j+1))/2.
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.	
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i+1,j))/2. 		
		a2=(suxiy(i,j)+suxiy(i+1,j))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.	
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i+1,j)+svetax(i+1,j+1))/2. 		
		a2=(svetay(i+1,j)+svetay(i+1,j+1))/2.
		a3=vp(i+1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		f1=vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta
		beta=(vp(i+1,j)+vp(i+1,j+1))/2.
		f1=f1 - vst_adjs(k,ii+1,jj+1)*svetax(i+1,j+1)/beta
		f2=f2 - vst_adjs(k,ii+1,jj+1)*svetay(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i+1,j)+suxix(i+2,j))/2. 		
		a2=(suxiy(i+1,j)+suxiy(i+2,j))/2. 
		a3=vp(i+1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta
		beta=(vp(i+1,j)+vp(i+2,j))/2.	
		f1=f1 - ust_adjs(k,ii+2,jj)*suxix(i+2,j)/beta
		f2=f2 - ust_adjs(k,ii+2,jj)*suxiy(i+2,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,ne,ex}, S^{z}_{i+1/2,j,k+1/2}S_{\eta,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i,j)+swz(i+1,j))/2.
		b1=(svretax(i,j)+svretax(i,j+1)+svretax(i+1,j)+svretax(i+1,j+1))/4.
		b2=(svretay(i,j)+svretay(i,j+1)+svretay(i+1,j)+svretay(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j)+swz(i+1,j))/2.
		a3=(vp(i,j)+vp(i+1,j))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k-1,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k-1,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 + ust_adjs(k+1,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 + ust_adjs(k+1,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\eta,i+1/2,j+1/2,k}S^{z}_{i+1/2,j,k+1/2}------------
		a1=(svretax(i,j)+svretax(i,j+1)+svretax(i+1,j)+svretax(i+1,j+1))/4.
		a2=(svretay(i,j)+svretay(i,j+1)+svretay(i+1,j)+svretay(i+1,j+1))/4.
		b3=(swz(i,j)+swz(i+1,j))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=suxix(i+1,j)
		a2=suxiy(i+1,j)
		a3=(vp(i,j)+vp(i+1,j))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i+1,j)
		f3=f3 - wst_adjs(k,ii+1,jj)*swz(i+1,j)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k+1,ii,jj)*swz(i,j)/beta
		beta=vp(i+1,j)
		f3=f3 + wst_adjs(k+1,ii+1,jj)*swz(i+1,j)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,se,ex}, S^{\eta}_{i+1/2,j-1/2,k}S_{\eta,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i,j-1)+svetax(i,j)+svetax(i+1,j-1)+svetax(i+1,j))/4.
		a2=(svetay(i,j-1)+svetay(i,j)+svetay(i+1,j-1)+svetay(i+1,j))/4.
		b1=(svretax(i,j-1)+svretax(i,j)+svretax(i+1,j-1)+svretax(i+1,j))/4.
		b2=(svretay(i,j-1)+svretay(i,j)+svretay(i+1,j-1)+svretay(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(suxix(i+1,j-2)+suxix(i+1,j-1))/2. 		
		a2=(suxiy(i+1,j-2)+suxiy(i+1,j-1))/2. 	
		a3=(vp(i,j-2)+vp(i,j-1)+vp(i+1,j-2)+vp(i+1,j-1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j-2)+vp(i,j-1))/2.
		f1=vst_adjs(k,ii,jj-1)*svetax(i,j-1)/beta
		f2=vst_adjs(k,ii,jj-1)*svetay(i,j-1)/beta
		beta=(vp(i+1,j-2)+vp(i+1,j-1))/2.
		f1=f1 - vst_adjs(k,ii+1,jj-1)*svetax(i+1,j-1)/beta
		f2=f2 - vst_adjs(k,ii+1,jj-1)*svetay(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i,j-1)+svetax(i+1,j-1))/2. 		
		a2=(svetay(i,j-1)+svetay(i+1,j-1))/2. 	
		a3=(vp(i,j-2)+vp(i,j-1)+vp(i+1,j-2)+vp(i+1,j-1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-2)+vp(i+1,j-2))/2.
		f1=ust_adjs(k,ii+1,jj-2)*suxix(i+1,j-2)/beta
		f2=ust_adjs(k,ii+1,jj-2)*suxiy(i+1,j-2)/beta
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(suxix(i+1,j-1)+suxix(i+1,j))/2. 		
		a2=(suxiy(i+1,j-1)+suxiy(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i,j)+vp(i+1,j-1)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		f1=f1 - vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=f2 - vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i+1,j))/2. 		
		a2=(svetay(i,j)+svetay(i+1,j))/2. 	
		a3=(vp(i,j-1)+vp(i,j)+vp(i+1,j-1)+vp(i+1,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,se,ex}, S^{\xi}_{i+1,j-1,k}S_{\eta,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=suxix(i+1,j-1)
		a2=suxiy(i+1,j-1)
		b1=(svretax(i,j-1)+svretax(i,j)+svretax(i+1,j-1)+svretax(i+1,j))/4.
		b2=(svretay(i,j-1)+svretay(i,j)+svretay(i+1,j-1)+svretay(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(svetax(i,j-1)+svetax(i,j))/2. 		
		a2=(svetay(i,j-1)+svetay(i,j))/2.
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-2)+vp(i,j-1))/2.
		f1=vst_adjs(k,ii,jj-1)*svetax(i,j-1)/beta
		f2=vst_adjs(k,ii,jj-1)*svetay(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i+1,j-1))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i+1,j-1))/2. 
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i+1,j-1)+svetax(i+1,j))/2. 		
		a2=(svetay(i+1,j-1)+svetay(i+1,j))/2.
		a3=vp(i+1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i+1,j-2)+vp(i+1,j-1))/2.
		f1=vst_adjs(k,ii+1,jj-1)*svetax(i+1,j-1)/beta
		f2=vst_adjs(k,ii+1,jj-1)*svetay(i+1,j-1)/beta
		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		f1=f1 - vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=f2 - vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i+1,j-1)+suxix(i+2,j-1))/2. 		
		a2=(suxiy(i+1,j-1)+suxiy(i+2,j-1))/2. 
		a3=vp(i+1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		beta=(vp(i+1,j-1)+vp(i+2,j-1))/2.
		f1=f1 - ust_adjs(k,ii+2,jj-1)*suxix(i+2,j-1)/beta
		f2=f2 - ust_adjs(k,ii+2,jj-1)*suxiy(i+2,j-1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,se,ex}, S^{z}_{i-1/2,j-1,k+1/2}S_{\eta,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i,j-1)+swz(i+1,j-1))/2.
		b1=(svretax(i,j-1)+svretax(i,j)+svretax(i+1,j-1)+svretax(i+1,j))/4.
		b2=(svretay(i,j-1)+svretay(i,j)+svretay(i+1,j-1)+svretay(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j-1)+swz(i+1,j-1))/2.
		a3=(vp(i,j-1)+vp(i+1,j-1))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k-1,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k-1,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		f1=f1 + ust_adjs(k+1,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 + ust_adjs(k+1,ii+1,jj-1)*suxiy(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\eta,i+1/2,j+1/2,k}S^{z}_{i-1/2,j,k+1/2}------------
		a1=(svretax(i,j-1)+svretax(i,j)+svretax(i+1,j-1)+svretax(i+1,j))/4.
		a2=(svretay(i,j-1)+svretay(i,j)+svretay(i+1,j-1)+svretay(i+1,j))/4.
		b3=(swz(i,j-1)+swz(i+1,j-1))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=suxix(i+1,j-1)
		a2=suxiy(i+1,j-1)
		a3=(vp(i,j-1)+vp(i+1,j-1))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i+1,j-1)
		f3=f3 - wst_adjs(k,ii+1,jj-1)*swz(i+1,j-1)/beta
		beta=vp(i,j-1)
		f3=f3 - wst_adjs(k+1,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i+1,j-1)
		f3=f3 + wst_adjs(k+1,ii+1,jj-1)*swz(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,nw,ex}, S^{\eta}_{i-1/2,j+1/2,k}S_{\eta,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i-1,j)+svetax(i-1,j+1)+svetax(i,j)+svetax(i,j+1))/4.
		a2=(svetay(i-1,j)+svetay(i-1,j+1)+svetay(i,j)+svetay(i,j+1))/4.
		b1=(svretax(i-1,j)+svretax(i-1,j+1)+svretax(i,j)+svretax(i,j+1))/4.
		b2=(svretay(i-1,j)+svretay(i-1,j+1)+svretay(i,j)+svretay(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(suxix(i,j-1)+suxix(i,j))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i-1,j)+vp(i,j-1)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i,j))/2. 		
		a2=(svetay(i-1,j)+svetay(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i-1,j)+vp(i,j-1)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(suxix(i,j)+suxix(i,j+1))/2. 		
		a2=(suxiy(i,j)+suxiy(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i-1,j+1)+vp(i,j)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j+1)+svetax(i,j+1))/2. 		
		a2=(svetay(i-1,j+1)+svetay(i,j+1))/2. 	
		a3=(vp(i-1,j)+vp(i-1,j+1)+vp(i,j)+vp(i,j+1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.		
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i-1,j+1)+vp(i,j+1))/2.	
		f1=f1 - ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=f2 - ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,nw,ex}, S^{\xi}_{i,j,k}S_{\eta,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=suxix(i,j)
		a2=suxiy(i,j)
		b1=(svretax(i-1,j)+svretax(i-1,j+1)+svretax(i,j)+svretax(i,j+1))/4.
		b2=(svretay(i-1,j)+svretay(i-1,j+1)+svretay(i,j)+svretay(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i-1,j+1))/2. 		
		a2=(svetay(i-1,j)+svetay(i-1,j+1))/2.
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j)+suxix(i,j))/2. 		
		a2=(suxiy(i-1,j)+suxiy(i,j))/2. 
		a3=vp(i-1,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		if (i-2 < 0) then
			beta=(vp(ngx-2,j)+vp(i-1,j))/2.
		else
			beta=(vp(i-2,j)+vp(i-1,j))/2.
		end if
		f1=ust_adjs(k,ii-1,jj)*suxix(i-1,j)/beta
		f2=ust_adjs(k,ii-1,jj)*suxiy(i-1,j)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j)+svetax(i,j+1))/2. 		
		a2=(svetay(i,j)+svetay(i,j+1))/2.
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k,ii,jj)*svetay(i,j)/beta
		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j)+suxix(i+1,j))/2. 		
		a2=(suxiy(i,j)+suxiy(i+1,j))/2. 
		a3=vp(i,j)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,nw,ex}, S^{z}_{i-1/2,j,k+1/2}S_{\eta,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i-1,j)+swz(i,j))/2.
		b1=(svretax(i-1,j)+svretax(i-1,j+1)+svretax(i,j)+svretax(i,j+1))/4.
		b2=(svretay(i-1,j)+svretay(i-1,j+1)+svretay(i,j)+svretay(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j)+swz(i,j))/2.
		a3=(vp(i-1,j)+vp(i,j))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k-1,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k-1,ii,jj)*suxiy(i,j)/beta
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta
		f1=f1 + ust_adjs(k+1,ii,jj)*suxix(i,j)/beta
		f2=f2 + ust_adjs(k+1,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\eta,i-1/2,j+1/2,k}S^{z}_{i-1/2,j,k+1/2}------------
		a1=(svretax(i-1,j)+svretax(i-1,j+1)+svretax(i,j)+svretax(i,j+1))/4.
		a2=(svretay(i-1,j)+svretay(i-1,j+1)+svretay(i,j)+svretay(i,j+1))/4.
		b3=(swz(i-1,j)+swz(i,j))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=suxix(i,j)
		a2=suxiy(i,j)
		a3=(vp(i-1,j)+vp(i,j))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i-1,j)
		f3=f3 - wst_adjs(k+1,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 + wst_adjs(k+1,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + (f3*e3)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,sw,ex}, S^{\eta}_{i-1/2,j-1/2,k}S_{\eta,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i-1,j-1)+svetax(i-1,j)+svetax(i,j-1)+svetax(i,j))/4.
		a2=(svetay(i-1,j-1)+svetay(i-1,j)+svetay(i,j-1)+svetay(i,j))/4.
		b1=(svretax(i-1,j-1)+svretax(i-1,j)+svretax(i,j-1)+svretax(i,j))/4.
		b2=(svretay(i-1,j-1)+svretay(i-1,j)+svretay(i,j-1)+svretay(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------		
		a1=(suxix(i,j-2)+suxix(i,j-1))/2. 		
		a2=(suxiy(i,j-2)+suxiy(i,j-1))/2. 	
		a3=(vp(i-1,j-2)+vp(i-1,j-1)+vp(i,j-2)+vp(i,j-1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j-2)+vp(i-1,j-1))/2.
		f1=vst_adjs(k,ii-1,jj-1)*svetax(i-1,j-1)/beta
		f2=vst_adjs(k,ii-1,jj-1)*svetay(i-1,j-1)/beta
		beta=(vp(i,j-2)+vp(i,j-1))/2.
		f1=f1 - vst_adjs(k,ii,jj-1)*svetax(i,j-1)/beta
		f2=f2 - vst_adjs(k,ii,jj-1)*svetay(i,j-1)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j-1)+svetax(i,j-1))/2. 		
		a2=(svetay(i-1,j-1)+svetay(i,j-1))/2. 	
		a3=(vp(i-1,j-2)+vp(i-1,j-1)+vp(i,j-2)+vp(i,j-1))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-2)+vp(i,j-2))/2.
		f1=ust_adjs(k,ii,jj-2)*suxix(i,j-2)/beta
		f2=ust_adjs(k,ii,jj-2)*suxiy(i,j-2)/beta
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------		
		a1=(suxix(i,j-1)+suxix(i,j))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i-1,j)+vp(i,j-1)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3			

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j)+svetax(i,j))/2. 		
		a2=(svetay(i-1,j)+svetay(i,j))/2. 	
		a3=(vp(i-1,j-1)+vp(i-1,j)+vp(i,j-1)+vp(i,j))/4.
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,sw,ex}, S^{\xi}_{i,j-1,k}S_{\eta,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=suxix(i,j-1)
		a2=suxiy(i,j-1)
		b1=(svretax(i-1,j-1)+svretax(i-1,j)+svretax(i,j-1)+svretax(i,j))/4.
		b2=(svretay(i-1,j-1)+svretay(i-1,j)+svretay(i,j-1)+svretay(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2
		c4=0.25*a2*b1
		c5=0.25*a2*b2

		!-----------------------------------------------------------------------
		a1=(svetax(i-1,j-1)+svetax(i-1,j))/2. 		
		a2=(svetay(i-1,j-1)+svetay(i-1,j))/2.
		a3=vp(i-1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i-1,j-2)+vp(i-1,j-1))/2.
		f1=vst_adjs(k,ii-1,jj-1)*svetax(i-1,j-1)/beta
		f2=vst_adjs(k,ii-1,jj-1)*svetay(i-1,j-1)/beta
		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i-1,j-1)+suxix(i,j-1))/2. 		
		a2=(suxiy(i-1,j-1)+suxiy(i,j-1))/2. 
		a3=vp(i-1,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		if (i-2<0) then
			beta=(vp(ngx-2,j-1)+vp(i-1,j-1))/2.
		else
			beta=(vp(i-2,j-1)+vp(i-1,j-1))/2.
		end if
		f1=ust_adjs(k,ii-1,jj-1)*suxix(i-1,j-1)/beta
		f2=ust_adjs(k,ii-1,jj-1)*suxiy(i-1,j-1)/beta
		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta	

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(svetax(i,j-1)+svetax(i,j))/2. 		
		a2=(svetay(i,j-1)+svetay(i,j))/2.
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3

		beta=(vp(i,j-2)+vp(i,j-1))/2.
		f1=vst_adjs(k,ii,jj-1)*svetax(i,j-1)/beta
		f2=vst_adjs(k,ii,jj-1)*svetay(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(suxix(i,j-1)+suxix(i+1,j-1))/2. 		
		a2=(suxiy(i,j-1)+suxiy(i+1,j-1))/2. 
		a3=vp(i,j-1)
		e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
		e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3	

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta	

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{\eta,sw,ex}, S^{z}_{i-1/2,j-1,k+1/2}S_{\eta,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(swz(i-1,j-1)+swz(i,j-1))/2.
		b1=(svretax(i-1,j-1)+svretax(i-1,j)+svretax(i,j-1)+svretax(i,j))/4.
		b2=(svretay(i-1,j-1)+svretay(i-1,j)+svretay(i,j-1)+svretay(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j-1)+swz(i,j-1))/2.
		a3=(vp(i-1,j-1)+vp(i,j-1))/2.
		e1=a1*c1*visc/a3  
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k-1,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k-1,ii,jj-1)*suxiy(i,j-1)/beta
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta
		f1=f1 + ust_adjs(k+1,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 + ust_adjs(k+1,ii,jj-1)*suxiy(i,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------S_{\eta,i-1/2,j-1/2,k}S^{z}_{i-1/2,j-1,k+1/2}------------
		a1=(svretax(i-1,j-1)+svretax(i-1,j)+svretax(i,j-1)+svretax(i,j))/4.
		a2=(svretay(i-1,j-1)+svretay(i-1,j)+svretay(i,j-1)+svretay(i,j))/4.
		b3=(swz(i-1,j-1)+swz(i,j-1))/2.
		c1=0.25*a1*b3
		c2=0.25*a2*b3

		a1=suxix(i,j-1)
		a2=suxiy(i,j-1)
		a3=(vp(i-1,j-1)+vp(i,j-1))/2.
		e3=(c1*a1+c2*a2)*visc/a3		

		beta=vp(i-1,j-1)
		f3=wst_adjs(k,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i,j-1)
		f3=f3 - wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i-1,j-1)
		f3=f3 - wst_adjs(k+1,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i,j-1)
		f3=f3 + wst_adjs(k+1,ii,jj-1)*swz(i,j-1)/beta

		diffex_adj = diffex_adj + (f3*e3)

		v_adj(k,ii,jj) = v_adj(k,ii,jj) + dt*diffex_adj

	end do
	end do
	end do

end 

subroutine advancez_Adjoint_diffex()
	use advance_module

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adjs(:, :, 0) = ust_adjs(:, :, 1)
       		vst_adjs(:, :, 1) =  0.0			
       		wst_adjs(:, :, 0) = wst_adjs(:, :, 1)		

       	      	!Extend
		vst_adjs(:, :, 0) = vst_adjs(:, :, 2)		! this halo is not included in following codes
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adjs(:, :, nlyb  ) = ust_adjs(:, :, nlyb-1)
              	vst_adjs(:, :, nlyb  ) =  0.0			
             	wst_adjs(:, :, nlyb  ) = wst_adjs(:, :, nlyb-1)	
		
		!Extend
		vst_adjs(:, :, nlyb+1) = vst_adjs(:, :, nlyb-1)	! this halo is not included in following codes
     		
	end if
 	
	!ccccccccccccccccccccc
	!   for w_adj equation
	!ccccccccccccccccccccc
	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! explicit diffusion term  D^{\dag}_{\eta,1,ex}
		!-----------------------------------------------------------------------
		! D^{\dag}_{z,P,ex}, S^{\xi}_{i+1/2,j,k-1/2}S_{z,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=swrz(i,j)
		b1=(suxix(i,j)+suxix(i+1,j))/2.
		b2=(suxiy(i,j)+suxiy(i+1,j))/2.
		c1=a1*b1
		c2=a1*b2

		a1=(swz(i-1,j)+swz(i,j))/2.
		a3=(vp(i-1,j)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k-1,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k-1,ii,jj)*suxiy(i,j)/beta
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = f1*e1 + f2*e2
		!-----------------------------------------------------------------------
		a1=(swz(i,j)+swz(i+1,j))/2.
		a3=(vp(i,j)+vp(i+1,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k-1,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k-1,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)		
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,P,ex}, S^{\eta}_{i,j+1/2,k-1/2}S_{z,i,j,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=swrz(i,j)
		b1=(svetax(i,j)+svetax(i,j+1))/2.
		b2=(svetay(i,j)+svetay(i,j+1))/2.
		c1=a1*b1
		c2=a1*b2

		a1=(swz(i,j-1)+swz(i,j))/2.
		a3=(vp(i,j-1)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k-1,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k-1,ii,jj)*svetay(i,j)/beta
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)	
		!-----------------------------------------------------------------------
		a1=(swz(i,j)+swz(i,j+1))/2.
		a3=(vp(i,j)+vp(i,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k-1,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k-1,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)	
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,ne,ex}, S^{\xi}_{i+1,j+1/2,k-1/2}S_{z,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i+1,j)+suxix(i+1,j+1))/2.
		a2=(suxiy(i+1,j)+suxiy(i+1,j+1))/2.
		b1=(swrz(i,j)+swrz(i,j+1)+swrz(i+1,j)+swrz(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=svetax(i,j+1)
		a2=svetay(i,j+1)
		a3=(vp(i,j)+vp(i,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i,j+1)
		f3=f3 - wst_adjs(k,ii,jj+1)*swz(i,j+1)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=svetax(i+1,j+1)
		a2=svetay(i+1,j+1)
		a3=(vp(i+1,j)+vp(i+1,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i+1,j)
		f3=wst_adjs(k,ii+1,jj)*swz(i+1,j)/beta
		beta=vp(i+1,j+1)
		f3=f3 - wst_adjs(k,ii+1,jj+1)*swz(i+1,j+1)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i+1/2,j+1/2,k}S^{\xi}_{i+1,j+1/2,k-1/2}-----------------
		a1=(swrz(i,j)+swrz(i,j+1)+swrz(i+1,j)+swrz(i+1,j+1))/4.
		b1=(suxix(i+1,j)+suxix(i+1,j+1))/2.
		b2=(suxiy(i+1,j)+suxiy(i+1,j+1))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j)+swz(i,j+1))/2.
		a3=(vp(i,j)+vp(i,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k-1,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k-1,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i+1,j)+swz(i+1,j+1))/2.
		a3=(vp(i+1,j)+vp(i+1,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i+1,j)+vp(i+1,j+1))/2.
		f1=vst_adjs(k-1,ii+1,jj+1)*svetax(i+1,j+1)/beta
		f2=vst_adjs(k-1,ii+1,jj+1)*svetay(i+1,j+1)/beta
		f1=f1 - vst_adjs(k,ii+1,jj+1)*svetax(i+1,j+1)/beta
		f2=f2 - vst_adjs(k,ii+1,jj+1)*svetay(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,ne,ex}, S^{\eta}_{i+1/2,j+1,k-1/2}S_{z,i+1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i,j+1)+svetax(i+1,j+1))/2.
		a2=(svetay(i,j+1)+svetay(i+1,j+1))/2.
		b1=(swrz(i,j)+swrz(i,j+1)+swrz(i+1,j)+swrz(i+1,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=suxix(i+1,j)
		a2=suxiy(i+1,j)
		a3=(vp(i,j)+vp(i+1,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i+1,j)
		f3=f3 - wst_adjs(k,ii+1,jj)*swz(i+1,j)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=suxix(i+1,j+1)
		a2=suxiy(i+1,j+1)
		a3=(vp(i,j+1)+vp(i+1,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j+1)
		f3=wst_adjs(k,ii,jj+1)*swz(i,j+1)/beta
		beta=vp(i+1,j+1)
		f3=f3 - wst_adjs(k,ii+1,jj+1)*swz(i+1,j+1)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i+1/2,j+1/2,k}S^{\eta}_{i+1/2,j+1,k-1/2}----------------
		a1=(swrz(i,j)+swrz(i,j+1)+swrz(i+1,j)+swrz(i+1,j+1))/4.
		b1=(svetax(i,j+1)+svetax(i+1,j+1))/2.
		b2=(svetay(i,j+1)+svetay(i+1,j+1))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j)+swz(i+1,j))/2.
		a3=(vp(i,j)+vp(i+1,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k-1,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k-1,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i,j+1)+swz(i+1,j+1))/2.
		a3=(vp(i,j+1)+vp(i+1,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j+1)+vp(i+1,j+1))/2.		
		f1=ust_adjs(k-1,ii+1,jj+1)*suxix(i+1,j+1)/beta
		f2=ust_adjs(k-1,ii+1,jj+1)*suxiy(i+1,j+1)/beta
		f1=f1 - ust_adjs(k,ii+1,jj+1)*suxix(i+1,j+1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj+1)*suxiy(i+1,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,nw,ex}, S^{\xi}_{i,j+1/2,k-1/2}S_{z,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i,j)+suxix(i,j+1))/2.
		a2=(suxiy(i,j)+suxiy(i,j+1))/2.
		b1=(swrz(i-1,j)+swrz(i-1,j+1)+swrz(i,j)+swrz(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=svetax(i-1,j+1)
		a2=svetay(i-1,j+1)
		a3=(vp(i-1,j)+vp(i-1,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i-1,j+1)
		f3=f3 - wst_adjs(k,ii-1,jj+1)*swz(i-1,j+1)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=svetax(i,j+1)
		a2=svetay(i,j+1)
		a3=(vp(i,j)+vp(i,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i,j+1)
		f3=f3 - wst_adjs(k,ii,jj+1)*swz(i,j+1)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i-1/2,j+1/2,k}S^{\xi}_{i,j+1/2,k-1/2}-----------------
		a1=(swrz(i-1,j)+swrz(i-1,j+1)+swrz(i,j)+swrz(i,j+1))/4.
		b1=(suxix(i,j)+suxix(i,j+1))/2.
		b2=(suxiy(i,j)+suxiy(i,j+1))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j)+swz(i-1,j+1))/2.
		a3=(vp(i-1,j)+vp(i-1,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i-1,j+1))/2.
		f1=vst_adjs(k-1,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=vst_adjs(k-1,ii-1,jj+1)*svetay(i-1,j+1)/beta
		f1=f1 - vst_adjs(k,ii-1,jj+1)*svetax(i-1,j+1)/beta
		f2=f2 - vst_adjs(k,ii-1,jj+1)*svetay(i-1,j+1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i,j)+swz(i,j+1))/2.
		a3=(vp(i,j)+vp(i,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i,j+1))/2.
		f1=vst_adjs(k-1,ii,jj+1)*svetax(i,j+1)/beta
		f2=vst_adjs(k-1,ii,jj+1)*svetay(i,j+1)/beta
		f1=f1 - vst_adjs(k,ii,jj+1)*svetax(i,j+1)/beta
		f2=f2 - vst_adjs(k,ii,jj+1)*svetay(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,nw,ex}, S^{\eta}_{i-1/2,j+1,k-1/2}S_{z,i-1/2,j+1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i-1,j+1)+svetax(i,j+1))/2.
		a2=(svetay(i-1,j+1)+svetay(i,j+1))/2.
		b1=(swrz(i-1,j)+swrz(i-1,j+1)+swrz(i,j)+swrz(i,j+1))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=suxix(i,j)
		a2=suxiy(i,j)
		a3=(vp(i-1,j)+vp(i,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=suxix(i,j+1)
		a2=suxiy(i,j+1)
		a3=(vp(i-1,j+1)+vp(i,j+1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j+1)
		f3=wst_adjs(k,ii-1,jj+1)*swz(i-1,j+1)/beta
		beta=vp(i,j+1)
		f3=f3 - wst_adjs(k,ii,jj+1)*swz(i,j+1)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i-1/2,j+1/2,k}S^{\eta}_{i-1/2,j+1,k-1/2}----------------
		a1=(swrz(i-1,j)+swrz(i-1,j+1)+swrz(i,j)+swrz(i,j+1))/4.
		b1=(svetax(i-1,j+1)+svetax(i,j+1))/2.
		b2=(svetay(i-1,j+1)+svetay(i,j+1))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j)+swz(i,j))/2.
		a3=(vp(i-1,j)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k-1,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k-1,ii,jj)*suxiy(i,j)/beta
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i-1,j+1)+swz(i,j+1))/2.
		a3=(vp(i-1,j+1)+vp(i,j+1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j+1)+vp(i,j+1))/2.
		f1=ust_adjs(k-1,ii,jj+1)*suxix(i,j+1)/beta
		f2=ust_adjs(k-1,ii,jj+1)*suxiy(i,j+1)/beta
		f1=f1 - ust_adjs(k,ii,jj+1)*suxix(i,j+1)/beta
		f2=f2 - ust_adjs(k,ii,jj+1)*suxiy(i,j+1)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,se,ex}, S^{\xi}_{i+1,j-1/2,k-1/2}S_{z,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i+1,j-1)+suxix(i+1,j))/2.
		a2=(suxiy(i+1,j-1)+suxiy(i+1,j))/2.
		b1=(swrz(i,j-1)+swrz(i,j)+swrz(i+1,j-1)+swrz(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=svetax(i,j)
		a2=svetay(i,j)
		a3=(vp(i,j-1)+vp(i,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=svetax(i+1,j)
		a2=svetay(i+1,j)
		a3=(vp(i+1,j-1)+vp(i+1,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i+1,j-1)
		f3=wst_adjs(k,ii+1,jj-1)*swz(i+1,j-1)/beta
		beta=vp(i+1,j)
		f3=f3 - wst_adjs(k,ii+1,jj)*swz(i+1,j)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i+1/2,j/2,k}S^{\xi}_{i+1,j/2,k-1/2}-----------------
		a1=(swrz(i,j-1)+swrz(i,j)+swrz(i+1,j-1)+swrz(i+1,j))/4.
		b1=(suxix(i+1,j-1)+suxix(i+1,j))/2.
		b2=(suxiy(i+1,j-1)+suxiy(i+1,j))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j-1)+swz(i,j))/2.
		a3=(vp(i,j-1)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k-1,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k-1,ii,jj)*svetay(i,j)/beta
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i+1,j-1)+swz(i+1,j))/2.
		a3=(vp(i+1,j-1)+vp(i+1,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i+1,j-1)+vp(i+1,j))/2.
		f1=vst_adjs(k-1,ii+1,jj)*svetax(i+1,j)/beta
		f2=vst_adjs(k-1,ii+1,jj)*svetay(i+1,j)/beta
		f1=f1 - vst_adjs(k,ii+1,jj)*svetax(i+1,j)/beta
		f2=f2 - vst_adjs(k,ii+1,jj)*svetay(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,se,ex}, S^{\eta}_{i+1/2,j,k-1/2}S_{z,i+1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i,j)+svetax(i+1,j))/2.
		a2=(svetay(i,j)+svetay(i+1,j))/2.
		b1=(swrz(i,j-1)+swrz(i,j)+swrz(i+1,j-1)+swrz(i+1,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=suxix(i+1,j-1)
		a2=suxiy(i+1,j-1)
		a3=(vp(i,j-1)+vp(i+1,j-1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i+1,j-1)
		f3=f3 - wst_adjs(k,ii+1,jj-1)*swz(i+1,j-1)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=suxix(i+1,j)
		a2=suxiy(i+1,j)
		a3=(vp(i,j)+vp(i+1,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j)
		f3=wst_adjs(k,ii,jj)*swz(i,j)/beta
		beta=vp(i+1,j)
		f3=f3 - wst_adjs(k,ii+1,jj)*swz(i+1,j)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i+1/2,j-1/2,k}S^{\eta}_{i+1/2,j,k-1/2}----------------
		a1=(swrz(i,j-1)+swrz(i,j)+swrz(i+1,j-1)+swrz(i+1,j))/4.
		b1=(svetax(i,j)+svetax(i+1,j))/2.
		b2=(svetay(i,j)+svetay(i+1,j))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i,j-1)+swz(i+1,j-1))/2.
		a3=(vp(i,j-1)+vp(i+1,j-1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i+1,j-1))/2.
		f1=ust_adjs(k-1,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=ust_adjs(k-1,ii+1,jj-1)*suxiy(i+1,j-1)/beta
		f1=f1 - ust_adjs(k,ii+1,jj-1)*suxix(i+1,j-1)/beta
		f2=f2 - ust_adjs(k,ii+1,jj-1)*suxiy(i+1,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i,j)+swz(i+1,j))/2.
		a3=(vp(i,j)+vp(i+1,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j)+vp(i+1,j))/2.
		f1=ust_adjs(k-1,ii+1,jj)*suxix(i+1,j)/beta
		f2=ust_adjs(k-1,ii+1,jj)*suxiy(i+1,j)/beta
		f1=f1 - ust_adjs(k,ii+1,jj)*suxix(i+1,j)/beta
		f2=f2 - ust_adjs(k,ii+1,jj)*suxiy(i+1,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,sw,ex}, S^{\xi}_{i,j-1/2,k-1/2}S_{z,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(suxix(i,j-1)+suxix(i,j))/2.
		a2=(suxiy(i,j-1)+suxiy(i,j))/2.
		b1=(swrz(i-1,j-1)+swrz(i-1,j)+swrz(i,j-1)+swrz(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=svetax(i-1,j)
		a2=svetay(i-1,j)
		a3=(vp(i-1,j-1)+vp(i-1,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j-1)
		f3=wst_adjs(k,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i-1,j)
		f3=f3 - wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=svetax(i,j)
		a2=svetay(i,j)
		a3=(vp(i,j-1)+vp(i,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i,j-1)
		f3=wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i-1/2,j-1/2,k}S^{\xi}_{i,j-1/2,k-1/2}-----------------
		a1=(swrz(i-1,j-1)+swrz(i-1,j)+swrz(i,j-1)+swrz(i,j))/4.
		b1=(suxix(i,j-1)+suxix(i,j))/2.
		b2=(suxiy(i,j-1)+suxiy(i,j))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j-1)+swz(i-1,j))/2.
		a3=(vp(i-1,j-1)+vp(i-1,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j-1)+vp(i-1,j))/2.
		f1=vst_adjs(k-1,ii-1,jj)*svetax(i-1,j)/beta
		f2=vst_adjs(k-1,ii-1,jj)*svetay(i-1,j)/beta
		f1=f1 - vst_adjs(k,ii-1,jj)*svetax(i-1,j)/beta
		f2=f2 - vst_adjs(k,ii-1,jj)*svetay(i-1,j)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i,j-1)+swz(i,j))/2.
		a3=(vp(i,j-1)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i,j-1)+vp(i,j))/2.
		f1=vst_adjs(k-1,ii,jj)*svetax(i,j)/beta
		f2=vst_adjs(k-1,ii,jj)*svetay(i,j)/beta
		f1=f1 - vst_adjs(k,ii,jj)*svetax(i,j)/beta
		f2=f2 - vst_adjs(k,ii,jj)*svetay(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! D^{\dag}_{z,sw,ex}, S^{\eta}_{i-1/2,j,k-1/2}S_{z,i-1/2,j-1/2,k}
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		a1=(svetax(i-1,j)+svetax(i,j))/2.
		a2=(svetay(i-1,j)+svetay(i,j))/2.
		b1=(swrz(i-1,j-1)+swrz(i-1,j)+swrz(i,j-1)+swrz(i,j))/4.
		c1=0.25*a1*b1
		c2=0.25*a2*b1

		a1=suxix(i,j-1)
		a2=suxiy(i,j-1)
		a3=(vp(i-1,j-1)+vp(i,j-1))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j-1)
		f3=wst_adjs(k,ii-1,jj-1)*swz(i-1,j-1)/beta
		beta=vp(i,j-1)
		f3=f3 - wst_adjs(k,ii,jj-1)*swz(i,j-1)/beta

		diffex_adj = diffex_adj + f3*e3
		!-----------------------------------------------------------------------
		a1=suxix(i,j)
		a2=suxiy(i,j)
		a3=(vp(i-1,j)+vp(i,j))/2.
		e3=(a1*c1+a2*c2)*visc/a3

		beta=vp(i-1,j)
		f3=wst_adjs(k,ii-1,jj)*swz(i-1,j)/beta
		beta=vp(i,j)
		f3=f3 - wst_adjs(k,ii,jj)*swz(i,j)/beta

		diffex_adj = diffex_adj - f3*e3
		!-----------S_{z,i-1/2,j-1/2,k}S^{\eta}_{i-1/2,j,k-1/2}----------------
		a1=(swrz(i-1,j-1)+swrz(i-1,j)+swrz(i,j-1)+swrz(i,j))/4.
		b1=(svetax(i-1,j)+svetax(i,j))/2.
		b2=(svetay(i-1,j)+svetay(i,j))/2.

		c1=0.25*a1*b1
		c2=0.25*a1*b2

		a1=(swz(i-1,j-1)+swz(i,j-1))/2.
		a3=(vp(i-1,j-1)+vp(i,j-1))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j-1)+vp(i,j-1))/2.
		f1=ust_adjs(k-1,ii,jj-1)*suxix(i,j-1)/beta
		f2=ust_adjs(k-1,ii,jj-1)*suxiy(i,j-1)/beta
		f1=f1 - ust_adjs(k,ii,jj-1)*suxix(i,j-1)/beta
		f2=f2 - ust_adjs(k,ii,jj-1)*suxiy(i,j-1)/beta

		diffex_adj = diffex_adj + (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------
		a1=(swz(i-1,j)+swz(i,j))/2.
		a3=(vp(i-1,j)+vp(i,j))/2.
		e1=a1*c1*visc/a3
		e2=a1*c2*visc/a3

		beta=(vp(i-1,j)+vp(i,j))/2.
		f1=ust_adjs(k-1,ii,jj)*suxix(i,j)/beta
		f2=ust_adjs(k-1,ii,jj)*suxiy(i,j)/beta
		f1=f1 - ust_adjs(k,ii,jj)*suxix(i,j)/beta
		f2=f2 - ust_adjs(k,ii,jj)*suxiy(i,j)/beta

		diffex_adj = diffex_adj - (f1*e1 + f2*e2)
		!-----------------------------------------------------------------------

		w_adj(k,ii,jj) = w_adj(k,ii,jj) + dt*diffex_adj

	end do
	end do
	end do


end

subroutine advancep_Adjoint()
	use advance_module
	integer :: flag_cost, my_jstart, my_jend, global_jndx, loc_jndx, jjout1

	!ccccccccccccccccccccccccccccccccccccccc
	!   calculate adjoint divergence of ust
	!ccccccccccccccccccccccccccccccccccccccc
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
			  	 ( - (ust_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (ust_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (ust_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  ust_adj(k,ii,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  ust_adj(k,ii,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  	&
			     -  ust_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 		&
			     -  ust_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (ust_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  ust_adj(k,ii+1,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  ust_adj(k,ii+1,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (ust_adj(k,ii+1  ,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  	&
			     -  ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 	&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
		! MZW: have considered that ust_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_en1   * &	
				( - (ust_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (ust_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (ust_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *	& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) ) *	&
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	 				&
			     - (ust_adj(k,ii+1,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *	& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) ) *	&
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )				
		! MZW: have considered that ust_adj(k,ii,jj) = 0 at wall 
		gradx_adj= gradx_adj + jflag_st1   * &	
				( - (ust_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (ust_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     + (ust_adj(k,ii  ,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *	& 
				(svetax(i-1,j+1  )+svetax(i,j+1  ))/2.*0.25				&
			     + (ust_adj(k,ii  ,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. ) ) *	& 
				(svetay(i-1,j+1  )+svetay(i,j+1  ))/2.*0.25				&
			     + (ust_adj(k,ii+1 ,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *	& 
				(svetax(i,j+1  )+svetax(i+1,j+1  ))/2.*0.25				&
			     + (ust_adj(k,ii+1 ,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. ) ) *	& 
				(svetay(i,j+1  )+svetay(i+1,j+1  ))/2.*0.25   )

		gradx_adj= gradx_adj + jflag_en2  * &
		   	    ( - (ust_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (ust_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (ust_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			      )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			      )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj-1)*suxix(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. )  	&
			     -  ust_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj-1)*suxiy(i,j-1)/( (vp(i,j-1)+vp(i-1,j-1))/2. ) 		&
			     -  ust_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (ust_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			      )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			      )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - (ust_adj(k,ii+1  ,jj-1)*suxix(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. )  	&
			     -  ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj-1)*suxiy(i+1,j-1)/( (vp(i+1,j-1)+vp(i,j-1))/2. ) 	&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )

		gradx_adj= gradx_adj + jflag_st2  * &
		   	    ( - (ust_adj(k,ii,jj)*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )  		&
			     -	ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxix(i,j)+suxix(i+1,j))/2.						&
			     - (ust_adj(k,ii,jj)*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. )      		&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. )   )*	& 
				(suxiy(i,j)+suxiy(i+1,j))/2.						&
			     - (ust_adj(k,ii  ,jj)*suxix(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  ust_adj(k,ii,jj+1)*suxix(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetax(i-1,j+1)+svetax(i,j+1))/2.*0.25					&
			     - (ust_adj(k,ii  ,jj)*suxiy(i  ,j  )/( (vp(i,j)+vp(i-1,j))/2. ) 		&
			     -  ust_adj(k,ii,jj+1)*suxiy(i,j+1)/( (vp(i,j+1)+vp(i-1,j+1))/2. )  )*	& 
				(svetay(i-1,j+1)+svetay(i,j+1))/2.*0.25					&
			     - (  	&
			     -  ust_adj(k,ii,jj  )*suxix(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetax(i-1,j  )+svetax(i,j  ))/2.*0.25					&
			     - (	&
			     -  ust_adj(k,ii,jj  )*suxiy(i,j  )/( (vp(i,j)+vp(i-1,j))/2. ))*		& 
				(svetay(i-1,j  )+svetay(i,j  ))/2.*0.25	    				&
			     - (ust_adj(k,ii+1,jj)*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  ust_adj(k,ii+1,jj+1)*suxix(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetax(i,j+1)+svetax(i+1,j+1))/2.*0.25					&
			     - (ust_adj(k,ii+1,jj)*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ) 		&
			     -  ust_adj(k,ii+1,jj+1)*suxiy(i+1,j+1)/( (vp(i+1,j+1)+vp(i,j+1))/2. )  )*	& 
				(svetay(i,j+1)+svetay(i+1,j+1))/2.*0.25					&
			     - ( 	&
			     -  ust_adj(k,ii+1,jj  )*suxix(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetax(i,j  )+svetax(i+1,j  ))/2.*0.25					&
			     - (	&
			     -  ust_adj(k,ii+1,jj  )*suxiy(i+1,j  )/( (vp(i+1,j)+vp(i,j))/2. ))*		& 
				(svetay(i,j  )+svetay(i+1,j  ))/2.*0.25	    )
				
		!-----------------------------------------------------------------------
		! Notice that p_adj is actually p_adj * dt
		p_adj(k,ii,jj) = p_adj(k,ii,jj) - dt * dt * gradx_adj

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
				( - (vst_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	vst_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (vst_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  vst_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (vst_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  vst_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (vst_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  vst_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (vst_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  vst_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (vst_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  vst_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (vst_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  vst_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (vst_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  vst_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (vst_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  vst_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (vst_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  vst_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )
		!MZW: periodic boundary condition has been considered here, i2_T+1 -> i1_T
		! at inlet and outlet, assume the grid is orthogonal
		grady_adj=    grady_adj + iflag_st1 * 	&
				( - (vst_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	vst_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (vst_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  vst_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (				&
			     -  vst_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (				&
			     -  vst_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25	    				&
			     - (		 		&
			     -  vst_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (		 		&
			     -  vst_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25		)

		grady_adj=    grady_adj + iflag_en1 * 	&
				( - (vst_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	vst_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (vst_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  vst_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (vst_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
				     )*				& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (vst_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  	&
			     	     )*				& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (vst_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     	     )*				& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (vst_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     	     )*				& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )

		grady_adj=   grady_adj + iflag_st2 *   &
				( - (vst_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	vst_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (vst_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  vst_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (vst_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  vst_adj(k,ii+1,jj)*svetax(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (vst_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     -  vst_adj(k,ii+1,jj)*svetay(i+1,j)/( (vp(i+1,j)+vp(i+1,j-1))/2. )  )*	& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (		  		&
			     -  vst_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (			  	&
			     -  vst_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (vst_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  vst_adj(k,ii+1,jj+1)*svetax(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (vst_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     -  vst_adj(k,ii+1,jj+1)*svetay(i+1,j+1)/( (vp(i+1,j+1)+vp(i+1,j))/2. )  )*	& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (			 	&
			     -  vst_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (			  	&
			     -  vst_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )	

		grady_adj=   grady_adj + iflag_en2 *   &
				( - (vst_adj(k,ii,jj)*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )  		&
			     -	vst_adj(k,ii,jj+1  )*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetax(i,j)+svetax(i,j+1))/2.						&
			     - (vst_adj(k,ii,jj)*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. )      		&
			     -  vst_adj(k,ii,jj+1  )*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. )   )*	& 
				(svetay(i,j)+svetay(i,j+1))/2.						&
			     - (vst_adj(k,ii  ,jj)*svetax(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     	  )*				& 
				(suxix(i+1,j-1)+suxix(i+1,j))/2.*0.25					&
			     - (vst_adj(k,ii  ,jj)*svetay(i  ,j  )/( (vp(i,j)+vp(i,j-1))/2. ) 		&
			     	  )*				& 
				(suxiy(i+1,j-1)+suxiy(i+1,j))/2.*0.25					&
			     - (vst_adj(k,ii-1,jj)*svetax(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  vst_adj(k,ii,jj  )*svetax(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxix(i,j-1  )+suxix(i,j  ))/2.*0.25					&
			     - (vst_adj(k,ii-1,jj)*svetay(i-1,j)/( (vp(i-1,j)+vp(i-1,j-1))/2. )  		&
			     -  vst_adj(k,ii,jj  )*svetay(i,j  )/( (vp(i,j)+vp(i,j-1))/2. ))*		& 
				(suxiy(i,j-1  )+suxiy(i,j  ))/2.*0.25	    				&
			     - (vst_adj(k,ii,jj+1)*svetax(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			    	  )*				& 
				(suxix(i+1,j)+suxix(i+1,j+1))/2.*0.25					&
			     - (vst_adj(k,ii,jj+1)*svetay(i,j+1  )/( (vp(i,j+1)+vp(i,j))/2. ) 		&
			     	  )*				& 
				(suxiy(i+1,j)+suxiy(i+1,j+1))/2.*0.25					&
			     - (vst_adj(k,ii-1  ,jj+1)*svetax(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  vst_adj(k,ii,jj+1)*svetax(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxix(i,j  )+suxix(i,j+1 ))/2.*0.25					&
			     - (vst_adj(k,ii-1  ,jj+1)*svetay(i-1,j+1)/( (vp(i-1,j+1)+vp(i-1,j))/2. )  	&
			     -  vst_adj(k,ii,jj+1)*svetay(i,j+1)/( (vp(i,j+1)+vp(i,j))/2. ))*		& 
				(suxiy(i,j  )+suxiy(i,j+1 ))/2.*0.25   )
			
		!-----------------------------------------------------------------------
		p_adj(k,ii,jj) = p_adj(k,ii,jj) - dt * dt * grady_adj

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
		gradz_adj= (- (wst_adj(k,ii,jj)*swz(i,j)/(vp(i,j))					&
			    -wst_adj(k+1,ii,jj)*swz(i,j)/(vp(i,j))	)*				&
			    swz(i,j)		)							
		!-----------------------------------------------------------------------
		p_adj(k,ii,jj) = p_adj(k,ii,jj) - dt * dt * gradz_adj

	end do
	end do
	end do
	!------------------------observation----------------------------------

        
	flag_cost = 0

	my_jstart = jTmin_1(jblock)
	my_jend   = jTmax_1(jblock)

	jjout1 = 1
	global_jndx = jjout1
	loc_jndx    = jmap_1(global_jndx)
	if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
		flag_cost = 1
	end if

	! Here p is p*dt !!!!!!!!!!!!!!!!!!!!	
	if (flag_cost == 1) then
		do ii=i1_T,i2_T
		do kk=1,ngzm
			k = kk
			p_adj(k,ii,loc_jndx) = p_adj(k,ii,loc_jndx) - (p_tmp(k,ii,1) - p_obs(k,ii,1))
		end do
		end do
	end if

	jjout1 = ngy-1
	global_jndx = jjout1
	loc_jndx    = jmap_1(global_jndx)
	if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
		flag_cost = 2
	end if
	
	if (flag_cost == 2) then	
		do ii=i1_T,i2_T
		do kk=1,ngzm
			k = kk
			p_adj(k,ii,loc_jndx) = p_adj(k,ii,loc_jndx) - (p_tmp(k,ii,2) - p_obs(k,ii,2))
		end do
		end do
	end if

	flag_cost = 0
	
end

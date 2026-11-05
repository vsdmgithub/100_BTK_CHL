! MZW 02/07/17
!! calcualte adjoint operator of A1, A2, A3 in u, v, w
!!-----------------------------------------
!! invertu_Adjoint()
!! invertv_Adjoint()
!! invertw_Adjoint()


subroutine invertu_Adjoint()
	use data_export
	real, dimension(:,:)  , allocatable :: aa, bb, cc, rrz, rrtemp
!	real, dimension(:,:)  , allocatable :: aa_g, bb_g, cc_g  ! for debug
	real, dimension(ngzm,i1_u:i2_u,j1_u:j2_u) :: rr 
	real :: dt2_inv
	integer :: iPrj

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MZW--calculate the rhs of adjoint equation (6.171)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	rhs = 0.0
	call readInt("iPrj", iPrj)     ! 1 for fixed mdot, 0 for fixed dpdx
!	------------ This part is put in Poisson solver------------
!	call readInt("iPrj", iPrj)
!	if (iPrj.eq.1) then
!		GPxi = 0.
!		call calculate_GPxi_Adjoint(GPxi_adj)
!		!if (irank.eq.1) then
!		!	write(6,*),'MASS FLOW RATE FIXED'
!		!end if
!	else
!		write(*,*) 'PRESSURE GRADIENT FIXED MODE HAS NOT BEEN DEVELOPED'
!		stop
!	end if

	dt2_inv = 1.0/dt**2
	do jj=j1_u,j2_u
	do ii=i1_u,i2_u
	do kk=1,ngzm
		! here GPxi_adj = F^{n+1}_x * dx / V, dx/V has been included in calculate_GPxi
		rhs(kk,ii,jj) =  u_adj(kk,ii,jj) + iPrj*GPxi_adj &
				+ (phi_adj(kk,ii-1,jj) - phi_adj(kk,ii,jj)) * dt2_inv
		!rhs(kk,ii,jj) = cos(real(ii+jj))      ! for debug
	end do
	end do
	end do

	rhsu(:,:,:) = rhs(:,:,:)
	! Check rhs
	!write(*,*) 'rhs of invertu_Adjoint at kk=1, jj=1'
	!write(*,*) rhs(1,i1_u:i2_u,1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SY--initialize the temporary rhs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr=0.0 
!	! set values for rhs, test tri-diagonal solver
!	! MZW ---------------
!	do jj=j1_u, j2_u
!		rhs(1:ngzm,i1_u:i2_u,jj) = jj !* 1.E-6
!	end do
!	! --------------------
	rr(:,:,:) = rhs(1:ngz-1, i1_u:i2_u, j1_u:j2_u)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     allocate(aa(i1_u:i2_u,1:ngzm),bb (i1_u:i2_u,1:ngzm)) 
      allocate(cc(i1_u:i2_u,1:ngzm),rrz(i1_u:i2_u,1:ngzm)) 
      allocate(rrtemp(1:ngzm,i1_u:i2_u)) 
      aa=0.0; bb=0.0; cc=0.0; rrz=0.0; rrtemp=0.0

      do 3000 jj=j1_u, j2_u
	rrtemp(1:ngzm, i1_u:i2_u) = rr(1:ngzm, i1_u:i2_u, jj)

      do ii=i1_u, i2_u
      do kk=1   , ngzm
        rrz(ii,kk) = rrtemp(kk,ii)
      end do
      end do
      do kk=1   , ngzm
      do ii=i1_u, i2_u
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)
        k    = kk

!---------------------------------------
!--------(eq 6.19) first term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{xi,P}+S_{xi,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{xi,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i-1,j))/2.
        tta=(vp(i,j)+vp(i-1,j))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
	beta = -dt/(2.*(vp(i,j)+vp(i-1,j))/2.)
        cc(ii,kk)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2) 

!---------------------------------------
!--------(eq 6.19) second term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{xi,P}+S_{xi,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{xi,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i-1,j))/2.
        tta=(vp(i,j)+vp(i-1,j))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{xi,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i-1,j))/2.
        tta=(vp(i,j)+vp(i-1,j))/2.
        e1 =e1+(d3*c7)*       visc/tta
        e2 =e2+(d3*c8)*       visc/tta

	beta = -dt/(2.*(vp(i,j)+vp(i-1,j))/2.)
        bb(ii,kk)=1.-beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

!---------------------------------------
!--------(eq 6.19) third term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{xi,P}+S_{xi,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{xi,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i-1,j))/2.
        tta=(vp(i,j)+vp(i-1,j))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
	beta = -dt/(2.*(vp(i,j)+vp(i-1,j))/2.)
        aa(ii,kk)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

      end do ! ii
      end do ! kk


!*****************************************************************************
      if(BC_front.eq.2) then
! MZW 02/14/17 check coefficients
	!if(jj.eq.2) then
	!write(*,*) '---------Start printing adjoint operator---------'
	!write(*,*) 'Start printing aa at (ii,jj)=', i1_u, jj
	!write(*,*) aa(i1_u,1:ngzm)
	!write(*,*) 'Start printing bb at (ii,jj)=', i1_u, jj
	!write(*,*) bb(i1_u,1:ngzm)
	!write(*,*) 'Start printing cc at (ii,jj)=', i1_u, jj
	!write(*,*) cc(i1_u,1:ngzm)
	!end if

      call triPeriodicM_(id_z, aa(i1_u,1), bb (i1_u,1),        &
                               cc(i1_u,1), rrz(i1_u,1), ngzm, niu)
	! MZW 02/14/17 check results
	!if (irank == iroot) then
	!	write(*,*) 'Writing TRI results of A_{xi,2} at (ii,jj)=',i1_u,jj
	!	write(*,*) rrz(i1_u,1:ngzm)
	!end if
      end if
      if(BC_front.ne.2) then
		write(*,*) "Not periodic boundary conditions!!!"
		stop
      end if
!*****************************************************************************

      do ii=i1_u, i2_u
      do kk=1   , ngzm
        rrtemp(kk,ii) = rrz(ii,kk)
      end do
      end do

      rr(1:ngzm, i1_u:i2_u, jj) = rrtemp(1:ngzm,i1_u:i2_u)
3000 continue
!*****************************************************************************
	deallocate(aa,bb,cc)
!*****************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! next solve for (1-A2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
allocate(aa(i1_u:i2_u,j1_u:j2_u),bb(i1_u:i2_u,j1_u:j2_u),cc(i1_u:i2_u,j1_u:j2_u))
aa=0.0; bb=0.0; cc=0.0
      do 83   jj=j1_u, j2_u
      do 83   ii=i1_u, i2_u
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!---------------------------------------
!-------( eq 6.18 )--first term---------
!---------------------------------------
!
!------S^{\eta}_P S_{\xi,P} + S_{\xi,P} S^{\eta}
        a1 =(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
        a2 =(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
        b1 = surxix(i,j)
        b2 = surxiy(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
 !-----S^{\eta}_n/(Re_n V_n)\bu
        d1 =(svetax(i,j+1)+svetax(i-1,j+1))/2.
        d2 =(svetay(i,j+1)+svetay(i-1,j+1))/2.
        tta=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.

        e1 =(d1*c1+d2*c4)*       visc/tta
        e2 =(d1*c2+d2*c5)*       visc/tta
!------dt/(2V_{\xi,N}) S^{\xi}_N \bu
	beta = -dt/(2.*(vp(i,j+1)+vp(i-1,j+1))/2.)
        cc(ii,jj)=beta*(suxix(i,j+1)*e1+suxiy(i,j+1)*e2) 


!---------------------------------------
!-------( eq 6.18 )--second term--------
!---------------------------------------
        a1 =(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
        a2 =(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
        b1 = surxix(i,j)
        b2 = surxiy(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
        
	d1 =(svetax(i,j+1)+svetax(i-1,j+1))/2.
        d2 =(svetay(i,j+1)+svetay(i-1,j+1))/2.
        tta=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
        e1 =(d1*c1+d2*c4)*       visc/tta
        e2 =(d1*c2+d2*c5)*       visc/tta
        
	d1 =(svetax(i,j)+svetax(i-1,j))/2.
        d2 =(svetay(i,j)+svetay(i-1,j))/2.
        tta=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.

        e1 =e1+(d1*c1+d2*c4)*       visc/tta
        e2 =e2+(d1*c2+d2*c5)*       visc/tta
	beta = -dt/(2.*(vp(i,j)+vp(i-1,j))/2.)
        bb(ii,jj)=1.-beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

!---------------------------------------
!-------( eq 6.18 )--third term--------
!---------------------------------------
        a1 =(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
        a2 =(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
        b1 = surxix(i,j)
        b2 = surxiy(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)

        d1 =(svetax(i,j)+svetax(i-1,j))/2.
        d2 =(svetay(i,j)+svetay(i-1,j))/2.
        tta=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
        e1 =(d1*c1+d2*c4)*       visc/tta
        e2 =(d1*c2+d2*c5)*       visc/tta
	beta = -dt/(2.*(vp(i,j-1)+vp(i-1,j-1))/2.)
        aa(ii,jj)=beta*(suxix(i,j-1)*e1+suxiy(i,j-1)*e2)

 83   continue

! Deal with B.Cs
      if (jblock.eq.1) then
         do 77 ii=i1_u,i2_u
           aa(ii,j1_u)=0.
 77      continue
      end if

      if (jblock.eq.npy) then
        do 91 ii=i1_u,i2_u
           cc(ii,j2_u)=0.
 91     continue
      end if
! Solve tri-diagnal system
!                              ORIGINAL
!*****************************************************************************
	! MZW 02/14/17 check coefficients
	!write(*,*) '---------Start printing adjoint operator---------'
	!ip=i2_u
	!write(*,*) 'Start printing aa at ii=', ip
	!write(*,*) aa(ip,j1_u:j2_u)
	!write(*,*) 'Start printing bb at ii=', ip
	!write(*,*) bb(ip,j1_u:j2_u)
	!write(*,*) 'Start printing cc at ii=', ip
	!write(*,*) cc(ip,j1_u:j2_u)
	! just for debug!!!!!
	!rr(:,:,:) = rhs(1:ngz-1, i1_u:i2_u, j1_u:j2_u)
	! just for debug!!!!
      call triDiagonalM_mrhs(id_y, aa(i1_u,j1_u), bb(  i1_u,j1_u),        &
                                   cc(i1_u,j1_u), rr(1,i1_u,j1_u), nju, niu, ngzm)

!temp = 0.0
!	do jj=j1_u,j2_u
!		temp = temp + rr(1,i1_u,jj) * rhs(1,i1_u,jj)
!	end do
!	!call globalSum(temp,1)
!	if (irank == 1) then
!	write(*,*) '----------adjoint operator---------'
!	write(*,*) temp
!	end if

	! MZW 02/14/17 check results
	!if (irank == iroot) then
	!	write(*,*) 'Writing TRI results of A_{xi,2}'
	!	write(*,*) rr(1,i1_u,:)
	!end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! next solve for (1-A1)(ust-u(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      aa=0.0; bb=0.0; cc=0.0 

      do 75   jj=j1_u, j2_u
      do 75   ii=i1_u, i2_u
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!--------------------------------------------
!------------(eq 6.17)--first term-----------
!--------------------------------------------
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
	beta = -dt/(2.*(vp(i+1,j)+vp(i,j))/2.)

	cc(ii,jj)=beta*(suxix(i+1,j)*e1+suxiy(i+1,j)*e2)

!--------------------------------------------
!-----------(equ 6.17)--second term----------
!--------------------------------------------	
	a1 = suxix (i,j)
	a2 = suxiy (i,j)
	b1 = surxix(i,j)
	b2 = surxiy(i,j)
	c1 = (a1*b1+b1*a1)
	c2 = (a1*b2+b1*a2)
	c4 = (a2*b1+b2*a1)
	c5 = (a2*b2+b2*a2)
	d1 = (suxix(i,j)+suxix(i+1,j))/2. * visc/vp(i,j) &
		& + (suxix(i,j)+suxix(i-1,j))/2. * visc/vp(i-1,j)
	d2 = (suxiy(i,j)+suxiy(i+1,j))/2. * visc/vp(i,j) &
		& + (suxiy(i,j)+suxiy(i-1,j))/2. * visc/vp(i-1,j)
	e1 = (d1*c1+d2*c4)
	e2 = (d1*c2+d2*c5)
	beta = -dt/(2.*(vp(i,j)+vp(i-1,j))/2.)

	bb(ii,jj)=1.-beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

!--------------------------------------------
!------------(eq 6.17)--third term-----------
!--------------------------------------------
	a1 = suxix (i,j)
	a2 = suxiy (i,j)
	b1 = surxix(i,j)
	b2 = surxiy(i,j)
	c1 = (a1*b1+b1*a1)
	c2 = (a1*b2+b1*a2)
	c4 = (a2*b1+b2*a1)
	c5 = (a2*b2+b2*a2)
	d1 = (suxix(i,j)+suxix(i-1,j))/2.
	d2 = (suxiy(i,j)+suxiy(i-1,j))/2. 
	e1 = (d1*c1+d2*c4)*   visc/vp(i-1,j)
	e2 = (d1*c2+d2*c5)*   visc/vp(i-1,j)
	beta = -dt/(2.*(vp(i-1,j)+vp(i-2,j))/2.)

	aa(ii,jj)=beta*(suxix(i-1,j)*e1+suxiy(i-1,j)*e2)

 75   continue


! Deal with boundary conditions
! Solve tri-diagonal system
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! QW: Solve periodic system
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (BC_left.eq.2) then

        !rr(:,:,:) = rhs(1:ngz-1, i1_u:i2_u, j1_u:j2_u)
	! MZW 02/14/17 check coefficients
	!write(*,*) '---------Start printing adjoint operator---------'
	!jp = j2_u-1
	!write(*,*) 'Start printing aa at jj=', jp
	!write(*,*) aa(i1_u:i2_u,jp)
	!write(*,*) 'Start printing bb at jj=', jp
	!write(*,*) bb(i1_u:i2_u,jp)
	!write(*,*) 'Start printing cc at jj=', jp
	!write(*,*) cc(i1_u:i2_u,jp)
        do jj=j1_u,j2_u
                call triPeriodic_mrhs(id_x, aa(i1_u,jj), bb(i1_u,jj),                &
                                            cc(i1_u,jj), rr(1,i1_u,jj), ngz-1, niu)
        enddo
	! MZW 02/14/17 check results
	!if (irank == iroot) then
	!	write(*,*) 'Writing TRI results of adjoint A_{xi,1}'
	!	write(*,*) rr(1,i1_u:i2_u,j1_u)
	!end if
else 
	write(*,*) 'Not periodic boundary condition!!!'
	stop
end if


! check adjoint relation
!temp = 0.0
!	do jj=j1_u,j2_u
!	do ii=i1_u,i2_u
!	do kk=1,ngzm
!		temp = temp + rr(kk,ii,jj) * rhs(kk,ii,jj)
!	end do
!	end do
!	end do
!	!call globalSum(temp,1)
!	if (irank == 1) then
!	write(*,*) '----------adjoint operator---------'
!	write(*,*) temp
!	end if

! update ust here  !!!!!!!!!!!!!!!!!
	ust_adj(1:ngzm,i1_u:i2_u,j1_u:j2_u) = rr(1:ngzm,i1_u:i2_u,j1_u:j2_u)
	! need speical treatment of ust_adj(1:ngzm,i1_u:i2_u,j1_u)

end


subroutine invertv_Adjoint()
	use data_export
	real, dimension(:,:), allocatable :: aa, bb, cc, rrz, rrtemp
	!real data_gather(5,5) ! just for debug
	!real, dimension(:,:)  , allocatable :: aa_g, bb_g, cc_g
	real, dimension(ngzm, i1_v:i2_v,j1_v:j2_v) :: rr
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MZW--calculate the rhs of adjoint equation (6.171) except GPxi_adj term
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	rhs = 0.0

	dt2_inv = 1.0/dt**2
	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
	do kk=1,ngzm
		! here GPxi_adj = F^{n+1}_x * dx / V, dx/V has been included in calculate_GPxi
		rhs(kk,ii,jj) =  v_adj(kk,ii,jj) + (phi_adj(kk,ii,jj-1) - phi_adj(kk,ii,jj)) * dt2_inv
		!rhs(kk,ii,jj) = cos(real(ii+jj+kk))  
	end do
	end do
	end do

	rhsv(:,:,:) = rhs(:,:,:)
	! Check rhs
	!write(*,*) 'rhs of invertv_Adjoint at kk=1, ii=1'
	!write(*,*) rhs(1,1,j1_v:j2_v)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MZW--initialize the temporary rhs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr=0.0
	! set values for rhs, test tri-diagonal solver
	! MZW ---------------
	!do jj=j1_v, j2_v
	!	rhs(1:ngzm,i1_v:i2_v,jj) = jj * 1.E-6
	!end do
	! --------------------
	rr(:,:,:) = rhs(1:ngz-1, i1_v:i2_v, j1_v:j2_v)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A3*)(1-A2*)(1-A1*)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A3*)
      allocate(aa(i1_v:i2_v,1:ngzm),bb (i1_v:i2_v,1:ngzm)) 
      allocate(cc(i1_v:i2_v,1:ngzm),rrz(i1_v:i2_v,1:ngzm)) 
      allocate(rrtemp(1:ngzm,i1_v:i2_v)) 
      aa=0.0; bb=0.0; cc=0.0; rrz=0.0; rrtemp=0.0

      do 3000 jj=j1_v, j2_v
	rrtemp(1:ngzm, i1_v:i2_v) = rr(1:ngzm, i1_v:i2_v, jj)

	do ii=i1_v, i2_v
	do kk=1   , ngzm
	        rrz(ii,kk) = rrtemp(kk,ii)
	end do
	end do

      do kk=1   , ngzm
      do ii=i1_v, i2_v
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)
        k    = kk

!---------------------------------------
!-------(eq 6.22 first term)-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{F}S_{eta,F}+S_{eta,F}S^{z}_{F})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{eta,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i,j-1))/2.
        tta=(vp(i,j)+vp(i,j-1))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
	beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)
        cc(ii,kk)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------(eq 6.22 second term)-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{eta,P}+S_{eta,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{eta,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i,j-1))/2.
        tta=(vp(i,j)+vp(i,j-1))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{eta,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i,j-1))/2.
        tta=(vp(i,j)+vp(i,j-1))/2.
        e1 =e1+(d3*c7)*       visc/tta
        e2 =e2+(d3*c8)*       visc/tta
	beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)
        bb(ii,kk)=1.-beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!--------(eq 6.22 third term)-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{F}S_{eta,F}+S_{eta,F}S^{z}_{F})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        c7=(a3*b1)
        c8=(a3*b2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{eta,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =(swz(i,j)+swz(i,j-1))/2.
        tta=(vp(i,j)+vp(i,j-1))/2.
        e1 =(d3*c7)*       visc/tta
        e2 =(d3*c8)*       visc/tta
	beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)
        aa(ii,kk)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

      end do ! ii
      end do ! kk
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(BC_front.ne.2) then
        do ii=i1_v,i2_v
          !bb(ii,1   )=bb(ii,1   )+aa(ii,1   )
          !bb(ii,ngzm)=bb(ii,ngzm)+cc(ii,ngzm)
          aa(ii,1   )=0.
          cc(ii,ngzm)=0.
        end do
      end if

!*****************************************************************************
      if(BC_front.eq.2) then
      call triPeriodicM_(id_z, aa(i1_v,1), bb (i1_v,1),        &
                               cc(i1_v,1), rrz(i1_v,1), ngzm, niv)
	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{eta,3} at (ii,jj)=',i1_v,jj
!		write(*,*) rrz(i1_v,1:ngzm)
!	end if
      end if
      if(BC_front.ne.2) then
      call triDiagonalM_(id_z, aa(i1_v,1), bb (i1_v,1),        &
                               cc(i1_v,1), rrz(i1_v,1), ngzm, niv)
      end if
!*****************************************************************************
      do ii=i1_v, i2_v
      do kk=1   , ngzm
        rrtemp(kk,ii) = rrz(ii,kk)
      end do
      end do

      rr(1:ngzm, i1_v:i2_v, jj) = rrtemp(1:ngzm,i1_v:i2_v)

3000 continue 

! (1-A2*)
	deallocate(aa,bb,cc)
      allocate(aa(i1_v:i2_v,j1_v:j2_v),bb(i1_v:i2_v,j1_v:j2_v),cc(i1_v:i2_v,j1_v:j2_v)) 
     ! allocate(aa_g(15, 15), bb_g(15, 15), cc_g(15, 15)) ! for printing all the aa

      aa=0.0; bb=0.0; cc=0.0 
      do 93   jj=j1_v, j2_v
      do 93   ii=i1_v, i2_v
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!---------------------------------------
!-------(eq 6.21) first term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =svetax (i,j)
        a2 =svetay (i,j)
        b1 =svretax(i,j)
        b2 =svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(svetax(i,j)+svetax(i,j+1))/2.
        d2 =(svetay(i,j)+svetay(i,j+1))/2.
        tta=vp(i,j)
        e1 =(d1*c1+d2*c4)*       visc/tta 
        e2 =(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i,j+1)+vp(i,j))/2.)
        cc(ii,jj)=beta*(svetax(i,j+1)*e1+svetay(i,j+1)*e2)

!---------------------------------------
!-------(eq 6.21) second term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =svetax (i,j)
        a2 =svetay (i,j)
        b1 =svretax(i,j)
        b2 =svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(svetax(i,j)+svetax(i,j+1))/2.
        d2 =(svetay(i,j)+svetay(i,j+1))/2.
        tta=vp(i,j)
        e1 =(d1*c1+d2*c4)*       visc/tta
        e2 =(d1*c2+d2*c5)*       visc/tta 
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(svetax(i,j)+svetax(i,j-1))/2.
        d2 =(svetay(i,j)+svetay(i,j-1))/2.
        tta=vp(i,j-1)
        e1 =e1+(d1*c1+d2*c4)*       visc/tta
        e2 =e2+(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)
        bb(ii,jj)=1.-beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------(eq 6.21) third term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =svetax (i,j)
        a2 =svetay (i,j)
        b1 =svretax(i,j)
        b2 =svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(svetax(i,j)+svetax(i,j-1))/2.
        d2 =(svetay(i,j)+svetay(i,j-1))/2.
        tta=vp(i,j-1)
        e1 =(d1*c1+d2*c4)*       visc/tta 
        e2 =(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i,j-1)+vp(i,j-2))/2.)
        aa(ii,jj)=beta*(svetax(i,j-1)*e1+svetay(i,j-1)*e2)

93 continue
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idvdy_1.eq.0) then
        do 1003 ii=i1_v,i2_v
          bb(ii,j1_v)=bb(ii,j1_v)+aa(ii,j1_v)
 1003   continue
      end if

      if(idvdy_2.eq.0) then
        do 1004 ii=i1_v,i2_v
          bb(ii,j2_v)=bb(ii,j2_v)+cc(ii,j2_v)
 1004   continue
      end if
! MZW in current adjoint code only this part is executed
      if (jblock.eq.1) then
        do 89 ii=i1_v,i2_v
          aa(ii,j1_v)=0.
 89     continue
      end if

      if (jblock.eq.npy) then
        do 95 ii=i1_v,i2_v
          cc(ii,j2_v)=0.
 95     continue
      end if

!*****************************************************************************
!      do i=i1_v,i2_v
!        call triDiagonal_mrhs(id_y, aa(i,j1_v), bb(  i,j1_v),             &
!                                    cc(i,j1_v), rr(1,i,j1_v), ngzm, njv)
!      end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
      call triDiagonalM_mrhs(id_y, aa(i1_v,j1_v), bb(  i1_v,j1_v),        &
                                   cc(i1_v,j1_v), rr(1,i1_v,j1_v), njv, niv, ngzm)

!temp = 0.0
!	do jj=j1_v,j2_v
!		temp = temp + rr(1,i1_v,jj) * rhs(1,i1_v,jj)
!	end do
!	call globalSum(temp,1)
!	if (irank == 1) then
!	write(*,*) '----------adjoint operator---------'
!	write(*,*) temp
!	end if

	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{eta,2}'
!		write(*,*) rr(1,i1_v,:)
!	end if
!*****************************************************************************

! (1-A1*)      
      aa=0.0; bb=0.0; cc=0.0 
      do 75   jj=j1_v, j2_v
      do 75   ii=i1_v, i2_v
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!---------------------------------------
!--------(eq 6.20) first term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(suxix  (i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2 =(suxiy  (i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1 = svretax(i,j)
        b2 = svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(suxix(i+1,j)+suxix(i+1,j-1))/2.
        d2 =(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        tta=(vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        e1 =(d1*c1+d2*c4)*       visc/tta
        e2 =(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i+1,j)+vp(i+1,j-1))/2.)
	cc(ii,jj)=beta*(svetax(i+1,j)*e1+svetay(i+1,j)*e2)

!---------------------------------------
!--------(eq 6.20) second term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(suxix  (i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2 =(suxiy  (i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1 = svretax(i,j)
        b2 = svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(suxix(i+1,j)+suxix(i+1,j-1))/2.
        d2 =(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        tta=(vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        e1 =(d1*c1+d2*c4)*       visc/tta 
        e2 =(d1*c2+d2*c5)*       visc/tta 
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(suxix(i,j)+suxix(i,j-1))/2.
        d2 =(suxiy(i,j)+suxiy(i,j-1))/2.
        tta=(vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        e1 =e1+(d1*c1+d2*c4)*       visc/tta 
        e2 =e2+(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)
	bb(ii,jj)=1.-beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!--------(eq 6.20) third term-----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(suxix  (i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2 =(suxiy  (i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1 = svretax(i,j)
        b2 = svretay(i,j)
        c1 =(a1*b1+b1*a1)
        c2 =(a1*b2+b1*a2)
        c4 =(a2*b1+b2*a1)
        c5 =(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =(suxix(i,j)+suxix(i,j-1))/2.
        d2 =(suxiy(i,j)+suxiy(i,j-1))/2.
        tta=(vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        e1 =(d1*c1+d2*c4)*       visc/tta 
        e2 =(d1*c2+d2*c5)*       visc/tta 
	beta=-dt/(2.*(vp(i-1,j)+vp(i-1,j-1))/2.)
        aa(ii,jj)=beta*(svetax(i-1,j)*e1+svetay(i-1,j)*e2)

75 continue 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (BC_left.ne.2) then
      if(idvdx_1.eq.0) then
        do 1001 jj=j1_v,j2_v
          bb(i1_v,jj)=bb(i1_v,jj)+aa(i1_v,jj)
 1001   continue
      end if

      if(idvdx_2.eq.0) then
        do 1002 jj=j1_v,j2_v
          bb(i2_v,jj)=bb(i2_v,jj)+cc(i2_v,jj)
 1002   continue
      end if

      if (iblock.eq.1) then
        do 81 jj=j1_v,j2_v
          aa(i1_v ,jj)=0.
 81     continue
      end if

      if (iblock.eq.npx) then
        do 98 jj=j1_v,j2_v
          cc(i2_v ,jj)=0.
 98     continue
      end if
end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! QW: Solve periodic system
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if(BC_left .eq. 2)then

        do jj=j1_v,j2_v
                call triPeriodic_mrhs(id_x, aa(i1_v,jj), bb(i1_v,jj),                &
                                            cc(i1_v,jj), rr(1,i1_v,jj), ngz-1, niv)
        enddo
	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{eta,1}'
!		write(*,*) rr(1,i1_v,:)
!	end if
end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! TJ: Solve non-periodic system
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (BC_left.ne.2) then

	rr(:,:,:) = rhs(1:ngz-1, i1_v:i2_v, j1_v:j2_v)

	do jj=j1_v,j2_v
		call triDiagonal_mrhs(id_x, aa(i1_v,jj), bb(i1_v,jj),                &
                                            cc(i1_v,jj), rr(1,i1_v,jj), ngz-1, niv)
        end do
end if

! check adjoint relation
!temp = 0.0
!	do jj=j1_v,j2_v
!	do ii=i1_v,i2_v
!	do kk=1,ngzm
!		temp = temp + rr(kk,ii,jj) * rhs(kk,ii,jj)
!	end do
!	end do
!	end do
!	call globalSum(temp,1)
!	if (irank == 1) then
!	write(*,*) '----------adjoint operator---------'
!	write(*,*) temp
!	end if

! update adjoint vst here !!!
	vst_adj(1:ngzm,i1_v:i2_v,j1_v:j2_v) = rr(1:ngzm,i1_v:i2_v,j1_v:j2_v)
	! Check ust_adj
	!write(*,*) 'After invert, ust at kk=1, ii=1'
	!write(*,*) vst_adj(1,1,j1_v:j2_v)

	! MZW 02/07/17 print coeff matrix to check
	!write(*,*) '---------Start printing adjoint operator---------'
	!do ii=i1_v, i2_v
	!	write(*,*) 'Start printing aa at ii=', ii
	!	write(*,*) aa(ii,j1_v:j2_v)
	!	write(*,*) 'Start printing bb at ii=', ii
	!	write(*,*) bb(ii,j1_v:j2_v)
	!	write(*,*) 'Start printing cc at ii=', ii
	!	write(*,*) cc(ii,j1_v:j2_v)
	!end do
	!data_gather = 0.0
	!data_gather(i1_v:i2_v,j1_v:j2_v) = aa(:,:)
	!call gather_abc(aa,aa_g,5,5)
	!data_gather(i1_v:i2_v,j1_v:j2_v) = bb(:,:)
	!call gather_abc(bb,bb_g,5,5)
	!data_gather(i1_v:i2_v,j1_v:j2_v) = cc(:,:)
	!call gather_abc(cc,cc_g,5,5)
	!if (irank == iroot) then
	!	write(*,*) '---------Start printing adjoint operator---------'
	!	write(*,*)'Start printing global aa at ii=1'
	!	write(*,*) aa_g(1,1:15)
	!	write(*,*)'Start printing global bb at ii=1'
	!	write(*,*) bb_g(1,1:15)
	!	write(*,*)'Start printing global cc at ii=1'
	!	write(*,*) cc_g(1,1:15)
	!	
	!end if
end


subroutine invertw_Adjoint()
      use data_export

      real, dimension(:,:)  , allocatable :: aa, bb, cc, rrz, rrtemp
      !real, dimension(:,:)  , allocatable :: aa_g, bb_g, cc_g
      real, dimension(ngzm,i1_w:i2_w,j1_w:j2_w) :: rr

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MZW--calculate the rhs of adjoint equation (6.171) except GPxi_adj term
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	rhs = 0.0

	dt2_inv = 1.0/dt**2
	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
	do kk=1,ngzm
		! here GPxi_adj = F^{n+1}_x * dx / V, dx/V has been included in calculate_GPxi
		rhs(kk,ii,jj) =  w_adj(kk,ii,jj) + (phi_adj(kk-1,ii,jj) - phi_adj(kk,ii,jj)) * dt2_inv
		!rhs(kk,ii,jj) = cos(real(ii+jj+kk))
	end do
	end do
	end do

	rhsw(:,:,:) = rhs(:,:,:)
	! Check rhs
	!write(*,*) 'rhs of invertw_Adjoint at ii=1, jj=1'
	!write(*,*) rhs(1:ngzm,1,1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SY--determine the starting index in the z direction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(BC_front.eq.2) then
        ka=1  ! periodic bc
      else
        ka=2
      end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SY--initialize the temporary rhs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr=0.0
	! set values for rhs, test tri-diagonal solver
	! MZW ---------------
	!do jj=j1_w, j2_w
	!	rhs(1:ngzm,i1_w:i2_w,jj) = jj * 1.E-6
	!end do
	! --------------------
	rr(ka:ngz-1, i1_w:i2_w, j1_w:j2_w) = rhs(ka:ngz-1, i1_w:i2_w, j1_w:j2_w)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A1*)(1-A2*)(1-A3*)(wstar-w(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A3*)
      allocate(aa(i1_w:i2_w,ka:ngzm),bb (i1_w:i2_w,ka:ngzm)) 
      allocate(cc(i1_w:i2_w,ka:ngzm),rrz(i1_w:i2_w,ka:ngzm)) 
      allocate(rrtemp(ka:ngzm,i1_w:i2_w)) 
      aa=0.0; bb=0.0; cc=0.0; rrz=0.0; rrtemp=0.0

      do 3000 jj=j1_w, j2_w

        rrtemp(ka:ngzm, i1_w:i2_w) = rr(ka:ngzm, i1_w:i2_w, jj)

      do ii=i1_w, i2_w
      do kk=ka  , ngzm
        rrz(ii,kk) = rrtemp(kk,ii)
      end do
      end do

      do kk=ka  , ngzm
      do ii=i1_w, i2_w
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)
        k    = kk

!---------------------------------------
!---------(eq 6.25)-first term----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{z,P}+S_{z,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=swz(i,j)
        b3=swrz(i,j)
        c9=(a3*b3+b3*a3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{z,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =swz(i,j)
        tta=vp(i,j)
        e3 =(d3*c9)*       visc/tta
	beta=-dt/(2.*vp(i,j))
        cc(ii,kk)=beta*swz(i,j)*e3

!---------------------------------------
!---------(eq 6.25)-second term---------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{z,P}+S_{z,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=swz(i,j)
        b3=swrz(i,j)
        c9=(a3*b3+b3*a3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{f}/(ReV_{f}) \bu d_{z,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =swz(i,j)
        tta=vp(i,j)
        e3 =(d3*c9)*       visc/tta
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{z,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =swz(i,j)
        tta=vp(i,j)
        e3 =e3+(d3*c9)*       visc/tta
	beta=-dt/(2.*vp(i,j))
        bb(ii,kk)=1.-beta*swz(i,j)*e3

!---------------------------------------
!---------(eq 6.25)-third term----------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{P}S_{z,P}+S_{z,P}S^{z}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=swz(i,j)
        b3=swrz(i,j)
        c9=(a3*b3+b3*a3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{z}_{b}/(ReV_{b}) \bu d_{z,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d3 =swz(i,j)
        tta=vp(i,j)
        e3 =(d3*c9)*       visc/tta
	beta=-dt/(2.*vp(i,j))
        aa(ii,kk)=beta*swz(i,j)*e3

      end do ! ii
      end do ! kk
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(BC_front.ne.2) then
        do ii=i1_w,i2_w
          !bb(ii,ka  )=bb(ii,ka  )+aa(ii,ka  )
          !bb(ii,ngzm)=bb(ii,ngzm)+cc(ii,ngzm)
          aa(ii,ka  )=0.
          cc(ii,ngzm)=0.
        end do
      end if

!*****************************************************************************
      if(BC_front.eq.2) then
      call triPeriodicM_(id_z, aa(i1_w,ka), bb (i1_w,ka),        &
                               cc(i1_w,ka), rrz(i1_w,ka), ngz-ka, niw)
      end if
	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{z,3} at (ii,jj)=',i1_w,jj
!		write(*,*) rrz(i1_w,1:ngzm)
!	end if

      if(BC_front.ne.2) then
      call triDiagonalM_(id_z, aa(i1_w,ka), bb (i1_w,ka),        &
                               cc(i1_w,ka), rrz(i1_w,ka), ngz-ka, niw)
      end if
!*****************************************************************************
      do ii=i1_w, i2_w
      do kk=ka  , ngzm
        rrtemp(kk,ii) = rrz(ii,kk)
      end do
      end do

      rr(ka:ngzm, i1_w:i2_w, jj) = rrtemp(ka:ngzm,i1_w:i2_w)

3000 continue

! (1-A2*)
      deallocate(aa,bb,cc) 
      allocate(aa(i1_w:i2_w,j1_w:j2_w),bb(i1_w:i2_w,j1_w:j2_w),cc(i1_w:i2_w,j1_w:j2_w)) 
      !allocate(aa_g(ngx-1, ngy-1), bb_g(ngx-1, ngy-1), cc_g(ngx-1, ngy-1)) ! for printing all the aa

aa=0.0;bb=0.0;cc=0.0

      do 96   jj=j1_w, j2_w
      do 96   ii=i1_w, i2_w
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!---------------------------------------
!----------(eq 6.24)--first term--------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(svetax(i,j)+svetax(i,j+1))/2.
        a2 =(svetay(i,j)+svetay(i,j+1))/2.
        b3 = swrz  (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{z,n,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1 =svetax(i,j+1)
        d2 =svetay(i,j+1)
        tta=(vp(i,j)+vp(i,j+1))/2. 
        e3 =(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i,j+1))
        cc(ii,jj)=beta*(swz(i,j+1)*e3) 

!---------------------------------------
!----------(eq 6.24)-second term--------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(svetax(i,j)+svetax(i,j+1))/2.
        a2 =(svetay(i,j)+svetay(i,j+1))/2.
        b3 = swrz  (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{z,n,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1 =svetax(i,j+1)
        d2 =svetay(i,j+1)
        tta=(vp(i,j)+vp(i,j+1))/2. 
        e3 =(d1*c3+d2*c6)*       visc/tta 
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =svetax(i,j)
        d2 =svetay(i,j)
        tta=(vp(i,j)+vp(i,j-1))/2.
        e3 =e3+(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i,j))
        bb(ii,jj)=1.-beta*(swz(i,j)*e3)

!---------------------------------------
!----------(eq 6.24)--third term--------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1 =(svetax(i,j)+svetax(i,j+1))/2.
        a2 =(svetay(i,j)+svetay(i,j+1))/2.
        b3 = swrz  (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1 =svetax(i,j)
        d2 =svetay(i,j)
        tta=(vp(i,j)+vp(i,j-1))/2.
        e3 =(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i,j-1))
        aa(ii,jj)=beta*swz(i,j-1)*e3


96 continue 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idwdy_1.eq.0) then
        do 1003 ii=i1_w,i2_w
          bb(ii,j1_w)=bb(ii,j1_w)+aa(ii,j1_w)
 1003   continue
      end if

      if(idwdy_2.eq.0) then
        do 1004 ii=i1_w,i2_w
          bb(ii,j2_w)=bb(ii,j2_w)+cc(ii,j2_w)
 1004   continue
      end if
! in current adjoint code (periodic), only this part is executed
      if (jblock.eq.1) then
        do 112 ii=i1_w,i2_w
          aa(ii,j1_w)=0.
 112    continue
      end if

      if (jblock.eq.npy) then
        do 121 ii=i1_w,i2_w
          cc(ii,j2_w)=0.
 121    continue
      end if

!*****************************************************************************
!   do i=i1_w,i2_w
!     call triDiagonal_mrhs(id_y, aa(i,j1_w), bb(  i,j1_w),               &
!                                 cc(i,j1_w), rr(1,i,j1_w), ngzm, njw)
!   end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
      call triDiagonalM_mrhs(id_y, aa(i1_w,j1_w), bb(   i1_w,j1_w),        &
                                   cc(i1_w,j1_w), rr(ka,i1_w,j1_w), njw, niw, ngz-ka)   
	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{z,2}'
!		write(*,*) rr(1,i1_w,:)
!	end if
!*****************************************************************************

! (1-A1*)
      aa=0.0; bb=0.0; cc=0.0 

      do 75   jj=j1_w, j2_w
      do 75   ii=i1_w, i2_w
        i    = ibmap_1(ii)
        j    = jbmap_1(jj)

!---------------------------------------
!--------(eq 6.23 first term)-----------
!---------------------------------------
        a1 =(suxix(i,j)+suxix(i+1,j))/2.
        a2 =(suxiy(i,j)+suxiy(i+1,j))/2.
        b3 = swrz (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{z,e,im}
!ccccccccccccccccccccccccccccccccccccccccccc
        d1 =suxix(i+1,j)
        d2 =suxiy(i+1,j)
        tta=(vp(i,j)+vp(i+1,j))/2. 

        e3=(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i+1,j))
        cc(ii,jj)=beta*(swz(i+1,j)*e3)

!---------------------------------------
!-------(eq 6.23 second term)-----------
!---------------------------------------
        a1 =(suxix(i,j)+suxix(i+1,j))/2.
        a2 =(suxiy(i,j)+suxiy(i+1,j))/2.
        b3 = swrz (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{z,e,im}
!ccccccccccccccccccccccccccccccccccccccccccc
        d1 = suxix(i+1,j)
        d2 = suxiy(i+1,j)
        tta=(vp(i,j)+vp(i+1,j))/2. 
        e3 =(d1*c3+d2*c6)*       visc/tta 
!cccccccccccccccccccccccccccccccccccccccccccc
!       -S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1 =suxix(i,j)
        d2 =suxiy(i,j)
        tta=(vp(i,j)+vp(i-1,j))/2.
        e3 =e3+(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i,j))
        bb(ii,jj)=1.-beta*(swz(i,j)*e3)

!---------------------------------------
!--------(eq 6.23 third term)-----------
!---------------------------------------
        a1 =(suxix(i,j)+suxix(i+1,j))/2.
        a2 =(suxiy(i,j)+suxiy(i+1,j))/2.
        b3 = swrz (i,j)
        c3 =(a1*b3)
        c6 =(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       -S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1 =suxix(i,j)
        d2 =suxiy(i,j)
        tta=(vp(i,j)+vp(i-1,j))/2.
        e3 =(d1*c3+d2*c6)*       visc/tta 
	beta=-dt/(2.*vp(i-1,j))
        aa(ii,jj)=beta*(swz(i-1,j)*e3)

75 continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (BC_left.ne.2) then
      if(idwdx_1.eq.0) then
        do 1001 jj=j1_w,j2_w
          bb(i1_w,jj)=bb(i1_w,jj)+aa(i1_w,jj)
 1001   continue
      end if

      if(idwdx_2.eq.0) then
        do 1002 jj=j1_w,j2_w
          bb(i2_w,jj)=bb(i2_w,jj)+cc(i2_w,jj)
 1002   continue
      end if

      if (iblock.eq.1) then
        do 92 jj=j1_w,j2_w
          aa(i1_w, jj)=0.
 92     continue
      end if

      if (iblock.eq.npx) then
        do 94 jj=j1_w,j2_w
          cc(i2_w,jj)=0.
 94     continue
      end if
end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! QW: Solve periodic system 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if(BC_left .eq. 2)then

        !rr(ka:ngz-1, i1_w:i2_w, j1_w:j2_w) = rhs(ka:ngz-1, i1_w:i2_w, j1_w:j2_w)

        do jj=j1_w,j2_w
                call triPeriodic_mrhs(id_x, aa(i1_w,jj), bb(i1_w,jj),                &
                                            cc(i1_w,jj), rr(ka,i1_w,jj), ngz-ka, niw)
        enddo
	! MZW 02/14/17 check results
!	if (irank == iroot) then
!		write(*,*) 'Writing TRI results of A_{z,1}'
!		write(*,*) rr(1,i1_w,:)
!	end if
end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! TJ: Solve non-periodic system
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if(BC_left.ne.2) then

	!rr(ka:ngz-1, i1_w:i2_w, j1_w:j2_w) = rhs(ka:ngz-1, i1_w:i2_w, j1_w:j2_w)

	do jj=j1_w,j2_w
		call triDiagonal_mrhs(id_x, aa(i1_w,jj), bb(i1_w,jj),                &
                                            cc(i1_w,jj), rr(ka,i1_w,jj), ngz-ka, niw)
	end do
end if

! check adjoint relation
!temp = 0.0
!	do jj=j1_w,j2_w
!	do ii=i1_w,i2_w
!	do kk=1,ngzm
!		temp = temp + rr(kk,ii,jj) * rhs(kk,ii,jj)
!	end do
!	end do
!	end do 
!	!call globalSum(temp,1)
!	if (irank == 1) then
!	write(*,*) '----------adjoint operator---------'
!	write(*,*) temp
!	end if
!! update adjoint wst here !!!
	wst_adj(1:ngzm,i1_w:i2_w,j1_w:j2_w) = rr(1:ngzm,i1_w:i2_w,j1_w:j2_w)
!	! Check ust_adj
!	write(*,*) 'After invert, wst at ii=1, jj=1'
!	write(*,*) wst_adj(1:ngzm,1,1)

	! MZW 02/07/17 print coeff matrix to check
	! write(*,*) '---------Start printing adjoint operator---------'
	! do ii=i1_w, i2_w
	! 	write(*,*) 'irank=', irank, ' : Start printing aa at ii=', ii
	!	write(*,*) aa(ii,j1_w:j2_w)
	!	write(*,*) 'irank=', irank, ' : Start printing bb at ii=', ii
	!	write(*,*) bb(ii,j1_w:j2_w)
	!	write(*,*) 'irank=', irank, ' : Start printing cc at ii=', ii
	!	write(*,*) cc(ii,j1_w:j2_w)
	!end do
	!call gather_abc(aa,aa_g,niw,njw)
	!call gather_abc(bb,bb_g,niw,njw)
	!call gather_abc(cc,cc_g,niw,njw)
	!if (irank == iroot) then
	!	write(*,*) '---------Start printing adjoint operator---------'
	!	write(*,*)'Start printing global aa at ii=1'
	!	write(*,*) aa_g(1,1:ngy-1)
	!	write(*,*)'Start printing global bb at ii=1'
	!	write(*,*) bb_g(1,1:ngy-1)
	!	write(*,*)'Start printing global cc at ii=1'
	!	write(*,*) cc_g(1,1:ngy-1)
	!	
	!end if
end







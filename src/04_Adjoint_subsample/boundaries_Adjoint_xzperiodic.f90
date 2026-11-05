!=======================================================================
! MZW 02/14/17
! Boundary conditions for channel flow adjoint variables

!=======================================================================
subroutine CartesianBC_Adjoint(deltaT) 
! originally copied from CartesianBC(deltaT)

      	use boundaries
     
      	!===================================================================
      	! 	NO SLIP B.C. 
      	!===================================================================
	!--------------------------------------------------------------
	!	Bottom B.C. 
	!--------------------------------------------------------------
	if (jblock.eq.1) then
		
		uc_adj(:, :, 0) = -uc_adj(:, :, 1) 
		vc_adj(:, :, 1) =  0.0            
		wc_adj(:, :, 0) = -wc_adj(:, :, 1)
	
		!Extend
		vc_adj(:, :, 0) = vc_adj(:, :, 2)        
        	
	end if

	!--------------------------------------------------------------
	!	 Top B.C.
	!-------------------------------------------------------------- 
	if (jblock.eq.npy) then
        	
        	uc_adj(:, :, nlyb  ) = -uc_adj(:, :, nlyb-1)
		vc_adj(:, :, nlyb  ) =  0.0
		wc_adj(:, :, nlyb  ) = -wc_adj(:, :, nlyb-1)
         		
		!Extend
		vc_adj(:, :, nlyb+1) = vc_adj(:, :, nlyb-1) 
        	
	end if

	!-----------------------------------------------------------------------
        !       Periodic in streamwise (x) direction
        !-----------------------------------------------------------------------
	if (iblock.eq.1 .or. iblock.eq.npx) then
		if(npx.eq.1) then 
		
			uleftbc_temp(:,:,1) = uc_adj(:,1,:)
                        uleftbc_temp(:,:,2) = uc_adj(:,2,:)
		       	vleftbc_temp(:,:,1) = vc_adj(:,1,:)
                        wleftbc_temp(:,:,1) = wc_adj(:,1,:)
	
			urightbc_temp(:,:,1) = uc_adj(:,nlxb-1,:)
			urightbc_temp(:,:,2) = uc_adj(:, nlxb, :)
		       	vrightbc_temp(:,:,1) = vc_adj(:,nlxb-1,:)
                        wrightbc_temp(:,:,1) = wc_adj(:,nlxb-1,:)
	
		else
			! ==> (iblock.eq.1  ) receives ucbce
			! ==> (iblock.eq.npx) receives ucbcw
			call updateBorder_Adjoint_xperiodic(1)
		end if
	end if

        if (iblock.eq.1) then
		uc_adj(:,0,:) = urightbc_temp(:,:,1)			   
		uc_adj(:,1,:) = urightbc_temp(:,:,2)
		vc_adj(:,0,:) = vrightbc_temp(:,:,1)			   
                wc_adj(:,0,:) = wrightbc_temp(:,:,1)           	  
	end if

        if (iblock.eq.npx) then
		
		vc_adj(:, nlxb, :) = vleftbc_temp(:,:,1)		   
                wc_adj(:, nlxb, :) = wleftbc_temp(:,:,1)		 
 	
		!Extend
		uc_adj(:,nlxb+1,:) = uleftbc_temp(:,:,2)		  
 
        end if

        !-----------------------------------------------------------------------
        !       Periodic in spanwise (z) direction
        !-----------------------------------------------------------------------
        uc_adj(0    , :, :) = uc_adj(ngz-1, :, :)
        uc_adj(ngz  , :, :) = uc_adj(1    , :, :)

        vc_adj(0    , :, :) = vc_adj(ngz-1, :, :)
        vc_adj(ngz  , :, :) = vc_adj(1    , :, :)

        wc_adj(0    , :, :) = wc_adj(ngz-1, :, :)
        wc_adj(ngz  , :, :) = wc_adj(1    , :, :)
        !DANGER ZONE: wc(  1  , :,:) = 0.5*(wc(1, :,:)+wc(ngz, :,:)) (wc only goes from 1:ngz-1 so wc(ngz,:,:) does not exist)
	!Extend
        wc_adj(ngz+1, :, :) = wc_adj(2    , :, :)

      	!----------------- Block boundaries -------------------
        !-- Exchange halos in x-direction , then in y-direction
        
	!-- x-direction
	call updateBorder_lim(uc_adj, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
        call updateBorder_lim(vc_adj, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
        call updateBorder_lim(wc_adj, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

        !-- y-direction
        call updateBorder_lim(uc_adj, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
        call updateBorder_lim(vc_adj, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
        call updateBorder_lim(wc_adj, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)
      	
	return
end
!============================================================================================

subroutine FluxBC_Adjoint()

!     Determine the BCs in terms of flux.     

      	use boundaries
      	real*8  :: vtempb(ngz+1), vtempe(ngz+1)
      	real*8  :: utempb(ngz+1), utempe(ngz+1)
      	real*8, dimension(2) :: Atmp

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		u_adj(:, :, 0) = -u_adj(:, :, 1)
       		v_adj(:, :, 1) =  0.0
       		w_adj(:, :, 0) = -w_adj(:, :, 1)

       	      	!Extend
		v_adj(:, :, 0) = v_adj(:, :, 2)	
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	u_adj(:, :, nlyb  ) = -u_adj(:, :, nlyb-1)
              	v_adj(:, :, nlyb  ) =  0.0
             	w_adj(:, :, nlyb  ) = -w_adj(:, :, nlyb-1)
		
		!Extend
		v_adj(:, :, nlyb+1) = v_adj(:, :, nlyb-1)	
     		
	end if 
	
	!-----------------------------------------------------------------------
        !	Periodic in streamwise (x) direction     
        !-----------------------------------------------------------------------
        if (iblock.eq.1 .or. iblock.eq.npx) then
                if (npx.eq.1) then

                        uleftbc_temp(:,:,1) = u_adj(:,1,:)
                        uleftbc_temp(:,:,2) = u_adj(:,2,:)
                        vleftbc_temp(:,:,1) = v_adj(:,1,:)
                        !SJL: necessary for diffusion term in u-equation
                        vleftbc_temp(:,:,2) = v_adj(:,2,:)
                        wleftbc_temp(:,:,1) = w_adj(:,1,:)

                        urightbc_temp(:,:,1) = u_adj(:,nlxb-1,:)
                        urightbc_temp(:,:,2) = u_adj(:, nlxb, :)
                        vrightbc_temp(:,:,1) = v_adj(:,nlxb-1,:)
                        wrightbc_temp(:,:,1) = w_adj(:,nlxb-1,:)

                else
                        call updateBorder_Adjoint_xperiodic(2)
                end if

        end if

        if (iblock.eq.1) then

                u_adj(:,0,:) = urightbc_temp(:,:,1)
                u_adj(:,1,:) = urightbc_temp(:,:,2)
                v_adj(:,0,:) = vrightbc_temp(:,:,1)
                w_adj(:,0,:) = wrightbc_temp(:,:,1)

        end if

        if (iblock.eq.npx) then

                v_adj(:, nlxb, :) = vleftbc_temp(:,:,1)
                w_adj(:, nlxb, :) = wleftbc_temp(:,:,1)

                !Extend
                u_adj(:,nlxb+1,:) = uleftbc_temp(:,:,2)
		!Extend by MZW
                v_adj(:,nlxb+1,:) = vleftbc_temp(:,:,2)

        end if

	!-----------------------------------------------------------------
        !	Periodic in spanwise (z) direction
        !-----------------------------------------------------------------
        u_adj( 0 , :,:) = u_adj(ngz-1, :,:)
        u_adj(ngz, :,:) = u_adj(  1  , :,:)
        
	v_adj( 0 , :,:)=v_adj(ngz-1, :,:)
        v_adj(ngz, :,:)=v_adj(  1  , :,:)

        w_adj(  0  , :,:) = w_adj(ngz-1, :,:)
        !TAZ DANGER ZONE: w(  1  , :,:) = 0.5*(w(1, :,:)+w(ngz, :,:))
        w_adj(ngz  , :,:) = w_adj(1, :,:)
	!Extend
        w_adj(ngz+1, :,:) = w_adj(2, :,:)

        !=======================================================================
        ! sum of mass flux along boundaries is zero
        !=======================================================================
        bmassin =0.0
        bmassout=0.0
        if (jblock.eq.1  ) then
        do i=i1_T,i2_T
           bmassin  = bmassin  + sum(v_adj(1:ngzm, i, 1   )) !TJ-- equal to zero?
        end do
        end if

        if (jblock.eq.npy) then
        do i=i1_T,i2_T
           bmassin  = bmassin  - sum(v_adj(1:ngzm, i, nlyb)) !TJ-- equal to zero?
        end do
        end if

        bmass = bmassin
        call globalSum(bmass,1)

        if (iblock.eq.1  ) then
        do j=j1_T,j2_T
           bmassin  = bmassin  + sum(u_adj(1:ngzm, 1   , j))
        end do
        end if

        if (iblock.eq.npx) then
        do j=j1_T,j2_T
           bmassout = bmassout + sum(u_adj(1:ngzm, nlxb, j))
        end do
        end if

        do 710 j=j1_T,j2_T
        do 710 i=i1_T,i2_T
           bmassin  = bmassin  + w_adj(1,i,j)-w_adj(ngz,i,j) !TJ-- equal to zero?
   710  continue 

        Atmp(1) = bmassin
        Atmp(2) = bmassout
        call globalSum(Atmp, 2)
        bmassin  = Atmp(1)
        bmassout = Atmp(2)

        if (irank.eq.1)   then
                !write(6,*) 'Adjoint mass inflow  =  ', bmassin
                !write(6,*) 'Adjoint mass outflow =  ', bmassout
!!! (for debug MZW)                 write(6,'(a26,es14.7,a12,es14.7)') 'Before Adjust Adjoint mass inflow = ', bmassin, ', outflow = ', bmassout !JL151026
        end if

        !if(abs(bmassout).ge.1.e-10) then 
        !   if (iblock.eq.npx) then
        !       u(1:ngzm,nlxb,1:nly-1)=u(1:ngzm,nlxb,1:nly-1)*bmassin/bmassout
        !   end if
        !end if

        bmassout=0.0
        if (iblock.eq.npx) then
           do j=j1_T,j2_T
              bmassout = bmassout + sum(u_adj(1:ngzm, nlxb, j))
           end do
        end if

        call globalSum(bmassout, 1)
      
	if (irank.eq.1) then
                !write(6,*) 'After adjustment'
                !write(6,*) 'Mass inflow  =  ', bmassin
                !write(6,*) 'Mass outflow =  ', bmassout
!!!                write(6,'(a26,es14.7,a12,es14.7)') 'After  Adjust Adjoint mass inflow = ', bmassin, ', outflow = ', bmassout !JL151026
!!! (for debug MZW)                 write(6,*) 'adjoint mass corrected: inflow - outflow =  ',bmassout-bmassin
        end if

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: this is incorrect for streamwise periodicity
	!
	!if (iblock.eq.npx) then
        !	!Extend
        !	u(:, nlxb+1, 0:nly) = 2.0*u(:, nlxb, 0:nly)-u(:, nlxb-1, 0:nly)
        !end if
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: these boundary conditions have been enforced already
	!
        !TAZ DANGER ZONE:  Xiahua has u-periodicity here:
	!if (BC_front.eq.2) then
        !	u( 0 , :, :) = u(ngz-1, :, :)
        !	u(ngz, :, :) = u(1    , :, :)
        !end if
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

!----------------- Block boundaries -------------------
      	!-- Exchange halos in x-direction , then in y-direction
      	!-- x-direction
	call updateBorder_lim(u_adj, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
      	call updateBorder_lim(v_adj, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
      	call updateBorder_lim(w_adj, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

      	!-- y-direction
      	call updateBorder_lim(u_adj, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
      	call updateBorder_lim(v_adj, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
      	call updateBorder_lim(w_adj, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

	 !print*, "u inside FluxBC"
         !do k=0,ngz
         !do j=0,nly
         	!print *, "u(", k, ",", nlxb+1, ",", j, ") :" ,u(k, nlxb+1, j)
         !end do
         !end do

return
end

subroutine FluxBC_Adjoint_st()

!     Determine the BCs in terms of flux.     

      	use boundaries
      	real*8  :: vtempb(ngz+1), vtempe(ngz+1)
      	real*8  :: utempb(ngz+1), utempe(ngz+1)
      	real*8, dimension(2) :: Atmp

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		ust_adj(:, :, 0) = -ust_adj(:, :, 1)
       		vst_adj(:, :, 1) =  0.0
       		wst_adj(:, :, 0) = -wst_adj(:, :, 1)

       	      	!Extend
		vst_adj(:, :, 0) = vst_adj(:, :, 2)	
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
         	
           	ust_adj(:, :, nlyb  ) = -ust_adj(:, :, nlyb-1)
              	vst_adj(:, :, nlyb  ) =  0.0
             	wst_adj(:, :, nlyb  ) = -wst_adj(:, :, nlyb-1)
		
		!Extend
		vst_adj(:, :, nlyb+1) = vst_adj(:, :, nlyb-1)	
     		
	end if 
	
	!-----------------------------------------------------------------------
        !	Periodic in streamwise (x) direction     
        !-----------------------------------------------------------------------
        if (iblock.eq.1 .or. iblock.eq.npx) then
                if (npx.eq.1) then

                        uleftbc_temp(:,:,1) = ust_adj(:,1,:)
                        uleftbc_temp(:,:,2) = ust_adj(:,2,:)
                        vleftbc_temp(:,:,1) = vst_adj(:,1,:)
                        !SJL: necessary for diffusion term in u-equation
                        vleftbc_temp(:,:,2) = vst_adj(:,2,:)
                        wleftbc_temp(:,:,1) = wst_adj(:,1,:)

                        urightbc_temp(:,:,1) = ust_adj(:,nlxb-1,:)
                        urightbc_temp(:,:,2) = ust_adj(:, nlxb, :)
                        vrightbc_temp(:,:,1) = vst_adj(:,nlxb-1,:)
                        wrightbc_temp(:,:,1) = wst_adj(:,nlxb-1,:)

                else
                        call updateBorder_Adjoint_xperiodic(3)	! 3 for ust_adj, vst_adj, wst_adj
                end if

        end if

        if (iblock.eq.1) then

                ust_adj(:,0,:) = urightbc_temp(:,:,1)
                ust_adj(:,1,:) = urightbc_temp(:,:,2)
                vst_adj(:,0,:) = vrightbc_temp(:,:,1)
                wst_adj(:,0,:) = wrightbc_temp(:,:,1)

        end if

        if (iblock.eq.npx) then

                vst_adj(:, nlxb, :) = vleftbc_temp(:,:,1)
                wst_adj(:, nlxb, :) = wleftbc_temp(:,:,1)

                !Extend
                ust_adj(:,nlxb+1,:) = uleftbc_temp(:,:,2)
		!Extend by MZW
                vst_adj(:,nlxb+1,:) = vleftbc_temp(:,:,2)

        end if

	!-----------------------------------------------------------------
        !	Periodic in spanwise (z) direction
        !-----------------------------------------------------------------
        ust_adj( 0 , :,:) = ust_adj(ngz-1, :,:)
        ust_adj(ngz, :,:) = ust_adj(  1  , :,:)
        
	vst_adj( 0 , :,:)=vst_adj(ngz-1, :,:)
        vst_adj(ngz, :,:)=vst_adj(  1  , :,:)

        wst_adj(  0  , :,:) = wst_adj(ngz-1, :,:)
        !TAZ DANGER ZONE: w(  1  , :,:) = 0.5*(w(1, :,:)+w(ngz, :,:))
        wst_adj(ngz  , :,:) = wst_adj(1, :,:)
	!Extend
        wst_adj(ngz+1, :,:) = wst_adj(2, :,:)
	

!----------------- Block boundaries -------------------
      	!-- Exchange halos in x-direction , then in y-direction
      	!-- x-direction
	call updateBorder_lim(ust_adj, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
      	call updateBorder_lim(vst_adj, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
      	call updateBorder_lim(wst_adj, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

      	!-- y-direction
      	call updateBorder_lim(ust_adj, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
      	call updateBorder_lim(vst_adj, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
      	call updateBorder_lim(wst_adj, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

	 !print*, "u inside FluxBC"
         !do k=0,ngz
         !do j=0,nly
         	!print *, "u(", k, ",", nlxb+1, ",", j, ") :" ,u(k, nlxb+1, j)
         !end do
         !end do

return
end

!=======================================================================
! MZW 02/14/17
! Message-passing routines for parallel processing adjoint variables using MPI
! updateBorder_Adjoint_xperiodic(flag)

!=======================================================================

subroutine updateBorder_Adjoint_xperiodic(flag)
! originally copied from updateBorder_xperiodic(flag)
! only the case flag = 1, 2  is changed for adjoint variables
! other cases need to be edited later

        use messenger
        include "mpif.h"

        integer flag, src_dst, icount
        real*8, allocatable :: buf1(:,:,:), buf2(:,:,:)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        if (iblock.eq.1  ) idims(1) = npx-1
        if (iblock.eq.npx) idims(1) = 0
        idims(2) = jblock-1
        idims(3) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        if (flag.eq.1) then

                allocate(buf1(1:ngzm,0:nly,4))
                allocate(buf2(1:ngzm,0:nly,4))
                icount = ngzm*(nly+1)*4

                if (iblock.eq.1) then

                        buf1(:,:,1) = uc_adj(1:ngzm,1,0:nly)
                        buf1(:,:,2) = uc_adj(1:ngzm,2,0:nly)
                        buf1(:,:,3) = vc_adj(1:ngzm,1,0:nly)
                        buf1(:,:,4) = wc_adj(1:ngzm,1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

                        urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
                        urightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
                        vrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
                        wrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

                endif

                if (iblock.eq.npx) then

                        buf1(:,:,1) = uc_adj(1:ngzm,nlxb-1,0:nly)
                        buf1(:,:,2) = uc_adj(1:ngzm, nlxb, 0:nly)
                        buf1(:,:,3) = vc_adj(1:ngzm,nlxb-1,0:nly)
                        buf1(:,:,4) = wc_adj(1:ngzm,nlxb-1,0:nly)


                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

                        uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
                        uleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
                        vleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
                        wleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

                endif

	else if (flag.eq.2) then
		!SJL: I need 5 slots instead of 4 because we need to know v beyond halos
		!     at the right global boundaries, which is available from the left 
		!     global boundaries by periodic boundary conditions.
		allocate(buf1(1:ngzm,0:nly,5))
		allocate(buf2(1:ngzm,0:nly,5))
		icount = ngzm*(nly+1)*5

		if (iblock.eq.1) then

			buf1(:,:,1) = u_adj(1:ngzm,1,0:nly)
			buf1(:,:,2) = u_adj(1:ngzm,2,0:nly)
			buf1(:,:,3) = v_adj(1:ngzm,1,0:nly)
			!SJL: 4th slot of buf1 is necessary for fixing "out-of-advance-bound"
			!     of v in explicit diffusion term in xi-momentum equation.
			!     (see u-advancement subroutine in "advancement_quick_..._.f90"
			!     This is only needed for periodic boundary conditions in xi-dir.
			buf1(:,:,4) = v_adj(1:ngzm,2,0:nly) 
			buf1(:,:,5) = w_adj(1:ngzm,1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			urightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
			!SJL: receive the dummy because it is not needed in iblock=1.
			vrightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,4)
			wrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,5)

                endif  

                if (iblock.eq.npx) then

			buf1(:,:,1) = u_adj(1:ngzm,nlxb-1,0:nly)
			buf1(:,:,2) = u_adj(1:ngzm, nlxb, 0:nly)
			buf1(:,:,3) = v_adj(1:ngzm,nlxb-1,0:nly)
			!SJL: dummy is sent to iblock=1 because it is not needed there.
			buf1(:,:,4) = 0.0
			buf1(:,:,5) = w_adj(1:ngzm,nlxb-1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			uleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
                        vleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,4)
			wleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,5)

                endif

	else if (flag.eq.3) then
		!SJL: I need 5 slots instead of 4 because we need to know v beyond halos
		!     at the right global boundaries, which is available from the left 
		!     global boundaries by periodic boundary conditions.
		allocate(buf1(1:ngzm,0:nly,5))
		allocate(buf2(1:ngzm,0:nly,5))
		icount = ngzm*(nly+1)*5

		if (iblock.eq.1) then

			buf1(:,:,1) = ust_adj(1:ngzm,1,0:nly)
			buf1(:,:,2) = ust_adj(1:ngzm,2,0:nly)
			buf1(:,:,3) = vst_adj(1:ngzm,1,0:nly)
			!SJL: 4th slot of buf1 is necessary for fixing "out-of-advance-bound"
			!     of v in explicit diffusion term in xi-momentum equation.
			!     (see u-advancement subroutine in "advancement_quick_..._.f90"
			!     This is only needed for periodic boundary conditions in xi-dir.
			buf1(:,:,4) = vst_adj(1:ngzm,2,0:nly) 
			buf1(:,:,5) = wst_adj(1:ngzm,1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			urightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
			!SJL: receive the dummy because it is not needed in iblock=1.
			vrightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,4)
			wrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,5)

                endif  

                if (iblock.eq.npx) then

			buf1(:,:,1) = ust_adj(1:ngzm,nlxb-1,0:nly)
			buf1(:,:,2) = ust_adj(1:ngzm, nlxb, 0:nly)
			buf1(:,:,3) = vst_adj(1:ngzm,nlxb-1,0:nly)
			!SJL: dummy is sent to iblock=1 because it is not needed there.
			buf1(:,:,4) = 0.0
			buf1(:,:,5) = wst_adj(1:ngzm,nlxb-1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			uleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
                        vleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,4)
			wleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,5)

                endif


	else if (flag.eq.4) then

			allocate(buf1(1:ngzm,0:nly,1))
			allocate(buf2(1:ngzm,0:nly,1))
			icount = ngzm*(nly+1)*1

		if (iblock.eq.1) then

			buf1(:,:,1) = p(1:ngzm,1,0:nly)	

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
		end if

		if (iblock.eq.npx) then
			buf1(:,:,1) = p(1:ngzm,nlxb-1,0:nly)

                   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                     buf2, icount, MPI_REAL8, src_dst, 0, &
                                     icomm_grid, istatus, ierr)
		
			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)	

		end if

        endif

        deallocate(buf1)
        deallocate(buf2)

    return
end

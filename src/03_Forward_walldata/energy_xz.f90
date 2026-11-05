! calculate perturbation kinetic energy
! averaged in x and z direction
! updated 09/16/17
! calculate energy at staggered grid first
! then average at cell center
subroutine energy_xz()
	use data_export
	use messenger
	include "mpif.h"
	real Uavg(0:nly), Vavg(0:nly+1), Wavg(0:nly)
	real :: umid, vmid, wmid, xL, yL, zL	! velocities at pressure node (center of cell C.V.)
	real :: TKEmid, TKEsum			! Total Kinetic Energy at pressure node (center of cell C.V.)

	! Average in the homogeneous x and z direction
	call readFloat("xL",xL)
	call readFloat("yL",yL)
	call readFloat("zL",zL)
	Uavg = 0.0
	Vavg = 0.0
	Wavg = 0.0
	if(iblock == 1) then
		istart = i1_u !- 1
	else 
		istart = i1_u
	end if
	do jj=j1_u,j2_u
		do ii=istart,i2_u
		do kk=1,ngzm
			Uavg(jj) = Uavg(jj) + uc(kk,ii,jj)
		end do
		end do
	end do
	do jj=j1_v,j2_v
		do ii=i1_v,i2_v
		do kk=1,ngzm
			Vavg(jj) = Vavg(jj) + vc(kk,ii,jj)
		end do
		end do
	end do
	do jj=j1_w,j2_w
		do ii=i1_w,i2_w
		do kk=1,ngzm
			Wavg(jj) = Wavg(jj) + wc(kk,ii,jj)
		end do
		end do
	end do
	call globalDirSum(Uavg,nly+1,1) ! sum in streamwise direction
	call globalDirSum(Vavg,nly+2,1)
	call globalDirSum(Wavg,nly+1,1)
	Uavg = Uavg/dble((ngx-1)*(ngz-1))
	Vavg = Vavg/dble((ngx-1)*(ngz-1))
	Wavg = Wavg/dble((ngx-1)*(ngz-1))
	! Kinetic energy
	TKEsum = 0.0
	do jj=j1_T,j2_T
	do ii=i1_T,i2_T
	do kk=1,ngzm
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		k = kk
		uu = 0.5*((uc(kk,ii,jj) - Uavg(jj))**2 + (uc(kk,ii+1,jj) - Uavg(jj))**2)
		vv = 0.5*((vc(kk,ii,jj) - Vavg(jj))**2 + (vc(kk,ii,jj+1) - Vavg(jj+1))**2)
		ww = 0.5*((wc(kk,ii,jj) - Wavg(jj))**2 + (wc(kk+1,ii,jj) - Wavg(jj))**2)
		TKEmid = 0.5*(uu + vv + ww)
		TKEsum = TKEsum + TKEmid * vp(i,j)
	end do
	end do
	end do
	call globalSum(TKEsum,1)
	TKEsum = TKEsum/(xL*yL*zL)
	if(irank==1) then
		write(111,*) stime, TKEsum
	end if
end

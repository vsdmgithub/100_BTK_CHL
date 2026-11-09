! MZW
! calculate the cost function and save the value
! this subroutine should be called :
! 1. before main time loop;
! 2. at the end of each time loop

subroutine cost_function()
	!use data_export
	use mesh_export
	use messenger
	include "mpif.h"
	integer :: flag_cost, my_jstart, my_jend, global_jndx, loc_jndx, jjout1
	real :: cost_single, cost_dudy, cost_dwdy, cost_p, cost_subdomain, weight_f, weight_p, weight_t

	flag_cost = 0
	cost_single = 0.0
        cost_dudy = 0.0
        cost_dwdy = 0.0
        cost_p = 0.0
        cost_subdomain = 0.0


        !--------------------------skin friction and pressure----------------
        !weight_f = visc**2
	!!weight_f = 0

	!my_jstart = jTmin_1(jblock)
	!my_jend   = jTmax_1(jblock)

	!jjout1 = 1
	!global_jndx = jjout1
	!loc_jndx    = jmap_1(global_jndx)
	!if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
	!	flag_cost = 1
	!end if

	!! Here p is p*dt !!!!!!!!!!!!!!!!!!!!	
	!if (flag_cost == 1) then
        !        if(ipre == 1) then
	!	do ii=i1_T,i2_T
	!	do kk=1,ngzm
	!		i = ibmap_1(ii)
	!		k = kk
	!		cost_p = cost_p + 0.5 *  (GPxi/dt*(xm(i)-xL) + &
	!					p(k,ii,loc_jndx)/dt - p_obs(k,ii,1)/dt)**2
	!		!cost_single = cost_single + 0.5  * (p(k,ii,loc_jndx)/dt - p_obs(k,ii,1)/dt)**2
	!	end do
	!	end do
        !        end if

        !        if(ifric == 1) then
	!	do ii=i1_u,i2_u
	!	do kk=1,ngzm
	!		k = kk
	!		cost_dudy = cost_dudy + 0.5 * weight_f * (dudy_tmp(kk,ii,1) - dudy_obs(kk,ii,1))**2
	!	end do
	!	end do

	!	do ii=i1_w,i2_w
	!	do kk=1,ngzm
	!		k = kk
	!		cost_dwdy = cost_dwdy + 0.5 * weight_f * (dwdy_tmp(kk,ii,1) - dwdy_obs(kk,ii,1))**2
	!	end do
	!	end do
        !        end if
	!end if

	!jjout1 = ngy-1
	!global_jndx = jjout1
	!loc_jndx    = jmap_1(global_jndx)
	!if ( global_jndx.ge.my_jstart  .and.  global_jndx.le.my_jend ) then	! MZW: should be ge and le!!!
	!	flag_cost = 2
	!end if
	!
	!if (flag_cost == 2) then
        !        if(ipre == 1) then
	!	do ii=i1_T,i2_T
	!	do kk=1,ngzm
	!		i = ibmap_1(ii)
	!		k = kk
	!		cost_p = cost_p + 0.5 * (GPxi/dt*(xm(i)-xL) + &
	!					p(k,ii,loc_jndx)/dt - p_obs(k,ii,2)/dt)**2
	!		!cost_single = cost_single + 0.5 * (p(k,ii,loc_jndx)/dt - p_obs(k,ii,2)/dt)**2
	!	end do
	!	end do
        !        end if

        !        if(ifric == 1) then
	!	do ii=i1_u,i2_u
	!	do kk=1,ngzm
	!		k = kk
	!		cost_dudy = cost_dudy + 0.5 * weight_f * (dudy_tmp(kk,ii,2) - dudy_obs(kk,ii,2))**2
	!	end do
	!	end do

	!	do ii=i1_w,i2_w
	!	do kk=1,ngzm
	!		k = kk
	!		cost_dwdy = cost_dwdy + 0.5 * weight_f * (dwdy_tmp(kk,ii,2) - dwdy_obs(kk,ii,2))**2
	!	end do
	!	end do
        !        end if
	!end if

        !--------------------------subdomain--------------------------
        if (isubdomain == 1) then
                !weight_t = exp(-2*0.021*180.0**2/2800.0*stime)
                !weight_t = 1.0
                !if(irank==1) write(*,*) 'weight_t=',weight_t
                do jj=j1_u,j2_u
                do ii=i1_u,i2_u
                do kk=1,ngzm
                        i = ibmap_1(ii)
                        j = jbmap_1(jj)
                        k = kk
                        !cost_subdomain = cost_subdomain + 0.5*uflag(k,ii,jj)*weight_t*surxix(i,j)**2*(u(k,ii,jj) - uob(k,ii,jj))**2
                        cost_subdomain = cost_subdomain + 0.5 * uflag(k,ii,jj)*surxix(i,j)**2*(u(k,ii,jj) - uob(k,ii,jj))**2
                end do
                end do
                end do
                do jj=j1_v,j2_v
                do ii=i1_v,i2_v
                do kk=1,ngzm
                        i = ibmap_1(ii)
                        j = jbmap_1(jj)
                        k = k
                        !cost_subdomain = cost_subdomain + 0.5 *vflag(k,ii,jj)*weight_t*svretay(i,j)**2*(v(k,ii,jj) - vob(k,ii,jj))**2
                        cost_subdomain = cost_subdomain + 0.5 * vflag(k,ii,jj)*svretay(i,j)**2*(v(k,ii,jj) - vob(k,ii,jj))**2
                end do
                end do
                end do
                do jj=j1_w,j2_w
                do ii=i1_w,i2_w
                do kk=1,ngzm
                        i = ibmap_1(ii)
                        j = jbmap_1(jj)
                        k = kk
                        !cost_subdomain = cost_subdomain + 0.5 * wflag(k,ii,jj)*weight_t*swrz(i,j)**2*(w(k,ii,jj) - wob(k,ii,jj))**2
                        cost_subdomain = cost_subdomain + 0.5 * wflag(k,ii,jj)*swrz(i,j)**2*(w(k,ii,jj) - wob(k,ii,jj))**2
                end do
                end do
                end do

        end if



	!call globalSum(cost_single,1)
        call globalSum(cost_p,1)
        call globalSum(cost_dudy,1)
        call globalSum(cost_dwdy,1)
        call globalSum(cost_subdomain,1)
        !write(33,'(4E16.9)') stime, cost_p, cost_dudy, cost_dwdy
	cost = cost + cost_p + cost_dudy + cost_dwdy + cost_subdomain
        write(33,*) stime, cost
        if(irank==1) write(*,*) 'ntime, cost=',ntime,cost
	flag_cost = 0
end

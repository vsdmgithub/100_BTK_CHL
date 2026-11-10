 subroutine calc_vort_mag
 use mvariable
 use msetup
 use mgeomdata
 use mindices
 use mopenmpcut
 implicit none
 real,dimension(3,3)::sym,asym,st2om2,eigenvec
 real,dimension(3  )::eigenval
 integer::l,m,ierr

 integer::thd
 integer::im,i,ip,jm,j,jp,km,k,kp
 real   ::dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
 real   ::dxm,dxp,dym,dyp
 real   ::uim,uip,ujm,uj,ujp,ukm,ukp
 real   ::vim,vip,vjm,vj,vjp,vkm,vkp
 real   ::wim,wip,wjm,wj,wjp,wkm,wkp

 write(*,*) 'Calculation of the vorticity magnitude...'
 ens(:,:, 1   ) = 0.0
 ens(:,:,nybox) = 0.0

!$omp parallel do
!$omp& private(
!$omp& thd,sym,asym,st2om2,eignvec
!$omp& ,im,i,ip,jm,j,jp,km,k,kp
!$omp& ,dxm,dxp,dym,dyp
!$omp& ,uim,uip,ujm,uj,ujp,ukm,ukp
!$omp& ,vim,vip,vjm,vj,vjp,vkm,vkp
!$omp& ,wim,wip,wjm,wj,wjp,wkm,wkp
!$omp& ,dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
!$omp& )
  do thd = 1 , nthds
   do j = jstart(thd) , jend(thd)
   if( j == 1 .or. j == nybox ) go to 200
    jp = j + 1
    jm = j - 1
    do i = 2 , nxboxm
     ip = i + 1
     im = i - 1
     do k = 1 , nzboxm
      kp = k + 1
      km = k - 1
      if(kp > nzboxm) kp = kp - nzboxm
      if(km < 1     ) km = km + nzboxm 

!     uij_{11} = dudxn = du / dx
      dxm = xpg(i_3d_full(i ),j_3d_full(j)) - xpg(i_3d_full(im),j_3d_full(j))
      dxp = xpg(i_3d_full(ip),j_3d_full(j)) - xpg(i_3d_full(i ),j_3d_full(j))

      uim = u(k,im,j)
      uip = u(k,ip,j)
      dudx1 = ( uip - uim ) / ( dxm + dxp )

!     uij_{21} = dvdxn = dv / dx
      vim = v(k,im,j)
      vip = v(k,ip,j)
      dvdx1 = ( vip - vim ) / ( dxm + dxp )

!     uij_{31} = dwdxn = dw / dx
      wim = w(k,im,j)
      wip = w(k,ip,j)
      dwdx1 = ( wip - wim ) / ( dxm + dxp )

!     uij_{12} = dudyn = du / dy
      dym = ypg(i_3d_full(i),j_3d_full(j )) - ypg(i_3d_full(i),j_3d_full(jm))
      dyp = ypg(i_3d_full(i),j_3d_full(jp)) - ypg(i_3d_full(i),j_3d_full(j ))

      ujm = u(k,i,jm)
      uj  = u(k,i,j )
      ujp = u(k,i,jp)
      dudy1 =  ( ujp * dym**2 - ujm * dyp**2 + uj * ( dyp**2 - dym**2 ) ) &
             / ( dyp * dym * ( dym + dyp ) )

!     uij_{22} = dvdyn = dv / dy
      vjm = v(k,i,jm)
      vj  = v(k,i,j )
      vjp = v(k,i,jp)
      dvdy1 =  ( vjp * dym**2 - vjm * dyp**2 + vj * ( dyp**2 - dym**2 ) ) &
             / ( dyp * dym * ( dym + dyp ) )

!     uij_{32} = dwdyn = dw / dy
      wjm = w(k,i,jm)
      wj  = w(k,i,j )
      wjp = w(k,i,jp)
      dwdy1 =  ( wjp * dym**2 - wjm * dyp**2 + wj * ( dyp**2 - dym**2 ) ) &
             / ( dyp * dym * ( dym + dyp ) )

!     uij_{13} = dudzn = du / dz
      ukm = u(km,i,j)
      ukp = u(kp,i,j)
      dudz1 = ( ukp - ukm ) / ( 2. * dzbox )      ! order of dz^2

!     uij_{23} = dvdzn = dv / dz  or  dv / (r dtheta) - w / r
      vkm = v(km,i,j)
      vkp = v(kp,i,j)
      dvdz1 = ( vkp - vkm ) / ( 2. * dzbox )

!     uij_{33} = dwdzn = dw / dz  or  dw / (r dtheta) + v / r
      wkm = w(km,i,j)
      wkp = w(kp,i,j)
      dwdz1 = ( wkp - wkm ) / ( 2. * dzbox )

!---------------Calculate RAMDA -----------------------
!   asym(i,j) = Omega_ij, sym(i,j)=S_ij
! st2om2(i,j) = sym(i,k)*sym(k,j)+asym(i,k)*asym(k,j)

      asym(1,1)= 0.0; asym(2,2)= 0.0; asym(3,3)=0.0
      asym(1,2)= 0.5 * ( dudy1 - dvdx1 )
      asym(1,3)= 0.5 * ( dudz1 - dwdx1 )
      asym(2,3)= 0.5 * ( dvdz1 - dwdy1 )
      asym(2,1)= -asym(1,2)
      asym(3,1)= -asym(1,3)
      asym(3,2)= -asym(2,3)
!      sym (1,1)= 1.0 * ( dudx1         )
!      sym (2,2)= 1.0 * ( dvdy1         )
!      sym (3,3)= 1.0 * ( dwdz1         )
!      sym (1,2)= 0.5 * ( dudy1 + dvdx1 )
!      sym (1,3)= 0.5 * ( dudz1 + dwdx1 )
!      sym (2,3)= 0.5 * ( dvdz1 + dwdy1 )
!      sym (2,1)= sym(1,2)
!      sym (3,1)= sym(1,3)
!      sym (3,2)= sym(2,3)

      ens(k,i,j) =  sqrt(               &
                           asym(1,2)**2 &
                         + asym(1,3)**2 &
                         + asym(2,3)**2 &
                                        )

     end do
    end do
200 continue
   end do
  end do

 return
 end subroutine calc_vort_mag

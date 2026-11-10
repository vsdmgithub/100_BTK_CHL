 subroutine calc_lamci
 use mvariable
 use msetup
 use mgeomdata
 use mindices
 use mopenmpcut
 implicit none
 integer::thd
 real   ::uij(3,3),residual_ci,resimax_ci,lam_ci2
 real,dimension(:),allocatable::resimax_ci_thd
 integer::im,i,ip,jm,j,jp,km,k,kp
 real   ::dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
 real   ::dxm,dxp,dym,dyp
 real   ::uim,uip,ujm,uj,ujp,ukm,ukp
 real   ::vim,vip,vjm,vj,vjp,vkm,vkp
 real   ::wim,wip,wjm,wj,wjp,wkm,wkp

!$omp parallel do
!$omp& private(i,j,k)
!$omp& shared(lam,nxbox,nybox,nzbox)
!$omp& default(none)
 do j = 1 , nybox
 do i = 1 , nxbox
 do k = 1 , nzbox
  lam(k,i,j) = 0.0
 end do
 end do
 end do

 allocate( resimax_ci_thd( nthds ) )

 write(*,*) 'Calculation of the swirling strength...'

!$omp parallel do
!$omp& private(
!$omp& thd,uij,lam_ci2,residual_ci
!$omp& ,im,i,ip,jm,j,jp,km,k,kp
!$omp& ,dxm,dxp,dym,dyp
!$omp& ,uim,uip,ujm,uj,ujp,ukm,ukp
!$omp& ,vim,vip,vjm,vj,vjp,vkm,vkp
!$omp& ,wim,wip,wjm,wj,wjp,wkm,wkp
!$omp& ,dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
!$omp& )
  do thd = 1 , nthds
   resimax_ci_thd(thd) = 0.0
   do j = jstart(thd) , jend(thd)
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


      uij(1,1) = dudx1
      uij(2,1) = dvdx1
      uij(3,1) = dwdx1

      uij(1,2) = dudy1
      uij(2,2) = dvdy1
      uij(3,2) = dwdy1

      uij(1,3) = dudz1
      uij(2,3) = dvdz1
      uij(3,3) = dwdz1

!!!!! Zhou, Adrian, Balachandra and Kendall (1999) !!!!!!
      call eigenval33_lambdaci(uij,lam_ci2,residual_ci)

      if( residual_ci > resimax_ci_thd(thd) ) then
       resimax_ci_thd(thd) = residual_ci
      end if

      lam(k,i,j) = sqrt( lam_ci2 )

     end do
    end do
   end do
  end do
!$omp end parallel do


  resimax_ci = 0.0
  do thd = 1 , nthds
   if ( resimax_ci_thd( thd ) > resimax_ci ) then
    resimax_ci = resimax_ci_thd( thd )
   end if
  end do

  write(*,*) 'Newton-Rhapson Method Residual'
  write(*,*) 'lambda ci residual = ', resimax_ci

  do j = 1 , nyboxm
   do i = 1 , nxboxm
    lam(nzbox,i,j) = lam(1,i,j)
   end do
  end do

 return
 end subroutine calc_lamci


!*******************************************************
!START SUBROUTINE EIGENVALUE CALCULATION for LAMBDA CI *
!*******************************************************
 subroutine eigenval33_lambdaci(co,x_ci2,y)
 implicit none
 integer::indx
 real::co(3,3),x_ci2,x_cr,a(0:3),x,y,yp,xl,xu,yl,yu

!     for incompressible flow
!     a(3)*x^3+a(2)*x^2+a(1)*x+a(0)=0
!     for real x
      a(3) = -1.d0
      a(2) =  0.d0 !+co(1,1) + co(2,2) + co(3,3)
      a(1) = -co(1,1) * co(2,2)            &
             -co(2,2) * co(3,3)            &
             -co(1,1) * co(3,3)            &
             +co(2,1) * co(1,2)            &
             +co(2,3) * co(3,2)            &
             +co(3,1) * co(1,3)
      a(0) = +co(1,1) * co(2,2) * co(3,3)  &
             +co(2,1) * co(3,2) * co(1,3)  &
             -co(1,1) * co(2,3) * co(3,2)  &
             +co(3,1) * co(1,2) * co(2,3)  &
             -co(3,1) * co(1,3) * co(2,2)  &
             -co(2,1) * co(1,2) * co(3,3)

!     Newton-Rhapson Method
      indx = 1
      x = 0.d0
      y = a(3) * x ** 3  &
        + a(2) * x ** 2  &
        + a(1) * x       &
        + a(0)
      do
       yp = 3.d0 *a(3) *x ** 2 &
          + 2.d0 *a(2) *x      &
          + a(1)
!       if(yp==0.d0) yp=1.d0
       x = x - y / yp
       y = a(3) * x ** 3       &
         + a(2) * x ** 2       &
         + a(1) * x            &
         + a(0)
       if ( ( abs( y ) <= 1.d-15 ) .or. ( indx == 100 ) ) exit
       indx = indx + 1
      end do
      y = abs( y )

!     If Newton-Rhapson Method is failed => False-Position Method
      if ( y >= 1.d-10 ) then
       x = 0.d0
       y = a(3) * x ** 3  &
         + a(2) * x ** 2  &
         + a(1) * x       &
         + a(0)
       if ( y > 0.d0 ) then
        xl = x
        yl = y
        xu = 0.d0
        fpm1 : do
        xu = xu + 1.d0
        yu = a(3) * xu ** 3  &
           + a(2) * xu ** 2  &
           + a(1) * xu       &
           + a(0)
        if ( yu < 0.d0 ) exit fpm1
        end do fpm1
       end if

       if ( y < 0.d0 ) then
        xu = x
        yu = y
        xl = 0.d0
        fpm2 : do
        xl = xl - 1.d0
        yl = a(3) * xl ** 3  &
           + a(2) * xl ** 2  &
           + a(1) * xl       &
           + a(0) 
        if ( yl > 0.d0 ) exit fpm2
        end do fpm2
       end if

       indx = 1
       fpm3 : do
        x = xu - yu * ( xl - xu ) / ( yl -yu )
        y = a(3) * x ** 3  &
          + a(2) * x ** 2  &
          + a(1) * x       &
          + a(0) 
        if ( abs( y ) <= 1.d-12 .or. indx == 1000000 ) exit fpm3
        if ( y > 0.d0 ) then
         xl = x
         yl = y
        else
         xu = x
         yu = y
        end if
       indx = indx + 1
       end do fpm3
      end if

      y = abs( y )

!     x+x_cr+i*x_ci+x_cr-i*x_ci = -a(2) / a(3)
      x_cr = -0.5d0 * x

!     x*(x_cr**2+x_ci**2) = -a(0) /a(3)
      x_ci2 = - a(0) / a(3) / x - x_cr ** 2

      if ( x_ci2 < 0.0 ) x_ci2 = 0.0

      return
 end subroutine eigenval33_lambdaci
!=======================================================
!END SUBROUTINE EIGENVALUE CALCULATION for LAMBDA CI ===
!=======================================================

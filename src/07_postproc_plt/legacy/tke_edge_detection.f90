subroutine tke_edge_detection
      use mgeomdata
      use mvariable
      use msetup
      use mindices
      use mreadavg
      implicit none
      integer::i,j,k,im,ip,jm,jp,imore,navg
      real::tke_jc,tke_jm,tke_jp,tke_im,tke_ip

!=============================================================
! 1) Calculation of TKE based on U_infty
!=============================================================
!$omp  parallel do
!$omp& private(i,j,k)
!$omp& shared(u,v,w,tke_calc,tke_calc_xy,Uy,i_3d_full,j_3d_full,nxbox,nybox,nzbox)
!$omp& default(none)
      do j = 1 , nybox
      do i = 1 , nxbox
      do k = 1 , nzbox
       tke_calc(k,i,j)                                            &
              = ( ( u(k,i,j)-1.0                               ) **2  & ! U_inf = 1.0
              +   ( v(k,i,j)-Uy(i_3d_full(i),j_3d_full(nybox)) ) **2  & 
              +     w(k,i,j)                                     **2 ) * 0.5
       tke_calc_xy(k,i,j)                                         &
              = ( ( u(k,i,j)-1.0                               ) **2  &  ! U_inf = 1.0
              +   ( v(k,i,j)-Uy(i_3d_full(i),j_3d_full(nybox)) ) **2 ) * 0.5
!       tke_calc(k,i,j) = (u(k,i,j)-1.0)**2 & ! U_inf = 1.0
!                       +  v(k,i,j)     **2 &
!                       +  w(k,i,j)     **2
!       tke_calc_xy(k,i,j) = (u(k,i,j)-1.0)**2 & ! U_inf = 1.0
!                         +  v(k,i,j)     **2  
      end do
      end do
      end do

!      call writexy(tke_calc,1,'tke_raw',imore)

!=============================================================
! 2) Adjusting exception at the boundary (to avoid zero)
!=============================================================
!$omp  parallel do
!$omp& private(i,k)
!$omp& shared(tke_calc,tke_calc_xy,nxbox,nzbox)
!$omp& default(none)
      do i = 1 , nxbox
      do k = 1 , nzbox
       tke_calc(k,i,1) = tke_calc(k,i,2)
       tke_calc_xy(k,i,1) = tke_calc_xy(k,i,2)
      end do
      end do

!$omp  parallel do
!$omp& private(j,k)
!$omp& shared(tke_calc,tke_calc_xy,nxbox,nybox,nzbox)
!$omp& default(none)
      do j = 1 , nybox
      do k = 1 , nzbox
       tke_calc(k,1    ,j) = tke_calc(k,2      ,j)
       tke_calc(k,nxbox,j) = tke_calc(k,nxbox-1,j)
       tke_calc_xy(k,1    ,j) = tke_calc_xy(k,2      ,j)
       tke_calc_xy(k,nxbox,j) = tke_calc_xy(k,nxbox-1,j)
      end do
      end do


!=============================================================
! 2.5) Spatio- and temporal-averaging
!=============================================================
      do j = 1 , nybox
      do i = 1 , nxbox
      do k = 1 , nzbox
       mean_tke(i,j) = mean_tke(i,j) + tke_calc(k,i,j)
       mean_tke_xy(i,j) = mean_tke_xy(i,j) + tke_calc_xy(k,i,j)
      end do
      end do
      end do

!=============================================================
! 3) Filtering of single-point area
!=============================================================
!$omp  parallel do
!$omp& private(i,j,k,im,ip,jm,jp,tke_jc,tke_jp,tke_jm,tke_ip,tke_im)
!$omp& shared(tke_calc,tke_th,nxbox,nybox,nzbox)
!$omp& default(none)
      do j = 2 , nybox - 1
      do i = 2 , nxbox - 1
      do k = 1 , nzbox
       im = i - 1
       ip = i + 1

       jm = j - 1
       jp = j + 1

       tke_jc = tke_calc(k,i,j )
       tke_jp = tke_calc(k,i,jp)
       tke_jm = tke_calc(k,i,jm)
       tke_ip = tke_calc(k,ip,j)
       tke_im = tke_calc(k,im,j)

       if(       tke_jc <= tke_th   &
           .and. tke_jp >  tke_th   &
           .and. tke_jm >  tke_th   &
           .and. tke_ip >  tke_th   &
           .and. tke_im >  tke_th ) &
        tke_calc(k,i,j) = 0.25 * ( tke_jp + tke_jm + tke_ip + tke_im )
       if(       tke_jc >= tke_th   &
           .and. tke_jp <  tke_th   &
           .and. tke_jm <  tke_th   &
           .and. tke_ip <  tke_th   &
           .and. tke_im <  tke_th ) &
        tke_calc(k,i,j) = 0.25 * ( tke_jp + tke_jm + tke_ip + tke_im )
      end do
      end do
      end do

!      call writexy(tke_calc,1,'tke_fsp',imore)

!=============================================================
! 4) Edge profile along x
!=============================================================
      do i = 1 , nxbox
      do k = 1 , nzbox

       Yedge(i,k) = ypg(i,j_3d_full(nybox))
       Jedge(i,k) =                 nybox

       do j = nybox - 1 , 2 , -1
        if( tke_calc(k,i,j) >= tke_th ) then
         Yedge(i,k) = ypg(i,j_3d_full(j))
         Jedge(i,k) =                 j
         go to 100
        end if
       end do ! j-loop
100    continue

       Jedgelow(i,k) = Jedge(i,k) - 10 ! temp
       Jedgeupp(i,k) = Jedge(i,k) +  5 ! temp

      end do ! k-loop
      end do ! i-loop

      return
end subroutine tke_edge_detection

subroutine tke_edge_thickness
      use mgeomdata
      use mvariable
      use msetup
      use mindices
      use mreadavg
      implicit none
      integer::i,j,k,im,ip,jm,jp,imore

      lam(:,:,:) = 0.0

!=============================================================
      do i = 1 , nxbox
      do k = 1 , nzbox

       do j = Jedgelow(i,k) , Jedgeupp(i,k)
        lam(k,i,j) = 1.0
       end do ! j-loop

      end do ! k-loop
      end do ! i-loop

      return
end subroutine tke_edge_thickness

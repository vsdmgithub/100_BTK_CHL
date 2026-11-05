!234567
      subroutine indices
      use mindices
      use msetup
      implicit none
      integer::i,j,k
      allocate( i_3d_full( nxbox ) )
      allocate( j_3d_full( nybox ) )
      allocate( k_3d_full( nzbox ) )

      do i = 1 , nxbox
       i_3d_full( i ) = i
      end do

      do j = 1 , nybox
       j_3d_full( j ) = j
      end do

      do k = 1 , nzbox
       k_3d_full( k ) = k
      end do


      return
      end subroutine indices

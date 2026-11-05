subroutine writeins_tke_edge_2d(imore)
      use mgeomdata
      use mvariable
      use msetup
      use mindices
      use mreadavg
      implicit none
      integer::i,j,k,ii,jj,kk,imore
      character(70)::outfname
      integer::TecIni111,TecDat,TecZne,TecEnd,TecNod,TecFil
      integer::O,OO,Debug,VIsDouble,FileType
      double precision, dimension(:,:,:),allocatable:: OutArray
      character(6)::cidx

      FileType = 2    ! 0 = full, 1 = grid, 2 = solution
      Debug = 0       ! 0 = no debug & 1 = debug
      VIsDouble = 0   ! 0 = single & 1 = double precision

      write(cidx,fmt='(i6.6)') imore

      !----------- TECPLOT FILE -----------
      outfname = './OUTPUT/ins_tke_edge_XY_'
      O=TecIni111('InsField'//char(0), & ! title
               'tke_edge'//char(0), &    ! variables
               trim(outfname)//cidx//'.plt'//char(0), & ! Fname
               '.'//char(0), &           ! ScratchDir
               FileType,Debug,VIsDouble);

      !Write the zone header information
      O=TecZne('SOULTION'//char(0), &
               nxbox,nybox,1,'BLOCK'//char(0),char(0));
      OO=nxbox*nybox

      allocate( OutArray(nxbox,nybox,1) )
      !--- Output velocity data ----
      do j = 1, nybox
      do i = 1, nxbox
        OutArray(i,j,1) = tke_calc(kref,i,j)
      end do
      end do
      O=TecDat(OO,OutArray(1,1,1),1);

      O=TecEnd();

      deallocate( OutArray )

      !----------- TECPLOT FILE -----------
      outfname = './OUTPUT/ins_tke_edge_XZ_'
      O=TecIni111('InsField'//char(0), & ! title
               'tke_edge'//char(0), &    ! variables
               trim(outfname)//cidx//'.plt'//char(0), & ! Fname
               '.'//char(0), &           ! ScratchDir
               FileType,Debug,VIsDouble);

      !Write the zone header information
      O=TecZne('SOULTION'//char(0), &
               nxbox,1,nzbox,'BLOCK'//char(0),char(0));
      OO=nxbox*nzbox

      allocate( OutArray(nxbox,1,nzbox) )
      !--- Output velocity data ----
      do k = 1, nzbox
      do i = 1, nxbox
        OutArray(i,1,k) = tke_calc(k,i,jref)
      end do
      end do
      O=TecDat(OO,OutArray(1,1,1),1);

      O=TecEnd();

      deallocate( OutArray )

      !----------- TECPLOT FILE -----------
      outfname = './OUTPUT/ins_tke_edge_YZ_'
      O=TecIni111('InsField'//char(0), & ! title
               'tke_edge'//char(0), &    ! variables
               trim(outfname)//cidx//'.plt'//char(0), & ! Fname
               '.'//char(0), &           ! ScratchDir
               FileType,Debug,VIsDouble);

      !Write the zone header information
      O=TecZne('SOULTION'//char(0), &
               1,nybox,nzbox,'BLOCK'//char(0),char(0));
      OO=nybox*nzbox

      allocate( OutArray(1,nybox,nzbox) )
      !--- Output velocity data ----
      do k = 1, nzbox
      do j = 1, nybox
        OutArray(1,j,k) = tke_calc(k,iref,j)
      end do
      end do
      O=TecDat(OO,OutArray(1,1,1),1);

      O=TecEnd();

      deallocate( OutArray )

      return
end subroutine writeins_tke_edge_2d

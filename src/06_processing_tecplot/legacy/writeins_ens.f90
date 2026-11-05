subroutine writeins_ens(imore)
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

      write(*,*) '    Writing the instantaneous flow field'
      !----------- TECPLOT FILE -----------
      outfname = './OUTPUT/ins_ens_'
      O=TecIni111('InsField'//char(0), & ! title
               'enstrophy'//char(0), &   ! variables
               trim(outfname)//cidx//'.plt'//char(0), & ! Fname
               '.'//char(0), &           ! ScratchDir
               FileType,Debug,VIsDouble);

      !Write the zone header information
      O=TecZne('SOULTION'//char(0), &
               nxbox,nybox,nzbox,'BLOCK'//char(0),char(0));
      OO=nxbox*nybox*nzbox

      allocate( OutArray(nxbox,nybox,nzbox) )
      !--- Output velocity data ----
      do k = 1, nzbox
      do j = 1, nybox
      do i = 1, nxbox
        OutArray(i,j,k) = ens(k,i,j)
      end do
      end do
      end do
      O=TecDat(OO,OutArray(1,1,1),1);

      O=TecEnd();

      deallocate( OutArray )
      return
end subroutine writeins_ens

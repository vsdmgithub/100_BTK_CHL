subroutine writexy(var,kref,filename,imore)
      use mgeomdata
      use msetup, only: nxbox,nybox
      implicit none
      double precision::var(m3,m1,m2)
      integer::i,j,imore,kref
      character( 7)::filename
      character(99)::outfname
      integer::TecIni111,TecDat,TecZne,TecEnd,TecNod,TecFil
      integer::O,OO,Debug,VIsDouble,FileType
      double precision, dimension(:,:,:),allocatable:: OutArray
      character(6)::cidx

      FileType = 2    ! 0 = full, 1 = grid, 2 = solution
      Debug = 0       ! 0 = no debug & 1 = debug
      VIsDouble = 0   ! 0 = single & 1 = double precision

      write(cidx,fmt='(i6.6)') imore

      !----------- TECPLOT FILE -----------
      outfname = './OUTPUT/ins_XY_'//trim(filename)//'_'//cidx//'.plt'
      write(*,*) 'write... '//trim(outfname)
      O=TecIni111('InsField'//char(0),  & ! title
               trim(filename)//char(0), & ! variables
               trim(outfname)//char(0), & ! Fname
               '.'//char(0), &            ! ScratchDir
               FileType,Debug,VIsDouble);

      !Write the zone header information
      O=TecZne('SOULTION'//char(0), &
               nxbox,nybox,1,'BLOCK'//char(0),char(0));
      OO=nxbox*nybox

      allocate( OutArray(nxbox,nybox,1) )
      !--- Output velocity data ----
      do j = 1, nybox
      do i = 1, nxbox
        OutArray(i,j,1) = var(kref,i,j)
      end do
      end do
      O=TecDat(OO,OutArray(1,1,1),1);

      O=TecEnd();

      deallocate( OutArray )

      return
end subroutine writexy

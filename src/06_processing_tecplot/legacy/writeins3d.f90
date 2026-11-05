subroutine writeins3d(dat,fname,varname,imore)
 use mgeomdata
 use msetup
 use mindices
 use mreadavg
 implicit none
 !double precision, dimension(nzboxp,nxboxp,nyboxp)::dat
 double precision, dimension(0:nzbox,0:nxbox,0:nybox)::dat
 integer::i,j,k,ii,jj,kk,imore
 character(70)::outfname
 !integer::TecIni111,TecDat,TecZne,TecEnd,TecNod,TecFil  ! MZW: Use TecIni
 integer::TecIni,TecDat,TecZne,TecEnd,TecNod,TecFil
 integer::O,OO,Debug,VIsDouble,FileType
 double precision, dimension(:,:,:),allocatable:: OutArray
 character(7)::cidx
 character(3)::fname
 character(1)::varname

 FileType = 2    ! 0 = full, 1 = grid, 2 = solution
 Debug = 0       ! 0 = no debug & 1 = debug
 VIsDouble = 0   ! 0 = single & 1 = double precision

 write(cidx,fmt='(i7.7)') imore

 write(*,*) '    Writing the instantaneous flow field'
 !----------- TECPLOT FILE -----------
outfname = trim(outdir)//fname//'_'

 O=TecIni('InsField'//char(0), &        ! title
          varname//char(0),       &        ! variables
          trim(outfname)//cidx//'.plt'//char(0), & ! Fname
          '.'//char(0), &                  ! ScratchDir
          FileType,Debug,VIsDouble);

 !Write the zone header information
 O=TecZne('FIELD'//char(0), &
          nxbox,nybox,nzbox,'BLOCK'//char(0),char(0));
 OO=nxbox*nybox*nzbox

 allocate( OutArray(nxbox,nybox,nzbox) )
 !--- Output velocity data ----
 do k = 1, nzbox
 do j = 1, nybox
 do i = 1, nxbox
   OutArray(i,j,k) = dat(k,i,j)
 end do
 end do
 end do
 O=TecDat(OO,OutArray(1,1,1),1);

 O=TecEnd();

 deallocate( OutArray )
 return
end subroutine writeins3d

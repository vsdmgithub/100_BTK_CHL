subroutine writeins3d_all(dat1,dat2,dat3,dat4,imore)
 use mgeomdata
 use msetup
 use mindices
 use mreadavg
 implicit none
 double precision, dimension(0:nzbox,0:nxbox,0:nybox)::dat1,dat2,dat3,dat4
 integer::i,j,k,ii,jj,kk,imore
 character(100)::outfname
 !integer::TecIni111,TecDat,TecZne,TecEnd,TecNod,TecFil  ! MZW: Use TecIni
 integer::TecIni,TecDat,TecZne,TecEnd,TecNod,TecFil
 integer::O,OO,Debug,VIsDouble,FileType
 double precision, dimension(:,:,:),allocatable:: Out1,Out2,Out3,Out4
 character(7)::cidx

 FileType = 2    ! 0 = full, 1 = grid, 2 = solution
 Debug = 0       ! 0 = no debug & 1 = debug
 VIsDouble = 0   ! 0 = single & 1 = double precision

 write(cidx,fmt='(i7.7)') imore

 write(*,*) '    Writing the instantaneous flow field'
 !----------- TECPLOT FILE -----------
 outfname = trim(outdir)//'sol_'//cidx//'.plt'

 O=TecIni('InsField'//char(0), &        ! title
          'u v w T'//char(0), &        ! variables
          outfname(1:LEN_TRIM(outfname))//char(0), & ! Fname
          '.'//char(0), &                  ! ScratchDir
          FileType,Debug,VIsDouble);

 !Write the zone header information
 O=TecZne('FIELD'//char(0), &
          nxbox,nybox,nzbox,'BLOCK'//char(0),char(0));
 OO=nxbox*nybox*nzbox

 allocate( Out1(nxbox,nybox,nzbox),Out2(nxbox,nybox,nzbox),Out3(nxbox,nybox,nzbox),Out4(nxbox,nybox,nzbox)  )
 !--- Output velocity data ----
 do k = 1, nzbox
 do j = 1, nybox
 do i = 1, nxbox
   Out1(i,j,k) = dat1(k,i,j)
   Out2(i,j,k) = dat2(k,i,j)
   Out3(i,j,k) = dat3(k,i,j)
   Out4(i,j,k) = dat4(k,i,j)
 end do
 end do
 end do
 O=TecDat(OO,Out1(1,1,1),1);
 O=TecDat(OO,Out2(1,1,1),1);
 O=TecDat(OO,Out3(1,1,1),1);
 O=TecDat(OO,Out4(1,1,1),1);
 O=TecEnd();
 deallocate( Out1 )
 deallocate( Out2 )
 deallocate( Out3 )
 deallocate( Out4 )
 return

end subroutine writeins3d_all

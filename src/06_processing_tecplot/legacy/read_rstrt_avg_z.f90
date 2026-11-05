!234567
  subroutine read_restart_avg_z(imore)
  use mvariable
  use msetup
  use mgeomdata
  implicit none
  integer(int_p)::i,j,k,imore,length,kk
  real(real_p)::dx_tmp,dy_tmp
  character(7)::cidx
  character(50)::cfnm
  real(real_p),dimension(:)  ,allocatable::uz,vz,wz,cf
  real(real_p),dimension(:,:),allocatable::um_stag,vm_stag,wm_stag,um,vm,wm
  integer(8)::read_count
!--Tecplot format--
  character(70)::outfname
  integer::TecIni111,TecDat,TecZne,TecEnd,TecNod,TecFil
  integer::O,OO,Debug,VIsDouble,FileType
  double precision, dimension(:,:,:),allocatable:: OutArray

!JL  allocate(   ut(0:ngz  ,0:ngx+1,0:ngy  ) )
!JL  allocate(   vt(0:ngz  ,0:ngx  ,0:ngy+1) )
!JL  allocate(   wt(0:ngz+1,0:ngx  ,0:ngy  ) )
  allocate( uz(0:ngz  ) )
  allocate( vz(0:ngz  ) )
  allocate( wz(0:ngz+1) )

  allocate( um_stag(0:ngx+1,0:ngy  ) )
  allocate( vm_stag(0:ngx  ,0:ngy+1) )
  allocate( wm_stag(0:ngx  ,0:ngy  ) )
     
  write(cidx,fmt='(i7.7)') imore

  !---------------------------------------------------------
  !       Read in velocity field
  !---------------------------------------------------------
  inquire(iolength=length) uz(1)
  open (18,file=trim(cfn)//'ucvcwc.eigr.'//cidx,form="unformatted", &
           access="direct",action='read',status='old',recl=length)
    read_count = 0
    ! U-velocity
    do j = 0, ngy
    do i = 0, ngx+1
      do k = 0, ngz
        read_count = read_count + 1
        read (18,rec=read_count) uz(k)
      end do
      um_stag(i,j) = sum( uz(1:ngzm) ) / real( ngzm, real_p )
    end do
    end do
    write(*,*) '   * Reading u-component... done'
  
    ! V-velocity
    do j = 0, ngy+1
    do i = 0, ngx
      do k = 0, ngz
        read_count = read_count + 1
        read (18,rec=read_count) vz(k)
      end do
      vm_stag(i,j) = sum( vz(1:ngzm) ) / real( ngzm, real_p )
    end do
    end do
    write(*,*) '   * Reading v-component... done'
  
    ! W-velocity
    do j = 0, ngy
    do i = 0, ngx
      do k = 0, ngz+1
        read_count = read_count + 1
        read (18,rec=read_count) wz(k)
      end do
      wm_stag(i,j) = sum( wz(1:ngzm) ) / real( ngzm, real_p )
    end do
    end do
    write(*,*) '   * Reading w-component... done'
    close(18)

    deallocate( uz, vz, wz )

  !---------------------------------------------------------
  !       Compute simple statistics
  !---------------------------------------------------------
  ! collocate to cell center
    allocate( um(ngxm, ngym) )
    allocate( vm(ngxm, ngym) )
    allocate( wm(ngxm, ngym) )
 
    do j = 1, ngym
    do i = 1, ngxm
     um(i,j) = 0.5 * ( um_stag(i,j) + um_stag(i+1,j) )
     vm(i,j) = 0.5 * ( vm_stag(i,j) + vm_stag(i,j+1) )
     wm(i,j) = 0.5 * ( wm_stag(i,j) + wm_stag(i,j  ) )
    end do
    end do
  
    deallocate( um_stag, vm_stag, wm_stag )
  
  ! compute wall-shear stress
    allocate( cf(ngxm) )
    do i = 1, ngxm
      dx_tmp = xpp(i,1)-xpv(i,1)
      dy_tmp = ypp(i,1)-ypv(i,1)

      cf(i) = ( um(i,1) / sqrt(dx_tmp**2 + dy_tmp**2) ) / Re
    end do


  !---------------------------------------------------------
  !       Output result
  !---------------------------------------------------------
    FileType = 0    ! 0 = full, 1 = grid, 2 = solution
    Debug = 0       ! 0 = no debug & 1 = debug
    VIsDouble = 0   ! 0 = single & 1 = double precision
   
    !----------- TECPLOT FILE -----------
    outfname = trim(outdir)//'avg_singleK_'
    !outfname = './OUTPUT/avg_singleK_'
    O=TecIni111('AvgField'//char(0), & ! title
             'x y um vm wm cf'//char(0), &          ! variables
             trim(outfname)//cidx//'.plt'//char(0), & ! Fname
             '.'//char(0), &           ! ScratchDir
             FileType,Debug,VIsDouble);
   
    !Write the zone header information
    O=TecZne('FIELD'//char(0), &
             ngxm,ngym,1,'BLOCK'//char(0),char(0));
    OO=ngxm * ngym
   
    allocate( OutArray(ngxm,ngym,1) )
    !--- Output grid node locations -----
    do j = 1, ngym
    do i = 1, ngxm 
      OutArray(i,j,1) = xpp(i,j)
    end do
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
    do j = 1, ngym
    do i = 1, ngxm 
      OutArray(i,j,1) = ypp(i,j)
    end do
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
   
    !--- Output velocity data ----
    do j = 1, ngym
    do i = 1, ngxm 
      OutArray(i,j,1) = um(i,j)
    end do
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
    do j = 1, ngym
    do i = 1, ngxm 
      OutArray(i,j,1) = vm(i,j)
    end do
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
    do j = 1, ngym
    do i = 1, ngxm 
      OutArray(i,j,1) = wm(i,j)
    end do
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
    do i = 1, ngxm 
      OutArray(i,:,1) = cf(i)
    end do
    O=TecDat(OO,OutArray(1,1,1),1);
   
    O=TecEnd();
   
    deallocate( OutArray )


  !---------------------------------------------------------
  !       Finalize
  !---------------------------------------------------------
    deallocate( cf )
    deallocate( um, vm, wm )

  return
  end subroutine

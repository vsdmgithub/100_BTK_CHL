!=======================================================================
! Adjusting time step to control CFL number
!
! CFLcontrol_init()
! CFLcontrol_read()
! CFLcontrol_write()
! CFLcontrol_adjust(timeStep)
! CFLcontrol_maxCFL()*
! CFLcontrol_CFLcomponents(c)
! CFLcontrol_status(heading, values)
! CFLcontrol_info(iunit)
!

module CFLcontrol
	use data_export
	use mesh_export

	integer iCFL          ! whether to run at constant CFL number
	real    CFLlimit      ! maximum allowable CFL number
	real    dtInitial     ! initial deltaT that was specified

	real    CFLmax        ! maximum CFL number
	real    dtMax         ! maximum time step
	real    dtMin         ! minimum time step
	integer nReduce       ! number of times deltaT reduced
	integer nIncrease     ! number of times deltaT increased

	real    CFL           ! current CFL number
	integer initFlag

end module

!=======================================================================
subroutine CFLcontrol_init()
	use CFLcontrol

	call readInt("iCFL", iCFL)
	call readFloat("CFLlimit", CFLlimit)
	if (iCFL == 0) then
		call readFloat("dt", dt)
		dtInitial = dt
		nReduce = 0
		nIncrease = 0
	else
		dt = 1.e-6
	end if

	!T I added the next two lines
	call readFloat("dt", dt)
	dtInitial = dt

	initFlag = 1
	CFLmax = 0.
	dtMin = 0.
	dtMax = 0.

	return
end

subroutine CFLcontrol_read()
	use CFLcontrol
	if (iCFL == 1) call readFloat("dt", dt)
	return
end

subroutine CFLcontrol_write()
	use CFLcontrol
	call writeInt("iCFL", iCFL)
	call writeFloat("CFLlimit", CFLlimit)
	if (iCFL == 0) then
		call writeFloat("dt", dtInitial)
	else
		call writeFloat("dt", dt)
	end if
	return
end

!=======================================================================
subroutine CFLcontrol_adjust()
	use CFLcontrol

	! Calculate maximum CFL number
	call CFLcontrol_maxCFL()
	if (CFL > CFLmax) CFLmax = CFL

	! Calculate the next time-step dt based on CFL number
	ratio = CFLlimit/CFL
	if (iCFL == 0) then
		! Run at fixed dt, if possible
		if (ratio < 1.) then
			! CFL too large; reduce it
			dt = dt*ratio
			nReduce = nReduce + 1
		else if (dt .ne. dtInitial) then
			! CFL ok; return to initial dt, but slowly
			! Increase by 10%
			dt = min(dt*1.1, dtInitial)
			nIncrease = nIncrease +1
		end if
	else
		! Run at fixed CFL
		dt = dt*ratio
	end if

	dtMin = min(dtMin, dt)
	dtMax = max(dtMax, dt)
	timeStep = dt

	if (initFlag == 1) then
		initFlag = 0
		CFLmax = 0.
		dtMin = dt
		dtMax = dt
	end if

	return
end

subroutine CFLcontrol_maxCFL()
!cccccccccccccccccccccccccccccccccccccccccccccc
! calculate CFL number if dt is specified
!cccccccccccccccccccccccccccccccccccccccccccccc
	use CFLcontrol

	real*8  :: vvmax, velmax, wcmax, wmax
	real*8  :: uavg, vavg, wavg
	real    :: Atmp(3)
	integer :: i,j,k, ii,jj
	real	:: LocalCFL
	integer :: iCFLmax, jCFLmax, kCFLmax

	real    :: Btmp(7)
	real	:: tmpVal
	real	:: LocalCFLx(7)
	integer :: iCFLmaxX(7)  , jCFLmaxX(7)  , kCFLmaxX(7)

	vvmax =-1.e9
	velmax=-1.e6
	vcmax =-1.e6
	wcmax =-1.e6
	vmax  =-1.e6
	wmax  =-1.e6

	tmpVal     =-1.e9
	LocalCFLx  =-1.e9

	do j=j1_T,j2_T   
	do i=i1_T,i2_T   
	do k=1,ngz-1
		ii = ibmap_1(i)
		jj = jbmap_1(j)
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!       vvmax=amax1(CFLmax,vvmax)
		!       CFLmax=(abs(uc(i,j,k))/(xpg(i+1,j)-xpg(i,j))+   &
		!               abs(vc(i,j,k))/(ypg(i,j+1)-ypg(i,j))+   &
		!               abs(wc(i,j,k))/(zpg(k+1)-zpg(k)))
		!       uavg=(uc(i+1,j,k)+uc(i,j,k))*.5
		!       vavg=(vc(i,j+1,k)+vc(i,j,k))*.5
		!       wavg=(wc(i,j,k+1)+wc(i,j,k))*.5
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		! vvmax=max(abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))+ &
  		! 	abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))+ &
  		! 	abs(wc(k,i,j))/(zpg(k+1)  -zpg(k)),vvmax)

		LocalCFL=abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))+ &
  			 abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))+ &
  			 abs(wc(k,i,j))/(zpg(k+1)  -zpg(k))

		if (LocalCFL.gt.vvmax) then
			vvmax = LocalCFL;
			iCFLmax = ii ;  jCFLmax = jj ;  kCFLmax = k ;
		end if

!JL		!---- Debug components -----
!JL		tmpVal = dt *	abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))
!JL		if (tmpVal.gt.LocalCFLx(1)  ) then
!JL			LocalCFLx(1)   = tmpVal
!JL			iCFLmaxX(1)   = ii ;  jCFLmaxX(1)   = jj ;  kCFLmaxX(1)   = k ;
!JL		end if
!JL
!JL		tmpVal = dt *	abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))
!JL		if (tmpVal.gt.LocalCFLx(2)  ) then
!JL			LocalCFLx(2)   = tmpVal
!JL			iCFLmaxX(2)   = ii ;  jCFLmaxX(2)   = jj ;  kCFLmaxX(2)   = k ;
!JL		end if
!JL
!JL		tmpVal = dt *	abs(wc(k,i,j))/(zpg(k+1)  -zpg(k))
!JL		if (tmpVal.gt.LocalCFLx(3)  ) then
!JL			LocalCFLx(3)   = tmpVal
!JL			iCFLmaxX(3)   = ii ;  jCFLmaxX(3)   = jj ;  kCFLmaxX(3)   = k ;
!JL		end if
!JL
!JL		tmpVal = dt * ( abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))+ &
!JL				abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))	)
!JL                if (tmpVal.gt.LocalCFLx(4) ) then
!JL			LocalCFLx(4)  = tmpVal
!JL			iCFLmaxX(4)  = ii ;  jCFLmaxX(4)  = jj ;  kCFLmaxX(4)  = k ;
!JL		end if
!JL
!JL		tmpVal = dt * ( abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))+ &
!JL				abs(wc(k,i,j))/(zpg(k+1)  -zpg(k))		)
!JL		if (tmpVal.gt.LocalCFLx(5) ) then
!JL			LocalCFLx(5)  = tmpVal
!JL			iCFLmaxX(5)  = ii ;  jCFLmaxX(5)  = jj ;  kCFLmaxX(5)  = k ;
!JL		end if
!JL
!JL		tmpVal = dt * ( abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))+ &
!JL				abs(wc(k,i,j))/(zpg(k+1)  -zpg(k))		)
!JL                if (tmpVal.gt.LocalCFLx(6) ) then
!JL			LocalCFLx(6)  = tmpVal
!JL			iCFLmaxX(6)  = ii ;  jCFLmaxX(6)  = jj ;  kCFLmaxX(6)  = k ;
!JL		end if
!JL
!JL		tmpVal = dt * ( abs(uc(k,i,j))/(xpg(ii+1,jj)-xpg(ii,jj))+ &
!JL				abs(vc(k,i,j))/(ypg(ii,jj+1)-ypg(ii,jj))+ &
!JL				abs(wc(k,i,j))/(zpg(k+1)  -zpg(k))		)
!JL                if (tmpVal.gt.LocalCFLx(7)) then
!JL			LocalCFLx(7) = tmpVal
!JL			iCFLmaxX(7) = ii ;  jCFLmaxX(7) = jj ;  kCFLmaxX(7) = k ;
!JL		end if

		uavg=(uc(k,i+1,j)+uc(k,i,j))/2.
		vavg=(vc(k,i,j+1)+vc(k,i,j))/2.
                wavg=(wc(k+1,i,j)+wc(k,i,j))/2.
		velmax=max(sqrt(uavg**2+vavg**2+wavg**2),velmax)

		vmax =max( v(k,i,j),vmax ) !TJ: determine vmax
		vcmax=max(vc(k,i,j),vcmax) !TJ: determine vcmax
                wmax =max( w(k,i,j),wmax )
		wcmax=max(wc(k,i,j),wcmax) !TJ: determine wcmax

	end do
	end do
	end do

	CFL=dt*vvmax

	!T Global CFL maximum
	Atmp(1) = CFL
	Atmp(2) = velmax
	Atmp(3) = wmax

	call globalMax(Atmp, 3)
	call globalMax(vcmax,1)
	call globalMax(vmax ,1)	!SJL: added to get correct global max
	call globalMax(wcmax,1)	!SJL: added to get correct global max

	!if (abs(CFL-Atmp(1))/CFL .le. 1e-10) then
        !SJL: allows CFL to have zero without error
        if (abs(CFL-Atmp(1)) .le. 1e-10*abs(CFL)) then
		print '(a,f12.6,3(a,i5))', 'Maximum CFL = '	, CFL, &
		                                      ' @ '	, iCFLmax ,',', jCFLmax ,',', kCFLmax
		print '(a,3f12.6)'      , 'CFL x,y,z   = '	, dt*abs( uc(kCFLmax,imap_1(iCFLmax),jmap_1(jCFLmax))/(xpg(iCFLmax,jCFLmax)-xpg(iCFLmax+1,jCFLmax)) )  &
								, dt*abs( vc(kCFLmax,imap_1(iCFLmax),jmap_1(jCFLmax))/(ypg(iCFLmax,jCFLmax)-ypg(iCFLmax,jCFLmax+1)) )  &
								, dt*abs( wc(kCFLmax,imap_1(iCFLmax),jmap_1(jCFLmax))/(zL/ngz) )
	end if


	CFL    = Atmp(1)
	velmax = Atmp(2)
	wmax   = Atmp(3)

	!JL if (irank.eq.1)  write(*,6) ntime,stime,dt,CFL,velmax
	if (irank.eq.1) then !JL151026 begin
                write(*,6) ntime,stime,dt,CFL,velmax  
                write(*,*) ' '
        end if               !JL151026 end
 6      format('n=',i7,'  stime=',1pe10.3,'  dt=',1pe10.3, &
     		'   CFL=',1pe10.3,'  velmax=',1pe10.3)

	!if (irank.eq.1)  write(6,*) 'vcmax=   ',vcmax
	!if (irank.eq.1)  write(6,*) 'wcmax=   ',wcmax
	!if (irank.eq.1)  write(6,*) 'vmax=   ',vmax
	!if (irank.eq.1)  write(6,*) 'wmax=   ',wmax
        if (irank.eq.1)  write(6,'(a7,es14.7,a9,es14.7)') 'vcmax= ', vcmax, ', wcmax= ', wcmax !JL151026
        if (irank.eq.1)  write(6,'(a7,es14.7,a9,es14.7)') 'v_max= ', vmax , ', w_max= ', wmax  !JL151026


	Btmp = LocalCFLx

	call globalMax(Btmp, 7)
!	do i=1,7
!		! if (abs(LocalCFLx(i)  -Btmp(i))/Btmp(i) .le. 1e-10) then
!                !SJL: allows Btmp(i) to have zero without error
!		if (abs(LocalCFLx(i)  -Btmp(i)) .le. 1e-10*abs(Btmp(i))) then
!			if (i==1) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max X   CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==2) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max  Y  CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==3) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max   Z CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==4) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max XY  CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==5) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max X Z CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==6) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max  YZ CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!			if (i==7) print '(a,i3,f12.6,a,3(x,i6,a))'    , 'Max XYZ CFL = '  , i, LocalCFLx(i), '@ (',iCFLmaxX(i),',', jCFLmaxX(i) ,',', kCFLmaxX(i),')'
!		end if
!	end do


	return
end


subroutine CFLcontrol_CFLcomponents(c)
	use CFLcontrol
	real c(6)

	! Components of CFL number
	c = 0.
	do k=k1_,k2_
	do j=j1_,j2_
	do i=i1_,i2_
		c(1) = max(c(1), abs(U(i,j,k))*dxi )
		c(2) = max(c(2), abs(V(i,j,k))*dymi(j) )
		c(3) = max(c(3), abs(W(i,j,k))*dzi )
		c(4) = max(c(4), abs(P(i,j,k))*4.*dxi**2 )
		c(5) = max(c(5), abs(P(i,j,k))*4.*dymi(j)**2 )
		c(6) = max(c(6), abs(P(i,j,k))*4.*dzi**2 )
	end do
	end do
	end do

	c = dt*c
	call globalMax(c, 6)

	return
end

!=======================================================================
subroutine CFLcontrol_status(heading, values)
	use CFLcontrol
	character*(*) heading, values

	write (heading, 11) "dt", "CFL"
	write (values, 12) dt, CFL

11  format (8a12)
12  format (8f12.5)

	return
end

subroutine CFLcontrol_info(iunit)
	use CFLcontrol

	call sectionTitle(iunit, "CFL control")
	write (iunit, 11) dtMin, dtMax, CFLmax
	if (iCFL == 0) then
		write (iunit, 12) nReduce, nIncrease
	end if

11  format (4x, "Minimum time step            : ", f11.4 / &
            4x, "Maximum time step            : ", f11.4 / &
	        4x, "Maximum CFL number           : ", f11.4)
12  format (4x, "Number of times dt reduced   : ", i11 / &
	        4x, "Number of times dt increased : ", i11)

	return
end


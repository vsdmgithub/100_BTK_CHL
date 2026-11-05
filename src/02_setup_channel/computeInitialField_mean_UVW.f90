subroutine ComputeInitialFieldStaggeredUVW_mean()
	use initialField
	use analyticField

	! specify input and out (dummy) variables
	!integer, intent(in) :: ngx, ngy, ngz
	!real, dimension(0:ngy), intent(in) :: U_b
	!real, dimension(0:ngz,0:ngx+1,0:ngy), intent(out) :: U_bmean

	!--------------------------------------------------------------
	!--- U ---
	do j=1,niym
	do i=1,nix
	do k=1,nizm
		!--- compute 3-D perturbation field ---
		U(k,i,j) = Ub(j)
	end do
	end do
	end do

	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!--- V ---
	!do i=0,ngx
	!do j=0,ngy+1
	!do k=0,ngz
	    
	!    !--- argument of complex exponential function ---
	!    xz_V = alp*x_vw(i) + beta*z_uv(k);	
	!    
	!    !--- compute 3-D perturbation field ---
	!     V(k,i,j) = eps*real( Voss(j)*exp(zi*xz_V) );
	!    !V(i,j,k) = eps*abs( Voss(k)*exp(zi*xy_V) );
	!    ! V(i,j,k) = -1.0*eps*imag(Voss(k)*sin(beta*y_v(j))*exp(zi*alp*x_v(i)));  
	!end do
	!end do
	!end do
    
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!--- W ---
	!do i=0,ngx
	!do j=0,ngy
	!do k=0,ngz+1

	!    %--- argument of complex exponential function ---
	!    xz_W = alp*x_vw(i) + beta*z_w(k);	
	!    
	!    %--- compute 3-D perturbation field ---
	!     W(k,i,j) = eps*real( Woss(j)*exp(zi*xz_W) );
	!     %W(i,j,k) = eps*abs( Woss(k)*exp(zi*xy_W) );
	!    % W(i,j,k) = eps*real(Woss(k)*cos(beta*y_w(j))*exp(zi*alp*x_w(i)));  
	!     
	!end do
	!end do
	!end do
    
	!--- Navier-slip at lower wall ---
	U(:,:, 0 ) = -U(:,:,  1 )

	!--- Navier-slip at upper wall ---
	U(:,:,niy) = -U(:,:,niym)

	!--- x-periodicity ---
	U(:,  1  ,:) = U(:,nix ,:)
	U(:,  0  ,:) = U(:,nixm,:)
	U(:,nix+1,:) = U(:,  2 ,:)

	!--- z-periodicity ---
	U(  0  ,:,:) = U(nizm,:,:)
	U( niz ,:,:) = U(  1 ,:,:)

	return
end

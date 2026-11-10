!234567
module mparameter
implicit none
include 'PARAM.H'
include 'param.inc'
end module mparameter
!===============================================
module msetup
use mparameter
implicit none
integer(int_p)::inssn,insen,insitv
integer(int_p)::wgrid,iwrite3d,iwrite2d,isingle
real(real_p)::Re
character(50)::cfn
character(50)::outdir
end module msetup
!===============================================
! Fields are stored in Tecplot order (i,j,k) = (x,y,z) so they can be
! handed to TecDat without any transposition copy.
module mvariable
use mparameter
implicit none
real(real_p),allocatable,dimension(:,:,:):: u,v,w
! Staggered work arrays, solver order (k,i,j), allocated once
real(real_p),allocatable,dimension(:,:,:):: ut,vt,wt
end module mvariable
!===============================================
module mgeomdata
use mparameter
implicit none
real(real_p):: xpg(ngx,ngy),ypg(ngx,ngy)
real(real_p):: dz,Lz
end module mgeomdata
!===============================================

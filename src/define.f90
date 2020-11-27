module definvar
    implicit none
	!$ integer :: num_threads 
    integer*8 :: nmo, nprims, nEDFprims, ncenter, charge, ix, iy, iz
    integer :: ifinish, lastcen, nx, ny, nz, max_edr_exponents=50
    real*8 :: ep, rr, sftx, sfty, sftz, sftx2, sfty2, sftz2, expterm
	real*8 :: dx, dy, dz, orgx, orgy, orgz
    real*8 :: expcutoff=-40D0, aug3D = 6.0D0
    real*8, allocatable :: MOocc(:), MOene(:), MOco(:,:)
	real*8,parameter :: pi=3.141592653589793D0
	character(len=7) :: cubname
	 
    type atoms
        character(len=2) :: name
        real*8 :: xcoo,ycoo,zcoo,charg 
    end type atoms 
    type(atoms), allocatable :: a(:)

    type prims
        integer*8 :: center, type
        real*8 :: expo
    end type prims
    type(prims), allocatable :: b(:)
	
	type EDFun
		integer*8 :: EDFcenter, EDFtype
        real*8 :: EDFexpo, EDFco
	end type EDFun
	type(EDFun), allocatable :: EDF(:)
	

    integer :: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,&
	&0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
    integer :: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,& 
	&0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
    integer :: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,& 
	&5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)

	!Convert shell type to the number of basis functions in the shell: 0=s,1=p,-1=sp,2=6d,-2=5d,3=10f,-3=7f,4=15g,-4=9g,5=21h,-5=11h
	integer :: shtype2nbas(-5:5)=(/ 11,9,7,5,4,1,3,6,10,15,21 /)     
    
end module definvar
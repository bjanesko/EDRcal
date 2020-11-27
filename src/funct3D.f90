! ======================================
! =========== 3D Grid Setting ==========
! ======================================
subroutine setgrid()
	use definvar
	implicit none
	integer :: grdqal
	integer*8 :: ntotal
	real*8 :: molxlen, molylen, molzlen
	
	write(*,*) " Please select the quality of grid"
	write(*,*) " 1 Low    ; about 125000 points in total to cover the whole system"
	write(*,*) " 2 Medium ; about 512000 points in total to cover the whole system"
	write(*,*) " 3 High   ; about 1728000 points in total to cover the whole system"
	write(*,*) " 4 Manually enter the number of points to cover the whole system"
	read(*,*) grdqal
	
	if (grdqal==1) then
		ntotal = 125000
	else if (grdqal==2) then 
		ntotal = 512000
	else if (grdqal==3) then
		ntotal = 1728000
	else if (grdqal==4) then
		write(*,*) " Enter the number of points:"
		read(*,*) ntotal
	else 
		write(*,*) " Option not supported. Select option 1 to 4"
		read(*,*) grdqal
		if (grdqal==4) then
			write(*,*) " Enter the number of points:"
			read(*,*) ntotal
		end if
	end if
	
	
	molxlen = (maxval(a%xcoo)-minval(a%xcoo))+2*aug3D
	molylen = (maxval(a%ycoo)-minval(a%ycoo))+2*aug3D
	molzlen = (maxval(a%zcoo)-minval(a%zcoo))+2*aug3D
	
	if (molxlen==0.0d0) then
		molxlen = 3.0d0
	else if (molylen==0.0d0) then
		molylen = 3.0d0
	else if (molzlen==0.0d0) then
		molzlen = 3.0d0
	end if
	
	dx = (molxlen*molylen*molzlen/dfloat(ntotal))**(1.0d0/3.0d0)
	dy = dx
	dz = dx
	
	nx = nint(molxlen/dx)+1
	ny = nint(molylen/dy)+1 
	nz = nint(molzlen/dz)+1

	orgx = minval(a%xcoo)-aug3D
	orgy = minval(a%ycoo)-aug3D
	orgz = minval(a%zcoo)-aug3D
	
end subroutine setgrid
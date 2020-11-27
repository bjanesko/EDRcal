program main
    use definvar
	use functions
	USE OMP_LIB
    implicit none
	integer*8 :: orbstart, orbend
	integer :: i,j,k,prop, nedr
	integer :: wtime1,wtime2
	real*8, allocatable :: cubmat(:,:,:), xarr(:), yarr(:), zarr(:)
	real*8 :: dedr, edrastart, edrainc 
	character :: c80tmp*80, cubfilname*200, propname*200
	
	! Prologue
	write(*,*)
	write(*,"(1x,A)") "EDRcal v1.0 (November 27, 2020)"
	write(*,"(1x,A)") "https://janeskoresearchgroup.wordpress.com/"
	write(*,"(1x,A)") "Please read the following work for details:"
	write(*,"(1x,A)") "Janesko, B. G.; Scalmani, G.; Frisch, M. J., J. Chem. Phys. 2014, 141, 144104."
	write(*,"(1x,A)") "Mehmood, A.; Janesko, B. G., Angew. Chem. Int. Ed. 2017, 56, 6878-6881."
	write(*,"(1x,A)") "Janesko, B. G.; et. al., J. Chem. Theory Comput. 2016, 12, 3185-3194."
	write(*,"(1x,A)") "Mehmood, A.; Jones, S. I.; Tao, P.; Janesko, B. G., J. Chem. Inf. Model. 2018, 58, 1836-1846."
	write(*,*)
	! End Prologue 
	
	!$ num_threads = omp_get_max_threads()
	!$ write(*,"((A),i5)") " Number of threads being used:", num_threads
    
	call inpnamtyp()
	
	write(*,*) " Select property: "
	write(*,*) " 1 Total electron density "
	write(*,*) " 2 Electron delocalization range function EDR(r;d) "
	write(*,*) " 3 Orbital overlap distance function D(r) "
	read(*,*) prop
	
	if (prop==2) then
		call EDR_interface(2,dedr,nedr,edrastart,edrainc)
		write(*,*)
	else if (prop==3) then
		call EDR_interface(3,dedr,nedr,edrastart,edrainc)
		write(*,*)
	end if
	
	ifinish=0
	
	call setgrid()
	
	allocate(xarr(nx))
	allocate(yarr(ny))
	allocate(zarr(nz))
	
	do k=1,nz
		write(c80tmp,"(D20.13)") orgz+(k-1)*dz
		read(c80tmp,*) zarr(k)
	end do
	do j=1,ny
		write(c80tmp,"(D20.13)") orgy+(j-1)*dy
		read(c80tmp,*) yarr(j)
	end do
	do i=1,nx
		write(c80tmp,"(D20.13)") orgx+(i-1)*dx
		read(c80tmp,*) xarr(i)
	end do
	
	allocate(cubmat(nx,ny,nz))
	cubmat = 0d0

    !$ call omp_set_num_threads(num_threads)
	
	if (prop==1) then
		orbstart = 1
		orbend = nmo
		call walltime(wtime1)
		write(*,*)
		!$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k) SCHEDULE(dynamic) NUM_THREADS(num_threads)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					cubmat(i,j,k)=rho(orbstart,orbend,xarr(i),yarr(j),zarr(k))				
				end do
			end do	
			ifinish=ifinish+1
			call showprog(ifinish,nz)
		end do
		!$OMP END PARALLEL DO
	
	else if (prop == 2) then
		orbstart = 1
		orbend = nmo
		call walltime(wtime1)
		write(*,*)
		!$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k) SCHEDULE(dynamic) NUM_THREADS(num_threads)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					cubmat(i,j,k)=edr(orbstart,orbend,dedr,xarr(i),yarr(j),zarr(k))					
				end do
			end do	
			ifinish=ifinish+1
			call showprog(ifinish,nz)
		end do
		!$OMP END PARALLEL DO
	
	else if (prop == 3) then
		orbstart = 1
		orbend = nmo
		call walltime(wtime1)
		write(*,*)
		!$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k) SCHEDULE(dynamic) NUM_THREADS(num_threads)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					cubmat(i,j,k)=dr(orbstart,orbend,nedr,edrastart,edrainc,xarr(i),yarr(j),zarr(k))				
				end do
			end do
			ifinish=ifinish+1
			call showprog(ifinish,nz)
		end do
		!$OMP END PARALLEL DO	
	end if
	
	deallocate(xarr)
	deallocate(yarr)
	deallocate(zarr)	
	
	if (prop==1) then
		write(cubfilname,"('rho.cub')")
		write(propname,"('Total density')")
		
	else if (prop==2) then
		write(cubfilname,"('EDR.cub')")
		write(propname,"('EDR(r;d)')")
		
	else if (prop==3) then
		write(cubfilname,"('Dr.cub')")
		write(propname,"('D(r)')")
	end if
	
	write(*,*)
	write(*,"(2x,'Writing ',(a),' data to a cube file. Please wait...')") trim(propname)
	
	call writcube(cubfilname,propname,cubmat,nx,ny,nz,dx,dy,dz,orgx,orgy,orgz)
	deallocate(cubmat)
	
	write(*,"(2x,'Cube file ',(a),' is generated in current folder.')") trim(cubfilname)
	
	call walltime(wtime2)
	write(*,*)
	write(*,"(2x,'The process took ',i6,' seconds of wall clock time.')") wtime2-wtime1
end program main

! ======================================
! Subroutine : EDR and Dr user interface
! ======================================
subroutine EDR_interface(prop,dedr,nedr,edrastart,edrainc)
	use definvar
	implicit none
	integer, intent(in) ::  prop
	integer, intent(out) :: nedr
	real*8, intent(out) :: dedr, edrastart, edrainc
	integer :: edrmaxpara, wrtnumedr
	real*8 :: wrtstart, wrtexpo(max_edr_exponents)

 if (prop==2) then 
	write(*,*) "Input length scale d (Bohr)   e.g. 0.85"
	read(*,*) dedr
else if (prop==3) then 
	write(*,*) "1 Manually input total number, start and increment in EDR exponents"
	write(*,*) "2 Use default values   i.e. 20,2.50,1.50"
	read(*,*) edrmaxpara
	if (edrmaxpara==1) then  
		write(*,*) "Please input in order: exponents start increment   e.g. 20 2.5 1.5"
		write(*,*) "Note: Max. allowed exponents are 50 and min. allowed increment is 1.01"
		read(*,*) nedr,edrastart,edrainc
		if (nedr<1) then
			write(*,*) "Error: Bad Number of EDR exponents. Should be between 1 to 50"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		else if (nedr>50) then
			write(*,*) "Error: Bad Number of EDR exponents. Should be between 1 to 50"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		end if
		if (edrainc<1.01d0) then
			write(*,*) "Error: Bad increment in EDR exponents. Should not be less than 1.01"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		end if
	else if (edrmaxpara==2) then
		nedr=20
		edrastart=2.5D0
		edrainc=1.5D0
	end if
	write(*,*) "The following EDR exponents will be used in calculation:"
	wrtstart=edrastart
	do wrtnumedr=1,nedr
		wrtexpo(wrtnumedr)=wrtstart
		wrtstart=wrtstart/edrainc
		write(*,"(E13.5)") wrtexpo(wrtnumedr) 
	end do
end if
end subroutine
! ======================================
! =====  Subroutine : Progress bar =====
! ======================================
subroutine showprog(inow,nall)
	implicit none
	integer :: inow,nall,iprog
	integer :: itmp=0
	character :: c80tmp*80
	
	iprog=int(dfloat(inow)/nall*50)
	c80tmp=' Progress: ['
	c80tmp(13:62)=repeat('+',iprog)
	c80tmp(13+iprog:62)=repeat('-',50-iprog)
	c80tmp(63:63)=']'
	
	write(c80tmp(64:),"(f8.2,' %')") dfloat(inow)/nall*100
	
	itmp=itmp+1
	
	if (itmp==1) c80tmp(79:79)='-'
	if (itmp==2) c80tmp(79:79)='\'
	if (itmp==3) c80tmp(79:79)='|'
	if (itmp==4) then
		c80tmp(79:79)='/'
		itmp=0
	end if
	write(*,"(2a$)") trim(c80tmp),char(13)
	if (inow>=nall) write(*,*)
end subroutine
! ======================================
! ======= Subroutine : Timestamp =======
! ======================================
subroutine walltime(now)
	implicit none
	character :: nowdate*20,nowtime*20
	integer :: now, nowhour, nowminute, nowsecond
	
	call date_and_time(nowdate,nowtime)
	
	read(nowtime(1:2),*) nowhour
	read(nowtime(3:4),*) nowminute
	read(nowtime(5:6),*) nowsecond
	
	now = nowhour*3600 + nowminute*60 + nowsecond
end subroutine
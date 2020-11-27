! ======================================
! =====  Read Input File Name Type  ====
! ======================================
subroutine inpnamtyp()
	use definvar
	implicit none
	integer :: filen
	character :: inpfile*200
	
	write(*,*)
	write(*,*) " Write the name of input file."
	write(*,*) " Supported file types are .fch, .fchk, .wfn, and .wfx"
	write(*,*)
	
	do while(.true.)
		read(*,"(a)") inpfile
		write(*,*)
		filen = len_trim(inpfile)
	
		if (inpfile(1:1) == '"' .OR. inpfile(1:1) == "'") inpfile(1:1) = " " 					     ! Remove first " or ' character from file name or path
		if (inpfile(filen:filen) == '"' .OR. inpfile(filen:filen) == "'") inpfile(filen:filen) = " " ! Remove last " or ' character from file name or path
	
		inpfile=adjustl(inpfile)

		if (index(inpfile, '.fch')/=0 .OR. index(inpfile, '.FCH')/=0 .OR. &
			index(inpfile, '.fchk')/=0 .OR. index(inpfile, '.FCHK')/=0 .OR. &
			index(inpfile, '.FChk')/=0 .OR. index(inpfile, '.FCHk')/=0) then
			write(*,*) " The input file type is .fch"
			write(*,*) " Loading the input file..."
			call read_fchk(inpfile)
			exit
		else if (index(inpfile, '.WFN')/=0 .OR. index(inpfile, '.wfn')/=0) then
			write(*,*) " The input file type is .wfn"
			write(*,*) " Loading the input file..."
			call read_wfn(inpfile)
			exit
		else if (index(inpfile, '.WFX')/=0 .OR. index(inpfile, '.wfx')/=0) then 
			write(*,*) " The input file type is .wfx"
			write(*,*) " Loading the input file..."
			call read_wfx(inpfile)
			exit
		else
			write(*,*) " The input file type is not supported. Try again."
			cycle
		end if
	end do
	write(*,*)

end subroutine inpnamtyp
! ======================================
! ===========  Read .wfn File  =========
! ======================================
subroutine read_wfn(filename) 
	use definvar
	implicit none
	character :: title*80, c80tmp*80, filename*200
	real*8 :: totalene, virlrato
	integer :: eqlsign1, eqlsign2,i,j,k

	open(1,file=filename,access="sequential",status="old")
	
	read(1,"(A)") title
	read(1,"(a8,i15,13x,i7,11x,i9,7x)") c80tmp, nmo, nprims, ncenter
	
	if (index(c80tmp,"SLATER")/=0) then
		write(*,"(1x,A)") " The type of orbitals in the input .wfn file are Slater type, which are currently not supported by EDRgen."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop
	end if

	allocate(a(ncenter))
	allocate(b(nprims))
	allocate(MOocc(nmo))
	allocate(MOene(nmo))
	allocate(MOco(nmo,nprims))

	do i = 1, ncenter
		read(1,"(a24,3f12.8,10x,f5.1)") c80tmp,a(i)%xcoo,a(i)%ycoo,a(i)%zcoo,a(i)%charg
		a(i)%name = c80tmp(3:4)
	end do

	read(1,"(20x,20i3)") (b(i)%center, i = 1,nprims)
	read(1,"(20x,20i3)") (b(i)%type, i = 1,nprims)
	read(1,"(10x,5E14.7)") (b(i)%expo, i = 1,nprims)

	do i = 1, nmo
		read(1,"(A)") c80tmp
		do j = 1, 80
			if (c80tmp(j:j)=='=') then
				read(c80tmp(j+1:),*) MOocc(i)
				exit
			end if
		end do
		do j = 80, 1, -1
			if (c80tmp(j:j)=='=') then
				read(c80tmp(j+1:),*) MOene(i)
				exit
			end if
		end do

		read(1,"(5f16.8)") (MOco(i,j), j = 1,nprims) 
	end do

	read(1,*)
	
	eqlsign1 = 0
	eqlsign2 = 0
	read(1,"(A80)") c80tmp
	do i = 1, 80
		if(c80tmp(i:i)=='=') then
			eqlsign1 = i
			exit
		end if
	end do 
	do i = 80, 1, -1
		if(c80tmp(i:i)=='=') then
			eqlsign2 = i
			exit
		end if
	end do
	if(eqlsign1/=0)	read(c80tmp(eqlsign1+1:),*) totalene
	if(eqlsign1==0) write(*,*) " Warning: Total Energy is not located in input .wfn file...!"
	if(eqlsign2/=0) read(c80tmp(eqlsign2+1:),*) virlrato
	if(eqlsign2==0) write(*,*) " Warning: Virial ratio is not located in input .wfn file...!"
	
	write (*,"(2x,A)") "File loaded successfully. "
		
end subroutine read_wfn
! ======================================
! ===========  Read .wfx File  =========
! ======================================
subroutine read_wfx(filnam)
	use definvar
	implicit none
	character(Len=*),intent(in) :: filnam
	integer :: i,j,k,l,ierror, found, tempfnd
	real*8 :: totalene, virlrato
	integer*8 :: ncorele 
	
	if (allocated(a)) deallocate(a)
	if (allocated(b)) deallocate(b)
	if (allocated(MOocc)) deallocate(MOocc)
	if (allocated(MOene)) deallocate(MOene)
	if (allocated(MOco)) deallocate(MOco)
	if (allocated(EDF)) deallocate(EDF)
	
	open(1,file=filnam,status='old',access='sequential',action='read',iostat=ierror)
	
	call strlocat(1,"<Number of Nuclei>",found)
	if (found==1) then
		read(1,*)
		read(1,*)  ncenter
	else if (found==0) then
		write(*,"(1x,A)") " Information about number of atoms not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 
	end if

	allocate(a(ncenter))
	
	call strlocat(1,"<Atomic Numbers>",found)
	if (found==1) then
		read(1,*)
		do i = 1, ncenter
			read(1,*)  a(i)%charg
		end do
	else if (found==0) then
		write(*,"(1x,A)") " Warning: Atomic numbers of atoms not found."
	end if
	
	call strlocat(1,"<Nuclear Cartesian Coordinates>",found)
	if (found==1) then
		read(1,*)
		do i = 1, ncenter
			read(1,*)  a(i)%xcoo,a(i)%ycoo,a(i)%zcoo
		end do
	else if (found==0) then
		write(*,"(1x,A)") " Nuclear coordinates not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 
	end if
	
	call strlocat(1,"<Number of Primitives>",found)
	if (found==1) then
		read(1,*)
		read(1,*)  nprims
	else if (found==0) then
		write(*,"(1x,A)") " Information about number of primitives not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 
	end if
	
	allocate(b(nprims))

	call strlocat(1,"<Primitive Centers>",found)
	read(1,*)
	read(1,*)  b%center
	
	call strlocat(1,"<Primitive Types>",found)
	read(1,*)
	read(1,*)  b%type

	call strlocat(1,"<Primitive Exponents>",found)
	read(1,*)
	read(1,*)  b%expo

	call strlocat(1,"<Number of Occupied Molecular Orbitals>",found)
	read(1,*)
	read(1,*)  nmo
	
	allocate(MOocc(nmo))
	allocate(MOene(nmo))
	allocate(MOco(nmo,nprims))
	
	call strlocat(1,"<Molecular Orbital Occupation Numbers>",found)
	read(1,*)	
	do i = 1, nmo
		read(1,*) MOocc(i)	
	end do
	
	call strlocat(1,"<Molecular Orbital Energies>",found)	
	read(1,*)
	do i = 1, nmo
		read(1,*) MOene(i)	
	end do
	
	call strlocat(1,"<Molecular Orbital Primitive Coefficients>",found)	
	read(1,*)
	do i = 1, nmo
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*) (MOco(i,j), j= 1,nprims)	
	end do

	call strlocat(1,"<Number of EDF Primitives>",found)
	if (found==1) then
		read(1,*)
		read(1,*) nEDFprims
		
		call strlocat(1,"<Number of Core Electrons>",tempfnd)
		read(1,*)
		read(1,*) ncorele
		
		write(*,*) " Note: The input file contains electron density functions (EDF)"
		write(*,"(1x,a,i7,a,i6,a)") " Loading ",nEDFprims, " EDFs for ",ncorele," core electrons."
		
		allocate(EDF(nEDFprims))
		
		call strlocat(1,"<EDF Primitive Centers>",tempfnd)
		read(1,*)
		read(1,*) EDF%EDFcenter

		call strlocat(1,"<EDF Primitive Types>",tempfnd)
		read(1,*)
		read(1,*) EDF%EDFtype
		if (maxval(EDF%EDFtype) .GT. 1) then
			write(*,*) " GTFs electron density function is not a S type"
			write(*,*) " GTFs must be S type"
			write(*,*) " Press ENTER button to exit."
			read(*,*)
			stop  
		end if
		
		call strlocat(1,"<EDF Primitive Exponents>",tempfnd)
		read(1,*)
		read(1,*) EDF%EDFexpo

		call strlocat(1,"<EDF Primitive Coefficients>",tempfnd)
		read(1,*)
		read(1,*) EDF%EDFco
	end if
	
	call strlocat(1,"<Energy = T + Vne + Vee + Vnn>",found)
	if (found==1) then
		read(1,*)
		read(1,*)  totalene
	else if (found==0) then
		write(*,"(1x,A)") " Warning: Total Energy is not located in input .wfn file...!"
	end if
	
	call strlocat(1,"<Virial Ratio (-V/T)>",found)
	if (found==1) then
		read(1,*)
		read(1,*)  virlrato
	else if (found==0) then
		write(*,"(1x,A)") " Warning: Virial ratio is not located in input .wfn file...!"
	end if
	
	write (*,"(2x,A)") "File loaded successfully. "	
	
end subroutine read_wfx
! ======================================
! ===========  Read .fchk File  =========
! ======================================
subroutine read_fchk(filnam)
	use definvar
	implicit none
	character(Len=*), intent(in) :: filnam
	integer :: i,j,k,l,ierror,found, wfntype, isspherical, ipos5D, ipos6D, ish, ishtyp5D, ishtyp6D, numshorb5D, numshorb6D
	integer :: iexp, imo, ibasis, s2f(-5:5,21)=0
	integer*8 :: nbasis, nbasis5D, nbasisCar, nindbasis, tempread, nshell, nprimshell
	integer*8, allocatable :: shelltype(:), shelltype5D(:), shelltype6D(:), shellcon(:), shell2atom(:), MOtype(:)
	integer, allocatable :: basshell(:), bascen(:), bastype(:), basstart(:), basend(:), primstart(:), primend(:) 
	real*8, allocatable :: primexp(:), concoeff(:), SPconcoeff(:), aMOcoeff(:,:), bMOcoeff(:,:), primconnorm(:)
	real*8, allocatable :: CObasa(:,:),CObasb(:,:), CObasa5D(:,:), CObasb5D(:,:)
	real*8 :: nelec, nalpha, nbeta, normgau, tnormgau, temp
	real*8 :: conv5d6d(6,5), conv7f10f(10,7), conv9g15g(15,9), conv11h21h(21,11)
	character :: c80tmp*80
	s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
	s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
	s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
	s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
	s2f(-1,1:4)=(/ 1,2,3,4 /)
	s2f(0,1)=1
	s2f(1,1:3)=(/ 2,3,4 /)
	s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
	s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
	s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
	s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)

	call gensphcartab(1,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

	if (allocated(a)) deallocate(a)
	if (allocated(b)) deallocate(b)
	if (allocated(MOocc)) deallocate(MOocc)
	if (allocated(MOene)) deallocate(MOene)
	if (allocated(MOco)) deallocate(MOco)
	
	open(1,file=filnam,status='old',access='sequential',action='read',iostat=ierror)
	
	read(1,*)
	read(1,"(A)") c80tmp
	
	if (c80tmp(11:11)=='R') then
		wfntype = 0
	else if (c80tmp(11:11)=='U') then
		wfntype = 1
	else if (c80tmp(11:11)=='RO') then
		wfntype = 2
		if (c80tmp(13:13)=='3') then
			wfntype = 0			!RO3LYP
		end if
	end if
	
	call strlocat(1,"Number of electrons",found)
	if (found==1) then
		read(1,"(49x,f12.0)") nelec
		read(1,"(49x,f12.0)") nalpha
		read(1,"(49x,f12.0)") nbeta
	else if (found==0) then
		write(*,"(1x,A)") " Information about number of electrons not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 
	end if

	call strlocat(1,"Number of basis functions",found)
	if (found==1) then
		read(1,"(45x,i16)") nbasis
		nindbasis = nbasis
	else if (found==0) then
		write(*,"(1x,A)") " Information about number of basis not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 	
	end if
	
	call strlocat(1,"Number of independent functions",found)
	if (found==1) then
		read(1,"(45x,i16)") nindbasis
	end if

	call strlocat(1,"Atomic numbers",found)
	if (found==1) then
		read(1,"(50x,i11)") ncenter
		allocate(a(ncenter))
		read(1,*) a%charg
	else if (found==0) then
		write(*,"(1x,A)") " Information about atoms not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 	
	end if
	
	call strlocat(1,"Current cartesian coordinates",found)
	if (found==1) then
		read(1,*) 
		read(1,"(5(1PE16.8))") (a(i)%xcoo, a(i)%ycoo, a(i)%zcoo, i=1,ncenter)
	else if (found==0) then
		write(*,"(1x,A)") " Atomic Coordinates not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 	
	end if
	
	call strlocat(1,"Shell types",found)
	if (found==1) then
		read(1,"(50x,i11)") nshell
		allocate(shelltype(nshell))
		read(1,*) shelltype
	else if (found==0) then
		write(*,"(1x,A)") " Information about shells not found."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 	
	end if
	
	isspherical=0
	if (any(shelltype<=-2)) isspherical=1
	if (any(abs(shelltype)>5)) then
		write(*,"(1x,A)") " Error: Input file contains GTFs with angular moment higher than h, which are not supported."
		write(*,"(1x,A)") " Press ENTER to exit:"
		read(*,*)
		stop 	
	end if

	allocate(shellcon(nshell))
	allocate(shell2atom(nshell))
	
	call strlocat(1,"Number of primitives per shell",found)
	if (found==1) then
		read(1,*)
		read(1,*) shellcon
	end if

	call strlocat(1,"Shell to atom map",found)
	if (found==1) then
		read(1,*)
		read(1,*) shell2atom
	end if

	call strlocat(1,"Primitive exponents",found)
	if (found==1) then
		read(1,"(50x,i11)") nprimshell
		allocate(primexp(nprimshell))
		read(1,*) primexp
	end if
	
	call strlocat(1,"Contraction coefficients",found)
	if (found==1) then
		read(1,*)
		allocate(concoeff(nprimshell))
		read(1,*) concoeff
		read(1,"(A)") c80tmp
		if (index(c80tmp,"P(S=P) Contraction coefficients")/=0) then
			allocate(SPconcoeff(nprimshell))
			read(1,*) SPconcoeff
		end if
	end if
	
	if (wfntype==0 .OR. wfntype==2) then ! R or RO
		nmo = nbasis
		allocate(MOocc(nmo))
		allocate(MOene(nmo))
		allocate(MOtype(nmo))
		allocate(aMOcoeff(nmo,nbasis))
		
		MOocc = 0d0
		MOene = 0d0
		MOtype = 0
		aMOcoeff = 0d0
		
		if (wfntype==0) then
			MOocc(1:nint(nelec/2))=2.0D0
		else if (wfntype==2) then
			MOocc(1:int(nbeta))=2.0D0
			MOocc(int(nbeta)+1:int(nalpha))=1.0D0
			MOtype(int(nbeta)+1:int(nalpha))=1
		end if
		
		call strlocat(1,"Alpha Orbital Energies",found)
		if (found==0) call strlocat(1,"orbital energies",found)
		read(1,*)
		read(1,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
		
		call strlocat(1,"Alpha MO coefficients",found)
		if (found==0) call strlocat(1,"MO coefficients",found)
		read(1,*)
		read(1,"(5(1PE16.8))") ((aMOcoeff(i,j),j=1,nbasis),i=1,nindbasis)	
	
	else if (wfntype==1) then
		nmo = 2*nbasis
		allocate(MOocc(nmo))
		allocate(MOene(nmo))
		allocate(MOtype(nmo))
		allocate(aMOcoeff(nmo,nbasis))
		allocate(bMOcoeff(nmo,nbasis))
		
		MOocc = 0d0
		MOene = 0d0
		MOtype = 0
		aMOcoeff = 0d0
		bMOcoeff = 0d0
		
		MOocc(1:int(nalpha))=1.0D0
		MOocc(nbasis+1:nbasis+int(nbeta))=1.0D0
		MOtype(1:nbasis) = 1
		MOtype(nbasis+1:nmo) = 2
		
		call strlocat(1,"Alpha Orbital Energies",found)
		if (found==0) call strlocat(1,"alpha orbital energies",found)
		read(1,*)
		read(1,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
		
		call strlocat(1,"Beta Orbital Energies",found)
		if (found==0) call strlocat(1,"beta orbital energies",found)
		if (found==0) then
			write(*,"(1x,A)") " Beta orbital energies not found."
			write(*,"(1x,A)") " Press ENTER to exit:"
			read(*,*)
			stop 	
		end if
		read(1,*)
		read(1,"(5(1PE16.8))") (MOene(i),i=nbasis+1,nbasis+nindbasis)
		
		call strlocat(1,"Alpha MO coefficients",found)
		if (found==0) call strlocat(1,"alpha MO coefficients",found)
		read(1,*)
		read(1,"(5(1PE16.8))") ((aMOcoeff(i,j),j=1,nbasis),i=1,nindbasis)
		
		call strlocat(1,"Beta MO coefficients",found)
		if (found==0) call strlocat(1,"beta MO coefficients",found)
		read(1,*)
		read(1,"(5(1PE16.8))") ((bMOcoeff(i,j),j=1,nbasis),i=1,nindbasis)
	end if
	
	if (isspherical==1) then
		allocate(shelltype5D(nshell))
		shelltype5D = shelltype
		where (shelltype<=-2) shelltype = -shelltype
		nbasis5D = nbasis
		nbasis = 0
		do i = 1, nshell
			nbasis = nbasis + shtype2nbas(shelltype(i))
		end do
	end if
	
	allocate(shelltype6D(nshell))
	shelltype6D = shelltype
	nbasisCar = nbasis
	
	nprims=0
	do i = 1,nshell
		nprims = nprims + shtype2nbas(shelltype(i))*shellcon(i)
	end do
	
	allocate(MOco(nmo,nprims))
	allocate(b(nprims))
	allocate(basshell(nbasis))
	allocate(bascen(nbasis))
	allocate(bastype(nbasis))
	allocate(primstart(nbasis))
	allocate(primend(nbasis))
	allocate(primconnorm(nprims))
	allocate(basstart(ncenter))
	allocate(basend(ncenter))
	
	if (isspherical==0) then 
		allocate(CObasa(nbasis,nbasis))
		CObasa = transpose(aMOcoeff)
		if (wfntype==1) then
			allocate(CObasb(nbasis,nbasis))
			CObasb = transpose(bMOcoeff)
		end if
	
	else if (isspherical==1) then
		allocate(CObasa(nbasis,nbasis))
		allocate(CObasa5D(nbasis5D,nbasis5D))
		CObasa5D = transpose(aMOcoeff)
		CObasa = 0
		if (wfntype==1) then
			allocate(CObasb(nbasis,nbasis))
			allocate(CObasb5D(nbasis5D,nbasis5D))
			CObasb5D = transpose(bMOcoeff)
			CObasb = 0
		end if
		
		ipos5D = 1
		ipos6D = 1
		do ish = 1, nshell
			ishtyp5D = shelltype5D(ish)
			ishtyp6D = shelltype(ish)
			numshorb5D = shtype2nbas(ishtyp5D)
			numshorb6D = shtype2nbas(ishtyp6D)
			if (ishtyp5D>=-1) then 
				CObasa(ipos6D:ipos6D + numshorb6D-1,1:nbasis5D) = CObasa5D(ipos5D:ipos5D + numshorb5D-1,:)
				if (wfntype==1) CObasb(ipos6D:ipos6D + numshorb6D-1,1:nbasis5D) = CObasb5D(ipos5D:ipos5D + numshorb5D-1,:)			
			else if (ishtyp5D==-2) then
				CObasa(ipos6D:ipos6D + 5,1:nbasis5D) = matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
				if (wfntype==1) CObasb(ipos6D:ipos6D+5,1:nbasis5D) = matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
			else if (ishtyp5D==-3) then
				CObasa(ipos6D:ipos6D+9,1:nbasis5D) = matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
				if (wfntype==1) CObasb(ipos6D:ipos6D+9,1:nbasis5D) = matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
			else if (ishtyp5D==-4) then
				CObasa(ipos6D:ipos6D+14,1:nbasis5D) = matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
				if (wfntype==1) CObasb(ipos6D:ipos6D+14,1:nbasis5D) = matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
			else if (ishtyp5D==-5) then
				CObasa(ipos6D:ipos6D+20,1:nbasis5D) = matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
				if (wfntype==1) CObasb(ipos6D:ipos6D+20,1:nbasis5D) = matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
			end if
			ipos5D = ipos5D + numshorb5D
			ipos6D = ipos6D + numshorb6D
		end do
	end if
	
	k = 1 
	iexp = 1
	ibasis = 1
	do i = 1,nshell
		b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center = shell2atom(i)
		basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1) = i
		bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1) = shell2atom(i)
		do j = 1,shtype2nbas(shelltype(i))
			b(k:k+shellcon(i)-1)%type = s2f(shelltype(i),j)
			bastype(ibasis) = s2f(shelltype(i),j)
			primstart(ibasis) = k 
			primend(ibasis) = k + shellcon(i)-1
			do l = 1, shellcon(i)
				b(k)%expo = primexp(iexp+l-1)
				tnormgau = normgau(b(k)%type,b(k)%expo)
				temp = concoeff(iexp+l-1)
				if (shelltype(i)==-1.and.j/=1) temp = SPconcoeff(iexp+l-1)
				primconnorm(k) = temp*tnormgau
				do imo = 1, nmo
					if (wfntype==0.or.wfntype==2) then
						MOco(imo,k) = CObasa(ibasis,imo)*temp*tnormgau
					else if (wfntype==1) then
						if (isspherical==1) then
							if (imo<=nbasis5D) MOco(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
							if (imo>nbasis5D) MOco(imo,k)=CObasb(ibasis,imo-nbasis5D)*temp*tnormgau
						else
							if (imo<=nbasis) MOco(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
							if (imo>nbasis) MOco(imo,k)=CObasb(ibasis,imo-nbasis)*temp*tnormgau
						end if
					end if
				end do
				k = k + 1
			end do
			ibasis = ibasis + 1
		end do
		iexp = iexp + shellcon(i)
	end do
	
	write (*,"(2x,A)") "File loaded successfully. "

end subroutine
! ======================================
! ===========  Write Cube File =========
! ======================================

subroutine writcube(filname,propname,cubmat,nox,noy,noz,dix,diy,diz,orgcx,orgcy,orgcz)
	use definvar
	implicit none
	integer :: nox,noy,noz
	integer :: i,j,k
	character :: filname*200, propname*200
	real*8 :: cubmat(nx,ny,nz),dix,diy,diz,orgcx,orgcy,orgcz
	
	open(2,file=filname,status="replace")
	
	write(2,"(a)") " Generated by EDRcal v1.0"
	write(2,"(1x,a,a,i8,a)") trim(propname)," projected on",nx*ny*nz, " grid points"
	write(2,"(i5,3f12.6)") ncenter, orgcx, orgcy, orgcz
	write(2,"(i5,3f12.6)") nox, dix, 0.0, 0.0
	write(2,"(i5,3f12.6)") noy, 0.0, diy, 0.0
	write(2,"(i5,3f12.6)") noz, 0.0, 0.0, diz
	
	do i = 1, ncenter
		write(2,"(i5,4f12.6)") int(a(i)%charg), a(i)%charg, a(i)%xcoo, a(i)%ycoo, a(i)%zcoo
	end do
	
	where (abs(cubmat)<=1D-99) cubmat=0D0 
	do i = 1, nox
		do j = 1, noy
			write(2,"(6(1PE14.5E3))",advance="no") cubmat(i,j,1:noz)
			write(2,*)
		end do
	end do
	close(2)
end subroutine writcube
! ======================================
! ===========  Search String  ==========
! ======================================
subroutine strlocat(filename,string,found)
	implicit none
	integer, intent(in) :: filename
	character(Len=*), intent (in) :: string
	integer, intent(out) :: found
	integer :: ierror
	character :: c200tmp*200
	
	rewind(filename)
	
	do while(.TRUE.)  
		read(filename,"(a)",iostat=ierror) c200tmp
		if (index(c200tmp,string)/=0) then
			backspace(filename)
			found = 1
			return
		end if
		if (ierror/=0) exit
	end do
	
	found = 0

end subroutine strlocat

! ======================================
! ===  Generate Spherical Harmonic  ====
! ======================================
subroutine gensphcartab(iprog,matd,matf,matg,math)
	implicit none
	real*8 :: matd(6,5), matf(10,7), matg(15,9), math(21,11)
	integer :: iprog
	matd = 0D0
	matf = 0D0
	matg = 0D0
	math = 0D0

	matd(1:3,1) = (/ -0.5D0,-0.5D0,1D0 /)
	matd(5,2) = 1D0
	matd(6,3) = 1D0
	matd(1:2,4) = (/ sqrt(3D0)/2D0,-sqrt(3D0)/2D0 /)
	matd(4,5) = 1D0

	matf(3,1)=1D0
	matf(6,1)=-1.5D0/sqrt(5D0)
	matf(9,1)=-1.5D0/sqrt(5D0)
	matf(1,2)=-sqrt(3D0/8D0)
	matf(4,2)=-sqrt(3D0/40D0)
	matf(7,2)=sqrt(6D0/5D0)
	matf(2,3)=-sqrt(3D0/8D0)
	matf(5,3)=-sqrt(3D0/40D0)
	matf(8,3)=sqrt(6D0/5D0)
	matf(6,4)=sqrt(3D0)/2D0
	matf(9,4)=-sqrt(3D0)/2D0
	matf(10,5)=1D0
	matf(1,6)=sqrt(5D0/8D0)
	matf(4,6)=-3D0/sqrt(8D0)
	matf(2,7)=-sqrt(5D0/8D0)
	matf(5,7)=3D0/sqrt(8D0)

	if (iprog==1) then !for .fch
		matg(1,1)=1D0
		matg(3,1)=-3D0*sqrt(3D0/35D0)
		matg(5,1)=3D0/8D0
		matg(10,1)=-3D0*sqrt(3D0/35D0)
		matg(12,1)=3D0/4D0*sqrt(3D0/35D0)
		matg(15,1)=3D0/8D0
		matg(6,2)=2D0*sqrt(5D0/14D0)
		matg(8,2)=-1.5D0/sqrt(14D0)
		matg(13,2)=-1.5D0*sqrt(5D0/14D0)
		matg(2,3)=2D0*sqrt(5D0/14D0)
		matg(4,3)=-1.5D0*sqrt(5D0/14D0)
		matg(11,3)=-1.5D0/sqrt(14D0)
		matg(3,4)=-3D0*sqrt(3D0/28D0)
		matg(5,4)=sqrt(5D0)/4D0
		matg(10,4)=3D0*sqrt(3D0/28D0)
		matg(15,4)=-sqrt(5D0)/4D0
		matg(7,5)=3D0/sqrt(7D0)
		matg(9,5)=-sqrt(5D0/28D0)
		matg(14,5)=-sqrt(5D0/28D0)
		matg(8,6)=-3D0/sqrt(8D0)
		matg(13,6)=sqrt(5D0/8D0)
		matg(4,7)=-sqrt(5D0/8D0)
		matg(11,7)=3D0/sqrt(8D0)
		matg(5,8)=sqrt(35D0)/8D0
		matg(12,8)=-3D0/4D0*sqrt(3D0)
		matg(15,8)=sqrt(35D0)/8D0
		matg(9,9)=-sqrt(5D0)/2D0
		matg(14,9)=sqrt(5D0)/2D0
	else if (iprog==2) then !For .molden
		matg(3,1)=1D0
		matg(1,1)=3D0/8D0
		matg(2,1)=3D0/8D0
		matg(11,1)=-3D0*sqrt(3D0/35D0)
		matg(12,1)=-3D0*sqrt(3D0/35D0)
		matg(10,1)=3D0/4D0*sqrt(3D0/35D0)
		matg(8,2)=2D0*sqrt(5D0/14D0)
		matg(5,2)=-1.5D0*sqrt(5D0/14D0)
		matg(14,2)=-1.5D0/sqrt(14D0)
		matg(9,3)=2D0*sqrt(5D0/14D0)
		matg(7,3)=-1.5D0*sqrt(5D0/14D0)
		matg(13,3)=-1.5D0/sqrt(14D0)
		matg(11,4)=3D0*sqrt(3D0/28D0)
		matg(12,4)=-3D0*sqrt(3D0/28D0)
		matg(1,4)=-sqrt(5D0)/4D0
		matg(2,4)=sqrt(5D0)/4D0
		matg(15,5)=3D0/sqrt(7D0)
		matg(4,5)=-sqrt(5D0/28D0)
		matg(6,5)=-sqrt(5D0/28D0)
		matg(5,6)=sqrt(5D0/8D0)
		matg(14,6)=-3D0/sqrt(8D0)
		matg(7,7)=-sqrt(5D0/8D0)
		matg(13,7)=3D0/sqrt(8D0)
		matg(1,8)=sqrt(35D0)/8D0
		matg(2,8)=sqrt(35D0)/8D0
		matg(10,8)=-3D0/4D0*sqrt(3D0)
		matg(4,9)=sqrt(5D0)/2D0
		matg(6,9)=-sqrt(5D0)/2D0
	end if

	math(1,1)=1D0
	math(12,1)=-5D0/sqrt(21D0)
	math(3,1)=-5D0/sqrt(21D0)
	math(19,1)=5D0/8D0
	math(5,1)=5D0/8D0
	math(14,1)=sqrt(15D0/7D0)/4D0
	math(7,2)=sqrt(5D0/3D0)
	math(16,2)=-3D0*sqrt(5D0/28D0)
	math(9,2)=-3D0/sqrt(28D0)
	math(21,2)=sqrt(15D0)/8D0
	math(11,2)=sqrt(5D0/3D0)/8D0
	math(18,2)=sqrt(5D0/7D0)/4D0
	math(2,3)=sqrt(5D0/3D0)
	math(4,3)=-3D0*sqrt(5D0/28D0)
	math(13,3)=-3D0/sqrt(28D0)
	math(6,3)=sqrt(15D0)/8D0
	math(20,3)=sqrt(5D0/3D0)/8D0
	math(15,3)=sqrt(5D0/7D0)/4D0
	math(12,4)=sqrt(5D0)/2D0
	math(3,4)=-sqrt(5D0)/2D0
	math(19,4)=-sqrt(35D0/3D0)/4D0
	math(5,4)=sqrt(35D0/3D0)/4D0
	math(8,5)=sqrt(5D0/3D0)
	math(17,5)=-sqrt(5D0/12D0)
	math(10,5)=-sqrt(5D0/12D0)
	math(16,6)=sqrt(5D0/6D0)
	math(9,6)=-sqrt(1.5D0)
	math(21,6)=-sqrt(17.5D0)/8D0
	math(11,6)=sqrt(17.5D0)/8D0
	math(18,6)=sqrt(5D0/6D0)/4D0	
	math(4,7)=-sqrt(5D0/6D0)
	math(13,7)=sqrt(1.5D0)
	math(20,7)=-sqrt(17.5D0)/8D0
	math(6,7)=sqrt(17.5D0)/8D0
	math(15,7)=-sqrt(5D0/6D0)/4D0	
	math(19,8)=sqrt(35D0)/8D0
	math(5,8)=sqrt(35D0)/8D0
	math(14,8)=-0.75D0*sqrt(3D0)	
	math(17,9)=sqrt(5D0)/2D0
	math(10,9)=-sqrt(5D0)/2D0	
	math(21,10)=3D0/8D0*sqrt(3.5D0)
	math(11,10)=5D0/8D0*sqrt(3.5D0)
	math(18,10)=-1.25D0*sqrt(1.5D0)	
	math(6,11)=3D0/8D0*sqrt(3.5D0)
	math(20,11)=5D0/8D0*sqrt(3.5D0)
	math(15,11)=-1.25D0*sqrt(1.5D0)
end subroutine
! ======================================
! ========  Calculate Factorial  =======
! ======================================
integer function ft(i)
	implicit none
	integer :: i, j
	
	ft = i
	if (i==0) ft = 1
	do j = i-1, 1, -1
		ft = ft * j
	end do
	
end function
! ======================================
! Normalize Coefficient of Cartesian GTF 
! ======================================

real*8 function normgau(itype,exp)
	use definvar
	implicit none
	integer*8 :: itype
	integer :: ft, iix, iiy, iiz
	real*8 :: exp
	
	iix=type2ix(itype)
	iiy=type2iy(itype)
	iiz=type2iz(itype)
	normgau = (2*exp/pi)**0.75D0*dsqrt((8*exp)**(iix+iiy+iiz)*ft(iix)*ft(iiy)*ft(iiz)/(ft(2*iix)*ft(2*iiy)*ft(2*iiz)))
	
end function

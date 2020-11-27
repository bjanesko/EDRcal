module functions
    use definvar
    implicit none
        
contains
! ======================================
! == Function: Total Electron Density ==
! ======================================
function rho(orbst,orbend,x,y,z)
	use definvar
	implicit none
	integer*8 :: orbs, orbst, orbend
	real*8 :: x, y, z, rho, wfnval(nmo), EDFrho
	
	call orbderv(orbst,orbend,x,y,z,wfnval)
	
	rho = 0d0
	do orbs = orbst,orbend
		rho = rho + MOocc(orbs)*wfnval(orbs)**2
	end do
	
	if (allocated(EDF)) then
		call EDFderv(x,y,z,EDFrho)
		rho = rho + EDFrho
	end if
	
end function
! ======================================
! = Subroutine: Total Electron Density =
! ======================================
subroutine orbderv(orbstrt,orbend,x,y,z,wfnval)
    use definvar
    implicit none
    real*8 :: x,y,z
    integer*8 :: orbstrt,orbend
	integer :: i,j,k, l
    real*8 :: wfnval(nmo)
    real*8 :: gtfval

    lastcen = -1
    wfnval = 0d0
    do j = 1, nprims
        ix = type2ix(b(j)%type)
        iy = type2iy(b(j)%type)
        iz = type2iz(b(j)%type)
        ep = b(j)%expo
        
        if (b(j)%center/=lastcen) then
          sftx = x-a(b(j)%center)%xcoo
          sfty = y-a(b(j)%center)%ycoo
          sftz = z-a(b(j)%center)%zcoo
          sftx2 = sftx*sftx
          sfty2 = sfty*sfty
          sftz2 = sftz*sftz
          rr = sftx2 +  sfty2 + sftz2  
        end if
        
        if (expcutoff>0 .OR. -ep*rr>expcutoff) then
            expterm = exp(-ep*rr)
        else 
            expterm = 0d0
        end if 
        
        lastcen = b(j)%center 

        if (expterm==0d0) cycle
        
        if (b(j)%type == 1) then
            gtfval = expterm
        else if (b(j)%type == 2) then
            gtfval = sftx*expterm
        else if (b(j)%type == 3) then
            gtfval = sfty*expterm
        else if (b(j)%type == 4) then
            gtfval = sftz*expterm 
        else if (b(j)%type == 5) then
            gtfval = sftx2*expterm
		else if (b(j)%type == 6) then
			gtfval = sfty2*expterm
		else if (b(j)%type == 7) then
			gtfval = sftz2*expterm
		else if (b(j)%type == 8) then
			gtfval = sftx*sfty*expterm
		else if (b(j)%type == 9) then
			gtfval = sftx*sftz*expterm
		else if (b(j)%type == 10) then
			gtfval = sfty*sftz*expterm
		else
			gtfval = sftx**ix * sfty**iy * sftz**iz * expterm
        end if

		do i = orbstrt,orbend
			wfnval(i) = wfnval(i) + MOco(i,j)*gtfval
		end do
    end do

end subroutine orbderv

! ======================================
! Subroutine: Electron Density Functions
! ======================================
subroutine EDFderv(x,y,z,EDFrho)
    use definvar
    implicit none
    real*8, intent(in) :: x,y,z
    real*8, intent(out) :: EDFrho
	integer :: i,j,k, l	
	real*8 :: EDFx,EDFy,EDFz,EDFx2,EDFy2,EDFz2,EDFrr,EDFep,EDFexpterm
	
	EDFrho = 0d0
	
	do i = 1, nEDFprims
		EDFx = x - a(EDF(i)%EDFcenter)%xcoo
		EDFy = y - a(EDF(i)%EDFcenter)%ycoo
		EDFz = z - a(EDF(i)%EDFcenter)%zcoo
		EDFx2 = EDFx * EDFx
		EDFy2 = EDFy * EDFy
		EDFz2 = EDFz * EDFz
		EDFrr = EDFx2 + EDFy2 + EDFz2
		EDFep = EDF(i)%EDFexpo
		EDFexpterm = exp(-EDFep*EDFrr)
		EDFrho = EDFrho + (EDF(i)%EDFco)*EDFexpterm
	end do
end subroutine EDFderv
! ======================================
! =========  Function: EDR(r;d) ========
! ======================================
function edr(orbst,orbend,dedr,x,y,z)
	implicit none
	real*8 :: edr,ed(max_edr_exponents),edrval(max_edr_exponents)
	real*8 :: dedr,edrdmaxval,x,y,z
	integer*8 :: orbst,orbend
	integer :: nedr
	nedr=1
	ed(1)=dedr**(-2.0d0)
	call EDRcal(2,orbst,orbend,x,y,z,nedr,ed,edrval,edrdmaxval)
	edr=edrval(1)
end function

! ======================================
! =========    Function: D(r)  =========
! ======================================

function dr(orbst,orbend,nedr,edrastart,edrainc,x,y,z)
	implicit none
	real*8 :: ed(max_edr_exponents),edrdmaxval,edrval(max_edr_exponents)
	real*8 :: edrexponent, dr, x, y, z, edrastart, edrainc
	integer*8 :: orbst,orbend
	integer :: iedr, nedr
	edrexponent = edrastart  
	do iedr=1,nedr
		ed(iedr)=edrexponent
		edrexponent=edrexponent/edrainc
	end do
	call EDRcal(3,orbst,orbend,x,y,z,nedr,ed,edrval,edrdmaxval)
	dr=edrdmaxval
end function

! ======================================
! = Subroutine: Main EDR(r;d) and D(r) =
! ======================================
subroutine EDRcal(prop,orbst,orbend,x,y,z,nedr,ed,edrval,edrdmaxval)
	use definvar
	implicit none
	real*8, intent(in) :: x,y,z,ed(max_edr_exponents)
	integer*8, intent(in) :: orbst, orbend
	integer, intent(in) :: nedr, prop
	real*8, intent(out) :: edrval(max_edr_exponents), edrdmaxval
	real*8 :: dens, psi(nmo), Bint(nmo,max_edr_exponents), amu0(max_edr_exponents)  
	real*8 :: xamu(3,max_edr_exponents), AMUVal(max_edr_exponents), sftval, gtfval
	real*8 :: dmaxdummy, edmax, shftx,shfty,shftz,shftx2,shfty2,shftz2,sqr,epp,exterm,gtval
	integer :: i, j, k, l, iedr, ixyz, ival, prevcen
	integer*8 :: tx, ty, tz
	
	prevcen = -1
	edrval = 0d0
	
	if (prop==3) then
		edrdmaxval = 0d0	
	end if
	
	dens=rho(orbst,orbend,x,y,z)

	if (dens.GT.1D-10) then
		psi = 0d0
		Bint=0d0	
	    do j = 1, nprims
			tx = type2ix(b(j)%type)
			ty = type2iy(b(j)%type)
			tz = type2iz(b(j)%type)
			epp = b(j)%expo
			
			if (b(j)%center/=prevcen) then 
				shftx = x-a(b(j)%center)%xcoo
				shfty = y-a(b(j)%center)%ycoo
				shftz = z-a(b(j)%center)%zcoo
				shftx2 = shftx*shftx
				shfty2 = shfty*shfty
				shftz2 = shftz*shftz
				sqr = shftx2 +  shfty2 + shftz2 
			end if
			
			exterm = 0.0
			amu0 = 0d0	
			if (expcutoff>0 .OR. -epp*sqr>expcutoff) then
				exterm = exp(-epp*sqr)
				do iedr=1,nedr
					amu0(iedr)=(2d0*ed(iedr)/pi)**(3d0/4d0) *(pi/(epp+ed(iedr)))**(3d0/2d0)&
					* exp(-epp*ed(iedr)/(epp+ed(iedr))*sqr)
				end do
			else 
				exterm = 0d0
			end if	
			
			prevcen = b(j)%center 
			
			if (exterm==0d0) cycle

			gtval = shftx**tx * shfty**ty * shftz**tz * exterm

			do iedr=1,nedr
				do ixyz=1,3
					ival=tx
					sftval=shftx
					if(ixyz.eq.2) then
						ival=ty
						sftval=shfty
					else if(ixyz.eq.3) then
						ival=tz
						sftval=shftz
					end if
					If(ival.eq.0) then
						xamu(ixyz,iedr)=1d0 
					else if(ival.eq.1) then 
						xamu(ixyz,iedr)=sftval*ed(iedr)/(ed(iedr)+epp)
					else If(ival.eq.2) then
						xamu(ixyz,iedr)=(sftval*ed(iedr)/(ed(iedr)+epp))**2d0 + 1d0/(2d0*(ed(iedr)+epp))
					else If(ival.eq.3) then
						xamu(ixyz,iedr)=(sftval*ed(iedr)/(ed(iedr)+epp))**3d0 + sftval*3d0*ed(iedr)&
						/(2d0*(ed(iedr)+epp)**2d0)
					else If(ival.eq.4) then
						xamu(ixyz,iedr)=sftval**4d0* (ed(iedr)/(ed(iedr)+epp))**4d0 + sftval**2d0 &
						*3d0*ed(iedr)**2d0/(ed(iedr)+epp)**3d0 + 3d0/(4d0*(ed(iedr)+epp)**2d0)
					else If(ival.eq.5) then
						xamu(ixyz,iedr)=sftval**5d0*  ed(iedr)**5d0/(ed(iedr)+epp)**5d0 + sftval**3d0* &
						5d0*ed(iedr)**3d0/(ed(iedr)+epp)**4d0 + sftval *15d0*ed(iedr)/(4d0*(ed(iedr)+epp)**3d0)
					else 
					write(*,*) "Angular momentum out of range"
					Call EXIT()
					end if
				end do
			end do
			do iedr=1,nedr
				AMUVal(iedr)=amu0(iedr)*xamu(1,iedr)*xamu(2,iedr)*xamu(3,iedr)
			end do
			
			do i = orbst,orbend
				if (nint(MOocc(i)).GE.1D0) then	
					psi(i) = psi(i) + MOco(i,j)*gtval
				end if
			end do
			
			do iedr = 1,nedr
				do i = orbst,orbend
					if (nint(MOocc(i)).GE.1D0) then	
						Bint(i,iedr) = Bint(i,iedr) + MOco(i,j)*AMUVal(iedr)
					end if
				end do
			end do	
		end do
		edrval = 0d0
		do i = orbst,orbend
			do iedr = 1,nedr
				edrval(iedr)=edrval(iedr)+psi(i)*Bint(i,iedr)
			end do
		end do
		
		if (prop==3) then
			call three_point_interpolation(nedr,ed,edrval,edmax,dmaxdummy)
			edrdmaxval=edmax
		else if (prop==2) then
			do iedr=1,nedr
				edrval(iedr)=edrval(iedr)*dens**(-0.5D0)
			end do
		else	
			write(*,*) " The selected option is not available."
			write(*,*) " Select option 2 or 3."
			call EXIT()
		end if
	end if
	
end subroutine EDRcal

! ======================================
! = Subroutine: 3 Points Interpolation =
! ======================================
subroutine three_point_interpolation(n,x,y,xmax,ymax)
	integer, intent(in) :: n
	real*8, intent(in) :: x(max_edr_exponents),y(max_edr_exponents)
	real*8, intent(out) :: xmax,ymax
	integer i , imax
	real*8  x1,x2,x3, y1,y2,y3, a,b
100 format ('XXX ',3F9.5)
	ymax = -1.0d0
	imax = -1 
	do i=1,n
		if(y(i) .gt. ymax) then 
			ymax = y(i)
			imax = i
		endif 
	end do 
	if(imax<1 .or. imax>n) then
		write(*,*) "Error: Bad imax"
		call EXIT()
	end if
	if(imax .eq. 1 .or. imax.eq.n) then
		xmax = x(imax)**(-0.5d0)
		return 
	endif
	x1 = x(imax-1)**(-0.5d0)
	x2 = x(imax  )**(-0.5d0)
	x3 = x(imax+1)**(-0.5d0)
	y1 = y(imax-1)
	y2 = y(imax  )
	y3 = y(imax+1)
	a = ( (y3-y2)/(x3-x2) -(y2-y1)/(x2-x1) )/(x3-x1)
	b = ( (y3-y2)/(x3-x2)*(x2-x1) + (y2-y1)/(x2-x1)*(x3-x2) )&
      /(x3-x1)
	xmax = x2 - b/(2d0*a)
	ymax = y2 - b**(2d0)/(4d0*a)
end subroutine 

end module functions
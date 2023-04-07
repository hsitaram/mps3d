module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=10
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=3	  !electron species
      integer, parameter :: nreac=15
      real*8, parameter  :: cm_to_m = 0.01
      character (LEN=10), dimension(nspecies) :: specnames
      integer :: reactants(nreac,nspecies)
      integer :: products(nreac,nspecies)
      real*8  :: k_arrh(3,nreac)
      real*8  :: elecenergy(nreac)
      real*8  :: gasenergy(nreac)
      logical :: isratearrh(nreac)
      real*8  :: molmass(nspecies)
      real*8  :: spec_charge(nspecies)
      
      integer :: ionspecmin
      integer :: ionspecmax
      integer :: neutralspecmin
      integer :: neutralspecmax
      integer :: no_of_ions
      integer :: no_of_neutrals

      contains
!====================================================================
subroutine assignreactions()

	integer :: H2,E,Hp,H2p,H
	integer :: O2,Op,O2p,Om,O
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

	H2  = 1
	O2  = 2
	E   = 3
	Hp  = 4
	H2p = 5
	Op  = 6
	O2p = 7
	Om  = 8
	H   = 9
	O   = 10

	!convention for energy is added is positive and
	!removed is negative.
	!eg: electron loses 19.8 eV when reaction 
	!1 happens

	rnum=1
	!***********************************
	!H + E --> H^+ + 2E
	!***********************************
	reactants(rnum,H) = 1
	reactants(rnum,E) = 1
	
	products(rnum,Hp)  = 1	
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 6.5023d-9*(cm_to_m**3)/(eVtoK**0.48931)
	k_arrh(2,rnum)   = 0.48931
	k_arrh(3,rnum)   = 149624.36
	!-----------------------------------
	elecenergy(rnum) = -12.89365
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> H^+ + H + 2E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E ) = 1

	products(rnum,Hp)  = 1
	products(rnum,H)   = 1
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.9962d-8*(cm_to_m**3)/(eVtoK**0.44456)
	k_arrh(2,rnum)   = 0.44456
	k_arrh(3,rnum)   = 437818.75
	!-----------------------------------
	elecenergy(rnum) = -37.72836
        gasenergy(rnum)  = 0.0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + E --> H^+ + H + E
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Hp)   = 1
	products(rnum,H)    = 1
	products(rnum,E)    = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.0702d-7*(cm_to_m**3)/(eVtoK**0.04876)
	k_arrh(2,rnum)   = 0.04876
	k_arrh(3,rnum)   = 112450.85
	!-----------------------------------
	elecenergy(rnum) = -9.69028
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + E --> H+ + H+ + 2E
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Hp)  = 2
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.1202d-9*(cm_to_m**3)/(eVtoK**0.31394)
	k_arrh(2,rnum)   = 0.31394
	k_arrh(3,rnum)   = 270371.5
	!-----------------------------------
	elecenergy(rnum) = -23.29885
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + H --> H_2 + H^+
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,H)   = 1
	
	products(rnum,H2)   = 1
	products(rnum,Hp)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 9.0d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 0.0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + H+ --> H_2+ + H
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,Hp) = 1

	products(rnum,H2p)  = 1
	products(rnum,H)    = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.19d-22*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> H_2^+ + 2E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E)  = 1

	products(rnum,H2p)  = 1
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 3.1228d-8*(cm_to_m**3)/(eVtoK**0.17156)
	k_arrh(2,rnum)   = 0.17156
	k_arrh(3,rnum)   = 232987.49
	!-----------------------------------
	elecenergy(rnum) = -20.07734
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> 2H + E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E)  = 1

	products(rnum,H)   = 2
	products(rnum,E)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.7527d-7*(cm_to_m**3)/(eVtoK**(-1.23668))
	k_arrh(2,rnum)   = -1.23668
	k_arrh(3,rnum)   = 146128.85
	!-----------------------------------
	elecenergy(rnum) = -12.59243
	gasenergy(rnum)  =   0.0
	!+++++++++++++++++++++++++++++++++++
	
	
	rnum=rnum+1
	!***********************************
	!O_2 + E --> O_2^+ + 2E
	!***********************************
	reactants(rnum,O2)  = 1
	reactants(rnum,E)   = 1
	
	products(rnum,O2p)  = 1	
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	!Arrhenius rates from Lee et al., Phys. Plasmas.,13, (057102) 2006
	k_arrh(1,rnum)   = 9.0d-10*(cm_to_m**3)/(eVtoK**0.5)
	k_arrh(2,rnum)   = 0.5
	k_arrh(3,rnum)   = 146216.7
	!-----------------------------------
	elecenergy(rnum) = -12.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O_2 + E --> 2 O + E
	!***********************************
	reactants(rnum,O2) = 1
	reactants(rnum,E ) = 1

	products(rnum,O)   = 2
	products(rnum,E)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 4.23d-9*(cm_to_m**3)/(eVtoK**0.0)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 64521.02
	!-----------------------------------
	!Energy from Lee et al., Phys. Plasmas.,13, (057102) 2006
	elecenergy(rnum) = -6.4
        gasenergy(rnum)  = 0.0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O_2 + E --> O + O-
	!***********************************
	reactants(rnum,O2)  = 1
	reactants(rnum,E)   = 1
	
	products(rnum,O)    = 1
	products(rnum,Om)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	!Arrhenius rates from Lee et al., Phys. Plasmas.,13, (057102) 2006
	k_arrh(1,rnum)   = 8.8d-11*(cm_to_m**3)/(eVtoK**0.0)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 51059.8
	!-----------------------------------
	elecenergy(rnum) = -3.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O + E --> O+ + 2E
	!***********************************
	reactants(rnum,O)   = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Op)   = 1
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 9.d-9*(cm_to_m**3)/(eVtoK**0.7)
	k_arrh(2,rnum)   = 0.7
	k_arrh(3,rnum)   = 157821.2
	!-----------------------------------
	elecenergy(rnum) = -13.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + O_2^+ --> O + O_2
	!***********************************
	reactants(rnum,Om)    = 1
	reactants(rnum,O2p)   = 1
	
	products(rnum,O)      = 1
	products(rnum,O2)     = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	!this reaction rate is calculated at 300 K gas 
	!temperature from Doug's thesis
	k_arrh(1,rnum)   =  1.99d-7*(cm_to_m**3)
	k_arrh(2,rnum)   =  0.0
	k_arrh(3,rnum)   =  0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 10.69
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + O+ --> O + O
	!***********************************
	reactants(rnum,Om) = 1
	reactants(rnum,Op) = 1

	products(rnum,O)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	!this reaction rate is calculated at 300 K gas 
	!temperature from Doug's thesis
	k_arrh(1,rnum)   = 2.66d-7*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 12.69
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + E --> O + 2E
	!***********************************
	reactants(rnum,Om) = 1
	reactants(rnum,E)  = 1

	products(rnum,O)   = 1
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	
	!Arrhenius rates from Lee et al., Phys. Plasmas.,13, (057102) 2006
	k_arrh(1,rnum)   = 2.0d-7*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 63824.75

	!Arrhenius rates from Doug's thesis
	!k_arrh(1,rnum)   = 2.10d-10*(cm_to_m**3)
	!k_arrh(2,rnum)   = 0.5
	!k_arrh(3,rnum)   = 39400.0
	!-----------------------------------
	
	!Energy from Doug's thesis
	!elecenergy(rnum) = -0.92

	!Energy from Lee et al., Phys. Plasmas.,13, (057102) 2006
	elecenergy(rnum) = -5.5
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1)  = 'H2'
	      specnames(2)  = 'O2'
	      specnames(3)  = 'E'
	      specnames(4)  = 'H+'
	      specnames(5)  = 'H2+'
	      specnames(6)  = 'O+'
	      specnames(7)  = 'O2+'
	      specnames(8)  = 'O-'
	      specnames(9)  = 'H'
	      specnames(10) = 'O'

	      molmass(1)  =  2.d0*mass_prot
	      molmass(2)  = 32.d0*mass_prot
	      molmass(3)  =       mass_elec
	      molmass(4)  =  1.d0*mass_prot
	      molmass(5)  =  2.d0*mass_prot
	      molmass(6)  = 16.d0*mass_prot
	      molmass(7)  = 32.d0*mass_prot
	      molmass(8)  = 16.d0*mass_prot
	      molmass(9)  =  1.d0*mass_prot
	      molmass(10) = 16.d0*mass_prot

	      spec_charge(1)  =  0.d0
	      spec_charge(2)  =  0.d0
	      spec_charge(3)  = -1.d0
	      spec_charge(4)  =  1.d0
	      spec_charge(5)  =  1.d0
	      spec_charge(6)  =  1.d0
	      spec_charge(7)  =  1.d0
	      spec_charge(8)  = -1.d0
	      spec_charge(9)  =  0.d0
	      spec_charge(10) =  0.d0

	      ionspecmin = 4
	      ionspecmax = 8

	      no_of_ions = 5

	      neutralspecmin = 9
	      neutralspecmax = 10

	      no_of_neutrals = 2

end subroutine setspecparams
!====================================================================
subroutine initializechemistry()

	call setspecparams()
	call assignreactions()

end subroutine initializechemistry
!====================================================================
subroutine getspecproduction(specnum,Te,Tg,specden,specprod,efld)

	integer, intent(in)  :: specnum
	real*8, intent(in)   :: Te,Tg
	real*8, intent(in)   :: specden(nspecies)
	real*8,intent(inout) :: specprod
	real*8,optional :: efld
	real*8 :: efield
	real*8 :: k
	integer :: i,j
	integer :: nmoles
	integer :: nreacmoles
	real*8 :: specmult

	if(present(efld)) then
		efield=efld
	else
		efield=1.0
	endif
	
	specprod = 0.d0
	do i=1,nreac

		nmoles = products(i,specnum)-reactants(i,specnum)

		if(nmoles .ne. 0) then

			if(isratearrh(i) .eqv. .true.) then
				k=getarrhrate(k_arrh(:,i),Te)
			else
				k=getcustomrate(i,k_arrh(:,i),Te,efield)
			endif

			specmult = 1.d0

			do j=1,nspecies
				nreacmoles = reactants(i,j)
				specmult = specmult*(specden(j)**nreacmoles)
			enddo

			specprod = specprod + nmoles*k*specmult
			!write(*,'(A,I5,A,I5,E20.10,E20.10,E20.10)')"reac:",i,"specprod:", &
			!		nmoles,k,specmult,nmoles*k*specmult
		endif
	enddo

end subroutine getspecproduction
!====================================================================
subroutine getelectroninelasticterm(Te,Tg,specden,inelterm,efld)
	
	real*8, intent(in)   :: Te,Tg
	real*8, intent(in)   :: specden(nspecies)
	real*8, intent(inout) :: inelterm
	real*8,optional :: efld
	integer :: i,j,nreacmoles
	real*8 :: k,specmult
	real*8 :: efield

	inelterm = 0.d0
	
	if(present(efld)) then
		efield=efld
	else
		efield=1.0
	endif
	
	do i=1,nreac

		if(elecenergy(i) .ne. 0) then
			
				if(isratearrh(i) .eqv. .true.) then
					k=getarrhrate(k_arrh(:,i),Te)
				else
					k=getcustomrate(i,k_arrh(:,i),Te,efield)
				endif

				specmult = 1.d0
				do j=1,nspecies
					nreacmoles = reactants(i,j)
					!print *,"nreacmoles at ",j,nreacmoles
					specmult = specmult*(specden(j)**nreacmoles)
				enddo

				inelterm = inelterm + k*specmult*elecenergy(i)
				!write(*,'(A,I5,A,E20.10,E20.10,E20.10,E20.10)')"reac:",i,"inelterm:", &
				!		k,specmult,elecenergy(i),k*specmult*elecenergy(i)

		endif
	enddo

end subroutine getelectroninelasticterm
!====================================================================
function getspecdcoeff(specnum,specarray,elecfield,Te,Tg,Pg)  result(dcoeff)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in) :: specarray(nspecies)
	real*8, intent(in) :: elecfield

	real*8 :: dcoeff
	real*8 :: Patm,mob
	
	integer :: E,Hp,H2p,Op,O2p,Om
	integer :: O,H
	
	E   = 3
	Hp  = 4
	H2p = 5
	Op  = 6
	O2p = 7
	Om  = 8
	H   = 9
	O   = 10

	Patm = 101325.d0

	if(specnum .eq. E) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff = k_B*Te*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. Hp) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. H2p) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. Op) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. O2p) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. Om) then
		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. H) then
		
		dcoeff=0.03
	
	else if(specnum .eq. O) then
		
		dcoeff=0.0075
	else
		write(*,*)"species does not exist"
		stop
	endif

end function getspecdcoeff
!====================================================================
function getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg)  result(mobility)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in)   :: specarray(nspecies)
	real*8, intent(in)   :: elecfield

	real*8 :: mobility,mob_inv_H2,mob_inv_O2
	real*8 :: Patm,Natm,Troom
	real*8 :: neutral_den

	integer :: E,Hp,H2p,Op,O2p,Om,H2,O2
	real*8 :: Te_in_eV
	real*8 :: log_Te
	
	H2  = 1
	O2  = 2
	E   = 3
	Hp  = 4
	H2p = 5
	Op  = 6
	O2p = 7
	Om  = 8

	neutral_den = Pg/k_B/Tg
	Te_in_eV = Te/eVtoK
	log_Te   = log(Te_in_eV)

	Patm = 101325.d0
	Troom = 300.d0
	Natm = Patm/k_B/Troom

	if(specnum .eq. E) then

		mob_inv_H2 = specarray(H2)/exp(0.061*log_Te**2 - 0.342*log_Te + 55.89)
		mob_inv_O2 = specarray(O2)/exp(0.031*log_Te**2 - 0.658*log_Te + 56.71)  
		mobility = -1.d0/(mob_inv_H2 + mob_inv_O2)

	else if(specnum .eq. Hp) then
		
		mob_inv_H2 = specarray(H2)/(15.d0*(cm_to_m**2) )
		mob_inv_O2 = specarray(O2)/(15.d0*(cm_to_m**2) )
		mobility   = Natm/(mob_inv_H2+mob_inv_O2)

	else if(specnum .eq. H2p) then
		
		mob_inv_H2 = specarray(H2)/(13.d0*(cm_to_m**2) )
		mob_inv_O2 = specarray(O2)/(13.d0*(cm_to_m**2) )
		mobility   = Natm/(mob_inv_H2+mob_inv_O2)

	else if(specnum .eq. Op) then

		mob_inv_H2 = specarray(H2)/(11.d0*(cm_to_m**2) )
		mob_inv_O2 = specarray(O2)/(3.5d0*(cm_to_m**2) )
		mobility   = Natm/(mob_inv_H2+mob_inv_O2)

	else if(specnum .eq. O2p) then
	
		mob_inv_H2 = specarray(H2)/(8.d0 *(cm_to_m**2) )
		mob_inv_O2 = specarray(O2)/(2.1d0*(cm_to_m**2) )
		mobility   = Natm/(mob_inv_H2+mob_inv_O2)
	
	else if(specnum .eq. Om) then
	
		mob_inv_H2 = specarray(H2)/(24.d0*(cm_to_m**2) )
		mob_inv_O2 = specarray(O2)/(3.8d0*(cm_to_m**2) )
		mobility   = -Natm/(mob_inv_H2+mob_inv_O2)
	else
		write(*,*)"species does not exist"
		stop
	endif

end function getspecmobility
!====================================================================
function getelectroncollisionfrequency(specarray,elecfield,Te,Tg,Pg) result(collfreq)
	
	real*8, intent(in) :: Te,Tg,Pg
	real*8, intent(in) :: specarray(nspecies)
	real*8, intent(in) :: elecfield
	real*8 :: collfreq
	real*8 :: mu
	real*8 :: m_e


	mu = getspecmobility(especnum,specarray,elecfield,Te,Tg,Pg)
	m_e = molmass(especnum)

	collfreq = -echarge/m_e/mu

end function getelectroncollisionfrequency
!====================================================================
function getarrhrate(k,T) result(rateconst)

	real*8, intent(in) :: k(3)
	real*8, intent(in) :: T
	real*8 :: rateconst

	real*8 :: A,alpha,Ea

	A     = k(1)
	alpha = k(2)
	Ea    = k(3)

	rateconst = A*(T**alpha)*exp(-Ea/T)

end function getarrhrate
!====================================================================
function getcustomrate(reacnum,kparams,T,efield) result(rateconst)
	
	integer,intent(in) :: reacnum
	real*8, intent(in) :: kparams(3)
	real*8, intent(in) :: T,efield
	real*8 :: rateconst
	real*8 :: expterm

	real*8 :: T_in_eV

	T_in_eV = T/(eVtoK)

	if(reacnum .eq. 9) then
		!Fit from shankar's thesis	
		rateconst = 1.043d-13
		expterm = -2.74d5/T - 2.0015d8/(T**2) + 5.33d13/(T**3) - 4.001d17/(T**4)
		rateconst=rateconst*exp(expterm)

	else if(reacnum .eq. 10) then
		
		!Fit from shankar's thesis	
		rateconst = 9.577d-16
		expterm = 1.246d5/T - 8.647d9/(T**2) + 1.381d14/(T**3) - 6.953d17/(T**4)
		rateconst=rateconst*exp(expterm)

	else if(reacnum .eq. 11) then
		!Fit from shankar's thesis	
		rateconst = 1.46d-14
		expterm = 6.21d4/T - 7.27d9/(T**2) + 1.25d14/(T**3) - 6.57d17/(T**4)

		!Fit from Lee et al.,J. Electrochem. Soc.,141,6,1994
		!rateconst = 4.6d-11*(cm_to_m**3)
		!expterm   = (-2.91/T_in_eV + 12.6/(T_in_eV**2) - 6.92/(T_in_eV**3))
		
		rateconst=rateconst*exp(expterm)
	
	else if(reacnum .eq. 15) then
		!Fit from Lee et al.,J. Electrochem. Soc.,141,6,1994
		rateconst = 1.73d-7*(cm_to_m**3)
		expterm = -5.56/T_in_eV + 7.3/(T_in_eV**2) - 3.48/(T_in_eV**3)
		rateconst=rateconst*exp(expterm)
	else
		print *,"No custom rate for reaction",reacnum
		stop
	endif

end function getcustomrate
!====================================================================

end module chem_module

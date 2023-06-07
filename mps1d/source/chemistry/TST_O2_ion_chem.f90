module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=6
      integer, parameter :: bgspecnum=6   !background gas
      integer, parameter :: especnum=1	  !electron species
      integer, parameter :: nreac=7
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
	  real*8 :: ydot(nspecies)
      
      integer :: ionspecmin
      integer :: ionspecmax
      integer :: neutralspecmin
      integer :: neutralspecmax
      integer :: no_of_ions
      integer :: no_of_neutrals

      contains
!====================================================================
subroutine assignreactions()

	integer :: E,Op,O2p,Om,O,O2
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

    E = 1
    Op = 2
    O2p = 3
    Om = 4
    O = 5
    O2 = 6

end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1) = 'E'
	      specnames(2) = 'O^+'
	      specnames(3) = 'O2^+'
	      specnames(4) = 'O^-'
	      specnames(5) = 'O'
	      specnames(6) = 'O2'

	      molmass(1) = mass_elec
	      molmass(2) = 16.d0*mass_prot
	      molmass(3) = 32.d0*mass_prot
	      molmass(4) = 16.d0*mass_prot
	      molmass(5) = 16.d0*mass_prot
	      molmass(6) = 32.d0*mass_prot

	      spec_charge(1) = -1.d0
	      spec_charge(2) = 1.d0
	      spec_charge(3) = 1.d0
	      spec_charge(4) = -1.d0
	      spec_charge(5) = 0.d0
	      spec_charge(6) = 0.d0          

	      ionspecmin = 2
	      ionspecmax = 4

	      no_of_ions = 3

	      neutralspecmin = 5
	      neutralspecmax = 5

	      no_of_neutrals = 1

end subroutine setspecparams
!====================================================================

subroutine read_elec_and_gas_energies()

	integer :: i, energyinpfptr

	energyinpfptr = 230
	open(unit=energyinpfptr,file="elec_gas_energy_O2_ion.inp")
	do i = 1,nreac
		read(energyinpfptr,*) elecenergy(i),gasenergy(i)
	enddo

end subroutine read_elec_and_gas_energies

subroutine initializechemistry()

	call setspecparams()
	call assignreactions()
	call read_elec_and_gas_energies()

end subroutine initializechemistry
!====================================================================

function zdp_getreactionrates(Te, Tgas) result(rrt)

        real*8 :: rrt(nreac)
        real*8, intent(in) :: Te, Tgas
 
        rrt(1) = (9.0D-10)*(TE**0.5)*EXP(-146216.7/TE) / (eVtoK**0.5)
        rrt(2) = (4.23D-9)*EXP(-64521.02/TE)
        rrt(3) = (8.8D-11)*EXP(-51059.8/TE)
        rrt(4) = (9.0D-09)*(TE**0.7)*EXP(-157821.2/TE) / (eVtoK**0.7)
        rrt(5) = 1.99D-07
        rrt(6) = 2.66D-07
        rrt(7) = 2.0D-07*EXP(-63824.75/TE)
  
        rrt = rrt*cm_to_m**3

end function zdp_getreactionrates

function zdp_getratesofprogress(Te,Tg,specden) result(rrt)

    	!integer, intent(inout) :: specnum
    	real*8, intent(in) :: Te, Tg
 		real*8, intent(in) :: specden(nspecies)
    	real*8 :: density(nspecies)
    	real*8 :: rrt(nreac)
		integer :: j
    	!real*8 :: ydot(nspecies)        

    	density = specden

		rrt = zdp_getreactionrates(Te,Tg)
        
        rrt(1) = rrt(1) * density(1) * density(6) 
        rrt(2) = rrt(2) * density(1) * density(6) 
        rrt(3) = rrt(3) * density(1) * density(6) 
        rrt(4) = rrt(4) * density(1) * density(5) 
        rrt(5) = rrt(5) * density(3) * density(4) 
        rrt(6) = rrt(6) * density(2) * density(4) 
        rrt(7) = rrt(7) * density(1) * density(4) 

end function zdp_getratesofprogress

subroutine zdp_getspecproduction(Te,Tg,specden) 

    	real*8, intent(in) :: Te, Tg
 		real*8, intent(in) :: specden(nspecies)
		real*8 :: rrt(nreac)

		rrt = zdp_getratesofprogress(Te,Tg,specden)

        ydot(1) = +rrt(1)-rrt(3)+rrt(4)+rrt(7) 
        ydot(2) = +rrt(4)-rrt(6) 
        ydot(3) = +rrt(1)-rrt(5) 
        ydot(4) = +rrt(3)-rrt(5)-rrt(6)-rrt(7) 
        ydot(5) = +2.d0 * rrt(2)+rrt(3)-rrt(4)+rrt(5)+2.d0 * rrt(6)+rrt(7) 
        ydot(6) = -rrt(1)-rrt(2)-rrt(3)+rrt(5) 

end subroutine zdp_getspecproduction

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

subroutine zdp_getelectroninelasticterm(Te,Tg,specden,inelterm,rrt)

	real*8 :: rrt(nreac)
	integer :: i
	real*8 :: Te, Tg, specden(nspecies)
	real*8 :: inelterm

	inelterm = 0.0

	do i = 1,nreac
		inelterm = inelterm + rrt(i)*elecenergy(i)
		!write(*,*) "elecenergy(",i,") = ", elecenergy(i)
	enddo

end subroutine zdp_getelectroninelasticterm


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

!	call zdp_getelectroninelasticterm()

end subroutine getelectroninelasticterm
!====================================================================
function getspecdcoeff(specnum,specarray,elecfield,Te,Tg,Pg)  result(dcoeff)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in) :: specarray(nspecies)
	real*8, intent(in) :: elecfield

	real*8 :: dcoeff
	real*8 :: Patm,mob
	
	integer :: E,Op,O2p,Om,O,O2
	
    E = 1
    Op = 2
    O2p = 3
    Om = 4
    O = 5
    O2 = 6

	Patm = 101325.d0

	if(specnum .eq. E) then
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff = k_B*Te*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. Op) then		
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. O2p) then
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. Om) then
		mob = getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg) 
		dcoeff= k_B*Tg*mob/(spec_charge(specnum)*echarge)

	else if(specnum .eq. O) then
		dcoeff=0.0075

	else
		write(*,*)"species does not exist"
	endif

end function getspecdcoeff
!====================================================================
function getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg)  result(mobility)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in)   :: specarray(nspecies)
	real*8, intent(in)   :: elecfield

	real*8 :: mobility,mob_inv_O2
	real*8 :: Patm,Natm,Troom
	real*8 :: neutral_den

	integer :: E,Op,O2p,Om,O,O2
	real*8 :: Te_in_eV
	real*8 :: log_Te
	
    E = 1
    Op = 2
    O2p = 3
    Om = 4
    O = 5
    O2 = 6

	neutral_den = Pg/k_B/Tg
	Te_in_eV = Te/eVtoK
	log_Te   = log(Te_in_eV)
	
	Patm = 101325.d0
	Troom = 300.d0
	Natm = Patm/k_B/Troom

	if(specnum .eq. E) then
		mob_inv_O2 = specarray(O2)/exp(0.031*log_Te**2 - 0.658*log_Te + 56.71)
		mobility = -1.d0/(mob_inv_O2)

	else if(specnum .eq. Op) then
		mob_inv_O2 = specarray(O2)/(3.5d0*(cm_to_m**2))
		mobility   = Natm/(mob_inv_O2)        

	else if(specnum .eq. O2p) then
		mob_inv_O2 = specarray(O2)/(2.1d0*(cm_to_m**2) )
		mobility   = Natm/(mob_inv_O2)    

	else if(specnum .eq. Om) then
		mob_inv_O2 = specarray(O2)/(3.8d0*(cm_to_m**2) )
		mobility   = -Natm/(mob_inv_O2)

    else
		write(*,*)"species does not exist"
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
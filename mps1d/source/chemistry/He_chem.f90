module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=6
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=2	  !electron species
      integer, parameter :: nreac=9
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

	integer :: HE,E,HEm,HEp,HE2m,HE2p
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = ZERO
	gasenergy  = ZERO
	k_arrh     = ZERO
	isratearrh = .true.

	HE   = 1
	E    = 2
	HEp  = 3
	HE2p = 4
        HEm  = 5
	HE2m = 6

	!convention for energy is added is positive and
	!removed is negative.
        !reaction rates obtained from Yuan and Raja,
        !IEEE. Trans. Plasma Sci.,31,4,2003

	rnum=1
	!***********************************
	!E + HE --> HEm + E
	!***********************************
	reactants(rnum,E)  = 1
	reactants(rnum,HE) = 1
	
	products(rnum,HEm)  = 1	
	products(rnum,E)    = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.308d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.31
	k_arrh(3,rnum)   = 2.297d5
	!-----------------------------------
	elecenergy(rnum) =  -19.8
	gasenergy(rnum)  =  ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!E + HEm --> HE + E
	!***********************************
	reactants(rnum,E )  = 1
	reactants(rnum,HEm) = 1

	products(rnum,HE)  = 1
	products(rnum,E)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.099d-11*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.31
	k_arrh(3,rnum)   = ZERO
	!-----------------------------------
	elecenergy(rnum) = 19.8
        gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!E + HE --> HE+ + 2E
	!***********************************
	reactants(rnum,E)    = 1
	reactants(rnum,HE)   = 1
	
	products(rnum,HEp)   = 1
	products(rnum,E)     = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.584d-12*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.68
	k_arrh(3,rnum)   = 2.854092d5
	!-----------------------------------
	elecenergy(rnum) = -24.6
	gasenergy(rnum)  =  ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!E + HEm --> HE+ + 2E
	!***********************************
	reactants(rnum,E)   = 1
	reactants(rnum,HEm) = 1
	
	products(rnum,HEp)  = 1
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 4.661d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.6
	k_arrh(3,rnum)   = 5.546d4
	!-----------------------------------
	elecenergy(rnum) = -4.78
	gasenergy(rnum)  =  ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!E + HE2m --> HE2p + 2E
	!***********************************
	reactants(rnum,E)    = 1
	reactants(rnum,HE2m) = 1
	
	products(rnum,HE2p)   = 1
	products(rnum,E)      = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.268d-12*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.71
	k_arrh(3,rnum)   = 3.945d4
	!-----------------------------------
	elecenergy(rnum) = -3.4
	gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!E + HE2p --> HEm + HE
	!***********************************
	reactants(rnum,E)    = 1
	reactants(rnum,HE2p) = 1

	products(rnum,HEm)  = 1
	products(rnum,HE)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 5.386d-7*(cm_to_m**3)
	k_arrh(2,rnum)   = -0.5
	k_arrh(3,rnum)   = ZERO
	!-----------------------------------
	elecenergy(rnum) = ZERO
	gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!2 HEm --> HE+ + HE + E
	!***********************************
	reactants(rnum,HEm) = 2

	products(rnum,HEp)  = 1
	products(rnum,HE)   = 1 
	products(rnum,E)    = 1 
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.7d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = ZERO
	k_arrh(3,rnum)   = ZERO
	!-----------------------------------
	elecenergy(rnum) = 15.d0
	gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!HEm + 2HE --> HE2m + HE
	!***********************************
	reactants(rnum,HEm) = 1
	reactants(rnum,HE)  = 2

	products(rnum,HE2m) = 1
	products(rnum,HE)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.3d-33*(cm_to_m**3)
	k_arrh(2,rnum)   = ZERO
	k_arrh(3,rnum)   = ZERO
	!-----------------------------------
	elecenergy(rnum) = ZERO
	gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++
	
        
        rnum=rnum+1
	!***********************************
	!HE+ + 2HE --> HE2+ + HE
	!***********************************
	reactants(rnum,HEp) = 1
	reactants(rnum,HE)  = 2

	products(rnum,HE2p) = 1
	products(rnum,HE)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.d-31*(cm_to_m**3)
	k_arrh(2,rnum)   = ZERO
	k_arrh(3,rnum)   = ZERO
	!-----------------------------------
	elecenergy(rnum) = ZERO
	gasenergy(rnum)  = ZERO
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1) = 'HE'
	      specnames(2) = 'E'
	      specnames(3) = 'HE+'
	      specnames(4) = 'HE2+'
	      specnames(5) = 'HEm'
	      specnames(6) = 'HE2m'

	      molmass(1) = 4.d0*mass_prot
	      molmass(2) = mass_elec
	      molmass(3) = 4.d0*mass_prot
	      molmass(4) = 8.d0*mass_prot
	      molmass(5) = 4.d0*mass_prot
	      molmass(6) = 8.d0*mass_prot

	      spec_charge(1) =  ZERO
	      spec_charge(2) = -ONE
	      spec_charge(3) =  ONE
	      spec_charge(4) =  ONE
	      spec_charge(5) =  ZERO
	      spec_charge(6) =  ZERO

	      ionspecmin = 3
	      ionspecmax = 4

	      no_of_ions = 2

	      neutralspecmin = 5
	      neutralspecmax = 6

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
				!print *,"reac:",i,"inelterm:", k,specmult,elecenergy(i),inelterm

		endif
	enddo

end subroutine getelectroninelasticterm
!====================================================================
function getspecdcoeff(specnum,specarray,elecfield,Te,Tg,Pg)  result(dcoeff)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in) :: specarray(nspecies)
	real*8, intent(in) :: elecfield

	integer :: HE,E,HEm,HEp,HE2m,HE2p
	real*8 :: dcoeff
	real*8 :: Patm,Pres_ratio

	Patm = ONE_ATM_IN_PA
        Pres_ratio = Patm/Pg

	HE   = 1
	E    = 2
	HEp  = 3
	HE2p = 4
        HEm  = 5
	HE2m = 6

	if(specnum .eq. E) then
		dcoeff=0.1737*(Te/17406.d0)*Pres_ratio
	else if(specnum .eq. HEp) then
		dcoeff=5.026d-5*Pres_ratio
	else if(specnum .eq. HE2p) then
		dcoeff=8.148d-5*Pres_ratio
	else if(specnum .eq. HEm) then
		dcoeff=4.116d-4*Pres_ratio
	else if(specnum .eq. HE2m) then
		dcoeff=2.029d-4*Pres_ratio
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

	integer :: HE,E,HEm,HEp,HE2m,HE2p
	real*8 :: mobility
	real*8 :: Patm,Pres_ratio
        real*8 :: EbyN,N, Te_in_eV, meanE
        real*8 :: one_td,Hepmob,He2pmob,elecmob

        one_td=1e-21 !V m2

	Patm = ONE_ATM_IN_PA

        N=Pg/k_B/Tg
        EbyN=abs(elecfield)/N/one_td

	HE   = 1
	E    = 2
	HEp  = 3
	HE2p = 4
        HEm  = 5
	HE2m = 6
        
        Pres_ratio = (Patm/Pg)
        Te_in_eV = Te/eVtoK
        meanE = Te_in_eV*(1.5d0)

        
        !mobility from 
        !Turner, Miles M., et al. "Simulation benchmarks for low-pressure
        !plasmas: Capacitive discharges." 
        !Physics of Plasmas 20.1 (2013): 013507.
        Hepmob=2.69*(1+1.2d-3*(EbyN**2)+4.2d-8*(EbyN)**4)**(-0.125)
        He2pmob=Hepmob

        !computed by Taaresh (U Minnesota) using Turner's cross sections
        elecmob = (-1.0)*exp(55.0 + 0.3942*log(meanE) + 2.134/meanE &
                   -0.6433/meanE**2 + (0.7112d-01)/meanE**3) / N

        !reduced mobility obtained from Yuan and Raja,
        !IEEE. Trans. Plasma Sci.,31,4,2003
        !they say these are 393K, may be temperature
        !scaling is necessary
        !elecmob = -0.1132*Pres_ratio
	!Hepmob = 0.001482*Pres_ratio
	!He2pmob = 0.002403*Pres_ratio

	if(specnum .eq. E) then
		mobility = elecmob
	else if(specnum .eq. HEp) then
		mobility = Hepmob
	else if(specnum .eq. He2p) then
		mobility = He2pmob
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

	real*8 :: alpha,E0

	alpha     = kparams(1)
	E0        = kparams(3)

	rateconst = alpha*abs(efield)*exp(-E0/abs(efield))

end function getcustomrate
!====================================================================

end module chem_module

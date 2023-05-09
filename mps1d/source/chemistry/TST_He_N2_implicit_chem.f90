module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=10
      integer, parameter :: bgspecnum=9   !background gas
      integer, parameter :: especnum=1	  !electron species
      integer, parameter :: nreac=20
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

	integer :: E,Hep,He2p,Np,N2p
	integer :: He_star,He2_star,N,He,N2
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

	E  = 1
	Hep  = 2
	He2p   = 3
	Np  = 4
	N2p = 5
	He_star  = 6
	He2_star = 7
	N  = 8
	He  = 9
	N2 = 10

end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1)  = 'E'
	      specnames(2)  = 'He+'
	      specnames(3)  = 'He2+'
	      specnames(4)  = 'N+'
	      specnames(5)  = 'N2+'
	      specnames(6)  = 'He_star'
	      specnames(7)  = 'He2_star'
	      specnames(8)  = 'N'
	      specnames(9)  = 'He'
	      specnames(10) = 'N2'

	      molmass(1)  =  mass_elec
	      molmass(2)  = 4.d0*mass_prot
	      molmass(3)  = 8.d0*mass_prot
	      molmass(4)  = 14.d0*mass_prot
	      molmass(5)  = 28.d0*mass_prot
	      molmass(6)  = 4.d0*mass_prot
	      molmass(7)  = 8.d0*mass_prot
	      molmass(8)  = 14.d0*mass_prot
	      molmass(9)  =  4.d0*mass_prot
	      molmass(10) = 28.d0*mass_prot

	      spec_charge(1)  = -1.d0
	      spec_charge(2)  =  1.d0
	      spec_charge(3)  =  1.d0
	      spec_charge(4)  =  1.d0
	      spec_charge(5)  =  1.d0
	      spec_charge(6)  =  0.d0
	      spec_charge(7)  =  0.d0
	      spec_charge(8)  =  0.d0
	      spec_charge(9)  =  0.d0
	      spec_charge(10) =  0.d0

	      ionspecmin = 2
	      ionspecmax = 5

	      no_of_ions = 4

	      neutralspecmin = 6
	      neutralspecmax = 10

	      no_of_neutrals = 5

end subroutine setspecparams
!====================================================================

subroutine read_elec_and_gas_energies()

	integer :: i, energyinpfptr

	energyinpfptr = 230
	open(unit=energyinpfptr,file="elec_gas_energy_He_N2.inp")
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

		rrt(01) = (2.308D-10)*(TE**0.31)*EXP(-(2.297D+05)/TE)
		rrt(02) = (1.099D-11)*(TE**0.31)
		rrt(03) = (2.584D-12)*(TE**0.68)*EXP(-(2.854092D+05)/TE)
		rrt(04) = (4.661D-10)*(TE**0.6)*EXP(-(5.546D+04)/TE)
		rrt(05) = (1.268D-12)*(TE**0.71)*EXP(-(3.945D+04)/TE)
		rrt(06) = (5.386D-07)*(TE**(-0.5))
		rrt(07) = 2.7D-10
		rrt(08) = 1.3D-33
		rrt(09) = 1.0D-31
		rrt(10) = 7.0D-11
		rrt(11) = 7.0D-11
		rrt(12) = 5.0D-10
		rrt(13) = 7.0D-10
		rrt(14) = 5.0D-10
		rrt(15) = 7.0D-10
		rrt(16) = 5.651D-27*(TE**(-0.8))
		rrt(17) = 2.540D-06*(TE**(-0.5))
		rrt(18) = 1.959D-06*(TE**(-0.7))*EXP(-(1.132D+05)/TE)
		rrt(19) = 8.401D-05*EXP(-(1.682D+05)/TE)
		rrt(20) = 4.483D-07*(TE**(-0.3))*EXP(-(1.81D+05)/TE)

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
        
		rrt(01) = rrt(01) * density(01) * density(09) 
  		rrt(02) = rrt(02) * density(01) * density(06) 
  		rrt(03) = rrt(03) * density(01) * density(09) 
		rrt(04) = rrt(04) * density(01) * density(06) 
		rrt(05) = rrt(05) * density(01) * density(07) 
		rrt(06) = rrt(06) * density(01) * density(03) 
		rrt(07) = rrt(07) * density(06)**2 
		rrt(08) = rrt(08) * density(06) * density(09)**2 
		rrt(09) = rrt(09) * density(02) * density(09)**2 
		rrt(10) = rrt(10) * density(06) * density(10) 
		rrt(11) = rrt(11) * density(07) * density(10) 
		rrt(12) = rrt(12) * density(02) * density(10) 
		rrt(13) = rrt(13) * density(02) * density(10) 
		rrt(14) = rrt(14) * density(03) * density(10) 
		rrt(15) = rrt(15) * density(03) * density(10) 
		rrt(16) = rrt(16) * density(01)**2 * density(05) 
		rrt(17) = rrt(17) * density(01) * density(05) 
		rrt(18) = rrt(18) * density(01) * density(10) 
		rrt(19) = rrt(19) * density(01) * density(08) 
		rrt(20) = rrt(20) * density(01) * density(10)
		
end function zdp_getratesofprogress

subroutine zdp_getspecproduction(Te,Tg,specden,j,ydot_j) 

    	real*8, intent(in) :: Te, Tg
 		real*8, intent(in) :: specden(nspecies)
		real*8 :: rrt(nreac)
        integer :: j
		real*8, intent(inout) :: ydot_j

		rrt = zdp_getratesofprogress(Te,Tg,specden)

        ydot_j = 0.0

		ydot(01) = +rrt(03)+rrt(04)+rrt(05)-rrt(06)+rrt(07)+rrt(10)+rrt(11)-rrt(16)-rrt(17)+rrt(19)+rrt(20) 
		ydot(02) = +rrt(03)+rrt(04)+rrt(07)-rrt(09)-rrt(12)-rrt(13) 
		ydot(03) = +rrt(05)-rrt(06)+rrt(09)-rrt(14)-rrt(15) 
		ydot(04) = +rrt(13)+rrt(15)+rrt(19) 
		ydot(05) = +rrt(10)+rrt(11)+rrt(12)+rrt(14)-rrt(16)-rrt(17)+rrt(20) 
		ydot(06) = +rrt(01)-rrt(02)-rrt(04)+rrt(06)- 2.d0 * rrt(07)-rrt(08)-rrt(10) 
		ydot(07) = -rrt(05)+rrt(08)-rrt(11) 
		ydot(08) = +rrt(13)+rrt(15)+ 2.d0 * rrt(17)+ 2.d0 * rrt(18)-rrt(19) 
		ydot(09) = -rrt(01)+rrt(02)-rrt(03)+rrt(06)+rrt(07)-rrt(08)-rrt(09)+rrt(10)+& 
					2.d0 * rrt(11)+rrt(12)+rrt(13)+ 2.d0 * rrt(14) + 2.d0 * rrt(15) 
		ydot(10) = -rrt(10)-rrt(11)-rrt(12)-rrt(13)-rrt(14)-rrt(15)+rrt(16)-rrt(18)-rrt(20) 

        ydot_j = ydot(j)

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
	
	integer :: E,Hep,He2p,Np,N2p
	integer :: He_star,He2_star,N,He,N2

	E  = 1
	Hep  = 2
	He2p   = 3
	Np  = 4
	N2p = 5
	He_star  = 6
	He2_star = 7
	N  = 8
	He  = 9
	N2 = 10
	
	Patm = 101325.d0

	if(specnum .eq. E) then
		
		dcoeff = 0.1737*(Te/17406.d0)*(Patm/Pg)

	else if(specnum .eq. Hep) then
		 
		dcoeff= (5.026d-05)*(Patm/Pg)

	else if(specnum .eq. He2p) then
		 
		dcoeff= (8.148d-05)*(Patm/Pg)

	else if(specnum .eq. Np) then
		
		dcoeff= (9.710d-05)*(Patm/Pg)

	else if(specnum .eq. N2p) then
		 
		dcoeff= (1.015d-04)*(Patm/Pg)

	else if(specnum .eq. He_star) then
		 
		dcoeff= (4.116d-04)*(Patm/Pg)

	else if(specnum .eq. He2_star) then
		
		dcoeff= (2.029d-04)*(Patm/Pg)
	
	else if(specnum .eq. N) then
		
		dcoeff= (1.955d-04)*(Patm/Pg)

	else if(specnum .eq. He) then
		
		dcoeff= (4.116d-04)*(Patm/Pg)

	else if(specnum .eq. N2) then
		
		dcoeff= (1.075d-04)*(Patm/Pg)

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

	integer :: E,Hep,He2p,Np,N2p
	integer :: He_star,He2_star,N,He,N2

	real*8 :: Te_in_eV
	real*8 :: log_Te
	
	E  = 1
	Hep = 2
	He2p = 3
	Np  = 4
	N2p = 5

	neutral_den = Pg/k_B/Tg
	Te_in_eV = Te/eVtoK
	log_Te   = log(Te_in_eV)
	
	Patm = 101325.d0
	Troom = 300.d0
	Natm = Patm/k_B/Troom

	if(specnum .eq. E) then
 
		mobility = -0.11320

	else if(specnum .eq. Hep) then
		
		mobility   = 1.482d-03

	else if(specnum .eq. He2p) then

		mobility   = 2.403d-03

	else if(specnum .eq. Np) then

		mobility   = 2.863d-03

	else if(specnum .eq. N2p) then
	
		mobility   = 2.993d-03

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

module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=3
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=2	  !electron species
      integer, parameter :: nreac=1
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
      integer :: solved_ions_num
      integer :: solved_neutrals_num
      integer :: no_of_ions
      integer :: no_of_neutrals
      
      real*8 :: k_BMAN  

      contains
!====================================================================
subroutine assignreactions()

	integer :: E,AR,ARp
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

        AR = 1
        E = 2
        ARp = 3        

	!convention for energy is added is positive and
	!removed is negative.
	!eg: electron loses 19.8 eV when reaction 
	!1 happens

	rnum=1
	!***********************************
	!E + AR --> ARp + 2E
	!***********************************
	reactants(rnum,AR) = 1
	reactants(rnum,E) = 1
	
	products(rnum,ARp)  = 1	
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.5d-6*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 278508.0
	!-----------------------------------
	elecenergy(rnum) = -15.578
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	real*8 :: M_PROT, M_ELEC
        M_PROT = mass_prot
        M_ELEC = mass_elec      

              specnames(1)  = 'AR'
	      specnames(2)  = 'E'
              specnames(3)  = 'ARp'

	      molmass(1)  = 40.d0*M_PROT
	      molmass(2)  = M_ELEC 
              molmass(3)  = 40.d0*M_PROT

	      spec_charge(1)  =  0.d0
	      spec_charge(2)  = -1.d0
	      spec_charge(3)  =  1.d0

	      ionspecmin = 3
	      ionspecmax = 3

	      solved_ions_num = 1

	      neutralspecmin = 1 
	      neutralspecmax = 1

	      solved_neutrals_num = 1

              no_of_ions = solved_ions_num
              no_of_neutrals = solved_neutrals_num

end subroutine setspecparams
!====================================================================
subroutine initializechemistry()

	call setspecparams()
	call assignreactions()

        
end subroutine initializechemistry
!====================================================================

! TST

function zdp_getreactionrates(Te, Tgas) result(rrt)

        real*8 :: rrt(nreac)
        real*8, intent(in) :: Te, Tgas

        !write(*,*) "zdp_getreactionrates() called"

        rrt(1) = 2.5D-06*(Te**0.0d0)*EXP(-278508.0/Te)

        rrt = rrt*cm_to_m**3
!        write(*,*) "rrt = ",  rrt(1)

end function

subroutine zdp_getspecproduction(Te,Tg,specden) 

        !integer, intent(inout) :: specnum
        real*8, intent(in) :: Te, Tg
        real*8, intent(in) :: specden(nspecies)
        real*8 :: density(nspecies)
        real*8 :: rrt(nreac)
        !real*8 :: ydot(nspecies)        

        rrt = zdp_getreactionrates(Te,Tg)

        density = specden

        rrt(1) = rrt(1) * density(1) * density(2)
        ydot(1) = -rrt(1)
        ydot(2) = +rrt(1)
        ydot(3) = +rrt(1)

       !write(*,*) "density = ", density

end subroutine zdp_getspecproduction

subroutine getspecproduction(specnum,Te,Tg,specden,specprod,efld)

	integer :: specnum
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

        if (specnum .eq. especnum) then
                write(*,*) "specprod in getspecproduction = ", specprod
        endif
        

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
        real*8 :: k_BMAN
	real*8 :: dcoeff
	real*8 :: Patm,mob
	
	integer :: E,AR,ARp
	
        AR  = 1
        E   = 2
        ARp = 3

	Patm = 101325.d0
        k_BMAN = k_B
	if(specnum .eq. E) then
		
                dcoeff=100.d0

	else if(specnum .eq. ARp) then
		
		dcoeff=0.01
        else if(specnum .eq. AR) then

                dcoeff=0.01
	else
		write(*,*)"dcoeff - species does not exist"
		stop
	endif

end function getspecdcoeff


!====================================================================
function getspecmobility(specnum,specarray,elecfield,Te,Tg,Pg)  result(mobility)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8, intent(in)   :: specarray(nspecies)
	real*8, intent(in)   :: elecfield

	real*8 :: mobility,mob_inv_AR
	real*8 :: Patm,Natm,Troom
	real*8 :: neutral_den

	integer :: E,AR,ARp
	real*8 :: Te_in_eV
	real*8 :: log_Te
	
        AR = 1
        E = 2
        ARp = 3

        k_BMAN = k_B
	neutral_den = Pg/k_BMAN/Tg
	Te_in_eV = Te/eVtoK
	log_Te   = log(Te_in_eV)

	Patm = 101325.d0
	Troom = 300.d0
	Natm = Patm/k_BMAN/Troom

	if(specnum .eq. E) then
                mobility = -20.d0
        else if(specnum .eq. ARp) then
                mobility   = 0.2
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

	collfreq = -ECHARGE/m_e/mu

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
!        write(*,*) "rateconstant = ", rateconst

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

        rateconst = 0.0

end function getcustomrate
!====================================================================

end module chem_module

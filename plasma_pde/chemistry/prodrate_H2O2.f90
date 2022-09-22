program chem_driver

      use chem_module
	implicit none
	
      real*8 :: specdenvec(nspecies);
      real*8 :: specprod,inelterm,emob
      real*8 :: efield,etemp,gtemp
      integer :: H2,E,Hp,H2p,H
      integer :: O2,Op,O2p,Om,O

      print *,"nspecies:",nspecies

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

      specdenvec=1.d12

!typical values atm. pressure streamer======
      specdenvec(H2)=1.63d25
      specdenvec(O2)=8.1d24
      specdenvec(E)=5.0d19
      
      specdenvec(Hp)=2.13d18
      specdenvec(H2p)=4.8d19
      specdenvec(Op)=4.8d14

      specdenvec(O2p) =6.39d18
      specdenvec(Om) =6.6d18
      specdenvec(H) =4.07d20
      specdenvec(O)=3.2d20
     
      efield = 2.77d5
      etemp  = 5.18d4
      gtemp  = 300.d0
!============================================

!typical values H2_O2 test case==============
      specdenvec(H2)=2.0d24
      specdenvec(O2)=1.0d24
      specdenvec(E)=3.24d14
      
      specdenvec(Hp)=3.17d13
      specdenvec(H2p)=2.27d14
      specdenvec(Op)=2.5d13

      specdenvec(O2p) = 7.4d13
      specdenvec(Om)  = 3.13d13
      specdenvec(H)   = 1.3d15
      specdenvec(O)   = 4.82d14

      efield = 437439.d0
      etemp  = 54292.1
      gtemp  = 300.d0
!============================================

      call initializechemistry()
      
      call getspecproduction(E,etemp,gtemp,specdenvec,specprod,efield);
      write(*,'(A E20.10)')"specprod E:",specprod
      
      call getspecproduction(H,etemp,gtemp,specdenvec,specprod,efield);
      write(*,'(A E20.10)')"specprod H:",specprod
      
      call getspecproduction(O,etemp,gtemp,specdenvec,specprod,efield);
      write(*,'(A E20.10)')"specprod O:",specprod

      !print *,"dcoeff :",getspecdcoeff(4,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(5,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(6,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(7,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(8,1.d0*11604.5,300.d0,13300.d0)

end program chem_driver

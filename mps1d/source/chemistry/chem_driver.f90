program chem_driver

      use chem_module
	implicit none
	
      real*8 :: specdenvec(nspecies);
      real*8 :: specprod,inelterm,emob

      print *,"nspecies:",nspecies

      specdenvec=1.d12

      specdenvec(1)=2.d22
      specdenvec(2)=1.d22
      specdenvec(3)=1.d14

      call initializechemistry()
      call getspecproduction(3,11604.d0,300.d0,specdenvec,specprod,2.d7);
      call getelectroninelasticterm(11604.d0,300.d0,specdenvec,inelterm,2.d7);

      emob=getspecmobility(2,specdenvec,2.d7,5.d0*11604.5,300.d0,13300.d0)
      print *," mobility is ",emob
      write(*,'(A E20.10)')"specprod:",specprod
      write(*,'(A E20.10)')"inelterm:",inelterm

      print *,"dcoeff :",getspecdcoeff(2,specdenvec,2.d7,&
		      5.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(4,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(5,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(6,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(7,1.d0*11604.5,300.d0,13300.d0)
      !print *,"dcoeff :",getspecdcoeff(8,1.d0*11604.5,300.d0,13300.d0)

end program chem_driver

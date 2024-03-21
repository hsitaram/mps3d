module convdiff
use vectorfunctions
use solvergmres_module
implicit none

integer,parameter :: minsize=5

contains
!===============================================================
subroutine find_fo_upwind_flux(cL,cR,uL,uR,flux)

	real*8,intent(in) :: cL,cR
	real*8,intent(in) :: uL,uR
	real*8,intent(out) :: flux

	real*8 :: c_half

	c_half=0.5*(cL+cR)
	
	if(c_half .ge. 0.d0) then
	    flux=cL*uL
	else
            flux=cR*uR
	endif
	     
end subroutine find_fo_upwind_flux
!===============================================================
subroutine find_so_WAF_flux(cLm1,cL,cR,cRp1,uLm1,uL,uR,uRp1,dt,dx,flux)

	real*8,intent(in) :: cLm1,cL,cR,cRp1
	real*8,intent(in) :: uLm1,uL,uR,uRp1
        real*8,intent(in) :: dt,dx
	real*8,intent(out) :: flux

	real*8 :: r, c_half,c
        real*8 :: eps,lim,sign_c

        eps=1.d-10
	c_half=0.5*(cL+cR)
        c=c_half*dt/dx

        if(c .le. 0) then
          sign_c=-1.d0
        else
          sign_c=1.d0
        endif

        if(c>=0.d0) then
           r=(uL-uLm1)/(uR-uL+eps)
        else
           r=(uRp1-uR)/(uR-uL+eps)
        endif

        lim=0.d0
        if(r<0.d0) then
           lim=1.d0
        else if((r .ge. 0.d0) .and. (r .lt. 1)) then
           lim=1.d0-(1.d0-abs(c))*r
        else
           lim=abs(c)
        endif
        
        !print *,"lim:",lim,r
        flux=0.5*(cL*uL+cR*uR)-0.5*sign_c*lim*(cR*uR-cL*uL)
        !flux=0.5*(cL*uL+cR*uR)

end subroutine find_so_WAF_flux
!===============================================================
subroutine slope_limited_reconstruct_lr(phi_im1,phi_i,phi_ip1,phil,phir)
   
   real*8, intent(in) :: phi_im1,phi_i,phi_ip1
   real*8, intent(out) :: phil,phir
   real*8 :: beta,del_imhalf,del_iphalf,del_i

   beta=1.0 !MINBEE
   del_imhalf=phi_i-phi_im1
   del_iphalf=phi_ip1-phi_i

   if(del_iphalf>0.d0) then
      del_i=max(0.d0,min(beta*del_imhalf,del_iphalf),min(del_imhalf,beta*del_iphalf))
   else
     del_i=min(0.d0,max(beta*del_imhalf,del_iphalf),max(del_imhalf,beta*del_iphalf))
   endif

   phil=phi_i-0.5d0*del_i
   phir=phi_i+0.5d0*del_i
   
end subroutine slope_limited_reconstruct_lr
!===============================================================
subroutine find_so_MH_flux(cim1,ci,cip1,cip2,uim1,ui,uip1,uip2,dt,dx,flux)

        !we are looking at face (i+1/2)
        !c is velocity, u is state
	real*8,intent(in) :: cim1,ci,cip1,cip2
	real*8,intent(in) :: uim1,ui,uip1,uip2
        real*8,intent(in) :: dt,dx
	real*8,intent(out) :: flux

        real*8 :: vel,u_i_L,u_i_R
        real*8 :: u_ip1_L,u_ip1_R
        real*8 :: c_i_L,c_i_R
        real*8 :: c_ip1_L,c_ip1_R
        real*8 :: uLbar,uRbar

        real*8 :: flx_i_L,flx_i_R
        real*8 :: flx_ip1_L,flx_ip1_R
        real*8 :: sign_vel

        vel=0.5d0*(ci+cip1)
        if(vel .le. 0) then
          sign_vel=-1.d0
        else
          sign_vel=1.d0
        endif
        call slope_limited_reconstruct_lr(cim1,ci,cip1,c_i_L,c_i_R)
        call slope_limited_reconstruct_lr(ci,cip1,cip2,c_ip1_L,c_ip1_R)

        call slope_limited_reconstruct_lr(uim1,ui,uip1,u_i_L,u_i_R)
        call slope_limited_reconstruct_lr(ui,uip1,uip2,u_ip1_L,u_ip1_R)
        call slope_limited_reconstruct_lr(uim1*cim1,ui*ci,uip1*cip1,flx_i_L,flx_i_R)
        call slope_limited_reconstruct_lr(ui*ci,uip1*cip1,uip2*cip2,flx_ip1_L,flx_ip1_R)

        !half time level reconstruction
        uLbar=u_i_R+0.5*dt/dx*(flx_i_L-flx_i_R)
        uRbar=u_ip1_L+0.5*dt/dx*(flx_ip1_L-flx_ip1_R)

        flux=uLbar*c_i_R*0.5*(1.d0+sign_vel)+&
                uRbar*c_ip1_L*0.5*(1.d0-sign_vel)


end subroutine find_so_MH_flux
!===============================================================
subroutine finddiffusiveflux(DL,DR,uL,uR,dx,flux)

	real*8,intent(in) :: DL,DR
	real*8,intent(in) :: uL,uR
	real*8,intent(in) :: dx
	real*8,intent(out) :: flux

	real*8 :: D_half
        real*8, parameter :: eps=1e-15

	!D_half = HALF*(DL+DR)
	D_half = 2.d0*DL*DR/(DL+DR+eps);
	flux = D_half*(uR-uL)/dx

end subroutine finddiffusiveflux
!===============================================================
subroutine findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n,hypfluxorder)

      integer, intent(in) :: n
      real*8, intent(inout) :: AX(n)
      real*8, intent(in)  :: X(n)
      integer, intent(in) :: hypfluxorder

      real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
      real*8,intent(in)  :: dx,dt
      logical,intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8, intent(in) :: timederivfactor
      
      integer :: i
      real*8 :: flux
      real*8 :: dx2

      dx2 = dx*dx

      AX(:) = timederivfactor*X(:)/dt
      
      !convection and diffusion terms
      do i=1,n-1

      if(hypfluxorder==1 .or. timederivfactor==0.d0) then
              call find_fo_upwind_flux(vel(i),vel(i+1),X(i),X(i+1),flux)
      else if(hypfluxorder==2) then
              print *,"calling waf"
              if(i==1) then
                      call find_so_WAF_flux(vel(i),vel(i),vel(i+1),vel(i+2),X(i),X(i),X(i+1),X(i+2),&
                      dt,dx,flux)
              else if(i==n-1) then
                      call find_so_WAF_flux(vel(i-1),vel(i),vel(i+1),vel(i+1),X(i-1),X(i),X(i+1),X(i+1),&
                      dt,dx,flux)
              else 
                      !call find_so_MH_flux(vel(i-1),vel(i),vel(i+1),vel(i+2),X(i-1),X(i),X(i+1),X(i+2),&
                      !        timederivfactor*dt,dx,flux)
                      call find_so_WAF_flux(vel(i-1),vel(i),vel(i+1),vel(i+2),X(i-1),X(i),X(i+1),X(i+2),&
                      dt,dx,flux)
              endif

      else
              print *,"Hyporder not implemented"
              exit
      endif

      AX(i)   = AX(i)   + flux/dx
      AX(i+1) = AX(i+1) - flux/dx
      call finddiffusiveflux(dcoeff(i),dcoeff(i+1),X(i),X(i+1),dx,flux)
      AX(i)   = AX(i)   - flux/dx
      AX(i+1) = AX(i+1) + flux/dx
      enddo

      !print *,"AX inside findAX:",AX

      !reaction term
      do i=1,n
      AX(i) = AX(i) - reac(i)*X(i)
      enddo

      !boundary conditions
      if(dirc_bc_flags(1) .eqv. .true.) then
              AX(1)=X(1)
      endif
      if(dirc_bc_flags(2) .eqv. .true.) then
              AX(n)=X(n)
      endif

end subroutine findAX
!===============================================================
subroutine findrhs(b,xold,timederivfactor,&
                source,dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n)

        integer, intent(in) :: n
        real*8, intent(inout) :: b(n)
        real*8, intent(in) :: xold(n),source(n)
        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8, intent(in) :: dx,dt
        real*8, intent(in) :: timederivfactor
        integer :: i

        do i=1,n
        b(i) = timederivfactor*xold(i)/dt + source(i)
        enddo	

        if(dirc_bc_flags(1) .eqv. .true.) then
                b(1) = dircvals(1)
        endif
        if(dirc_bc_flags(2) .eqv. .true.) then
                b(n) = dircvals(2)
        endif
        if(flux_bc_flags(1) .eqv. .true.) then
                b(1) = b(1) + fluxvals(1)/(dx)
        endif
        if(flux_bc_flags(2) .eqv. .true.) then
                b(n) = b(n) - fluxvals(2)/(dx)
        endif

end subroutine findrhs
!===============================================================
subroutine noprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,&
                dx,dt,n)
        integer, intent(in) :: n
        real*8, intent(inout) :: MinvX(n)
        real*8, intent(in)  :: X(n)

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
        real*8,intent(in)  :: timederivfactor

        MinvX = X

end subroutine noprecond
!===============================================================
subroutine gauss_seidel_smoothing(res,b,X,timederivfactor,vel,dcoeff,reac,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n,maxiter,tol,hypfluxorder)

        integer, intent(in)   :: n
        real*8, intent(in)    :: b(n)
        real*8, intent(inout) :: X(n)
        real*8, intent(inout) :: res(n)
        real*8, intent(in)    :: tol
        integer, intent(in) :: maxiter
        integer, intent(in) :: hypfluxorder

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
        real*8,intent(in)  :: timederivfactor

        integer :: i,it

        real*8 :: diag
        real*8 :: cL,cR,chalf
        real*8 :: dL,dR,dhalf
        real*8 :: offdiag
        real*8 :: AX(n)
        real*8 :: dx2
        real*8 :: resnorm

        dx2 = dx*dx
        AX  = ZERO

        !gauss seidel iterations
        do it=1,maxiter

        do i=1,n

        diag = -reac(i) + timederivfactor*ONE/dt

        offdiag=ZERO

        !right face
        if(i .lt. n) then
                !convection term
                cL = vel(i)
                cR = vel(i+1)
                chalf = HALF*(cL+cR)
                if(chalf .ge. 0) then
                        diag = diag + chalf/dx
                else
                        offdiag = offdiag + chalf*X(i+1)/dx
                endif

                !diffusion term
                dL = dcoeff(i)
                dR = dcoeff(i+1)
                dhalf = HALF*(dL + dR)

                diag = diag + dhalf/dx2
                offdiag = offdiag - X(i+1)*dhalf/dx2
        endif

        !left face
        if(i .gt. 1) then

                !convection term
                cL = vel(i-1)
                cR = vel(i)
                chalf = HALF*(cL+cR)
                if(chalf .ge. 0) then
                        offdiag = offdiag - chalf*X(i-1)/dx
                else
                        diag = diag - chalf/dx
                endif

                !diffusion term
                dL = dcoeff(i-1)
                dR = dcoeff(i)
                dhalf = HALF*(dL+dR)

                diag = diag + dhalf/dx2
                offdiag = offdiag - X(i-1)*dhalf/dx2
        endif

        if(i .eq. 1) then
                if(dirc_bc_flags(1) .eqv. .true.) then
                        diag    = ONE
                        offdiag = ZERO
                endif
        endif	

        if(i .eq. n) then
                if(dirc_bc_flags(2) .eqv. .true.) then
                        diag    = ONE
                        offdiag = ZERO
                endif
        endif

        X(i) = (b(i)-offdiag)/diag

        enddo

        call  findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,n,hypfluxorder)

        res = b - AX

        call findnorm(resnorm,res,n)

        if(resnorm .lt. tol) then
                !print *,"it:",it
                exit
        endif 

        enddo

        call  findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,dx,dt,n,hypfluxorder)

        res = b - AX 

end subroutine gauss_seidel_smoothing
!===============================================================
subroutine restriction(Xh,X2h,n2h) 

        integer, intent(in) :: n2h
        real*8, intent(inout) :: Xh(2*n2h-1)
        real*8, intent(inout) :: X2h(n2h)

        integer :: i

        X2h(1)     = Xh(1)
        X2h(n2h) = Xh(2*n2h-1)

        do i=2,n2h-1
        X2h(i)=ONEFOURTH*(Xh(2*i-2) + TWO*Xh(2*i-1) + Xh(2*i))
        enddo


end subroutine restriction
!===============================================================
subroutine prolong(Xh,X2h,n2h)

        integer, intent(in) :: n2h
        real*8, intent(inout) :: Xh(2*n2h-1)
        real*8, intent(inout) :: X2h(n2h)

        integer :: i

        do i=1,n2h-1
        Xh(2*i-1) = X2h(i)
        Xh(2*i)   = HALF*(X2h(i) + X2h(i+1))
        enddo

        Xh(2*n2h-1) = X2h(n2h)

end subroutine prolong
!===============================================================
recursive subroutine dovcycle(X,b,timederivfactor,vel,dcoeff,reac,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n,hypfluxorder)

        integer, intent(in)   :: n
        real*8, intent(inout) :: X(n)
        real*8, intent(inout) :: b(n)
        integer, intent(in)   :: hypfluxorder

        logical, intent(in)  :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in)   :: dircvals(2),fluxvals(2)
        real*8, intent(in)   :: vel(n),dcoeff(n),reac(n)
        real*8 ,intent(in)   :: dx,dt
        real*8 ,intent(in)   :: timederivfactor

        real*8  :: resh(n)
        real*8 :: eh(n)
        real*8 :: velh(n)
        real*8 :: dcoeffh(n)
        real*8 :: reach(n)

        real*8 :: res2h(n/2+1)
        real*8 :: e2h(n/2+1)
        real*8 :: vel2h(n/2+1)
        real*8 :: dcoeff2h(n/2+1)
        real*8 :: reac2h(n/2+1)
        real*8 :: Xmin(minsize)
        real*8 :: AX(minsize)
        logical :: success
        logical :: printflag
        real*8 :: initial_res

        success = .false.
        printflag = .false.


        velh    = vel
        dcoeffh = dcoeff
        reach   = reac

        resh     = ZERO
        eh       = ZERO

        res2h    = ZERO
        e2h      = ZERO
        vel2h    = ZERO
        dcoeff2h = ZERO
        reac2h   = ZERO
        Xmin     = ZERO

        if(n .le. minsize) then
                !call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
                !		dirc_bc_flags,flux_bc_flags,dircvals,&
                !		fluxvals,dx,dt,n,10000,NEARZERO)

                call performgmres(b,Xmin,X,timederivfactor,&
                        velh,dcoeffh,reach,dirc_bc_flags,&
                        flux_bc_flags,dircvals,fluxvals,dx,dt,&
                        n,n,1,findAX,noprecond,&
                        1.d-12,success,printflag,initial_res,hypfluxorder)

                AX = ZERO

                call  findAX(AX,X,timederivfactor,velh,dcoeffh,reach,dirc_bc_flags,&
                        flux_bc_flags,dircvals,fluxvals,dx,dt,n,hypfluxorder)

                resh = b - AX 

        else

                !initial smoothing
                call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
                        dirc_bc_flags,flux_bc_flags,dircvals,&
                        fluxvals,dx,dt,n,1,NEARZERO,hypfluxorder)

                !restriction of residual from fine to coarse grid
                call restriction(resh,res2h,n/2+1)

                call restriction(reach,reac2h,n/2+1)
                call restriction(dcoeffh,dcoeff2h,n/2+1)
                call restriction(velh,vel2h,n/2+1)

                call dovcycle(e2h,res2h,timederivfactor,vel2h,dcoeff2h,reac2h,&
                        dirc_bc_flags,flux_bc_flags,dircvals,&
                        fluxvals,TWO*dx,dt,n/2+1,hypfluxorder)

                !prolong error from coarse to fine grid
                call prolong(eh,e2h,n/2+1)

                !update
                X(:) = X(:) + eh(:)

                !post smooth
                call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
                        dirc_bc_flags,flux_bc_flags,dircvals,&
                        fluxvals,dx,dt,n,1,NEARZERO,hypfluxorder)

        endif

end subroutine dovcycle
!===============================================================
subroutine mgridprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
                flux_bc_flags,dircvals,fluxvals,&
                dx,dt,n,hypfluxorder)

        integer, intent(in)    :: n
        real*8, intent(inout)    :: MinvX(n)
        real*8, intent(inout)  :: X(n)
        integer, intent(in)    :: hypfluxorder

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
        real*8,intent(in)  :: timederivfactor

        integer :: nvcycles,i
        !real*8 :: AX(n),res(n)
        !real*8 :: norm
        !integer :: j

        !norm=0.d0
        nvcycles=1

        do i=1,nvcycles
        call dovcycle(MinvX,X,timederivfactor,vel,dcoeff,reac,&
                dirc_bc_flags,flux_bc_flags,dircvals,&
                fluxvals,dx,dt,n,hypfluxorder)
        enddo


end subroutine mgridprecond
!===============================================================

end module convdiff

	Program mcmc
	Implicit none
	!-------------------------------------------------------------------------------------------------------
	! Markov Chain Monte carlo code to do parameter estimation.
	! Gopolang Mohlabeng
	! Computational Physics, Project 5
	!-------------------------------------------------------------------------------------------------------
	Integer			:: i, j, k, ndat, npar, n
	parameter(ndat=101, npar=2, n=50000)
	Real*8			:: x(ndat), y(ndat), errd(ndat), t(ndat)
	Real*8, Dimension(npar)	:: A, Aold, Anew, chain
	Real, Dimension(npar)	:: z1, z2 
	Integer			:: idun
	Logical			:: accept = .False., reject = .False.
	Real			:: randnum, ran2
	Real*8			:: oldchi, newchi, RR, alpha, stepsize, likeli, pi, P, eqold, eqnew
	Real*8			:: lnlike(n), L(n), cor1, cor2, r, A1(n), A2(n), likelihood(n), prob

	External prob
	pi = 4.d0*Atan(1.d0)

! stepsize is very important because it determines the acceptance ratio or rather is determined by the ratio.

	OPEN(UNIT=54,FILE='measurements.dat',STATUS='UNKNOWN',form='formatted') 
 	do i=1,ndat
 		read(54,*) x(i), y(i), errd(i)	
 	enddo  
 	close(54)

! These are the values of the initial parameter values
	chain(1) = 2.0d0
	chain(2) = 0.1d0
	    idun = -2

! Choose the initial parameter vector
	do i = 1, npar
	   A(i) = chain(i)
	enddo
	stepsize = 0.001d0
	oldchi = newchi

	OPEN(UNIT=38,FILE='posterior.dat',STATUS='unknown')
! Doing the Monte Carlo loop with large n
	do i = 1, n
	   do j = 1, npar
	      call GRNF(z1(j),z2(j))   ! Gaussian Random numbers
	      Anew(j) = A(j) + stepsize*z1(j)
	   enddo
	   oldchi = 0.d0	! chi^2 = -2log(Likelihhod)
	   lnlike(i) = 0.d0
	   do k = 1, ndat
	      t(k) = Anew(1) + Anew(2)*(x(k)*x(k)) !this is the model
	      oldchi = oldchi + ((y(k) - t(k))**2)/(errd(k)*errd(k))
              likeli = prob(y(k), t(k), errd(k)) ! Posterior density which is the likelihood
	      lnlike(i) = lnlike(i) + Log(likeli) ! Log P
	   enddo
	   if(i .ne. 1) then 
	     RR = Exp(oldchi - newchi)
	     alpha = Exp(lnlike(i) - lnlike(i-1)) ! the acceptance ratio
	     randnum = ran2(idun)
	     ! The acceptance/ rejection criterion
	     if(alpha .ge. 1.d0) then
	       A = Anew
	     else if (randnum .lt. alpha) then  
	       A = Anew
	     else
	       Anew = Anew 
	    endif
	  endif
	  L(i)  = Exp(lnlike(i))/(2.5d+105) !Likelihood with normalization
	  write(38,*), Anew, L(i)
	enddo	
	print*, 'The values of the parameters are, a1 = ', Anew(1), 'and a2 = ', Anew(2)
	print*, 'The value of the chi^2 is', oldchi
	print*, 'The acceptance ratio is', alpha
	close(38)	
	

	End Program mcmc

	Real*8 Function prob(dat, model, errors) ! Normal pdf P(D|x) centered on the model parameters
	Implicit none
	Real*8			:: dat, model, errors, pi
	pi = 4.d0*Atan(1.d0)
	prob = (1.d0/sqrt(2.d0*pi*(errors*errors)))*exp((-1.d0/2.d0)*((dat - model)**2)/(errors**2))	
	End Function


!	Computing the Gaussian Random Numbers
	SUBROUTINE GRNF (X,Y)
! 	Two Gaussian random numbers generated from two uniform random
! 	numbers. Copyright (c) Tao Pang 1997.
!
  	Implicit none
  	Real, Intent (out) :: X,Y
	Integer :: idun
  	Real :: Pi,r1,r2,r, ran2
	External ran2
!
  	Pi = 4.0*Atan(1.0)
  	r1 = -Alog(1.0-ran2(idun))
  	r2 = 2.0*Pi*ran2(idun)
  	r1 = SQRT(2.0*r1)
  	X  = r1*COS(r2)
  	Y  = r1*SIN(r2)
	END SUBROUTINE GRNF

	FUNCTION ran2(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL ran2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
           idum=max(-idum,1)
           idum2=idum
           do 11 j=NTAB+8,1,-1
              k=idum/IQ1
              idum=IA1*(idum-k*IQ1)-k*IR1
              if (idum.lt.0) idum=idum+IM1
                if (j.le.NTAB) iv(j)=idum
11         continue
           iy=iv(1)
         endif
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
            k=idum2/IQ2
            idum2=IA2*(idum2-k*IQ2)-k*IR2
         if (idum2.lt.0) idum2=idum2+IM2
            j=1+iy/NDIV
            iy=iv(j)-idum2
            iv(j)=idum
         if(iy.lt.1)iy=iy+IMM1
            ran2=min(AM*iy,RNMX)
         return
         END







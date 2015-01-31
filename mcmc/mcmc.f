c-------------------------------------------------------------------------------------------------c
      program mcmc
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This program implements a simple monte carlo markov chain to find the best fit to a random     c
c  set of data.                                                                                   d
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c parameters used in this program only
      integer ndata, nsteps, naccept, nreject, i, j
      double precision x(101), y(101), err(101), test(101)
      double precision chisq, old_chisq, min_chisq, param(2), new_param(2), best_param(2)
      double precision stepsize, ratio, prob, R, L
      double precision x_ave, y_ave, xsq_ave, xy_ave, m_exp, b_exp
      
c initialize various parameters
      idum = -42
      ndata = 101
      nsteps = 100000
      naccept = 1
      nreject = 0
      stepsize = 0.01d0
      
c read in the data from file
      open(unit=42, file='measurements.txt')
      do i=1,ndata
        read(42,*) x(i), y(i), err(i)
        err(i) = 0.1d0
      enddo
      close(42)
      
c open up the file that will have the posterior distribution
      open (unit=42,file='likelihood.txt')

c set the initial parameters the we are trying to find
      param(1) = 2.d0
      param(2) = 2.d0
      
c calculate the new chisq at this point
      old_chisq = 0.d0
      do j=1,ndata
        test(j) = param(1) + param(2)*x(j)**2
        old_chisq = old_chisq + (y(j)-test(j))**2/(err(j)**2)
      enddo
      min_chisq = old_chisq
        
c loop over the steps of the markov chain
      do i=1,nsteps

c get a new set of parameters with a random gaussian step
        new_param(1) = randgauss(idum,param(1),stepsize)
        new_param(2) = randgauss(idum,param(2),stepsize)
        
c calculate the new chisq at this point
        chisq = 0.d0
        do j=1,ndata
          test(j) = new_param(1) + new_param(2)*x(j)**2
          chisq = chisq + (y(j)-test(j))**2/(err(j)**2)
        enddo
        if (chisq.lt.min_chisq) then
          min_chisq = chisq
          best_param(1) = new_param(1)
          best_param(2) = new_param(2)
        endif

c calculate the acceptance ratio and liklihood        
        R = dexp(0.5d0*(old_chisq - chisq))
        L = dexp(-0.5d0*chisq)
        prob = random(idum)
        
        if (prob .le. R) then
          write(42,*) param(1), param(2), L, R
c          write(*,*) chisq/dble(ndata-2),param(1),param(2)
          param(1) = new_param(1)
          param(2) = new_param(2)
          old_chisq = chisq
          naccept = naccept + 1
        else
          nreject = nreject + 1
        endif
      enddo
      close(42)
      
      ratio = dble(naccept)/dble(nreject)

      write(*,*) 'best chisq = ',min_chisq
      write(*,*) 'theta 1 = ',best_param(1),' theta 2 = ',best_param(2)
      write(*,*) 'acceptance ratio = ',ratio

c calculate the expected value
      x_ave = 0.d0
      y_ave = 0.d0
      xsq_ave = 0.d0
      xy_ave = 0.d0
      do i=1,ndata
        x_ave = x_ave + x(i)
        y_ave = y_ave + y(i)
        xsq_ave = xsq_ave + x(i)**2
        xy_ave = xy_ave + x(i)*y(i)
      enddo
      x_ave = x_ave/dble(ndata)
      y_ave = y_ave/dble(ndata)
      xsq_ave = xsq_ave/dble(ndata)
      xy_ave = xy_ave/dble(ndata)

      m_exp = (xy_ave - x_ave*y_ave)/(xsq_ave - x_ave**2)
      b_exp = (xsq_ave*y_ave - xy_ave*x_ave)/(xsq_ave - x_ave**2)
      
      write(*,*) 'expected value for y-intercept and slope'
      write(*,*) b_exp, m_exp
      end
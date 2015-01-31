c-------------------------------------------------------------------------------------------------c
      subroutine binomial_fit(input_data, n, best_p, func_min, fileunit)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine will take an input set of discrete, normalized probabilities and find the      c
c  value, p, for the binomial distribution that best fits the data via a least squares analysis.  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer*8 n, fileunit
      double precision input_data(101,2), best_p, func_min
      
c parameters used in this subroutine only
      double precision binomial_function
      external binomial_function

c common block for passing the data to the function
      double precision pass_data(0:100,2)
      common/datapass/ pass_data
      integer*8 n_pass
      common/npass/ n_pass
      
c parameters for brent
      double precision ax, bx, cx, tol

c check to see if the data was properly sent
      write(*,*) n
      do x1=1,n
        write(*,*) x1, input_data(x1,1), input_data(x1,2)
      enddo
      read(*,*)

c reset the array that's passed to the minimization function
      do x1=0,100
        pass_data(x1) = 0.d0
      enddo
      
      n_pass = n-1
      do x1=0,n_pass
        pass_data(x1) = input_data(x1+1,2)
      enddo

c initalize parameters for the function minimizing subroutine brent      
      ax = -1.d0
      bx = 0.5d0
      cx = 2.d0
      tol = 1.0d-10
      func_min = brent(ax,bx,cx,binomial_function,tol,best_p)

c if the fileunit is not equal to zero then write out the ourput file
      do x1=1,n
        write(fileunit,*) x1, input_data(x1,1), input_data(x1,2)
      enddo
      
      return
      end
      


c-------------------------------------------------------------------------------------------------c
      function binomial_function(p_test)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function is the least squares function to find the best p in the binomial distribution.   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters and function definition
      double precision p_test, binomial_function
      
c common block with the data
      double precision pass_data(0:100)
      common/datapass/ pass_data
      integer*8 n_pass
      common/npass/ n_pass
      
c parameters used in this function only
      double precision binomial_probability
      
c initialize the function
      binomial_function = 0.d0
      
      do x1=0,n_pass-1
        binomial_function = binomial_function + (binomial_probability(p_test,n_pass,x1) - pass_data(x1))**2
      enddo
      
      return
      end



c-------------------------------------------------------------------------------------------------c
      function binomial_probability(p, n, r)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function calculates the binomial distribution probability given the p value,              c
c  the n number of things taken r at a time.                                                      c 
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      
c input parameters and function definition
      double precision p, binomial_probability
      integer*8 n, r

c parameters used in this function only
      integer*8 r_use, x1
      
      binomial_probability = 1.d0
      
      if (r.le.n/2) then
        r_use = r
      else
        r_use = n-r
      endif
      
      do x1=1,r_use
        binomial_probability = binomial_probability*dble(n-x1+1)/dble(x1)
      enddo
      
      binomial_probability = binomial_probability*(p**r)*((1.d0-p)**(n-r))

      return
      end
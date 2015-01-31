c-------------------------------------------------------------------------------------------------c
      program binomial_driver
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This program tests the binomial fit subroutine.                                                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c parameters used in this program only
      double precision datas(16), p_test, p_result, binomial_probability
      integer*8 n, nm1, counter

c parameters for testing
      integer*8 test_r, test_n
      double precision test_p
      
c-------------------------------------------------------------------------------------------------c
c test the binomial probability function with variables
      test_r = 3
      test_n = 10
      test_p = 0.5d0
      write(*,*) binomial_probability(test_p,test_n,test_r), 120.d0/2.d0**10
      read(*,*)
      
C c test the binomial probability function without variables
C       write(*,*) binomial_probability(0.5d0,10,3), 120.d0/2.d0**10
C       read(*,*)
C OK for this one regular integers such as "10" and "3" aren't the same as integer*8 numbers... blarg...
c-------------------------------------------------------------------------------------------------c



c-------------------------------------------------------------------------------------------------c
c test the binomial_fit subroutine with fabricated data from the binomial probability function      
c initailze the datas
      n = 20
      nm1 = 19
      p_test = 0.75d0

      do counter=0,19
        datas(counter+1) = binomial_probability(p_test,nm1,counter)
        write(*,*) counter, datas(counter+1)
      enddo

      call binomial_fit(datas,n,p_result)
      write(*,*) 'The result from binomial fit is ',p_result
      read(*,*)
c-------------------------------------------------------------------------------------------------c


c-------------------------------------------------------------------------------------------------c
c test the binomial_fit subroutine with actual data that I want to fit
      datas(1) = 0.7716049383d-3
      datas(2) = 0.3086419753d-2
      datas(3) = 0.7716049383d-2
      datas(4) = 0.1620370370d-1
      datas(5) = 0.2932098765d-1
      datas(6) = 0.4783950617d-1
      datas(7) = 0.7021604938d-1
      datas(8) = 0.9413580247d-1
      datas(9) = 0.1141975309d0
      datas(10) = 0.1288580247d0
      datas(11) = 0.1327160494d0
      datas(12) = 0.1234567901d0
      datas(13) = 0.1010802469d0
      datas(14) = 0.7253086420d-1
      datas(15) = 0.4166666667d-1
      datas(16) = 0.1620370370d-1
      
      n = 16
      nm1 = 15
      p_result = 0.5d0

C       do counter=0,15
C         write(*,*) counter+3, datas(counter+1), binomial_probability(p_result,nm1,counter)
C       enddo
C       read(*,*)
            
      call binomial_fit(datas,n,p_result)
      
      write(*,*) 'The result from binomial fit is ',p_result
      do counter=0,15
        write(*,*) counter+3, datas(counter+1), binomial_probability(p_result,nm1,counter)
      enddo
      
      end
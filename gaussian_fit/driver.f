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
      double precision datas(20,2), p_test, p_result, binomial_probability
      integer*8 n, nm1, counter, fileunit

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
        datas(counter+1,1) = counter + 1
        datas(counter+1,2) = binomial_probability(p_test,nm1,counter)
        write(*,*) datas(counter+1,1), datas(counter+1,2)
      enddo

      fileunit = 0
      call binomial_fit(datas,n,p_result,finalunit)
      write(*,*) 'The result from binomial fit is ',p_result
      read(*,*)
c-------------------------------------------------------------------------------------------------c


c-------------------------------------------------------------------------------------------------c
c test the binomial_fit subroutine with actual data that I want to fit
      do counter=1,16
        datas(x1,1) = counter
      enddo
      datas(1,2) = 0.7716049383d-3
      datas(2,2) = 0.3086419753d-2
      datas(3,2) = 0.7716049383d-2
      datas(4,2) = 0.1620370370d-1
      datas(5,2) = 0.2932098765d-1
      datas(6,2) = 0.4783950617d-1
      datas(7,2) = 0.7021604938d-1
      datas(8,2) = 0.9413580247d-1
      datas(9,2) = 0.1141975309d0
      datas(10,2) = 0.1288580247d0
      datas(11,2) = 0.1327160494d0
      datas(12,2) = 0.1234567901d0
      datas(13,2) = 0.1010802469d0
      datas(14,2) = 0.7253086420d-1
      datas(15,2) = 0.4166666667d-1
      datas(16,2) = 0.1620370370d-1
      
      n = 16
      nm1 = 15
      p_result = 0.5d0

C       do counter=0,15
C         write(*,*) counter+3, datas(counter+1), binomial_probability(p_result,nm1,counter)
C       enddo
C       read(*,*)
      
      fileunit = 42
      open(unit=fileunit,file='binomial_output.txt')
      call binomial_fit(datas,n,p_result,fileunit)
      close(fileunit)
      
      open(unit=fileunit,file='output.txt')
      write(*,*) 'The result from binomial fit is ',p_result
      do counter=0,15
        write(*,*) counter+3, datas(counter+1,2), binomial_probability(p_result,nm1,counter)
        write(fileunit,fmt='(I4,2X,2(SE20.10,2X))') counter+3, datas(counter+1,2), 
     .        binomial_probability(p_result,nm1,counter)
      enddo
      close(fileunit)
      
      end
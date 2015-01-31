c-------------------------------------------------------------------------------------------------c
      program stats
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This program will analyze a set of data with noise in both the x and y direction.              c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      include 'noise.inc'

c parameters for finding the expected b and m values given input noise
      double precision m, b, sigx, sigy, x(100), y(100), x_new(100), y_new(100)
      double precision x_ave, y_ave, xsq_ave, ysq_ave, xy_ave, mx_ave, mxsq_ave, my_ave, mysq_ave
      double precision bx_ave, bxsq_ave, by_ave, bysq_ave, sig_mx, sig_bx, sig_my, sig_by
      double precision mx, bx, my, by
      double precision mx_old, bx_old, my_old, by_old, sigmx_old, sigbx_old, sigmy_old, sigby_old
      integer ndata, ntrials, m_hist, b_hist
      integer mx_hist, bx_hist, my_hist, by_hist
      logical found_answer

c initialize the idum variable
      idum = -42

C c parameters for checking the gaussian number generator
C       double precision x_val, y_val, gauss
C       integer gausshist
C 
C c test the random gaussian number generator to see if it works (mean = 5, sd = 1)
C       open(unit=42,file='output/gaussian.txt')
C       do x1=0,1000
C         x_val = dble(x1)/100.d0
C         y_val = dexp((-1.d0*(x_val-5.d0)**2)/(2.d0))/dsqrt(2.d0*pi)
C         write(42,*) x_val,y_val
C       enddo
C       close(42)
C       
C c open up a histogram
C       call histinit
C       call hbook(gausshist,'gauss_hist',100,0.d0,10.d0,0.d0)
C       do x1=1,1000000
C         gauss = randgauss(idum,5.d0,1.d0)
C         call hfill(gausshist,gauss,1.d0)
C       enddo
C       call hist2file(gausshist,1.d-5)


c initialize some histograms so we can look at the final distributions
      call histinit
      call hbook(mx_hist,'mx',200,-10.d0,10.d0,0.d0)
      call hbook(bx_hist,'bx',200,-10.d0,10.d0,0.d0)
      call hbook(my_hist,'my',200,-10.d0,10.d0,0.d0)
      call hbook(by_hist,'by',200,-10.d0,10.d0,0.d0)
      
c we start off with a straight line described by y = mx + b
      m = 2.d0
      b = -3.d0
      ndata = 100
      
      do x1=1,ndata
        x(x1) = dble(x1)
        y(x1) = m*x(x1) + b
      enddo
      
c here are the user defined erorrs in the x and y direction
      sigx = 2.0d1
      sigy = 1.0d0

c loop over sigma x
      do x1=0,100
        sigy = 10.d0**(3.d0*dble(x1)/100.d0)

c initialize the averages of the slopes and y-intercepts
        mx_ave = 0.d0
        mxsq_ave = 0.d0
        bx_ave = 0.d0
        bxsq_ave = 0.d0
        my_ave = 0.d0
        mysq_ave = 0.d0
        by_ave = 0.d0
        bysq_ave = 0.d0
        mx_old = 0.d0
        bx_old = 0.d0
        my_old = 0.d0
        by_old = 0.d0
        sigmx_old = 0.d0
        sigbx_old = 0.d0
        sigmy_old = 0.d0
        sigby_old = 0.d0
        found_answer = .false.
        
        ntrials = 0
        do while(.not.found_answer)
          
          ntrials = ntrials + 1
          do x2=1,ndata
            x_new(x2) = randgauss(idum,x(x2),sigx)
            y_new(x2) = randgauss(idum,y(x2),sigy)
          enddo

c initialize all the averages        
          x_ave = 0
          y_ave = 0
          xsq_ave = 0
          ysq_ave = 0
          xy_ave = 0

c loop through the data and tally up the needed averages
          do x2=1,ndata
            x_ave = x_ave + x_new(x2)
            y_ave = y_ave + y_new(x2)
            xsq_ave = xsq_ave + x_new(x2)**2
            ysq_ave = ysq_ave + y_new(x2)**2
            xy_ave = xy_ave + x_new(x2)*y_new(x2)
          enddo

c calculate the data averages
          x_ave = x_ave/dble(ndata)
          y_ave = y_ave/dble(ndata)
          xsq_ave = xsq_ave/dble(ndata)
          ysq_ave = ysq_ave/dble(ndata)
          xy_ave = xy_ave/dble(ndata)
        
c use these to calculate the slopes and their respective averages
          mx = (xy_ave - x_ave*y_ave)/(xsq_ave - x_ave**2)
          bx = (xsq_ave*y_ave - xy_ave*x_ave)/(xsq_ave - x_ave**2)
          my = (xy_ave - x_ave*y_ave)/(ysq_ave - y_ave**2)
          by = (ysq_ave*x_ave - xy_ave*y_ave)/(ysq_ave - y_ave**2)

c add these slopes to the averages 
          mx_ave = mx_ave + mx
          mxsq_ave = mxsq_ave + mx**2
          bx_ave = bx_ave + bx
          bxsq_ave = bxsq_ave + bx**2
          my_ave = my_ave + my
          mysq_ave = mysq_ave + my**2
          by_ave = by_ave + by
          bysq_ave = bysq_ave + by**2

c check every 1000 trials 
          if (1000*(ntrials/1000).eq.ntrials) then
            sig_mx = dsqrt(ntrials*mxsq_ave - mx_ave**2)/dble(ntrials)
            sig_bx = dsqrt(ntrials*bxsq_ave - bx_ave**2)/dble(ntrials)
            sig_my = dsqrt(ntrials*mysq_ave - my_ave**2)/dble(ntrials)
            sig_by = dsqrt(ntrials*bysq_ave - by_ave**2)/dble(ntrials)
            
            if ((((mx_ave/dble(ntrials)-mx_old)/mx_old).lt.1.d-5).and.
     .        ((dabs(bx_ave/dble(ntrials)-bx_old)/bx_old).lt.1.d-5).and.
     .        ((dabs(my_ave/dble(ntrials)-my_old)/my_old).lt.1.d-5).and.
     .        ((dabs(by_ave/dble(ntrials)-by_old)/by_old).lt.1.d-5).and.
     .        ((dabs(sig_mx-sigmx_old)/sigmx_old).lt.1.d-5).and.
     .        ((dabs(sig_bx-sigbx_old)/sigbx_old).lt.1.d-5).and.
     .        ((dabs(sig_my-sigmy_old)/sigmy_old).lt.1.d-5).and.
     .        ((dabs(sig_by-sigby_old)/sigby_old).lt.1.d-5)) found_answer = .true.
     
            mx_old = mx_ave/dble(ntrials)
            bx_old = bx_ave/dble(ntrials)
            my_old = my_ave/dble(ntrials)
            by_old = by_ave/dble(ntrials)
            sigmx_old = sig_mx
            sigbx_old = sig_bx
            sigmy_old = sig_my
            sigby_old = sig_by
          endif
          
C c add all the results to the respective histograms
C           call hfill(mx_hist,mx,1.d0)
C           call hfill(bx_hist,bx,1.d0)
C           call hfill(my_hist,my,1.d0)
C           call hfill(by_hist,by,1.d0)
        enddo

c find the averages and standard deviations of the slopes and intercepts      
        mx_ave = mx_ave/dble(ntrials)
        mxsq_ave = mxsq_ave/dble(ntrials)
        bx_ave = bx_ave/dble(ntrials)
        bxsq_ave = bxsq_ave/dble(ntrials)
        my_ave = my_ave/dble(ntrials)
        mysq_ave = mysq_ave/dble(ntrials)
        by_ave = by_ave/dble(ntrials)
        bysq_ave = bysq_ave/dble(ntrials)

c find the standard deviations of the slopes and intercepts      
        sig_mx = dsqrt(mxsq_ave - mx_ave**2)
        sig_bx = dsqrt(bxsq_ave - bx_ave**2)
        sig_my = dsqrt(mysq_ave - my_ave**2)
        sig_by = dsqrt(bysq_ave - by_ave**2)

        write(*,fmt='(11(SE16.8,1X))') sigx, sigy, dble(ntrials), mx_ave, sig_mx, bx_ave, sig_bx, my_ave, sig_my,
     .        by_ave, sig_by
      enddo

C c output the histograms
C       call hist2file(mx_hist,1.d0)
C       call hist2file(bx_hist,1.d0)
C       call hist2file(my_hist,1.d0)
C       call hist2file(by_hist,1.d0)
C 
C       write(*,*) 'noise'
C       write(*,*) sigx, sigy
C       write(*,*) 'y versus x'
C       write(*,*) m, b
C       write(*,*) 'calculated'
C       write(*,*) mx, sig_mx
C       write(*,*) bx, sig_bx
C       write(*,*) 'x versus y'
C       write(*,*) 1.d0/m, -b/m
C       write(*,*) 'calculated'
C       write(*,*) my, sig_my
C       write(*,*) by, sig_by
      
      end  
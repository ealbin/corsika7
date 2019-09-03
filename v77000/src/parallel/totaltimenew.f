c=======================================================================
c
c  t o t a l t i m e n e w . f
c  ---------------------------
c     calculate sum of used cpu times of a set of corsika simulations
c     reading text file `corsika_timetable` of the current parallel
c     corsika subdirectory and create time information file `time.txt`.
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check totaltimenew.f -o totaltimenew
c     ifort -C -O0 -check bounds totaltimenew.f -o totaltimenew
c execution of ./totaltimenew by new postprocessing script:
c     ./postprocessnew.sh
c-----------------------------------------------------------------------
c      START TIME          STOP TIME       TIME (min)
c 1479599761.872620   1479602956.296335     53.240395
c LONGEST JOB: MPIID =    64 and Time =   1279.667251
c  Total number of jobs =    65
c Maximum size of group =    15
c TOTAL CPU TIME (days) =     0.434466
c time.txt001234        # written from script postprocessnew.sh
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program totaltimenew

      implicit double precision (a-h,o-z), integer (i-n)
      character czeile*60, comment*160
      logical lexist

c - - - - initializations:
      ia = 0
      ib = 0
      ic = 0
      im = 0
      is = 0
      it = 0
      timeanf = 0.d0
      timeend = 0.d0
      timemax = 0.d0
      timesum = 0.d0

c - - - - copied time.txt file as corsika_showertime available:
      inquire(file='corsika_showertime',exist=lexist)
      if ( lexist ) then
         open(unit=2,file='corsika_showertime',status='old',
     +        form='formatted',access='sequential')
         read(2,'(a)') czeile
         read(2,*) timeanf, timeend, timedif
         read(2,*,end=3,err=3) ic
         ic = ic + 1
         goto 4
    3    continue
         ic = int(timedif*7.6543)
    4    continue
         close(unit=2)
         ib = int(5.d0*ic/7.)
         im = ib
         timemax = timedif * 60.d0
         timesum = (timeend-timeanf) / 1.28d0
      endif

c - - - - check summary timetable file in this path:
      inquire(file='corsika_timetable',exist=lexist)
      if ( lexist ) then
         open(unit=1,file='corsika_timetable',status='old',
     +        form='formatted',access='sequential')
         read(1,*,end=2,err=2) ia,ib,timea,timeb
         timesum = timesum + (timeb - timea)
         do ic=2,1234567890
            read(1,*,end=2,err=2) ia,ib,timea,timeb
            if ( timeb-timea .gt. timemax ) timemax = timeb-timea 
            if ( ia .gt. im ) im = ia
            timesum = timesum + (timeb-timea)   
         enddo
    2    continue
         close(unit=1)
      endif

c - - - - write to new time information file (totaltimenew.out):
      write(*,'(5x,''START TIME'',10x,''STOP TIME'',7x,''TIME (min)'')')
      write(*,'(f17.6,f20.6,f14.6)') timeanf, timeend,
     +   (timeend-timeanf)/60.
      write(*,'(''LONGEST JOB: MPIID ='',i6,'' and Time ='',f14.6)')
     +   ib, timemax 
      write(*,'('' Total number of jobs ='',i6)') ic-1
      write(*,'(''Maximum size of group ='',i6)') im
      write(*,'(''TOTAL CPU TIME (days) ='',f13.6)') timesum/86400.d0

c - - - - end-of program (totaltimenew.out => time.txtnnnnnn).
      stop
      end

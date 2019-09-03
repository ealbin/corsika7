c=======================================================================
c
c  t o t a l t i m e . f
c  ---------------------
c     calculate sum of used cpu times of a set of corsika simulations
c     by checking all `.lst` files (start time and end of run time);
c     the line `head -107` in the script `totaltime.sh` must give the
c     actual line of all simulation protocol files `DAT*.lst`   
c     where the start time like `PRESENT TIME : 30.09.2012  04:40:41`
c     is printed;
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check totaltime.f -o totaltime
c     f77 -fbounds-check -m32 totaltime.f -o totaltime
c     ifort -C totaltime.f -o totaltime
c execution by script:
c     ./postprocess.sh
c-----------------------------------------------------------------------
cDAT401200.lst
c PRESENT TIME : 30.09.2012  04:40:41
c PRESENT TIME : 30.09.2012  07:22:52
c
c ========== END OF RUN ================================================
cDAT401201.lst
c PRESENT TIME : 30.09.2012  04:42:49
c PRESENT TIME : 30.09.2012  06:38:34
c
c ========== END OF RUN ================================================
c-----------------------------------------------------------------------
c     START TIME          STOP TIME         TIME (min)
c 1352728629.908327   1352729950.580470   22.011202
c LONGEST JOB: MPIID = 16 and Time = 671.813086
c  Total number of jobs = 55 
c Maximum size of group = 16
c TOTAL CPU TIME (days) = 0.186547
c-----------------------------------------------------------------------

      program totaltime

      implicit double precision (a-h,o-z), integer (i-n)

      dimension idate(6), ldate(6), idats(6)

      character clines(5)*40, czeile0*1, czeile1*50,
     +    czeile3*34, czeile4*23, czeile5*23, czeile6*23

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - initializations:
      data czeile1/
     +   '    START TIME          STOP TIME       TIME (min)'/
      data czeile3/'LONGEST JOB: MPIID =    and Time ='/
      data czeile4/' Total number of jobs ='/
      data czeile5/'Maximum size of group ='/
      data czeile6/'TOTAL CPU TIME (days) ='/
      chdsek = 1.d0 / 3600.d0
      hoursum = 0.d0
      ignored = 0
      jobs = 13
      idats(1) = 1970
      idats(2) = 1
      idats(3) = 1
      idats(4) = 0
      idats(5) = 0
      idats(6) = 0
      mjob = 13

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - loop over all DAT....lst files:
      do  lfiles=1,12345678

c - - - - - - read 5 lines (2,3 containing times):
         read(*,'(a)',end=4,err=4) clines(1)
         read(*,'(a)') clines(2)
         read(*,'(a)') clines(3)
         read(*,'(a)') clines(4)
         read(*,'(a)') clines(5)
         if ( lfiles .eq. 1 ) then
            write(*,'(19x,a)') clines(1)
            read(clines(1)(4:9),'(i6)') lrunnr
         endif 

c - - - - - - copy start date and time to integer vector:
         if ( index(clines(2),'PRESENT TIME') .le. 0 ) then
            write(*,*) lfiles,'th file ignored: ',clines(1)
            ignored = ignored + 1
            if ( ignored .lt. 6 ) goto 3
            write(*,'('' Check command '', 
     +         '' `head -... $i | tail -1 >> $timeinpt` '',
     +         '' in the script!'')')
            goto 4
         endif
         read(clines(2)(17:18),'(i2)') idate(3)
         read(clines(2)(20:21),'(i2)') idate(2)
         read(clines(2)(23:26),'(i4)') idate(1)
         read(clines(2)(29:30),'(i2)') idate(4)
         read(clines(2)(32:33),'(i2)') idate(5)
         read(clines(2)(35:36),'(i2)') idate(6)

c - - - - - - copy end date and time to integer vector:
         if ( index(clines(3),'PRESENT TIME') .le. 0 ) then
            write(*,*) lfiles,'th file ignored: ',clines(1)
            goto 3
         endif
         read(clines(3)(17:18),'(i2)') ldate(3)
         read(clines(3)(20:21),'(i2)') ldate(2)
         read(clines(3)(23:26),'(i4)') ldate(1)
         read(clines(3)(29:30),'(i2)') ldate(4) 
         read(clines(3)(32:33),'(i2)') ldate(5)
         read(clines(3)(35:36),'(i2)') ldate(6)

c - - - - - - sum up time difference in hours:
         call datsek(idate,idats,idsek)
         call datsek(ldate,idats,ldsek)
         hoursd = chdsek * (ldsek - idsek)
         hoursum = hoursum + hoursd

c - - - - - - being ignored a set of 5 lines:
    3    continue
**       if ( mod(lfiles,50) .eq. 0 )
**          write(*,'(i18,'' DAT'',i6.6,''*.lst done,'')') lfiles,lrunnr
      enddo
    4 continue
**    write(*,'(i18,'' DAT'',i6.6,''*.lst total.'')') lfiles-1,lrunnr

c - - - - print sum of all simulation times:
      if ( lfiles .gt. 3 ) write(*,'(19x,6(''. ''))')
      if ( lfiles .gt. 1 ) write(*,'(i18,1x,a)') lfiles-1,clines(1)
      write(*,'(f36.6,f16.6,''   hours   days'')') hoursum, hoursum/24.

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - read existing lines from file time.txt:
      open(unit=2,file='time.txt',form='formatted',access='sequential',
     +     status='old')
      read(2,'(a)') czeile0
      read(2,*,end=7,err=7) startime, stoptime, umaxtime
      read(2,*,end=7,err=7) jobs
    7 continue
      close(unit=2)

c - - - - add statistic lines to file time.txt incl. total time:
      open(unit=3,file='time.txt',form='formatted',access='sequential',
     +     status='unknown')
      write(3,'(a)') czeile1
      write(3,'(f17.6,f20.6,f13.6)') startime,stoptime,umaxtime
      write(3,'(a,i3,a,f13.6)') czeile3(1:20), mod(mjob,100),
     +    czeile3(24:34), umaxtime*24.
      write(3,'(a,i5)') czeile4, jobs
      write(3,'(a,i5)') czeile5, mjob
      write(3,'(a,f12.6)') czeile6, hoursum/24.
      close(unit=3) 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - end of program.
    9 continue
      stop
      end
c=======================================================================
 
      subroutine datsek(idate,idat0,idsek)
 
c-----------------------------------------------------------------------
c
c  calculate total number of elapsed seconds on a date given by an
c  integer vector idate(6) with six elements (year, month, day, hours,
c  minutes, seconds) corresponding to a start date idat0(6).
c
c-----------------------------------------------------------------------
 
      dimension idate(6),idat0(6)
      integer imons(13,4)
      data imons/0,31,59,90,120,151,181,212,243,273,304,334,365,
     +         0,31,60,91,121,152,182,213,244,274,305,335,366,
     +         0,31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     +         0,31,29,31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
 
c--------- count leap days of bygone years -----------------------------
      idsek = 0
      if (idate(1).lt.100) idate(1)=1900+idate(1)
      if (idat0(1).lt.100) idat0(1)=1900+idat0(1)
      il0 = 1
      if (mod(idat0(1),4).eq.0.and.(mod(idat0(1),100).ne.0.or
     +   .mod(idat0(1),400).eq.0)) il0 = 2
      ileap = 0
      ila = 1
      do  1  i=idat0(1),idate(1)
         if (mod(i,4).eq.0.and.(mod(i,100).ne.0.or.mod(i,400).eq.0))then
            ileap = ileap + 1
            if (i.eq.idate(1)) then
               ila = 2
               ileap = ileap - 1
            endif
         endif
    1 continue
 
c--------- calculate total number of seconds ---------------------------
      idsek = (365 * (idate(1)-idat0(1)) + ileap + imons(idate(2),ila) -
     - imons(idat0(2),il0) + idate(3) - idat0(3) ) * 86400 +
     + ( (idate(4)-idat0(4)) * 60 + idate(5) - idat0(5) ) * 60 +
     + idate(6) - idat0(6)
 
      return
      end

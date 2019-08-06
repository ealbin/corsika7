c=======================================================================
c
c  r e a d t i m e s n e w . f
c  ---------------------------
c     calculate sum of used cpu times of a set of corsika simulations
c     by checking all `.lst` files (start time and end of run time). 
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check readtimesnew.f -o readtimesnew
c     ifort -C -check bounds readtimesnew.f -o readtimesnew
c execution by script:
c     ./readtimesnew.sh
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

      program readtimesnew

      implicit double precision (a-h,o-z), integer (i-n)

      character czeile(5)*40, listfile*36

      dimension idate(6), ldate(6), idats(6)

      chdsek = 1.d0 / 3600.d0
      hoursum = 0.d0

      idats(1) = 1970
      idats(2) = 1
      idats(3) = 1
      idats(4) = 0
      idats(5) = 0
      idats(6) = 0

      do  lfiles=1,12345678

c - - - - - - read 5 lines (2,3 containing times):
         read(*,'(a)',end=8,err=8) czeile(1)
         read(*,'(a)') czeile(2)
         read(*,'(a)') czeile(3)
         read(*,'(a)') czeile(4)
         read(*,'(a)') czeile(5)
         if ( lfiles .eq. 1 ) write(*,'(19x,a)') czeile(1)
         read(czeile(1)(4:9),'(i6)') jrunnr

c - - - - - - copy start date and time to integer vector:
         if ( index(czeile(2),'PRESENT TIME') .le. 0 ) then
            write(*,*) lfiles,'th file ignored: ',czeile(1)
            goto 7
         endif
         read(czeile(2)(17:18),'(i2)') idate(3)
         read(czeile(2)(20:21),'(i2)') idate(2)
         read(czeile(2)(23:26),'(i4)') idate(1)
         read(czeile(2)(29:30),'(i2)') idate(4)
         read(czeile(2)(32:33),'(i2)') idate(5)
         read(czeile(2)(35:36),'(i2)') idate(6)

c - - - - - - copy end date and time to integer vector:
         if ( index(czeile(3),'PRESENT TIME') .le. 0 ) then
            write(*,*) lfiles,'th file ignored: ',czeile(1)
            goto 7
         endif
         read(czeile(3)(17:18),'(i2)') ldate(3)
         read(czeile(3)(20:21),'(i2)') ldate(2)
         read(czeile(3)(23:26),'(i4)') ldate(1)
         read(czeile(3)(29:30),'(i2)') ldate(4) 
         read(czeile(3)(32:33),'(i2)') ldate(5)
         read(czeile(3)(35:36),'(i2)') ldate(6)

c - - - - - - sum up time difference in hours:
         call datsek(idate,idats,idsek)
         call datsek(ldate,idats,ldsek)
         hoursd = chdsek * (ldsek - idsek)
         hoursum = hoursum + hoursd

c - - - - - - being ignored a set of 5 lines:
    7    continue

      enddo

    8 continue

c - - - - print sum of run times:
      if ( lfiles .gt. 3 ) write(*,'(19x,6(''. ''))')
      if ( lfiles .gt. 1 ) write(*,'(19x,a)') czeile(1)
      write(*,'(f42.4,'' hours ='',f10.4,'' days'')')
     +   hoursum,hoursum/24.
      open(unit=9,file='readtimesnew.out',access='sequential',
     +     form='formatted',status='unknown')
      write(9,*) hoursum/24., jrunnr
      close(unit=9)

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

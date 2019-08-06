c=======================================================================
c
c  s u m l o n g i f i l e s . f
c  -----------------------------
c     sum up content of a number of long files for one shower;
c     print sum of energy quantities,
c-----------------------------------------------------------------------
c  gfortran -O0 -fbounds-check sumlongifiles.f -o sumlongifiles
c  ifort -C -O0 sumlongifiles.f -o sumlongifiles
c  ls -1 DAT*.long > sumlongifiles.002543
c  sumlongifiles < sumlongifiles.i002543 > sumlongifiles.out002543
c  /bin/cp -p sumlongifiles.sum002543 DAT002543-999989999.long
c-----------------------------------------------------------------------
c LONGITUDINAL DISTRIBUTION IN   jvs VERTICAL STEPS OF    fgr. G/CM**2
c DEPTH    GAMMAS   POSITRONS   ELECTRONS         MU+         MU-     HADRONS  
c  5. 1.15077E+04 2.80446E+01 2.52274E+03 0.00000E+00 0.00000E+00 1.00000E+00
c 10. 1.19640E+04 1.03913E+02 1.57065E+02 0.00000E+00 0.00000E+00 1.00000E+00
c 15. 1.89921E+04 1.02062E+02 6.02637E+02 0.00000E+00 0.00000E+00 1.00000E+00
c 20. 2.21814E+04 6.25936E+02 7.22297E+02 0.00000E+00 0.00000E+00 1.00000E+00
c 25. 2.33861E+04 2.67479E+02 6.90617E+02 0.00000E+00 0.00000E+00 1.00000E+00
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program sumlongifiles

      implicit double precision (a-h,o-z), integer (i-n)

      character clongname(50000)*120, czeilong(4)*120, czeile*120
      character chlongsum*23, chpartic*29

      dimension qdistrb(0:9,1040), qenergy(0:9,1040), qzeile(0:9)
      dimension qxdistr(1040), qxenerg(1040)
      dimension longlen(50000)
      real pdata(400)
      
c - - - - - - read all long files and keep names in a character array - -
      do  long=1,12543
         read(*,'(a)',err=101,end=102) clongname(long)
         jl = 120 + 1 
  100    continue 
         jl = jl - 1
         if ( clongname(long)(jl:jl) .eq. ' ' ) goto 100
         longlen(long) = jl       
      enddo
      goto 102
  101 continue
      write(*,*) ' e r r o r   reading clongname ',clongname(long)
  102 continue  
      long = long - 1
      read(clongname(1)(4:9),'(i6)') lrunnr      
      chlongsum = 'sumlongifiles.sum000000'
      write(chlongsum(18:23),'(i6.6)') lrunnr

c - - - - - - read first particle data file (single precision) - - - - -
      chpartic = clongname(1)(1:longlen(1)-5) 
      open(unit=2,file=chpartic,form='unformatted',status='old')
      read(2) (pdata(il),il=1,400)
      close(unit=2)
      lsubblo = 312
      if ( 217433.0 .lt. pdata(273+1) .and. pdata(273+1) .lt. 217433.2 )
     +   lsubblo = 273
      parengy = pdata(lsubblo+4)

c - - - - - - read first long file - - - - - - - - - - - - - - - - - - -
      open(unit=1,file=clongname(1)(1:longlen(1)),
     +            form='formatted',status='old')
      ! - - - - - clear all arrays:
      do  is=1,1040
      do  il=0,9
         qdistrb(il,is) = 0.d0
         qenergy(il,is) = 0.d0
      enddo
      qxdistr(is) = 0.d0
      qxenerg(is) = 0.d0
      enddo
      ! - - - - - read first longi table:
      read(1,'(a)',end=119,err=118) czeilong(1)
      read(1,'(a)',end=119,err=118) czeilong(2)
      read(czeilong(1)(31:35),'(i5)') lsteps
      read(czeilong(1)(55:59),'(f5.0)') grstep
      write(czeilong(1)(80:86),'(i7.6)') lrunnr 
      do  is=1,lsteps
         read(1,*) (qzeile(il),il=0,9)
         do  il=1,9
            qdistrb(il,is) = qzeile(il)
         enddo 
         qxdistr(is) = qzeile(0)
      enddo
      ! - - - - - read second longi table:
      read(1,'(a)',end=119,err=118) czeilong(3)
      read(1,'(a)',end=119,err=118) czeilong(4)
      ! read(czeilong(3)(33:37),'(i5)') lsteps
      ! read(czeilong(3)(57:61),'(f5.0)') grstep
      write(czeilong(3)(82:88),'(i7.6)') lrunnr 
      do  is=1,lsteps
         read(1,*) (qzeile(i),i=0,9) 
         do  il=1,9
            qenergy(il,is) = qzeile(il)
         enddo 
         qxenerg(is) = qzeile(0)
      enddo
      close(unit=1)
      write(*,'(3x,''sumlongifiles.out'',i6.6,''   lsteps ='',i4,3x,
     +   ''grstep ='',i4)') lrunnr,lsteps,int(grstep)

c - - - - - - read following long files and sum content - - - - - - - -
      do  ifile=2,long
         if ( mod(ifile,100) .eq. 0 ) write(*,*) '       ifile =',ifile
         open(unit=1,file=clongname(ifile)(1:longlen(ifile)),
     +               form='formatted',status='old')
         read(1,'(a)',end=119,err=118) czeile
         read(1,'(a)',end=119,err=118) czeile
         do  is=1,lsteps
            read(1,*,end=129,err=128) (qzeile(il),il=0,9)
            do  il=1,9
               qdistrb(il,is) = qdistrb(il,is) + qzeile(il)
            enddo
         enddo
         read(1,'(a)',end=119,err=118) czeile
         read(1,'(a)',end=119,err=118) czeile
         do  is=1,lsteps
            read(1,*,end=129,err=128) (qzeile(i),i=0,9)
            do  il=1,9
               qenergy(il,is) = qenergy(il,is) + qzeile(il)
            enddo
         enddo
  128    continue
  129    continue         
         close(unit=1)
      enddo
      if ( mod(ifile-1,100) .ne. 0 ) write(*,*) '       ifile =',ifile-1

c - - - - - - write sum tables - - - - - - - - - - - - - - - - - - - - -
      open(unit=8,file=chlongsum,form='formatted',status='unknown')
      write(8,'(a)') czeilong(1)
      write(8,'(a)') czeilong(2)
      do  is=1,lsteps
         write(8,'(f6.0,1p,9e12.5)') qxdistr(is),(qdistrb(il,is),il=1,9)
      enddo
      write(8,'(a)') czeilong(3)
      write(8,'(a)') czeilong(4)
      engysum = 0.d0
      do  is=1,lsteps
         engysum = engysum + qenergy(9,is) 
         write(8,'(f6.1,1p,9e12.5)') qxenerg(is),(qenergy(il,is),il=1,9)
      enddo
      close(unit=8)
c - - - - - - calculate energy sum of all levels - - - - - - - - - - - -
      write(*,*) '     parengy =',parengy,' GeV'
      write(*,*) '     engysum =',engysum
  118 continue
  119 continue
      stop
      end

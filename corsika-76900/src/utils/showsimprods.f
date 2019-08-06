c=======================================================================
c
c  s h o w s i m p r o d s . f
c  ---------------------------
c     reading corsika particle data files from a (one name per line)
c     tabular each simulation ([DAT,CER]iiiiii) will be printed out
c     in one line the following quantities:
c          Primary, lg(E), theta, phi, nsh,  runnr, size,
c              obslvme, h1stme, thilev, thiwmax, lg(thirad),
c              verspgm, models, rundate, Xmagn, Zmagn;
c     on 32bit machines all 64bit simulations will be detected
c     correctly and in the list marked by `_64`, Cherenkov files by
c     `_CE`, and both by `_C4`; particle data files must be available
c     as [DAT,CER]iiiiii, protocol files as [DAT,CER]iiiiii.lst;
c     otherwise additional conditions must be implemented.
c-----------------------------------------------------------------------
c compilation:
c           gfortran -fbounds-check showsimprods.f -o showsimprods
c           f77 -fbounds-check showsimprods.f -o showsimprods
c           ifort -C -check bounds showsimprods.f -o showsimprods 
c execution:
c           ls -l DAT?????? | ./showsimprods
c           ls -l [c,f,p]?e??m???_* | ./showsimprods
c           ls -l [C,D]???????? | ./showsimprods > showsimprods.corsika
c-----------------------------------------------------------------------
c     output-files:
c           unit=*: tabular output of simulations.
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                             author: juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
c         Primary   lg(E)  theta    phi  nsh  runnr    size  obslvme  h1stme   .....
c Iron       5626   16.00    0.0    0.0   1  199079     3.9  1413.82   18350.  .....
c _stck_in_     4   15.09    0.0    0.0   1  169051     2.4   194.00   18765.  .....
c Fluorine   1909   16.00    0.0    0.0   1  199070     3.9  1416.51   22222.  .....
c proton       14   16.00   30.0   -3.3   1  199080     3.4  1429.25   22224.  .....
c Manganese  5525   16.00   30.0   -3.3   1  199082    10.1  1428.13   22222.  .....
c proton _64   14   17.50   37.9 -138.0   1  000044    90.8  -500.00  -24880.  .....
c=======================================================================

      program showsimprods

      implicit double precision (a-h,o-z), integer(i-n) 
      parameter (nfmx=50000,nchx=160)

      character cdata(nfmx)*160,cdat*160,csklst*160,czeile*80,crunh*4
      character qpatext(0:200)*19,chemical(0:200)*12,chfmti*5,chaspec*1

      dimension fsize(nfmx),qrunh(312),qevth(312)
      dimension nflen(nfmx),nfsho(nfmx),ntskf(nfmx),nbits(nfmx)
      real pdata(5733)

      logical lstexi

      equivalence (crunh,prunh)
      data crunh/'RUNH'/, chaspec/' '/

      data chemical/  ' unknown    ',
     c ' Hydrogen   ',' Helium     ',' Lithium    ',' Beryllium  ',
     c ' Boron      ',' Carbon     ',' Nitrogen   ',' Oxygen     ',
     c ' Fluorine   ',' Neon       ',' Sodium     ',' Magnesium  ',
     c ' Aluminium  ',' Silicon    ',' Phosphorus ',' Sulfur     ',
     c ' Chlorine   ',' Argon      ',' Potassium  ',' Calcium    ',
     c ' Scandium   ',' Titanium   ',' Vanadium   ',' Chromium   ',
     c ' Manganese  ',' Iron       ',' Cobalt     ',' Nickel     ',
     c ' Copper     ',' Zinc       ',' Gallium    ',' Germanium  ',
     c ' Arsenic    ',' Selenium   ',' Bromine    ',' Krypton    ',
     c ' Rubidium   ',' Strontium  ',' Yttrium    ',' Zirconium  ',
     c ' Niobium    ',' Molybdenum ',' Technetium ',' Ruthenium  ',
     c ' Rhodium    ',' Palladium  ',' Silver     ',' Cadmium    ',
     c ' Indium     ',' Tin        ',' Antimony   ',' Tellurium  ',
     c ' Iodine     ',' Xenon      ',' Caesium    ',' Barium     ',
     c ' Lanthanum  ',' Cerium     ',' Praseodym. ',' Neodymium  ',
     c ' Promethium ',' Samarium   ',' Europium   ',' Gadolinium ',
     c ' Terbium    ',' Dysprosium ',' Holmium    ',' Erbium     ',
     c ' Thulium    ',' Ytterbium  ',' Lutetium   ',' Hafnium    ',
     c ' Tantalum   ',' Tungsten   ',' Rhenium    ',' Osmium     ',
     c ' Iridium    ',' Platinum   ',' Gold       ',' Mercury    ',
     c ' Thallium   ',' Lead       ',' Bismuth    ',' Polonium   ',
     c ' Astatine   ',' Radon      ',' Francium   ',' Radium     ',
     c ' Actinium   ',' Thorium    ',' Protactin. ',' Uranium    ',
     c ' Neptunium  ',' Plutonium  ',' Americium  ',' Curium     ',
     c ' Berkelium  ',' Californium',' Einsteinium',101*'            '/

      data (qpatext(i),i=0,60)/ '               ',
     c' gamma             ',' positron          ',' electron          ',
     c' _stck_in_         ',' muon+             ',' muon-             ',
     c' pi0               ',' pi+               ',' pi-               ',
     c' K0long            ',' K+                ',' K-                ',
     c' neutron           ',' proton            ',' anti proton       ',
     c' K0short           ','                   ',' Lambda            ',
     c' Sigma+            ',' Sigma0            ',' Sigma-            ',
     c' Xi0               ',' Xi-               ',' Omega-            ',
     c' anti neutron      ',' anti Lambda       ',' anti Sigma-       ',
     c' anti Sigma0       ',' anti Sigma+       ',' anti Xi0          ',
     c' anti Xi+          ',' anti Omega+       ',16*'                ',
     c'                   ',' omega             ',' rho0              ',
     c' rho+              ',' rho-              ',' Delta++           ',
     c' Delta+            ',' Delta0            ',' Delta-            ',
     c' anti Delta--      ',' anti Delta-       ',' anti Delta0       '/

      data (qpatext(i),i=61,129)/
     c' anti Delta+       ',' K*0               ',' K*+               ',
     c' K*-               ',' anti K*0          ',' electron neutrino ',
     c' anti elec neutrino',' muon neutrino     ',' anti muon neutrino',
     c'                   ',' eta=>2*gamma      ',' eta=>3*pi0        ',
     c' eta=>pi+pi-pi0    ',' eta=>pi+pi-gamma  ',40*'                ',
     c'                   ',' D0                ',' D+                ',
     c' anti D-           ',' anti D0           ',' Ds+               ',
     c' anti Ds-          ',' eta c             ',' D*0               ',
     c' D*+               ',' anti D*-          ',' anti D*0          ',
     c' D*s+              ',' anti D*s-         ','                   '/

      data (qpatext(i),i=130,200)/
     c' J/psi             ',' tau +             ',' tau -             ',
     c' tau neutrino      ',' anti tau neutrino ','                   ',
     c'                   ',' Lambda c +        ',' Xi c +            ',
     c' Xi c 0            ',' Sigma c ++        ',' Sigma c +         ',
     c' Sigma c 0         ',' Xi c prime +      ',' Xi c prime 0      ',
     c' Omega c 0         ','                   ','                   ',
     c'                   ',' anti Lambda c -   ',' anti Xi c -       ',
     c' anti Xi c 0       ',' anti Sigma c --   ',' anti Sigma c -    ',
     c' anti Sigma c 0    ',' anti Xi c prime - ',' anti Xi c prime 0 ',
     c' anti Omega c 0    ','                   ','                   ',
     c'                   ',' Sigma c * ++      ',' Sigma c * +       ',
     c' Sigma c * 0       ',7*'                 ',' anti Sigma c * -- ',
     c' anti Sigma c * -  ',' anti Sigma c * 0  ',25*'                ',
     c' Cherenkov photon  ','                   '/

      qpatext(176) = ' Stnd              '
      qpatext(196) = ' Thin              '
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      obslev = 0.1
      lthi = 0

c--read first file name and fix position of file size-------------------
      ifil = 1
      read(*,'(a)',end=411,err=411) cdat
      if ( index(cdat,'-rw') .le. 0 ) then
         write(*,'(/,14x,''Use of `ls -1 ...` not supported, but'',
     +      '' `ls -l ...` is ok.'',/)')
         goto 445
      endif
      ibl1 = index(cdat,' ')
      ibl2 = index(cdat(ibl1+2:nchx),' ') ! after -rwx-`s 2 blanks.
      ibl2 = ibl2 + ibl1 + 1
      ibl3 = index(cdat(ibl2+1:nchx),' ')
      ibl3 = ibl3 + ibl2
      ibl4 = index(cdat(ibl3+1:nchx),' ')
      ibl4 = ibl4 + ibl3
      imon = index(cdat,'Jan') ! test string of date and time.
      if ( imon .eq. 0 ) imon = index(cdat,'Feb')
      if ( imon .eq. 0 ) imon = index(cdat,'Mar')
      if ( imon .eq. 0 ) imon = index(cdat,'Mär')
      if ( imon .eq. 0 ) imon = index(cdat,'Apr')
      if ( imon .eq. 0 ) imon = index(cdat,'Mai')
      if ( imon .eq. 0 ) imon = index(cdat,'May')
      if ( imon .eq. 0 ) imon = index(cdat,'Jun')
      if ( imon .eq. 0 ) imon = index(cdat,'Jul')
      if ( imon .eq. 0 ) imon = index(cdat,'Aug')
      if ( imon .eq. 0 ) imon = index(cdat,'Sep')
      if ( imon .eq. 0 ) imon = index(cdat,'Oct')
      if ( imon .eq. 0 ) imon = index(cdat,'Okt')
      if ( imon .eq. 0 ) imon = index(cdat,'Nov')
      if ( imon .eq. 0 ) imon = index(cdat,'Dec')
      if ( imon .eq. 0 ) imon = index(cdat,'Dez')
      if ( imon .eq. 0 ) then
         imon = ibl4 + index(cdat(ibl4+1:nchx),'-') - 4
      else
         if ( index(cdat,'. ') .eq. imon-2 ) imon = imon - 4
      endif
      iright= imon - 2
      ileft = ibl4 + 1
      write(chfmti(1:5),'(''(i'',i2,'')'')') iright-ileft+1
      read(cdat(ileft:iright-3),chfmti) isize
      fsize(ifil) = 1.d-3 * isize ! now in Mbytes.
      if ( fsize(ifil) .lt. 0.1 ) fsize(ifil) = 0.1
      itaskcer = index(cdat,'CER')
      itaskdat = index(cdat,'DAT')
      isimprod = 0
      if ( itaskdat .le. 0 .and. itaskcer .le. 0 ) then
         isimprod = index(cdat,'coe')
         if ( isimprod .le. 0 ) isimprod = index(cdat,'fee')
         if ( isimprod .le. 0 ) isimprod = index(cdat,'hee')
         if ( isimprod .le. 0 ) isimprod = index(cdat,'fee')
         if ( isimprod .le. 0 ) isimprod = index(cdat,'pre')
         if ( isimprod .le. 0 ) isimprod = index(cdat,'sie')
         itask = isimprod 
         ilen = nchx
  400    continue
         ilen = ilen - 1
         if ( cdat(ilen:ilen) .eq. ' ' ) goto 400
         cdata(ifil) = cdat(isimprod:ilen)
         nflen(ifil) = ilen - isimprod + 1
         ntskf(ifil) = 900000 + ifil
      else
         itask = -8
         if ( itaskcer .gt. 0 ) then
            itask = itaskcer
         elseif ( itaskdat .gt. 0 ) then
            itask = itaskdat
         endif 
         read(cdat(itask+3:itask+8),'(i6)') ntask
         ntskf(ifil) = ntask
      endif
c--check structure of file names (with / or without):
      iddot = index(cdat,':')
      isla = index(cdat,'/')
      islb = 0
      if ( isla .gt. 0 ) then
         if ( iddot .gt. 0 ) isla = iddot + 4
         islb = nchx
  401    continue
         islb = islb - 1
         if ( cdat(islb:islb) .ne. '/' ) goto 401
      elseif ( itask .gt. 0 ) then
         isla = itask
      endif
c--check length of first line in the file list:
      ilen = nchx
  402 continue
      ilen = ilen - 1
      if ( cdat(ilen:ilen) .eq. ' ' ) goto 402
      cdata(ifil) = cdat(itask:ilen)
      nflen(ifil) = ilen - itask + 1
      nfsho(ifil) = 1
c--check `[DAT,CER]iiiiii.lst` file to get number of simulated showers:
      if ( index(cdat,'.part') .eq. 0 ) then
         csklst = cdata(ifil)(1:nflen(ifil))//'.lst'
      elseif ( index(cdat,'.part') .gt. ilen-5 ) then
         csklst = cdata(ifil)(1:nflen(ifil)-5)//'.lst'
      endif
      inquire(file=csklst,exist=lstexi)
      if ( lstexi ) then
         open(unit=2,file=csklst,form='formatted',status='old')
         do  iq=1,123456
            read(2,'(a)',end=404,err=404) czeile
            if ( index(czeile,'NUMBER OF GENERATED').ge.1 ) then
               close(2)
               read(czeile(30:40),'(i11)') nfsho(ifil)
               if ( nfsho(ifil) .ge. 100000 ) nfsho(ifil) = 99999
               goto 404
            endif
         enddo
      endif
  404 continue
 
c--read run parameters including file names (more than 1 file)----------
      do  ifil=2,nfmx
         read(*,'(a)',end=411,err=411) cdat
         read(cdat(ileft:iright-3),chfmti) isize
         fsize(ifil) = 1.d-3 * isize ! now in Mbytes.
         if ( fsize(ifil) .lt. 0.1 ) fsize(ifil) = 0.1
         itaskcer = index(cdat,'CER')
         itaskdat = index(cdat,'DAT')
         isimprod = 0
         if ( itaskdat .le. 0 .and. itaskcer .le. 0 ) then
           isimprod = index(cdat,'coe')
           if ( isimprod .le. 0 ) isimprod = index(cdat,'fee')
           if ( isimprod .le. 0 ) isimprod = index(cdat,'hee')
           if ( isimprod .le. 0 ) isimprod = index(cdat,'fee')
           if ( isimprod .le. 0 ) isimprod = index(cdat,'pre')
           if ( isimprod .le. 0 ) isimprod = index(cdat,'sie')
           itask = isimprod
           ntskf(ifil) = 900000 + ifil
         else
           itask = -8
           if ( itaskcer .gt. 0 ) then
             itask = itaskcer
           elseif ( itaskdat .gt. 0 ) then
             itask = itaskdat
           endif
           read(cdat(itask+3:itask+8),'(i6)') ntask
           ntskf(ifil) = ntask
         endif
         if ( isimprod .gt. 0 ) itask = isimprod
c----check structure of file names (with / or without):
         iddot = index(cdat,':')
         isla = index(cdat,'/')
         islb = 0
         if ( isla .gt. 0 ) then
            if ( iddot .gt. 0 ) isla = iddot + 4
            islb = nchx
  406       continue
            islb = islb - 1
            if ( cdat(islb:islb) .ne. '/' ) goto 406
         elseif ( itask .gt. 0 ) then
            isla = itask
         endif
c----check length of first line in the file list:
         ilen = nchx
  407    continue
         ilen = ilen - 1
         if ( cdat(ilen:ilen) .eq. ' ' ) goto 407
         write(*,*) '   ifil = ',ifil,
     +      '   ilen = ',ilen,'  itask = ',itask
         cdata(ifil) = cdat(itask:ilen)
         nflen(ifil) = ilen - itask + 1
         nfsho(ifil) = 1
c----check protocol file `.lst` file to get number of simulated showers:
         if ( index(cdat,'.part') .eq. 0 ) then
            csklst = cdata(ifil)(1:nflen(ifil))//'.lst'
         elseif ( index(cdat,'.part') .gt. ilen-5 ) then
            csklst = cdata(ifil)(1:nflen(ifil)-5)//'.lst'
         endif
         inquire(file=csklst,exist=lstexi)
         if ( lstexi ) then
            open(unit=2,file=csklst,form='formatted',status='old')
            do  iq=1,123456
               read(2,'(a)',end=409,err=409) czeile
               if ( index(czeile,'NUMBER OF GENERATED').ge.1 ) then
                  close(2)
                  read(czeile(30:40),'(i11)') nfsho(ifil) 
                  if ( nfsho(ifil) .ge. 100000 ) nfsho(ifil) = 99999
                  goto 409
               endif
            enddo
         endif
  409    continue
      enddo ! end-of loop ifil=2,nfmx.
  411 continue
      nfil = ifil - 1

c--print title line-----------------------------------------------------
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   2x,'' runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')
      mcodprev = 0
      engyprev = 0.
      thetprev = 0.
      phiaprev = 0.
      obslprev = 0.
      lthiprev = 0
      nplotsh = 0
      ntotal = 0

c--work on all particle data files--------------------------------------
      do  444  ifil=1,nfil
      if ( ifil .gt. 1 ) close(unit=3)
      itaskcer = index(cdata(ifil),'CER')
      itaskdat = index(cdata(ifil),'DAT')
      ntotal = ntotal + 1
      open(unit=3,file=cdata(ifil),status='unknown',form='unformatted')
      read(unit=3,end=424,err=422) (pdata(i),i=1,5733) ! 5733 elements.
      goto 424 
  422 continue
      write(*,*) '     ERROR in reading pdata vector, ifil=',ifil
  424 continue

c - - - - - check on 64 or 32 bit and standard or thinning corsika:
      ibit = 32
      if ( abs(pdata(1)).lt.1.e-6 .and.
     +   ( pdata(2).ge.211284. .and. pdata(2).le.211286. ) ) ibit = 64
      if ( pdata(273).ge.217432. .and. pdata(273).le.217434.) ibit = -64
      if ( pdata(312).ge.217432. .and. pdata(312).le.217434.) ibit = -64
      nbits(ifil) = ibit
      ibit = ibit / 64
      isubr = 312
      if ( pdata(ibit+274).ge.217432. .and. pdata(ibit+274).le.217434.)
     +   isubr=273
      qrunh(1) = prunh
      do  i=2,isubr
         qrunh(i) = pdata(ibit+i)
      enddo
      do  i=isubr+1,isubr*2
         qevth(i-isubr) = pdata(ibit+i)
      enddo
      ntask = ntskf(ifil)
      if ( qevth(11) .lt. 0.01 ) qevth(12) = 0.0

c - - - - - check on primary particle:
      if ( 0 .lt. qevth(3) .and. qevth(3) .lt. 200. ) then
         ip = int(qevth(3))
      else if ( qevth(3) .lt. 6000. ) then
         ip = 200
         qpatext(200) = chemical(int(mod(qevth(3),100.)))//'       '
         if ( qevth(3).eq.201 ) qpatext(200) = ' Deuteron          '
         if ( qevth(3).eq.301 ) qpatext(200) = ' Tritium           '
      else
         write(*,*) '       invalid particle id ',qevth(3)
      endif
      czeile(1:11) = qpatext(ip)(1:11) ! text of primary particle.

c - - - - - get codes for models of intervals of energy - - - - -
      models = 0
      do  i=73,80
         models = models + 10**(80-i) * int(qevth(i))
      enddo

c - - - - - check date of simulation - - - - -
      imont = mod(qevth(45),1.e4) ! including number of day.
      ijahr = qevth(45) * 1.00001e-4
      if ( ijahr < 30 ) ijahr = ijahr + 2000
      if ( ijahr < 100) ijahr = ijahr + 1900
      idate = 10000.*ijahr + imont ! yyyymmdd as 8 digits.
      if ( mod(idate,100) .gt. 31 ) idate = 31 + 100 * int(idate/100)
      if ( mod(idate,100) .eq.  0 ) idate =  1 + idate 

c - - - - - check on 64bit simulation and/or Cerenkov output - - - - -
      if ( nbits(ifil) .eq. 64 ) then
         if ( itaskcer .gt. 0 ) then
            czeile(9:11) = '_C4' 
         else
            czeile(9:11) = '_64'
         endif 
      elseif ( itaskcer .gt. 0 ) then
         czeile(9:11) = '_CE'
      endif

c - - - - - check on thinning parameters for qevth(150)>0 - - - - -
      if ( qevth(150) .gt. 0. ) then ! or qevth(148) .lt. 1. !!!
         if ( qevth(152) .eq. 0. ) qevth(152) = 1.e4
         if ( qevth(151) .lt. 1. ) then
            if ( qevth(148) .eq. 0. ) qevth(148) = 1.
            isla = log10(qevth(4)) + log10(qevth(148)+4.67735141e-34)
            islb = 10 
            qevth(151) = 0.
            do  i=1,isla
              qevth(151) = qevth(151) + 10.d0**(isla-i) * (islb-i) 
            enddo
         endif
         if ( qevth(152) .gt. 1. ) then
            qevth(152) = 0.01234 + log10(qevth(152)*1.e-2)
         else
            qevth(152) = 0. ! qevth(151) = 0.
         endif
      endif
      if ( 0. lt. qevth(148) .and. qevth(148) .lt. 0.1234 ) then
         qevth(148) = log10(qevth(148)+4.67735141e-34)
      else
         qevth(148) = 0.
      endif

c - - - - - print out tabular quantities - - - - -
      if ( ip .gt. 0 ) then
        if ( int(qevth(44)) .ne. ntask ) qevth(44) = ntask
        mcode = int(qevth(3))
        ! - - - - - - mcode=0 for stackin simulations - - - - - - - - -
        if ( mcode .eq. 0 ) mcode = 4
        energy = qevth(4) 
        if ( qevth(4) .gt. 0. ) energy = 9.+log10(qevth(4))
        theta = c180pi*qevth(11)
        phia = c180pi*qevth(12)
        if ( mcode .ne. mcodprev .or. energy .ne. engyprev .or.
     +     theta .ne. thetprev .or. phia .ne. phiaprev .or.
     +    obslev .ne. obslprev .or. lthi .ne. lthiprev ) then
          ! write(*,'(1x)') ! check extra blank lines (if wanted).
        endif
        if ( qevth(159) .lt. 1.e-3 ) then
          if ( qevth(44) .lt. 80000. .or. 80009. .lt. qevth(44) ) then
           ! no PLOTSH simulation, do not print energy cuts: 
           if ( fsize(ifil) .lt. 10000. ) then
            write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.1,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), qevth(148),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
           else
            write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.0,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), qevth(148),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
           endif
          else 
           ! assumed to be a PLOTSH simulation, print energy cuts: 
           nplotsh = nplotsh + 1
           if ( fsize(ifil) .lt. 10000. ) then
            write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.1,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,1p,4e9.1)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), qevth(148),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
     +      ,qevth(61), qevth(62), qevth(63), qevth(64)
           else
            write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.0,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,1p,4e9.1)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), qevth(148),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
     +      ,qevth(61), qevth(62), qevth(63), qevth(64)
           endif 
          endif 
        else
         if ( fsize(ifil) .lt. 10000. ) then
          write(*,'(a11,i5,a1,f7.2,f6.1,f9.3,i4,i8.6,f8.1,
     +    f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f11.6)')
     +    czeile(1:11),int(qevth(3)),chaspec,energy,
     +    c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +    int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +    1.d-2*qevth(7), qevth(148),
     +    qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +    idate, qevth(71), qevth(72)
     +    ,qevth(159)
         else
          write(*,'(a11,i5,a1,f7.2,f6.1,f9.3,i4,i8.6,f8.0,
     +    f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f11.6)')
     +    czeile(1:11),int(qevth(3)),chaspec,energy,
     +    c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +    int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +    1.d-2*qevth(7), qevth(148),
     +    qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +    idate, qevth(71), qevth(72)
     +    ,qevth(159)
         endif 
        endif
      endif
      qpatext(200) = '                   '
      mcodprev = mcode
      engyprev = energy
      thetprev = theta
      phiaprev = phia
      obslprev = obslev
      lthiprev = lthi

c - - - - - end-of loop ifil=1,nfil.
  444 continue 
      close(unit=3)

c--print closing comment lines------------------------------------------
      if ( nplotsh .gt. 0 ) then
        write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +    2x,'' runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +    2x,''verspgm    models   rundate  Xmagn  Zmagn'',3x,
     +    ''hadron-'',2x,''muon-'',2x,''electr-'',1x,''photoncut'',/)')
      else
        write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +    2x,'' runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +    2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')
      endif
  445 continue
      write(*,'(14x,''Total number of files:'',i7)') ntotal

c--print explanation of model digits------------------------------------
      write(*,'(14x,''Appendix `_64`,`_CE`,`_C4`: by 64bit executable'',
     + '', Cerenkov output read, or both was detected, otherwise no'',
     + '' such infos;'')')
      write(*,'(14x,''explanation of column `models`:'',
     + 3x,''(10^7): EGS flag;   (10^6): NKG flag;   (10^5):'',
     + '' lowEnergy flag, 1=gheisha, 2=urqmd,'')')
      write(*,'(14x,''3=fluka;  (10^4): highEnergy, 0=hdpm,'',
     + '' 1=venus, 2=sibyll, 3=qgsjet, 4=dpmjet, 5=nexus, 6=epos;'',
     + 2x,''(10^3): Cerenkov flag;'')')
      write(*,'(14x,''(10^2): Neutrino flag;  (10^1): Curved'',
     + '' flag, 0=standard, 1=opt.Aires, 2=curved;  (10^0):'',
     + '' Computer, 3=unix, 4=macint.;'')')
      write(*,'(14x,''thirad: radial thinning to 10^thirad meter.'')')
      write(*,'(''  '',2(/,''  ''))')

c--end of data----------------------------------------------------------
      stop
      end

c=======================================================================
c
c  s h o w s i m u l i s t . f
c  ---------------------------
c     reading corsika particle data files from a (one name per line)
c     tabular each simulation ([DAT,CER]iiiiii) will be printed out
c     in one line the following quantities:
c         Primary, lg(E), theta, phi, nsh,  runnr, size,
c             obslvme, h1stme, thilev, thiwmax, lg(thirad),
c                 verspgm, models, rundate, Xmagn, Zmagn;
c     on 32bit machines all original 64bit simulations will be detected
c     correctly, Cherenkov files will be marked by `_ce`; particle data
c     files must be available as [DAT,CER]iiiiii[.cher], protocol files
c     as [DAT,CER]iiiiii.lst; otherwise additional conditions must be
c     implemented in this source code.
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check showsimulist.f -o showsimulist
c     ifort -C -O0 -check bounds showsimulist.f -o showsimulist 
c execution:
c     ls -l DAT??????[.cher] | ./showsimulist
c     ls -l CER?????? | ./showsimulist > showsimulist.cherenkov
c-----------------------------------------------------------------------
c output-file:
c     unit=*: tabular output of simulations.
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c         Primary   lg(E)  theta    phi  nsh  runnr    size  obslvme  h1stme   .....
c Iron       5626   16.00    0.0    0.0   1  199079     3.9  1413.82   18350.  .....
c _stck_in_     4   15.09    0.0    0.0   1  169051     2.4   194.00   18765.  .....
c Fluorine   1909   16.00    0.0    0.0   1  199070     3.9  1416.51   22222.  .....
c proton       14   16.00   30.0   -3.3   1  199080     3.4  1429.25   22224.  .....
c Manganese  5525   16.00   30.0   -3.3   1  199082    10.1  1428.13   22222.  .....
c proton _ce   14   17.50   37.9 -138.0   1  000044    90.8  -500.00  -24880.  .....
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program showsimulist

      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (nfmx=50000,nchx=250)

      character cdata(nfmx)*250,cdat*250,csklst*250,czeile*80
      character qpatext(0:200)*19,chemical(0:100)*12,chfmti*5,chaspec*1

      dimension fsize(nfmx),qrunh(312),qevth(312)
      dimension nfsho(nfmx),ntskf(nfmx),nbits(nfmx),nflen(nfmx)
      real pdata(5733)

      logical lstexi

      data chaspec/' '/

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
     c ' Berkelium  ',' Californium',' Einsteinium','            '/

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
     c' K*-               ',' anti K*0          ',' elec neutrino     ',
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
      nfilstop = 0
      obslev = 0.1
      lthi = 0

c--read first file name and fix position of file size-------------------
      ifil = 1
      read(*,'(a)',end=411,err=411) cdat
      ldat = nchx
  401 continue
      ldat = ldat - 1
      if ( cdat(ldat:ldat) .eq. ' ' ) goto 401     
      if ( index(cdat(1:ldat),'-rw') .le. 0 .and.
     +     index(cdat(1:ldat),'-r') .le. 0 ) then
         write(*,'(/,14x,''Use of `ls -1 ...` not supported, but'',
     +      '' `ls -l ...` will be ok;'')')
         write(*,'(14x,''or file access rights do not begin'',
     +      '' with `-rw`; first use `chmod u+w *`;'')')
         goto 445
      endif
      ibl1 = index(cdat,' ')
      ibl2 = index(cdat(ibl1+2:nchx),' ') ! after -rwx-... 2 blanks.
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
      if ( iright-ileft+1 .gt. 12 ) ileft = iright - 11
      write(chfmti(1:5),'(''(i'',i2,'')'')') iright-ileft+1
      read(cdat(ileft:iright-3),chfmti) isize
      fsize(ifil) = 1.d-3 * isize ! now in Mbytes.
      if ( fsize(ifil) .lt. 0.1 ) fsize(ifil) = 0.1
      itaskcer = index(cdat,'CER')
      itaskdat = index(cdat,'DAT')
      itaskchv = index(cdat,'cher')
      itask = -8
      if ( itaskcer .gt. 0 ) then
         itask = itaskcer
      else if ( itaskdat .gt. 0 ) then
         itask = itaskdat
      endif 
      read(cdat(itask+3:itask+8),'(i6)') ntask ! irunnr
      ntskf(ifil) = ntask
c--check structure of file names (with : or / or without):
      ilen = nchx
  402 continue
      ilen = ilen - 1
      if ( cdat(ilen:ilen) .eq. ' ' ) goto 402
      iddot = index(cdat,':')
      isla = ilen - 8
      if ( iddot .gt. 10 ) then
         isla = iddot + 4
      else
         isla = ilen
  403    continue
         isla = isla - 1
         if ( cdat(isla:isla) .ne. ' ' ) goto 403
         isla = isla + 1
      endif 
c--check length of first line in the file list:
      cdata(ifil) = cdat(isla:ilen)
      nflen(ifil) = ilen - isla + 1
      nfsho(ifil) = 1
c--check `DATiiiiii.lst` file to get number of simulated showers:
      if ( index(cdat,'.part') .eq. 0 ) then
         csklst = cdata(ifil)(1:nflen(ifil))//'.lst'
      else if ( index(cdat,'.part') .gt. ilen-5 ) then
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
               if ( nfsho(ifil) .ge. 10000 ) nfsho(ifil) = 9999
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
         itaskchv = index(cdat,'cher')
         itask = -8
         if ( itaskcer .gt. 0 ) then
            itask = itaskcer
         else if ( itaskdat .gt. 0 ) then
            itask = itaskdat
         endif
         read(cdat(itask+3:itask+8),'(i6)') ntask
         ntskf(ifil) = ntask
         ilen = nchx
  406    continue
         ilen = ilen - 1
         if ( cdat(ilen:ilen) .eq. ' ' ) goto 406
         iddot = index(cdat,':')
         isla = ilen - 8
         if ( iddot .gt. 10 ) then
            isla = iddot + 4
         else
            isla = ilen
  407       continue
            isla = isla - 1
            if ( cdat(isla:isla) .ne. ' ' ) goto 407
            isla = isla + 1
         endif
         cdata(ifil) = cdat(isla:ilen)
         nflen(ifil) = ilen - isla + 1
         nfsho(ifil) = 1
         if ( index(cdat,'.part') .eq. 0 ) then
            csklst = cdata(ifil)(1:nflen(ifil))//'.lst'
         else if ( index(cdat,'.part') .gt. ilen-5 ) then
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
                  if ( nfsho(ifil) .ge. 10000 ) nfsho(ifil) = 9999
                  goto 409
               endif
            enddo
         endif
  409    continue
      enddo ! end-of loop ifil=2,nfmx.
  411 continue
      nfil = ifil - 1

c--print title lines----------------------------------------------------
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   2x,'' runnr   sizeM  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')
      mcodprev = 0
      engyprev = 0.
      thetprev = 0.
      phiaprev = 0.
      obslprev = 0.
      lthiprev = 0
      nfilstop = 0

c--work on all particle data files--------------------------------------
      do  444  ifil=1,nfil
      if ( ifil .gt. 1 ) close(unit=3)
      do  i=1,5733
         pdata(i) = -1.
      enddo  
      itaskcer = index(cdata(ifil),'CER')
      itaskdat = index(cdata(ifil),'DAT')
      itaskchv = index(cdat,'cher')
      nfilstop = nfilstop + 1
      open(unit=3,file=cdata(ifil)(1:nflen(ifil)),status='old',
     +     form='unformatted',access='sequential')
      read(unit=3,end=422,err=422) (pdata(i),i=1,5733) ! 5733 elements.
      goto 424 
  422 continue
      ! particle data file empty, continue with next file:
      if ( pdata(1) .eq. -1. ) goto 444
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
      qrunh(1) = 211285.2812500000000000d0
      do  i=2,isubr
         qrunh(i) = pdata(ibit+i)
      enddo
      do  i=1,isubr
         qevth(i) = pdata(ibit+isubr+i)
      enddo
      ntask = ntskf(ifil)
      if ( qevth(11) .lt. 0.01 ) qevth(12) = 0.0
      ip = nint(qevth(3))

c - - - - - check on primary particle:
      if ( 0 .lt. qevth(3) .and. qevth(3) .lt. 200. ) then
         ip = int(qevth(3))
      else if ( qevth(3) .le. 5656. ) then
         ip = 200
         qpatext(200) = chemical(int(mod(qevth(3),100.)))//'       '
         if ( qevth(3).eq.201 ) qpatext(200) = ' Deuteron          '
         if ( qevth(3).eq.301 ) qpatext(200) = ' Tritium           '
      else
         write(*,*) '       invalid particle id ',qevth(3)
      endif

c - - - - - check models and date of simulation - - - - -
      models = 0
      do  i=73,80
         models = models + 10**(80-i) * int(qevth(i))
      enddo
      imont = mod(qevth(45),1.e4)
      ijahr = qevth(45) * 1.00001e-4
      if ( ijahr < 30 ) ijahr = ijahr + 2000
      if ( ijahr < 100) ijahr = ijahr + 1900
      idate = 10000. * ijahr + imont
      if ( mod(idate,100) .gt. 31 ) idate = 31 + 100 * int(idate/100)
      if ( mod(idate,100) .eq.  0 ) idate =  1 + idate 
      if ( qevth(148) .eq. 0. ) qevth(148) = 1.
      czeile(1:11) = qpatext(ip)(1:11)

c - - - - - check on original 64bit simulation and Cherenkov output - - -
      if ( nbits(ifil) .eq. 64 ) then
         if ( itaskcer .gt. 0 ) then
            czeile(9:11) = '_c4' ! should not occur. 
         else
            czeile(9:11) = '_64' ! should not occur.
         endif 
      else if ( itaskcer .gt. 0 .or. itaskchv .gt. 0 ) then
         czeile(9:11) = '_ce'
      endif

c - - - - - check on thinning parameters - - - - -
      if ( qevth(151) .lt. 1. ) then
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
         qevth(152) = 0.
         qevth(151) = 0.
      endif

c - - - - - print out tabular of shower quantities - - - - -
      if ( ip .gt. 0 ) then
         if ( int(qevth(44)) .ne. ntask ) qevth(44) = ntask
         mcode = int(qevth(3))
         ! - - - - - - mcode=0 for stackin simulations - - -
         if ( mcode .eq. 0 ) mcode = 4
         energy = qevth(4) 
         if ( qevth(4) .gt. 0. ) energy = 9.+log10(qevth(4))
         theta = c180pi*qevth(11)
         phia = c180pi*qevth(12)
         if ( mcode .ne. mcodprev .or. energy .ne. engyprev .or.
     +      theta .ne. thetprev .or. phia   .ne. phiaprev .or.
     +      obslev .ne. obslprev .or. lthi   .ne. lthiprev ) then
            ! write(*,'(1x)') ! check extra blank lines.
         endif

         if ( fsize(ifil) .lt. 10000. ) then ! < 10 GBytes
            write(*,'(a11,i5,a1,f7.2,f6.1,f8.2,i5,i8.6,f8.1,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,i10.8,i10,2f7.2)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
         else
            write(*,'(a11,i5,a1,f7.2,f6.1,f8.2,i5,i8.6,f8.0,
     +      f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,i10.8,i10,2f7.2)')
     +      czeile(1:11),int(qevth(3)),chaspec,energy,
     +      c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +      int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +      1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +      qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +      idate, qevth(71), qevth(72)
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
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   2x,'' runnr   sizeM  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')
  445 continue
      write(*,'(14x,''Total number of files:'',i7)') nfilstop

c--print explanation of model digits------------------------------------
      write(*,'(14x,''Appendix `_ce`: Cherenkov output read,'',
     + '' otherwise no such extra info (i.e. only blanks).'')')
      write(*,'(14x,''Explanation of column `models`:'',
     + 2x,''(10^7): EGS flag;   (10^6): NKG flag;   (10^5):'',
     + '' lowEnergy flag, 1=gheisha, 2=urqmd,'')')
      write(*,'(14x,''3=fluka;  (10^4): highEnergy, 0=hdpm,'',
     + '' 1=venus, 2=sibyll, 3=qgsjet, 4=dpmjet, 5=nexus, 6=epos;'',
     + 2x,''(10^3): Cherenkov flag;'')')
      write(*,'(14x,''(10^2): Neutrino flag;   (10^1): Curved'',
     + '' flag, 0=standard, 2=curved;   (10^0):'',
     + '' Computer, 3=unix, 4=macintosh.'')')
      write(*,'(14x,''thirad: radial thinning to 10^thirad meter.'')')
      write(*,'(''  '',2(/,''  ''))')

c--end of data----------------------------------------------------------
      stop
      end

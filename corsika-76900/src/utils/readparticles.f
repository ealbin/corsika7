c=======================================================================
c
c  r e a d p a r t i c l e s . f
c  ----------------------------- 
c  read a number of particle data files of a corsika air shower
c  (parallel) simulation also incl. multi-thinning option and
c  print tabular of counted particles;
c-----------------------------------------------------------------------
c CompLink:
c      gfortran -O0 -fbounds-check readparticles.f -o readparticles 
c      ifort -C -O0 -check bounds readparticles.f -o readparticles
c RunProg:
c      ./readparticles.f < readparticles.i > readparticles.out111888
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c      input-files:
c           unit=*: number of showers and file name(s):
c           ------------------------------------------------
c                       1          'total_number_of_showers'
c                       1          'total_number_of_files'
c            '/lxdata/d2lx14/joe/DAT045216'
c                       1          '_showers_of_this_file_'
c           ------------------------------------------------
c           unit=3: current corsika particle data file.
c     output-files: 
c           unit=*: protocol output.
c           unit=9: ascii file named DATnnnnnn.ascithin or ....ascistnd.
c Remark: former true 64-bit corsika simulations will be detected.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readparticles
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cout*250,crunh*4,cevte*4,qpatext(1:200)*19
      character cdata(5000)*250,cdat*250,cblk*250
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,irec,isho,isub,idat,ilen,ifil,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      common /ctext/cdata,qpatext
      data crunh/'RUNH'/,cevte/'EVTE'/

c--initialize names of particles----------------------------------------
      data (qpatext(i),i=1,60)/
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
     c' anti Xi+          ',' anti Omega+       ',8*'                 ',
     c' Deuteron          ',' Tritium           ',' Helium3           ',
     c' Helium4           ',' addi muon+        ',' addi muon-        ',
     c 3*'                ',' omega             ',' rho0              ',
     c' rho+              ',' rho-              ',' Delta++           ',
     c' Delta+            ',' Delta0            ',' Delta-            ',
     c' anti Delta--      ',' anti Delta-       ',' anti Delta0       '/
      data (qpatext(i),i=61,129)/
     c' anti Delta+       ',' K*0               ',' K*+               ',
     c' K*-               ',' anti K*0          ',' electron neutrino ',
     c' anti elec neutrino',' muon neutrino     ',' anti muon neutrino',
     c'                   ',' eta=>2*gamma      ',' eta=>3*pi0        ',
     c' eta=>pi+pi-pi0    ',' eta=>pi+pi-gamma  ',' addi muon+        ',
     c' addi muon-        ','                   ',37*'                ',
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

c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      isho = 0
      iend = 0
      irec = 0
      irun = 0
 
c--read run parameters including file names-----------------------------
      read(*,*,end=498,err=498) lsho
      read(*,*,end=498,err=498) nfil
      ! write(*,'(22x,''total number of showers:'',i7)') lsho
      ! write(*,'(22x,''total number of files  :'',i7)') nfil
      if ( nfil .gt. 100 ) goto 498
      do  ifil=1,nfil
         read(*,*,end=498,err=498) cdat
         read(*,*,end=498,err=498) nfsh
         cdata(ifil) = cdat
      enddo
      idat = index(cdata(1),'DAT')
      read(cdata(1)(idat+3:idat+8),'(i6)') irun
 
c----------one or more showers in big disk file-------------------------
      write(*,'(1x,6(''_ ''),''readparticles.f '',24('' _''))')
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      ilost = 0
      iswit = 0
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
* - - - - - - read data record with lenstnd words - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      ishift = 2
      itype = 2
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenstnd)
      close(unit=3)
* - - - - - - check on reading 32bit simulation on 64bit machine:
      if ( 211285.2 .lt. pdata(1) .and.
     +     pdata(1) .lt. 211285.4 ) then
         ishift = 0
      elseif ( 211285.2 .lt. pdata(2) .and.
     +     pdata(2) .lt. 211285.4 ) then
         ishift = -1
         do  i=1,936-1
            pdata(i) = pdata(i+1)
         enddo
      elseif ( 2.0202 .lt. pdata(3) .and. pdata(3) .lt. 9.9999 ) then 
         ! check version number instead of testing (273+1) and (312+1).
         ishift = 1
         do  i=936,2,-1
            pdata(i) = pdata(i-1)
         enddo
         pdata(1) = prunh ! 211285.281 
      endif
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         cout = cdata(ifil)(idat:ilen)//'.ascistnd'
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
         cout = cdata(ifil)(idat:ilen)//'.ascithin'
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 498
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      qdate = dble(pdata(3))
      if ( qdate .lt. 765432.d0 ) then
         qdate = 20000000.d0 + qdate
      else
         qdate = 19000000.d0 + qdate
      endif
      write(*,'(/,6x,''runnr ='',i7.6,4x,''log10(E)='',f10.4,
     +   4x,''theta='',f8.2,3x,''phi='',f9.2)') irun,
     +   log10(1.d9*pdata(lenblk+4)),57.29578*pdata(lenblk+11),
     +   57.29578*pdata(lenblk+12)
      write(*,'(/,6x,''primary ='',i5,4x,''obslev(m)='',f9.2,4x,
     +   ''date='',i9,3x,''vers='',f8.5)') nint(pdata(3+lenblk)),
     +   pdata(48+lenblk)*1.d-2,20000000+nint(pdata(3)),pdata(4)
      do  i=1,936 ! keep three subblocks including some particles.
         qdata(i) = pdata(i)
      enddo
      mthin = 0
      do  i=575,936,7
         if ( int(qdata(i)/1000.) .eq. 8888. ) mthin = mthin + 1
      enddo
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,12x,''`thinning run`'',5x,a,/)')
     +      cdata(ifil)(1:ilen)
      else
         if ( mthin .eq. 0 ) then
            write(*,'(/,12x''`standard run`'',5x,a,/)')
     +         cdata(ifil)(1:ilen)
         else
            write(*,'(/,12x,''`multi-thin run`'',5x,a,/)')
     +         cdata(ifil)(1:ilen)
         endif
      endif
      if ( ishift .eq. -1 ) goto 443 ! should not occur. 
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
* - - - - - - check on original 64bit simulation:
      pdatone = pdata(1)
      if ( irec .gt. 1 .and. ( 3.21d-41 .lt. pdatone
     +                 .and. pdatone .lt. 3.23d-41 ) ) goto 443
      if (mod(irec,8000) .eq. 0) write(*,'(27x,''irec ='',i8)') irec
* - - - - - - - - - - -
      if ( ishift .eq. 1 ) then
* - - - - - - shift data corresp. to a 32bit simul. on a 64bit machine:
         do  i=lenrec,2,-1
            pdata(i) = pdata(i-1)
         enddo
         if ( irec .gt. 1 ) then
            pdata(1) = 177177.1
            ipart = int(pdata(lenpar+1)*1.000001d-3)
            if ( ipart .ge. 1 .and. ipart .le. 3 )
     +         pdata(1) = 1000. + mod(pdata(lenpar+1),1000.)
            if ( ipart .eq. 0 ) pdata(1) = 1091. ! fixed generation.
            if ( ipart .eq. 5 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 6 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 75 ) pdata(1) = 177177.1
            if ( ipart .eq. 76 ) pdata(1) = 177177.1
            if ( 3301.2 .lt. pdata(lenblk+1) .and.
     +       pdata(lenblk+1) .lt. 3301.4 ) pdata(1) = pevte ! 3397.39185
            if ( pdata(1) .eq. 177177.1 ) then
               ilost = ilost + 1
            else
               iswit = iswit + 1
            endif
         else ! irec = 1.
            pdata(1) = prunh ! 211285.281
         endif
      elseif ( ishift .eq. -1 ) then
* - - - - - - shift data corresp. to a 64bit simul. on a 32bit machine:
         do  i=1,lenrec-1
            pdata(i) = pdata(i+1)
         enddo
         pdata(lenrec) = -pdata(lenrec-lenpar) 
         ! check negative elements carefully.
      endif
* - - - - - - - - - - -
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
* - - - - - - end of corsika particle data file reached: 
      ! write(*,'(27x,''irec ='',i8)') irec
  443 continue
* - - - - - - distinguish contents of pdata vector:
      if ( ishift .eq. 1 ) then
         write(*,*) '    ______________________________________________'
         write(*,*) '    WARNING: will be better using a 32bit machine.'
         if ( ilost+iswit .gt. 0 ) then
            if ( ilost .ge. 0 )
     +         write(*,*) '    minimal lost particles =',ilost
            if ( iswit .ge. 0 )
     +         write(*,*) '    max switched particles =',iswit
         endif
      elseif ( ishift .eq. -1 ) then
**      write(*,'(5x,''irec=1.01'',f11.1,5f15.1)') (pdata(i),i=1,6)
**      write(*,'(10x,''1.02'',f11.1,5f15.1)') 
**   +       (pdata(i),i=lenblk+1,lenblk+6)
**      write(*,'(10x,''1.03'',f11.1,5f15.1)')
**   +       (pdata(i),i=lenblk*2+1,lenblk*2+5)
**      write(*,*) '             ',(pdata(i),i=lenblk*2+lenpar+1,
**   +                                         lenblk*2+lenpar+5)
        write(*,*) '    irec=1.01',(pdata(i),i=1,6)
        write(*,*) '         1.02',(pdata(i),i=lenblk+1,lenblk+6)
        write(*,*) '         1.03',(pdata(i),i=lenblk*2+1,lenblk*2+5)
        write(*,*) '             ',(pdata(i),i=lenblk*2+lenpar+1,
     +                                         lenblk*2+lenpar+5)
        write(*,*) '    _______________________________________________'
        write(*,*) '    fortran: impossible to read a 64bit simulation;'
        write(*,*) '        use new 32bit particle data files.'
        open(unit=9,file=cout,status='unknown',form='formatted')
        close(unit=9)
      elseif ( ishift .eq. 0 ) then
**       write(*,*) '    ______________________________________________'
**       write(*,*) '    WARNING: last record of the ascii output file '
**       write(*,*) '        may be incompletely filled with zeroes as '
**       write(*,*) '        compared to the exactly last record of the'
**       write(*,*) '        corsika particle data file.'
**       write(*,*) '    INFO: check record length 22940',
**   +             ' of `standard` simulation'
**       write(*,*) '        and other record length 26216',
**   +             ' of `thinning` simulation.'
      endif
  444 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '       irec =',irec
      inull = 0
      do  l=1,5733,7
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,5733,273
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      write(*,*)
     +'    ____________________________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      write(*,'(9x,''estimated number of particles'',f15.0,
     +   '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull  
      write(*,*) '        subblocks EVTE and RUNE missing!'
      goto 499 
  498 continue
      write(*,*) '    _________________________________________________'
      write(*,*) '    ERROR: missing some input parameters.'
  499 continue
c--print out tabulars---------------------------------------------------
      npot = 16 ! logarithmic bins per decade.
      fpot = 10.d0**(1.d0/npot) ! factor per bin.
      if ( npot .lt. 0 ) then
      write(*,'(10x,''photons'')')
      write(*,'(10x,''-------'')')
      do  jp=-4,4
         ja = npot * (jp+4)
         write(*,'(/,1p,e7.0,1x,0p,15f8.5)') 10.d0**jp,(fpot**i,i=1,15)
         write(*,'(16i8)') (jpartic(i,1),i=ja,ja+npot-1)
      enddo
      write(*,'(10x,''positrons'')')
      write(*,'(10x,''---------'')')
      do  jp=-4,4
         ja = npot * (jp+4)
         write(*,'(/,1p,e7.0,1x,0p,15f8.5)') 10.d0**jp,(fpot**i,i=1,15)
         write(*,'(16i8)') (jpartic(i,2),i=ja,ja+npot-1)
      enddo
      write(*,'(10x,''electrons'')')
      write(*,'(10x,''---------'')')
      do  jp=-4,4
         ja = npot * (jp+4)
         write(*,'(/,1p,e7.0,1x,0p,15f8.5)') 10.d0**jp,(fpot**i,i=1,15)
         write(*,'(16i8)') (jpartic(i,3),i=ja,ja+npot-1)
      enddo
      endif
* - - - - - - print total counts of particles (also with weights):
*     write(*,'(/,6x,''code name'',18x,''particles '',a)')
*    +   cdata(nfil)(idat:ilen)
*     write(*,'(6x,''---- ----'',18x,''--------- '',16(''-''))')
*     do  l=1,46 ! 50
*        if ( qpartic(l) .gt. 0. ) write(*,'(i10,a19,f13.0)')
*    +      l,qpatext(l),qpartic(l)
*     enddo
*     if ( lenblk .eq. 312 ) then ! prints particles incl. weights.
*       write(*,'(36x,''counts without weights.'')')
*       write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
*       do  l=1,46 ! 50
*         if ( qpartic(l+100) .gt. 0. ) write(*,'(i10,a19,f14.1)')
*    +       l,qpatext(l),qpartic(l+100)
*       enddo
*       write(*,'(36x,''counts including weights.'')')
*     endif
*     write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
*     write(*,'(2x)')
      stop
      end
c=======================================================================
c
c     block data initialization
c
c=======================================================================

      block data blinit

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)

      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      character cdata(5000)*250,qpatext(1:200)*19
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,irec,isho,isub,idat,ilen,ifil,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      common /ctext/cdata,qpatext
      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data parmas/
     + 2*0.0000d0 , 0.000511d0 , 0.000511d0 , 0.000000d0 , 0.105658d0,
     + 0.105658d0 , 0.134973d0 , 0.139568d0 , 0.139568d0 , 0.497671d0,
     + 0.493646d0 , 0.493646d0 , 0.939566d0 , 0.938272d0 , 0.938272d0,
     + 0.497671d0 , 0.5488d0   , 1.11563d0  , 1.18937d0  , 1.192550d0,
     + 1.19743d0  , 1.3149d0   , 1.32132d0  , 1.67243d0  , 0.939566d0,
     + 1.11563d0  , 1.18937d0  , 1.19255d0  , 1.19743d0  , 1.3149d0,
     + 1.32132d0  , 1.67243d0  , 1.7841d0   , 1.7841d0   , 1.8693d0,
     + 1.8693d0   , 1.8645d0   , 1.8645d0   , 1.9693d0   , 1.9693d0,
     + 2.2852d0   , 80.6d0     , 80.6d0     , 91.161d0   , 1.8770d0,
     + 2.817d0    , 3.755d0    , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.7669d0   , 0.7681d0   , 0.7681d0   , 1.2309d0   , 1.2323d0,
     + 1.2336d0   , 1.2349d0   , 1.2309d0   , 1.2323d0   , 1.2336d0,
     + 1.2349d0   , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.105658d0,
     + 0.105658d0 , 0.0d0      , 0.0d0      , 0.0d0      , 22*0.0d0/
      end
c=======================================================================
c
c     analyze the contents of all 21 subblocks in a record
c
c=======================================================================
c
c        primary    energy    runnumb   simnrsh     #hadr    #muons
c        #gammas     #elec   #elecnkg    obslv    theta    phi
c         h1wkm    h1wgcmq  ( emymax    ehasum    ehamax    nha )
c
c=======================================================================
 
      subroutine blwork(iend,lenrec)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)
 
      double precision parmas(0:101),phead(30)
      character cdata(5000)*250,qpatext(1:200)*19
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin)
      common /integ/lobs,lsho,irec,isho,isub,idat,ilen,ifil,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      common /ctext/cdata,qpatext
      isub = 0
      lenblk = lenrec / 21
      lenpar = lenblk / 39

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      isub = isub + 1
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------
         lobs = int(pdata(lia+4))
         if (isho.le.0) then
            phead(6) = lobs
            do  908  i=1,lobs
               phead(6+i) = 1.d-2 * pdata(lia+4+i)
  908       continue
            if (lobs.ge.1.and.lobs.le.10) then
               if (lobs.gt.1) write(*,
     +            '(6x,''observation level:'',f10.2,'' meter'')')
     +            phead(6+lobs)
            else
               write(*,'(6x,''nr of observation level undefined.'')')
               iend = iend + 9876543
               goto 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.d-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
      elseif ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
         if (iend.lt.0) goto 948
         isho = isho + 1
c- - - - - - - - - - simulated shower number - - - - - - - - - - - - -
         phead(1) = pdata(lia+1)
c - - - - - - - - - - height of first interaction (km)  - - - - - - - -
         phead(2) = 1.d-5 * pdata(lia+6)
c - - - - - - - - - - height of first interaction (grams/cm2)- - - - -
         pdata(lia+4) = real(thickgr(dble(pdata(lia+6))))
         phead(3) = pdata(lia+4)
c - - - - - - - - - - primary particle code - - - - - - - - - - - - - -
         phead(4) = pdata(lia+2)
c - - - - - - - - - - primary particle energy in gev  - - - - - - - - -
         phead(5) = pdata(lia+3)
c - - - - - - - - - - - phi angle this shower in radian - - - - - - - -
         phead(19) = pdata(lia+11)
c - - - - - - - - - - theta angle this shower in radian - - - - - - - -
         phead(20) = pdata(lia+10)
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = real(timev0)
c----------------subblock longi information-----------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
      ! if longi info and tabular in particle data subblocks.
c----------------subblock event end-------------------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then 
         write(*,'(12x,''isub='',i3,'' (evte)    irec='','//
     +      'i9,''    isho='',i7)') isub,irec,isho
         if (iend.lt.0) then ! ignore first showers.
            iend = iend + 1
            isho = iend
            goto 948
         endif
         write(*,'(/,6x,''code name'',18x,''particles '',a)')
     +      cdata(ifil)(idat:ilen)
         write(*,'(6x,''---- ----'',18x,''--------- '',16(''-''))')
         do  l=1,46 ! 50
            if ( qpartic(l) .gt. 0. ) write(*,'(i10,a19,f13.0)')
     +         l,qpatext(l),qpartic(l)
         enddo
         if ( lenblk .eq. 312 ) then ! prints particles incl. weights.
           write(*,'(36x,''counts without weights.'')')
           write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
           do  l=1,46 ! 50
             if ( qpartic(l+100) .gt. 0. ) write(*,'(i10,a19,f14.1)')
     +          l,qpatext(l),qpartic(l+100)
           enddo
           write(*,'(36x,''counts including weights.'')')
         endif
         write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
         write(*,'(2x)')
         do  l=0,400
            qpartic(l) = 0.d0
         enddo
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + lenpar + isho
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
            icode = int(1.d-3*(pdata(l)*1.0000001))
            if ( 0 .lt. icode .and. icode .le. 32 ) then
               equad = parmas(icode)*parmas(icode)
     +               + pdata(l+1)*pdata(l+1)
     +               + pdata(l+2)*pdata(l+2) + pdata(l+3)*pdata(l+3)
               if ( equad .gt. 0 ) then 
                  epart = sqrt(equad) 
                  eplog = log10(epart) + 4. ! multiply by 10^4
                  jlog = int(eplog*16.d0) ! 16 log bins per decade
                  if ( jlog .lt. 0 ) jlog = 0
                  if ( jlog .gt. 201 ) jlog = 201 
                  jpartic(jlog,icode) = jpartic(jlog,icode) + 1 
               endif  
            endif
            if ( icode .gt. 50 ) then
               if ( icode .eq. 75 ) then
                  icode = 45
               else if ( icode .eq. 76 ) then
                  icode = 46
               else if ( icode .eq. 201 ) then
                  icode = 41
               else if ( icode .eq. 301 ) then
                  icode = 42
               else if ( icode .eq. 302 ) then
                  icode = 43
               else if ( icode .eq. 402 ) then
                  icode = 44
               else if ( icode .gt. 200 .and. icode .le. 987 ) then
                  icode = 4 ! count `exotic` nuclei.
               else
                  if ( icode .ge. 85 .and. icode .le. 96 ) then
                     icode = 47 ! ignore 85,86,95,96.
                  endif
                  if ( icode .eq. 8888 ) icode = 50
               endif
            endif
            if ( icode .gt. 0 ) then
               qpartic(icode) = qpartic(icode) + 1.
               qpartic(icode+100) = qpartic(icode+100)+pdata(l+lenpar-1)
            endif
  917    continue
      endif
      if (iend.gt.0) goto 949
  948 continue
 
c-----------all records analyzed or end of showers reached--------------
  949 continue
      return
      end
c=======================================================================
 
      double precision function thickgr( arg )
 
c-----------------------------------------------------------------------
c  calculates thickness (g/cm**2) of atmosphere depending on height (cm)
c  argument:    arg    = height in cm
c-----------------------------------------------------------------------
 
      double precision aatm(5),batm(5),catm(5),arg
      common /atmos/aatm,batm,catm
 
      if     ( arg .lt. 4.d5 ) then
         thickgr = aatm(1) + batm(1) * exp( -arg/catm(1) )
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg/catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg/catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg/catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif
 
      return
      end

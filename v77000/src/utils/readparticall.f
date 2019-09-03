c=======================================================================
c
c  r e a d p a r t i c a l l . f
c  ----------------------------- 
c     read all particle data files of a (parallel) corsika simulation 
c     and print tabular of counted particles. 
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check readparticall.f -o readparticall 
c     ifort -C -O0 -check bounds readparticall.f -o readparticall
c-----------------------------------------------------------------------
c RunProg:
c     ls -1 DAT012345 > readparticall.i012345  
c     ls -1 DAT000182-0????? > readparticall.i000182
c     ./readparticall < readparticall.i012345
c     ./readparticall < readparticall.i000182 > readparticall.out000182   
c-----------------------------------------------------------------------
c input-files:
c           unit=*: names of particle data files.
c           unit=3: current corsika particle data file.
c output-file: 
c           unit=*: protocol output.
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readparticall
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cdata(8000)*250,cdat*250
      character crunh*4,cevte*4,qpatext(1:200)*19
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,nfil,ifil,irec,isho,isub,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
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

c--initialize some more quantities--------------------------------------
      cdat='                                                           '
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      isho = 0
      iend = 0
      irec = 0
 
c--read all names of particle data files--------------------------------
      write(*,'(1x,6(''_ ''),''readparticall.f '',24('' _''))')
      do  ifil=1,8000
         read(*,'(a)',end=408,err=407) cdat
         ilen = index(cdat,' ') - 1
         cdata(ifil) = cdat(1:ilen)
      enddo
      write(*,'(6x,''maximum number of particle data files reached.'')')
      goto 409
  407 continue
      write(*,'(6x,''ERROR in reading names of particle data files.'')')
      goto 409
  408 continue
      ! write(*,'(6x,''END condition reading particle data files.'')')
  409 continue
      nfil = ifil - 1
      ifil = 1

c--check and fix length of file name, test simulation-------------------
      lfil = 250+1
  414 continue
      lfil = lfil - 1
      if ( cdata(1)(lfil:lfil) .eq. ' ' ) goto 414
      ilen = lfil
      idat = index(cdata(1)(1:ilen),'DAT')
      if ( idat .ge. 1 ) then
         read(cdata(1)(idat+3:idat+8),'(i6)') irun
      else
         if ( index(cdata(1)(1:ilen),'_') .eq. ilen-2 ) then
            irun = 987654
         else
            irun = 987987
         endif
      endif

c--get first record of first particle data file:
      open(unit=3,file=cdata(1)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(unit=3,err=496,end=496) (pdata(i),i=1,lenstnd)
      close(unit=3)
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39
      do  i=1,936 ! keep three subblocks including some particles.
         qdata(i) = pdata(i)
      enddo
      mthin = 0
      do  i=575,936,7
         if ( int(qdata(i)/1000.) .eq. 8888. ) mthin = mthin + 1
      enddo
      qdate = dble(pdata(3))
      if ( qdate .lt. 765432.d0 ) then
         qdate = 20000000.d0 + qdate
      else
         qdate = 19000000.d0 + qdate
      endif
      write(*,'(/,6x,''runnr ='',i7.6,4x,''log10(E)='',f10.4,
     +   4x,''theta='',f8.2,3x,''phi='',f9.2)') int(pdata(2)),
     +   log10(1.d9*pdata(lenblk+4)),57.29578*pdata(lenblk+11),
     +   57.29578*pdata(lenblk+12)
      write(*,'(/,6x,''primary ='',i5,4x,''obslev(m)='',f9.2,4x,
     +   ''date='',i9,3x,''vers='',f8.5)') nint(pdata(3+lenblk)),
     +   pdata(48+lenblk)*1.d-2,20000000+nint(pdata(3)),pdata(4)
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,12x,''`thinning run`'',5x,a,/)')
     +      cdata(1)(1:ilen)
      else
         if ( mthin .eq. 0 ) then
            write(*,'(/,12x''`standard run`'',5x,a,/)')
     +         cdata(1)(1:ilen)
         else
            write(*,'(/,12x,''`multi-thin run`'',5x,a,/)')
     +         cdata(1)(1:ilen)
         endif
      endif
 
c----------many big particle data files for one shower------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      isho = 0
      iend = 0
      irec = 0
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      if ( ifil .eq. nfil .and. mod(irec,4000) .eq. 0 )
     +   write(*,'(27x,''irec ='',i8)') irec
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
* - - - - - - end of corsika particle data file reached: 
      write(*,'(27x,''irec ='',i8)') irec
  444 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      inull = 0
      do  l=1,lenstnd,lenpar
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,lenstnd,lenblk
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      write(*,*)
     +'    ____________________________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      write(*,'(9x,''estimated number of particles'',f15.0,
     +   '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull  
      write(*,*) '        subblocks EVTE and RUNE missing!'
  499 continue
c--print out tabulars---------------------------------------------------
      npot = 16 ! logarithmic bins per decade.
      fpot = 10.d0**(1.d0/npot) ! factor per bin.
      if ( fpot .lt. 0. ) then ! optionally print out.
         ! write(*,'(2x)')
         write(*,'(10x,''photons'')')
         write(*,'(10x,''-------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
           write(*,'(16i8)') (jpartic(i,1),i=ja,ja+npot-1)
         enddo
         write(*,'(10x,''positrons'')')
         write(*,'(10x,''---------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
            write(*,'(16i8)') (jpartic(i,2),i=ja,ja+npot-1)
         enddo
         write(*,'(10x,''electrons'')')
         write(*,'(10x,''---------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
            write(*,'(16i8)') (jpartic(i,3),i=ja,ja+npot-1)
         enddo
      endif
      ! print total counts of particles (and also with weights):
      write(*,'(/,6x,''code name'',18x,''particles '',a)')
     +   cdata(nfil)(1:ilen)
      write(*,'(6x,''---- ----'',18x,''--------- '',16(''-''))')
      do  l=1,46
         if ( qpartic(l) .gt. 0. ) write(*,'(i10,a19,f13.0)')
     +      l,qpatext(l),qpartic(l)
      enddo
      write(*,'(36x,''counts without weights.'')')
      if ( lenblk .eq. 312 ) then
        write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
        do  l=1,46
          if ( qpartic(l+200) .gt. 0. ) write(*,'(i10,a19,f14.1)')
     +       l,qpatext(l),qpartic(l+200)
        enddo
        write(*,'(36x,''counts including weights.'')')
      endif
      write(*,'(6x,''---- ---- ----'',13x,''--------- '')')
      write(*,'(2x)')
c--end-of program.
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
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,nfil,ifil,irec,isho,isub,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
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
      dimension qpartic(0:400),jpartic(0:201,32)
      real pdata(lenthin)
      common /integ/lobs,lsho,nfil,ifil,irec,isho,isub,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
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
               phead(6+i) = 1.e-2 * pdata(lia+4+i)
  908       continue
            if (lobs.ge.1.and.lobs.le.10) then
               if (lobs.gt.1) write(*,
     +            '(6x,''observation level:'',f10.2,'' meter'')')
     +            phead(6+lobs)
            else
               write(*,'(6x,''nr of observation level undefined.'')')
               iend = iend + 987654321
               goto 948 ! 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.e-2 * pdata(lia+lobs)
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
         phead(2) = 1.e-5 * pdata(lia+6)
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
         write(*,*) '       isub=',isub,' (evte)    irec=',irec
     +      ,'    isho =',isho
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         ! if ( ifil .eq. nfil ) write(*,'(27x,''irec='',i9)') irec
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
           icode = int(1.d-3*(pdata(l)*1.0000001))
           if ( 0 .lt. icode .and. icode .le. 32 ) then
             equad = parmas(icode)*parmas(icode) + pdata(l+1)*pdata(l+1)
     +             + pdata(l+2)*pdata(l+2) + pdata(l+3)*pdata(l+3)
             if ( equad .gt. 0. ) then 
               epart = sqrt(equad) 
               eplog = log10(epart) + 4. ! multiply by 10000.
               jlog = int(eplog*16.d0) ! 16 log bins per decade
               if ( jlog .le. 0 ) jlog = 0
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
              else if ( icode .gt. 200 .and. icode .lt. 987 ) then
                 write(*,*) '   a    icode = ',icode
                 icode = 4 ! count `exotic` nuclei.
              else
                 if ( icode .ge. 85 .and. icode .le. 96 ) then
                    icode = 47 ! ignore 85,86,95,96.
                 else 
                    if ( icode .eq. 8888 ) then
                       icode = 50
                    else
                       write(*,*) '   b    icode = ',icode
                    endif
                 endif
              endif
           endif
           if ( icode .gt. 0 .and. icode .le. 46 ) then
              qpartic(icode) = qpartic(icode) + 1.
              qpartic(icode+200) = qpartic(icode+200)+pdata(l+lenpar-1)
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

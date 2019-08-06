c=======================================================================
c
c  r e a d m t h i n p r t c l s . f
c  --------------------------------- 
c  read a number of particle data files of a parallel simulation 
c  with multi-thinning option and split them into thinning levels; 
c-----------------------------------------------------------------------
c compile+link:
c     gfortran -O0 -fbounds-check readmthinprtcls.f -o readmthinprtcls 
c     ifort -C -O0 -check bounds readmthinprtcls.f -o readmthinprtcls
c-----------------------------------------------------------------------
c select-files:
c  echo "DAT151503" > readmthinmuons.inpmthi
c  ls -1 DAT001515-* | grep t -v | grep n -v > readmthinprtcls.inpmthi
c  ls -1 DAT001522-0????? > readmthinprtcls.inpmthi
c running:
c     ./readmthinprtcls.sh DAT151503
c     ./readmthinprtcls < readmthinprtcls.inpmthi
c-----------------------------------------------------------------------
c input-files:
c              unit=*: names of particle data files.
c              unit=3: current corsika particle data file.
c output-file: 
c              unit=*: protocol output.
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readmthinprtcls
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cdata(25000)*120,cdat*120,cblk*120
      character qpatext(200)*19,cfilemuons*120
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(40),qparhit(6,20,0:50)
      double precision qpartic(0:200,0:6)
      double precision qparmth(0:200,6),qparsum(0:200,0:6)
      double precision xcore(20),ycore(20),qparh3m(6,20,0:50)
      double precision qnpl1h(0:65),qnpl0h(0:65),rcore(20)
      integer ldata(25000)
      real pdata(lenthin),qdata(936),bdata(624),edata(624)
      common /atmos/aatm,batm,catm
      common /integ/mtype8,irec,isho,isub,ifil,nfil,mcores
     +  ,mlevel,mthins,mprt,qindepar(6),qweight1(6),qwei8888(6)
     +  ,qnpl1h,qnpl0h
      common /utabl/pdata,parmas,phead,qpartic,bdata,edata,qparh3m
     +  ,qparmth,qparsum,qparhit,xmcoord,ymcoord,dradius,qradius
     +  ,xshcore,yshcore,xdistm,ydistm,hdistm,xcore,ycore,rcore

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
      cblk='                                                  '
      cdat=cblk
      xdistm =  750.d0 ! instead of 1500.d0
      hdistm = 0.5d0 * xdistm
      ydistm = sqrt(3.d0) * hdistm
      qwei8888(1) = 0.d0 ! count all lines with 8888 codes.
      qweight1(1) = 0.d0 ! count all lines with at least one weight >= 1.
      qindepar(1) = 0.d0 ! total number of weights (of each mlv) >= 1.
      mthins = 0
 
c--read all names of particle data files--------------------------------
      write(*,'(1x,6(''_ ''),''readmthinprtcls.f '',23('' _''))')
      do  ifil=1,25000
         read(*,'(a29)',end=428,err=427) cdat
         ilen = index(cdat,' ') - 1
         cdata(ifil) = cdat(1:ilen)
         ldata(ifil) = ilen
      enddo
      write(*,'(6x,''maximum number of particle data files reached.'')')
      goto 429
  427 continue
      write(*,'(6x,''ERROR in reading names of particle data files.'')')
      goto 429
  428 continue
      ! write(*,'(6x,''END condition reading particle data files.'')')
  429 continue
      nfil = ifil - 1

c--get first record of first particle data file:
      idat = index(cdata(1),'DAT')
      ilen = ldata(1)
      cfilemuons=cdata(1)(idat:ilen)//'.muons'
      read(cdata(1)(idat+3:idat+8),'(i6)') irun
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
         write(*,'(/,12x,''`thinning run`'')')
         mthin = -1234
      else
         if ( mthin .eq. 0 ) then
            write(*,'(/,12x,''`standard run`'')')
         else if ( mthin .gt. 0 ) then
            write(*,'(/,12x,''`multi-thin run`'')')
         endif
      endif
      mwarn = 0

c--open muons tabular file (ascii):
      open(unit=8,file=cfilemuons,access='sequential',
     +     form='formatted',status='unknown')
      
c--check particle data subblocks on multi-thinning identification:
      mtype8 = 0
      do  l=1+2*lenblk,5*lenblk,lenpar
         icode = int(1.d-3*(pdata(l)))
         if ( icode .eq. 8888 ) mtype8 = mtype8 + 1 
      enddo

c--get total number of multi-thinning levels:
      write(*,'(2x)')
      mthins = int(qdata(450))
      if ( mthins .gt. 0 ) then
        write(*,'(12x,''mthlev'',16x,''maxweight'',15x,''maxweight'')')
        do  mlv=1,mthins
          write(*,'(12x,''mthi='',i1,1x,1p,4e12.2)')
     +      mlv,(qdata(450+mlv+i*6),i=0,3)
        enddo
        write(*,'(24x,''thinlev'',3x,''hadr,muon'',5x,
     +         ''thinlev'',3x,''elec,gamm'')')
      else 
        mwarn = mwarn + 1
        write(*,'(12x,''WARNING: actually no multi-thinning levels'',
     +    '' defined.'')') 
      endif

c--get augerhit core positions (in meter):
      write(*,'(2x)')
      mcores = int(qdata(lenblk+98)) ! copied from pdata vector.
      if ( mcores .gt. 0 ) then
        dradius = qdata(lenblk+174)*1.d-2 ! switch from cm to meter. 
        qradius = dradius * dradius ! [m^2]
        qdata(lenblk+176) = 0. ! current number of scattered core.
        do  lcor=1,mcores
          xshcore = qdata(lenblk+98+lcor) * 1.d-2
          yshcore = qdata(lenblk+118+lcor)* 1.d-2
          xcore(lcor) = xshcore
          ycore(lcor) = yshcore
          rcore(lcor) = sqrt( xshcore*xshcore + yshcore*yshcore )
          write(*,'(14x,3f12.4)') xshcore, yshcore, rcore(lcor)
          write(*,'(40x,''  x='',f8.2,''  y='',f8.2,''  r='',f8.2)')
     +      xcore(lcor), ycore(lcor), rcore(lcor)
        enddo
        write(*,'(2x)')
        write(*,'(12x,''number of core positions ='',i7)') mcores
        write(*,'(12x,''radius around a detector ='',f7.2,'' m'')')
     +    dradius
      else 
        mwarn = mwarn + 1
        write(*,'(12x,
     +          ''WARNING: actually no core positions defined.'',/)')
        if ( mwarn .gt. 1 ) then
          write(*,'(12x,
     +            ''Simply count all types of particles.'',/)')
        endif
      endif

c--work on only one thinning level:
      mlevel = mthins ! fix number of mthins level, loop missing.
      do  mlv=1,mthins 
        qwei8888(mlv) = 0.d0 ! count all lines with 8888 codes.
        qweight1(mlv) = 0.d0 ! count all lines with at least one weight >= 1.
        qindepar(mlv) = 0.d0 ! total number of weights (of each mlv) >= 1.
      enddo

c----------many big particle data files for one shower------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      isho = 0
      iend = 0
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
      if ( ifil .eq. nfil ) write(*,'(20x,a)') cdata(ifil)(1:ilen)
      read(cdata(ifil)(idat+10:idat+15),'(i6)') mprt 
* - - - - - - create names of individual particle data files:
      ! cpllpart = 'DAT000000-000000.mthiasci'
      ! write(cpllpart(4:9),'(i6.6)') irun 
      ! write(cpllpart(11:16),'(i6.6)') mprt
      ! write(cpllpart(21:21),'(i1)') mlevel
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil)(idat:ilen),status='old',
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
      if ( ifil .eq. nfil ) then 
        if ( mthins .ge. 1 ) then
          write(*,'(38('' _''))')
          write(*,'(2x)')
          write(*,'(6x,''number of particles not taking weights '',
     +      ''into account:'')')
          write(*,'(6x,''code'',7('' ------ mlv='',i1,'' ----'',:))')
     +      (mlv,mlv=0,mthins)
          do  l=1,50
            if ( qparsum(l,0) .gt. 0. ) write(*,'(i8,7f18.1)')
     +         l,qpartic(l,0),(qparsum(l,mlv),mlv=1,mthins)
          enddo
          write(*,'(2x)')
          write(*,'(6x,''number of particles taking weights '',
     +      ''into account:'')')
          write(*,'(6x,''code'',7('' ------ mlv='',i1,'' ----'',:))')
     +      (mlv,mlv=0,mthins)
          do  l=1,50
            if ( qparsum(l,0) .gt. 0. ) write(*,'(i8,7f18.1)')
     +         l,qpartic(l,0),(qparmth(l,mlv),mlv=1,mthins) ! +100
          enddo
        endif
      endif
  444 continue
* - - - - - - end-of loop over particle data files.
      close(unit=3)

c--end of current thinning level reached--------------------------------
      if ( mthins .ge. 1 ) then
        write(*,'(38('' _''))')
        write(*,'(2x)')
        write(*,'(34x,''maxweight'',15x,''maxweight'')')
        do  mlv=1,mthins
          write(*,'(12x,''mthi='',i1,'':'',1p,4e12.2)')
     +      mlv,(qdata(450+mlv+i*6),i=0,3)
        enddo
        write(*,'(24x,''thinlev'',3x,''hadr,muon'',5x,
     +         ''thinlev'',3x,''elec,gamm'')')
        write(*,'(2x)')
        write(*,'(9x,''qwei8888 = '',
     +    ''count total number of particles after thinning.'')')
        write(*,'(9x,''qweight1 = '',
     +    ''count all lines with at least one weight >= 1.'')')
        write(*,'(9x,''qindepar = '',
     +    ''total number of weights (of each mlv) >= 1.'')')
        if ( qweight1(1) .gt. 0. ) then
          write(*,'(2x)')
          write(*,'(6x,''average weight per thinned particle'')')
          qsum8888 = 0.d0
          do  mlv=1,mthins
            qsum8888 = qsum8888 + qwei8888(mlv)
            write(*,'(12x,''qindepar / qwei8888 ='',1p,e16.4,0p,
     +        ''   mlv ='',i2)') qindepar(mlv)/qwei8888(mlv),mlv
          enddo 
          write(*,'(2x)')
          write(*,'(6x,''average active weights per'',
     +      '' active codes 8888'')')
          write(*,'(12x,''qsum8888 / qweight1'',1p,e18.4)') 
     +      qsum8888/qweight1(1)
        endif
        write(*,'(2x)')
        write(*,'(6x,''counted sum of unthinned particles'')')
        qprtcls = 0.d0
        do  l=1,50
          qprtcls = qprtcls + qpartic(l,0)
        enddo
        write(*,'(37x,1p,e12.4)') qprtcls
        write(*,'(2x)')
        if ( qsum8888 .lt. 0.01d0 ) qsum8888 = 0.d0 
        write(*,'(6x,''sum of thinned particles qsum8888'')')
        write(*,'(37x,1p,e12.4)') qsum8888 
        write(*,'(2x)')
        write(*,'(6x,''sum of particles with at least one '',
     +    ''weight >= 1'')')
        write(*,'(37x,1p,e12.4)') qweight1(1)
        write(*,'(2x)')
        write(*,'(6x,''ratio of particles with at least one '',
     +    ''weight >= 1 to all particles'')')
        write(*,'(12x,''qweight1 / qprtcls = '',1p,e16.4)')
     +    qweight1(1)/max(1.d0,qprtcls)
        qnplq0sum = 0.d0
        do  l=0,64
           qnplq0sum = qnplq0sum + qnpl1h(l)+qnpl0h(l)
        enddo
        if ( qnplq0sum .gt. 0.d0 ) then
           write(*,'(2x)')
           write(*,'(6x,''particles hit tanks: code'',6x,
     +        ''mthinned'',7x,''only circ'')')
           do  l=0,64
              if ( qnpl1h(l)+qnpl0h(l) .gt. 0. )
     +           write(*,'(27x,i2.2,2f16.1)') l,qnpl1h(l),qnpl0h(l)
           enddo
           write(*,'(27x,''code'',6x,
     +        ''mthinned'',7x,''only circ'')')
        write(*,*) '                        qnplq0sum = ',qnplq0sum
        endif
        write(*,'(2x)')
        write(*,'(38('' _''))')
      endif

c--print total counts of particles--------------------------------------
      write(*,'(/,6x,''code name'',12x,''total particles '',a)')
     +    cdata(nfil)(1:ilen)
      write(*,'(6x,''---- ----'',12x,''--------------- '',16(''-''))')
      do  l=1,46
         if ( qpartic(l,0) .gt. 0. ) write(*,'(i10,a19,f13.0)')
     +        l,qpatext(l),qpartic(l,0)
      enddo
      write(*,'(36x,''counts without weights.'')') 
      if ( lenblk .eq. 312 ) then
         write(*,'(6x,''---- ----'',12x,''--------------- '')')
         do  l=1,46
            if ( qpartic(l,0) .gt. 0. ) write(*,'(i10,a19,f14.1)')
     +         l,qpatext(l),qpartic(l+100,0)
         enddo
         write(*,'(36x,''counts including weights.'')') 
      endif
      write(*,'(6x,''---- ----'',12x,''--------------- '')')

c--print total counts of innermost circle at tank positions-------------
      if ( mcores .ge. 1 ) then
        write(*,'(2x)')
        write(*,'(38('' _''))')
        write(*,'(2x)')
        write(*,'(''  icode  lcor      countall        count3m'',
     +    9x,''(3m/all)'')')
        do  lcor=1,mcores
          qlsumhit = 0
          qlsumh3m = 0
          do  l=1,50
            if ( qparhit(1,lcor,l) .gt. 0. ) then
              qlsumhit = qlsumhit + qparhit(1,lcor,l)
              qlsumh3m = qlsumh3m + qparh3m(1,lcor,l)
              write(*,'(2i6,1p,3e16.4)') l,lcor,qparhit(1,lcor,l),
     +          qparh3m(1,lcor,l),qparh3m(1,lcor,l)/qparhit(1,lcor,l)
            endif
          enddo
          if ( qlsumhit .gt. 0 ) then
            write(*,'(i12,1p,3e16.4)') lcor,qlsumhit,
     +        qlsumh3m,qlsumh3m/qlsumhit
        write(*,'(''  icode  lcor      countall        count3m'',
     +    9x,''(3m/all)'')')
          endif
        enddo
        if ( mthins. eq. 0 ) then
          write(*,'(2x)')
          write(*,'(38('' _''))')
          write(*,'(/,6x,''code name'',12x,''total particles '',a)')
     +      cdata(nfil)(1:ilen)
        write(*,'(6x,''---- ----'',12x,''--------------- '',16(''-''))')
          do  l=1,46
            if ( qpartic(l,0) .gt. 0. ) write(*,'(i10,a19,f13.0)')
     +         l,qpatext(l),qpartic(l,0)
          enddo
          write(*,'(6x,''---- ----------'',6x,''--------------- '')')
          write(*,'(2x)')
        endif
      endif
      if ( mtype8 .gt. 0 ) then
         write(*,'(38('' _''))')
         write(*,'(2x)')
      endif
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      ! - - - - error exit reading pdata vector:
      inull = 0
      do  l=1,lenrec,lenpar
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,lenrec,lenblk
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      write(*,*)
     +'    ____________________________________________________________'
      if ( lenpar .eq. 7 ) then
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      elseif ( lenpar .eq. 8 ) then
      write(*,*) '    ERROR: simulation type of corsika is `thinning`.'
      endif
      write(*,'(9x,''estimated number of particles'',f15.0,
     +   '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull  
      write(*,*) '        subblocks EVTE and RUNE missing!'
      goto 499 
      ! - - - - error exit input parameters:
      write(*,*) '    _________________________________________________'
      write(*,*) '    ERROR: missing some input parameters.'
c--end-of program or print out tabulars---------------------------------
  499 continue
      write(*,'(2x)')
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
      double precision parmas(0:101),phead(40),qparhit(6,20,0:50)
      double precision qpartic(0:200,0:6)
      double precision qparmth(0:200,6),qparsum(0:200,0:6)
      double precision xcore(20),ycore(20),qparh3m(6,20,0:50)
      double precision qnpl1h(0:65),qnpl0h(0:65),rcore(20)
      real pdata(lenthin),bdata(624),edata(624)

      common /atmos/aatm,batm,catm
      common /integ/mtype8,irec,isho,isub,ifil,nfil,mcores
     +  ,mlevel,mthins,mprt,qindepar(6),qweight1(6),qwei8888(6)
     +  ,qnpl1h,qnpl0h
      common /utabl/pdata,parmas,phead,qpartic,bdata,edata,qparh3m
     +  ,qparmth,qparsum,qparhit,xmcoord,ymcoord,dradius,qradius
     +  ,xshcore,yshcore,xdistm,ydistm,hdistm,xcore,ycore,rcore

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
 
      double precision parmas(0:101),phead(40),qparhit(6,20,0:50)
      double precision qpartic(0:200,0:6)
      double precision qparmth(0:200,6),qparsum(0:200,0:6)
      double precision xcore(20),ycore(20),qparh3m(6,20,0:50)
      double precision qnpl1h(0:65),qnpl0h(0:65),rcore(20)
      real pdata(lenthin),bdata(624),edata(624)
      dimension indeprec(6),nbin(0:9)
      common /integ/mtype8,irec,isho,isub,ifil,nfil,mcores
     +  ,mlevel,mthins,mprt,qindepar(6),qweight1(6),qwei8888(6)
     +  ,qnpl1h,qnpl0h
      common /utabl/pdata,parmas,phead,qpartic,bdata,edata,qparh3m
     +  ,qparmth,qparsum,qparhit,xmcoord,ymcoord,dradius,qradius
     +  ,xshcore,yshcore,xdistm,ydistm,hdistm,xcore,ycore,rcore

      lobs = 1
      isub = 0
      lenblk = lenrec / 21
      lenpar = lenblk / 39
      lblthin = lenthin / 21

c-----------keep last two particle info sets (because of MUADDI)--------
      do  i=0,2*lenpar-1
         phead(i+23) = pdata(i+lenrec-2*lenpar+1)
      enddo

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      isub = isub + 1
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------
         lobs = nint(pdata(lia+4))
         if ( isho .le. 0 ) then
            phead(6) = lobs
            do  908  i=1,lobs
               phead(6+i) = 1.e-2 * pdata(lia+4+i)
  908       continue
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
         phead(1) = pdata(lia+1) ! simulated shower number.
         phead(2) = 1.e-5 * pdata(lia+6) ! first interaction (km).
         pdata(lia+4) = real(thickgr(dble(pdata(lia+6))))
         phead(3) = pdata(lia+4)
         phead(4) = pdata(lia+2) ! primary particle code.
         phead(5) = pdata(lia+3) ! primary particle energy in GeV.
         phead(19) = pdata(lia+11) ! phi angle in radian.
         phead(20) = pdata(lia+10) ! theta angle in radian.
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100.d0 ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = real(timev0)
c----------------subblock longi information-----------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         ! contents ignored and/or never within particle data files!
c----------------subblock event end-------------------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then 
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         if ( ifil .eq. nfil ) write(*,'(27x,''irec ='',i8)') irec
         goto 949
      else
c-----------subblock with particle data---------------------------------
      do  917  l=lia,lia+lenblk-1,lenpar
         icode = int(1.d-3*(pdata(l)))
         if ( icode .eq. 8888 ) then
           jplev = mod(int(pdata(l)),1000)
c - - - - - - jplev = 8888100: hit tank, but only standard simulated.
           ! 8888101 <= jplev <= 8888163: hit tank as thinned particle.
           ! 8888001 <= jplev <= 8888063: hit selection circle (2nd 
           !         parameter of keyword AUGSCT) as thinned particle.
           ! jplev = 8888000: should not occur.
           if ( jplev .ge. 100 ) then
             jpdec = mod(jplev,100)
             do  lb=0,9
               nbin(lb) = 0
             enddo
             jdec = jpdec
             lb = -1
  909        continue
             lb = lb + 1
             nbin(lb) = mod(jdec,2)
             jdec = jdec / 2
             if ( jdec .gt. 0 ) goto 909
             qnpl1h(64) = qnpl1h(64) + 1.
             do  i=0,63
                if ( jpdec .eq. i ) qnpl1h(i) = qnpl1h(i) + 1.
             enddo  
           else ! if ( jplev .lt. 100 ) then
             qnpl0h(64) = qnpl0h(64) + 1.
             do  i=0,63
                if ( jplev .eq. i ) qnpl0h(i) = qnpl0h(i) + 1.
             enddo  
           endif  
c - - - - - - icode = 8888: count all particles with any weight:
           jcode = 0
           if ( l .gt. lenpar ) then
             xmcoord = 1.d-2 * pdata(4+l-lenpar)
             ymcoord = 1.d-2 * pdata(5+l-lenpar)
             jcode = int(1.d-3*(pdata(l-lenpar)))
             if ( jcode .eq. 5 .or. jcode .eq. 6 ) then 
               write(8,'(1p,8e14.6)')
     +         (pdata(l-lenpar+i),i=0,lenpar-1)
               write(8,'(1p,8e14.6)')
     +         (pdata(l+i),i=0,lenpar-1)
             endif
           else if ( l .eq. 1 ) then
             xmcoord = 1.d-2 * phead(28)
             ymcoord = 1.d-2 * phead(29)
             jcode = int(1.d-3*(phead(30)))
             if ( jcode .eq. 5 .or. jcode .eq. 6 ) then 
               write(8,'(1p,8e14.6)')
     +         (phead(i),i=30,30+lenpar-1)
               write(8,'(1p,8e14.6)')
     +         (pdata(l+i),i=0,lenpar-1)
             endif
           endif
           if ( jcode .gt. 50 ) then
              if ( jcode .eq. 75 ) then
                 jcode = 45
              else if ( jcode .eq. 76 ) then
                 jcode = 46   
              else if ( jcode .eq. 201 ) then
                 jcode = 41
              else if ( jcode .eq. 301 ) then
                 jcode = 42
              else if ( jcode .eq. 302 ) then
                 jcode = 43
              else if ( jcode .eq. 402 ) then
                 jcode = 44
              else if ( jcode .gt. 200 .and. jcode .lt. 987 ) then
                 jcode = 4
              else
                 if ( jcode .ge. 85 .and. jcode .le. 96 ) then
                    jcode = 47 ! ignore 85,86,95,96.
                 else
                    jcode = 50 
                 endif
              endif
           endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           if ( mthins .gt. 0 ) then
             do  mlv=1,mthins
               indeprec(mlv) = 0 ! only mlv=1 is used of it.
             enddo 
             do  914  mlv=1,mthins ! loop over all levels:
              if ( pdata(l+mlv) .ge. 1. ) then
               qwei8888(mlv) = qwei8888(mlv) + 1. 
              ! count number of particles with weights:
               qparsum(jcode,0) = qparsum(jcode,0) + 1.
               qparsum(jcode,mlv) = qparsum(jcode,mlv) + 1.
               qparmth(jcode,mlv) = qparmth(jcode,mlv) + pdata(l+mlv)
               qpartic(jcode+100,0)=qpartic(jcode+100,0) + pdata(l+mlv)
              ! total number of weights (of all mthins) >= 1:
               qindepar(mlv) = qindepar(mlv) + pdata(l+mlv)
              ! count number of weights >= 1 for all particles:
               indeprec(1) = indeprec(1) + 1
              endif 
  914        continue ! end-of loop over all levels. 
           endif 
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           if ( mcores .ge. 1 ) then
             ! this scope only for a given augerhit core position:
             do  912  lcor=1,mcores 
               xshcore = xcore(lcor)
               yshcore = ycore(lcor) 
               ! test particle coordinates on horizontal stripes:
               yp = ymcoord - yshcore
               iytst = nint(yp / ydistm)
               ytest = yp - ydistm*iytst
               if (.not.(-dradius.le.ytest .and. ytest.le.dradius) )
     +           goto 912
               ! test coordinates on 60 degrees stripes:
               xp = xmcoord - xshcore
               xtest = xp + hdistm * mod(iytst,2)
               ixtst = nint(xtest / xdistm)
               xtest = xtest - xdistm*ixtst
               if (.not.(-dradius.le.xtest .and. xtest.le.dradius) )
     +           goto 912
               pdist = xtest*xtest + ytest*ytest
               if ( pdist .le. qradius ) then
                 ! particle belongs to lcor-th core position:
                 qparhit(1,lcor,jcode) = qparhit(1,lcor,jcode) + 1.
                 if ( pdist .le. 9. ) then
                   qparh3m(1,lcor,jcode) = qparh3m(1,lcor,jcode) + 1.
                 endif 
               endif
  912        continue ! end-of test on augerhit core positions.
           endif
           ! do  mlv=1,mthins ! actually not used.
           if ( indeprec(1) .ge. 1 ) qweight1(1) = qweight1(1) + 1.
           ! enddo
           if ( jcode .gt. 0 .and. jcode .lt. 50 ) then
              qpartic(jcode,0) = qpartic(jcode,0) + 1.
           endif
         else
c - - - - - - other codes for particles instead of 8888 - - - - - - - -
           jcode = int(1.d-3*(pdata(l)))  
           if ( jcode .eq. 75 ) write(8,'(1p,8e14.6)')
     +         (pdata(l+i),i=0,lenpar-1)
           if ( jcode .eq. 76 ) write(8,'(1p,8e14.6)')
     +         (pdata(l+i),i=0,lenpar-1)
           if ( mtype8 .gt. 0 ) goto 917 ! mthins codes already done. 
           xmcoord = 1.d-2 * pdata(4+l)
           ymcoord = 1.d-2 * pdata(5+l)
           if ( jcode .gt. 50 ) then
              if ( jcode .eq. 75 ) then
                 jcode = 45
              else if ( jcode .eq. 76 ) then
                 jcode = 46   
              else if ( jcode .eq. 201 ) then
                 jcode = 41
              else if ( jcode .eq. 301 ) then
                 jcode = 42
              else if ( jcode .eq. 302 ) then
                 jcode = 43
              else if ( jcode .eq. 402 ) then
                 jcode = 44
              else if ( jcode .gt. 200 .and. jcode .lt. 987 ) then
                 jcode = 4
              else
                 if ( jcode .ge. 85 .and. jcode .le. 96 ) then
                    jcode = 47 ! ignore 85,86,95,96.
                 else
                    jcode = 50
                 endif
              endif
           endif
           if ( jcode .gt. 0 .and. jcode .lt. 50 )
     +        qpartic(jcode,0) = qpartic(jcode,0) + 1.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           if ( mcores .ge. 1 ) then
            ! test on augerhit core positions:
            do  916  lcor=1,mcores
             xshcore = xcore(lcor)
             yshcore = ycore(lcor)
c - - - - - - test particle coordinates on horizontal stripes:
             yp = ymcoord - yshcore
             iytst = nint(yp / ydistm)
             ytest = yp - ydistm*iytst
             if (.not.(-dradius.le.ytest .and. ytest.le.dradius) )
     +         goto 916
c - - - - - - test coordinates on 60 degrees stripes:
             xp = xmcoord - xshcore
             xtest = xp + hdistm * mod(iytst,2)
             ixtst = nint(xtest / xdistm)
             xtest = xtest - xdistm*ixtst
             if (.not.(-dradius.le.xtest .and. xtest.le.dradius) )
     +         goto 916
             pdist = xtest*xtest + ytest*ytest
             if ( pdist .le. qradius ) then
               ! particle belongs to lcor-th core position:  
               qparhit(1,lcor,jcode) = qparhit(1,lcor,jcode) + 1.
               if ( pdist .le. 9. ) then
               qparh3m(1,lcor,jcode) = qparh3m(1,lcor,jcode) + 1.
               endif 
             endif
  916       continue
           endif ! end-of mcores >= 1.
         endif ! end-of other codes for particles instead of 8888.
         if ( lenblk .eq. 312 ) then
           qpartic(jcode+100,0) = qpartic(jcode+100,0)+pdata(l+lenpar-1)
         endif
  917 continue
      endif ! end-of subblock with particle data.
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

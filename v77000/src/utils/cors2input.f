c=======================================================================
c
c  c o r s 2 i n p u t . f
c  =======================
c  creates a (possible) corsika steering (i.e. input) file using the
c  existing particle data output file, and also already available for 
c  corsika simulations of pure 64 bit architecture and a particle
c  data file of a parallel simulation (MPI used); data files of the
c  former `simprod` simulation sets are also correctly read;
c  but: some steering quantities may be defined additionally and
c  individually;
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check cors2input.f -o cors2inputfor
c     ifort -C -O0 -check bounds cors2input.f -o cors2inputfor
c execution:
c     ./cors2inputfor < cors2input.name
c     ./cors2inputfor 
c        < then enter particle data file name >
c-----------------------------------------------------------------------
c "standard" Corsika:
c     output format for particle output (blocklength = 22932+8 fixed)
c     each record consists of 21 subblocks of 273 words (5733).
c "thinning" Corsika:
c     output format for particle output (blocklength = 26208+8 fixed)
c     each record consists of 21 subblocks of 312 words (6552).
c-----------------------------------------------------------------------
c     Transfer indices of thinning elements:
c         evth(149) = thin(1); evth(151) = thin(2); evth(152) = thin(3);
c         evth(148) = thin(1)/thinh(1); evth(150) = thin(2)/thinh(2);
c     Reading DAT-file: if (pdata(l) .gt. 9.9e6) bs = mod(pdata(l),1.e5)
c     Reading CER-file: bs = pdata(l) ! bunchsize
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program cors2input

      implicit double precision (a-h,o-z), integer (i-n)

      character chfile*240,chflag(0:1)*7,czeile*240,cinput*240,chuser*7
      character qpatext(0:200)*20,qnuclei(0:29)*2,csklst*240
      character cfmtint*5,crunh*4
      dimension drunh(312),devth(312)
      dimension nevent(16),ncount(16),nclast(16),nbit(0:2)
      double precision aatm(5),batm(5),catm(5),thickgr
      double precision aatmax(5),batmax(5),catmax(5),atmlay(5)
      real pdata(5733),prunh 

      equivalence (crunh,prunh) 
      logical lexist,lstexi
      common /atmos/aatm,batm,catm

      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data aatmax
     + /-129.86987d0, -13.912578d0, 1.1379052d0, -4.5502d-4, 1.12829d-2/
      data batmax
     + /1170.07784d0, 1310.69613d0, 1490.6966d0, 503.613568d0, 1.0d0/
      data catmax
     + /971950.04d0, 682326.92d0, 615751.06d0, 795110.76d0, 1.d-9/
      data atmlay /10.7d5, 14.6d5, 36.6d5, 100.d5, 112.8292d5/
      data chflag/'      F','      T'/,cfmtint/'(i10)'/,crunh/'RUNH'/ 
c - - - - - - number of showers of simprod simulation series:
c     1.0000E+14....1.7783E+14      28117 sh  =  18 * 1500 sh  + 1117 sh
c     1.7783E+14....3.1623E+14      15811 sh  =  15 * 1000 sh  +  811 sh
c     3.1623E+14....5.6234E+14       8891 sh  =  22 *  400 sh  +   91 sh
c     5.6234E+14....1.0000E+15       5000 sh  =  25 *  200 sh
c     1.0000E+15....1.7783E+15       2812 sh  =  28 *  100 sh  +   12 sh
c     1.7783E+15....3.1623E+15       1581 sh  =  26 *   60 sh  +   21 sh
c     3.1623E+15....5.6234E+15        889 sh  =  29 *   30 sh  +   19 sh
c     5.6234E+15....1.0000E+16        500 sh  =  25 *   20 sh
c     1.0000E+16....1.7783E+16        281 sh  =  25 *   11 sh  +    6 sh
c     1.7783E+16....3.1623E+16        158 sh  =  31 *    5 sh  +    3 sh
c     3.1623E+16....5.6234E+16         89 sh  =  44 *    2 sh  +    1 sh
c     5.6234E+16....1.0000E+17         50 sh  =  50 *    1 sh
c     1.0000E+17....1.7783E+17         28 sh  =  28 *    1 sh
c     1.7783E+17....3.1623E+17         16 sh  =  16 *    1 sh
c     3.1623E+17....5.6234E+17          9 sh  =   9 *    1 sh
c     5.6234E+17....1.0000E+18          5 sh  =   5 *    1 sh
      data nevent/
     +   1500, 1000, 400, 200, 100, 60, 30, 20, 11, 5, 2, 1, 1, 1, 1, 1/
      data ncount/
     +    18,   15,  22,  25,  28, 26, 29, 25, 25, 31,44,50,28,16, 9, 5/
      data nclast/
     +   1117,  811,  91,   0, 12, 21, 19,  0,  6,  3, 1, 0, 0, 0, 0, 0/
      data nbit/32,32,64/
      data qnuclei/'  ',
     +             'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     +             'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +             'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu'/
      data (qpatext(ip),ip=0,60)/'                   ',
     c' gamma             ',' positron          ',' electron          ',
     c' (stackin simulat.)',' muon +            ',' muon -            ',
     c' pion 0            ',' pion +            ',' pion -            ',
     c' Kaon 0 long       ',' Kaon +            ',' Kaon -            ',
     c' neutron           ',' proton            ',' anti proton       ',
     c' Kaon 0 short      ','                   ',' Lambda            ',
     c' Sigma +           ',' Sigma 0           ',' Sigma -           ',
     c' Xi 0              ',' Xi -              ',' Omega -           ',
     c' anti neutron      ',' anti Lambda       ',' anti Sigma -      ',
     c' anti Sigma 0      ',' anti Sigma +      ',' anti Xi 0         ',
     c' anti Xi +         ',' anti Omega +      ',16*'                ',
     c'                   ',' omega             ',' rho 0             ',
     c' rho +             ',' rho -             ',' Delta ++          ',
     c' Delta +           ',' Delta 0           ',' Delta -           ',
     c' anti Delta --     ',' anti Delta -      ',' anti Delta 0      '/
      data (qpatext(ip),ip=61,129)/
     c' anti Delta +      ',' Kaon * 0          ',' Kaon * +          ',
     c' Kaon * -          ',' anti Kaon * 0     ',' electron neutrino ',
     c' anti elec neutrino',' muon neutrino     ',' anti muon neutrino',
     c'                   ',' eta=>2*gamma      ',' eta=>3*pi0        ',
     c' eta=>pi+ pi- pi0  ',' eta=>pi+ pi- gamma',40*'                ',
     c'                   ',' D 0               ',' D +               ',
     c' anti D -          ',' anti D 0          ',' D s +             ',
     c' anti D s -        ',' eta c             ',' D * 0             ',
     c' D * +             ',' anti D * -        ',' anti D * 0        ',
     c' D * s +           ',' anti D * s -      ','                   '/
      data (qpatext(ip),ip=130,200)/
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
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - read name of a corsika data file (without quotes):
      lenchf = 9
      iusc = 0
      read(*,'(a)',err=28,end=29) chfile
      i = 240 + 1
      ! get real length of file name:
   11 continue
      i = i - 1
      if ( chfile(i:i) .eq. ' ' ) goto 11
      lenchf = i
      write(*,'(/,12x,a)') chfile(1:lenchf)
      isla = 0
      itpa = index(chfile(1:lenchf),'.part')
      iche = index(chfile(1:lenchf),'CER')
      idat = index(chfile(1:lenchf),'DAT')
      if ( chfile(1:1).eq.'/' ) then
         i = lenchf
   12    continue
         i = i - 1
         if ( chfile(i:i) .ne. '/' ) goto 12
         isla = i ! position of last slash of file name.
      else if ( chfile(1:2).eq.'./' ) then
         isla = 2
      endif
      idat = isla + 1
      iusc = index(chfile(1:lenchf),'_') ! simprod run.
      if ( idat .eq. 0 .and. iusc .eq. 0 ) then
         if ( iche .gt. 0 ) idat = iche
      endif
      chuser = 'prktkm '
      if ( idat .ge. 0 ) then
         cinput = chfile(idat:lenchf)//'.inp'
      else if ( itpa .gt. 3 ) then
         cinput = chfile(idat:lenchf-5)//'.inp'
      else if ( iusc .gt. 9 ) then ! simprod simulation found.
         ! /lxdata/d13lx24/simprod/epos_stat-A/coe16m158_01
         chuser = 'simprod'
         write(cfmtint(3:4),'(i2)') lenchf-iusc
         read(chfile(iusc+1:lenchf),cfmtint) msimu
         if ( chfile(iusc-4:iusc-4).eq.'M' .or.
     +        chfile(iusc-4:iusc-4).eq.'m'
     +     .or. chfile(iusc-7:iusc-7).eq.'E' .or.
     +        chfile(iusc-7:iusc-7).eq.'e' ) then
            cinput = chfile(iusc-9:lenchf)//'.inp'
         else
            cinput = chfile(1:lenchf)//'.inp'
         endif
      endif

c - - - - - - - - binary corsika or cherenkov particle data file:
      open(3,file=chfile(1:lenchf),status='old',form='unformatted')
      read(3,err=27,end=27) (pdata(i),i=1,5733)
      close(3)

c - - - - - check 32 or 64 bit and standard or thinning simulation:
      if ( pdata(1).ge.211284. .and. pdata(1).le.211286. ) then 
         ibit = 0
         if ( pdata(313).ge.217432. .and. pdata(313).le.217434. ) then
            isubr = 312
         else if (pdata(274).ge.217432..and.pdata(274).le.217434.) then
            isubr = 273
         endif
      else if ( abs(pdata(1)).lt.1.e-6 ) then 
         ibit = 1
         if ( pdata(2).ge.211284. .and. pdata(2).le.211286. ) then
           if ( pdata(314).ge.217432. .and. pdata(314).le.217434.) then
               isubr = 312
           else if (pdata(275).ge.217432..and.pdata(275).le.217434.)then
               isubr = 273
           endif
         endif
      else if ( pdata(3).ge.2.0202 .and. pdata(3).lt.9.9999 ) then
         ibit = -1
         if ( pdata(312).ge.217432. .and. pdata(312).le.217434.) then
            isubr = 312
         else if (pdata(273).ge.217432..and.pdata(273).le.217434.) then
            isubr = 273
         endif
      else
         write(*,*) '     pdata(1) =',pdata(1),' case should not occur!'
      endif
      ! copy run header to double array drunh:
      drunh(1) = prunh ! 211285.281
      do  i=2,isubr
         drunh(i) = dble(pdata(ibit+i))
      enddo
      ! exclude invalid default dates, set date to 01.Jan.2011:
      if ( 22. .lt. drunh(3)/10000. .and. drunh(3)/10000. .lt. 99. )
     +   drunh(3) = 110101.
      ! copy event header to double array devth:
      do  i=isubr+1,isubr*2
         devth(i-isubr) = dble(pdata(ibit+i))
      enddo
c - - - - - - check subblocks on multi-thinning identification:
      memadd = 0
      mtype8 = 0
      do  i=1+2*isubr,5*isubr,isubr/39
         icode = int(1.d-3*(pdata(i)*1.0000001d0))
         if ( icode .eq. 8888 ) mtype8 = mtype8 + 1
         if ( icode .eq. 85 .or. icode .eq. 86 .or.
     +        icode .eq. 95 .or. icode .eq. 96 ) memadd = memadd + 1
      enddo
c - - - - - - use quantities of runh end evth subblocks:
      if ( devth(7) .lt. 0. ) devth(7) = -devth(7) ! curved version
      iday = mod(int(drunh(3)),100)
      imon = mod(int(drunh(3)),10000) / 100
      iyer = int(drunh(3)) / 10000 + 1900
      if ( iyer .lt. 1988 ) iyer = iyer + 100
      if ( drunh(4) .gt. 1.e6 ) drunh(4) = drunh(4) * 1.e-6
      write(*,'(/,4x,''runh'',9x,''runnumb'',6x,''date'',9x,'//
     +   '''versprog'',5x,''nobslev'',6x,''obslev1'',6x,''obslev2'','//
     +   '7x,''eslope'')')
      write(*,'(1p,e13.5,0p,i10.6,''.'',2x,1p,5e13.5,0p,f11.4,/)')
     +   drunh(1),int(drunh(2)),(drunh(l),l=3,7),devth(58)
      write(*,'(4x,''engymin'',6x,''engymax'',6x,''flagEGS4'','//
     +   '5x,''flagNKG'',6x,''ecut.hadr'',4x,''ecut.muon'','//
     +   '4x,''ecut.elec'',4x,''ecut.phot'')')
      write(*,'(1p,e14.7,e15.7,0p,f9.4,f12.4,2x,1p,4e13.5,/)')
     +   (drunh(l),l=17,24)
      write(*,'(4x,''evth'',9x,''evtnumb'',6x,''particle'',5x,'//
     +   '''energy'',7x,''altit/gr'',5x,''nrtarget'',5x,'//
     +   '''height/cm'')')
      write(*,'(1p,6e13.5,e16.8,/)') (devth(l),l=1,7)
      write(*,'(4x,''theta/deg'',4x,'//
     +   '''phi/deg'',7x,''hilow'',7x,''nr.seeds'',7x,''seed1'')')
      write(*,'(0p,f13.7,f14.7,f11.4,f12.4,1p,e18.8,/)')
     +(devth(l)*c180pi,l=11,12),devth(155),(devth(l),l=13,14)
      if ( isubr .eq. 312 ) then ! kind of corsika simulation. 
        if ( nbit(1+ibit) .eq. 32 ) then
          write(*,'(19x,''thinning simulation. lsubrec ='',i4)') isubr
        else
          write(*,'(i14,'' bit thinning simulation. lsubrec ='',i4)')
     +      nbit(1+ibit),isubr ! should not occur.
        endif
      else if ( isubr .eq. 273 ) then
        if ( nbit(1+ibit) .eq. 32 ) then
         if ( mtype8 .eq. 0 ) then
          write(*,'(19x,''standard simulation. lsubrec ='',i4)') isubr 
         else if ( mtype8 .gt. 0 ) then
          write(*,'(19x,''multi-thinning simu. lsubrec ='',i4)') isubr
         endif
        else
          write(*,'(i14,'' bit standard simulation. lsubrec ='',i4)') 
     +      nbit(1+ibit),isubr ! should not occur.
        endif
      endif
      write(*,'(/,12x,a,/)') cinput(1:lenchf+4)

c - - - - - test on DATnnnnnn.lst file to get original steering infos:
      csklst = chfile(1:lenchf)//'.lst'
      inquire(file=csklst,exist=lstexi)
      iprm = int(devth(3))
      open(7,file=cinput,form='formatted',
     +     access='sequential',status='unknown')

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( lstexi ) then
         ! use protocol file DATnnnnnn.lst (add accurate angles):
         open(unit=2,file=csklst,form='formatted',
     +        access='sequential',status='old')
   19    continue     
         read(2,'(a)',end=22,err=23) czeile(1:80)
         if ( index(czeile,'RUNNR ').le.0 ) goto 19
         ic = 80
   20    continue
         ic = ic - 1
         if ( czeile(ic:ic) .eq. ' ' ) goto 20
         write(7,'(a)') czeile(1:ic)
         do  iq=2,1234
            read(2,'(a)',end=22,err=23) czeile(1:80)
            if (czeile(1:6) .eq. 'LONGI ' ) then
              gramms = thickgr(dble(devth(48)))
              if ( 2838.0e2 .le. devth(48).and.devth(48) .lt. 2859.9e2 )
     +          gramms = thick_south_oct(dble(devth(48)))
              write(7,'(''FIXHEI '',f16.5,''E2 '',i7,'//
     +          'f15.7,'' g/cm^2'')') devth(7)/100.,int(devth(6)),gramms
              write(7,'(''* ERANGE '',1p,2e16.7)') (devth(4),i=1,2)
              write(7,'(''* ESLOPE '', f16.7)') devth(58)
              write(7,'(''* THETAP '',2f16.7)') (devth(11)*c180pi,i=1,2)
              write(7,'(''* PHIP   '',2f16.7)') (devth(12)*c180pi,i=1,2)
              write(7,'(''* OUTFILE DAT'',i6.6,''.firstint'')')
     +          int(drunh(2))
              write(7,'(''HILOW'',f15.2)') devth(155)
            endif
            ic = 80
   21       continue
            ic = ic - 1
            if ( czeile(ic:ic) .eq. ' ' ) goto 21
            write(7,'(a)') czeile(1:ic)
            if ( index(czeile,'EXIT').ge.1 ) then
               close(2)
               goto 23
            endif
         enddo
   22    continue
   23    continue
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
      ! use quantities of first data record of the particle data file:
      write(7,'(''RUNNR  '',i10,4x,''version'',f9.6,'//
     +   '4x,''rundate'',i3.2,''.'',i2.2,''.'',i4)')
     +   int(drunh(2)),drunh(4),iday,imon,iyer
      if ( drunh(2) .lt. 10000. .and. devth(13) .ge. 6. ) then
        if ( devth(4) .gt. 1.d7 ) then
          icut = int(2.d0 * log10(drunh(17)*1.000001) - 14.d0) 
          drunh(80) = 4.5432d5 * 1.95583d0**icut
          write(7,'(''* PARALLEL'',f10.0,f12.0,''   1   F'')')
     +      drunh(80)/1.d3,drunh(80)
        endif
      endif
      if ( chuser .ne. 'simprod' ) then
         write(7,'(''EVTNR  '',i10)') int(devth(2))
      else
         write(7,'(''EVTNR  '',i10,24x,a)')
     +      int(devth(2)),cinput(1:lenchf+4)
      endif
      is = int(devth(13))
      do  ii=1,is
         write(7,'(''SEED   '',3i10)') (int(devth(13+3*(ii-1)+i)),i=1,3)
      enddo
      if ( iprm .eq. 4 ) then
         write(7,'(''* should not occur: iprm == 4 '')')
      endif
c - - - - - - - primary, energy, angles - - - - - - - - - - - - - - - -
      if ( ( 1 .le. iprm .and. iprm .le. 74 ) .or.
     +   ( 116 .le. iprm .and. iprm .le. 194 ) ) then
         write(7,'(''PRMPAR '',i10,23x,a)') iprm,qpatext(iprm)
      else if ( 201 .le. iprm .and. iprm .le. 5629 ) then
         write(7,'(''PRMPAR '',i10,22x,i4,''-'',a)')
     +      iprm,int(iprm/100),qnuclei(mod(iprm,100))
      else if ( ( mod(iprm,100) .le. 0 .or. 29 .lt. mod(iprm,100) ) .or.
     +   ( drunh(3) .lt. 100. ) ) then
         write(*,*) '  Warning: invalid drunh(3). ',
     +   'Stop! drunh(3)=',drunh(3)
         write(*,*) '  '
         goto 29
      endif
      write(7,'(''ERANGE '',1p,2e13.4)') devth(59),devth(60)
      write(7,'(''ESLOPE '',f13.2)') devth(58)
      write(7,'(''THETAP '',2f14.3)') devth(81),devth(82)
      write(7,'(''PHIP   '',2f14.3)') devth(83),devth(84)
      isho = 1
      iobs = int(devth(47))
      if ( devth(79) .eq. 1. .and. iobs .lt. 10 ) then
         if ( drunh(15) .gt. 0 ) isho = int(drunh(15))  
      else if ( devth(79) .eq. 0. .or. devth(79) .eq. 2. ) then 
         if ( chuser .ne. 'simprod' ) then
            ! check lst file to get number of simulated showers:
            csklst = chfile(1:lenchf)//'.lst'
            inquire(file=csklst,exist=lstexi)
            if ( lstexi ) then 
               open(unit=2,file=csklst,form='formatted',status='old')
               do  iq=1,12345 
                  read(2,'(a)',end=24,err=24) czeile
                  if ( index(czeile,'NUMBER OF GENERATED').ge.1 ) then
                     close(2)
                     read(czeile(30:40),'(i11)') isho
                     goto 24
                  endif
               enddo
            endif
   24       continue
         else ! user 'simprod' special selection of shower numbers.
            is = int(4.*log10(devth(4)) - 19.)
            isho = nevent(is)
            if ( lenchf-iusc .eq. 2 ) then
               if ( msimu .gt. ncount(is) ) isho = nclast(is)
            else 
               isho = 1
            endif
         endif
      endif
      write(7,'(''NSHOW  '',i10)') isho
      iobs = 1
      do  i=1,iobs
         gramms = thickgr(dble(devth(47+i)))
         if ( 2838.0e2 .le. devth(48).and.devth(48) .lt. 2859.9e2 )
     +      gramms = thick_south_oct(dble(devth(48)))
         if ( gramms .lt. 900. ) gramms = gramms + 0.00123456
         write(7,'(''OBSLEV '',f13.2,''E2'',f17.7,'' g/cm^2'')')
     +      devth(47+i)/100.,gramms
      enddo
c - - - - - - - cuts, prints, flags - - - - - - - - - - - - - - - - - -
      write(7,'(''ECUTS     '',2f11.3,1p,2e11.2)') (devth(i),i=61,64)
      write(7,'(''ECTMAP          1.E11'')')
      write(7,'(''RADNKG '',f13.2,''E2'')') devth(147)/100.
      write(7,'(''MAXPRT '',i10)') min(1,int(devth(2)))
      write(7,'(''HADFLG '',6i5)') (int(devth(i)),i=65,70)
c - - - - - - - model quantities - - - - - - - - - - - - - - - - - - - -
      write(7,'(''ELMFLG '',2(3x,a7))')
     +   chflag(int(devth(73))),chflag(int(devth(74)))
      if (devth(76).eq.0.) then
         write(7,'(''* HDPM            T         0'')')
      endif
      if (devth(76).eq.1.) then
         write(7,'(''VENUS           T         0'')')
         if (devth(145).eq.1.) then
            write(7,'(''VENSIG          T'')')
         else if (devth(145).eq.2.) then
            write(7,'(''NEXSIG          T'')')
         else
            write(*,*) ' check elements `evth` 140,142,144 '
         endif
         write(7,'(''VENPAR ''''      ''''         0.'')') ! parch paval
      endif
      if (devth(76).eq.2.) then
         write(7,'(''SIBYLL          T         0'')')
         if (devth(140).ge.1.) then
            write(7,'(''SIBSIG          T'')')
         else
            write(*,*) ' check elements 142,144,145 '
         endif
      endif
      if (devth(76).eq.3.) then
         write(7,'(''QGSJET          T         0'')')
         if (devth(142).ge.1.) then
            write(7,'(''QGSSIG          T'')')
         else
            write(*,*) ' check elements 140,144,145 '
         endif
      endif
      if (devth(76).eq.4.) then
         write(7,'(''DPMJET          T         0'')')
         if (devth(144).ge.1.) then
            write(7,'(''DPJSIG          T'')')
         else
            write(*,*) ' check elements 140,142,145 '
         endif
      endif
      if (devth(76).eq.5.) then
         write(7,'(''NEXUS           T         0'')')
         if (devth(145).eq.1.) then
            write(7,'(''VENSIG          T'')')
         else if (devth(145).eq.2.) then
            write(7,'(''NEXSIG          T'')')
            ! init input files for nexus use:
            write(7,'(''NEXPAR fname inics nexus/nexus.inics'')')
            write(7,'(''NEXPAR fname iniev nexus/nexus.iniev'')')
            write(7,'(''NEXPAR fname initl nexus/nexus.initl'')')
            write(7,'(''NEXPAR fname inirj nexus/nexus.inirj'')')
            ! dummy out files for epos (debug case):
            write(7,'(''NEXPAR fname check none'')')
            write(7,'(''NEXPAR fname histo none'')')
            write(7,'(''NEXPAR fname data  none'')')
            write(7,'(''NEXPAR fname copy  none'')')
         else
            write(*,*) ' check elements 140,142,144 '
         endif
      endif
      if (devth(76).eq.6.) then
         write(7,'(''EPOS            T         0'')')
         ! init input files for epos use:
         write(7,'(''EPOPAR input epos/epos.param      '')')
         write(7,'(''EPOPAR fname inics epos/epos.inics'')')
         write(7,'(''EPOPAR fname iniev epos/epos.iniev'')')
         write(7,'(''EPOPAR fname initl epos/epos.initl'')')
         write(7,'(''EPOPAR fname inirj epos/epos.inirj'')')
         write(7,'(''EPOPAR fname inihy epos/epos.ini1b'')')
         ! dummy out files for epos (debug case):
         write(7,'(''EPOPAR fname check none'')')
         write(7,'(''EPOPAR fname histo none'')')
         write(7,'(''EPOPAR fname data  none'')')
         write(7,'(''EPOPAR fname copy  none'')')
      endif 
c - - - - - - - logicals, hilow, longi, magnet - - - - - - - - - - - - -
      tradius = 0.d0 ! [cm]
      if ( devth(4) .ge. 1.0d+08 ) tradius = 1.d0  
      if ( devth(4) .ge. 3.2d+08 ) tradius = 3.d0
      if ( devth(4) .ge. 1.0d+09 ) tradius = 10.d0  
      if ( devth(4) .ge. 3.2d+09 ) tradius = 30.d0
      if ( devth(4) .ge. 1.0d+10 ) tradius = 100.d0
      if ( devth(4) .ge. 3.2d+10 ) tradius = 120.d0
      if ( devth(4) .ge. 1.0d+11 ) tradius = 140.d0 ! [m]
      write(7,'(''MUMULT          T'')')
      if ( chflag(int(devth(94))) .eq. 'T' )
     +   write(7,'(''MUADDI    '',a7)') chflag(int(devth(94)))
      if ( memadd .gt. 0 ) write(7,'(''EMADDI          T'')')
      if ( tradius .eq. 0.d0 ) then 
         write(7,'(''CORECUT '',f13.3'' cm'')') devth(152)
      else
         write(7,'(''CORECUT '',f13.3'' cm      # thinrad'',
     +      f9.4,''e2'')') devth(152),tradius
      endif
      if ( devth(95) .lt. 1. ) devth(95) = 1.
      write(7,'(''STEPFC '',f14.3)') devth(95)
      if ( devth(93) .ne. 0. ) write(7,'(''ARRANG '',f14.3)') devth(93)
      if ( isho .eq. 1 ) then
         gramms = thickgr(dble(devth(7)))
         if ( 2838.0e2 .le. devth(48).and.devth(48) .lt. 2859.9e2 )
     +      gramms = thick_south_oct(dble(devth(7)))
         write(7,'(''FIXHEI '',f16.5,''E2 '',i7,'//
     +      'f16.7,'' g/cm^2'')') devth(7)/100.,int(devth(6)),gramms
         write(7,'(''* ERANGE '',1p,2e16.7)') (devth(4),i=1,2)
         write(7,'(''* THETAP '',2f16.7)') (devth(11)*c180pi,i=1,2)
         write(7,'(''* PHIP   '',2f16.7)') (devth(12)*c180pi,i=1,2)
         write(7,'(''* OUTFILE DAT'',i6.6,''.firstint'')')
     +      int(drunh(2))
      endif 
c - - - - - - - low energy model and hilow - - - - - - - - - - - - - - -
      if ( devth(75) .ge. 3. ) then
         hilow = 200.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
         write(7,'(''* low energy model fluka was used.'')')
      else if ( devth(75) .eq. 2. ) then
         hilow = 80.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
         write(7,'(''URQMD             T'')')
      else if ( devth(75) .le. 1. ) then
         hilow = 80.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
      endif
      write(7,'(''HILOW'',f15.2)') devth(155)
c - - - - - - - check longitudinal step size - - - - - - - - - - - - - -
      gramlong = 5.
      inquire(file=chfile(1:lenchf)//'.long',exist=lexist)
      if ( lexist ) then
         gramlong = 12.34 
         open(4,file=chfile(1:lenchf)//'.long',form='formatted',
     +        access='sequential',status='old')
         read(4,'(30x,i5,19x,f5.0)',err=25,end=25) longstep,gramlong
   25    continue
         close(4)
      endif
      write(7,'(''LONGI'',11x,''T'',f9.1,2(''      T''))') gramlong
c - - - - - - - check quantities of magnetic field - - - - - - - - - - -   
      if ( devth(72) .gt. -17. .and. devth(72) .lt. -11. .and.
     +     devth(71) .gt.  18. .and. devth(71) .lt.  22. ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''Auger'')')
     +      devth(71),devth(72)
      else if ( devth(72) .gt. 41. .and. devth(72) .lt. 45. .and.
     +          devth(71) .gt. 18. .and. devth(71) .lt. 22. ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''Karlsruhe'')')
     +      devth(71),devth(72)
      else if ( abs(devth(71)) .lt. 1.e-2 .and.
     +          abs(devth(72)) .lt. 1.e-2 ) then
         write(7,'(''MAGNET  '',1p,2e13.3,''    NoMag'')')
     +      devth(71),devth(72)
      else if ( ( 26. .lt. devth(71) .and. devth(71) .lt.  29. ) .and.
     +         ( -46. .lt. devth(72) .and. devth(72) .lt. -50. ) ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''SKA (West-Australia)'')')
     +      devth(71),devth(72)
      else if ( ( 14. .lt. devth(71) .and. devth(71) .lt.  18. ) .and.
     +         ( -55. .lt. devth(72) .and. devth(72) .lt. -51. ) ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''South Pole'')')
     +      devth(71),devth(72)
      else if ( ( 18. .lt. devth(71) .and. devth(71) .lt.  20. ) .and.
     +         (  57. .lt. devth(72) .and. devth(72) .lt.  59. ) ) then
         write(7,'(''MAGNET  '',2f12.4,6x,''Tunka-Rex'')')
     +      devth(71),devth(72)
      else
         write(7,'(''MAGNET  '',2f12.4)') devth(71),devth(72)
      endif
c - - - - - - - thinning quantities - - - - - - - - - - - - - - - - - -
      if (devth(150).gt.0.) then
         devth(149) = 1.00001d0 * devth(149)
         write(7,'(''THIN     '',1p,3e12.4)')
     +      devth(149),devth(151),devth(152)
         write(7,'(''THINH '',2f12.0)') 1.,100.
         !  devth(149)/devth(148),devth(151)/devth(150)
      endif
c - - - - - - - cherenkov quantities - - - - - - - - - - - - - - - - - -
      if (devth(77).gt.0.) then
         if (devth(85).gt.0.) then
            write(7,'(''CERSIZ  '',f10.0)') devth(85)
            write(7,'(''CERFIL    '',i7)') int(devth(92)) 
            ! telescope coord. x,y,z, and radius r (IACT or not):
            write(7,'(''* TELESCOPE  -3.17411e+06  1.50956e+06'',
     +         2x,''21054.8    15.e2   Heat    `IACT`'')')
         endif
         if (devth(96).gt.0.)
     +      write(7,'(''CWAVLG '',2f10.0)') devth(96),devth(97)
         if (devth(86).gt.0. .or. devth(87).gt.0.) then
            write(7,'(''CERARY '',2i5,1p,4e11.2,''  `not IACT`'')')
     +         int(devth(86)),int(devth(87)),(devth(88+i),i=0,3)
         endif
         if (devth(98).gt.1.) 
     +      write(7,'(''CSCAT   '',i9,2x,2(f8.0,''e2''))')
     +         int(devth(98)),drunh(248)/100.,drunh(249)/100. 
      endif
c - - - - - - - core positions - - - - - - - - - - - - - - - - - - - - -
      if (devth(98).gt.0.) then
         mpos = int(devth(98))
         do  is=1,mpos
            write(*,'(''COREPOS   '',1p,2e14.6,8x,''(pos.'',
     +         i2.2,'')'')') devth(98+is),devth(118+is),is 
            write(7,'(''COREPOS   '',1p,2e14.6,8x,''(pos.'',
     +         i2.2,'')'')') devth(98+is),devth(118+is),is 
         enddo
      endif
c - - - - - - - multithin parameters - - - - - - - - - - - - - - - - - -
      if ( devth(177) .gt. 0 ) then
         multhi = int(devth(177))
         write(7,'(''MTHINR  '',f13.3'' cm'')') devth(152)
         do  is=1,multhi
            write(7,'(''MTHINH  '',1p,4(e10.1))') devth(177+is),
     +         devth(183+is),devth(189+is)/devth(177+is),
     +         devth(195+is)/devth(183+is)
            write(7,'(''MSEED   '',3i10)')
     +         int(devth(201+is)),int(devth(207+is)),int(devth(213+is))
         enddo
      endif
c - - - - - - - direct, host, user - - - - - - - - - - - - - - - - - - -
      if ( chuser.eq.'maximo' ) then
         write(7,'(''ATMOD       0'')')
         write(7,'(''ATMA  '',1p,4e15.7,e13.5)') (aatmax(i),i=1,5)
         write(7,'(''ATMB  '',1p,4e15.7,e13.5)') (batmax(i),i=1,5)
         write(7,'(''ATMC  '',1p,4e15.7,e13.5)') (catmax(i),i=1,5)
         write(7,'(''ATMLAY'',5(f9.3,''e5''))') (atmlay(i)*1.d-5,i=1,5)
         write(7,'(''DIRECT /cr/data02/joe/corsika-run/'')')
      endif
      if ( chuser.eq.'simprod' ) then
         write(7,'(''DIRECT ./'')')
      endif
      if ( devth(13) .ge. 6.d0 .and. devth(4) .gt. 1.d7 ) then 
         if ( int(devth(48)*1.00001234d-2) .eq. 2838 ) 
     +      write(7,'(''ATMOD        13'',12x,''South Pole'')')
         write(7,'(''DIRECT csk'',i6.6,''/'')') int(drunh(2))
      else 
        if ( int(devth(48)*1.00001234d-2) .eq. 2838 ) then
          write(7,'(''ATMOD        13'',12x,''South Pole'')')
        else 
          ! write(7,'(''* ATMOSPHERE    20    F     (needs file '',
          !         + ''atmprof20.dat)'')')
          write(7,'(''* ATMOD       0'',12x,''monthly-atmosphere-04'',
     +      4x,''(optional)'')')
          write(7,'(''* ATMA     -1.33894966E+02 -4.70898991E+01'',
     +      ''  8.32783065E-01 -2.10543673E-04'')')
          write(7,'(''* ATMB      1.17413853E+03  1.21573550E+03'',
     +      ''  1.42569186E+03  5.03931697E+02'')')
          write(7,'(''* ATMC      9.67164959E+05  7.68155103E+05'',
     +      ''  6.21852254E+05  7.85601620E+05'')')
          write(7,'(''* ATMLAY    9.90000000E+05  1.21000000E+06'',
     +      ''  3.80000000E+06  1.00000000E+07'')')
        endif
        write(7,'(''DIRECT ./'')')
      endif
      if ( devth(13) .ge. 6.d0 .and. devth(4) .gt. 1.d7 ) then
         write(7,'(''HOST   uc1n996'')')
         write(7,'(''USER   jl5949'')')
      else
         write(7,'(''HOST   iklx288'')')
         write(7,'(''USER   joe'')')
      endif
      write(7,'(''EXIT'')')
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endif
      close(7)
      goto 29 ! only one file will be processed (do not: goto 10). 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   27 continue
      write(*,'(''    err: file name: '',a,''!'')') chfile(1:lenchf)
      goto 29
   28 continue
      write(*,'(''    err: file name: '',a,''!'')') chfile(1:lenchf)
   29 continue
      stop
      end
c=======================================================================

      double precision function heightcm( arg )

c-----------------------------------------------------------------------
c     height (cm) above sea level as function of gramms per cm^2
c-----------------------------------------------------------------------

      implicit double precision (a-h)
      double precision aatm(5),batm(5),catm(5),arg
      common /atmos/aatm,batm,catm

      if ( arg .gt. 631.1d0 ) then
        heightcm = catm(1) * log( batm(1) / (arg - aatm(1)) )
      else if ( arg .gt. 271.7d0 ) then
        heightcm = catm(2) * log( batm(2) / (arg - aatm(2)) )
      else if ( arg .gt. 3.0395d0 ) then
        heightcm = catm(3) * log( batm(3) / (arg - aatm(3)) )
      else if ( arg .gt. 1.28292d-3 ) then
        heightcm = catm(4) * log( batm(4) / (arg - aatm(4)) )
      else
        heightcm = (aatm(5) - arg) / catm(5)
      endif

      return
      end
c=======================================================================

      double precision function height_south_oct(arg)

c-----------------------------------------------------------------------
c  calculates the height [cm] as function of thickness [g/cm**2]
c  for South Pole October atmosphere
c-----------------------------------------------------------------------

      double precision aatm(5), batm(5), catm(5), thickl(5)
      data aatm/-142.801d0,-70.1538d0,1.14855d0,9.10269d-4,1.52236d-3/
      data batm/1177.19d0,1125.11d0,1304.77d0,433.823d0,1.d0/
      data catm/861745.d0,765925.d0,581351.d0,775155.d0,7.4095699d0/
      data thickl/1.03500000d+03, 724.15535998265d0, 522.47650689988d0,
     *            154.56131091416714d0, 2.1803354127318884d-3/
      double precision arg,heigh

      if ( arg .gt. thickl(2) ) then
        heigh = catm(1) * dlog( batm(1) / (arg - aatm(1)) )
      else if ( arg .gt. thickl(3) ) then
        heigh = catm(2) * dlog( batm(2) / (arg - aatm(2)) )
      else if ( arg .gt. thickl(4) ) then
        heigh = catm(3) * dlog( batm(3) / (arg - aatm(3)) )
      else if ( arg .gt. thickl(5) ) then
        heigh = catm(4) * dlog( batm(4) / (arg - aatm(4)) )
      else
        heigh = (aatm(5) - arg) * catm(5)
      endif
      height_south_oct = heigh

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
         thickgr = aatm(1) + batm(1) * exp( -arg / catm(1) )
      else if ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg / catm(2) )
      else if ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg / catm(3) )
      else if ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg / catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif

      return
      end
c=======================================================================

      double precision function thick_south_oct(arg)

c-----------------------------------------------------------------------
c  calculates the thickness [g/cm**2] as function of height [cm]
c  for the south pole atmosphere of October 01.
c-----------------------------------------------------------------------

      double precision aatm(5), batm(5), catm(5), hlay(5)
      data aatm/-142.801d0,-70.1538d0,1.14855d0,9.10269d-4,1.52236d-3/
      data batm/1177.19d0,1125.11d0,1304.77d0,433.823d0,1.d0/
      data catm/861745.d0,765925.d0,581351.d0,775155.d0,7.4095699d0/
      data hlay/ -5.D5, 2.66667D5, 5.33333D5, 8.D5, 9.8765D6/ !SPolJan
      double precision arg,thick

      IF     ( ARG .LT. HLAY(2) ) THEN
        THICK = AATM(1) + BATM(1) * dEXP( (-ARG) / cATM(1) )
      ELSEIF ( ARG .LT. HLAY(3) ) THEN
        THICK = AATM(2) + BATM(2) * dEXP( (-ARG) / cATM(2) )
      ELSEIF ( ARG .LT. HLAY(4) ) THEN
        THICK = AATM(3) + BATM(3) * dEXP( (-ARG) / cATM(3) )
      ELSEIF ( ARG .LT. HLAY(5) ) THEN
        THICK = AATM(4) + BATM(4) * dEXP( (-ARG) / cATM(4) )
      ELSE
        THICK = AATM(5) - ARG / cATM(5)
      ENDIF
      thick_south_oct = thick

      return
      end

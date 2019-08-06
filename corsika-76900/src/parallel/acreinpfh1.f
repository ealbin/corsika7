c=======================================================================
c
c  a c r e i n p f h 1 . f
c  -----------------------
c     Automatic creation of successive steering files to run corsika
c        simulations and of the corresponding shell script files.
c               (qgsjet, gheisha, single shower files)
c-----------------------------------------------------------------------
c     #!/bin/bash
c     # gfortran -O0 -fbounds-check acreinpfh1.f -o acreinpfh1
c     ifort -C -O0 -check bounds acreinpfh1.f -o acreinpfh1
c     ./acreinpfh1 > acreinpfh1.tabout
c     cat acreinpfh1.tabout
c     chmod +x jobfh1*
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program acreinpfh1
  
      implicit double precision (a-h,o-z), integer (i-n)
      character cpinput*15,cplldir*10,cstar0*40,cstar2*13,corsexec*32
      character cjobfh1*13,cquota*8,cquotb*8,crunmpi*26
      dimension mfanz(16),mshif(16)
      dimension ap(16,12),mprim(10),qengy(0:16),theta(10)
      double precision aatm(5),batm(5),catm(5),thickgr
      common /atmos/aatm,batm,catm
      DATA aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      DATA batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      DATA catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data ((ap(i,l),i=1,16),l=1,12)/
     + -129.86987635d0, -13.912578027d0, 1.137905219d0,-4.5501979957d-4,
     + 1170.0778448d0, 1310.6961298d0, 1490.6966654d0, 503.61356785d0,
     + 971950.03798d0, 682326.92144d0, 615751.05849d0, 795110.76d0,
     + 1070000.d0, 1460000.d0, 3660000.d0, 10000000.d0,
     + -140.02703037d0, -32.154019677d0,1.3324218584d0,-3.4759783013d-4,
     + 1177.9948918d0, 1226.57938d0, 1544.3942757d0, 529.86754375d0,
     + 990841.49847d0, 746882.89502d0, 608546.43122d0, 787929.8796d0,
     + 980000.d0, 1430000.d0, 3550000.d0, 10000000.d0,
     + -140.79368856d0, -36.408818366d0,1.4042651261d0,-2.4583883619d-4,
     + 1177.4425334d0, 1189.02147d0, 1556.0070888d0, 539.33390781d0,
     + 984028.35984d0, 766300.21844d0, 604631.06014d0, 782862.57672d0,
     + 900000.d0, 1450000.d0, 3500000.d0, 10000000.d0,
     + -133.89496569d0, -47.089899137d0,0.83278306486d0,-2.105436733d-4,
     + 1174.1385323d0, 1215.7355043d0, 1425.6918596d0, 503.93169671d0,
     + 967164.95864d0, 768155.10252d0, 621852.25389d0, 785601.62044d0,
     + 990000.d0, 1210000.d0, 3800000.d0, 10000000.d0,
     + -152.40981237d0,-12.095554799d0,0.85136995339d0,-4.5928490669d-5,
     + 1190.3622315d0, 1265.7764387d0, 1437.9375049d0, 524.79662686d0,
     + 985619.00774d0, 681017.97599d0, 617989.22267d0, 776008.70875d0,
     + 950000.d0, 1500000.d0, 3750000.d0, 10000000.d0,
     + -187.01662033d0, -25.280960312d0, 0.89432403304d0, 3.47276315d-5,
     + 1227.6177212d0, 1178.8761669d0, 1445.4602411d0, 535.28753109d0,
     + 1012299.9934d0, 732521.03651d0, 614783.70338d0, 771077.62332d0,
     + 750000.d0, 1450000.d0, 3700000.d0, 10000000.d0,
     + -136.34855518d0, -20.914318979d0, 0.87280369409d0,3.732759787d-5,
     + 1182.905474d0, 1213.2533581d0, 1413.9801198d0, 547.58686206d0,
     + 947431.21788d0, 710893.29507d0, 619412.05748d0, 769605.79586d0,
     + 890000.d0, 1400000.d0, 3700000.d0, 10000000.d0,
     + -139.77721439d0,-25.948986634d0,0.70618984827d0,-9.6699826826d-6,
     + 1184.5463072d0, 1205.4324091d0, 1366.057468d0, 542.40761127d0,
     + 953640.44565d0, 723050.32076d0, 629158.2012d0, 772372.2386d0,
     + 880000.d0, 1280000.d0, 3830000.d0, 10000000.d0,
     + -158.54644384d0,-95.110327312d0,0.76408408672d0,-1.5130518952d-4,
     + 1200.3725926d0, 1165.8522743d0, 1373.2445748d0, 521.58648708d0,
     + 986103.78389d0, 881190.80912d0, 628689.50007d0, 781005.59304d0,
     + 690000.d0, 1100000.d0, 3820000.d0, 10000000.d0,
     + -150.1396913d0, -58.140415154d0,0.87315735073d0,-3.6221952628d-4,
     + 1189.9866918d0, 1193.1675932d0, 1422.7442647d0, 490.3837468d0,
     + 990579.96621d0, 800597.35355d0, 623333.46715d0, 793328.53093d0,
     + 900000.d0, 1200000.d0, 3800000.d0, 10000000.d0,
     + -194.83705256d0, -68.148895317d0,0.91697530761d0,-5.928015954d-4,
     + 1231.0108397d0, 1164.4301252d0, 1412.7462127d0, 460.92737873d0,
     + 1042817.5709d0, 832823.86466d0, 624149.1426d0, 805671.77828d0,
     + 700000.d0, 1200000.d0, 3810000.d0, 10000000.d0,
     + -152.92783842d0, -37.51920314d0, 1.5800882678d0,-5.9051575509d-4,
     + 1190.2311741d0, 1189.4576063d0, 1575.714709d0, 484.1030213d0,
     + 997989.25368d0, 766606.85869d0, 599559.45014d0, 802421.45312d0,
     + 870000.d0, 1450000.d0, 3480000.d0, 10000000.d0/

c - - - - - - set number of files, showers, values of energy - - - - - -
      data mprim/ 14, 1206, 5626, 402, 1407, 1608, 2814, 1, 1, 1/
      data mfanz/ 200, 160, 130, 100,  80,  20,  20,  20, 8*20/
      data mshif/   0, 200, 400, 600, 700, 800, 900, 950, 8*990/
      data theta/ 5., 22., 38., 0., 12., 25., 36., 45., 53., 60./
      data qengy/1.e07,1.78e7,3.16e7,5.62e7,1.00e8,1.78e8,3.16e8,5.62e8,
     +      1.00e9,1.78e9,3.16e9,5.62e9,1.00e10,1.78e10,3.16e10,5.62e10,
     +      1.00e11/
      data multiple/2/ ! must be a divisor of isnmax.

c - - - - - - some initialisations for corsika files - - - - - - - - - -
      iatm = 1
      ishow = 1
      obslev = 1452.d2
      obsgrm = thickgr(obslev) 
      isnmax = 12
      corsexec = 'corsika74045_stnd_QGSII4_gheisha'   

c - - - - - - set of 10 showers for each primary, energy, angle, atmos:
      do  29  iset=1,1 ! qgs4 + gheisha.

c - - - - - - some initialisations for corsika files - - - - - - - - - -
      ecut1 = 0.1d0 
      ecut2 = 0.05d0
      ecut3 = 2.5d-4
      ecut4 = 2.5d-4 
      iauge = 1
      ifluk = 0
      eslop = -1.0 ! -2.0
      themin = 14.0
      themax = 14.0
      phideg = 180.d0-153.435d0 ! 26.565
      ectmax = 1.d5 
      ectcut = ectmax*1.d-3

c - - - - - - particle energy loop - - - - - - - - - - - - - - - - - - -
      do  28  iegy=0,0
      engya = qengy(iegy)
      engyb = qengy(iegy)
      write(*,'(''  '')')
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      if ( engya .lt. engyb ) then
         write(*,'(''   energ='',1p,e9.2,'' ...'',e9.2)') 
     +      engya+6.,engyb+6.
      else
         write(*,'(''   energ='',1p,e9.2)') engya+6.
      endif
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      tweight = 1.d-6 * engya ! * 5.d0/8.d0
      tradius = 100.d0 ! 1.d-2 * tweight ! [m]
      thadron = 100.d0
      ! isnmax = mfanz(iegy)
      isnmax = 5
      mstrt0 = 1320 ! + 10*iegy + mshif(iegy)
**    if ( iegy .eq. 3 ) tradius = 40.d0
**    if ( iegy .ge. 5 ) then
**       tradius = 100. ! [m]
**    else
**       thadron = tweight / 100.d0
**       if ( thadron .lt. 1.d0 ) thadron = 1.d0
**       if ( tweight .ge. 100.d0 ) thadron = 100.d0
**    endif

c - - - - - - primary particle loop - - - - - - - -
      do  27  iprm=1,1
      mstrt = mstrt0  ! + 10*ithe + 100*(iprm-1)
      iprim = mprim(iprm)
      write(*,'(10x,''prim='',i6)') mprim(iprm)

c - - - - - - theta angle loop - - - - - - - -
      do  26  ithe=1,1
      themin = 14.d0
      themax = themin
      if ( themax .gt. themin ) then
         write(*,'(16x,''theta='',f7.2,'' ...'',f7.2)') themin,themax
      else
         write(*,'(16x,''theta='',f7.2)') themin
      endif

c - - - - - - (single) shower loop - - - - - - - -
      do  25  isnr=1,isnmax
         msrun = mstrt + isnr - 1
         if ( isnr .eq. 1 .or. isnr .eq. isnmax) 
     +      write(*,'(23x,''msrun='',i8)') msrun
         if ( isnr .eq. 2 .and. isnr .ne. isnmax)
     +      write(*,*) '                      . . . . . . . . .'
         irun = msrun
         im = mod(irun,1000000)

c - - - - - - some more quantities for shell scripts - - - - - - - - - -
         cstar0 = 'corsika75715_stnd_QGSII4_gheisha'
         cstar2 = 'DAT000000.lst'
         write(cstar2(4:9),'(I6.6)') im
         cpinput = 'parallel-000000'
         write(cpinput(10:15),'(I6.6)') im
         cplldir = 'csk000000/' 
         write(cplldir(4:9),'(I6.6)') im
         cjobfh1 = 'jobfh1-000000'
         write(cjobfh1(8:13),'(I6.6)') im
         crunmpi = 'runmpi000000.txt'
         write(crunmpi(7:12),'(I6.6)') im
         crunmpi = cplldir//crunmpi(1:16)

c - - - - - - - - create corsika steering file - - - - - - - - - - - - -
         iseed = msrun*3
         open(unit=9,FILE=cpinput,STATUS='UNKNOWN')
         write(9,'(''RUNNR '',I11)') irun
         if ( index(cstar0,'runner') .gt. 0 ) then
           write(9,'(''PARALLEL'',f9.0,f12.0,''   1  F'')')ectcut,ectmax
         endif
         write(9,'(''NSHOW '',I10)') 1
         write(9,'(''EVTNR '',I10)') 1
         write(9,'(''SEED  '',3I10)') iseed,0,0
         write(9,'(''SEED  '',3I10)') iseed+1,0,0
         write(9,'(''SEED  '',3I10)') iseed+2,0,0
         write(9,'(''SEED  '',3I10)') iseed+3,0,0
         write(9,'(''SEED  '',3I10)') iseed+4,0,0
         write(9,'(''SEED  '',3I10)') iseed+5,0,0
         write(9,'(''PRMPAR'',I10)') iprim
         if ( engya .lt. engyb ) write(9,'(''ESLOPE'',F10.2)') eslop
         write(9,'(''ERANGE'',6x,1P,2E12.4)') engya,engyb
         write(9,'(''THETAP'',2F13.3)') themin,themax
         write(9,'(''PHIP  '',2F13.3)') phideg,phideg ! 0.,360.
       ! write(9,'(''OBSLEV     4000.E2      631.101 g/cm^2'')')
       ! write(9,'(''OBSLEV     3500.E2      673.273 g/cm^2'')')
       ! write(9,'(''OBSLEV     3000.E2      717.622 g/cm^2'')')
       ! write(9,'(''OBSLEV     2500.E2      764.258 g/cm^2'')')
       ! write(9,'(''OBSLEV     2000.E2      813.300 g/cm^2'')')
       ! write(9,'(''OBSLEV     1500.E2      864.871 g/cm^2'')')
       ! write(9,'(''OBSLEV     1000.E2      919.102 g/cm^2'')')
       ! write(9,'(''OBSLEV      500.E2      976.130 g/cm^2'')')
         if ( iauge .eq. 1 ) then
            write(9,'(''OBSLEV   1452.e2        870.000 g/cm^2'')')
            write(9,'(''MAGNET     19.51        -14.18     Auger'')')
         else
            write(9,'(''OBSLEV    110.e2       1022.647 g/cm^2'')')
            write(9,'(''MAGNET     20.40         43.20   Karlsruhe'')')
         endif
         write(9,'(''MAXPRT'',I10)') 1
         write(9,'(''ECTMAP     1.E11'')')
         write(9,'(''ECUTS    '',2F9.4,1p,2e11.2)')
     +      ecut1,ecut2,ecut3,ecut4 ! 0.1,0.05,2.5d-4,2.5d-4
         write(9,'(''RADNKG      5.E5'')')
         write(9,'(''HADFLG'',6I5)') (0,i=1,5),2
         write(9,'(''ELMFLG         T         T'')')
         if ( index(cstar0,'SIB') .gt. 0 ) then
           write(9,'(''SIBYLL         T'',I10)') 0
           write(9,'(''SIBSIG         T'')')
         endif
         if ( index(cstar0,'QGS') .gt. 0 ) then
           write(9,'(''QGSJET         T'',I10)') 0
           write(9,'(''QGSSIG         T'')')
         endif
         if ( index(cstar0,'EPOS') .gt. 0 ) then
           write(9,'(''EPOPAR input epos/epos.param'')')
           write(9,'(''EPOPAR fname inics epos/epos.inics'')')
           write(9,'(''EPOPAR fname iniev epos/epos.iniev'')')
           write(9,'(''EPOPAR fname initl epos/epos.initl'')')
           write(9,'(''EPOPAR fname inirj epos/epos.inirj'')')
           write(9,'(''EPOPAR fname inihy epos/epos.ini1b'')')
           write(9,'(''EPOPAR fname check none'')')
           write(9,'(''EPOPAR fname histo none'')')
           write(9,'(''EPOPAR fname data  none'')')
           write(9,'(''EPOPAR fname copy  none'')')
         endif
         write(9,'(''MUMULT         T'')')
         write(9,'(''MUADDI         T'')')
         write(9,'(''STEPFC        1.'')')
         write(9,'(''ATMOD  0'',14x,''monthly-atmosphere-'',i2.2)') iatm
         write(9,'(''ATMA    '',1p,4e16.8)') (ap(i,iatm),i=1,4)
         write(9,'(''ATMB    '',1p,4e16.8)') (ap(i,iatm),i=5,8)
         write(9,'(''ATMC    '',1p,4e16.8)') (ap(i,iatm),i=9,12)
         write(9,'(''ATMLAY  '',1p,4e16.8)') (ap(i,iatm),i=13,16)
         write(9,'(''LONGI          T      5.     T      T'')')
         if ( index(cstar0,'thin') .gt. 0 ) then
           write(9,'(''THIN      1.e-06'',f11.0,f8.0,''e2'')')
     +        tweight, tradius
           write(9,'(''THINH         1.'',f11.0)') thadron
         endif
       ! write(9,'(''* AUGSCT        20     45.   750.    T     T'')')
         write(9,'(''HILOW         1.e2'')')
         write(9,'(''OUTFILE '',a,''DAT'',i6.6''.firstint'')')
     +      cplldir,irun
         write(9,'(''DIRECT  '',a)') cplldir
         write(9,'(''HOST   forhlr1'')')
         write(9,'(''USER   jl5949'')')
         write(9,'(''EXIT  '')')
         close(unit=9)

c - - - - - - - - create shell script - - - - - - - - - - - - - - - - -
         itimec = 1000
         if ( engya .gt. 3.3333e10 ) itimec = 68 * 60 ! max 4080 min.
         if ( multiple .gt. 1 ) then
      ! - - - - - - - - multiple > 1:
         open(unit=9,FILE=cjobfh1,STATUS='UNKNOWN')
         write(9,'(''#!/bin/bash'')')
**       write(9,'(''# script to run a parallel corsika simulation'','//
**   +      '1x,''on `fh1.scc.kit.edu`'')')
         write(9,'(''#'')')
         write(9,'(''#MSUB -l nodes=2:ppn=20'')')
         write(9,'(''#MSUB -l walltime= 1:29:00'')')
         write(9,'(''#MSUB -l pmem=2000mb'')')
         write(9,'(''#MSUB -q multinode'')')
         write(9,'(''#MSUB -o job_fh1_jobfh1-'',i6.6,''_%j.out'')') im
         write(9,'(''#MSUB -e job_fh1_jobfh1-'',i6.6,''_%j.err'')') im
         write(9,'(''#MSUB -m n'')')
         write(9,'(''#'')')
         write(9,'(''module load mpi 2>&1'')')
         write(9,'(''if [ ! -e qgsdat-II-04 ] ; then'')')
         write(9,'(''   echo "         WARNING:'',
     +      '' no corsika submit path."'')')
         write(9,'(''   exit ;'')')
         write(9,'(''fi'')')
         write(9,'(''echo $PWD'')')
         write(9,'(''cd $PWD'')')
         write(9,'(''rm -rf '',a)') cplldir
         write(9,'(''mkdir '',a)') cplldir
         if ( index(corsexec,'coreas') .gt. 0 ) then
            write(9,'(''cp SIM'',i6.6,''.inp parallel-'',i6.6)') im,im
            write(9,'(''cp SIM'',i6.6,''.sh jobfh1-'',i6.6)') im,im
            write(9,'(''cp SIM'',i6.6,''.* '',a)') im,cplldir
         endif
         write(9,'(''cp '',a,1x,a)') cjobfh1,cplldir
         write(9,'(''cp '',a,1x,a)') cpinput,cplldir
         write(9,'(''cp summ* '',a)') cplldir
         write(9,'(''cp readpart* '',a)') cplldir
         if ( index(corsexec,'mthi') .gt. 0 )
     +      write(9,'(''cp readmthin* '',a)') cplldir
         if ( index(corsexec,'hit') .gt. 0 )
     +      write(9,'(''cp sortaugerhit* '',a)') cplldir
         write(9,'(''cp readcsk2asci* '',a)') cplldir
         write(9,'(''cp totaltimenew* '',a)') cplldir
         write(9,'(''cp postprocessnew* '',a)') cplldir
         write(9,'(''mpirun '',
     +   ''./mpi_corsika75703_stnd_QGSII4_gheisha_runner '',
     +   ''parallel-'',i6.6)') im
         close(unit=9)
         ! write(20,'(''msub '',A)') cjobfh1
         elseif ( multiple .eq. 1 ) then
      ! - - - - - - - - multiple = 1:
         open(unit=9,FILE=cjobfh1,STATUS='UNKNOWN')
         write(9,'(''#!/bin/bash'')')
         write(9,'(''#MSUB -l nodes=1:ppn=1'')')
         write(9,'(''#MSUB -l walltime=09:30:00'')')
         write(9,'(''#MSUB -l pmem=2000mb'')')
         write(9,'(''#MSUB -q singlenode'')')
         write(9,'(''#MSUB -e '',A,''_%j.err'')') cjobfh1
         write(9,'(''#MSUB -o '',A,''_%j.out'')') cjobfh1
         write(9,'(''#MSUB -m n'')')
         write(9,'(''module load mpi 2>&1'')')
         write(9,'(''# cd $WORK'',A)') ! cstar6 
         write(9,'(''./'',A,'' < '',A,'' > '',A)') cstar0,cpinput,cstar2
         write(9,'(''# echo $PWD'')')
         close(unit=9)
         ! write(20,'(''msub '',A)') cjobfh1
         else
      ! - - - - - - - error case multiple <= 0:
         write(*,'(''invalid value `multiple`'',i3)') multiple
         endif
   25 continue
   26 continue
   27 continue
   28 continue
   29 continue
c - - - - - - end of all loops - - - - - - - - - - - - - - - - - - - - -

      stop
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
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg / catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg / catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg / catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif

      return
      end
c=======================================================================

      double precision function heightus(arg)

c-----------------------------------------------------------------------
c  calculates the height [cm] as function of thickness [g/cm**2] 
c  for the u.s. standard atmosphere
c-----------------------------------------------------------------------
      double precision arg,heigh
      double precision aatm(5), batm(5), catm(5), thickl(5)
      data aatm/-186.5562D0,-94.919D0 ,.61289D0, 0.D0 , .01128292D0/
      data batm/1222.6562D0,1144.9069D0,1305.5948D0,540.1778D0,1.D0/
      data catm/ 994186.38D0,878153.55D0,636143.04D0,772170.16D0,1.D9/
      data thickl/1036.1d0,631.1008D0,271.7001D0,
     *                                     3.039500D0,1.2829201D-03/

      IF     ( ARG .GT. THICKL(2) ) THEN
        HEIGH = CATM(1) * dlog( BATM(1) / (ARG - AATM(1)) )
      ELSEIF ( ARG .GT. THICKL(3) ) THEN
        HEIGH = CATM(2) * dlog( BATM(2) / (ARG - AATM(2)) )
      ELSEIF ( ARG .GT. THICKL(4) ) THEN
        HEIGH = CATM(3) * dlog( BATM(3) / (ARG - AATM(3)) )
      ELSEIF ( ARG .GT. THICKL(5) ) THEN
        HEIGH = CATM(4) * dlog( BATM(4) / (ARG - AATM(4)) )
      ELSE
        HEIGH = (AATM(5) - ARG) * CATM(5)
      ENDIF
      heightus = heigh
  
      return
      end

c=======================================================================
c
c  b c r e i n p c o n t . f
c  -------------------------
c     Automatic creation of successive steering files to run corsika
c     simulations and all corresponding shell script files,
c     qgsjetII04 or epos-lhc, fluka, thinning, single shower files,
c     energy intervals 10^17...10^18...10^19...10^20.
c     (monthly atmospheres, M. Roth, A. Bridgeman, ...)
c usage:
c       ./bcr
c-----------------------------------------------------------------------
c     #!/bin/bash
c     # creation script `./bcr`
c     gfortran -O0 -fbounds-check bcreinpcont.f -o bcreinpcont
c     ./bcreinpcont > bcreinpcont.tabout
c     cat bcreinpcont.tabout
c     mv fort.20 sauger500000
c     chmod 750 saug*
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program bcreinpcont
      implicit double precision (a-h,o-z), integer (i-n)
      character cstar1*9,cstar2*48,cstar6*27,cstar7*14,cscript*10
      character cstar0*43,cinpst*10,cinput*9,csterr*37,cstout*38
      dimension ap(16,12),mprim(10),qengy(8),theta(8)
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
      data theta/ 0., 12., 22., 32., 38., 48., 56., 65./
      data mprim/ 14, 5626, 402, 1407, 1206, 2814, 703, 904, 1105, 1/
      data qengy/ 1.000001e7, 1.00000e8, 1.000000e9, 1.000000e10,
     +            1.00000e11,  3.1623e11, 1.000000e11, 3.1623e11/
 
c - - - - - - initialize random generator - - - - - - - - - - - - - - - 
      double precision c(8),u(97,8),uni,cd,cint,cm,twom24,twom48
      integer modcns,ijkl(8),i97(8),j97(8),ntot(8),ntot2(8),jseq
      common /crranma3/cd,cint,cm,twom24,twom48,modcns
      common /crranma4/c,u,ijkl,i97,j97,ntot,ntot2,jseq
      dimension iseed(103,8),rdnr(1)
      iseed(1,1) = 1
      iseed(2,1) = 0
      iseed(3,1) = 0
      call rmmaqd( iseed(1,1),1,'s' )

c - - - - - - some initialisations for corsika files - - - - - - - - - -
      ishow = 1
      ! cstar0 = 'corsika75061_thin_EPOS_fluka_ehist_muprod'
      cstar0 = 'corsika76401_thin_QGSII4_fluka_ehist_muprod'
      cstar6 = '/cr/data02/joe/corsika-run/' ! length=27
      cstar7 = '/home/tmp/joe/' ! length=14
      obslev = 1452.d2
      obsgrm = thickgr(obslev) 
      eslope = -1.d0
      thadron= 100.d0
      isnmax = 3

c - - - - - - set of showers for each model combination,
c - - - - - - energy, primary, angle, atmosphere:
c - - - - - - iset=1: qgsjetII04 + fluka, iset=2: epos-lhc + fluka.
      do  26  iset=1,1

c - - - - - - particle energy loop - - - - - - - - - - - - - - - - - - -
      do  25  iegy=1,1
      engya = qengy(iegy)
      engyb = qengy(iegy+1)
      write(*,'(''  '')')
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      if ( engya .lt. engyb ) then
         write(*,'(''   energ='',1p,e9.2,'' ...'',e9.2)') 
     +      engya+6.,engyb+6.
      else
         write(*,'(''   energ='',1p,e9.2)') engya+6.
      endif

c - - - - - - primary particle loop - - - - - - - -
      do  24  iprm=1,1 
      iprim = mprim(iprm)
      write(*,*) '         prim= ',mprim(iprm)

c - - - - - - angle theta - - - - - - - -
      themin =  0.d0 ! theta(ithe)
      themax = 70.d0 ! themin
      write(*,'(16x,''theta='',f6.2,'' ...'',f6.2)') themin,themax
      mstrt0 = 100000 + 500*(iprm-1) + 5000*(iegy-1)
      do  ir=1,isnmax*iprm
         call rmmard( rdnr(1),1,1 )
      enddo
      iatm = 0

c - - - - - - (single) shower loop - - - - - - - -
      do  22  isnr=1,isnmax
      msrun = mstrt0 + isnr - 1
      call rmmard( rdnr(1),1,1 )
      if ( eslope .ne. -1.d0 ) then
         ELL = engya**(eslope+1.D0)
         EUL = engyb**(eslope+1.D0)
         SLEX = 1.D0 / (eslope+1.D0)
         engyrndm = ( rdnr(1)*EUL + (1.D0-rdnr(1))*ELL )**SLEX
      else
         ELL = engyb / engya
         engyrndm = engya * (ELL**rdnr(1))
      endif
      tweight = 1.d-6 * engyrndm
      tradius = 1.d0 ! [m]
      if ( engyrndm .ge. 1.d+09 ) tradius = 10.d0
      if ( engyrndm .ge. 3.d+09 ) tradius = 30.d0
      if ( engyrndm .ge. 1.d+10 ) tradius = 100.d0 
      if ( engyrndm .ge. 3.d+10 ) tradius = 120.d0 
      if ( engyrndm .ge. 1.d+11 ) tradius = 140.d0 
      if ( mod(isnr,10) .eq. 1 ) iatm = iatm + 1 ! change atm after 10 sh.
      if ( iatm .gt. 12 ) iatm = 1
      if ( isnr .eq. 1 .or. isnr .eq. isnmax) 
     +   write(*,*) '                     msrun=',msrun

c - - - - - - some more quantities for shell scripts - - - - - - - -
      cinput = 'aug000000'
      cscript = 's'//cinput
      cstar1 = 'DAT000000'
      cstar2 = cstar1//'.lst'  
      itext = 0

c - - - - - - - - create corsika steering files - - - - - - - - - - 
      irun = msrun
      im = mod(irun,1000000)
      write(cinput(4:9),'(I6.6)') im
      ised1 = irun*3
      ised2 = ised1 + 1
      ised3 = ised1 + 2
      open(unit=9,FILE=cinput,STATUS='UNKNOWN')
      write(9,'(''RUNNR '',I10)') irun
      write(9,'(''EVTNR '',I10)') 1
      write(9,'(''SEED  '',3I10)') ised1,0,0
      write(9,'(''SEED  '',3I10)') ised2,0,0
      write(9,'(''SEED  '',3I10)') ised3,0,0
      write(9,'(''PRMPAR'',I10)') iprim
      write(9,'(''ERANGE   '',1P,2E15.4)') engyrndm,engyrndm
      write(9,'(''ESLOPE   '',f9.2)') eslope
      write(9,'(''THETAP  '',2F10.2)') themin,themax
      if ( themax .gt. 0.1234 ) then 
         write(9,'(''PHIP    '',2F10.2)') 0.,360.
      else
         write(9,'(''PHIP    '',2F10.2)') 0.,0.
      endif
      write(9,'(''NSHOW '',I10)') ishow
      write(9,'(''OBSLEV'',f12.2,''E2'',f14.3,'' g/cm^2'')')
     +      obslev*1.e-2, obsgrm
      write(9,'(''ECUTS    '',2F7.2,1P,2E9.1)') 0.1,0.01,5.d-5,5.d-5
      write(9,'(''MAXPRT'',I10)') 1
      write(9,'(''ECTMAP     1.E11'')')
      write(9,'(''RADNKG    200.E2'')')
      write(9,'(''HADFLG'',6I5)') (0,i=1,5),2
      write(9,'(''ELMFLG         T         T'')')
      if ( cstar0(19:19) .eq. 'Q' ) then
         write(9,'(''QGSJET         T'',I10)') 0
         write(9,'(''QGSSIG         T'')')
      endif
      if ( cstar0(19:19) .eq. 'E' ) then
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
      ! write(9,'(''EMADDI         T'')')
      write(9,'(''CORECUT '',f10.2)') tradius
      write(9,'(''STEPFC        1.'')')
      write(9,'(''HILOW       100.'')')
      write(9,'(''THIN      1.e-06'',f11.0,f8.0,''e2'')')
     +      tweight, tradius
      write(9,'(''THINH         1.'',f11.0)') thadron
      write(9,'(''ATMOD  0              monthly-atmosphere-'',i2.2)')
     +      iatm
      write(9,'(''ATMA    '',1p,4e16.8)') (ap(i,iatm),i=1,4)
      write(9,'(''ATMB    '',1p,4e16.8)') (ap(i,iatm),i=5,8)
      write(9,'(''ATMC    '',1p,4e16.8)') (ap(i,iatm),i=9,12)
      write(9,'(''ATMLAY  '',1p,4e16.8)') (ap(i,iatm),i=13,16)
      write(9,'(''LONGI          T      5.     T      T'')')
      write(9,'(''MAGNET       19.45     -14.15'')')
      write(9,'(''DIRECT '',A)') '/home/tmp/joe/' ! './'
      write(9,'(''HOST   iklx282'')')
      write(9,'(''USER   joe'')')
      write(9,'(''EXIT  '')')
      close(unit=9)
c - - - - - - - - create shell scripts - - - - - - - - - - - - -
      it = itext + 9
      write(cscript(5:10),'(I6.6)') im
      write(cstar2(it-5:it),'(I6.6)') irun
      open(unit=9,FILE=cscript,STATUS='UNKNOWN')
      write(9,'(''#!/bin/bash'')')
      write(9,'(''#'')')
      write(9,'(''#$ -cwd'')')
      write(9,'(''#$ -j y'')')
      write(9,'(''#$ -o '',A)') cstar6
      write(9,'(''#'')')
      write(9,'(''cd '',A)') cstar6
      write(9,'(''mkdir -p /home/tmp/joe/'')')
      write(9,'(''rm -rf /home/tmp/joe/DAT'',i6.6,
     +   ''* > /dev/null 2>&1'')') irun
      write(9,'(''echo "       hostname =" $HOSTNAME > '',
     +   ''/home/tmp/joe/DAT'',i6.6,''.tmp'')') irun
      write(9,'(''cp /home/tmp/joe/DAT'',i6.6,''.tmp '',
     +   ''/cr/data02/joe/corsika-run/'')') irun
      write(9,'(''./'',A,'' < '',A,'' > '',A)')
     +   cstar0,cinput,cstar2(1:13)
      write(9,'(''#'')')
      write(9,'(''rm -rf /home/tmp/joe/DAT'',i6.6,''.fl*'')') irun
      write(9,'(''ls -l /home/tmp/joe/DAT'',i6.6,
     +   ''* >> /home/tmp/joe/DAT'',i6.6,''.tmp'')') irun,irun
      write(9,'(''date >> /home/tmp/joe/DAT'',i6.6,''.tmp'')') irun
      write(9,'(''mv /home/tmp/joe/DAT'',i6.6,''.long '',a)')
     +   irun,cstar6
      write(9,'(''mv /home/tmp/joe/DAT'',i6.6,1x,a)') irun,cstar6
      write(9,'(''date >> /home/tmp/joe/DAT'',i6.6,''.tmp'')') irun
      write(9,'(''mv /home/tmp/joe/DAT'',i6.6,''.tmp '',a)')
     +   irun,cstar6
      close(unit=9)
      write(20,'(''qsub -l q=augx.q '',A)') cscript 
      write(20,'(''sleep 12'')')
   22 continue
   24 continue
   25 continue
   26 continue
c - - - - - - end of all loops - - - - - - - - - - - - - - - - - - - - -

      stop
      end
c=======================================================================

      subroutine rmmaqd( iseed,iseq,chopt )

c-----------------------------------------------------------------------
c  subroutine for initialization of rmmard
c  arguments:
c   iseed  = seed to initialize a sequence (3 integers)
c   iseq   = # of random sequence
c   chopt  = character to steer initialize options
c            ' '  sequence 1 is initialized with default seed
c            'r'  get status of generator by 3 seeds
c            'rv' complete status of generator is dumped (103 words)
c            's'  set random generator by 3 seeds
c            'sv' set random generator by array with 103 words
c            'v'  vector option set/get status using 103 words
c-----------------------------------------------------------------------

      implicit none
      integer          kseq
      parameter        (kseq = 8)
      common /crranma3/cd,cint,cm,twom24,twom48,modcns
      double precision cd,cint,cm,twom24,twom48
      integer          modcns
      common /crranma4/c,u,ijkl,i97,j97,ntot,ntot2,jseq
      double precision c(kseq),u(97,kseq),uni
      integer          ijkl(kseq),i97(kseq),j97(kseq),
     *                 ntot(kseq),ntot2(kseq),jseq
      double precision cc,s,t,uu(97)
      integer          iseed(3),i,idum,ii,ii97,ij,ij97,iorndm,
     *                 iseq,j,jj,k,kl,l,loop2,m,niter
      character        chopt*(*), cchopt*12
      logical          first
      save
      data             first / .true. /, iorndm / 11 /, jseq / 1 /
c-----------------------------------------------------------------------

      if ( first ) then
        twom24 = 2.d0**(-24)
        twom48 = 2.d0**(-48)
        cd     = 7654321.d0*twom24
        cm     = 16777213.d0*twom24
        cint   = 362436.d0*twom24
        modcns = 1000000000
        first  = .false.
      endif

      cchopt = chopt
      if ( cchopt .eq. ' ' ) then
        iseed(1) = 54217137
        iseed(2) = 0
        iseed(3) = 0
        cchopt   = 's'
        jseq     = 1
      endif

      if     ( index(cchopt,'s') .ne. 0 ) then
        if ( iseq .gt. 0  .and.  iseq .le. kseq ) jseq = iseq
        if ( index(cchopt,'v') .ne. 0 ) then
          read(iorndm,'(3z8)') ijkl(jseq),ntot(jseq),ntot2(jseq)
          read(iorndm,'(2z8,z16)') i97(jseq),j97(jseq),c(jseq)
          read(iorndm,'(24(4z16,/),z16)') u
          ij = ijkl(jseq)/30082
          kl = ijkl(jseq) - 30082 * ij
          i  = mod(ij/177, 177) + 2
          j  = mod(ij, 177)     + 2
          k  = mod(kl/169, 178) + 1
          l  = mod(kl, 169)
          cd =  7654321.d0 * twom24
          cm = 16777213.d0 * twom24
        else
          ijkl(jseq)  = iseed(1)
          ntot(jseq)  = iseed(2)
          ntot2(jseq) = iseed(3)
          ij = ijkl(jseq) / 30082
          kl = ijkl(jseq) - 30082*ij
          i  = mod(ij/177, 177) + 2
          j  = mod(ij, 177)     + 2
          k  = mod(kl/169, 178) + 1
          l  = mod(kl, 169)
          do  ii = 1, 97
            s = 0.d0
            t = 0.5d0
            do  jj = 1, 48
              m = mod(mod(i*j,179)*k, 179)
              i = j
              j = k
              k = m
              l = mod(53*l+1, 169)
              if ( mod(l*m,64) .ge. 32 ) s = s + t
              t = 0.5d0 * t
            enddo
            uu(ii) = s
          enddo
          cc    = cint
          ii97  = 97
          ij97  = 33
          ! complete init by skipping (ntot2*modcns+ntot) randomnumbers
          niter = modcns
          do  loop2 = 1, ntot2(jseq)+1
            if ( loop2 .gt. ntot2(jseq) ) niter = ntot(jseq)
            do  idum = 1, niter
              uni = uu(ii97) - uu(ij97)
              if ( uni .lt. 0.d0 ) uni = uni + 1.d0
              uu(ii97) = uni
              ii97     = ii97 - 1
              if ( ii97 .eq. 0 ) ii97 = 97
              ij97     = ij97 - 1
              if ( ij97 .eq. 0 ) ij97 = 97
              cc       = cc - cd
              if ( cc .lt. 0.d0 ) cc  = cc + cm
            enddo
          enddo
          i97(jseq) = ii97
          j97(jseq) = ij97
          c(jseq)   = cc
          do  jj = 1, 97
            u(jj,jseq) = uu(jj)
          enddo
        endif
      elseif ( index(cchopt,'r') .ne. 0 ) then
        if ( iseq .gt. 0 ) then
          jseq = iseq
        else
          iseq = jseq
        endif
        if ( index(cchopt,'v') .ne. 0 ) then
          write(iorndm,'(3z8)') ijkl(jseq),ntot(jseq),ntot2(jseq)
          write(iorndm,'(2z8,z16)') i97(jseq),j97(jseq),c(jseq)
          write(iorndm,'(24(4z16,/),z16)') u
        else
          iseed(1) = ijkl(jseq)
          iseed(2) = ntot(jseq)
          iseed(3) = ntot2(jseq)
        endif
      endif

      return
      end
c=======================================================================

      subroutine rmmard( rvec,lenv,iseq )

c-----------------------------------------------------------------------
c  r(ando)m (number generator of) mar(saglia type) d(ouble precision)
c  these routines (rmmard,rmmaqd) are modified versions of routines
c  from the cern libraries.
c  arguments:
c   rvec   = double prec. vector field to be filled with random numbers
c   lenv   = length of vector (# of randnumbers to be generated)
c   iseq   = # of random sequence
c-----------------------------------------------------------------------

      implicit none
      integer          kseq
      parameter        (kseq = 8)
      common /crranma3/cd,cint,cm,twom24,twom48,modcns
      double precision cd,cint,cm,twom24,twom48
      integer          modcns
      common /crranma4/c,u,ijkl,i97,j97,ntot,ntot2,jseq
      double precision c(kseq),u(97,kseq),uni
      integer          ijkl(kseq),i97(kseq),j97(kseq),
     *                 ntot(kseq),ntot2(kseq),jseq
      double precision rvec(*)
      integer          iseq,ivec,lenv
      save
c-----------------------------------------------------------------------

      if ( iseq .gt. 0 .and. iseq .le. kseq ) jseq = iseq

      do  ivec = 1, lenv
        uni = u(i97(jseq),jseq) - u(j97(jseq),jseq)
        if ( uni .lt. 0.d0 ) uni = uni + 1.d0
        u(i97(jseq),jseq) = uni
        i97(jseq)  = i97(jseq) - 1
        if ( i97(jseq) .eq. 0 ) i97(jseq) = 97
        j97(jseq)  = j97(jseq) - 1
        if ( j97(jseq) .eq. 0 ) j97(jseq) = 97
        c(jseq)    = c(jseq) - cd
        if ( c(jseq) .lt. 0.d0 ) c(jseq)  = c(jseq) + cm
        uni        = uni - c(jseq)
        if ( uni .lt. 0.d0 ) uni = uni + 1.d0
        if ( uni .eq. 0.d0 ) uni = twom48
        rvec(ivec) = uni
      enddo

      ntot(jseq) = ntot(jseq) + lenv
      if ( ntot(jseq) .ge. modcns )  then
        ntot2(jseq) = ntot2(jseq) + 1
        ntot(jseq)  = ntot(jseq) - modcns
      endif

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

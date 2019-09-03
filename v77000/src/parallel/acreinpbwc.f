c=======================================================================
c
c  a c r e i n p b w c . f
c  -----------------------
c     Automatic creation of successive steering files to run corsika
c        simulations and of the corresponding shell script files.
c               (qgsjet, gheisha, single shower files)
c-----------------------------------------------------------------------
c #!/bin/bash
c #
c # acr.sh:
c # =======
c # compile and run code `acreinpbwc.f` to create examples
c # of parallel corsika steering files and submit scripts;
c #
c   ifort -C -O0 -check bounds acreinpbwc.f -o acreinpbwc
c   ./acreinpbwc > acreinpbwc.tabout
c   cat acreinpbwc.tabout
c   chmod +x jobwcl*
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program acreinpbwc
      implicit double precision (a-h,o-z), integer (i-n)
      character cpinput*15,cplldir*10
      character cjobwcl*13,corsexec*43
      dimension mprim(8),mfanz(16),mshif(16)
      dimension qengy(0:16)

c - - - - - - set number of files, showers, values of energy - - - - - -
      data mprim/ 14, 402, 1206, 2814, 5626, 1, 1, 1/
      data mfanz/ 200, 160, 130, 100,  80,  20,  20,  20, 8*20/
      data mshif/   0, 200, 400, 600, 700, 800, 900, 950, 8*990/
      data qengy/1.e08,1.78e7,3.16e7,5.62e7,1.00e8,1.78e8,3.16e8,5.62e8,
     +      1.00e9,1.78e9,3.16e9,5.62e9,1.00e10,1.78e10,3.16e10,5.62e10,
     +      1.00e11/

c - - - - - - some initialisations for corsika files - - - - - - - - - -
      ecut1 = 0.1d0 
      ecut2 = 0.05d0
      ecut3 = 2.5d-4
      ecut4 = 2.5d-4 
      iauge = 1
      ifluk = 0
      eslop = -2.0
      themin = 20.0
      themax = 20.0
      phideg = 180.d0-153.435d0
      ectmax = 600000. 
      ectcut = ectmax*1.d-3
      corsexec = 'mpi_corsika76400_stnd_QGSII4_gheiatm_runner'
      isnmax = 3

c - - - - - - particle energy loop - - - - - - - - - - - - - - - - - - -
      do  29  iegy=0,0
      engya = qengy(iegy)
      engyb = engya 
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      write(*,'(''   energ='',1p,e9.2,'' ...'',e9.2)')
     +    engya+6.,engyb+6.
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      ! isnmax = mfanz(iegy)
      mstrt0 = 720 ! + 5 * (iegy-1) + mshif(iegy)

c - - - - - - theta angle loop - - - - - - - -
      do  27  ithe=0,0
      themin = 20.d0
      themax = themin

c - - - - - - primary particle loop - - - - - - - -
      do  26  iprm=1,1
      mstrt = mstrt0  ! + 10*ithe + 100*(iprm-1)
      iprim = mprim(iprm)
      write(*,'(10x,''prim='',i6)') mprim(iprm)

c - - - - - - (single) shower loop - - - - - - - -
      do  25  isnr=1,isnmax
         msrun = mstrt + isnr - 1
         if ( isnr .eq. 1 .or. isnr .eq. isnmax)
     +      write(*,*) '                     msrun=',msrun
         if ( isnr .eq. 2 .and. isnr .ne. isnmax)
     +      write(*,*) '                      . . . . . . . . .'
         irun = msrun
         im = mod(irun,1000000)

c - - - - - - some more quantities for shell scripts - - - - - - - -
         cpinput = 'parallel-000000'
         write(cpinput(10:15),'(I6.6)') im
         cplldir = 'csk000000/' 
         write(cplldir(4:9),'(I6.6)') im
         cjobwcl = 'jobwcl-000000'
         write(cjobwcl(8:13),'(I6.6)') im

c - - - - - - - - create corsika steering file - - - - - - - - - - - -
         iseed = msrun*3
         open(unit=9,FILE=cpinput,STATUS='UNKNOWN')
         write(9,'(''RUNNR '',I10)') irun
         write(9,'(''PARALLEL'',f9.0,f12.0,''   1  F'')') ectcut,ectmax
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
         if ( iauge .eq. 1 ) then
            write(9,'(''OBSLEV   1452.e2        870.000 g/cm^2'')')
            write(9,'(''MAGNET     19.46        -14.16     Auger'')')
         else
            write(9,'(''OBSLEV    110.e2       1022.647 g/cm^2'')')
            write(9,'(''MAGNET     20.40         43.20   Karlsruhe'')')
            ! write(9,'(''OBSLEV    370.e2        991.434 g/cm^2'')')
            ! write(9,'(''MAGNET     27.38   -48.34   SKA (West-Australia)')
         endif
       ! write(9,'(''OBSLEV     4000.E2      631.101 g/cm^2'')')
       ! write(9,'(''OBSLEV     3500.E2      673.273 g/cm^2'')')
       ! write(9,'(''OBSLEV     3000.E2      717.622 g/cm^2'')')
       ! write(9,'(''OBSLEV     2500.E2      764.258 g/cm^2'')')
       ! write(9,'(''OBSLEV     2000.E2      813.300 g/cm^2'')')
       ! write(9,'(''OBSLEV     1500.E2      864.871 g/cm^2'')')
       ! write(9,'(''OBSLEV     1000.E2      919.102 g/cm^2'')')
       ! write(9,'(''OBSLEV      500.E2      976.130 g/cm^2'')')
         write(9,'(''MAXPRT'',I10)') 1
         write(9,'(''ECTMAP     1.E11'')')
         if ( index(corsexec,'urqmd') .gt. 0 ) ecut1 = 0.3d0
         write(9,'(''ECUTS    '',2F9.4,1p,2e11.2)')
     +      ecut1,ecut2,ecut3,ecut4
         write(9,'(''RADNKG      5.E5'')') ! 5.E5 in coreas simul.
         write(9,'(''HADFLG'',6I5)') (0,i=1,5),2 ! 0,1,0,1,0,2 ! IceCube
         write(9,'(''ELMFLG         T         T'')')
       ! write(9,'(''SIBYLL         T'',I10)') 0
       ! write(9,'(''SIBSIG         T'')')
       ! write(9,'(''QGSJET         T'',I10)') 0
       ! write(9,'(''QGSSIG         T'')')
         write(9,'(''MUADDI         T'')')
         write(9,'(''MUMULT         T'')')
         write(9,'(''STEPFC        1.'')')
         write(9,'(''LONGI          T     2.0     T     T'')')
         write(9,'(''HILOW       100.00'')')
       ! write(9,'(''* AUGSCT        20     55.   750.    T     T'')')
         write(9,'(''OUTFILE '',a,''DAT'',i6.6,''.firstint'')')
     +      cplldir,im
         write(9,'(''DIRECT '',a)') cplldir
         write(9,'(''HOST   bwunicl'')')
         write(9,'(''USER   jl5949'')')
         write(9,'(''EXIT  '')')
         close(unit=9)

c - - - - - - - - create shell script - - - - - - - - - - - - - - - - -
         itimec = 1000
         if ( engya .gt. 3.3333e10 ) itimec = 68 * 60 ! max 4080 min.
         open(unit=9,FILE=cjobwcl,STATUS='UNKNOWN')
         write(9,'(''#!/bin/bash'')')
         write(9,'(''#'')')
*        write(9,'(''# parallel corsika air simulation '',
*    +      ''on `bwunicluster.scc.kit.edu`:'')')
         write(9,'(''#MSUB -l nodes=2:ppn=28'')')
         write(9,'(''#MSUB -l walltime=00:40:00'')')
         write(9,'(''#MSUB -q multinode'')')
         write(9,'(''#MSUB -l pmem=3000mb'')')
         write(9,'(''#MSUB -o job_uc1_jobwcl-'',i6.6,''_%j.out'')') im
         write(9,'(''#MSUB -e job_uc1_jobwcl-'',i6.6,''_%j.err'')') im
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
         write(9,'(''#'')')
         if ( index(corsexec,'coreas') .gt. 0 ) then
            write(9,'(''cp SIM'',i6.6,''.inp parallel-'',i6.6)') im,im
            write(9,'(''cp SIM'',i6.6,''.sh jobwcl-'',i6.6)') im,im
            write(9,'(''cp SIM'',i6.6,''.* '',a)') im,cplldir
         endif
         write(9,'(''cp '',a,1x,a)') cjobwcl,cplldir 
         write(9,'(''cp '',a,1x,a)') cpinput,cplldir
         write(9,'(''cp summ* '',a)') cplldir
         write(9,'(''cp readpart* '',a)') cplldir
         if ( index(corsexec,'mthi') .gt. 0 )
     +      write(9,'(''cp readmthin* '',a)') cplldir
         if ( index(corsexec,'hit') .gt. 0 )
     +      write(9,'(''cp sortaugerhit* '',a)') cplldir
         write(9,'(''cp readcsk2asci* '',a)') cplldir
         write(9,'(''cp concatcorsika* '',a)') cplldir
         write(9,'(''cp totaltimenew* '',a)') cplldir
         write(9,'(''cp postprocessnew* '',a)') cplldir
         write(9,'(''#'')')
         write(9,'(''mpirun ./'',a,'' parallel-'',i6.6)') corsexec,im
         close(unit=9)
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

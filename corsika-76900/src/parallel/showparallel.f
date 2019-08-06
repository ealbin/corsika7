c=======================================================================
c
c  s h o w p a r a l l e l . f
c  ---------------------------
c     writing a tabular of simulation quantities to get an overview of
c     all available parallel corsika simulation paths named csk002141/
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check showparallel.f -o showparallel
c     ifort -C -O0 -check bounds showparallel.f -o showparallel
c running:
c     # - - - - - create list of Job*.out files:
c     ls -1 csk*/Job*.out > showparallel.jobinfos
c     # - - - - - list paths and Job info files:
c     ./showparallel < showparallel.jobinfos > showparallel.bwc-work
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c Primary   lg(E)  theta    phi   runnr  sizeGBy  procs
c    T(days)  ecutmax  t(min)  files  RATIO  obslev  Xmagn  Zmagn
c        (a)        Xmagn  Zmagn  _corsika_executable_
c        (b)       Antens  gamma  _corsika_executable_
c        < ecutha  ecutmu  ecutel  ecutga  thilev  wmax  lg(thirad) >
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program showparallel

      implicit double precision (a-h,o-z), integer (i-n)

      character(len=240) chjobname(10000), cherrname(10000),
     +  chzeile, cskvers, chtext 
      character(len=24) chtime, chstck
      character(len=39) chpfile
      character(len=35) chdays
      character(len=7) cfmtflt
      character(len=4) cfmtint
      real pdata(1000)
      dimension ecut(4), ecmnt(4), jcexp(4), jobnlen(10000), qmthin(18),
     +          jseedq(6)
      logical lexist
      data cfmtint/'(i2)'/, cfmtflt/'(f 6.2)'/
      data nthinh/0/, nthi/0/, lthi/0/, jseedq/1,2,3,4,5,6/

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read all Job files and keep names in a character array - -
      do  job=1,10000
         read(*,'(a)',err=101,end=102) chjobname(job)
         jl = 240 + 1 
  100    continue
         jl = jl - 1
         if ( chjobname(job)(jl:jl) .eq. ' ' ) goto 100
         jobnlen(job) = jl       
      enddo
      goto 102
  101 continue
      write(*,'(8x,''e r r o r  reading chjobname'',a)') chjobname(job)
  102 continue
      job = job - 1
      mcodprev = 0
      engyprev = 0.
      thetprev = 0.
      phiaprev = 0.
      obslprev = 0.
      lthiprev = 0
      lenvers = 26
      lversion = 1 ! printing of current parallel executable. 
      mrun = 0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read and check content of each file Job*.out - - - - - -
      do  150  jj=1,job

      read(chjobname(jj)(4:9),'(i6)') mrun
      energy = 0.001
      ectmin = 2.e2
      ectmax = 2.e5
      sumgiga = 0.001
      do  j=1,4
         ecut(j) = 0.
      enddo   
      deratio = 0.
      thinlev = 0.
      thiwmax = 0.
      thirad  = 0.
      obslev = 0. 
      bxmag = 1.e-5
      bzmag = 1.e-5
      jbxmag = 0  ! variant for coreas.
      gbzmag = 0. ! variant for coreas.
      gbymag = 0. ! variant for coreas.
      theta = 0.
      phia = 0.
      phib = 0. 
      jprocs = 0 
      jseeds = 0
      mfiles = 0 
      mthinh = 0
      mseeds = 0
      mprocs = 0 
      mcode = 0
      lthi = 0 
      inquire(file=chjobname(jj),exist=lexist)
      if ( .not. lexist ) goto 149 

c - - - - get current number of subshowers from Job*.err:
      jbl1 = index(chjobname(jj),'.out')
      cherrname(jj) = chjobname(jj)(1:jbl1)//'err'
      open(unit=1,file=cherrname(jj)(1:jobnlen(jj)),
     +     form='formatted',access='sequential',status='old')
      ifiles = 99000
      jfiles = 99000
      read(1,'(2i7)',end=103,err=103) ifiles,jfiles
      if ( mrun .ge. 11000 .or. mrun .lt. 100 ) then ! found coreas:
         jbxmag = 0
         gbzmag = 0. ! = 2.0
         gbymag = 0. ! = 1.e40
         read(1,'(i7,2f7.1)',end=103,err=103) jbxmag,gbzmag,gbymag
       ! read(1,*,end=103,err=103) jbxmag,gbzmag,gbymag
       ! read(1,'(i7,f7.1,e8.1)',end=103,err=103) jbxmag,gbzmag,gbymag
      endif
  103 continue
      close(unit=1)
      gbxmag = jbxmag ! found coreas simulation.
      tminuts = 0.987654
      if ( jfiles .eq. 0 ) jfiles = ifiles

c - - - - implicit loop on lines in file Job*.out:
      open(unit=2,file=chjobname(jj)(1:jobnlen(jj)),
     +            form='formatted',status='old')
      lin = 0
  104 continue
      lin = lin + 1
      read(2,'(a)',err=144,end=144) chzeile 

      ! - - - comment lines in steering file - - - - - - - - - - - - - -

      if ( chzeile(1:1) .eq. '*' .or. chzeile(1:2) .eq. 'c ' 
     +   .or. chzeile(1:2) .eq. 'C ' ) then
            jeng = index( chzeile, 'ERANGE')
            if ( jeng .ge. 2 ) then
               ! - - - - original primary energy (stck_in simulation)
               jbl1 = index( chzeile(jeng:240), ' ') + jeng - 1
  105          continue
               jbl1 = jbl1 + 1
               if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 105
               jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
               write(cfmtflt(3:4),'(i2)') jbl2-jbl1
               read(chzeile(jbl1:jbl2-1),cfmtflt) energy
               deratio = 1.
               if ( ectmax .gt. 0. ) deratio = energy/ectmax
               energy = 9.+log10(energy)
            endif
            goto 104
      endif
      if ( chzeile(1:6) .eq. '      ' ) goto 104

      ! - - - version of corsika mpi runner executable - - - - - - - - -
         if ( index(chzeile,'corsika-initial') .gt. 0 ) then
            jbl1 = index( chzeile(1:240), 'parallel') + 23
            jbl2 = index( chzeile(1:240), '_pll') + 3
            lenvers = jbl2-jbl1+1
            cskvers(1:jbl2-jbl1+1) = chzeile(jbl1:jbl2)
         endif
         if ( chzeile(1:10) .eq. 'job_submit' ) then ! old hc3 syntax.
            jbl1 = index( chzeile(1:240), 'mpi_') + 11
            jbl2 = index( chzeile(1:240), 'parallel') - 8
            lenvers = jbl2-jbl1+1
            cskvers(1:jbl2-jbl1+1) = chzeile(jbl1:jbl2)  
         endif
         if ( chzeile(1:6) .eq. 'mpirun' ) then
            jbl1 = index( chzeile(1:240), 'mpi_') + 11
            jbl2 = index( chzeile(1:240), '_runner')
            lenvers = jbl2-jbl1+1
            cskvers(1:jbl2-jbl1+1) = chzeile(jbl1:jbl2)
         endif 

      ! - - - number of processors at stuttgart simulation - - - - - - -
         if ( chzeile(1:16) .eq. '#PBS -l mppwidth' ) then
            jbl1 = 18
            jbl2 = jbl1
  106       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .ne. ' ' ) goto 106
            jbl2 = jbl2 - 1 
            write(cfmtint(3:3),'(i1)') jbl2 - jbl1 + 1
            read(chzeile(jbl1:jbl2),cfmtint) mprocs
         endif

      ! - - - name of executable at stuttgart simulation - - - - - - - -
         if ( index(chzeile,'aprun') .gt. 0 ) then
            jbl1 = index(chzeile,'mpi_') + 4 
            jbl2 = index(chzeile,'_runner')
            lenvers = jbl2 - jbl1 + 1
            cskvers(1:lenvers) = chzeile(jbl1:jbl2)
         endif

      ! - - - Files: - - - - - - - - - - - - - - - - - - - ik3 cluster -
         if ( index( chzeile, 'Real-time') .gt. 0 ) then
            read(chzeile(1:10),'(i10)') jtimbeg
            lin = lin + 1
            read(2,*,err=144,end=144) jtimend
            lin = lin + 1
            read(2,*,err=144,end=144) jfiles
            lin = lin + 1
            read(2,*,err=144,end=144) jsecnds
            mfiles = jfiles
            mprocs = jfiles
         endif

      ! - - - Exit and read sum of bytes of DAT files - - - - - - - - -
         if ( index( chzeile, 'EXIT') .ge. 1 ) then
            read(2,*,err=144,end=144) sumgiga
            goto 144
         endif

      ! - - - Command: - - - - - - - - - - - - - - - - - - - - - - - - - 
         if ( index( chzeile, 'Command:') .gt. 0 ) then
            if ( index( chzeile, ' -m') .le. 0 ) goto 149 
            ! - - - - number of processors - - - - - - - - - - - - - - -
            jpcs = index( chzeile, ' -p') + 3 
            if ( jpcs .gt. 9 ) then
               do  j=0,240-jpcs
                  if ( chzeile(jpcs+j:jpcs+j) .ne. ' ' ) goto 107
               enddo
  107          continue  
               jpcs = jpcs + j
               jblk = index( chzeile(jpcs:240), ' ')
               write(cfmtint(3:3),'(i1)') jblk-1
               read(chzeile(jpcs:jpcs+jblk-2),cfmtint) mprocs
            endif
         endif

      ! - - - number of processors of the old hc3 simulation - - - - - -
         if ( chzeile(1:10) .eq. 'job_submit' ) then
            jbl1 = index( chzeile(1:240), ' -p') + 3
            jbl2 = index( chzeile(jbl1:240), ' -') + jbl1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1+1
            read(chzeile(jbl1:jbl2-1),cfmtint) mprocs
         endif

      ! - - - number of processors of the bwcluster-simulation - - - - -
         if ( chzeile(1:14) .eq. '#MSUB -l nodes' ) then
            jbl1 = index( chzeile(1:240), 's=') + 2
            jbl2 = index( chzeile(1:240), ':') - 1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1+1
            read(chzeile(jbl1:jbl2),cfmtint) jnodes
            jbl3 = index( chzeile(1:240), 'n=') + 2
            jbl4 = jbl3 + 1
            if ( chzeile(jbl4:jbl4) .eq. ' ' ) jbl4 = jbl3
            write(cfmtint(3:3),'(i1)') jbl4-jbl3+1
            read(chzeile(jbl3:jbl4),cfmtint) jprocs
            mprocs = jnodes * jprocs
         endif

      ! - - - parallel parameters - - - - - - - - - - - - - - - - - - -
         jpar = index( chzeile, 'PARALLEL')
         if ( jpar .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  108       continue  
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 108 
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1  
            read(chzeile(jbl1:jbl2-1),cfmtflt) ectmin
  109       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 109
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 -1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ectmax
         endif 

      ! - - - read and count seeds for random numbers - - - - - - - - -
         jsed = index( chzeile, 'SEED ')
         if ( jsed .eq. 1 ) then
            jseeds = jseeds + 1
            jbl1 = index( chzeile(jsed:240), ' ')
  110       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 110
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtint) jseedq(jseeds)
         endif

      ! - - - primary particle code - - - - - - - - - - - - - - - - - - 
         jcod = index( chzeile, 'PRMPAR')
         if ( jcod .eq. 1 ) then
            jbl1 = index( chzeile(jcod:240), ' ') ! + jcod
  111       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 111
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtint) mcode
         endif

      ! - - - primary particle energy - - - - - - - - - - - - - - - - - 
         jeng = index( chzeile, 'ERANGE')
         if ( jeng .ge. 1 ) then
            jbl1 = index( chzeile(jeng:240), ' ')
  112       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 112
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) energy
            deratio = 1.
            if ( ectmax .gt. 0. ) deratio = energy/ectmax
            energy = 9.+log10(energy)
         endif

      ! - - - theta angle in degrees - - - - - - - - - - - - - - - - - -
         jthe = index( chzeile, 'THETAP')
         if ( jthe .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  113       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 113
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) theta
  114       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 114
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) thetb
            if ( theta .lt. thetb ) then
               ! get theta angle from first particle data file:
               chpfile = chjobname(jj)(1:19)//'-000001'
               write(chpfile(11:13),'(a3)') 'DAT'
               inquire(file=chpfile,exist=lexist)
               if (.not.lexist) then
                  chpfile = chjobname(jj)(1:19)//'-000000000-000000001'
                  write(chpfile(11:13),'(a3)') 'DAT'
                  write(chpfile(21:29),'(i9.9)') jseedq(6) 
               endif
               open(unit=4,file=chpfile,
     +              form='unformatted',status='old')
               read(4,err=115,end=115) (pdata(j),j=1,936)
  115          continue
               close(unit=4)
               lblock = 273
               if ( 217433.0 .lt. pdata(312+1) .and.
     +            pdata(312+1) .lt. 217433.2 ) lblock = 312
               theta = 57.2957795 * pdata(11+lblock)
            endif
         endif

      ! - - - azimuth angle in degrees - - - - - - - - - - - - - - - - -
         jphi = index( chzeile, 'PHIP')
         if ( jphi .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  116       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 116
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) phia
  117       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 117
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) phib
            if ( phia .lt. phib ) then
               ! get phi angle from first particle data file:
               chpfile = chjobname(jj)(1:19)//'-000001'
               write(chpfile(11:13),'(a3)') 'DAT'
               inquire(file=chpfile,exist=lexist)
               if (.not.lexist) then
                  chpfile = chjobname(jj)(1:19)//'-000000000-000000001'
                  write(chpfile(11:13),'(a3)') 'DAT'
                  write(chpfile(21:29),'(i9.9)') jseedq(6)
               endif
               open(unit=4,file=chpfile,
     +              form='unformatted',status='old')
               read(4,err=118,end=118) (pdata(j),j=1,936)
  118          continue
               close(unit=4)
               lblock = 273
               if ( 217433.0 .lt. pdata(312+1) .and.
     +            pdata(312+1) .lt. 217433.2 ) lblock = 312 
               phia = 57.2957795 * pdata(12+lblock) 
            endif
         endif

      ! - - - observation level - - - - - - - - - - - - - - - - - - - -
         jlev = index( chzeile, 'OBSLEV')
         if ( jlev .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  119       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 119
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) obslev
            obslev = obslev / 100.
         endif

      ! - - - magnetic field components - - - - - - - - - - - - - - - -
         jmag = index( chzeile, 'MAGNET')
         if ( jmag .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  120       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 120
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) bxmag
  121       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 121
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) bzmag
         endif 

      ! - - - energy cuts of four particle groups - - - - - - - - - - -
         ject = index( chzeile, 'ECUTS ')
         if ( ject .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  122       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 122
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) ecut(1)
  123       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 123
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ecut(2)
  124       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 124
            jbl4 = index( chzeile(jbl3:240), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) ecut(3)
  125       continue
            jbl4 = jbl4 + 1
            if ( chzeile(jbl4:jbl4) .eq. ' ' ) goto 125
            jbl5 = index( chzeile(jbl4:240), ' ') + jbl4 - 1
            write(cfmtflt(3:4),'(i2)') jbl5-jbl4
            read(chzeile(jbl4:jbl5-1),cfmtflt) ecut(4)
            do  j=1,4
               if ( ecut(j) .lt. 1. ) then
                  jcexp(j) = int(log10(ecut(j))*1.000003-0.999)
               else
                  jcexp(j) = int(log10(ecut(j))*1.000003)    
               endif  
               ecmnt(j) = ecut(j)/10.d0**jcexp(j)
            enddo  
         endif 

      ! - - - corecut/mthinr specification - - - - - - - - - - - - - - -
         corecut = 0.d0
         jcor = index( chzeile, 'CORECUT')
         if ( jcor .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  126       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 126
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) corecut
         endif
         jmtr = index( chzeile, 'MTHINR')
         if ( jcor .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  127       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 127
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) corecut
         endif

      ! - - - multi thinning specifications- - - - - - - - - - - - - - -
         jthm = index( chzeile, 'MTHINH')
         if ( jthm .ge. 1 ) then
            mthinh = mthinh + 1
            if ( mthinh .eq. 1 ) nthinh = nthinh + 1
            jbl1 = index( chzeile, ' ')
  128       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 128
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) qthina
            qmthin(1+3*(mthinh-1)) = qthina
  129       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 129
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) qthinb
            qmthin(2+3*(mthinh-1)) = qthinb
  130       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 130
            jbl4 = index( chzeile(jbl3:240), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) qthinc
            qmthin(3*mthinh) = qthinc
         endif
         jsed = index( chzeile, 'MSEED')
         if ( jsed .ge. 1 ) mseeds = mseeds + 1

      ! - - - thinning specification - - - - - - - - - - - - - - - - - -
         jthi = index( chzeile, 'THIN ')
         if ( jthi .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  140       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 140
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) thinlev
            thinlev = log10(thinlev+4.67735141e-34)
  141       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 141
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) thiwmax
  142       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 142
            jbl4 = index( chzeile(jbl3:240), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) thirad
            thirad = log10(thirad+4.67735141e-34) - 2. ! meters
            nthi = nthi + 1 
            lthi = 1 
         endif       

      if ( lin .lt. 100 ) goto 104

c - - - end-of implicit loop on lines - - - - - - - - - - - - - - - - -
  144 continue 
      close(unit=2)

      ! - - - if energy still 0.0: check DAT00nnnn.stck file - - - - - - 
      if ( energy .lt. 1.d0 ) then
         chstck = 'DAT000000.stck'
         write(chstck(4:9),'(i6.6)') mrun
         inquire(file=chstck,exist=lexist)
         if ( .not. lexist ) then 
            ! - - - - - - check stck file in subdirectory :
            chstck = 'csk000000/DAT000000.stck'
            write(chstck(4:9),'(i6.6)') mrun
            write(chstck(14:19),'(i6.6)') mrun
            inquire(file=chstck,exist=lexist)
         endif
         if ( lexist ) then
            open(unit=2,file=chstck,form='formatted',status='old')
            read(2,*) jstck, energy, jprim, h1stck 
            close(unit=2)
            deratio = 1.
            if ( ectmax .gt. 0. ) deratio = energy/ectmax
            energy = 9.+log10(energy)          
         endif 
      endif

c - - - - - - time analysis from file time.txt - - - - - - - hc3 - - - -
      if ( tminuts .lt. 1. ) then
      iest = 0
      chtime = chjobname(jj)(1:10)//'time.txt'
  145 continue
      inquire(file=chtime,exist=lexist)
      if ( lexist ) then
         open(unit=3,file=chtime,form='formatted',status='old')
         read(3,'(a)',err=148,end=148) chtext
         read(3,*,err=148,end=148) secst,secfi,tminuts
         read(3,'(a)',err=148,end=148) chtext
         read(3,'(a)',err=148,end=148) chtext
         read(3,'(a)',err=148,end=148) chdays
         read(3,'(a)',err=148,end=148) chdays
         read(chdays(24:35),'(f12.6)') ectmin ! totaldays 
         jbl1 = index(chtext,'jobs =')
         jbl2 = jbl1+6
         ! count blanks before digits of mfiles:
  146    continue
         jbl2 = jbl2+1
         if ( chtext(jbl2:jbl2) .eq. ' ' ) goto 146       
         ! count valid digits of mfiles:
  147    continue
         jbl2 = jbl2+1            
         if ( chtext(jbl2:jbl2) .ne. ' ' ) goto 147
         write(cfmtint(3:3),'(i1)') jbl2-(jbl1+6)
         read(chtext(jbl1+6:jbl2-1),cfmtint) mfiles
  148    continue
         close(unit=3)
      else
         iest = iest + 1
         if ( iest .eq. 2 ) 
     +      chtime = chjobname(jj)(1:10)//'time.txt'//chjobname(jj)(4:9)
         if ( iest .eq. 3 ) chtime 
     +      = chjobname(jj)(1:10)//'time'//chjobname(jj)(4:9)//'.txt'
         if ( iest .le. 3 ) goto 145
      endif
      else
         ! tminuts already available from Job0*_scc.out (see above).
      endif

c - - - - - - - print tabular - - - - - - - - - - - - - - - - - - - - -
      if ( jj .eq. 1 ) then  
      ! - - - - - - print first title line - - - - - - - - - - - - - - -
      write(*,'(54x,''Total'',13x,''Wall'',9x,''E/ecutmax'')')
      if ( mrun .ge. 7000 .or. mrun .lt. 100 ) then
        write(*,'(''primary  lg(E)  theta   phi   runnr  sizeGBy'',
     +      2x,''procs  T(days)  ecutmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',2x,''Antens'',3x,''gamma'',
     +      2x,''_corsika_executable_'')')  
      else
        write(*,'(''primary  lg(E)  theta   phi   runnr  sizeGBy'',
     +      2x,''procs  T(days)  ecutmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''_corsika_executable_'')')  
      endif
      endif
      ! if ( 599 .le. mrun .and. mrun .le. 599 ) goto 149 ! no print.
      ! - - - - - - mcode=0 for stackin simulations - - - - - - - - - -
      if ( mcode .eq. 0 ) mcode = 4
      if ( mfiles .eq. 1 ) mprocs = 1
      ! - - - - - - distinguish `thin` and `stnd` - - - - - - - - - - -
      if ( sumgiga .lt. 1.e-2 ) sumgiga = 1.e-2
      if ( mcode .ne. mcodprev .or. energy .ne. engyprev .or.
     +     theta .ne. thetprev .or. phia   .ne. phiaprev .or.
     +    obslev .ne. obslprev .or. lthi   .ne. lthiprev ) then
        write(*,'(1x)')
      endif
      ! - - - - - - print infos for each simulation - - - - - - - - - -
      if ( jbxmag .gt. 0 ) then ! copy special coreas quantities:
         bxmag = gbxmag
         bzmag = gbzmag 
      endif
      if ( mfiles .ge. 100000 ) mfiles = 99999
      if ( tminuts .ge. 10000. ) tminuts = 9999.901234
      if ( deratio .ge. 1000000.d0 ) deratio = 999999.901234
      ! if ( ifiles .lt. mfiles ) then ! almost all DAT files deleted:
      if ( ifiles .lt. mfiles-1 .or. jfiles .lt. jprocs-1 ) then
       if ( lversion .eq. 0 ) then
        if ( lthi .ge. 1 ) then
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2), 
     +      ss,f6.0,1p,e9.1,0p,f6.1)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4),thinlev,thiwmax,thirad
        else
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2))')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4)
        endif
       else ! ( lversion .ne. 0 )
         if ( bxmag .ge. 100. ) then
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,f8.0,f8.2,''  _'',a)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      cskvers(1:lenvers)
         else
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,''  _'',a)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      cskvers(1:lenvers)
         endif
       endif 
      else ! all particle data files available:
       if ( lversion .eq. 0 ) then
        if ( lthi .ge. 1 ) then
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2), 
     +      ss,f6.0,1p,e9.1,0p,f6.1)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4),thinlev,thiwmax,thirad
        else
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2))')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4)
        endif
       else ! ( lversion .ne. 0 )
         if ( bxmag .ge. 100. ) then
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,f8.0,f8.2,''  _'',a)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      cskvers(1:lenvers) 
         else
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,f10.4,1p,e9.1,
     +      0p,f7.1,i7,f9.1,f8.1,2f8.2,''  _'',a)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectmin,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      cskvers(1:lenvers) 
         endif  
       endif
      endif
      ! - - - - - - print multi thinning parameters (optional):
      if ( 10599 .le. mrun .and. mrun .le. 10899 ) then
         if ( mthinh .lt. -1 ) then
            write(*,'(29x,''thicut'',f7.1,''cm  mthi'',
     +         i2,''  msee'',i2,''   thilev'',1p,6e9.1)')
     +         corecut,mthinh,mseeds,(qmthin(i),i=1,mthinh*3,3)
            write(*,'(63x,''weight'',1p,6e9.1)')
     +         (qmthin(i),i=2,mthinh*3,3)
         endif
      endif
      ! - - - - - - print ratio of ectmax/ectmin (run number dependent):
      !  if ( 899 .le. mrun .and. mrun .le. 899 )
      !     write(*,'(51x,2(i8,''*''))') int(0.5+ectmax/ectmin)

c - - - - - - - copy new quantities to previous ones - - - - - - - - - -
      mcodprev = mcode 
      engyprev = energy
      thetprev = theta
      phiaprev = phia
      obslprev = obslev
      lthiprev = lthi

  149 continue 

  150 continue ! end-of-loop jj=1,job

c - - - - - - print last title line of tabular - - - - - - - - - - - - -
      if ( mrun .ge. 7000 .or. mrun .lt. 100 ) then ! found coreas: 
        write(*,'(/,''primary  lg(E)  theta   phi   runnr  sizeGBy'',
     +      2x,''procs  T(days)  ecutmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',2x,''Antens'',3x,''gamma'',
     +      2x,''_corsika_executable_'')')  
      else
        write(*,'(/,''primary  lg(E)  theta   phi   runnr  sizeGBy'',
     +      2x,''procs  T(days)  ecutmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''_corsika_executable_'')')  
      endif
      write(*,'(54x,''Total'',13x,''Wall'',9x,''E/ecutmax'')')
      if ( nthi .gt. 0 )
     +   write(*,'(/,50x,''Number of thin simulations:'',i5)') nthi
      if ( nthinh .gt. 0 )
     +   write(*,'(/,50x,''Number of multithin simuls:'',i5)') nthinh

c - - - - - - end-of-program showparallel - - - - - - - - - - - - - -
      stop
      end

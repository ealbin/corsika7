c=======================================================================
c
c  t o t a l t i m e i k 3. f
c  --------------------------
c     calculate sum of used cpu times of a set of corsika simulations
c     by checking first lines of Job00ijkl_ik3.out;
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check totaltimeik3.f -o totaltimeik3
c execution in script:
c     ./postprocess-ik3.sh
c-----------------------------------------------------------------------
c example of time info:
c 1416835247 0:21.80 Real-time 0.23 TotalCPUseconds ./DAT000323-000000000-000000001.lst
c 1416835389 1:15.74 Real-time 0.53 TotalCPUseconds ./DAT000323-642460627-000000044.lst
c-----------------------------------------------------------------------
c      START TIME          STOP TIME       TIME (min)
c 1417007478.000336   1417007781.000336     63.705833
c LONGEST JOB: MPIID =  1771 and Time =  75932.855908
c  Total number of jobs =   69
c Maximum size of group =   13
c TOTAL CPU TIME (days) =    0.044240
c time.txt000336
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program totaltimeik3

      implicit double precision (a-h,o-z), integer (i-n)

      character chlinea*160, chlineb*160, chlinec*160, cjoberr*17,
     +    ctimtxt*14, czeile1*51, czeile3*37, czeile4*23,
     +    czeile5*23, czeile6*23, cfmtflt*7, cfmtint*4

      data cfmtflt/'(f 6.2)'/, cfmtint/'(i2)'/

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - initializations:
      data cjoberr/'Job000000_ik3.err'/ 
      data ctimtxt/'time.txt000000'/
      data czeile1/
     +   '     START TIME          STOP TIME       TIME (min)'/
      data czeile3/'LONGEST JOB: MPIID =       and Time ='/
      data czeile4/' Total number of jobs ='/
      data czeile5/'Maximum size of group ='/
      data czeile6/'TOTAL CPU TIME (days) ='/
      chdsek = 1.d0 / 3600.d0
      hoursum = 0.d0
      ignored = 0
      mjobs = 17
      ljob = 17

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - read Real-time lines from job-file or Job00ijkl_ik3.out: 
      read(*,'(a160)') chlinea
      read(*,'(a160)') chlineb
      do  il=3,33
         read(*,'(a160)') chlinec
         if ( index(chlinec,'initial-ik3') .gt. 0 ) goto 3
      enddo
    3 continue
      if ( index(chlinea,'error') .gt. 0 .or. 
     +     index(chlinea,'exited') .gt. 0 .or.
     +     index(chlinea,'runtime') .gt. 0 ) then
         write(*,*) ' error message in `Job00ijkl_ik3.out`. exit'
         goto 9
      endif
      read(chlinea(1:10),'(i10)') jobseca
      read(chlineb(1:10),'(i10)') jobsecb
      ib = index(chlinea,' ')
      id = index(chlinea,':')
      ia = index(chlinea,'Real-time')
      jd = index(chlinea(id+1:id+5),':') ! possible second :
      if ( jd .gt. 0 ) chlinea(id+jd:id+jd) = '.'
      write(cfmtint(3:3),'(i1)') id-ib
      jmin = 0
      if ( id-1 .ge. ib ) read(chlinea(ib:id-1),cfmtint) jmin
      read(chlinea(id+1:ia-1),cfmtflt) xsec
      jobseca = jobseca - jmin*60 - int(xsec+0.7)
      jrun = index(chlinec,'parallel') + 9
      read(chlinec(jrun:jrun+5),'(i6)') mrunnr
      write(cjoberr(4:9),'(i6.6)') mrunnr 
      write(ctimtxt(9:14),'(i6.6)') mrunnr
 
c - - - - get number of parts of the simulation from Job00ijkl_ik3.err:
      open(unit=1,file=cjoberr,form='formatted',
     +     access='sequential',status='unknown')
      read(1,*,end=4,err=4) mjobs
    4 continue
      close(unit=1)

c - - - - get totaltime in days and run number:
      open(unit=2,file='totaljobfile.out',access='sequential',
     +     form='formatted',status='old')
      read(2,*,end=5,err=5) timedays, jwallsec, wallmax, jobmax
    5 continue
      close(unit=2)
      if ( mjobs .eq. 1 ) ljob = mjobs

c - - - - write time infos to file `time.txt00ijkl`:
      open(unit=3,file=ctimtxt,form='formatted',
     +     access='sequential',status='unknown')
      write(3,'(a)') czeile1
      write(3,'(i10,''.'',i6.6,i13,''.'',i6.6,f14.6)')
     +   jobseca,mrunnr,jobsecb,mrunnr,1.d0/60.d0*jwallsec
      write(3,'(a,i6,a,f14.6)') czeile3(1:20), jobmax,
     +    czeile3(27:37), wallmax ! timedays*2.6543d3
      write(3,'(a,i5)') czeile4, mjobs
      write(3,'(a,i5)') czeile5, ljob
      write(3,'(a,f12.6)') czeile6, timedays
      write(3,'(''time.txt'',i6.6)') mrunnr
      close(unit=3)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - end of program.
    9 continue
      stop
      end

c=======================================================================
c
c  t o t a l j o b f i l e . f
c  ---------------------------
c     auxiliary program after a parallel corsika simulation by scripts
c     on the IKP computing cluster in B425 at KIT-CN campus north;
c     count number of files, sum up all times from file `job-file`.
c-----------------------------------------------------------------------
c compilation:
c     gfortran -O0 -fbounds-check totaljobfile.f -o totaljobfile 
c execution:
c     ./totaljobfile > totaljobfile.out
c-----------------------------------------------------------------------
c input `job-file`:
c 1526551423 0:17.84 Real-time 0.74 TotalCPUseconds ./DAT001409-000000000-000000001.lst
c 1526551483 0:26.79 Real-time 0.99 TotalCPUseconds ./DAT001409-000014115-000000003.lst
c 1526551490 0:27.13 Real-time 0.60 TotalCPUseconds ./DAT001409-000014115-000000004.lst
c 1526551499 0:25.03 Real-time 1.46 TotalCPUseconds ./DAT001409-000014115-000000006.lst
c 1526551512 0:27.60 Real-time 0.58 TotalCPUseconds ./DAT001409-451172026-000000007.lst
c 1526551555 0:40.76 Real-time 0.55 TotalCPUseconds ./DAT001409-572494754-000000010.lst
c 1526551591 2:38.02 Real-time 0.70 TotalCPUseconds ./DAT001409-000014115-000000002.lst
c 1526551602 1:15.66 Real-time 1.41 TotalCPUseconds ./DAT001409-042414867-000000012.lst
c-----------------------------------------------------------------------
c output `totaljobfile.out`:
c    1.024870E+03      437613      276922.0       629       874
c      time[days]   wall[sec]   maxjob[sec]     maxid    jfiles
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program totaljobfile

      implicit double precision (a-h,o-z), integer (i-n)

      character chzeile*100, cfmtflt*7, cfmtint*4

      data cfmtflt/'(f 6.2)'/, cfmtint/'(i4)'/, ifail/0/

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - get time data from time protocol file `job-file`:
      open(unit=1,file='job-file',form='formatted',status='old')
      totimin = 0. ! total time in minutes.
      wallmax = 0.
      jobmax = 0
      do  ii=1,1234567
         read(1,'(a)',end=9,err=9) chzeile
         if ( chzeile(11:11) .eq. ' ' ) then
            ib = index( chzeile, ' ' )
            id = index( chzeile, ':' )
            ia = index( chzeile, 'Real-time' )
            chzeile(id:id) = ' '
            iz = index( chzeile, ':' )
            if ( iz .gt. 0 ) then ! twice found doubledot character:
               write(cfmtint(3:3),'(i1)') id-ib
               read(chzeile(ib:id-1),cfmtint) jhrs
               write(cfmtint(3:3),'(i1)') iz-id+1
               read(chzeile(id+1:iz-1),cfmtint) jmin
               write(cfmtint(3:3),'(i1)') ia-iz+1
               read(chzeile(iz+1:ia-1),cfmtint) jsec
               walltim = 60.*jhrs + jmin + 1.d0/60.d0*jsec
            else ! once found doubledot character, i.e. time < 1 hour:
               write(cfmtint(3:3),'(i1)') id-ib
               read(chzeile(ib:id-1),cfmtint) jmin
               read(chzeile(id+1:ia-1),cfmtflt) xsec
               walltim = xsec/60. + jmin
            endif
            if ( walltim .gt. wallmax ) then
               wallmax = walltim 
               jobmax = ii
            endif
            totimin = totimin + walltim
            if ( ii .gt. 1 ) then
               read(chzeile(1:10),'(i10)') jobendt
            else
               read(chzeile(1:10),'(i10)') jobanft
               jobanft = jobanft - int(60.*walltim+1.) 
            endif
         else
            ifail = ifail + 1
         endif
      enddo
    9 continue
      close(unit=1)
      jfiles = ii - 1 ! number of files
      if ( jfiles .gt. 1 ) then
         jwallsec = jobendt - jobanft
      else
         jwallsec = int(60.*walltim)
      endif

c - - - - - writing total minutes and wall seconds (totaljobfile.out):
      write(*,'(1p,e17.6,0p,i12,f14.1,i10,i10)')
     +   totimin/1440.d0, jwallsec, 60.*wallmax, jobmax, jfiles
      write(*,'(7x,''time[days]   wall[sec]   maxjob[sec]'',
     +   ''     maxid    jfiles'')')

c - - - - - print in error case:
      if ( ifail .gt. 0 ) then 
         write(*,'(i28,''. Files in all,'',i5,'' failed.'')') ii,ifail
      endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - end-of program totaljobfile.
      stop
      end

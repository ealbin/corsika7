c=======================================================================
c
c   s h d a t a t a n k s . f
c   -------------------------
c     reduce particle data output by writing particles only at 
c     all existing SD tanks (radius few meters);
c     check all SD tank positions on matching particles,
c     shower directed to the middle tank position (i.e. tank nr. 805).
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check shdatatanks.f -o shdatatanks
c     f77 -fbounds-check shdatatanks.f -o shdatatanks
c     ifort -C -check bounds shdatatanks.f -o shdatatanks
c execution:
c     ls -1 DAT000070* | grep t -v | grep n -v > shdatatanks.i000070
c     ./shdatatanks < shdatatanks.i000070 > shdatatanks.out000070
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c          shdatatanks.pmsoutab:
c   2006 07 03      72.2   6115227.70    450615.30    1568.47
c   2006 07 05      79.3   6113927.51    452868.44    1542.87
c   2006 07 04      80.2   6115231.60    452122.20    1554.13
c   2006 07 07      89.2   6113928.98    451367.26    1555.97
c   2006 07 07      90.2   6112631.00    450619.10    1553.46
c   2006 07 07      91.2   6112631.90    452118.70    1541.04
c-----------------------------------------------------------------------
c          shdatatanks.pattern:
c
c 912 865 859 864-220-218 853 850 545 846 848-357-449-374-519-396-522-533-524
c
c   916 863 866-222-219-211 849 827 830 833-371-372-525-527-528-523-305-354-375
c
c-737-740-743-221-223-217 843 832 839 836 841-526-419-420-424-498-456-463-448-452
c
c  -747-741-745-818-797-210-213 840 818 834-421-425-426-428-422-460-491-464-393
c
c-729-735-746-819-791-817-798-209-215 835 837-429-432-433-430-515-492-473-807-803
c
c  -724-721-751-752-695-719-717 820 809 805-445-436-446-440-412-499-534-345-521
c
c-727-722-757-793-794-795-796-232 829 826 838-447-821-443-285-517-415-404-416-288
c
c   779-749-748-754-758-834-828-234-235-233-441-444-437 -12-414-467-720-718-238
c
c 586 593 749-799-750-753   0 828 813 825 811 814 810 -13 -14-529-459-714-318-236
c
c   587 551 795 791 797 816-800 823 831 509 511 503 505 -19 -15-713-716 692 721
c
c 563 568 560 558 565 798 796 817 822 502 523 510 506 -16 -17-715-320 696 691 697
c
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program shdatatanks

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (nfmx=10000,lenthin=6552,lenstnd=5733) 

      character cfilename(nfmx)*120, cfilenout(20)*120, chposit*25

      dimension qtnumb(2000),qynort(2000),qxeast(2000),qtalti(2000)

      dimension qtanks(100),xcore(20),ycore(20),rcore(20)

      real pdata(lenthin), rdata(lenthin), qdata(1560), ddata(312)

      logical lexist

      common /urdat/pdata,rdata,ddata,qdata,qtanks
      common /utabl/qtnumb,qynort,qxeast,qtalti,dradius,qradius,
     + xshcore,yshcore,tpoi,nfil,ifil,irec,isho,lrec,lsub,lpoi,lcor,
     + ltks,ntks

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - define radius and position of the shower center: 
      ! number of the tank as center of the array.
      midl = 805
      ! used radius [m] for data around a detector tank.
      dradius =  3.0
      qradius = dradius * dradius ! [m^2]

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read tabular of built tanks (from pms data base):
      open(unit=3,file='shdatatanks.pmsoutab',form='formatted',
     +     status='old',access='sequential')
      do  ltks=1,2000
         read(3,*,end=91,err=91) iy,im,id,
     +      qtnumb(ltks),qynort(ltks),qxeast(ltks),qtalti(ltks)   
         if ( int(qtnumb(ltks)) .eq. midl ) mtks = ltks
***      write(8,'(5x,''+'',f7.1,''d0,'',f11.2,''d0,'',f10.2,''d0,'',
***  +      f8.2,''d0,'')') qtnumb(ltks),
***  +      qynort(ltks),qxeast(ltks),qtalti(ltks)
      enddo
   91 continue
      close(unit=3)
      ltks = ltks - 1
      xmiddle = qxeast(mtks)
      ymiddle = qynort(mtks)
      do  l=1,ltks
         qxeast(l) = qxeast(l) - xmiddle
         qynort(l) = qynort(l) - ymiddle
      enddo

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read all names of particle data files (up to nfmx):
      do  ifil=1,nfmx
         read(*,'(a)',end=102,err=101) cfilename(ifil) 
      enddo
  101 continue
      write(*,*) '    ERROR: reading list of file names.'
  102 continue
      nfil = ifil - 1
      lfil = 120
  103 continue
      lfil = lfil - 1
      if ( cfilename(nfil)(lfil:lfil) .eq. ' ' ) goto 103  
      read(cfilename(1)(4:9),'(i6)') lrun ! from name DAT000407-000001

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      chposit = 'shdatapositions.tab000000'
      write(chposit(20:25),'(i6.6)') lrun
      inquire(file=chposit,exist=lexist)
      if ( .not.lexist ) then
         ! x-distance [m] of shower core relatively to tank `midl`.
         xshcore = 142.
         ! y-distance [m] of shower core relatively to tank `midl`.
         yshcore = 142.
         ! number of core positions.
         xcore(1) = xshcore
         ycore(1) = yshcore
         rcore(1) = sqrt(xshcore*xshcore+yshcore*yshcore)
         lcor = 1
      else
         open(unit=2,file=chposit,form='formatted',
     +        status='old',access='sequential')
         do  lcor=1,20
            read(2,'(24x,f9.3,7x,f9.3,7x,f9.3)')
     +         xcore(lcor),ycore(lcor),rcore(lcor)
         enddo 
         close(unit=2)
         lcor = lcor - 1
      endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read first record of first file:
      open(unit=1,file=cfilename(1)(1:lfil),form='unformatted',
     +     status='old',access='sequential')
      read(unit=1,err=496,end=495) (pdata(i),i=1,lenstnd)
      close(unit=1)
c - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         stop ! goto 498
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39
      do  i=1,936 ! keep first two subblocks and first particle subblock.
         qdata(i) = pdata(i)
      enddo
      qdata(lenblk+268) = dradius * 1.d2 ! radius of circle in cm.
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(12x,
     +      ''corsika particle data reduction `thinning run`'')')
      else
         write(*,'(12x,
     +      ''corsika particle data reduction `standard run`'')')
      endif
      write(*,'(12x,''radius around a detector ='',
     +   f7.2,'' m'')') dradius
      ncores = lcor
      if ( ncores .gt. 20 ) ncores = 20
      if ( ncores .gt. 1 ) 
     +   write(*,'(12x,''number of core positions ='',i7)') ncores
c - - - - - - calculate core dependent file names:    
      do  lcor=1,ncores
         xshcore = xcore(lcor)
         yshcore = ycore(lcor)
         lxshc = min( int(abs(xshcore)), 999)
         lyshc = min( int(abs(yshcore)), 999)
         if ( xshcore .ge. 0. ) then
            if ( yshcore .ge. 0. ) then
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''+'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            else
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''-'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            endif
         else
            if ( yshcore .ge. 0. ) then
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''-'',i3.3,''+'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            else
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''-'',i3.3,''-'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            endif
         endif
         if ( lcor .ge. 0 ) write(*,'(12x,a,''  x='',f8.2,''  y='',
     +      f8.2,''  r='',f8.2)') cfilenout(lcor)(1:28),
     +      xshcore, yshcore, sqrt(xshcore*xshcore+yshcore*yshcore)
      enddo
      if ( nfil .gt. 1 ) write(*,'(30x,''nfiles ='',i7)') nfil

c - - - - - loop on known core positions to sort matching particles:
      do  448  lcor=1,ncores
      write(*,'(1x,44(''-''))')
      write(*,'(3x,a,''  r='',f8.2)') cfilenout(lcor)(1:28),rcore(lcor)
      xshcore = xcore(lcor)
      yshcore = ycore(lcor)

c - - - - - - open new (shorter) particle data file:
      open(unit=4,file=cfilenout(lcor)(1:28),form='unformatted',
     +     status='unknown',access='sequential')
      tpoi = 0.d0
      tsum = 0.d0 
      nrec = 0
      ntks = 0
      isho = 0
      lrec = 0
      lpoi = 0
      do  ir=1,lenblk*2 ! copy runh and evth subblocks.
         rdata(ir) = qdata(ir)
      enddo
      lsub = 2
      lp = 0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - loop on names of particle data files:
      do  444  ifil=1,nfil
      if ( ifil .gt. 1 ) close(unit=1)
      if ( nfil .gt. 1 )
     +   write(*,'(i6,''. '',a)') ifil,cfilename(ifil)(1:lfil)
      iend = 0
      irec = 0
      open(unit=1,file=cfilename(ifil)(1:lfil),form='unformatted',
     +     status='old',access='sequential')
c - - - - - - read data record with lenrec words - - - - -
  431 continue
      read(unit=1,err=432,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      ! if (mod(irec,200) .eq. 0) write(*,'(12x,''irec ='',i11)') irec
      call blwork(iend,lenrec,lp)
      if ( iend .le. 0 ) goto 431
  432 continue
c - - - - - - end-of particle data file - - - -
      nrec = nrec + irec
      tsum = tsum + 819.d0*irec
      if ( irec .gt. 1 .or. ifil .eq. nfil ) then
         write(*,'(10x,f14.0,'' particles of'',f14.0,
     +      '' (ratio'',f9.6,'')'')') tpoi,tsum,tpoi/tsum
         if ( qdata(lenblk+152) .gt. 0. )
     +      write(*,'(25x,''particles outside'',f4.0,
     +      '' meter of the shower center.'')') qdata(lenblk+152)*1.d-2
         write(*,'(12x,i12,'' old records'')') nrec
         write(*,'(12x,i12,'' new records'')') lrec
         call sortpair(qtanks,ntks)
         write(*,'(25x,''tank_nrs'',/,8(25x,10i5,:,/))')
     +      (int(qtanks(i)),i=1,ntks)     
      endif
c - - - - - - end of corsika particle data file reached: 
  444 continue
      close(unit=1)

c - - - - - - end of loop on known core positions:
  448 continue
      goto 499
 
c--end of data----------------------------------------------------------
  495 continue
  496 continue
  497 continue
  498 continue
  499 continue
      stop
      end
c=======================================================================
c
c     sorting vector by paired transpositions
c
c-----------------------------------------------------------------------

      subroutine sortpair(qtanks,ntks)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension qtanks(ntks)

      lz = ntks
      do  6  la=1,lz
      lc = 0
      do  5  lb=2,lz-la+1
         if ( qtanks(lb-1) .gt. qtanks(lb) ) then
            qauxil = qtanks(lb)
            qtanks(lb) = qtanks(lb-1)                         
            qtanks(lb-1) = qauxil
         else
            lc = lc + 1
         endif
         if ( lc .eq. lz-la ) goto 7
    5 continue
    6 continue
    7 continue

      return
      end
c=======================================================================
c
c     analyze the contents of all 21 subblocks in a record
c
c-----------------------------------------------------------------------
 
      subroutine blwork(iend,lenrec,lp)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)

      dimension qtnumb(2000),qynort(2000),qxeast(2000),qtalti(2000)

      dimension qtanks(100)

      real pdata(lenthin), rdata(lenthin), qdata(1560), ddata(312)

      common /urdat/pdata,rdata,ddata,qdata,qtanks
      common /utabl/qtnumb,qynort,qxeast,qtalti,dradius,qradius,
     + xshcore,yshcore,tpoi,nfil,ifil,irec,isho,lrec,lsub,lpoi,lcor,
     + ltks,ntks

      lenblk = lenrec / 21
      lenpar = lenblk / 39

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------
      elseif ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
c----------------subblock longi information-----------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then
c----------------subblock event end-------------------------------------

         if ( mod(lp,lenblk) .eq. 1 ) then 
            ip = mod(int(1.*(lp-lia)/lenpar),39)
            if ( ip .eq. 0 ) ip = 39
***         write(*,'(12x,''irec ='',i11,'' last;  _particle='',i2.2,
***  +      ''   isub='',i2.2,i10,a)') irec,ip,
***  +         int((1.*lenblk+lia)/lenblk),lsub,'=lsub' 
         endif

c - - - - - - - first data file - - - - - - - - - - - - - - - - - - - -
         if ( ifil .eq. 1 ) then
            ! - - - - keep evte and rune subblocks of first file:
            do  i=1,lenblk
               qdata(lenblk*3+i) = pdata(lia+i-1)
            enddo
            do  i=1,lenblk
               qdata(lenblk*4+i) = pdata(lia+i-1+lenblk)
            enddo
            qdata(lenblk*4+3) = -1.*nfil
            qdata(lenblk*4+4) = dradius * 1.d2 ! radius of circle in cm.
         endif

c - - - - - - - last data file - - - - - - - - - - - - - - - - - - - - -
         if ( ifil .eq. nfil ) then
            write(*,'(43x,''last_particle='',i2.2,i10,''=lsub'')')
     +         lpoi,lsub+1
            if ( lsub .ge. 21 ) then
               ! - - - - possibly write full record to file:
               write(unit=4) (rdata(i),i=1,lenrec)
               lrec = lrec + 1
               do  i=1,lenrec
                 rdata(i) = 0.
               enddo
               lsub = 0
            endif 
            if ( lpoi .gt. 0 ) then
               ! - - - - copy end of particle data to buffer:
               lsub = lsub + 1
               is = lenblk * (lsub-1)
               do  ir=1,lenpar*lpoi
                  rdata(is+ir) = ddata(ir)
               enddo
               lpoi = 0
            endif
            if ( lsub .ge. 21 ) then
               ! - - - - possibly write full record to file:
               write(unit=4) (rdata(i),i=1,lenrec)
               lrec = lrec + 1
               do  i=1,lenrec
                 rdata(i) = 0.
               enddo
               lsub = 0
            endif
            ! - - - - write evte subblock:
            lsub = lsub + 1
            is = lenblk * (lsub-1)
            do  ir=1,lenblk  
               rdata(is+ir) = qdata(ir+lenblk*3)
            enddo
            if ( lsub .ge. 21 ) then
               ! - - - - possibly write full record to file:
               write(unit=4) (rdata(i),i=1,lenrec)
               lrec = lrec + 1
               do  i=1,lenrec
                 rdata(i) = 0.
               enddo
               lsub = 0
            endif
            ! - - - - write rune subblock:
            lsub = lsub + 1
            is = lenblk * (lsub-1)
            do  ir=1,lenblk
               rdata(is+ir) = qdata(ir+lenblk*4)
            enddo
            ! - - - - write last record to particle data file:
            write(unit=4) (rdata(i),i=1,lenrec)
            lrec = lrec + 1
            ! - - - - increase iend to mark end of data.
            iend = iend + nfil
         endif

c - - - - - - - loop 948 finished - - - - - - - - - - - - - - - - - - -
         goto 949

c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + irec
         goto 949
      else
c-----------------subblock with particle data---------------------------
c - - - - - - - - test on matching positions around tanks:
        do  944  lp=lia,lia+lenblk-1,lenpar  
          if ( pdata(lp).ge.75000. .and. pdata(lp).lt. 77000. ) goto 944
          if ( pdata(lp).gt.    0. .and. pdata(lp).lt.196000. ) then
c - - - - - - - - shift particle coordinates by core position:
            xp = 1.d-2 * pdata(lp+4) - xshcore  
            yp = 1.d-2 * pdata(lp+5) - yshcore 
c - - - - - - - - test all current tank positions:
            do  942  l=1,ltks
              pdist = (xp-qxeast(l))*(xp-qxeast(l)) +
     +                (yp-qynort(l))*(yp-qynort(l))
              if ( pdist .le. qradius ) then
c - - - - - - - - - - matching tank position:
                if ( tpoi .eq. 0. ) then
                  if ( qdata(lenblk+152) .gt. 0. ) then
                    write(*,'('' icode'',6x,''px'',12x,''py'',12x,
     +               ''pz'',11x,''x[m]'',6x,''y[m]'',7x,''t[nsec]'',
     +                4x,''weight'',3x,''dist[m]'',4x,''tank'')')  
                  else
                    write(*,'('' icode'',6x,''px'',12x,''py'',12x,
     +               ''pz'',11x,''x[m]'',6x,''y[m]'',7x,''t[nsec]'',
     +                3x,''dist[m]'',4x,''tank'')')  
                  endif   
                  ntks = 1
                  qtanks(ntks) = qtnumb(l)
                endif
                if ( tpoi .le. 11. ) then
                  icode = int(pdata(lp)*1.d-3)
                  if ( lcor .ge. 1 .and. pdata(lp+lenpar-1) .ge. 1.)
     +              pdata(lp+lenpar-1) = 1.d-6*pdata(lp+lenpar-1) + lcor
                  pdata(lp+4) = pdata(lp+4) * 1.d-2 
                  pdata(lp+5) = pdata(lp+5) * 1.d-2 
                  write(*,'(i5,''.'',1p,3e14.6,0p,2f10.2,f13.7,
     +              f10.2,f8.0)') icode,
     +              (pdata(i),i=lp+1,lp+lenpar-1),sqrt(pdist),qtnumb(l)
                endif
                do  i=1,ntks ! test on additional tank number.
                  if ( qtnumb(l) .eq. qtanks(i) ) goto 111  
                enddo
                ntks = ntks + 1
                qtanks(ntks) = qtnumb(l)
  111           continue       
                tpoi = tpoi + 1. ! count points.
                lpoi = lpoi + 1
                is = lenpar * (lpoi-1) 
                do  i=1,lenpar
                  ddata(is+i) = pdata(i+lp-1)
                enddo 
                if ( lpoi .eq. 39 ) then
                  lsub = lsub + 1
                  is = lenblk * (lsub-1) ! element index in subblock.
                  do  ir=1,lenblk
                    rdata(is+ir) = ddata(ir)
                  enddo
                  lpoi = 0
                  if ( lsub .eq. 21 ) then
                    ! - - - - write to new particle data file:
                    write(unit=4) (rdata(i),i=1,lenrec)
                    lrec = lrec + 1
                    do  i=1,lenrec
                      rdata(i) = 0.
                    enddo
                    lsub = 0
                  elseif ( lsub .lt. 21 ) then
                  ! write(*,*) ' record not full. '
                  endif
                elseif ( lpoi .lt. 39 ) then ! subblock not full.
                ! write(*,*) ' subblock not full. '
                endif
              endif
  942       continue
          elseif ( pdata(lp) .eq. 0. ) then
            ! end of particle data found within this subblock: 
c           write(*,'(12x,''irec ='',i11,'' last;   particle='',i2.2,
c    +         ''   isub='',i2.2,i10,a)') irec,int(1.*(lp-lia)/lenpar),
c    +         int((1.*lenblk+lia)/lenblk),lsub,'=lsub'
            goto 946 
          elseif ( pdata(lp) .ge. 201000. ) then
          ! write(*,*) '        irec =',irec,'   icod =',pdata(lp)*1.e-3
          endif
  944   continue
      endif
  946 continue
      if ( iend .gt. 0 ) goto 949
  948 continue ! end-of-loop over subblocks (lia=1,lenrec,lenblk).
 
c-----------end of record or end of particle data reached---------------
  949 continue
      return
      end

c=======================================================================
c
c  s o r t a u g e r h i t . f
c  ---------------------------
c     sort reduced particle data output (air shower simulation using
c     keyword AUGSCT) by writing particles separately sorted to  
c     max. 20 matching core positions.
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
c     ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
c RunProg:
c     # (a) write file names:
c     ls -1 DAT000070-?????????-????????? > sortaugerhit.i000070
c     ls -1 DAT000070-?????? > sortaugerhit.i000070
c     ls -1 DAT000070-* | grep t -v | grep n -v | grep a -v > sortaugerhit.i000070
c     # (b) run augerhit sorting program:
c     ./sortaugerhit < sortaugerhit.i000070 > sortaugerhit.out000070
c     # mv fort.20 DAT000070.prtcls20 # former test case.
c out-files:
c     ./DAT000070-augerhit*
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program sortaugerhit

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (nfmx=10000,lenthin=6552,lenstnd=5733) 
      character cfilename(nfmx)*120, cfilenout(20)*120
      double precision qdata(1560), xcore(20), ycore(20), rnd(2)
      real pdata(lenthin), rdata(lenthin), ddata(312)
      logical lexist, first

      common /urdat/pdata,rdata,ddata,qdata,xdistanz,ydistanz,hdistanz
      common /utabl/dradius,qradius,xshc,yshc,tpoi,
     + nfiles,ifil,irec,isho,lcor,lrec,lsub,lpoi,ltks,ncores,
     + muaddi,mzeros,mtoobig
      data first/.true./

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read all names of particle data files (up to nfmx):
      write(*,'(1x,6(''_ ''),''sortaugerhit.f'',25('' _''))')
      ifil = 1
      read(*,'(a)',end=108,err=107) cfilename(1)
      ! - - - check length of file name:
      lfil = 120 
  102 continue
      lfil = lfil - 1
      if ( cfilename(1)(lfil:lfil) .eq. ' ' ) goto 102  
      ! - - - read following names of particle data files:
  104 continue
      ifil = ifil + 1
      read(*,'(a)',end=108,err=107) cfilename(ifil)
      inquire(file=cfilename(ifil)(1:lfil),exist=lexist)
      if ( lexist ) goto 104
      write(*,*) '    ERROR: wrong name in the list of file names:'
      write(*,*) cfilename(ifil)
      ifil = ifil - 1
      goto 104
  107 continue
      write(*,*) '    ERROR: reading names of particle data files.'
  108 continue
      idat = index(cfilename(1),'DAT')
      if ( lfil .ne. 9 .and. lfil .ne. 14 .and. lfil .ne. 16 ) then
         if ( index(cfilename(1),'_') .eq. 10 ) irun = 987654 
      else
         irun = 999999
         if ( idat .ge. 1 ) 
     +      read(cfilename(1)(idat+3:idat+8),'(i6)') irun
      endif
      nfiles = ifil - 1

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read first record of first file to test simulation:
      open(unit=1,file=cfilename(1)(1:lfil),status='old',
     +     form='unformatted',access='sequential')
      read(unit=1,err=496,end=495) (pdata(i),i=1,lenstnd)
      close(unit=1)
      if ( index(cfilename(1),'_') .eq. 10 ) irun = int(pdata(2))
      if ( lfil .eq. 9 .or. lfil .eq. 14 ) cfilename(1)(10:10) = '-'
c - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      else if ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 498
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
      do  i=1,936 ! keep first three subblocks.
         qdata(i) = pdata(i)
      enddo
      do  i=937,1560 ! init for evte and rune subblocks.
         qdata(i) = 0.
      enddo
      if ( qdata(274)+qdata(547) .lt. 1. ) then
       write(*,'(/,12x,''corsika particle data as `thinning run`'')')
      else
       write(*,'(/,12x,''corsika particle data as `standard run`'')')
      endif
      write(*,'(12x,''height of 1st interaction [m] ='',f12.3)')
     +   abs(qdata(lenblk+7)*1.d-2)
      write(*,'(12x,''number of particle data files ='',i7)') nfiles
c - - - - - - test/copy augerhit parameters and calculate file names:
      if ( qdata(lenblk+98) .le. 0.d0 )
     + then ! number of scattered cores originally not set.
         write(*,'(12x,''INFO: no core positions defined;'',
     +      '' calculate new ones:'')')
         ncores = 20
         qdata(lenblk+98) = -ncores 
         qdata(lenblk+176) = ncores
         dradius = 35.d0 ! meter
         qdata(lenblk+174) = dradius*1.d2 ! in cm.
         xdistanz = 1500.d0 ! meter
         qdata(lenblk+175) = xdistanz*1.d2 ! in cm.
      else ! use number of scattered cores set at the simulation:
         ncores = int(qdata(lenblk+98))
         dradius = qdata(lenblk+174)*1.d-2 ! switch to meter. 
         if ( dradius .lt. 5.d0 ) dradius = 35.d0
         xdistanz = qdata(lenblk+175)*1.d-2 ! switch to meter.
         if ( xdistanz .lt. 5.d0 ) xdistanz = 1500.d0 ! Auger dist. 
      endif
      qradius = dradius * dradius ! [m^2]
      hdistanz = 0.5d0 * xdistanz
      ydistanz = sqrt(3.d0) * hdistanz
      xydist = hdistanz * 0.5d0
c - - - - - - calculate new positions of scattered cores:
      write(*,'(2x)')
      if ( qdata(lenblk+98) .lt. 0.d0 ) then
        if ( first ) then
          ! initialize sobol number generator SOBSEQ:
          first = .false.
          call sobseq( -2,rnd )
          do  i=1,irun ! iseed(1,3) or irunnr
            call sobseq( 2,rnd )
          enddo
        endif
        do  lcor=1,ncores
    2     continue
          call sobseq( 2,rnd ) ! get pair of quasi random numbers.
          xshc = hdistanz * (2.d0*rnd(1) - 1.d0)
          yshc = ydistanz * (2.d0*rnd(2) - 1.d0)
          ! test distance to central detector position
          if ( xshc*xshc + yshc*yshc .lt. 2.000d4 ) goto 2 ! r>200m
          if ( xshc*xshc + yshc*yshc .gt. 5.678d5 ) goto 2 ! r<753m
          ! now point is within the hexagon, outside the inner circle.
          xcore(lcor) = xshc
          ycore(lcor) = yshc 
          if ( sqrt(xshc*xshc+yshc*yshc) .gt. xydist )
     +       xydist = sqrt(xshc*xshc+yshc*yshc)
          lxshc = min( int(abs(xshc)), 999)
          lyshc = min( int(abs(yshc)), 999)
          if ( xshc .ge. 0. ) then
            if ( yshc .ge. 0. ) then
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''+'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            else
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''-'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            endif
          else
            if ( yshc .ge. 0. ) then
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
     +      xshc, yshc, sqrt(xshc*xshc+yshc*yshc)
        enddo
        write(*,'(2x)')
c - - - - - - use available scattered cores set before the simulation:
      else
        do  lcor=1,min(ncores,20)
         xshc = qdata(lenblk+98+lcor) * 1.d-2
         yshc = qdata(lenblk+118+lcor)* 1.d-2
         xcore(lcor) = xshc 
         ycore(lcor) = yshc
         if ( sqrt(xshc*xshc+yshc*yshc) .gt. xydist )
     +      xydist = sqrt(xshc*xshc+yshc*yshc)
         lxshc = min( int(abs(xshc)), 999)
         lyshc = min( int(abs(yshc)), 999)
         if ( xshc .ge. 0. ) then
            if ( yshc .ge. 0. ) then
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''+'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            else
               write(cfilenout(lcor)(1:28),
     +            '(a,''augerhit'',i2.2,''+'',i3.3,''-'',i3.3)')
     +            cfilename(1)(1:10), lcor, lxshc, lyshc
            endif
         else
            if ( yshc .ge. 0. ) then
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
     +      xshc, yshc, sqrt(xshc*xshc+yshc*yshc)
        enddo
        write(*,'(2x)')
      endif
      write(*,'(12x,''number of core positions ='',i7)') ncores
      write(*,'(12x,''radius around a detector ='',f7.2,'' m'')')dradius
      write(*,'(12x,''   maximum core distance ='',f7.2,'' m'')') xydist
      tpoisum = 0.d0
      trecsum = 0.d0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - plot core positions in x-y-grid (printer version):
      call plotcorepos(ncores,xcore,ycore)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - sort particles to matching core positions: 
      do  449  lcor=1,min(ncores,20)
      xshc = xcore(lcor)
      yshc = ycore(lcor)
      qdata(lenblk+176) = 1.d0*lcor ! current number of core positions.
      if ( lcor .gt. 1 ) close(unit=4)
      write(*,'(/,'' _________'',i3,''. core position   x='',
     +   f8.2,'' m   y='',f8.2,'' m   r='',f8.2,'' m  _________'')')
     +   lcor,xshc,yshc,sqrt(xshc*xshc+yshc*yshc)
      write(*,'(8x,a)') cfilenout(lcor)(1:28)
      muaddi = 0
      mzeros = 0
      mtoobig= 0

c - - - - - - open new core dependent particle data file:
      open(unit=4,file=cfilenout(lcor)(1:28),form='unformatted',
     +     status='unknown',access='sequential')
      tsum = 0.d0 
      tpoi = 0.d0
      nrec = 0
      lrec = 0
      lsub = 0
      lpoi = 0
      isho = 0
      do  ir=1,lenblk*2 ! copy runh and evth subblocks.
         rdata(ir) = real(qdata(ir))
      enddo
      lsub = 2
      lp = 0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - loop on names of particle data files:
      do  444  ifil=1,nfiles
      if ( ifil .gt. 1 ) close(unit=1)
      iend = 0
      irec = 0
      open(unit=1,file=cfilename(ifil)(1:lfil),form='unformatted',
     +     status='old',access='sequential')
c - - - - - - read data record with lenrec words - - - - -
  431 continue
      read(unit=1,err=432,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      ! if ( mod(irec,4000).eq.0 ) write(*,'(12x,''irec ='',i11)') irec
      call blwork(iend,lenrec,lp)
      if ( iend .le. 0 ) goto 431
  432 continue
c - - - - - - end-of corsika particle data file reached - - - -
      nrec = nrec + irec
      tsum = tsum + 819.d0*(-0.5d0+irec) - 80. ! estim. total particles.
      if ( ifil .eq. nfiles ) then
         if ( int((tpoi+818.d0)/819.d0) .gt. 1 ) then
            write(*,'(10x,f14.0,'' particles of'',f14.0,i12,
     +         '' records'')') tpoi,tsum,int((tpoi+818.d0)/819.d0)
         else
            if ( tpoi .ge. 1.d0 ) 
     +         write(*,'(10x,f14.0,'' particles of'',f14.0,i12,
     +            '' record '')') tpoi,tsum,int((tpoi+818.d0)/819.d0)
         endif
         tpoisum = tpoisum + tpoi
         trecsum = trecsum + int((tpoi+818.d0)/819.d0)
      endif
  444 continue
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - end-of loop on names of particle data files.
      close(unit=1) ! last close of read unit.
      ! if ( lcor .eq. ncores ) write(*,'(/,16x,''muaddi,mzeros,'',
      ! +   ''mtoobig='',3i11)') muaddi,mzeros,mtoobig
  449 continue
      write(*,'(/,10x,f14.0,'' particles total'',i23,
     +   '' records total'',/)') tpoisum,int(trecsum)
      goto 499
 
c--end of data----------------------------------------------------------
  495 continue
  496 continue
  498 continue
  499 continue
      write(*,'(1x,38(''_ ''),/)')
      stop
      end
c=======================================================================
c
c     analyze the contents of all 21 subblocks in a record
c
c-----------------------------------------------------------------------
 
      subroutine blwork(iend,lenrec,lp)

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (lenthin=6552)
      double precision qdata(1560)
      real pdata(lenthin), rdata(lenthin), ddata(312)
      common /urdat/pdata,rdata,ddata,qdata,xdistanz,ydistanz,hdistanz
      common /utabl/dradius,qradius,xshc,yshc,tpoi,
     + nfiles,ifil,irec,isho,lcor,lrec,lsub,lpoi,ltks,ncores,
     + muaddi,mzeros,mtoobig
      lenblk = lenrec / 21
      lenpar = lenblk / 39

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------

      else if ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------

      else if (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
c----------------subblock longi information-----------------------------

      else if ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then
c----------------subblock event end-------------------------------------
         if ( mod(lp,lenblk) .eq. 1 ) then 
            ip = mod(int(1.*(lp-lia)/lenpar),39)
            if ( ip .eq. 0 ) ip = 39
         endif
c - - - - - - - first data file - - - - - - - - - - - - - - - - - - - -
         if ( ifil .eq. 1 ) then
            ! - - - - keep evte subblock of first file:
            do  i=1,lenblk
               qdata(lenblk*3+i) = pdata(lia+i-1)
            enddo
            ! - - - - set new rune subblock:
            qdata(lenblk*4+1) = 3301.3325195312500000d0
            qdata(lenblk*4+2) = qdata(2)
            qdata(lenblk*4+3) = 1.*lcor
            qdata(lenblk*4+4) = dradius*1.d2
            qdata(lenblk*4+5) = qdata(lenblk+176)
            qdata(lenblk*4+6) = xshc
            qdata(lenblk*4+7) = yshc
         endif
c - - - - - - - last data file - - - - - - - - - - - - - - - - - - - - -
         if ( ifil .eq. nfiles ) then
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
               rdata(is+ir) = real(qdata(ir+lenblk*3))
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
               rdata(is+ir) = real(qdata(ir+lenblk*4))
            enddo
            ! - - - - write last record to particle data file:
            write(unit=4) (rdata(i),i=1,lenrec)
            lrec = lrec + 1
            ! - - - - increase iend to mark end of data.
            iend = iend + nfiles
         endif
c - - - - - - - loop 948 finished - - - - - - - - - - - - - - - - - - -
         goto 949

c----------------subblock run end---------------------------------------
      else if ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + irec
         goto 949
      else
c-----------------subblock with particle data---------------------------
c - - - - - - - - incl. test on matching positions around tanks:
        do  944  lp=lia,lia+lenblk-1,lenpar  
          if ( pdata(lp).le. 0. ) then
             mzeros = mzeros + 1
             goto 944
          endif
          if ( pdata(lp).ge. 75000. .and. pdata(lp).lt. 97000. ) then
             muaddi = muaddi + 1
             goto 944
          endif
          if ( pdata(lp).gt. 0. .and. pdata(lp).lt. 196000. ) then
             icode = int(pdata(lp)*1.d-3)
c - - - - - - - - switch particle coordinates to meter (!!):
             pdata(lp+4) = real(1.d-2 * pdata(lp+4))
             pdata(lp+5) = real(1.d-2 * pdata(lp+5)) 
c - - - - - - test particle coordinates on horizontal stripes:
             yp = pdata(lp+5) - yshc
             iytst = nint(yp / ydistanz)
             ytest = yp - ydistanz*iytst
             if (.not.(-dradius .le. ytest .and. ytest .le. dradius) )
     +          goto 944
c - - - - - - test coordinates on 60 degrees stripes:
             xp = pdata(lp+4) - xshc
             xtest = xp + hdistanz * mod(iytst,2)
             ixtst = nint(xtest / xdistanz) 
             xtest = xtest - xdistanz*ixtst
             if (.not.(-dradius .le. xtest .and. xtest .le. dradius) )
     +          goto 944
             pdist = xtest*xtest + ytest*ytest
             if ( pdist .le. qradius ) then
c - - - - - - - - - - particle belongs to lcor-th core position:  
                if ( tpoi .eq. 0. ) then
                  if ( lenpar .eq. 7 ) then
                    write(*,'('' icode'',5x,''px'',10x,''py'',10x,
     +               ''pz'',10x,''x[cm]'',5x,''y[cm]'',7x,''t[psec]'',
     +                2x,''dist[m]'',4x,''del(x)'',4x,''del(y)'')')
                  else ! lenpar=8:
                    write(*,'('' icode'',5x,''px'',10x,''py'',10x,
     +               ''pz'',10x,''x[cm]'',5x,''y[cm]'',7x,''t[psec]'',
     +                3x,''weight'',2x,''dist[m]'',4x,''del(x)'',
     +                4x,''del(y)'')')
                  endif
                endif
                ! test print of first lcor-th particles: 
                if ( tpoi .le. 3. ) then
                  pdata(lp+6) = real(1.d-6*pdata(lp+6)+lcor)
                  if ( lenpar .eq. 7 ) then
                    write(*,'(i5,''.'',1p,3e12.4,
     +                0p,2f10.2,f14.9,f9.2,2f10.2)') icode,
     +                (pdata(i),i=lp+1,lp+lenpar-1),sqrt(pdist)
     +                ,pdata(lp+4)-xshc,pdata(lp+5)-yshc
                  else ! lenpar=8:
                    write(*,'(i5,''.'',1p,3e12.4,
     +                0p,2f10.2,f14.9,2f9.2,2f10.2)') icode,
     +                (pdata(i),i=lp+1,lp+lenpar-1),sqrt(pdist)
     +                ,pdata(lp+4)-xshc,pdata(lp+5)-yshc
                  endif
                endif
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
                  else if ( lsub .lt. 21 ) then
                  ! write(*,*) ' record not full. '
                  endif
                else if ( lpoi .lt. 39 ) then ! subblock not full.
                ! write(*,*) ' subblock not full. '
                endif
             else
                mtoobig = mtoobig + 1
             endif
          else if ( pdata(lp) .eq. 0. ) then
            ! end of particle data found within this subblock: 
            write(*,*) '        irec =',irec,'   null =',pdata(lp)*1.e-3
            goto 946 
          else if ( pdata(lp) .ge. 201000. ) then
            if ( lcor .eq. ncores) write(21,*) '        ifil =',ifil,
     +         '   irec =',irec,'   icod =',pdata(lp)*1.e-3
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

c=======================================================================
c
c     plot core positions of scattered shower centers (800m*800m)
c
c-----------------------------------------------------------------------

      subroutine plotcorepos(mcores,xcore,ycore)

      implicit double precision (a-h,o-z), integer (i-n)
      double precision xcore(20), ycore(20)
      character charray(-80:80,-40:40)*1
      data charray/13041*' '/

c - - - - - - define constants and initialize grafic array:
      do  k=20,-20,-1
         charray(0,k) = '|'
         if ( mod(k,5) .eq. 0) charray(0,k) = '+'
         charray(40,k) = '|'
         if ( mod(k,5) .eq. 0) charray(40,k) = '+'
         charray(-40,k) = '|'
         if ( mod(k,5) .eq. 0) charray(-40,k) = '+'
      enddo 
      do  i=-40,40
         charray(i,0) = '-'
         if ( mod(i,10) .eq. 0) charray(i,0) = '+'
         charray(i,40) = '-'
         if ( mod(i,10) .eq. 0) charray(i,40) = '+'
         charray(i,-40) = '-'
         if ( mod(i,10) .eq. 0) charray(i,-40) = '+'
      enddo
      charray(0,0) = 'T'
      write(*,'(2x)')

c - - - - - - calculate mcores centers of showers: 
      do  i=1,mcores
         lx = int((xcore(i)+10.)*5.0d-2) ! 20 units / char. 
         ly = int((ycore(i)+20.)*2.5d-2) ! 40 units / line.
         charray(lx,ly) = '@' 
      enddo

c - - - - - - display printer plot of core positions:      
      write(*,'(42x,''meter'')')
      write(*,'(i5,8i10)') (100*i,i=-8,8,2)
      write(*,'(4x,81a1)') (charray(i,40),i=-40,40)
      do  k=19,-19,-1
         if ( mod(k,5) .eq. 0 ) then
             write(*,'(i4,81a1)') k*40,(charray(i,k),i=-40,40)   
         else 
            write(*,'(4x,81a1)') (charray(i,k),i=-40,40)
         endif
      enddo
      write(*,'(4x,81a1)') (charray(i,40),i=-40,40)
      write(*,'(i5,8i10)') (100*i,i=-8,8,2)
      write(*,'(42x,''meter'')')

      return
      end
c=======================================================================

      SUBROUTINE SOBSEQ( N,XX )

c-----------------------------------------------------------------------
c  SOBOL QUASI RANDOM NUMBER GENERATOR
c  REFERENCE : NUMERICAL RECIPES, W.H. PRESS ET AL.,
c              CAMBRIDGE UNIVERSITY PRESS, 1992  ISBN 0 521 43064 X
c  THIS SUBROUTINE IS CALLED FROM PLLCOR and SELCOR.
c  ARGUMENTS:
c   N      = NUMBER OF QUASI-RANDOM NUMBERS
c   XX     = ARRAY CONTAINING THE RANDOM NUMBERS
c  THIS ROUTINE USES 'LOGICAL AND' AND 'EXCLUSIVE OR' SYSTEM FUNCTIONS
c  'IAND' AND 'IEOR' WHICH ARE NON-STANDARD FORTRAN FUNCTIONS !!
c-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION XX(*),FAC
      INTEGER     N,MAXBIT,MAXDIM
      PARAMETER   ( MAXBIT = 30, MAXDIM = 6 )
      INTEGER     I,IM,IK,IPP,J,K,L,IP(MAXDIM),IU(MAXDIM,MAXBIT),
     *            IV(MAXBIT*MAXDIM),IX(MAXDIM),MDEG(MAXDIM)
      EQUIVALENCE (IV,IU)
      SAVE
      DATA        IP /0,1,1,2,1,4/, MDEG /1,2,3,3,4,4/, IX /6*0/,
     *            IV /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
c-----------------------------------------------------------------------
      IF ( N .LT. 0 ) THEN
        DO  K = 1, MAXDIM
          DO  J = 1, MDEG(K)
            IU(K,J) = IU(K,J) * 2**(MAXBIT-J)
          ENDDO
          DO  J = MDEG(K)+1, MAXBIT
            IPP = IP(K)
            I   = IU(K,J-MDEG(K))
c  IEOR IS A NON-STANDARD FORTRAN SYSTEM FUNCTION MAKING 'EXCLUSIVE OR'
            I   = IEOR(I,I/2**MDEG(K))
            DO  L = MDEG(K)-1, 1, -1
c  IAND IS A NON-STANDARD FORTRAN SYSTEM FUNCTION MAKING 'LOGICAL AND'
              IF ( IAND(IPP,1) .NE. 0 ) I = IEOR(I,IU(K,J-L))
              IPP = IPP/2
            ENDDO
            IU(K,J) = I
          ENDDO
        ENDDO
        FAC = 1.D0/(2.D0**MAXBIT)
        IK  = 0
      ELSE
        IM  = IK
        DO  J = 1, MAXBIT
c  IAND IS A NON-STANDARD FORTRAN SYSTEM FUNCTION MAKING 'LOGICAL AND'
          IF ( IAND(IM,1) .EQ. 0 ) GOTO 1
          IM = IM/2
        ENDDO
        WRITE(*,*)'SOBSEQ: MAXBIT =',MAXBIT,' TOO SMALL IN SOBSEQ'
        STOP
    1   CONTINUE
        IM = (J-1) * MAXDIM
        DO  K = 1, MIN( N, MAXDIM )
c  IEOR IS A NON-STANDARD FORTRAN SYSTEM FUNCTION MAKING 'EXCLUSIVE OR'
          IX(K) = IEOR(IX(K),IV(IM+K))
          XX(K) = IX(K) * FAC
        ENDDO
        IK = IK+1
      ENDIF
      RETURN
      END

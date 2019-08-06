c=======================================================================
c
c   s h d a t a r e d u c t i o n . f
c   ---------------------------------
c       check indices of 2-dim array of length 750 m in x-direction
c       and of lenght 1299 m in y-direction on (theoretical) 
c       positions of SD tanks. 
c-----------------------------------------------------------------------
c compilation:
c   gfortran -fbounds-check shdatareduction.f -o shdatareduction
c   f77 -fbounds-check shdatareduction.f -o shdatareduction
c   ifort -C -check bounds shdatareduction.f -o shdatareduction
c execution:
c   ls -1 DAT000071-* | grep t -v | grep n -v > shdatareduction.i000071
c   shdatareduction < shdatareduction.i000071 > shdatareduction.out000071
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;    
c-----------------------------------------------------------------------
c           rows  -22 <= ix <= +22, columns -24 <= iy <= 24
c           distx = 1500., disty = 1299.03810567
c - - - - - - - - - - test print on matching positions outside:
c                 if ( abs(pdata(lp+4))+abs(pdata(lp+5)) .gt. distnear )
c    +            write(7,'(f10.3,1x,1p,3e14.6,0p,3f12.2,1p,e14.6)')
c    +            pdata(lp)*1.e-3,pdata(lp+1),pdata(lp+2),pdata(lp+3),
c    +            pdata(lp+4)*1.e-2,pdata(lp+5)*1.e-2,pdata(lp+lenpar-1)
c-----------------------------------------------------------------------
c#!/bin/bash
c#
c# job_submit -p1 -cp -t40 -m1000 shdatareduction.sh000351
c# 
c# create file list and run `shdatareduction` program:
c# ---------------------------------------------------------------------
c# cd csk000351/
c  ls -1 DAT000351-* | grep t -v | grep n -v > shdatareduction.i000351
c# names of sub paths csk??????;
c# gfortran -fbounds-check shdatareduction.f -o shdatareduction
c# f77 -fbounds-check shdatareduction.f -o shdatareduction
c# ifort -C -check bounds shdatareduction.f -o shdatareduction
c  ./shdatareduction < shdatareduction.i000351 > shdatareduction.out000351
c-----------------------------------------------------------------------
c                                   juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program shdatareduction

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (nfmx=50000,lenthin=6552,lenstnd=5733) 

      character cfilename(nfmx)*120, cfilenout(20)*120, chposit*25

      real pdata(lenthin), rdata(lenthin), qdata(1560), ddata(312)

      dimension xcore(20),ycore(20),rcore(20)

      logical lexist
 
      common /urdat/pdata,rdata,ddata,qdata
      common /utabl/xshcore,yshcore,dradius,qradius,
     + tpoi,nfil,ifil,irec,isho,lrec,lsub,lpoi,lcor

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - initialize constants:
      ! used radius [m] for data around a detector tank.
      dradius = 55.d0 
      qradius = dradius * dradius ! [m^2]

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
      do  i=1,936 ! keep first three subblocks.
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

c - - - - - - open new (shorter) particle data file:
*     write(cfilenout(lcor)(1:lfil+8),'(a,''.rout-'',i2.2)')
*    +      cfilename(nfil)(1:lfil),lcor  
*     write(*,'(20x,a)') cfilenout(lcor)(1:lfil+8)
*     if ( nfil .gt. 1 ) write(*,'(12x,''nfiles ='',i11)') nfil
*     open(unit=4,file=cfilenout(lcor)(1:lfil+8),form='unformatted',
*    +     status='unknown',access='sequential')

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
     +   write(*,'(i9,''. '',a)') ifil,cfilename(ifil)(1:lfil)
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
      endif
c - - - - - - end of corsika particle data file reached: 
  444 continue
      close(unit=1)

c - - - - - - end of loop on known core positions.
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
c     analyze the contents of all 21 subblocks in a record
c
c-----------------------------------------------------------------------
 
      subroutine blwork(iend,lenrec,lp)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)

      real pdata(lenthin), rdata(lenthin), qdata(1560), ddata(312)

      common /urdat/pdata,rdata,ddata,qdata
      common /utabl/xshcore,yshcore,dradius,qradius,
     + tpoi,nfil,ifil,irec,isho,lrec,lsub,lpoi,lcor

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
            isub = 1 + int(lia/lenblk)
            write(*,'(12x,''irec ='',i11,''   isub ='',i3.2)') irec,isub
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
               if ( tpoi .eq. 0. ) tpoi = tpoi + lpoi ! NEW
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
        do  945  lp=lia,lia+lenblk-1,lenpar  
          if ( pdata(lp).ge.75000. .and. pdata(lp).lt. 97000. ) goto 945
          if ( pdata(lp).gt.    0. .and. pdata(lp).lt.196000. ) then
            yp = 1.d-2 * pdata(lp+5) - yshcore ! meter instead of cm.
            myp = int(abs(yp)+dradius)
c - - - - - - - - test on y-coordinate:
            if ( mod(myp,1299) .le. int(dradius*2.) ) then
              xp = 1.d-2 * pdata(lp+4) - xshcore ! meter instead of cm.
              mxp = int(abs(xp)+dradius)
c - - - - - - - - - test on x-coordinate:
              if ( mod(mxp,750) .le. int(dradius*2.) ) then
                ixh = int(1.3337d-3*mxp) * sign(1.d0,xp)
                iyy = int(7.7000d-4*myp) * sign(1.d0,yp)
c - - - - - - - - - - matching tank position for ixh=iyy:
                if ( mod(ixh,2) .eq. mod(iyy,2) ) then
                  pdist = (abs(xp)-750.d0*ixh)*(abs(xp)-750.d0*ixh) +
     +                    (abs(yp)-1299.d0*iyy)*(abs(yp)-1299.d0*iyy)
                  if ( pdist .gt. 1.005432d0*qradius ) goto 945
                  if ( tpoi .eq. 0. ) then
                  if ( qdata(lenblk+152) .gt. 0. ) then
                    write(*,'('' icode'',6x,''px'',12x,''py'',12x,
     +               ''pz'',11x,''x[m]'',6x,''y[m]'',7x,''t[nsec]'',
     +                4x,''weight'',3x,''dist[m]'')')
                  else
                    write(*,'('' icode'',6x,''px'',12x,''py'',12x,
     +               ''pz'',11x,''x[m]'',6x,''y[m]'',7x,''t[nsec]'',
     +                3x,''dist[m]'',
     +                3x,''ixh'',3x,''iyy'',5x,''xp'',4x,
     +                ''mxp'',5x,''yp'',4x,''myp'')')
                  endif
                  endif
                  if ( tpoi .le. 210000. ) then
                   if ( iyy .eq. 1 ) then
                    if ( lcor .ge. 1 .and. pdata(lp+lenpar-1) .ge. 1.)
     +              pdata(lp+lenpar-1) = 1.d-6*pdata(lp+lenpar-1) + lcor
                    icode = int(pdata(lp)*1.d-3)
                    pdata(lp+4) = pdata(lp+4) * 1.d-2
                    pdata(lp+5) = pdata(lp+5) * 1.d-2
                    write(*,'(i5,''.'',1p,3e14.6,0p,2f10.2,f13.7,
     +                f10.2,2i6,f9.2,i5,f9.2,i5)') icode,
     +                (pdata(i),i=lp+1,lp+lenpar-1),sqrt(pdist),
     +                ixh,iyy,xp,mxp,yp,myp
                   endif
                  endif
                  tpoi = tpoi + 1.
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
                    endif
                  endif
              ! elseif ( mod(ixh,2) .ne. mod(iyy,2) ) then ! dummy case.
                  ! write(*,*) ' nothing to do. '
                endif
c - - - - - - - - - - end-of matching tank position for ixh=iyy.
            ! elseif ( mod(mxp,750) .gt. int(dradius) ) then ! dummy case.
                ! write(*,*) ' nothing to do. '
              endif
c - - - - - - - - - end-of test on x-coordinate.
          ! elseif ( mod(myp,1299) .gt. int(dradius) ) then ! dummy case.
              ! write(*,*) ' nothing to do. '
            endif
c - - - - - - - - end-of test on y-coordinate.
          elseif ( pdata(lp) .le. 0. ) then
            ! end of particle data found within this subblock: 
            goto 946 
          elseif ( pdata(lp) .ge. 201000. ) then
            ! write(*,*) '        irec =',irec,'   icod =',pdata(lp)*1.e-3
          endif
  945   continue
      endif
  946 continue
      if ( iend .gt. 0 ) goto 949
  948 continue ! end-of-loop over subblocks (lia=1,lenrec,lenblk).
 
c-----------end of record or end of particle data reached---------------
  949 continue
      return
      end

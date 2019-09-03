c=======================================================================
c
c  r e a d c s k 2 p r t c l s . f
c  ------------------------------- 
c  write new particle data files of selected kind of particles,
c  three additional files, extension `.muons` for codes 5 and 6,
c  `.elect` for codes 2 and 3, and `.gamma` for particle code 1.
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check readcsk2prtcls.f -o readcsk2prtcls 
c     ifort -C -O0 -check bounds readcsk2prtcls.f -o readcsk2prtcls 
c RunProg:
c     ./readcsk2prtcls.sh DAT145216
c     ./readcsk2prtcls < readcsk2prtcls.inp
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c input-files:
c     unit=*: read file name(s):
c           /cr/corsdat1/joe/e16m100prot/DAT045216
c     unit=1,2,3: current corsika particle data file.
c output-files: 
c     unit=*: protocol output.
c     unit=11,12,13: sorted shorter particle data files.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readcsk2prtcls
 
      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (lenthin=6552,lenstnd=5733) 
      character cdata(8000)*250,cdat*250,coutpart*250
      
      double precision aatm(5),batm(5),catm(5)
      double precision qdata(936),parmas(0:101),phead(30)

      real pdata(lenthin),cskbuff(lenthin),cskblock(312)

      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
     +       ,ntypa,ntypb,npart,jpart,nfil,kfil
      common /utabl/pdata,cskblock,parmas,phead,cpi180,c180pi
      common /buffer/cskbuff,lenblk,jflg,jsub,jrec,junit,jseed

c--initialize some quantities-------------------------------------------
      cdat='                                                           '
      lrecmax = 234567890 ! up to 5.38 TBytes.
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      do  i=1,6552
         cskbuff(i) = 0.
      enddo
      do  i=1,312
         cskblock(i) = 0.
      enddo
 
c--read particle data file name-----------------------------------------
      do  ifil=1,8000
      read(*,'(a)',end=402,err=498) cdat
      cdata(ifil) = cdat
      enddo
  402 continue
      nfil = ifil - 1
      inew = nfil + 1
      ifil = 1

c--check and fix length of file name, test simulation-------------------
      isel = 1
      isho = 0
      iend = 0
      irec = 0
      lfil = 250+1
  404 continue
      lfil = lfil - 1
      if ( cdata(1)(lfil:lfil) .eq. ' ' ) goto 404
      ilen = lfil
      idat = index(cdata(1)(1:ilen),'DAT')
      if ( lfil .ne. 9 .and. lfil .ne. 14 .and. lfil .ne. 16 ) then
         if ( index(cdata(1),'_') .eq. 10 ) then
            irun = 987654
            idat = 1
          endif
      else
         irun = 999999
         if ( idat .ge. 1 )
     +      read(cdata(1)(idat+3:idat+8),'(i6)') irun
      endif
      cdata(inew) = cdata(nfil)(idat:ilen)
      write(cdata(inew)(ilen-5:ilen),'(i6.6)') inew
      jout = ilen - idat + 1
      if ( nfil .eq. 1 ) then
         coutpart(1:jout) = cdata(1)(idat:ilen)
      else
         coutpart(1:jout) = cdata(inew)(idat:ilen)
      endif
      ! - - - read one data record with lenstnd words:
      open(unit=isel,file=cdata(1)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(unit=isel,err=496,end=496) (pdata(i),i=1,lenstnd)
      close(unit=isel)
      qdate = dble(pdata(3))
      if ( qdate .lt. 765432.d0 ) then
         qdate = 20000000.d0 + qdate
      else
         qdate = 19000000.d0 + qdate
      endif
      write(*,'(/,13x,f12.3,i8.6,''.'',f12.0,f10.4,f9.4,f12.2)')
     +   pdata(1),int(pdata(2)),qdate,(pdata(i),i=4,6)
      write(*,'(18x,''RUNH'',6x,''runnr'',7x,''date'',5x,''version'',
     +   4x,''nobsl'',4x,''obslev'')')
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
      mthin = 0
      do  i=575,936,7
         if ( int(qdata(i)*1.d-3) .eq. 8888 ) mthin = mthin + 1
      enddo
      ! - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      else if ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      else
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 498
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      jseed = int(pdata(lenblk+13))

c--loop on selecting muons, electrons, gammas--------------------------- 
      do  439  isel=1,3
      if ( isel .eq. 1 ) then
         ntypa = 5
         ntypb = 6
         coutpart(jout+1:jout+6) = '.muons'
      else if ( isel .eq. 2 ) then
         ntypa = 2
         ntypb = 3
         coutpart(jout+1:jout+6) = '.elect'
      else if ( isel .eq. 3 ) then
         ntypa = 1
         ntypb = 1
         coutpart(jout+1:jout+6) = '.gamma'
      endif
      ! open new particle selection file:
      write(*,'(/,8x,a)') coutpart(1:jout+6)
      junit = 10 + isel
      open(unit=junit,file=coutpart(1:jout+6),status='unknown',
     +     form='unformatted',access='sequential')

c--loop on all particle data files from input----------------------------
      do  436  kfil=1,nfil

      isho = 0
      iend = 0
      irec = 0
      npart = 0
      jpart = 0
      jflg = 0
      jrec = 0
      jsub = 0

      open(unit=isel,file=cdata(kfil)(1:ilen),access='sequential',
     +     status='old',form='unformatted')
  431 continue
      ! - - - - - - read data record with lenrec words:
      read(unit=isel,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      if (mod(irec,4000) .eq. 0) write(*,*) '       irec =',irec
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
      close(unit=isel)

  436 continue

      ! - - - - - - end of corsika particle data file reached: 
      do  l=1,lenthin
         pdata(l) = 0.
      enddo
      close(unit=junit)
  439 continue ! end-of loop isel=1,3

      ! - - - - - - print kind of simulation:
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,8x,''`thinning run`'',5x,a,/)')
     +      cdata(nfil)(1:ilen)
      else
         if ( mthin .eq. 0 ) then
            write(*,'(/,8x,''`standard run`'',5x,a,/)')
     +      cdata(nfil)(1:ilen)
         else
            write(*,'(/,8x,''`multithin run`'',5x,a,/)')
     +      cdata(nfil)(1:ilen)
         endif
      endif
      goto 499
 
c--end of data----error exits-------------------------------------------
  496 continue
      write(*,'(27x,''irec ='',i8)') irec
      inull = 0
      do  l=1,lenblk,lenpar
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,lenblk,lenpar*39
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,*)
     +'    ____________________________________________________________'
         write(*,'(5x,
     +      ''ERROR: simulation type of corsika is `thinning`.'')')
         write(*,'(9x,''estimated number of particles'',f15.0,
     +      '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull
         write(*,*) '        subblocks EVTE and RUNE missing!'
      else
         if ( mthin .eq. 0 ) then
            write(*,*)
     +'    ____________________________________________________________'
            write(*,'(5x,
     +         ''ERROR: simulation type of corsika is `standard`.'')')
            write(*,'(9x,''estimated number of particles'',f15.0,
     +         '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull
            write(*,*) '        subblocks EVTE and RUNE missing!'
         else
            write(*,'(8x,14(''_''),/,8x,''`multithin run`'',/)')
            write(*,*)
     +'    ____________________________________________________________'
            write(*,'(5x,
     +         ''ERROR: simulation type of corsika is `multithin`.'')')
            write(*,'(9x,''estimated number of particles'',f15.0,
     +         '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull
            write(*,*) '        subblocks EVTE and RUNE missing!'
         endif  
      endif
      goto 499 
  498 continue
      write(*,*) '    _________________________________________________'
      write(*,*) '    ERROR: missing file name input parameter.'
  499 continue
      stop
      end
c=======================================================================
c
c     block data initialization
c
c-----------------------------------------------------------------------

      block data blinit

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (lenthin=6552)
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      real pdata(lenthin),cskbuff(lenthin),cskblock(312)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
     +       ,ntypa,ntypb,npart,jpart,nfil,kfil
      common /utabl/pdata,cskblock,parmas,phead,cpi180,c180pi
      common /buffer/cskbuff,lenblk,jflg,jsub,jrec,junit,jseed
      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data parmas/
     + 2*0.0000d0 , 0.000511d0 , 0.000511d0 , 0.000000d0 , 0.105658d0,
     + 0.105658d0 , 0.134973d0 , 0.139568d0 , 0.139568d0 , 0.497671d0,
     + 0.493646d0 , 0.493646d0 , 0.939566d0 , 0.938272d0 , 0.938272d0,
     + 0.497671d0 , 0.5488d0   , 1.11563d0  , 1.18937d0  , 1.192550d0,
     + 1.19743d0  , 1.3149d0   , 1.32132d0  , 1.67243d0  , 0.939566d0,
     + 1.11563d0  , 1.18937d0  , 1.19255d0  , 1.19743d0  , 1.3149d0,
     + 1.32132d0  , 1.67243d0  , 1.7841d0   , 1.7841d0   , 1.8693d0,
     + 1.8693d0   , 1.8645d0   , 1.8645d0   , 1.9693d0   , 1.9693d0,
     + 2.2852d0   , 80.6d0     , 80.6d0     , 91.161d0   , 1.8770d0,
     + 2.817d0    , 3.755d0    , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.7669d0   , 0.7681d0   , 0.7681d0   , 1.2309d0   , 1.2323d0,
     + 1.2336d0   , 1.2349d0   , 1.2309d0   , 1.2323d0   , 1.2336d0,
     + 1.2349d0   , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.105658d0,
     + 0.105658d0 , 0.0d0      , 0.0d0      , 0.0d0      , 22*0.0d0/
      end
c=======================================================================
c
c     analyze the contents of all 21 subblocks in a record
c
c-----------------------------------------------------------------------
c
c        primary    energy    runnumb   simnrsh     #hadr    #muons
c        #gammas     #elec   #elecnkg    obslv    theta    phi
c         h1wkm    h1wgcmq  ( emymax    ehasum    ehamax    nha )
c
c-----------------------------------------------------------------------
 
      subroutine blwork(iend,lenrec)

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (lenthin=6552)
      double precision parmas(0:101),phead(30)
      real pdata(lenthin),cskbuff(lenthin),cskblock(312)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
     +       ,ntypa,ntypb,npart,jpart,nfil,kfil
      common /utabl/pdata,cskblock,parmas,phead,cpi180,c180pi
      common /buffer/cskbuff,lenblk,jflg,jsub,jrec,junit,jseed

c-----------initialize counters-----------------------------------------
      isub = 0
      lenblk = lenrec / 21
      lenpar = lenblk / 39

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      isub = isub + 1
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------
         lobs = int(pdata(lia+4))
         if (isho.le.0) then
            phead(6) = lobs
            do  908  i=1,lobs
               phead(6+i) = 1.d-2 * pdata(lia+4+i)
  908       continue
            if (lobs.ge.1.and.lobs.le.10) then
               if (lobs.gt.1) write(*,
     +            '(6x,''observation level:'',f10.2,'' meter'')')
     +            phead(6+lobs)
            else
               write(*,'(6x,''nr of observation level undefined.'')')
               iend = iend + 9876543
               goto 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.d-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
         jflg = 0
         call tobuff(pdata(lia))
      else if ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
         if (iend.lt.0) goto 948
         isho = isho + 1
         phead(1) = pdata(lia+1) ! run number
         phead(2) = 1.d-5 * pdata(lia+6) ! h1st in km.
         pdata(lia+4) = sngl(thickgr(dble(pdata(lia+6))))
         phead(3) = pdata(lia+4) ! h1st in gr. 
         phead(4) = pdata(lia+2) ! primary particle code.
         phead(5) = pdata(lia+3) ! primary energy in GeV.
         phead(19)= pdata(lia+11) ! phi angle in rad.
         phead(20)= pdata(lia+10) ! theta angle in rad.
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = sngl(timev0)
         phead(22) = 1.234567d0
         if ( pdata(lia+148) .ne. 0. ) phead(22) = 1.d0/pdata(lia+148)
         jflg = 0
         call tobuff(pdata(lia))
c----------------subblock longi information (if any)--------------------
      else if (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         jflg = 0
         call tobuff(pdata(lia))
c----------------subblock event end-------------------------------------
      else if ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then
         if ( nfil .gt. 1 ) then
            if ( mod(kfil,20) .eq.1 .or. kfil .eq. nfil ) then
               write(*,*) '       isub =',isub,' (evte)    irec =',irec
     +            ,'    kfil = ',kfil
            endif
         else
            write(*,*) '       isub =',isub,' (evte)    irec =',irec
         endif
         jflg = 0
         if ( jpart .gt. 0 ) call tobuff(cskblock)
         call tobuff(pdata(lia))
c----------------subblock run end---------------------------------------
      else if ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + lenpar ! isho
         jflg = 1
         call tobuff(pdata(lia))
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
            icode = int(1.d-3*(pdata(l)*1.0000001)) ! particle code.
            if ( ntypa .le. icode .and. icode .le. ntypb ) then
               ! copy selected particle to new subblock buffer:
               do  i=1,lenpar
                  cskblock(jpart*lenpar+i) = pdata(l+i-1)
               enddo
               jpart = jpart + 1
               npart = npart + 1
               if ( jpart .eq. 39 ) then
                  ! copy subblock to output current record:
                  jflg = 0
                  call tobuff(cskblock)
                  do  i=1,312 ! lenblk
                     cskblock(i) = 0.
                  enddo
                  jpart = 0
*              else if ( jseed .ge. 6 ) then ! concatenated pll's.
*                 write(*,*) '     jseed = ',jseed,jpart
*                 if ( jpart .lt. 39 ) then
*                    ! copy subblock to output current record:
*                    jflg = 0
*                    call tobuff(cskblock)
*                    do  i=1,312 ! lenblk
*                       cskblock(i) = 0.
*                    enddo
*                    jpart = 0
*                  endif
               endif
            endif   
  917    continue
      endif
      ! if (iend.gt.0) goto 949
  948 continue
 
c-----------all records analyzed or end of showers reached--------------
  949 continue
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
         thickgr = aatm(1) + batm(1) * exp( -arg/catm(1) )
      else if ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg/catm(2) )
      else if ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg/catm(3) )
      else if ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg/catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif
 
      return
      end
c=======================================================================
c
c  tobuff(cskblock)
c
c     copy subblock data to output record and write to binary file
c
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c
c-----------------------------------------------------------------------

      subroutine tobuff(cskblock)

      implicit double precision (a-h,o-z), integer (i-n)
      real cskbuff(6552),cskblock(312)
      integer j, jflg, jsub, jrec, junit, jseed
      common /buffer/cskbuff,lenblk,jflg,jsub,jrec,junit,jseed

c - - - - - - copy new subblock data to next empty buffer subblock:
      if ( 0 .le. jflg .and. jflg .le. 1 ) then
         do  j=1,lenblk
            cskbuff(jsub*lenblk+j) = cskblock(j)
         enddo
         jsub = jsub + 1
      endif

c - - - - - - write next record of new particle data file:
      if ( jflg .eq. 1 .or. jsub .eq. 21 ) then
         jrec = jrec + 1
         write(junit) (cskbuff(j),j=1,lenblk*21)
         do  j=1,6552 ! lenblk*21
            cskbuff(j) = 0.
         enddo
         if ( jsub .eq. 21 ) jsub = 0
      endif
 
c - - - - - - end-of subroutine tobuff.
      return
      end

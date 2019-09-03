c=======================================================================
c
c  r e a d c s k 2 b e o k . f
c  --------------------------- 
c  write corsika shower data as readable ascii file, needs factor 3.5
c  of disk space as of binary particle data file; the protocol output
c  and the ascii particle data file will be written to current
c  directory: 2 records, if 'total_number_of_files'>0, or up to
c  234567890 records (i.e. 5.38 TBytes), if 'total_number_of_files'<0.
c-----------------------------------------------------------------------
c Complink:
c  gfortran -O0 -fbounds-check readcsk2beok.f -o readcsk2beok 
c  ifort -C -O0 -check bounds readcsk2beok.f -o readcsk2beok
c RunProg:
c  ls -1 /lsdf/corsika/e16m100prot-qgs4/th00/DAT?????? > readcsk2beok.names
c  ./readcsk2beok < readcsk2beok.names [ > readcsk2beok.protocol ]
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c  input-files:
c           unit=*: text file with names of particle data files.
c           unit=3: current corsika particle data file.
c  output-files: 
c           unit=*: protocol output to display or text file.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readcsk2beok
 
      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (lenthin=6552,lenstnd=5733) 
      character cdata(10000)*250,cdat*250,cblk*250
      character crunh*4,cevte*4
      double precision aatm(5),batm(5),catm(5)
      double precision qdata(936),parmas(0:101),phead(30)
      real pdata(lenthin),prunh,pevte
      logical lexist
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      data crunh/'RUNH'/,cevte/'EVTE'/,lrecmax/2/,ishift/0/

c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
 
c--read run parameters including file names-----------------------------
      do   lsho=1,2000 
         read(*,'(a)',end=420,err=420) cdat
         cdata(lsho) = cdat
      enddo
  420 continue
      lsho = lsho - 1
      nfil = lsho
      ibeg = 1
 
c----------one or more showers in big disk file-------------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      ! if (ifil.gt.1) close(unit=9)
      isho = 0
      iend = 0
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
      write(*,'(8x,a)') cdata(ifil)(1:ilen)
* - - - - - - read data record with lenstnd words - - - -
      inquire(file=cdata(ifil)(1:ilen),exist=lexist)
      if ( .not. lexist ) then
         write(*,'(29x,''file does not exist!'')')
         goto 444
      endif
      open(unit=3,file=cdata(ifil),status='old',form='unformatted')
      itype = 2
      read(unit=3,err=441,end=440) (pdata(i),i=1,lenstnd)
      close(unit=3)
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
      mthin = 0
      do  i=575,936,7
         if ( int(qdata(i)/1000.) .eq. 8888. ) mthin = mthin + 1
      enddo
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         lenrec = lenstnd
      else if ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 499
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      if ( ishift .eq. -1 ) goto 443 
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
  431 continue
      ! if ( irec .ge. 84040 ) write(*,*) '   next irec =',irec+1
      ! if ( irec .eq. 84040 ) write(*,*) (pdata(lenblk*i+1),i=0,20)
      read(unit=3,err=441,end=440) (pdata(i),i=1,lenrec)
      irec = irec + 1
      if ( pdata(1) .eq . -100000. ) goto 440
      if ( mod(irec,10000) .eq. 0) write(*,*) '       irec =',irec
      if ( i .lt. lenrec .or. pdata(lenrec) .eq. -100000. ) goto 440
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
      goto 443
  440 continue
      write(*,'(22x,''end of data: irec ='',i11,6x,a)')
     +   irec,cdata(ifil)(1:ilen)
      goto 444
  441 continue
      write(*,'(20x,''error in data: irec ='',i11,6x,a)')
     +   irec,cdata(ifil)(1:ilen)
      goto 444
  443 continue
* - - - - - - distinguish contents of pdata vector:
*     if ( qdata(274)+qdata(547) .lt. 1. ) then
*        write(*,'(29x,''thinning '',i8,''. file'')') ifil
*     else
*        if ( mthin .eq. 0 ) then
*           write(*,'(29x,''standard '',i8,''. file'')') ifil
*        else
*           write(*,'(29x,''multithin'',i8,''. file'')') ifil
*        endif  
*     endif
      do  l=1,lenthin
         pdata(l) = -100000.
      enddo
  444 continue
      close(unit=3)
 
c--end of data----------------------------------------------------------
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
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
      common /utabl/pdata,parmas,phead,cpi180,c180pi
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
      real pdata(lenthin)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lrecmax
      common /utabl/pdata,parmas,phead,cpi180,c180pi
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
               phead(6+i) = 1.e-2 * pdata(lia+4+i)
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
         phead(17) = 1.e-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
      else if ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
         if (iend.lt.0) goto 948
         isho = isho + 1
         phead(1) = pdata(lia+1) ! run number
         phead(2) = 1.e-5 * pdata(lia+6) ! h1st in km.
         pdata(lia+4) = sngl(thickgr(dble(pdata(lia+6))))
         phead(3) = pdata(lia+4) ! h1st in gr. 
         phead(4) = pdata(lia+2) ! primary particle code.
         phead(5) = pdata(lia+3) ! primary energy in GeV.
         phead(19) = pdata(lia+11) ! phi angle in rad.
         phead(20) = pdata(lia+10) ! theta angle in rad.
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = sngl(timev0)
c----------------subblock longi information (if any)--------------------
      else if (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
c----------------subblock event end-------------------------------------
      else if ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then
         if (iend.lt.0) then ! ignore first showers.
            iend = iend + 1
            isho = isho + 1
            goto 948
         endif
         write(*,*) '       isub =',isub,' (evte)    irec =',irec,
     +      '    isho =',isho
c----------------subblock run end---------------------------------------
      else if ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         write(9,*) '       isub =',isub-1,' (evte)    !'
* - - - - - - end of corsika particle data file reached: 
         write(9,*) '       isub =',isub,' (rune)    irec =',irec,
     +      '    irwc =',819*irec,'      ',int(phead(18))
         write(*,*) '       isub =',isub,' (rune)    irec =',irec,
     +      '    irwc =',819*irec,'      ',int(phead(18))
         iend = iend + iend + isho
         goto 949
      else
c-----------subblock with particle data---------------------------------
         ! icode = int(1.e-3*(pdata(l)*1.0000001)) ! particle code.
      endif
      if (iend.gt.0) goto 949
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
 
      if      ( arg .lt. 4.d5 ) then
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

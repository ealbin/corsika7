c=======================================================================
c
c  r e a d c s k 2 a s c i . f
c  --------------------------- 
c  write corsika shower data as readable ascii text file; it needs a
c  factor 3.5234 of disk space as of the binary particle data file;
c  the protocol output and the big ascii particle data file will be
c  written to the current directory: 2 records written, if 
c  'total_number_of_files' > 0, or up to maximum 234567890 records 
c  (i.e. 5.38 TBytes !!), if 'total_number_of_files' < 0;
c  ==> increase `lrecmax` data statement to print more records.
c-----------------------------------------------------------------------
c CompLink:
c           gfortran -O0 -fbounds-check readcsk2asci.f -o readcsk2asci 
c           ifort -C -O0 -check bounds readcsk2asci.f -o readcsk2asci
c RunProg:
c           ./readcsk2asci < readcsk2asci.i
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c input-files:
c           unit=*: number of showers and file name(s):
c           ------------------------------------------------
c                       1          'total_number_of_showers'
c                       1          'total_number_of_files'
c            '/lxdata/d2lx14/joe/DAT045216'
c                       1          '_showers_of_this_file_'
c           ------------------------------------------------
c           unit=3: current corsika particle data file.
c output-files: 
c           unit=*: protocol output.
c           unit=9: ascii file named DATnnnnnn.ascithin
c                   or ....ascistnd or ....ascimthi.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readcsk2asci
 
      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (lenthin=6552,lenstnd=5733) 
      character cdata(100)*250,cdat*250,cblk*250
      character cout*250,crunh*4,cevte*4
      double precision aatm(5),batm(5),catm(5)
      double precision qdata(936),parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,lrecmax
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      data lrecmax/234567890/ ! up to 5.38 TBytes.
      data crunh/'RUNH'/,cevte/'EVTE'/

c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
 
c--read run parameters including file names-----------------------------
      read(*,*,end=498,err=498) lsho
      read(*,*,end=498,err=498) nfil
      if ( nfil .eq. 0 ) then
         goto 499
      else if ( nfil .lt. 0 ) then
         nfil = -nfil
         lrecmax = 234567890 ! write all records to asci.   
      endif
      lrecmax = 234567890 ! write all records to asci.   
      ! write(*,'(22x,''total number of showers:'',i7)') lsho
      ! write(*,'(22x,''total number of files  :'',i7)') nfil
      if ( nfil .gt. 100 ) goto 498
      do  ifil=1,nfil
         read(*,*,end=498,err=498) cdat
         read(*,*,end=498,err=498) nfsh
         cdata(ifil) = cdat
      enddo
 
c----------one or more showers in big disk file-------------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      if (ifil.gt.1) close(unit=9)
      ibeg = 1
      isho = 0
      iend = 0
      irec = 0
      idat = index(cdata(ifil),'DAT')
      if ( idat .ge. 1 ) then
         ibeg = idat
      else 
         icer = index(cdata(ifil),'CER') 
         if ( icer .ge. 1 ) then
            ibeg = icer
         else 
            ibeg = 1
         endif
      endif
      ilen = index(cdata(ifil),' ') - 1
      write(*,'(/,8x,a)') cdata(ifil)(1:ilen)
* - - - - - - read data record with lenstnd words - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      ishift = 2
      itype = 2
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenstnd)
      close(unit=3)
      qdate = dble(pdata(3))
      if ( qdate .lt. 765432.d0 ) then
         qdate = 20000000.d0 + qdate
      else
         qdate = 19000000.d0 + qdate
      endif
      write(*,'(/,14x,f11.2,i8.6,''.'',f13.2,
     +   2f9.4,f12.2)') pdata(1),int(pdata(2)),qdate,(pdata(i),i=4,6)
      write(*,'(18x,''RUNH'',6x,''runnr'',7x,''date'',5x,''version'',
     +   4x,''nobsl'',4x,''obslev'')')
* - - - - - - check on reading 32bit simulation on 64bit machine:
      if ( 211285.2 .lt. pdata(1) .and.
     +     pdata(1) .lt. 211285.4 ) then
         ishift = 0
      else if ( 211285.2 .lt. pdata(2) .and.
     +     pdata(2) .lt. 211285.4 ) then
         ishift = -1
         do  i=1,936-1
            pdata(i) = pdata(i+1)
         enddo
      else if ( 2.0202 .lt. pdata(3) .and. pdata(3) .lt. 9.9999 ) then 
         ! check version number instead of testing (273+1) and (312+1).
         ishift = 1
         do  i=936,2,-1
            pdata(i) = pdata(i-1)
         enddo
         pdata(1) = prunh ! 211285.2812500000000000; 
      endif
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
      mthin = 0
      do  i=575,936,7
         if ( int(qdata(i)/1000.) .eq. 8888 ) mthin = mthin + 1
      enddo
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         if ( mthin .eq. 0 ) then
            cout = cdata(ifil)(ibeg:ilen)//'.ascistnd'
         else
            cout = cdata(ifil)(ibeg:ilen)//'.ascimthi'
         endif
         lenrec = lenstnd
      else if ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
         cout = cdata(ifil)(ibeg:ilen)//'.ascithin'
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 498
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      do  i=0,21
         lpdat(i) = i * lenblk
      enddo
      write(*,'(/,8x,a,/)') cout(1:ilen+9) ! name of asci protocol file.
      if ( ishift .eq. -1 ) goto 443 ! should not occur. 
      do  i=ilen+10,200
         cout(i:i) = ' '
      enddo
* - - - - - - read data record with lenrec words - - - -
      open(unit=9,file=cout(1:ilen+9),status='unknown',form='formatted')
      open(unit=3,file=cdata(ifil),status='old',form='unformatted')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
* - - - - - - check on original 64bit simulation:
      pdatone = pdata(1)
      if ( irec .gt. 1 .and. ( 3.21d-41 .lt. pdatone
     +                 .and. pdatone .lt. 3.23d-41 ) ) goto 443
      if (mod(irec,2000) .eq. 0) write(*,*) '       irec =',irec
* - - - - - - - - - - -
      if ( ishift .eq. 1 ) then
* - - - - - - shift data corresp. to a 32bit simul. on a 64bit machine:
         do  i=lenrec,2,-1
            pdata(i) = pdata(i-1)
         enddo
         if ( irec .gt. 1 ) then
            pdata(1) = 177177.1
            ipart = int(pdata(lenpar+1) * 1.000001d-3)
            if ( ipart .ge. 1 .and. ipart .le. 3 )
     +         pdata(1) = 1000. + mod(pdata(lenpar+1),1000.)
            if ( ipart .eq. 0 ) pdata(1) = 1091. ! fixed generation.
            if ( ipart .eq. 5 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 6 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 75 ) pdata(1) = 177177.1
            if ( ipart .eq. 76 ) pdata(1) = 177177.1
            if ( 3301.2 .lt. pdata(lenblk+1) .and.
     +         pdata(lenblk+1) .lt. 3301.4 ) pdata(1) = pevte ! 3397.39185
         else ! irec = 1.
            pdata(1) = prunh ! 211285.2812500000000000;
         endif
      else if ( ishift .eq. -1 ) then
* - - - - - - shift data corresp. to a 64bit simul. on a 32bit machine:
         do  i=1,lenrec-1
            pdata(i) = pdata(i+1)
         enddo
         pdata(lenrec) = -pdata(lenrec-lenpar) 
         ! check negative elements carefully.
      endif
* - - - - - - - - - - -
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
* - - - - - - end of corsika particle data file reached: 
      write(*,*) '       isub =',isub,' (rune)    irec =',irec,
     +      '    irwc =',819*irec
      if ( irec .le. lrecmax .or. lrecmax .ge. 10 ) then 
         if ( 0 .lt. isub .and. isub .lt. 21 ) then
            do  i=1+isub*39,819
               write(9,'(1p,8e14.6)') (0.,l=0,lenpar-1) 
            enddo
         endif
      endif
  443 continue
* - - - - - - distinguish contents of pdata vector:
      if ( ishift .eq. 1 ) then
        write(*,*) '    ______________________________________________'
        write(*,*) '    WARNING: will be better using new 32bit data. '
      else if ( ishift .eq. -1 ) then
        write(*,*) '    irec=1.01',(pdata(i),i=1,6)
        write(*,*) '         1.02',(pdata(i),i=lenblk+1,lenblk+6)
        write(*,*) '         1.03',(pdata(i),i=lenblk*2+1,lenblk*2+5)
        write(*,*) '             ',(pdata(i),i=lenblk*2+lenpar+1,
     +                                         lenblk*2+lenpar+5)
        write(*,*) '    _______________________________________________'
        write(*,*) '    fortran: impossible to read a 64bit simulation;'
        write(*,*) '        use new 32bit runs to get full ascii output'
        write(*,*) '        (or use a program written in C/C++).'
        open(unit=9,file=cout,status='unknown',form='formatted')
        do  l=1,lenblk*3,lenpar
           write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
        enddo
        close(unit=9)
      else if ( ishift .eq. 0 ) then
      !  write(*,*) '    ______________________________________________'
      !  write(*,*) '    WARNING: last record of the ascii output file '
      !  write(*,*) '        may be incompletely filled with zeroes as '
      !  write(*,*) '        compared to the exactly last record of the'
      !  write(*,*) '        corsika particle data file.'
      endif
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,8x,''`thinning run`'',5x,a,/)')
     +      cdata(ifil)(1:ilen)
      else
         if ( mthin .eq. 0 ) then
            write(*,'(/,8x,''`standard run`'',5x,a,/)')
     +      cdata(ifil)(1:ilen)
         else
            write(*,'(/,8x,''`multithin run`'',5x,a,/)')
     +      cdata(ifil)(1:ilen)
         endif  
      endif
      do  l=1,lenthin
         pdata(l) = 0.
      enddo
  444 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '       irec =',irec
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
      write(*,*) '    ERROR: missing some input parameters.'
  499 continue
      stop
      end
c=======================================================================
c
c     block data initialization
c
c=======================================================================

      block data blinit

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)

      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,lrecmax
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
c=======================================================================
c
c        primary    energy    runnumb   simnrsh     #hadr    #muons
c        #gammas     #elec   #elecnkg    obslv    theta    phi
c         h1wkm    h1wgcmq  ( emymax    ehasum    ehamax    nha )
c
c=======================================================================
 
      subroutine blwork(iend,lenrec)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)
 
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,lrecmax
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
               phead(6+i) = 1.d-2 * pdata(lia+4+i)
  908       continue
            if (lobs.ge.1.and.lobs.le.10) then
               if (lobs.gt.1) write(*,
     +            '(6x,''observation level:'',f10.2,'' meter'')')
     +            phead(6+lobs)
            else
               write(*,'(6x,''nr of observation level undefined.'')')
               iend = iend + 987654321
               goto 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.d-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
         do  909  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  909    continue
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
         phead(19) = pdata(lia+11) ! phi angle in rad.
         phead(20) = pdata(lia+10) ! theta angle in rad.
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = sngl(timev0)
         do  911  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  911    continue
         phead(22) = 0.987654d0
         if ( pdata(lia+148) .gt. 0. ) 
     +      phead(22) = 1.d0/pdata(lia+148) ! 1/thinlev
c----------------subblock longi information (if any)--------------------
      else if (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         do  912  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  912    continue
c----------------subblock event end-------------------------------------
      else if ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then
         if (iend.lt.0) then ! ignore first showers.
            iend = iend + 1
            isho = isho + 1
            goto 948
         endif
         do  913  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  913    continue
         write(*,*) '       isub =',isub,' (evte)    irec =',irec,
     +      '    isho =',isho
c----------------subblock run end---------------------------------------
      else if ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + lenpar + isho
         do  914  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  914    continue
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
            icode = int(1.d-3*(pdata(l)*1.0000001)) ! particle code.
            if ( 0. .lt. pdata(l) .and. pdata(l) .lt. 1000. ) then
               if ( pdata(l)*phead(22) .gt. 9.9e6 ) then
                  ! calculate cherenkov bunch size:
                  mbunch = int(dmod(phead(22)*pdata(l)-1.d0,1.d5)/10.d0)
               endif                
            endif
            if ( irec .le. lrecmax )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  917    continue
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

c=======================================================================
c
c  c o n c a t c o r s i k a . f
c  =============================
c create one big particle data file of many small ones of a parallel
c simulation (containing then many `zero` particles, about 20*jfiles).
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c CompLink:
c     gfortran -O0 -fbounds-check concatcorsika.f -o concatcorsika
c     ifort -C -O0 -check bounds concatcorsika.f -o concatcorsika
c RunProgr:
c     ./concatcorsika < concatcorsika.jfiles [ > concatcorsika.jprint ]
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c          RUNH = 211285.2812500000000000;
c          EVTH = 217433.0781250000000000;
c          LONG =  52815.2968750000000000;
c          EVTE =   3397.3918457031250000;
c          RUNE =   3301.3325195312500000;
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c input-files:
c     unit=*: list of particle data files 
c     unit=14: /lxdata/d2lx104/joe/csk106...
c output:
c     unit=24: new corsika particle data file.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------

      program concatcorsika

      implicit double precision (a-h,q-z), real (p), integer (i-n) 

      parameter (lenthin=6552, lenstnd=5733, maxrec=214000000)

      character chfiles(100000)*29, chtotal*16

      dimension pdata(lenthin)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read and initialize run parameters:
      do  ifil=1,100000
         read(*,'(a)',end=402,err=402) chfiles(ifil)
      enddo
  402 continue
      jfiles = ifil - 1
      read(chfiles(1)(4:9),'(i6)') jrunnr
      ! - - - - get length of file names (i.e. 16 or 29):
      ilen = 30
  404 continue
      ilen = ilen - 1
      if ( chfiles(1)(ilen:ilen) .eq. ' ' ) goto 404
      write(*,'(/,9x,''concatcorsika.sh:'')')
      write(*,'(22x,''1st: '',a)') chfiles(1)(1:ilen)
      ! - - - - detect `standard` simulation or `thinning`:
      open(unit=14,file=chfiles(1)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(14,end=429,err=429) (pdata(i),i=1,lenstnd)
      close(unit=14)
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.1 ) then
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.1 ) then
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 429
      endif 
      chtotal = chfiles(1)(1:10)//'999999'
      ! - - - - open concatenated particle data file:
      open(unit=24,file=chtotal,status='unknown',form='unformatted',
     +     access='sequential')
      lenblk = lenrec / 21
      lenpar = lenblk / 39
      nrec = 0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read first record of first particle data file:
      ifil = 1
      open(unit=14,file=chfiles(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(14,end=408,err=408) (pdata(i),i=1,lenrec)
      iend = 0
      do  isub=1,lenrec,lenblk
         if ( 3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4 ) then
            iend = iend + isub
            do  is=isub,lenrec
               pdata(is) = 0.e0
            enddo
            write(24) (pdata(i),i=1,lenrec)
            nrec = nrec + 1
         endif
      enddo
      if ( iend .le. 0 ) then
         write(24) (pdata(i),i=1,lenrec)
         nrec = nrec + 1
         ! - - - - read following records of first particle data file:
         do  406  irec=2,maxrec
            read(14,end=408,err=408) (pdata(i),i=1,lenrec) 
            do  isub=1,lenrec,lenblk
               if ( 3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4 )
     +         then
                  ! - - - - - subblock EVTE found:
                  if ( isub .eq. 1 ) goto 408
                  do  is=isub,lenrec
                     pdata(is) = 0.e0
                  enddo 
               endif
            enddo
            write(24) (pdata(i),i=1,lenrec)
            nrec = nrec + 1
  406    continue
      endif
  408 continue
      close(unit=14)
      write(*,'(12x,''nrec ='',i8,4x,''ifil ='',i6)') nrec,ifil

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read following particle data files:
      do  418  ifil=2,jfiles-1
        open(unit=14,file=chfiles(ifil)(1:ilen),status='old',
     +       form='unformatted',access='sequential')
        ! - - - - first record of following particle data file:
        read(14,end=412,err=412) (pdata(i),i=1,lenrec)
        iend = 0
        do  isub=1,lenrec,lenblk
          if (211285.2.lt.pdata(isub) .and. pdata(isub).lt.211285.3)then
            ! - - - - - subblock RUNH (and EVTH) found:
            do  is=isub,isub+2*lenblk-1
               pdata(is) = 0.e0
            enddo 
          endif
          if ( 3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4 ) then
            iend = iend + isub
            if ( isub .eq. 3 ) goto 412
            do  is=isub,lenrec
              pdata(is) = 0.e0
            enddo
            write(24) (pdata(i),i=1,lenrec)
            nrec = nrec + 1
          endif
        enddo
        ! - - - - file consists of >= 2 records:
        if ( iend .le. 0 ) then
          write(24) (pdata(i),i=1,lenrec)
          nrec = nrec + 1
          ! - - - - following records of current particle data file:
          do  410  irec=2,maxrec
            read(14,end=412,err=412) (pdata(i),i=1,lenrec)
            do  isub=1,lenrec,lenblk
              if (3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4)then
                ! - - - - - subblock EVTE found:
                if ( isub .eq. 1 ) goto 412
                do  is=isub,lenrec
                  pdata(is) = 0.e0
                enddo
              endif
            enddo
            write(24) (pdata(i),i=1,lenrec)
            nrec = nrec + 1
  410     continue
        endif
  412   continue
        close(unit=14)
        if ( mod(ifil,10) .eq. 1 )      
     +    write(*,'(12x,''nrec ='',i8,4x,''ifil ='',i6)') nrec,ifil
  418 continue
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read last particle data file:
      ifil = jfiles
      open(unit=14,file=chfiles(ifil)(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      ! - - - - read first record of last particle data file:
      read(14,end=422,err=422) (pdata(i),i=1,lenrec)
      iend = 0
      do  isub=1,lenrec,lenblk
        if (211285.2.lt.pdata(isub) .and. pdata(isub).lt.211285.3) then
          ! - - - - - subblock RUNH (and EVTH) found:
          do  is=isub,isub+2*lenblk-1
             pdata(is) = 0.e0
          enddo
        endif
        if ( 3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4 ) then
             iend = iend + isub
             do  is=isub,lenrec
                pdata(is) = 0.d0
             enddo
             write(24) (pdata(i),i=1,lenrec)
             nrec = nrec + 1
        endif
      enddo
      if ( iend .le. 0 ) then
        ! - - - - file consists of >= 2 records:
        write(24) (pdata(i),i=1,lenrec)
        nrec = nrec + 1
        ! - - - - following records of last particle data file:
        do  420  irec=2,maxrec
          read(14,end=422,err=422) (pdata(i),i=1,lenrec)
          do  isub=1,lenrec,lenblk
            if (3397.3.lt.pdata(isub) .and. pdata(isub).lt.3397.4) then
               ! - - - - - subblock EVTE found:
               if ( isub .le. 20 ) then
                  write(24) (pdata(i),i=1,lenrec)
                  nrec = nrec + 1
                  goto 422
               else if ( isub .eq. 21 ) then
                  write(24) (pdata(i),i=1,lenrec)
                  nrec = nrec + 1
                  read(14,end=422,err=422) (pdata(i),i=1,lenrec)
                  write(24) (pdata(i),i=1,lenrec)
                  nrec = nrec + 1
                  goto 422
               endif
            endif
          enddo
  420   continue
      endif
  422 continue
      close(unit=14)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - end-of program concatcorsika.
  428 continue
      write(*,'(12x,''nrec ='',i8,4x,''ifil ='',i6)') nrec,ifil
      close(unit=24)
      write(*,'(9x,''concatcorsika completed.'',/)')
  429 continue
      stop
      end

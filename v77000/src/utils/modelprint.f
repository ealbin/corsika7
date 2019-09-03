c=======================================================================
c
c  m o d e l p r i n t . f
c  ----------------------- 
c     display high energy and low energy models as text for the given
c     corsika particle data file or model id.
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check modelprint.f -o modelprint 
c     ifort -C -O0 -check bounds modelprint.f -o modelprint
c RunProg:
c     ./modelprint 
c     < enter particle data file name > 
c     /cr/auger02/joe/showtest/DAT054032
c     < enter model code of showsimulist tabular >
c     11121003
c-----------------------------------------------------------------------
c input-files:
c       unit=*: name of particle data file.
c       unit=3: current corsika particle data file.
c output-files: 
c       unit=*: protocol output.
c-----------------------------------------------------------------------
c         RUNH = 211285.2812500000000000;
c         EVTH = 217433.0781250000000000;
c         LONG =  52815.2968750000000000;
c         EVTE =   3397.3918457031250000;
c         RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program modelprint
 
      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (lenthin=6552,lenstnd=5733) 
      character cintext*120,cdat*120,cblk*120
      character cmodlow(3)*8,cmodhig(0:7)*8,cmodqgs(4)*8,cmodsib(4)*8
      character crunh*4,cevte*4,chartabl(64)*1,cfmtint*8
      character chmodels(0:7,0:7)*8
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      data crunh/'RUNH'/,cevte/'EVTE'/,cfmtint/'(i10)'/
      data cmodqgs/'old     ','01c     ','-II-4   ','-II-4   '/
      data cmodsib/'-1.6    ','-2.1    ','-2.3    ','        '/
      data cmodlow/'gheisha ',' urqmd  ',' fluka  '/
      data cmodhig/'  hdpm  ',' venus  ',' sibyll ','  qgsjet',
     +             '  dpmjet',' nexus  ','  epos  ','        '/
      data chartabl/'1','2','3','4','5','6','7','8','9','0',
     +  'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +  'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
     +  'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +  'n','o','p','q','r','s','t','u','v','w','x','y','z',
     +  ' ',' '/
      data chmodels/
     +    3*'      ','  unix  ','macinto ',3*' ',
     +    '        ','        ',' curved ',5*' ',
     +    '        ','neutrino',6*' ',
     +    '        ',7*'cherenko',
     +    '  hdpm  ',' venus  ',' sibyll ',' qgsjet ',' dpmjet ',
     +               ' nexus  ',' epos   ',' sibyll ',
     +    '        ','gheisha ',' urqmd  ',' fluka  ',4*' ',
     +    '        ','  NKG   ',6*' ',
     +    '        ','  EGS   ',6*' '/

c-----------------------------------------------------------------------
c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      mshcode = 0
      isho = 0
      iend = 0
      irec = 0
      ilet = 0
 
c-----------------------------------------------------------------------
c--read run parameters or file name---------------------------
c - - - - - - - read name of a corsika partcile data file or code:
      read(*,'(a120)',err=498,end=498) cintext
      i = 120 + 1
  410 continue
      i = i - 1
      if ( cintext(i:i) .eq. ' ' ) goto 410
      ilen = i

c-----------------------------------------------------------------------
c--check on file name `DAT` and/or with beginning `/` or model code:
      idat = index(cintext(1:ilen),'DAT')
      if ( idat .gt. 0 ) then
         ! data file name `DAT??????`:
         cdat(1:ilen) = cintext(1:ilen)
         lets = 3
      else if ( cdat(1:1).eq.'/' ) then
         ! data file name begins with `/`:
         write(*,*) ' position of last slash in the file name? '
         i = ilen
  412    continue
         i = i - 1
         if ( cdat(i:i) .ne. '/' ) goto 412
         isla = i ! position of last slash in the file name.
         cdat(1:ilen) = cintext(1:ilen)
         ! test on letters in file name (except `DAT`):
         lets = 0
         do  ib=11,62
            ilet = index(cintext(1:ilen),chartabl(ib))
            if ( ilet .gt. 0 ) lets = lets + 1
         enddo
      else if ( idat .le. 0 ) then
         ! test on letters in file name (except `DAT`):
         cdat(1:ilen) = cintext(1:ilen)
         lets = 0
         do  ib=11,62
            ilet = index(cintext(1:ilen),chartabl(ib))
            if ( ilet .gt. 0 ) lets = lets + 1
         enddo
      endif
      models = 0

c-----------------------------------------------------------------------
c----------one or more showers in big disk file-------------------------
      if ( lets .gt. 0 ) then
      itpa = index(cdat(1:ilen),'.part')
      iusc = index(cdat(1:ilen),'_')
      irec = 0
      write(*,'(1x,6(''_ ''),''modelprint.f'',23('' _''))')
      write(*,'(/,13x,a)') cdat(1:ilen)
      ! - - - read data record with lenstnd words - - - -
      open(unit=3,file=cdat(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(unit=3,err=496,end=496) (pdata(i),i=1,lenstnd)
      close(unit=3)
      ! - - - check on reading 32bit simulation:
      if ( 211285.2 .lt. pdata(1) .and.
     +     pdata(1) .lt. 211285.4 ) then
      elseif ( 211285.2 .lt. pdata(2) .and.
     +     pdata(2) .lt. 211285.4 ) then
         do  i=1,936-1
            pdata(i) = pdata(i+1)
         enddo
      elseif ( 2.0202 .lt. pdata(3) .and. pdata(3) .lt. 9.9999 ) then 
         ! check version number instead of testing (273+1) and (312+1).
         do  i=936,2,-1
            pdata(i) = pdata(i-1)
         enddo
         pdata(1) = prunh ! 211285.2812500000000000;
      endif
      ! - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 499
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      if ( lenblk .eq. 312 ) then
         write(*,'(23x,''`thinning run`'')')
      else if ( lenblk .eq. 273 ) then
         write(*,'(23x,''`standard run`'')')
      endif
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
      models = 0
      do  i=73,80
         models = models + 10**(80-i) * int(qdata(lenblk+i))
      enddo
      write(*,'(16x,''lenblk'',i5)') lenblk
      write(*,'(13x,''modelcode'',i10,/)') models
      modlow = int(pdata(lenblk+75))
      modhig = int(pdata(lenblk+76))
      modsib = int(pdata(lenblk+139))
      modqgs = int(pdata(lenblk+141))
      moddpm = int(pdata(lenblk+143))
      endif ! end-of lets>0.

c-----------------------------------------------------------------------
c--or model code from showsimulist tabular found as input string:
      if ( lets .eq. 0 .or. models .gt. 0 ) then
         if ( models .eq. 0 ) then
            write(cfmtint(3:4),'(i2)') ilen
            read(cintext(1:ilen),cfmtint) mshcode
         else
            mshcode = models
         endif
         ishcode = mshcode
         ilen = int(1.d0+log10(1.d0*ishcode))
         do  l=ilen-1,0,-1
            i = ishcode / 10**l
            if ( l .ne. 4 .and. i .gt. 0 ) then
               write(*,'(19x,''pos.'',i2,4x,''code'',i2,
     +            4x,a,i12)') l+1,i,chmodels(i,l),ishcode
            else if ( l .eq. 4 ) then ! high energy model:
              if ( modhig .eq. 3 ) then
               write(*,'(19x,''pos.'',i2,4x,''code'',i2,
     +            4x,a,i12)') l+1,i,chmodels(i,l),ishcode
               write(*,'(41x,a)') cmodqgs(modqgs)
              else if ( modhig .eq. 2 ) then
               write(*,'(19x,''pos.'',i2,4x,''code'',i2,
     +            4x,a,i12)') l+1,i,chmodels(i,l),ishcode
               write(*,'(41x,a)') cmodsib(modsib)
              else 
               write(*,'(19x,''pos.'',i2,4x,''code'',i2,
     +            4x,a,i12)') l+1,i,chmodels(i,l),ishcode
              endif
            endif
            ishcode = ishcode - i*10**l
         enddo
         write(*,'(2x)')
      endif
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '        irec = 1  pdata(273:275)',(pdata(i),i=273,275)
      goto 499 
  498 continue
      write(*,*) '    ERROR: missing some input parameters.'
  499 continue
c--end of program-------------------------------------------------------
      stop
      end

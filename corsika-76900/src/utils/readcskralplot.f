c=======================================================================
c
c  r e a d c s k r a l p l o t . f
c  ------------------------------- 
c      write corsika shower data as readable ascii file, needs 
c      factor 3.5 of disk space as of binary particle data file.
c      This is impossible for 64bit simulations on all machines.
c      End of ascii print out may be incomplete because of
c      differencies in reading particle data of 32bit simulations
c      on 64bit machines. For 64bit corsika simulations only three
c      subblocks will be printed out (switch to 32bit simulations).
c-----------------------------------------------------------------------
c      gfortran -fbounds-check readcskralplot.f -o readcskralplot 
c      f77 -fbounds-check readcskralplot.f -o readcskralplot
c      ifort -C -check bounds readcskralplot.f -o readcskralplot
c      ------------------------------------------------------------
c           runh=211285.281   evth=217433.078
c           long=52815.2969   evte=3397.39185   rune=3301.33252
c      ------------------------------------------------------------
c      input-files:
c           unit=*: file names:
c           /lxdata/d2lx14/joe/DAT045216
c           unit=3: current corsika particle data file.
c     output-files: 
c           unit=*: protocol and table output.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program readcskralplot

      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cout*120,crunh*4,cevte*4
      character cdata(12000)*29,cdat*29,cblk*120
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      double precision genhisalpa(100),genhismuon(100)
      double precision genhiselec(100),genhisgamm(100)
      double precision genradalpa(100),genradmuon(100)
      double precision genradelec(100),genradgamm(100)
      dimension lpdat(0:21)
      real pdata(lenthin),qdata(936),prunh,pevte

      character chtitle(480)*40
      character chline(0:101),chaxis(0:101),cformdec*32
      double precision qhistos(10020,480),qintvalog(0:101)
      common /histch/ chtitle
      common /histos/ qhistos

      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      common /pgens/genhisalpa,genhisgamm,genhiselec,genhismuon
      common /pgrad/genradalpa,genradgamm,genradelec,genradmuon
      data crunh/'RUNH'/,cevte/'EVTE'/

c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      isho = 0
      iend = 0
      irec = 0
 
c--read run parameters including file names-----------------------------
      do  ifil=1,12345
         read(*,'(a29)',end=427,err=427) cdat(1:29)
         cdata(ifil) = cdat(1:29)
      enddo
  427 continue
      nfiles = ifil - 1
      do i=1,100
         genhisalpa(i) = 0.D0
         genhismuon(i) = 0.D0
         genhiselec(i) = 0.D0
         genhisgamm(i) = 0.D0
         genradalpa(i) = 0.D0
         genradmuon(i) = 0.D0
         genradelec(i) = 0.D0
         genradgamm(i) = 0.D0
      enddo
 
c----------one or more showers in big disk file-------------------------
      do  444  ifil=1,nfiles
      if (ifil.gt.1) close(unit=3)
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
      if ( ilen .le. 0 ) ilen = 29
*     write(*,'(/,8x,''shower file '',a)') cdata(ifil)(1:ilen)
* - - - - - - read data record with lenstnd words - - - -
      open(unit=3,file=cdata(ifil),status='old',form='unformatted')
      itype = 2
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenstnd)
      close(unit=3)
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         cout = cdata(ifil)(1:ilen)//'.ascistnd'
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
         cout = cdata(ifil)(1:ilen)//'.ascithin'
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 499
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      do  i=0,21
         lpdat(i) = i * lenblk
      enddo
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil),status='old',form='unformatted')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      call blwork(iend,lenrec)
      goto 431
  432 continue
  444 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '        irec = 1  pdata(273:275)',(pdata(i),i=273,275)
      write(*,*) '    ______________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      goto 499 
  497 continue
      write(*,*) '        irec = 1  pdata(312:314)',(pdata(i),i=312,314)
      write(*,*) '    ______________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `thinning`.'
      goto 499 
  499 continue
c--end-of program and closing printouts---------------------------------

c-----------------------------------------------------------------------
c        hist vector header:
c     qhistos(1,ih) = ident number of the histogram
c     qhistos(2,ih) = number of bins (1st dim), nbin
c     qhistos(3,ih) = lower bound of 1st dim
c     qhistos(4,ih) = upper bound of 1st dim
c     qhistos(5,ih) = marker for all logarithmic axes
c     qhistos(6,ih) = number of bins (2nd dim), nbi2
c     qhistos(7,ih) = lower bound of 2nd dim
c     qhistos(8,ih) = upper bound of 2nd dim
c     qhistos(9,ih) = -99.8877
c     qhistos(10,ih) = not used and not written to unit 9
c     qhistos(11,ih) = first element of histogram contents
c        hist vector trailer 1-dim:                not written to unit 9
c     qhistos(10+nbin,ih) = last element of histogram contents
c     qhistos(10+nbin+3,ih) = min index of 1st dim with entries
c     qhistos(10+nbin+4,ih) = max index of 1st dim with entries
c     qhistos(10+nbin+5,ih) = minimum of number of entries
c     qhistos(10+nbin+6,ih) = maximum of number of entries
c        hist vector trailer 2-dim:                not written to unit 9
c     qhistos(10+nbin*nbi2,ih) = last element of histogram contents
c     qhistos(10+nbin*nbi2+3,ih) = min index of 1st dim with entries
c     qhistos(10+nbin*nbi2+4,ih) = max index of 1st dim with entries
c     qhistos(10+nbin*nbi2+5,ih) = minimum of number of entries
c     qhistos(10+nbin*nbi2+6,ih) = maximum of number of entries
c     qhistos(10+nbin*nbi2+7,ih) = min index of 2nd dim with entries    
c     qhistos(10+nbin*nbi2+8,ih) = max index of 2nd dim with entries
c-----------------------------------------------------------------------

c - - - - - - book histograms 81 to 88:
      nbin = 100
      xlow = 0.
      xupr = 100.
      ulog = 1.
      do  ihist=81,88  
        qhistos(1,ihist) = 1.d0*ihist
        qhistos(2,ihist) = 1.d0*nbin
        qhistos(3,ihist) = xlow
        qhistos(4,ihist) = xupr
        qhistos(5,ihist) = ulog
        qhistos(6,ihist) = 0.d0
        qhistos(7,ihist) = 0.d0
        qhistos(8,ihist) = 0.d0
        qhistos(9,ihist) = -99.8877d0
        do  i=10,nbin+10+1
           qhistos(i,ihist) = 0.
        enddo
      enddo
      ! fill histogram elements:
      ihist = 80
      chtitle(ihist+1) = 'generation of all particles (1,2,3,5,6)'
      chtitle(ihist+2) = 'generation of photons'
      chtitle(ihist+3) = 'generation of electr/positrons'
      chtitle(ihist+4) = 'generation of all muons'
      chtitle(ihist+5) = 'generation in 1000m of all particles'
      chtitle(ihist+6) = 'generation in 1000m of photons'
      chtitle(ihist+7) = 'generation in 1000m of electr/positrons'
      chtitle(ihist+8) = 'generation in 1000m of all muons'
      do  i=1,nbin
         qhistos(i+10,ihist+1) = genhisalpa(i)
         qhistos(i+10,ihist+2) = genhisgamm(i)
         qhistos(i+10,ihist+3) = genhiselec(i)
         qhistos(i+10,ihist+4) = genhismuon(i)
         qhistos(i+10,ihist+5) = genradalpa(i)
         qhistos(i+10,ihist+6) = genradgamm(i)
         qhistos(i+10,ihist+7) = genradelec(i)
         qhistos(i+10,ihist+8) = genradmuon(i)
      enddo

      do  ihist=81,88
 
c - - - - - - plot header of the histogram:
         nbin = qhistos(2,ihist)
         nbi2 = max(1,int(qhistos(6,ihist)))
         write(*,'(55('' =''))')
         write(*,'(1x,a40)') chtitle(ihist)
         write(19,'(1x,a40)') chtitle(ihist)
         if ( qhistos(6,ihist) .eq. 0. ) then ! 1-dim.
           write(*,'(2(0p,2f8.0,1p,2g12.4))')
     +         (qhistos(i,ihist),i=1,5)
           write(19,'(2(0p,2f8.0,1p,2g12.4),0p,f14.4)')
     +         (qhistos(i,ihist),i=1,5), 0.,-99.,-99., -99.8877
           write(19,'(1p,10e12.4)') (qhistos(i,ihist),i=11,10+nbin)
         else ! 2-dim.
           write(*,'(2(0p,2f8.0,1p,2g12.4))')
     +         (qhistos(i,ihist),i=1,8)
           write(19,'(2(0p,2f8.0,1p,2g12.4),0p,f14.4)')
     +         (qhistos(i,ihist),i=1,8), -99.8877
           write(19,'(1p,10e12.4)') (qhistos(i,ihist),i=11,10+nbin*nbi2)
         endif
c - - - - - - calculate minimum and maximum of the histogram:        
         qhistos(10+nbin*nbi2+3,ihist) = -1.
         qhistos(10+nbin*nbi2+4,ihist) = -1.
         qhmin = qhistos(11,ihist)
         qhmax = qhistos(11,ihist)
         do  ib=11,10+nbin*nbi2
            if ( qhistos(ib,ihist) .gt. qhmax )
     +          qhmax = qhistos(ib,ihist)
            if ( qhistos(ib,ihist) .lt. qhmin )
     +          qhmin = qhistos(ib,ihist)
            if ( qhistos(10+nbin*nbi2+3,ihist) .eq. -1. .and.
     +            qhistos(ib,ihist) .gt. 0. ) then
                qhistos(10+nbin*nbi2+3,ihist) = -10 + ib
            endif
            if ( qhistos(10+nbin*nbi2+4,ihist) .eq. -1. .and.
     +            qhistos(21+nbin*nbi2-ib,ihist) .gt. 0. ) then
                qhistos(10+nbin*nbi2+4,ihist) = 11. + nbin*nbi2 - ib
            endif
         enddo
         qhistos(10+nbin*nbi2+5,ihist) = qhmin
         qhistos(10+nbin*nbi2+6,ihist) = qhmax
         if ( qhmin .eq. qhmax ) goto 509

c - - - - - - - calculate stretching factor (5 or 2):
             qhlog = log10(qhmax)
             if ( qhlog .gt. 0. ) then
                ihlog = int(qhlog)
                qhlog = qhlog - ihlog
                qscal = 1. - ihlog
             else
                ihlog = 1. + int(-qhlog)
                qhlog = qhlog + ihlog
                qscal = 1.*ihlog
             endif
             if ( qhlog .gt. 0.69897d0 ) then
                qfact = 1.
             elseif ( qhlog .gt. 0.30103d0 ) then
                qfact = 2.
             else
                qfact = 5.
             endif
             qscal = 10.d0**qscal
             write(*,'(42x,''hmin='',1p,g12.5,''   hmax='',g12.5,
     +          9x,''(content *'',g9.2,'')'')') qhmin,qhmax,1.d0/qscal

c - - - - - - - plot title axis and scale quantities:
             do  i=1,100
                chaxis(i) = '-'
             enddo
             if ( qfact .eq. 5. ) then
                do  i=10,90,10
                   chaxis(i) = 'o'
                enddo
                chaxis(25) = '+'
                chaxis(50) = '|'
                chaxis(75) = '+'
                chaxis(100) = '|'
             elseif ( qfact .eq. 2. ) then
                do  i=10,90,20
                   chaxis(i) = '+'
                   chaxis(i+10) = '|'
                enddo
             else
                do  i=10,100,10
                   chaxis(i) = '|'
                enddo
             endif
             if ( qfact .eq. 5. ) then
                write(*,'(i8,4i25)') 0,(i,i=5,20,5)
             elseif ( qfact .eq. 2. ) then
                write(*,'(i8,5i20)') 0,(10*i,i=1,5)
             else
                write(*,'(i8,10i10)') 0,(10*i,i=1,10)
             endif
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)

c - - - - - - - plot histogram contents (to the right):
            do  ib=11,10+nbin
               do  is=1,100
                  chline(is) = ' '
               enddo
               if ( qfact .eq. 5. ) then
                  do  is=25,100,25
                     chline(is) = '|'
                  enddo
               elseif ( qfact .eq. 2. ) then
                  do  is=20,100,20
                     chline(is) = '|'
                  enddo
               else
                  do  is=10,100,10
                     chline(is) = '|'
                  enddo
               endif
               is = qscal * qfact * qhistos(ib,ihist)
               if ( is .eq. 0 .and. qhistos(ib,ihist) .gt. 0. ) is = 1
               chline(is) = '*'
               write(*,'(f6.2,'' |'',100a1)')
     +            qhistos(3,ihist)+(qhistos(4,ihist)-qhistos(3,ihist))/
     +            qhistos(2,ihist)*(ib-11),(chline(i),i=1,is)
            enddo

c - - - - - - - plot closing axis and scale quantities:
            write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
            if ( qfact .eq. 5. ) then
               write(*,'(i8,4i25)') 0,(i,i=5,20,5)
            elseif ( qfact .eq. 2. ) then
               write(*,'(i8,5i20)') 0,(10*i,i=1,5)
            else
               write(*,'(i8,10i10)') 0,(10*i,i=1,10)
            endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - - - - plot histogram in logarithmic scale:
            qhmin = 1.d0 ! logar. minimum fixed to 10^0.
            qhlog = log10(qhmax)
            idecads = qhlog + 1.d0
         mdecade = int(100.d0/idecads) ! characters per decade.
         ! - - - - - initialize logarithmic intervals:
         factor = 10.d0**(1.d0/mdecade)
         qintvalog(0) = 1.d0/sqrt(factor)
         do  mdec=1,mdecade
            qintvalog(mdec) = sqrt(factor) * factor**mdec
         enddo
         write(*,'(1x,a40)') chtitle(ihist)
         write(*,'(2(0p,2f8.0,1p,2g12.4))') (qhistos(i,ihist),i=1,5)
         write(*,'(42x,''hmin=10^'',i2.2,10x,''hmax=10^'',i2.2)')
     +      int(log10(qhmin)),idecads
         ! - - - - - mark logarithmic ordinate quantities:
         do  i=1,100
            chaxis(i) = '-'
         enddo
         do  i=mdecade,100,mdecade
            chaxis(i) = '|'
         enddo
         if ( mdecade .le. 20 ) then
            if ( mdecade .gt. 7 ) then
               do  i=mdecade/2,100,mdecade
                  chaxis(i) = '+'
               enddo
            endif
         else ! mdecade > 20: 
            do  i=mdecade/3,100,mdecade
               chaxis(i) = '+'
            enddo 
            do  i=2*mdecade/3,100,mdecade
               chaxis(i) = '+'
            enddo 
         endif
         ! - - - - - plot starting logarithmic ordinate:
         cformdec = '("    10^",i1,10i10)'
         write(cformdec(15:16),'(i2)') idecads
         write(cformdec(18:19),'(i2.2)') mdecade 
         write(*,cformdec) 0,(i,i=1,idecads)
         write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
         ! - - - - - calculate histogram markers:
            do  ib=11,nbin+10
               do  is=1,100
                  chline(is) = ' '
               enddo
               do  is=mdecade,100,mdecade
                  chline(is) = '|'
               enddo
               if ( mdecade .le. 20 ) then
                  if ( mdecade .gt. 7 ) then
                     do  i=mdecade/2,100,mdecade
                        chline(i) = '+'
                     enddo
                  endif
               else ! mdecade > 20: 
                  do  i=mdecade/3,100,mdecade
                     chline(i) = '+'
                  enddo
                  do  i=2*mdecade/3,100,mdecade
                     chline(i) = '+'
                  enddo
               endif
               zlog = log10(qhistos(ib,ihist))
               mlog = zlog ! lowest mark 10^0.
               do  mdec=1,mdecade
                  qvalog = qintvalog(mdec) * 10.d0**mlog
                  if ( qhistos(ib,ihist) .lt. qvalog ) goto 502
               enddo
  502          continue
               is = 0 
               if ( qhistos(ib,ihist) .gt. 1. ) then 
                  is = mlog * mdecade + mdec
               elseif ( qhistos(ib,ihist) .eq. 1. ) then
                  is = 1
               endif
               chline(is) = '*'
               write(*,'(f6.0,'' |'',100a1)')
     +            qhistos(3,ihist)+(qhistos(4,ihist)-qhistos(3,ihist))/
     +            qhistos(2,ihist)*(ib-11),(chline(i),i=1,is)
            enddo
         ! - - - - - plot closing logarithmic ordinate:
         write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
         write(*,cformdec) 0,(i,i=1,idecads)
c - - - - - - - - end-of plot histogram in logarithmic scale. 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      enddo ! end-of loop ihist=81,88.

c--end-of program and closing printouts---------------------------------
  509 continue
      write(*,*) '       genhisalpa '
      write(*,'(5f14.0)') (genhisalpa(i),i=1,100)  
      write(*,*) '       genhisgamm '
      write(*,'(5f14.0)') (genhisgamm(i),i=1,100)  
      write(*,*) '       genhiselec '
      write(*,'(5f14.0)') (genhiselec(i),i=1,100)  
      write(*,*) '       genhismuon '
      write(*,'(5f14.0)') (genhismuon(i),i=1,100)  
      write(*,*) '       genradalpa '
      write(*,'(5f14.0)') (genradalpa(i),i=1,100)  
      write(*,*) '       genradgamm '
      write(*,'(5f14.0)') (genradgamm(i),i=1,100)  
      write(*,*) '       genradelec '
      write(*,'(5f14.0)') (genradelec(i),i=1,100)  
      write(*,*) '       genradmuon '
      write(*,'(5f14.0)') (genradmuon(i),i=1,100)  
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(8x,''thinning run'')')
      else
         write(*,'(8x,''standard run'')')
      endif
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
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
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
 
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      double precision genhisalpa(100),genhismuon(100)
      double precision genhiselec(100),genhisgamm(100)
      double precision genradalpa(100),genradmuon(100)
      double precision genradelec(100),genradgamm(100)
      dimension lpdat(0:21)
      real pdata(lenthin)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      common /pgens/genhisalpa,genhisgamm,genhiselec,genhismuon
      common /pgrad/genradalpa,genradgamm,genradelec,genradmuon
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
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.e-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
      elseif ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
         isho = isho + 1
c- - - - - - - - - - simulated shower number - - - - - - - - - - - - -
         phead(1) = pdata(lia+1)
c - - - - - - - - - - height of first interaction (km)  - - - - - - - -
         phead(2) = 1.e-5 * pdata(lia+6)
c - - - - - - - - - - height of first interaction (grams/cm2)- - - - -
         pdata(lia+4) = thickgr(dble(pdata(lia+6)))
         phead(3) = pdata(lia+4)
c - - - - - - - - - - primary particle code - - - - - - - - - - - - - -
         phead(4) = pdata(lia+2)
c - - - - - - - - - - primary particle energy in gev  - - - - - - - - -
         phead(5) = pdata(lia+3)
c - - - - - - - - - - - phi angle this shower in radian - - - - - - - -
         phead(19) = pdata(lia+11)
c - - - - - - - - - - theta angle this shower in radian - - - - - - - -
         phead(20) = pdata(lia+10)
         timev0 = ( phead(2) * 1.e+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458e0
         pdata(lia+37) = timev0
c----------------subblock event end-------------------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         ! longi subblock ignored.
c----------------subblock event end-------------------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then 
         iend = iend + 1
         goto 949
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
           icode = int(1.e-3*(pdata(l)*1.0000001))
           idgen = ( int(pdata(l)) - icode*1000 ) / 10
           radqua = 1.d-4*pdata(l+4)*pdata(l+4) +
     +              1.d-4*pdata(l+5)*pdata(l+5)
           if ( icode.eq.1 ) then
             genhisgamm(idgen) = genhisgamm(idgen) + 1.d0 
             genhisalpa(idgen) = genhisalpa(idgen) + 1.d0 
           endif
           if ( icode.eq.2 .or. icode.eq.3 ) then 
             genhiselec(idgen) = genhiselec(idgen) + 1.d0 
             genhisalpa(idgen) = genhisalpa(idgen) + 1.d0 
           endif
           if ( icode.eq.5 .or. icode.eq.6 ) then
             genhismuon(idgen) = genhismuon(idgen) + 1.d0
             genhisalpa(idgen) = genhisalpa(idgen) + 1.d0 
           endif
           if ( 810000. .le. radqua .and. radqua .le. 1210000. ) then
             if ( icode.eq.1 ) then
               genradgamm(idgen) = genradgamm(idgen) + 1.d0
               genradalpa(idgen) = genradalpa(idgen) + 1.d0
             endif
             if ( icode.eq.2 .or. icode.eq.3 ) then
               genradelec(idgen) = genradelec(idgen) + 1.d0
               genradalpa(idgen) = genradalpa(idgen) + 1.d0
             endif
             if ( icode.eq.5 .or. icode.eq.6 ) then
               genradmuon(idgen) = genradmuon(idgen) + 1.d0
               genradalpa(idgen) = genradalpa(idgen) + 1.d0
             endif
           endif
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
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg/catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg/catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg/catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif
 
      return
      end

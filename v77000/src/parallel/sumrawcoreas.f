c=======================================================================
c
c  s u m r a w c o r e a s . f
c  ---------------------------
c  sum up antenna data files of all subdirectories `*_coreas/`
c  of a parallel corsika-coreas-simulation by scripts (not by MPI)
c  to one new file for each antenna (KIT_CN, IKP, B.425).
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check sumrawcoreas.f -o sumrawcoreas
c Runprog: being in path cskiiiiii
c     ./sumrawcoreas.sh
c     ./sumrawcoreas > sumrawcoreas.outrunnr
c-----------------------------------------------------------------------
c SIMcoreas.rawinfo:
c            072072    SIMcoreas.rawinfo
c     7056   total number of raw files
c       24   files per subdirectory = n_Antennas
c      882   entries per raw file
c  .../corsika.trunk/run/csk072072
c  SIM072072-893463520-000000201_coreas.bins
c-----------------------------------------------------------------------
c  structure of `raw*.dat` files:
c  4.105e-07   8.6428653514112e-15   2.4548122439862e-15   2.5393380524264e-15
c  4.11e-07   -6.0580621846899e-16   3.482588930098e-15   1.0644610143572e-15
c  4.115e-07   -5.9825251228488e-15   -4.327863108646e-15   -4.6702315909869e-15
c  4.12e-07   1.2313474984855e-16   -4.1473129068949e-16   1.33954102071e-14
c  4.125e-07   -2.6180367333925e-15   -5.1368214099552e-15   -8.9741182699155e-15
c  4.13e-07   -2.5609892790167e-15   -1.6169107445055e-15   2.5831183388761e-15
c  4.135e-07   -8.4475358210681e-16   -2.5303027545107e-15   5.4730484201067e-15
c  4.14e-07   -1.7262725075204e-14   2.2542545341602e-15   3.6295720400874e-15
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program sumrawcoreas

      implicit double precision (a-h,o-z), integer (i-n)

      character crawname(30000)*100, chrawline*100
      character cnewname(3000)*100, charpath*100, crawbins*100
      character chfilerawlst*17, tabkey*1
      dimension qdatraw(0:3,30000,80), lrawname(30000)
      logical lexist
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - some initializations:
      tabkey = char(9) ! editor vi: `ctrl-v-tab`
      tfact = 1.00000003d0
      tstep = 5.d-10*tfact ! seconds

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - get steering quantities from file `SIMcoreas.rawinfo`:
      open(unit=1,file='SIMcoreas.rawinfo',
     +     form='formatted',access='sequential',status='old')
      read(1,*) nrunnr
      read(1,*) nrawtotal
      read(1,*) nrawsubdir 
      read(1,*) nrawdat
      read(1,'(a)') charpath
      read(1,'(a)') crawbins
      close(unit=1)
      nrawfildir = nrawtotal / nrawsubdir
      ! - - - - print steering quantities:
      write(chfilerawlst,'(''SIM'',i6.6,''.raw.lst'')') nrunnr
      write(*,'(12(''_ ''),''sumrawcoreas.f'',20('' _''),/)')
      write(*,'(25x,''run number:'',i7.6)') nrunnr
      write(*,'(10x,''total number of raw files:'',i7)') nrawtotal
      write(*,'(10x,''raw files in SIM*_coreas/:'',i7)') nrawsubdir
      write(*,'(15x,''entries per raw file:'',i7)') nrawdat
      write(*,'(16x,a)') charpath(1:80)
      write(*,'(16x,a)') crawbins(1:80)
      write(*,'(25x,''nrawfildir:'',i7)') nrawfildir

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - read all names of SIM*_coreas/raw*.dat`-files:
      open(unit=2,file=chfilerawlst,
     +     form='formatted',access='sequential',status='old')
      do  ifil=1,nrawtotal
         read(2,'(a)',end=71,err=71) crawname(ifil)
         if ( mod(ifil,20*nrawsubdir) .eq. 1 )
     +     write(*,'(i7,''.  '',a80)') ifil,crawname(ifil)
         ic = 100 + 1 ! max length of crawname.
   70    continue
         ic = ic - 1 
         if ( crawname(ifil)(ic:ic) .eq. ' ' ) goto 70
         lrawname(ifil) = ic
      enddo
   71 continue
      close(unit=2)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - get current value of time step of Antenna files:
      ifil = 1
      open(unit=3,file=crawname(ifil)(1:lrawname(ifil)),
     +       form='formatted',access='sequential',status='old')
      read(3,'(a)') chrawline 
          lraw = len(chrawline)
          do  il=lraw,1,-1
            if ( chrawline(il:il) .ne. ' ' ) goto 72
          enddo
   72     continue
          lraw = il
          do  il=1,lraw
            if ( chrawline(il:il) .eq. tabkey ) chrawline(il:il) = ' '
          enddo
      is = 0 ! first position of testing/converting characters.
      call dtcdbl(chrawline,is,tval,1)
      read(3,'(a)') chrawline 
          lraw = len(chrawline)
          do  il=lraw,1,-1
            if ( chrawline(il:il) .ne. ' ' ) goto 73
          enddo
   73     continue
          lraw = il
          do  il=1,lraw
            if ( chrawline(il:il) .eq. tabkey ) chrawline(il:il) = ' '
          enddo
      is = 0 ! first position of testing/converting characters.
      call dtcdbl(chrawline,is,tstep,1)
      close(unit=3)
      tstep = tstep - tval
      tstep = tstep*tfact
      write(*,'(13x,''time-step in raw files:'',1p,e12.4)') tstep

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - main loop over all `nrawsubdir` Antenna positions:
c - - - - - - nrawtotal  = total number of raw files in all subdirect.
c - - - - - - nrawsubdir = number of raw files in subdirectory = nAnts.
c - - - - - - nrawdat    = entries per raw file.
c - - - - - - nrawfildir = nrawtotal / nrawsubdir
      do  84  nbin=1,nrawsubdir

      do  iv=1,30000
        qdatraw(0,iv,nbin) = 0.d0
        qdatraw(1,iv,nbin) = 0.d0
        qdatraw(2,iv,nbin) = 0.d0
        qdatraw(3,iv,nbin) = 0.d0
      enddo
      write(*,*) '===> nbin ',nbin,' of Antenna positions:'
 
c - - - - - run through all files of the same Antenna position:
      do  80  ifil=nbin,nrawtotal,nrawsubdir
        if ( ifil .lt. 3*nrawsubdir ) 
     +     write(*,*) '          ifil = ',ifil 
        if ( ifil .ge. nrawtotal-nrawsubdir ) 
     +     write(*,*) '          ifil = ',ifil
      ! if ( ifil.eq.nbin ) write(*,*) crawname(ifil)(1:lrawname(ifil))
        open(unit=3,file=crawname(ifil)(1:lrawname(ifil)),
     +       form='formatted',access='sequential',status='old')

        do  76  it=1,nrawdat ! all lines in current raw*.dat. 
          read(3,'(a)',end=78,err=78) chrawline
          lraw = len(chrawline)
          do  il=lraw,1,-1
            if ( chrawline(il:il) .ne. ' ' ) goto 74 
          enddo
   74     continue
          lraw = il
          do  il=1,lraw
            if ( chrawline(il:il) .eq. tabkey ) chrawline(il:il) = ' '
          enddo 
          is = 0 ! first position of testing characters.
          call dtcdbl(chrawline,is,tval,1)
          call dtcdbl(chrawline,is,xval,2)
          call dtcdbl(chrawline,is,yval,3)
          call dtcdbl(chrawline,is,zval,4)
          if ( it .eq. 1 ) then
            tmins = tval*tfact
            tmaxs = tmins + tstep*(nrawdat-1)
            ! write(*,'(2x,''it=1'')')
**          write(*,'(5x,''tmins'',6x,''tmaxs'',6x,''tstep'',3x,
**   +        ''nrawdat'',/,1p,e11.3,e11.3,e11.3,0p,i9)')
**   +          tmins, tmaxs, tstep, nrawdat
          endif
          tval = tmins + tstep*(it-1) 
          if ( it .eq. 1 ) tval = tval/tfact  
          tquot = (tval*tfact-tmins)/tstep
          if ( tval .lt. 0.d0 ) tquot = tquot*(1.000003d0)
          indx = 1 + int(tquot)
          if ( abs(tval) .lt. 3.d-16 ) tval = 0.d0          
*         if ( it .ge. nrawdat-55 ) then 
*         if ( ifil .gt. 1 ) then
*            write(*,'(''@'',1p,e11.4,3e24.15)')
*    +         (qdatraw(i,indx,nbin),i=0,3)
*         endif
*         endif
          qdatraw(0,indx,nbin) = tval
          qdatraw(1,indx,nbin) = qdatraw(1,indx,nbin) + xval
          qdatraw(2,indx,nbin) = qdatraw(2,indx,nbin) + yval
          qdatraw(3,indx,nbin) = qdatraw(3,indx,nbin) + zval
*         if ( it .ge. nrawdat-55 ) then 
*         if ( ifil .gt. 1 ) then
*           write(*,'(1p,e12.4,3e24.15)') tval, xval, yval, zval
*           write(*,'(12x,1p,3e24.15,i10)')
*    +        (qdatraw(i,indx,nbin),i=1,3),ifil
*         endif
*         endif
   76   continue ! end-of loop of all lines in current raw*.dat-file.

   78   continue ! error-exit reading current raw*.dat-file.
        close(unit=3)
   80 continue ! end-of loop over all files of the same Antenna posit.

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - now write back currently summed Antenna data:    
      lfil = ifil - nrawsubdir
      cnewname(nbin) = crawname(lfil)(1:lrawname(lfil)) 
      cnewname(nbin)(21:29) = cnewname(nbin)(11:19)
      write(*,*) cnewname(nbin)(1:lrawname(nbin))
      inquire(file=cnewname(nbin)(1:lrawname(nbin)),exist=lexist)
      if ( lexist ) then
        open(unit=8,file=cnewname(nbin),
     +       form='formatted',access='sequential',status='unknown')
      else
        open(unit=8,file=cnewname(nbin),
     +       form='formatted',access='sequential',status='new')
      endif
      do  ix=1,nrawdat
         tval = tmins + tstep*(ix-1)
         tval = 1.0000007d0 * tval
         if ( abs(tval) .lt. 3.d-16 ) tval = 0.d0          
         qdatraw(0,ix,nbin) = tval
         write(8,'(1p,e12.4,3e24.15)') (qdatraw(i,ix,nbin),i=0,3)  
      enddo
      close(unit=8)
   84 continue ! end-of loop over all Antenna positions.
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - program finished incl. error exit:
   89 continue
      write(*,'(39(''_ ''),/)')
 
      stop
      end
c=======================================================================

      subroutine dtcdbl( chline,is,dval,ik )

c-----------------------------------------------------------------------
c  read double precision parameter from data card character string
c  errors are indicated by writing a '!' in chline(1:1)
c  this subroutine is called from datac.
c   chline = character string of input chline
c   is   = pointer for start of next interpretation of chline
c   dval = double precision variable to be returned
c   ik   = number of argument after keyword
c-----------------------------------------------------------------------

      implicit none
      double precision dval
      integer ic,ie,ik,is,lc
      character cfmtreal*7, chline*(*)
      save

      lc = len(chline)
      do  ic=is+1,lc
         if ( chline(ic:ic) .ne. ' ' ) goto 7
      enddo
    7 continue
      if ( ic .gt. lc .or. chline(ic:ic) .eq. '!'
     +                .or. chline(ic:ic) .eq. ' ' ) then
         dval = 0.d0
         chline(1:1) = '!'
         return
      endif
      is = ic
      do  ic=is+1,lc
         if ( chline(ic:ic).eq.' ' .or. chline(ic:ic).eq.'!' ) goto 8
      enddo
    8 continue
      if ( chline(ic:ic).eq.' ' .or. chline(ic:ic).eq.'!' ) then
         ie = ic-1
      else
         ie = ic
      endif
      cfmtreal = '(f00.0)'
      write(cfmtreal(3:4),'(i2)') ie-is+1
      read(chline(is:ie),cfmtreal,err=9) dval
      is = ie
      return
    9 continue
      write(*,'(10x,''parameter'',i2,'' invalid: '',a)')
     +   ik,chline(is:ie)
      chline(1:1) = '!'
      dval = 0.d0
      is = ie

      return
      end

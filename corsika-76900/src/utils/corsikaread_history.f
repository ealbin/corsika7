C=======================================================================

      PROGRAM CORSIKAREAD_HISTORY

C-----------------------------------------------------------------------
C
C                            CORSIKAREAD_HISTORY  
C
C         ===========================================================
C           READ AND PRINT CORSIKA SHOWER DATA WITH EHISTORY OPTION
C         ===========================================================
C
C     output format for particle output (blocklength = 22932+8 fixed)
C     each block consists of 21 subblocks of 273 words.
C----------------------------------------------------------------------
C     compilation:
C      gfortran -fbounds-check -frecord-marker=4 corsikaread_history.f \
C                              -o corsikaread_history
C      f77 -fbounds-check -m32 corsikaread_history.f \
C                              -o corsikaread_history
C      ifort -C corsikaread_history.f -o corsikaread_history
C----------------------------------------------------------------------
C     How to use this program:
C     1) Execute the executable program defining a file 'output' where 
C        the output is written:
C              ./corsikaread_history  >output
C     2) The program asks: ' CORSIKA data file name:' and you type
C        the name of your PARTICLE DATA FILE by the keybordR 
C        (e.g. DATnnnnnn).
C     3) The file 'output' will contain only the muons with their 
C        extended muon additional information suppressing all
C        other particles arriving at the observation level.
C     4) If the THIN option is used, the LENPAR parameter has to be
C        enlarged from 7 to 8 to correspond with the different data 
C        structure of the PARTICLE DATA FILE in the THIN option.
C     5) This program may be used as a basis for extensions including
C        the production of histogram for plotting different features
C        of muons correlated with their production history.
C
C-----------------------------------------------------------------------
*-- Author :    Ralph Engel   04/09/2002
C-----------------------------------------------------------------------

      IMPLICIT NONE

*KEEP,iolun.
      COMMON /IOLUN/ lundat
      integer        lundat
*KEEP,iopar.
      integer        lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEND.

      real*4         curpar(lenpar)
      character*70   file_dat,filename
      integer        i,iret

      SAVE
C-----------------------------------------------------------------------

C  initialize particle masses
      call PAMAF

C  logical input unit
      lundat = 21

C  open data file
      write(0,'(a,$)') ' CORSIKA data file name: '
      read(5,'(a)') file_dat
      filename = file_dat
      write(6,'(/,2a)') ' Opening file ',filename
      open (lundat, file=filename, status='OLD', form='UNFORMATTED',
     &  access='SEQUENTIAL', err=911)

C  read file header
      call get_run_header(iret)
      if (iret.ne.0) goto 110

C  initialize histogramming
      call analysis(-1,curpar)

C  main loop over showers
      i = 0
 50   continue

C  read shower header
      i = i+1
      call get_shower_header(iret)
      write(6,11)
 11   format(1h ,'ID |gen.count.|  ENERGY  |   P_x   |   P_y    |',
     &           '   P_z    |      x     |      y     |      t',/,
     &  1h ,70x,'|    chi     |      z')
      if (iret.eq.1) then
C  finalize shower histogramming: last call to analysis
        call analysis(-2,curpar)
        write(6,'(/,a,/)') ' CORSIKA_READ: analysis of file completed.'
        close(lundat)
        stop
      elseif (iret.eq.2) then
        write(6,'(/,a,i6,/)')
     &' CORSIKA_READ: error reading shower header of shower:',i
        stop
      endif

C  loop over particles in one shower
 100  continue

        call get_particle(curpar,iret)

        if (iret.eq.0) then
C  valid particle entry
          call analysis(1,curpar)
        elseif (iret.eq.-1) then
C  end of particle list of one shower reached
          call analysis(2,curpar)
          goto 50
        else
          write(6,'(/,a,i4,/)')
     &      ' CORSIKA_READ: error reading next particle.'
          stop
        endif

      goto 100


C  error conditions
 110  continue
      write(6,'(/,a,i4,/)') ' CORSIKA_READ: error reading file ',iret
      stop

 911  continue
      write(6,'(/,2a,/)') ' CORSIKA_READ: error opening file: ',filename
      stop

      END

*-- Author :    Ralph Engel   04/09/2002
C=======================================================================

      SUBROUTINE get_block(datblk,iret)

C-----------------------------------------------------------------------
C  read one complete block from CORSIKA data file
C
C  input:  lundat     logical unit of CORSIKA data file
C
C  output: datblk     data block
C          iret       0  block filled with new data
C                     1  error due to end of file
C                     2  error reading file
C-----------------------------------------------------------------------

      IMPLICIT NONE

*KEEP,iolun.
      COMMON /IOLUN/ lundat
      integer        lundat
*KEEP,iopar.
      integer        lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEND.

      real*4         datblk(lenblk)
      real*4         buffer(lenblk,21)
      integer        iret,i,k,iblock
      data           iblock /22/

      SAVE
C-----------------------------------------------------------------------

      if (iblock.eq.22) then
        read(lundat,err=911,end=900) ((buffer(i,k),i=1,lenblk),k=1,21)
        iblock = 1
      endif

      do i=1,lenblk
        datblk(i) = buffer(i,iblock)
      enddo
      iblock = iblock+1

      iret = 0
      return

 900  continue
      write(6,'(/,a,/)') ' GET_BLOCK: error: reading past end of file.'
      iret = 1
      return

 911  continue
      write(6,'(/,a,/)') ' GET_BLOCK: error reading data file.'
      iret = 2

      return
      END

*-- Author :    Ralph Engel   04/09/2002
C=======================================================================

      SUBROUTINE get_run_header(iret)

C-----------------------------------------------------------------------
C  read initial header of CORSIKA data file
C-----------------------------------------------------------------------

      IMPLICIT NONE

*KEEP,iopar.
      integer        lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEEP,rhead.
      COMMON /RHEAD/ runhd(lenblk)
      real*4         runhd
      character*4    runlb
      equivalence    (runlb,runhd(1))
*KEEP,curblk.
      COMMON /CURBLK/ datblk(lenblk)
      real*4         datblk
      character*4    label
      equivalence    (label,datblk(1))
*KEND.

      integer        iret,k,idate

      SAVE
C-----------------------------------------------------------------------

      call get_block(datblk,iret)
      if (iret.ne.0) return

C  copy header to common block
      do k=1,lenblk
        runhd(k) = datblk(k)
      enddo

      if (runlb.ne.'RUNH') then
        write(6,'(/,a,a4,/)')
     &   ' GET_RUN_HEADER: error: no run header found ',runlb
        iret = 2
        return
      endif

      write(6,'(a)') ' CORSIKA data file header:'
      write(6,'(5x,a,i8)') 'run number:         ',nint(runhd(2))
      idate = nint(runhd(3)) + 20000000
      write(6,'(5x,a,i8)') 'date of start:      ',idate
      write(6,'(5x,a,f8.3)') 'program version:    ',(runhd(4))
      iret = 0

      return
      END

*-- Author :    Ralph Engel   04/09/2002
C=======================================================================

      SUBROUTINE get_shower_header(iret)

C-----------------------------------------------------------------------
C  read header of a single shower in CORSIKA data file
C-----------------------------------------------------------------------

      IMPLICIT NONE

*KEEP,iopar.
      integer        lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEEP,shead.
      COMMON /SHEAD/ shwhd(lenblk)
      real*4         shwhd
      character*4    shwlb
      equivalence    (shwlb,shwhd(1))
*KEEP,curblk.
      COMMON /CURBLK/ datblk(lenblk)
      real*4         datblk
      character*4    label
      equivalence    (label,datblk(1))
*KEEP,RTRAIL.
      COMMON /RTRAIL/ runtr(lenblk)
      real*4         runtr
      character*4    runtb
      equivalence    (runtb,runtr(1))
*KEND.

      integer        iret,k

      SAVE
C-----------------------------------------------------------------------

      iret = 0
      call get_block(datblk,iret)
      if (iret.ne.0) return

      if (label.ne.'EVTH') then
        if (label.eq.'RUNE') then
C  copy trailer to common block
          do k=1,lenblk
            runtr(k) = datblk(k)
          enddo
          iret = 1
          return
        else
          write(6,'(/,a,a4,2i8/)')
     &      ' GET_SHOWER_HEADER: error: no shower header found ',
     &      label,nint(datblk(2)),nint(datblk(3))
          iret = 2
          return
        endif
      endif

C  copy header to common block
      do k=1,lenblk
        shwhd(k) = datblk(k)
      enddo

      write(6,'(a,i6)')
     &  ' GET_SHOWER_HEADER: new shower: ',nint(shwhd(2))

      return
      END

*-- Author :    Ralph Engel   04/09/2002
C=======================================================================

      SUBROUTINE get_particle(curpar,iret)

C-----------------------------------------------------------------------
C  read one complete block from CORSIKA data file
C
C  input:  lundat     logical unit of CORSIKA data file
C
C  output: datblk     data block
C          iret       0  block filled with new data
C                    -1  end of shower
C                     1  end of file
C                     2  error reading file
C-----------------------------------------------------------------------

      IMPLICIT NONE

*KEEP,iopar.
      integer         lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEEP,curblk.
      COMMON /CURBLK/ datblk(lenblk)
      real*4          datblk
      character*4     label
      equivalence     (label,datblk(1))
*KEND.

      real*4          curpar(lenpar)
      integer         iret,k,ipart
      data            ipart /40/

      SAVE
C-----------------------------------------------------------------------

      iret = 0
 100  continue

      if (ipart.gt.39) then

        call get_block(datblk,iret)
C  check for errors
        if (iret.ne.0) then
          write(6,'(a,2i3)')
     &  ' GET_PARTICLE: error reading new block (ipart,iret)',ipart,iret
          return
        endif
        if (label.eq.'EVTE') then
C  end of shower reached (EVTE block currently not used)
          iret = -1
          return
        endif
        ipart = 1

      endif

      do  k=1,lenpar
        curpar(k) = datblk(k+lenpar*(ipart-1))
      enddo

      ipart = ipart+1

C  skip empty entries
      if (nint(curpar(1)).eq.0) then
*       write(6,'(a,i4)') ' GET_PARTICLE: empty particle entry',ipart-1
        goto 100
      endif

      return
      END

*-- Author :    Ralph Engel   26/11/2002
C=======================================================================

      SUBROUTINE analysis(Imode,curpar)

C-----------------------------------------------------------------------
C  simple example for analysis
C-----------------------------------------------------------------------

      IMPLICIT NONE
*KEEP,PAM.
      COMMON /PAM/     PAMA,SIGNUM,RESTMS
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000)
*KEEP,iopar.
      integer          lenblk,lenpar
      parameter (lenpar = 7)
c      parameter (lenpar = 8)              ! for thinning
      parameter (lenblk = 39*lenpar)
*KEND.

      real*4           curpar(lenpar)
      integer          Imode

      double precision xm,Elab,px,py,pz,xx,yy,hh,tt,zz,ehis(1000),gen
      integer          i,ID_par,Ishower,Nmax,igen,igenhad
      logical          lmu,lall

      SAVE
C-----------------------------------------------------------------------

      if (Imode.eq.-1) then

        Ishower = 0
        Nmax = 100
        write(6,'(a,i6)') ' ANALYSIS: initialization of histograms',Nmax
        do i=1,Nmax
          ehis(i) = 0.D0
        enddo
        lmu=.false.     !flag to print muons only
        gen=0d0
        lall=.true. !.false.    !set this to ".true" to get all particles

      else if (Imode.eq.1) then

C  particle four momentum
        ID_par = nint(CURPAR(1))/1000
        igen = mod(nint(abs(CURPAR(1))),1000)
        xm = PAMA(abs(ID_par))   !  Attention, neg. ID_par possible
        px = CURPAR(2)
        py = CURPAR(3)
        pz = CURPAR(4)
        Elab = sqrt(xm**2+px**2+py**2+pz**2)
        xx = CURPAR(5)
        yy = CURPAR(6)
        tt = CURPAR(7)

C  print muon/mother/grandmother entries
c  use the fact that we have 
c  "add info" -> "mother" -> "grandmother" -> "muon"
c  to recover full information on muons

        if ((ID_par.eq.5).or.(ID_par.eq.6).or.lmu
     &                 .or.(ID_par.eq.75).or.(ID_par.eq.76)) then
          if (ID_par.ge.70) then
            write(6,'(a)') 'muon production point entry found:'
            lmu=.true.
            hh = CURPAR(7)
            write(6,'(i4,i10,1p,4e11.3,3e13.5)')ID_par,igen,Elab,px,py
     *                                                  ,pz,xx,yy,hh
         else if (lmu.and.ID_par.lt.0.and.tt.gt.0d0) then
            write(6,'(a)') 'muon mother entry found:'
            igenhad=igen
            hh = CURPAR(7)
            write(6,'(i4,i10,1p,4e11.3,3e13.5)')ID_par,igen,Elab,px,py
     *                                                  ,pz,xx,yy,hh
          elseif (lmu.and.ID_par.lt.0.and.tt.lt.0d0) then
            write(6,'(a)') 'muon grandmother entry found:'
            gen = 100 * CURPAR(5)
            zz = CURPAR(6)
            hh = abs(CURPAR(7))
            igen=mod(igen,100)*1000+exp(float(igen/100))-1
            if(igenhad.gt.500)then
              igen=100000000+igen*1000+igenhad-500
            else
              igen=igen*1000+igenhad
            endif
            igenhad=0
            write(6,'(i4,i10,1p,4e11.3,13x,2e13.5)')ID_par,igen,Elab,px
     *                                                  ,py,pz,zz,hh
          elseif (lmu) then
            write(6,'(a)') 'muon entry found:'
            igen=nint(gen)+igen/10
            gen=0d0
            lmu=.false.
            write(6,'(i4,i10,1p,4e11.3,3e13.5,/)') ID_par,igen,Elab,px
     *                                                  ,py,pz,xx,yy,tt
          endif
        elseif(lall)then
          if (ID_par.lt.0.and.tt.gt.0d0) then
            write(6,'(a)') 'particle mother entry found:'
            hh = CURPAR(7)
            igenhad=igen
            write(6,'(i4,i10,1p,4e11.3,3e13.5)')ID_par,igen,Elab,px,py
     *                                                  ,pz,xx,yy,hh
          elseif (ID_par.lt.0.and.tt.lt.0d0) then
            write(6,'(a)') 'particle grandmother entry found:'
            gen = 100 * CURPAR(5)
            zz = CURPAR(6)
            hh = abs(CURPAR(7))
            igen=igenhad+igen*10000
            igenhad=0
            write(6,'(i4,i10,1p,4e11.3,13x,2e13.5)')ID_par,igen,Elab,px
     *                                                  ,py,pz,zz,hh
          else
            write(6,'(a)') 'particle entry found:'
            igen=nint(gen)+igen/10
            gen=0d0
            write(6,'(i4,i10,1p,4e11.3,3e13.5,/)') ID_par,igen,Elab,px
     *                                                  ,py,pz,xx,yy,tt
          endif
        endif

      else if (Imode.eq.2) then
        Ishower = Ishower+1
        write(6,'(a,i6)') ' ANALYSIS: one complete shower read',Ishower
      else if (Imode.eq.-2) then
        write(6,'(a,i6)') ' ANALYSIS: total number of showers:',Ishower
      endif

      return
      END

*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      SUBROUTINE PAMAF

C-----------------------------------------------------------------------
C  PA(RTICLE) MA(SS) F(ILLING)
C
C  FILLS PARTICLE MASS FOR PARTICLE IP IN ARRAY PAMA
C  RESONANCES AND STRANGE BARYONS INCLUDED
C  PARTICLE MASSES ACCORDING TO GEANT TABLE,
C  TAKEN FROM THE PERIODIC TABLE
C  OR CALCULATED WITH THE MASS FORMULA OF WEIZSAECKER
C-----------------------------------------------------------------------

      IMPLICIT NONE
*KEEP,CONSTA.
C*    COMMON /CONSTA/  PI,PI2,OB3,TB3,ENEPER
C*    DOUBLE PRECISION PI,PI2,OB3,TB3,ENEPER
*KEEP,PAM.
      COMMON /PAM/     PAMA,SIGNUM,RESTMS
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000)
*KEND.

      DOUBLE PRECISION CHARGE(75),MASSES(75)
      DOUBLE PRECISION CHARGE2(100),MASSES2(100)
C*    DOUBLE PRECISION AMUS(59,14),BIND,B1,B2,B3,B4,B5,SS
      INTEGER          IA,IC,IN,IP
C*    INTEGER          I,L
      SAVE
C-----------------------------------------------------------------------

C  MASSES REVISED NOV  2004 BY D. HECK
      DATA MASSES /
     * 0.D0       ,.51099892D-3,.51099892D-3,  0.D0     ,.105658369D0,
     * .105658369D0, .1349766D0, .13957018D0,.13957018D0, 0.497648D0 ,!10
     * 0.493677D0 , 0.493677D0 ,.93956536D0 ,.93827203D0,.93827203D0 ,
     * 0.497648D0 , 0.54775D0  , 1.115683D0 , 1.18937D0 , 1.192642D0 ,!20
     * 1.197449D0 , 1.31483D0  , 1.32131D0  , 1.67245D0 ,.93956536D0 ,
     * 1.115683D0 , 1.18937D0  , 1.192642D0 , 1.197449D0, 1.31483D0  ,!30
     * 1.32131D0  , 1.67245D0  , 0.D0       , 0.D0      , 0.D0       ,
     * 0.D0       , 0.D0       , 0.D0       , 0.D0      , 0.D0       ,!40
     * 0.D0       , 0.D0       , 0.D0       , 0.D0      , 0.D0       ,
     * 0.D0       , 0.D0       , 0.D0       , 0.D0      , 0.78259D0  ,!50
     * 0.7690D0   , 0.7665D0   , 0.7665D0   , 1.2305D0  , 1.2318D0   ,
     * 1.2331D0   , 1.2344D0   , 1.2309D0   , 1.2323D0  , 1.2336D0   ,!60
     * 1.2349D0   , 0.89610D0  , 0.89166D0  , 0.89166D0 , 0.89610D0  ,
     * 0.D0       , 0.D0       , 0.D0       , 0.D0      , 0.D0       ,!70
     * 0.54775D0  , 0.54775D0  , 0.54775D0  , 0.54775D0 , 0.D0       /

      DATA CHARGE /
     *  0.D0,+1.D0,-1.D0, 0.D0,+1.D0,-1.D0, 0.D0,+1.D0,-1.D0, 0.D0,
     * +1.D0,-1.D0, 0.D0,+1.D0,-1.D0, 0.D0, 0.D0, 0.D0,+1.D0, 0.D0,
     * -1.D0, 0.D0,-1.D0,-1.D0, 0.D0, 0.D0,-1.D0, 0.D0,+1.D0, 0.D0,
     * +1.D0,+1.D0,+1.D0,-1.D0,+1.D0,-1.D0, 0.D0, 0.D0,+1.D0,-1.D0,
     * +1.D0,+1.D0,-1.D0, 0.D0,+1.D0,+1.D0,+2.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0,+1.D0,-1.D0,+2.D0,+1.D0, 0.D0,-1.D0,-2.D0,-1.D0, 0.D0,
     * +1.D0, 0.D0,+1.D0,-1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /

C  CHARGE2,  MASSES2 AND DECTME2  RUN FROM PARTICLE CODE 101 TO 200
      DATA MASSES2 /
     * 15*0.D0,
     * 1.8645D0   , 1.8697D0   , 1.8697D0   , 1.8645D0   , 1.9682D0   , !120
     * 1.9682D0   , 2.9804D0   , 2.0067D0   , 2.0100D0   , 2.0100D0   ,
     * 2.0067D0   , 2.1121D0   , 2.1121D0   , 0.0D0      , 3.096916D0 , !130
     * 1.77699D0  , 1.77699D0  , 0.D0       , 0.D0       , 0.D0       ,
     * 0.D0       , 2.28646D0  , 2.4679D0   , 2.4710D0   , 2.45402D0  , !140
     * 2.4529D0   , 2.45376D0  , 2.5757D0   , 2.5780D0   , 2.6975D0   ,
     * 0.D0       , 0.D0       , 0.D0       , 2.28646D0  , 2.4679D0   , !150
     * 2.4710D0   , 2.45402D0  , 2.4529D0   , 2.45376D0  , 2.5757D0   ,
     * 2.5780D0   , 2.6975D0   , 0.D0       , 0.D0       , 0.D0       , !160
     * 2.5184D0   , 2.5175D0   , 2.5180D0   , 0.D0       , 0.D0       ,
     * 5*0.D0     ,                                                     !170
     * 2.5184D0   , 2.5175D0   , 2.5180D0   , 0.D0       , 0.D0       ,
     * 5*0.D0     ,                                                     !180
     * 20*0.D0/

      DATA CHARGE2 /
     * 10*0.D0,
     *  5*0.D0,                       0.D0,+1.D0,-1.D0, 0.D0,+1.D0,   !120
     * -1.D0, 0.D0, 0.D0,+1.D0,-1.D0, 0.D0,+1.D0,-1.D0, 0.D0, 0.D0,
     * -1.D0,+1.D0, 0.D0, 0.D0, 0.D0, 0.D0,+1.D0,+1.D0, 0.D0,+2.D0,   !140
     * +1.D0, 0.D0,+1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,-1.D0,-1.D0,
     *  0.D0,-2.D0,-1.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,   !160
     * +2.D0,+1.D0, 0.D0, 6*0.D0,                            +2.D0,
     * -2.D0,-1.D0, 0.D0, 6*0.D0,                            -2.D0,   !180
     * 20*0.D0/

C  ISOTOPE MASSES CALCULATED FROM: ATOMIC DATA AND NUCL.DATA TABLES 39
C  (1988) 289, (WAPSTRA'S VALUES, CORRECTED FOR ELECTRON MASSES)
C*    DATA ((AMUS(I,L),I=1,59),L=1,7) /
C*   * 1.8756D0,  2.8089D0,                                    57*0.D0,
C*   * 2.8083D0,  3.7273D0,  4.6678D0,  5.6054D0,  6.5454D0,   54*0.D0,
C*   * 2*0.D0  ,  5.6014D0,  6.5337D0,  7.4712D0,  8.4067D0,
C*   *                       9.3471D0, 10.2856D0,              51*0.D0,
C*   * 2*0.D0  ,  6.5341D0,  7.4547D0,  8.3926D0,  9.3253D0,
C*   *                      10.2644D0, 11.2008D0,              51*0.D0,
C*   * 2*0.D0  ,  7.4722D0,  8.3932D0,  9.3243D0, 10.2524D0,
C*   *           11.1886D0, 12.1232D0, 13.0618D0, 13.9986D0,   49*0.D0,
C*   * 2*0.D0  ,  8.4091D0,  9.3274D0, 10.2538D0, 11.1747D0, 12.1093D0,
C*   *           13.0406D0, 13.9790D0, 14.9143D0, 15.8531D0,   48*0.D0,
C*   * 4*0.D0  , 11.1915D0, 12.1110D0, 13.0400D0, 13.9687D0, 14.9057D0,
C*   *           15.8394D0, 16.7761D0, 17.7104D0,              47*0.D0/
C*    DATA ((AMUS(I,L),I=1,59),L=8,14) /
C*   * 4*0.D0, 12.1282D0, 13.0446D0, 13.9709D0, 14.8948D0, 15.8302D0,
C*   *             16.7617D0, 17.6973D0, 18.6293D0, 19.5650D0, 46*0.D0,
C*   * 7*0.D0, 15.8325D0, 16.7629D0, 17.6920D0, 18.6429D0, 19.5564D0,
C*   *             20.4907D0, 21.4227D0, 22.3587D0,            44*0.D0,
C*   * 6*0.D0, 15.8464D0, 16.7668D0, 17.6947D0, 18.6174D0, 19.5502D0,
C*   *  20.4794D0, 21.4137D0, 22.3444D0, 23.2839D0, 24.2138D0, 43*0.D0,
C*   * 8*0.D0, 18.6308D0, 19.5532D0, 20.4817D0, 21.4088D0, 22.3414D0,
C*   *  23.2720D0, 24.2059D0, 25.1387D0, 26.0746D0, 27.0099D0,
C*   *  27.9469D0, 28.8820D0, 29.8173D0, 30.7546D0, 31.6913D0, 36*0.D0,
C*   * 7*0.D0, 18.6410D0, 19.5658D0, 20.4860D0, 21.4124D0, 22.3354D0,
C*   *  23.2676D0, 24.1961D0, 25.1292D0, 26.0602D0, 26.9961D0,
C*   *  27.9291D0, 28.8660D0, 29.7994D0, 30.7376D0,            38*0.D0,
C*   * 9*0.D0, 21.4241D0, 22.3488D0, 23.2714D0, 24.1996D0, 25.1261D0,
C*   *  26.0579D0, 26.9880D0, 27.9218D0, 28.8541D0, 29.7894D0,
C*   *  30.7233D0, 31.6599D0, 32.5944D0, 33.5316D0,            36*0.D0,
C*   * 9*0.D0, 22.3591D0, 23.2836D0, 24.2041D0, 25.1304D0, 26.0527D0,
C*   *  26.9838D0, 27.9128D0, 28.8457D0, 29.7761D0, 30.7111D0,
C*   *  31.6431D0, 32.5803D0, 33.5128D0, 34.4505D0, 35.3837D0, 35*0.D0/
C-----------------------------------------------------------------------

C  GEANT PARTICLES  INCLUDING RHO, K*, AND DELTA
      DO  IP = 1,75
        PAMA  (IP) = MASSES(IP)
        SIGNUM(IP) = CHARGE(IP)
*       DECTIM(IP) = DECTME(IP)
      ENDDO

C  RESET REST OF THE ARRAY
      DO  IP = 76,6000
        PAMA  (IP) = 0.D0
        SIGNUM(IP) = 0.D0
      ENDDO

C  NOW FILL IN CHARMED PARTICLES AND OTHER EXOTICS
      DO  IP = 1, 99
        PAMA  (IP+100) = MASSES2(IP)
        SIGNUM(IP+100) = CHARGE2(IP)
*       DECTIM(IP+100) = DECTME2(IP)
      ENDDO

C  LIGHTEST NUCLEUS IS DEUTERON (IA=2, IC=1)
      DO  IA = 2, 59
        DO  IC = 1, IA
          IN = IA - IC
          IP = IA * 100 + IC
C*        IF ( IC .LE. 14 ) THEN
C  MASSES FROM MASS TABLE FOR ISOTOPES
C*          IF ( IN .EQ. 0 ) THEN
C*            PAMA(IP) = IC * PAMA(14)
C*          ELSE
C*            PAMA(IP) = AMUS(IN,IC)
C*          ENDIF
C  SIMPLE SUM OF PROTON AND NEUTRON MASSES
C*          IF ( PAMA(IP) .EQ. 0.D0 )
C*   *                 PAMA(IP) = IC * PAMA(14) + IN * PAMA(13)
C*        ELSE
C  WEIZSAECKERS MASS FORMULA GIVES BINDING ENERGY IN MEV
C*          B1 = 14.1D0 * IA
C*          B2 = (-13.D0) * IA**TB3
C*          B3 = (-0.595D0) * IC**2 / IA**OB3
C*          B4 = (-19.D0) * (IC-IN)**2 / IA
C*          B5 = 33.5D0 / IA**0.75D0
C*          IF     ( MOD(IC,2) .EQ. 0  .AND.  MOD(IN,2) .EQ. 0 ) THEN
C*            SS =  1.D0
C*          ELSEIF ( MOD(IC,2) .EQ. 1  .AND.  MOD(IN,2) .EQ. 1 ) THEN
C*            SS = -1.D0
C*          ELSE
C*            SS =  0.D0
C*          ENDIF
C*          BIND = (B1 + B2 + B3 + B4 + SS*B5)* 1.D-3
C*          BIND = MAX( 0.D0, BIND )
C*          PAMA(IP) = IN * MASSES(13) + IC * MASSES(14) - BIND
C*        ENDIF

C  DO NOT USE BINDING ENERGY EFFECTS
          PAMA(IP)   = IN * MASSES(13) + IC * MASSES(14)
          RESTMS(IP) = PAMA(IP)

C  NUCLEI ARE ASSUMED TO BE FULLY IONIZED
          SIGNUM(IP) = +IC
        ENDDO
      ENDDO

C  MASSES OF MULTINEUTRON CLUSTERS
      DO  IN = 1,59
        IP = 100 * IN
        PAMA  (IP) = IN * PAMA(13)
        SIGNUM(IP) = 0.D0
      ENDDO
C  REST MASS OF LIGHT NUCLEI (DEUTERIUM, TRITIUM, ALPHA)
      RESTMS(201) =        RESTMS(13) +        RESTMS(14)
      RESTMS(301) = 2.D0 * RESTMS(13) +        RESTMS(14)
      RESTMS(402) = 2.D0 * RESTMS(13) + 2.D0 * RESTMS(14)

      RETURN
      END

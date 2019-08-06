C=======================================================================
C
C  c o r s i k a r e a d . f  (without THINNING)
C           ====================================================
C                 READ  AND  PRINT  CORSIKA  SHOWER  DATA
C           ====================================================
C     Output format for particle output (blocklength = 22932+8 fixed)
C     each block consists of 21 subblocks of 273 words.
C----------------------------------------------------------------------
C     compilation:
C         gfortran -fbounds-check -frecord-marker=4 corsikahisto.f -o corsikahisto
C         f77 -fbounds-check -m32 corsikahisto.f -o corsikahisto
C         ifort -C corsikahisto.f -o corsikahisto
C----------------------------------------------------------------------
C     How to use this program:
C     1) Generate a file 'input' containing the path and name of the 
C        DATnnnnnn file to be analyzed by this program.
C        The name should not contain leading blank.
C        Next line should have an integer from 0 to 999999 for the number
C        if file to read (with name DATnnnnnn-xxxxxx from parallel run)
C     2) Execute this program with the file 'input' as standard input:
C              ./corsikaread <input >output
C     3) The file 'output' will contain a short overview of the 
C        content of the DATnnnnnn file to be analyzed.
C     4) The file DATnnnnnn.histo will contain histograms.
C----------------------------------------------------------------------
C     J.Oehlschlaeger, D. Heck, T. Pierog 16 June 2015 
C=======================================================================
      PROGRAM CORSIKAHISTO
      character cdat*70,hname*70,cname*63,cnum*7,cblk*70
      integer ihis
      parameter(ihis=2)

C--Initalize plots-------------------------------------------------------

      call Initialize
      CBLK='                                                  '
      CDAT=CBLK

C--READ FILE NAME-------------------------------------------------------
 429  CONTINUE
      write(*,*)"File name :"
      READ(*,428,END=440,ERR=439) cname
 428  FORMAT(A)
      write(*,*)"Number of files (0 if standart sim) :"
      READ(*,*) ndat
      if(ndat.gt.999999)stop 'Number of files should be < 1e6 !'


C--OPEN output file-----------------------------------------------------
      name=index(cname,' ')
      hname=cname(1:name-1)//'.histo '
      open(ihis,file=hname,status='unknown')

C--Analyse DAT file-----------------------------------------------------
      cnum='-000000'
      do n=min(ndat,1),ndat
        if(n.gt.0)then
          write(cnum(2:7),'(i6.6)')n
        else
          cnum='       '
        endif
        cdat=cname(1:name-1)//cnum//' '
        call corsread(cdat)
      enddo

C--Write output file-----------------------------------------------------

      call xShower(ihis)
      close(ihis)

      stop

  439 CONTINUE
      WRITE(*,*)'         READ ERROR ON STANDARD INPUT'
      GOTO 429
  440 CONTINUE
      WRITE(*,*)'         READ END ON STANDARD INPUT'
      STOP
      end

C=======================================================================
      subroutine corsread(CDAT)
      implicit double precision (a-h,o-z)
      CHARACTER CHV(5733)*4,CIDENT*4,CDAT*70
      DIMENSION PDATA(5733)
      real PDATA
      EQUIVALENCE (CHV(1),PDATA(1))
      common/obslev/obslvl(10),eprima,tfront(10)
      integer maxo,maxmm,musmm,iimom,i1mom,i2mom,i3mom
      parameter(maxo=2,maxmm=maxo)
      common/cxmom1/musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx,2)
      common/countgen/cntgen
      double precision zamin,zamax,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,numie
      common/cxtime/ tamin,tamax,ctime
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden

      parameter(maxiz=10,maxjz=maxiz,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz,maxie)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)
      

      IREC=0
      tmin=1d10
      tmax=0d0

      WRITE(*,430) CDAT
 430  FORMAT(1H ,'READ DATA FROM FILE = ',A)
      OPEN(UNIT=3,FILE=CDAT,STATUS='OLD',FORM='UNFORMATTED')
* - - - - - - read data records with 5733 words - - - -
  431 CONTINUE
      IREC = IREC + 1
      READ(UNIT=3,ERR=434,END=433) PDATA
      if ( mod(irec,100) .eq. 0 ) 
     +   WRITE(*,*)'         HAVE READ RECORD NR.',IREC
C-----------loop over subblocks-----------------------------------------
      DO    LIA=1,5733,273
        CIDENT(1:1) = CHV(LIA)(1:1)
        CIDENT(2:2) = CHV(LIA)(2:2)
        CIDENT(3:3) = CHV(LIA)(3:3)
        CIDENT(4:4) = CHV(LIA)(4:4)
        IF (PDATA(LIA).GE.211284.0.AND.
     +      PDATA(LIA).LE.211286.0) THEN
          CIDENT = 'RUNH'
          WRITE(*,*)'RUNH'
        ENDIF
        IF (PDATA(LIA).GE.217432.0.AND.
     +      PDATA(LIA).LE.217434.0) THEN
          CIDENT = 'EVTH'
          WRITE(*,*)'EVTH'
        ENDIF
        IF (PDATA(LIA).GE. 52814.0.AND.
     +      PDATA(LIA).LE. 52816.0) THEN
          CIDENT = 'LONG'
          WRITE(*,*)'LONG'
        ENDIF
        IF (PDATA(LIA).GE.  3396.0.AND.
     +      PDATA(LIA).LE.  3398.0) THEN
          CIDENT = 'EVTE'
          WRITE(*,*)'EVTE'
        ENDIF
        IF (PDATA(LIA).GE.  3300.0.AND.
     +      PDATA(LIA).LE.  3302.0) THEN
          CIDENT = 'RUNE'
          WRITE(*,*)'RUNE'
        ENDIF
C-----------which kind of block is it?----------------------------------
        IF ( CIDENT.EQ.'RUNH' .OR. CIDENT.EQ.'RUNE' .OR. 
     +       CIDENT.EQ.'LONG' .OR. CIDENT.EQ.'EVTH' .OR. 
     +                             CIDENT.EQ.'EVTE' ) THEN
          CHV(LIA) = CIDENT
          IF     ( CIDENT .EQ. 'RUNH' ) THEN
C----------------subblock run header------------------------------------
             PDATA(LIA) = 11111111.
c             DO    IL=LIA,LIA+272,7
c                WRITE(7,'(1P,7E13.5)') (PDATA(II+IL),II=0,6)
c             ENDDO
             numiz=nint(PDATA(LIA+4))
             print *,numiz,' observation level'
             do IL=0,numiz-1
               obslvl(IL+1)=PDATA(LIA+5+IL)/100d0     !height in m
               print *,IL+1,' at ',obslvl(IL+1),' m'
               zamin=min(zamin,obslvl(IL+1))
               zamax=max(zamax,obslvl(IL+1))
             enddo
             zamin=zamin*0.9d0
             zamax=zamax*1.1d0
             hatm=PDATA(LIA+254)/100d0     !border of atmosphere in m
          ELSEIF ( CIDENT .EQ. 'EVTH' ) THEN
C----------------subblock event header----------------------------------
             PDATA(LIA) = 33333333.
c             DO    IL=LIA,LIA+272,7
c                WRITE(7,'(1P,7E13.5)') (PDATA(II+IL),II=0,6)
c             ENDDO
             eprima=PDATA(LIA+3)
             costh=cos(PDATA(LIA+10))
             height0=PDATA(LIA+6)/100d0
c             height0=hatm
             if(height0.lt.0d0)height0=hatm
             do iz=1,numiz
               if(abs(costh).gt.0.0001)then
                 tfront(iz)=(height0-obslvl(iz))/0.299792458D0
               else
                 tfront(iz)=1d30
               endif
             enddo
C----------------subblock longitudinal data-----------------------------
          ELSEIF ( CIDENT .EQ. 'LONG' ) THEN
             PDATA(LIA) = 55555555.
c             DO    IL=LIA,LIA+272,7
c                WRITE(7,'(1P,7E13.5)') (PDATA(II+IL),II=0,6)
c             ENDDO
C----------------subblock event end-------------------------------------
          ELSEIF ( CIDENT .EQ. 'EVTE' ) THEN
             PDATA(LIA) = 77777777.
c             DO    IL=LIA,LIA+272,7
c                WRITE(7,'(1P,7E13.5)') (PDATA(II+IL),II=0,6)
c             ENDDO
C----------------subblock run end---------------------------------------
          ELSEIF ( CIDENT .EQ. 'RUNE' ) THEN
             PDATA(LIA) = 99999999.
c             DO    IL=LIA,LIA+272,7
c                WRITE(7,'(1P,7E13.5)') (PDATA(II+IL),II=0,6)
c             ENDDO
             GOTO 929
          ENDIF
        ELSE
C-----------subblock with particle data---------------------------------
           DO    IL=LIA,LIA+272,7
             idc=nint(PDATA(IL))
             if(idc.ne.8888000.and.idc.gt.0)then
               id=idc/1000
               if(id.lt.70)then
                 gen=dble(mod(idc,1000)/10)
                 if(gen.lt.99)then
                   if(gen.gt.50d0)then
                     gen=gen-50d0
                   elseif(gen.gt.30d0)then
                     gen=gen-30d0
                   endif
                 endif
                 lvl=mod(idc,10)
                 h1=obslvl(lvl)
                 px=PDATA(IL+1) !px momentum GeV/c
                 py=PDATA(IL+2) !py momentum GeV/c
                 pz=PDATA(IL+3) !pz momentum GeV/c
                 x1=PDATA(IL+4)/100d0 !x in m
                 y1=PDATA(IL+5)/100d0 !y in m
                 t1=PDATA(IL+6) !time in ns
                 wt=1d0         !weight
                 iimode=0       !particle
                 idi=idtrafo('cor','nxs',id)
                 call idmass(idi,am)
                 E1=sqrt(px*px+py*py+pz*pz+am*am)
                 tmin=min(tmin,t1)
                 tmax=max(tmax,t1)
c                 if(tfront(lvl).lt.0d0)
c     .           tfront(lvl)=0.25d0*(t1-tfront(lvl))
               else
                 iimode=1       !additional info
                 idi=0
                 E1=0
               endif
             elseif(idc.gt.0)then
               iimode=2      !multithin
             else
               iimode=-1      !no particle
             endif
             call cana(h1,x1,y1,t1,E1,px,py,pz,am,wt,gen,idi,iimode)
           ENDDO
         ENDIF
      ENDDO
  929 CONTINUE
      GOTO 431
 
C--END OF TEST----------------------------------------------------------
  433 CONTINUE
      WRITE(*,*)'         LAST RECORD ',irec-1
      CLOSE(UNIT=3)
      print *,"tmin",tmin,tmax,tfront(1),tmin-tfront(1)
      return
  434 CONTINUE
      WRITE(*,*)'         READ ERROR ON UNIT 3'
*     CLOSE(UNIT=3)
      GOTO 431
      END


c----------------------------------------------------------------------
      subroutine cana(h1,x1,y1,t1,E1,px,py,pz,am,wt,gen,id,iimode)
c-----------------------------------------------------------------------
c previous version with all moments
c     analyzes a particle arriving at (h1,x1,y1,t1) 
c     id                 =  nexus/isajet particle id
c     x1,y1  =  position in meter (x,y in the obs frame)
c     h1                 =  height (m)
c     t1                 =  time in ns
c     E1                 =  total energy in GeV,GeV/c^2
c     px,py,pz,am        =  momentum,mass in GeV,GeV/c^2 at the end in shower frame
c     wt                 =  weight
c     gen                =  generation
c     iimode             =  particle at gound (0) or additional info (1)
c
c     in case of abs(id)<=1 we have EGS4 which means:
c         0 = photon
c        -1 = electron
c         1 = positron
c        px,py,pz are the directional vectors
c        (not necessarily normalized to 1, but almost)
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      integer maxo,maxmm,musmm,iimom,i1mom,i2mom,i3mom
      parameter(maxo=2,maxmm=maxo)
      common/cxmom1/musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx,2)
      common/countgen/cntgen
      double precision zamin,zamax,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,numie
      common/cxtime/ tamin,tamax,ctime
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden

      parameter(maxiz=10,maxjz=maxiz,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz,maxie)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)
      common/obslev/obslvl(10),eprima,tfront(10)

      if(iimode.ne.0)return
      ee=E1
      amm=am
      ida=abs(id)
      if(ida.le.12)then    ! if (ida.le.1) means EGS4 
        k=1
      else
        k=2
      endif
      ee=ee-amm                 !  use kinetic energy !


c        added----------vc300304---from cana
      in=0
      inti=0
      if(id.eq.10)then
        in=1  !photon
      elseif(id.eq.12)then
        in=2  !electron
        inti=-1
      elseif(id.eq. -12)then
        in=3  !positron
        inti=1
      elseif(ida.eq.1120) then
        in=4  !nucleon
        inti=sign(1,id)
      elseif(ida.eq.1220) then
        in=4  !nucleon
      elseif(ida.eq.120) then
        in=5  !pi+pi-
        inti=sign(1,id)
      elseif(ida.eq.130) then
        in=7  !K+K-
        inti=sign(1,id)
      elseif(id.eq.-20)  then
        in=8  !Klong
      elseif(id.eq. 20)  then
        in=9  !Kshort
      elseif(ida.eq.14)  then
        inti=-sign(1,id)
        if(inti.gt.0)then
          in=11     !muon +
        else
          in=12     !muon -
        endif
      elseif(id.eq.-10)then
        in=6  !gamma from hadron
      endif
c        added----------vc300304---from cana

      if(in.eq.0)return

      igen=int(gen)
      if(igen.le.ngenmx)then
        cntgen(0,k)=cntgen(0,k)+1d0
        cntgen(igen,k)=cntgen(igen,k)+1d0
      endif

      iz=1
      do while(.not.abs(h1-obslvl(iz)).lt.0.0001d0)
        iz=iz+1
      enddo
      yieldz(in,iz)=yieldz(in,iz)+wt
      xx=x1
      yy=y1
      tt=t1-tfront(iz)

      it=max(1,int(log(tt/tamin)/log(ctime)+1.5d0))
      yiex(in,iz,it)=yiex(in,iz,it)+wt

      if(tt.lt.0)print *,'tt,it,ctime',tt,it,ctime

      if(ee-eamin(k).lt.-1d-8.or.ee-eamax(k).gt.1d-8)return
c eamin and eamax are the middle of the first and the last bin
      cdec=(eamax(k)/eamin(k))**(1.d0/dble(numie))


c        spectra for selected heightes
      ie=max(1,int(log(ee/eamin(k))/log(cdec)+1.5d0))
      spec(0,in,ie,iz)=spec(0,in,ie,iz)+wt



c  spectra for selected heightes

      pt2=px**2+py**2
      pp2=pt2+pz**2
      sin2t=pt2/pp2
      pp=sqrt(pp2)
      cost=pz/pp
      aa=px/pp
      bb=py/pp

c        write(6,*) 'pp px py pz aa bb sin2t cost= ',
c     &  pp, px, py, pz, aa, bb, sin2t, cost
      
      r2=xx**2+yy**2
      r=sqrt(r2)                !radial distance (m)
      axby=aa*xx+bb*yy
      h=h1    !so240903
c          fct=rhoair(h)/radlth                           !so020703
      rr1=sin2t
      rr2=r2                    !*fct*fct             !radial distance squared (unit of radiation lenght squared)
      rr3=axby                  !*fct

      do ii=0,musmm

        spec(ii,in,ie,iz)=spec(ii,in,ie,iz)
     &      +wt *rr1**i1mom(ii) *rr2**i2mom(ii) *rr3**i3mom(ii)


c        write(6,*) 'rr1, rr2, rr3 = ',
c     &  rr1, rr2, rr3, ii, i1mom(ii), i2mom(ii),i3mom(ii)
        
        spex(ii,in,ie,iz)=spex(ii,in,ie,iz)
     & +wt*cost *rr1**i1mom(ii) *rr2**i2mom(ii) *rr3**i3mom(ii)
      enddo                     !moments
      if(in.gt.10)then
        do ii=0,musmm
          spec(ii,10,ie,iz)=spec(ii,10,ie,iz)
     &      +wt *rr1**i1mom(ii) *rr2**i2mom(ii) *rr3**i3mom(ii)
          spex(ii,10,ie,iz)=spex(ii,10,ie,iz)
     & +wt*cost *rr1**i1mom(ii) *rr2**i2mom(ii) *rr3**i3mom(ii)
        enddo                   !moments
      endif

      if(r.gt.ramin.and.r.lt.ramax)then
        ir=int(log(r/ramin)/log(ramax/ramin)*numir)+1
        yieldr1(in,iz,ir)=yieldr1(in,iz,ir)+wt

        if(in.gt.10)then
          yieldr1(10,iz,ir)=yieldr1(10,iz,ir)+wt
        endif


        jr=1+(ir-irfirst)/modr
        if(jr.ge.1.and.jr.le.maxjr)then
          specr(in,ie,iz,jr)=specr(in,ie,iz,jr)+wt
          if(in.gt.10)then
            specr(10,ie,iz,jr)=specr(10,ie,iz,jr)+wt
          endif
          do i=1,4
            xx=-1.d100
            if(i.eq.1.and.pt2.gt.0.d0)then
c    cos(phi) distribution for a given height and radius
                 xx=px/pt2
               elseif(i.eq.2)then
c     sin^2(theta) distribution for a given height and radius
                 xx=sin2t
               elseif(i.eq.3)then
c     time distribution for a given height and radius
                 xx=tt
               elseif(i.eq.4)then
c     angle difference between particle position and particle direction
                 xx=atan2(yy,xx)-atan2(py,px)
               endif
               if(xx.gt.xamin(i).and.xx.lt.xamax(i))then
              ix=int((xx-xamin(i))/(xamax(i)-xamin(i))*dble(numix(i)))+1
               yieldx(i,in,iz,jr,ix)=yieldx(i,in,iz,jr,ix)+wt
               if(in.gt.10)then
               yieldx(i,10,iz,jr,ix)=yieldx(i,10,iz,jr,ix)+wt
               endif
               endif
             enddo
           endif
         endif

         if(mod(ie-iefirst,moden).eq.0.and.ie.ge.iefirst
     &        .and.1+(ie-iefirst)/moden.le.maxjr)then
           je=1+(ie-iefirst)/moden
           i=5
c          sin^2(theta) distribution for a given height and energy
           xx=log10(sin2t)
           if(xx.gt.xamin(i).and.xx.lt.xamax(i))then
             speca(in,iz,je)=speca(in,iz,je)+wt
             if(in.gt.10)then
               speca(10,iz,je)=speca(10,iz,je)+wt
             endif
             ix=int((xx-xamin(i))/(xamax(i)-xamin(i))*dble(numix(i)))+1
             yieldx(i,in,iz,je,ix)=yieldx(i,in,iz,je,ix)+wt
             if(in.gt.10)then
               yieldx(i,10,iz,je,ix)=yieldx(i,10,iz,je,ix)+wt
             endif
           endif
         endif

      end


c-------------------------------------------------------------------------
      subroutine xShower(ifho)
c-------------------------------------------------------------------------
c  plot particle type k1 to k2
c-------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer maxo,maxmm,musmm,iimom,i1mom,i2mom,i3mom
      parameter(maxo=2,maxmm=maxo)
      common/cxmom1/musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx,2)
      common/countgen/cntgen
      double precision zamin,zamax,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,numie
      common/cxtime/ tamin,tamax,ctime
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden
      parameter(maxiz=10,maxjz=maxiz,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz,maxie)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)
      common/obslev/obslvl(10),eprima,tfront(10)
      dimension stot(maxjz)
      integer kmn(2),kmx(2)
      character*6 c(maxjz)
      character*11 ptyp(maxin)
      character*3 nyx
      character*17 tx(maxjr,maxjz),ty(maxjr,maxjz,maxin)
      dimension emoy(maxjz),sum(maxjz)



      k1=1
      k2=12
      if(k1.gt.k2)then
        write(*,*)'error in xShower k1 < k2 !'
        return
      elseif(k2.le.3)then
        kkk=1
        kmin=2
        kmax=3
        nk3=0
      elseif(k1.ge.4)then
        kkk=2
        kmin=4
        kmax=9
        nk3=0
      else
        kkk=3
        kmin=4
        kmax=9
        nk3=1
      endif
      fev=1.d9
      a=1d-9
      ptyp(1)='    [g]    '
      ptyp(2)='   e^-!    '
      ptyp(3)='   e^+!    '
      ptyp(4)='  nucleons '
      ptyp(5)='  [p]^+/-! '
      ptyp(6)='  d[g]/dZ  '
      ptyp(7)='   K^+/-!  '
      ptyp(8)='   K?l!    '
      ptyp(9)='   K?s!    '
      ptyp(10)='  [m]^+/-! '
      ptyp(11)='  [m]^+! '
      ptyp(12)='  [m]^-! '

      nshower=1
      call NormalizeTables(1.d0/dble(max(1,nshower)))
      do jz=1,numiz
        iz=jz
        zi=obslvl(iz)
        write(c(jz),'(i6)')nint(zi)
      enddo

      write(ifho,'(a)')'orientation landscape'
      write(ifho,'(a)')'newpage'

c Generation Profile **********************************************
      write(ifho,'(a,a)')'!************************'
     &           ,' Generation profile  ************************'
      write(ifho,'(a)')'zone 1 1 1'
      do k=1,2
        if(cntgen(0,k).gt.0d0)then
          genmn=0d0
          genmx=0d0
          do i=1,ngenmx
            cntgen(i,k)=cntgen(i,k)/ cntgen(0,k)
            if(cntgen(i,k).le.0d0.and.genmx.le.0d0)genmn=dble(i)
            if(cntgen(i,k).gt.0d0)genmx=dble(i)
          enddo

          if(k.eq.1)write(ifho,'(a)')  'openhisto name GenerationEM'
          if(k.eq.2)write(ifho,'(a)')  'openhisto name GenerationHad'
          write(ifho,'(a)')  'htyp his'
          write(ifho,'(a)')  'xmod lin ymod lin'
          write(ifho,*)      'xrange ',genmn,genmx
          write(ifho,'(a)')  'yrange auto auto '
          write(ifho,'(a,1p,e10.3,a)')
     &         'txt "title Prim. Engy (eV) =',eprima*fev,'"'
          write(ifho,'(a)')  'txt  "xaxis number of Had. generations"'
          write(ifho,'(a)')  'txt  "yaxis Probability"'
          write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
c column names
          write(ifho,'(a)')'! Gener nbr  !   Probability '
c end column names
          write(ifho,*)'array 2'
          do ix=1,ngenmx
            write(ifho,'(i3,1p,e13.5)')ix,cntgen(ix,k)
          end do
          write(ifho,'(a)')       '  endarray'
          write(ifho,'(a)')       'closehisto'
          if(k.eq.1)write(ifho,'(a,a)')  ' plot 0-'
          if(k.eq.2)write(ifho,'(a,a)')  ' plot 0'
        endif
      enddo

      if(numiz.gt.1)then
      write(ifho,'(a)')'!---------------yieldz------------'
      write(ifho,'(a)')  'zone 2 2 1 openhisto name xyi'
      write(ifho,'(a)')  'htyp lin'
      write(ifho,'(a)')  'xmod lin ymod lin'
      write(ifho,*)  'xrange ',zamin,zamax
      write(ifho,'(a)')  'yrange auto auto '
      write(ifho,'(a,1p,e10.3,a)')  'text 0.1 0.9 "Prim. Engy (eV) ='
     &                           ,eprima*fev,'"'
      write(ifho,'(a)')  '- txt  "xaxis height H (m)"'
      if(kkk.eq.1)then
        write(ifho,'(a)')  '+ txt  "yaxis number of charged"'
      else
        write(ifho,'(a)')  '+ txt  "yaxis number of charged hadrons"'
      endif
      do ip=k1,k2
      write(ifho,'(a)')  '+ txt  "yaxis number of '//ptyp(ip)//'"'
      enddo
      if(kkk.eq.3)then
        write(ifho,'(a)')  '+ txt  "yaxis number of e^+/-!"'
      endif
      write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
      write(ifho,*)  'array ',-1-(k2-k1+1)-1-nk3
      do iz=1,numiz
        zi=obslvl(iz)
        stot(iz)=0.d0
        do ip=kmin,kmax
          if(ip.ne.6.and.ip.ne.1)stot(iz)=stot(iz)+yieldz(ip,iz)
        enddo
        if(kkk.ne.3)then
          write(ifho,'(e13.5,90e11.3)')zi,stot(iz)
     & ,(yieldz(ip,iz),ip=k1,k2)
        else
          write(ifho,'(e13.5,90e11.3)')zi,stot(iz)
     & ,(yieldz(ip,iz),ip=k1,k2),yieldz(2,iz)+yieldz(3,iz)
        endif
      end do
      write(ifho,'(a)')       '  endarray'
      write(ifho,'(a)')       'closehisto'
      do ip=1,k2-k1+2+nk3
       if(ip.le.9)write(ifho,'(a,i1,$)')  ' plot xyi+',ip
       if(ip.gt.9)write(ifho,'(a,i2,$)')  ' plot xyi+',ip
      enddo
      write(ifho,'(a)')' '
      endif

        k2m=min(k2,10)
      write(ifho,'(a)')'!---------------yiex------------'
      write(ifho,'(a)')  'zone 2 2 1 openhisto name xti'
      write(ifho,'(a)')  'htyp lin'
      write(ifho,'(a)')  'xmod log ymod lin'
      write(ifho,*)   'xrange ',tamin,tamax
      write(ifho,'(a)')  'yrange auto auto '
      write(ifho,'(a,1p,e10.3,a)')  'text 0.1 0.9 "Prim. Engy (eV) ='
     &                           ,eprima*fev,'"'
      write(ifho,'(a)')  '- txt  "xaxis time t (ns)"'
      if(kkk.eq.1)then
        write(ifho,'(a)')  '+ txt  "yaxis number of e+/e-"'
      else
        write(ifho,'(a)')  '+ txt  "yaxis number of hadrons"'
      endif
      do iz=1,numiz
      do ip=k1,k2m
      write(ifho,'(a)')  '+ txt  "yaxis number of '//ptyp(ip)//
     &' at '//c(iz)//'m"'
      enddo
      write(ifho,'(a)')  '+ txt  "yaxis charge"'
      if(kkk.eq.3)then
        write(ifho,'(a)')  '+ txt  "yaxis number of e^+/-!"'
      endif
      enddo
      write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
      write(ifho,*)  'array ',-1-numiz*((k2m-k1+1)+2+nk3)
      do it=1,numie
        ti=tamin*ctime**dble(it-1)
        do iz=1,numiz
        stot(iz)=0.d0
        do ip=kmin,kmax
          if(ip.ne.6)stot(iz)=stot(iz)+yiex(ip,iz,it)
        enddo
        enddo
        if(kkk.ne.3)then
          write(ifho,'(e13.5,900e16.8)')ti,(stot(iz)
     & ,(yiex(ip,iz,it),ip=k1,k2m),yiex(11,iz,it),iz=1,numiz)
        else
          write(ifho,'(e13.5,900e16.8)')ti,(stot(iz)
     & ,(yiex(ip,iz,it),ip=k1,k2m),yiex(11,iz,it)
     & ,yiex(2,iz,it)+yiex(3,iz,it),iz=1,numiz)
        endif
      end do
      write(ifho,'(a)')       '  endarray'
      write(ifho,'(a)')       'closehisto'
      do ip=1,numiz*(k2m-k1+3+nk3)
       if(ip.le.9)write(ifho,'(a,i1,$)')  ' plot xti+',ip
       if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i2,$)')  ' plot xti+',ip
       if(ip.gt.99)write(ifho,'(a,i3,$)')  ' plot xti+',ip
      enddo
      write(ifho,'(a)')' '

      anorm=dble(max(1,nshower))
      write(ifho,'(a)')'!--------------yieldr----------------------'
      write(ifho,'(a)')  'zone 3 6 1 openhisto name xyr'
      write(ifho,'(a)')  'htyp pnt'
      write(ifho,'(a)')  'xmod log ymod log'
      write(ifho,'(a,2e11.3)')  'xrange ',ramin,ramax
      write(ifho,'(a)')  'yrange auto auto '
      write(ifho,'(a)')'- txt "xaxis R (m) "'
      do ip=k1,k2
      do jj=1,numiz        !radial density normalized
      write(ifho,'(a)')
     *'++ txt "yaxis '//ptyp(ip)//'?norm! (h='//c(jj)//')"'
      enddo
      enddo
      do ip=k1,k2       !radial density
      do jj=1,numiz
      write(ifho,'(a)')'++ txt "yaxis '//ptyp(ip)//' (h='//c(jj)//')"'
      enddo
      enddo
      write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
      write(ifho,*)'array ',-1-4*numiz*(k2-k1+1)
      do ir=1,numir
        rr= ramin*(ramax/ramin)**((dble(ir)-0.5d0)/dble(numir))
        ra= (ramin*(ramax/ramin)**((dble(ir)-1.d0)/dble(numir)))**2.d0
        rb= (ramin*(ramax/ramin)**((dble(ir)-0.d0)/dble(numir)))**2.d0
        d=(rb-ra)*3.14159
        write(ifho,'(900e11.3)')rr
     & ,((abs(yieldr(ip,jz,ir))/d/max(a,yieldz(ip,jz))
     &   ,sqrt(max(0.d0,yieldr2(ip,jz,ir)-yieldr(ip,jz,ir)**2)
     &   /max(anorm-1.d0,1.d0))/d/max(a,yieldz(ip,jz))
     & ,jz=1,numiz),ip=k1,k2),((abs(yieldr(ip,jz,ir))/d
     &   ,sqrt(max(0.d0,yieldr2(ip,jz,ir)-yieldr(ip,jz,ir)**2)
     &   /max(anorm-1.d0,1.d0))/d,jz=1,numiz),ip=k1,k2)
      end do
      write(ifho,'(a)')       '  endarray'
      write(ifho,'(a)')       'closehisto'
      do ip=1,2*numiz*(k2-k1+1)
       if(ip.le.9)write(ifho,'(a,i1,$)')  ' plot xyr+',ip
       if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i2,$)')  ' plot xyr+',ip
       if(ip.gt.99)write(ifho,'(a,i3,$)')  ' plot xyr+',ip
      enddo
      write(ifho,'(a)')' '





      write(ifho,'(a)')'!--------------spec----------------------'

        do jz=1,numiz
          iz=jz
          zi=obslvl(iz)
          do jr=1,maxjr
            ir=irfirst+(jr-1)*modr
            ri= ramin*(ramax/ramin)**(dble(ir)/dble(numir))
            tx(jr,jz)(7:7)=','
            write(tx(jr,jz)(1:6),'(i6)')nint(zi)
            write(tx(jr,jz)(8:17),'(f9.3)')ri
            ie=iefirst+(jr-1)*moden
            do ip=k1,k2
              if(ip.le.3)then
                k=1
              else
                k=2
              endif
              ei=eamin(k)*(eamax(k)/eamin(k))**(dble(ie)/dble(numie))
              ty(jr,jz,ip)(7:7)=','
              write(ty(jr,jz,ip)(1:6),'(i6)')nint(zi)
              write(ty(jr,jz,ip)(8:17),'(f9.3)')ei
            enddo
          enddo
        enddo


      if(kkk.ne.3)then
        kk1=kkk
        kk2=kkk
        kmn(kkk)=k1
        kmx(kkk)=k2
      else
        kk1=1
        kk2=2
        kmn(1)=k1
        kmx(1)=3
        kmn(2)=4
        kmx(2)=k2
      endif


      do k=kk1,kk2              !loop particle type (EM, hadron or both)

        if(eamin(k).lt.eamax(k))then
        dle=log10(eamax(k)/eamin(k))/numie

        write(ifho,'(a,i1)')  'zone 2 3 1 openhisto name xsp',k
        write(ifho,'(a)')  'htyp lin'
        write(ifho,'(a)')  'xmod log ymod log'
        write(ifho,'(a,2e11.3)')  'xrange ',eamin(k),eamax(k)
        write(ifho,'(a)')  'yrange auto auto '
        write(ifho,'(a)')'- txt "xaxis energy (GeV)" '
        do 10 ip=kmn(k),kmx(k)
        do 10 j=1,numiz
 10       write(ifho,'(a)')
     &   '+ txt "yaxis EdN/dE '//ptyp(ip)(1:11)//' (h='//c(j)//')"'
        write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
        write(ifho,*)'array ',-1-numiz*(kmx(k)-kmn(k)+1)
        do nem=1,numiz
          emoy(nem)=0.d0
          sum(nem)=0.d0
        enddo
        do i=1,numie+1
          ee= eamin(k)*(eamax(k)/eamin(k))**(dble(i-1)/dble(numie))
          d=dle
          if(i.eq.1)d=0.5d0*dle       !first bin is half size because of the cutoff
          do nem=1,numiz
          emoy(nem)=emoy(nem)+ee*(spec(0,2,i,nem)
     &             +spec(0,3,i,nem))/d !mean energy of e+/e-
          sum(nem)=sum(nem)+(spec(0,2,i,nem)+spec(0,3,i,nem))/d
          enddo
          write(ifho,'(300e11.3)')ee
     &         ,((spec(0,ip,i,j)/d,j=1,numiz),ip=kmn(k),kmx(k))
        end do
        write(ifho,'(a)')'  endarray'
        do nem=1,numiz
          if(sum(nem).gt.0.d0)then
            emoy(nem)=emoy(nem)/sum(nem)
          else
            emoy(nem)=0.d0
          endif
        enddo
        do j=1,numiz
      write(ifho,'(a,i1,a,1p,e10.4,a)')'text 0.5 0.9 "Eel(',j-1,')= '
     &                                               ,emoy(j),' GeV"'
        enddo
        write(ifho,'(a)')'closehisto'
        do ip=1,numiz*(kmx(k)-kmn(k)+1)
          if(ip.le.9)write(ifho,'(a,i1,a,i1,$)')  ' plot xsp',k,'+',ip
          if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i1,a,i2,$)')
     & ' plot xsp',k,'+',ip
          if(ip.gt.99)write(ifho,'(a,i1,a,i3,$)')  ' plot xsp',k,'+',ip
        if(mod(ip,10).eq.0.or.ip.eq.3*(kmx(k)-kmn(k)+1))write(ifho,*)' '
        enddo


        if(musmm.ge.0)then
        write(ifho,'(a)')
        write(ifho,'(a)')  'resethisto'
        write(ifho,'(a,i1)')  'openhisto name xmo',k
        write(ifho,'(a)')  'htyp lin'
        write(ifho,'(a)')  'xmod log ymod log'
        write(ifho,'(a,2e11.3)')  'xrange ',eamin(k),eamax(k)
        write(ifho,'(a)')  'yrange auto auto '
        write(ifho,'(a)')'- txt "xaxis energy (GeV)" '
        do 1 ii=1,2
        do 1 ip=kmn(k),kmx(k)
        do 1 m=0,musmm
        do 1 j=1,numiz
  1     write(ifho,'(a,i2,1x,a)')
     &   '+ txt "yaxis M',m,ptyp(ip)(1:11)//' (h='//c(j)//')"'
        write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
        write(ifho,*)'array ',-1-2*numiz*(musmm+1)*(kmx(k)-kmn(k)+1)
        do i=1,numie+1
          ee= eamin(k)*(eamax(k)/eamin(k))**(dble(i-1)/dble(numie))
          d=dle
          if(i.eq.1)d=0.5d0*dle       !first bin is half size because of the cutoff
        write(ifho,'(800e11.3)')ee
     & ,((spec(0,ip,i,j)/d,j=1,numiz)
     &  ,((spec(m,ip,i,j)/max(a,spec(0,ip,i,j)),j=1,numiz),m=1,musmm)
     & ,ip=kmn(k),kmx(k))
        write(ifho,'(800e11.3)')
     &  ((1d0-max(a,spex(0,ip,i,j))/max(a,spec(0,ip,i,j)),j=1,numiz)
     &  ,((1d0-max(a,spex(m,ip,i,j))/max(a,spec(m,ip,i,j)),j=1,numiz),
     &        m=1,musmm)
     &      ,ip=kmn(k),kmx(k))
        end do
        write(ifho,'(a/a)')'  endarray','closehisto'
        do ip=1,(musmm+1)*numiz*(kmx(k)-kmn(k)+1)
          if(ip.le.9)write(ifho,'(a,i1,a,i1,$)')  ' plot xmo',k,'+',ip
          if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i1,a,i2,$)')
     &                                            ' plot xmo',k,'+',ip
          if(ip.gt.99.and.ip.le.999)write(ifho,'(a,i1,a,i3,$)')
     &                                            ' plot xmo',k,'+',ip
          if(ip.gt.999)write(ifho,'(a,i1,a,i4,$)') ' plot xmo',k,'+',ip
        if(mod(ip,10).eq.0.or.ip.eq.3*(kmx(k)-kmn(k)+1))write(ifho,*)' '
       if(mod(ip,10).eq.0.or.ip.eq.12*(kmx(k)-kmn(k)+1))write(ifho,*)' '
        enddo
        do ip=(musmm+1)*numiz*(kmx(k)-kmn(k)+1)+1,
     *        (musmm+1)*2*numiz*(kmx(k)-kmn(k)+1)
       if(ip.le.9)write(ifho,'(a,i1,a,i1,$)')' plot -ymod lin xmo',k,'+'
     *,ip
       if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i1,a,i2,$)')
     *                                   ' plot -ymod lin xmo',k,'+',ip
       if(ip.gt.99.and.ip.le.999)write(ifho,'(a,i1,a,i3,$)')
     *                                   ' plot -ymod lin xmo',k,'+',ip
       if(ip.gt.999)write(ifho,'(a,i1,a,i4,$)')' plot -ymod lin xmo',
     *  k,'+',ip
       if(mod(ip,10).eq.0.or.ip.eq.24*(kmx(k)-kmn(k)+1))write(ifho,*)' '
        enddo

        endif
        write(ifho,'(a)')
        write(ifho,'(a)')  'resethisto'


        write(ifho,'(a)')'!--------------specr----------------------'
        write(ifho,'(a,i1)')  'zone 3 6 1 openhisto name xsr',k
        write(ifho,'(a)')  'htyp lin'
        write(ifho,'(a)')  'xmod log ymod log'
        write(ifho,'(a,2e11.3)')  'xrange ',eamin(k),eamax(k)
        write(ifho,'(a)')  'yrange auto auto '
        write(ifho,'(a)')'- txt "xaxis energy (GeV)"'
        do 2 ip=kmn(k),kmx(k)
        do 2 jz=1,numiz
        do 2 jr=1,maxjr
 2        write(ifho,'(a)')
     &    '+ txt "yaxis '//ptyp(ip)//' ('//tx(jr,jz)//')"'
        write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
        write(ifho,*)'array ',-1-maxjr*numiz*(kmx(k)-kmn(k)+1)
        do ie=1,numie+1
          ee= eamin(k)*(eamax(k)/eamin(k))**(dble(ie-1)/dble(numie))
          d=dle
          if(ie.eq.1)d=0.5d0*dle !first bin is half size because of the cutoff
          write(ifho,'(900e11.3)')ee
     &   ,(((specr(ip,ie,jz,jr)/d,jr=1,maxjr)
     &   ,jz=1,numiz),ip=kmn(k),kmx(k))
        end do
        write(ifho,'(a/a)')'  endarray','closehisto'
        do ip=1,maxjr*numiz*(kmx(k)-kmn(k)+1)
          if(ip.le.9)write(ifho,'(a,i1,a,i1,$)')  ' plot xsr',k,'+',ip
          if(ip.gt.9.and.ip.le.99)write(ifho,'(a,i1,a,i2,$)')
     &                                       ' plot xsr',k,'+',ip
         if(ip.gt.99)write(ifho,'(a,i1,a,i3,$)')  ' plot xsr',k,'+',ip
        if(mod(ip,10).eq.0.or.ip.eq.9*(kmx(k)-kmn(k)+1))write(ifho,*)' '
        enddo
        write(ifho,'(a)')
        write(ifho,'(a)')  'resethisto'


      endif

      enddo   !---> end loop k (EM and hadron)



      do iyx=1,maxiex

        write(nyx,'(a,i1)')'xx',iyx
        write(ifho,'(3a)')'!-------------- ',nyx,' --------------------'
        write(ifho,'(a)')  'zone 3 6 1 openhisto name '//nyx
        write(ifho,'(a)')  'htyp lin'
        write(ifho,'(a)')  'xmod lin ymod log'
        write(ifho,'(a,2e11.3)')  'xrange ',xamin(iyx),xamax(iyx)
        write(ifho,'(a)')  'yrange auto auto '
        if(iyx.eq.1)then
          write(ifho,'(a)')'- txt "xaxis cos([f])"'
        elseif(iyx.eq.2)then
          write(ifho,'(a)')'- txt "xaxis sin^2!([q])"'
        elseif(iyx.eq.3)then
          write(ifho,'(a)')'- txt "xaxis t (ns)"'
        elseif(iyx.eq.4)then
          write(ifho,'(a)')'- txt "xaxis [D][f]"'
        elseif(iyx.eq.5)then
          write(ifho,'(a)')'- txt "xaxis log?10!(sin^2!([q])"'
        endif
        do ip=k1,k2
          do jz=1,numiz
            do jr=1,maxjr
              if(iyx.eq.5)then
                write(ifho,'(a)')
     &      '+ txt "yaxis '//ptyp(ip)//' ('//ty(jr,jz,ip)//')"'
              else
                write(ifho,'(a)')
     &      '+ txt "yaxis '//ptyp(ip)//' ('//tx(jr,jz)//')"'
              endif
            enddo
          enddo
        enddo
        write(ifho,'(a,d22.14)')  'histoweight ',dble(nshower)
        write(ifho,*)'array ',-1-maxjr*numiz*(k2-k1+1)
        d=(xamax(iyx)-xamin(iyx))/numix(iyx)
        do ix=1,numix(iyx)
          xx= xamin(iyx)+(xamax(iyx)-xamin(iyx))
     &                  *((dble(ix)-0.5d0)/dble(numix(iyx)))
          if(iyx.eq.5)then !log10(sin2(theta)) is normalized
            write(ifho,'(900e11.3)')xx
     & ,(((yieldx(iyx,ip,jz,je,ix)/d/max(1d-9,speca(ip,jz,je))
     & ,je=1,maxjr),jz=1,numiz),ip=k1,k2)
          else
            write(ifho,'(900e11.3)')xx
     & ,(((yieldx(iyx,ip,jz,jr,ix)/d,jr=1,maxjr),jz=1,numiz),ip=k1,k2)
          endif
        end do
        write(ifho,'(a/a)')'  endarray','closehisto'
        do ip=1,maxjr*numiz*(k2-k1+1)
          if(ip.le.9)write(ifho,'(a,a3,a1,i1,$)')  ' plot ',nyx,'+',ip
          if(ip.gt.9.and.ip.le.99)write(ifho,'(a,a3,a1,i2,$)')
     &  ' plot ',nyx,'+',ip
          if(ip.gt.99)write(ifho,'(a,a3,a1,i3,$)')  ' plot ',nyx,'+',ip
          if(mod(ip,10).eq.0.or.ip.eq.maxjr*numiz*(k2-k1+1))
     &    write(ifho,*)' '
        enddo

      enddo


      end

c----------------------------------------------------------------------
      subroutine NormalizeTables(value)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer maxo,maxmm,musmm,iimom,i1mom,i2mom,i3mom
      parameter(maxo=2,maxmm=maxo)
      common/cxmom1/musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx,2)
      common/countgen/cntgen
      double precision zamin,zamax,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,numie
      common/cxtime/ tamin,tamax,ctime
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden
      parameter(maxiz=10,maxjz=maxiz,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz,maxie)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)

      do in=1,maxin
        do iz=1,numiz
          yieldz(in,iz)=yieldz(in,iz)*value
          do it=1,numie
            yiex(in,iz,it)=yiex(in,iz,it)*value
          enddo
        enddo
        do jz=1,numiz
          do ir=1,numir
            yieldr(in,jz,ir)=yieldr1(in,jz,ir)*value
            yieldr2(in,jz,ir)=yieldr1(in,jz,ir)**2*value
          enddo
          do ie=1,numie+1
            do mm=0,maxmm
              spec(mm,in,ie,jz)=spec(mm,in,ie,jz)*value
              spex(mm,in,ie,jz)=spex(mm,in,ie,jz)*value
            enddo
            do jr=1,maxjr
              specr(in,ie,jz,jr)=specr(in,ie,jz,jr)*value
            enddo
          enddo
          do i=1,maxiex
            do ix=1,numix(i)
              do jr=1,maxjr
                 yieldx(i,in,jz,jr,ix)=yieldx(i,in,jz,jr,ix)*value
              enddo
            enddo
          enddo
          do jr=1,maxjr
            speca(in,jz,jr)=speca(in,jz,jr)*value
          enddo
        enddo
      enddo
      end

c                  end analysis part

c------------------------------------------------------------------------------
      integer function idtrafo(code1,code2,idi)
c------------------------------------------------------------------------------
c.....tranforms id of code1 (=idi) into id of code2 (=idtrafocx)
c.....supported codes:
c.....'nxs' = epos
c.....'pdg' = PDG 1996
c.....'qgs' = QGSJet
c.....'ghe' = Gheisha
c.....'sib' = Sibyll
c.....'cor' = Corsika (GEANT)
c.....'flk' = Fluka

C --- ighenex(I)=EPOS CODE CORRESPONDING TO GHEISHA CODE I ---

      common /ighnx/ ighenex(35)
      data ighenex/
     $               10,   11,   -12,    12,   -14,    14,   120,   110,
     $             -120,  130,    20,   -20,  -130,  1120, -1120,  1220,
     $            -1220, 2130, -2130,  1130,  1230,  2230, -1130, -1230,
     $            -2230, 1330,  2330, -1330, -2330,    17,    18,    19,
     $            3331, -3331,  30/

C --- DATA STMTS. FOR GEANT/GHEISHA PARTICLE CODE CONVERSIONS ---
C --- KIPART(I)=GHEISHA CODE CORRESPONDING TO GEANT   CODE I ---
C --- IKPART(I)=GEANT   CODE CORRESPONDING TO GHEISHA CODE I ---
      DIMENSION        KIPART(48)!,IKPART(35)
      DATA KIPART/
     $               1,   3,   4,   2,   5,   6,   8,   7,
     $               9,  12,  10,  13,  16,  14,  15,  11,
     $              35,  18,  20,  21,  22,  26,  27,  33,
     $              17,  19,  23,  24,  25,  28,  29,  34,
     $              35,  35,  35,  35,  35,  35,  35,  35,
     $              35,  35,  35,  35,  30,  31,  32,  35/

c      DATA IKPART/
c     $               1,   4,   2,   3,   5,   6,   8,   7,
c     $               9,  11,  16,  10,  12,  14,  15,  13,
c     $              25,  18,  26,  19,  20,  21,  27,  28,
c     $              29,  22,  23,  30,  31,  45,  46,  47,
c     $              24,  32,  48/
      INTEGER          ICFTABL(200),IFCTABL(-6:100)
C  ICTABL CONVERTS CORSIKA PARTICLES INTO FLUKA PARTICLES
C  FIRST TABLE ONLY IF CHARMED PARTICLES CAN BE TREATED
      DATA ICFTABL/
     *   7,   4,   3,   0,  10,  11,  23,  13,  14,  12,  ! 10
     *  15,  16,   8,   1,   2,  19,   0,  17,  21,  22,  ! 20
     *  20,  34,  36,  38,   9,  18,  31,  32,  33,  34,  ! 30
     *  37,  39,  24,  25, 6*0,
     *  0,    0,   0,   0,  -3,  -4,  -6,  -5,   0,   0,  ! 50
     *  10*0,
     *   0,   0,   0,   0,   0,   5,   6,  27,  28,   0,  ! 70
     *  10*0,
     *  10*0,
     *  10*0,                                             !100
     *  10*0,
     *   0,   0,   0,   0,   0,  47,  45,  46,  48,  49,  !120
     *  50,   0,   0,   0,   0,   0,   0,   0,   0,   0,  !130
     *  41,  42,  43,  44,   0,   0,  51,  52,  53,   0,  !140
     *   0,   0,  54,  55,  56,   0,   0,   0,  57,  58,  !150
     *  59,   0,   0,   0,  60,  61,  62,   0,   0,   0,  !160
     *  40*0/
C  IFCTABL CONVERTS FLUKA PARTICLES INTO CORSIKA PARTICLES
      DATA IFCTABL/
     *                402, 302, 301, 201,   0,   0,   0,
     *  14,  15,   3,   2,  66,  67,   1,  13,  25,   5,
     *   6,  10,   8,   9,  11,  12,  18,  26,  16,  21,
     *  19,  20,   7,  33,  34,   0,  68,  69,   0,   0,
     *  27,  28,  29,  22,  30,  23,  31,  24,  32,   0,
     * 131, 132, 133, 134, 117, 118, 116, 119, 120, 121,
     * 137, 138, 139, 143, 144, 145, 149, 150, 151, 155,
     * 156, 157,   0,   0,   36*0/
c-------------------------------------------------------------------------------

      character*3 code1,code2
      parameter (ncode=5,nidt=365)
      integer idt(ncode,nidt)

c            nxs|pdg|qgs|cor|sib
      data ((idt(i,j),i=1,ncode),j= 1,68)/
     *          1,2,99,99,99             !u quark
     *     ,    2,1,99,99,99             !d
     *     ,    3,3,99,99,99             !s
     *     ,    4,4,99,99,99             !c
     *     ,    5,5,99,99,99             !b
     *     ,    6,6,99,99,99             !t
     *     ,   10,22,99,1,1              !gamma
     *     ,   9 ,21,99,99,99            !gluon
     *     ,   12,11,11,3,3              !e-
     *     ,  -12,-11,-11,2,2            !e+
     *     ,   11,12,99,66,15            !nu_e-
     *     ,  -11,-12,99,67,16           !nu_e+
     *     ,   14,13,99,6,5              !mu-
     *     ,  -14,-13,99,5,4             !mu+
     *     ,   13,14,99,68,17            !nu_mu-
     *     ,  -13,-14,99,69,18           !nu_mu+
     *     ,   16,15,99,132,19           !tau-
     *     ,  -16,-15,99,131,-19         !tau+
     *     ,   15,16,99,133,20           !nu_tau-
     *     ,  -15,-16,99,134,-20         !nu_tau+
     *     ,  110,111,0,7,6              !pi0
     *     ,  120,211,1,8,7              !pi+
     *     , -120,-211,-1,9,8            !pi-
     *     ,  220,221,10,17,23           !eta
     *     ,  130,321,4,11,9             !k+
     *     , -130,-321,-4,12,10          !k-
     *     ,  230,311,5,33,21            !k0
     *     , -230,-311,-5,34,22          !k0b
     *     ,   20,310,5,16,12            !kshort
     *     ,  -20,130,-5,10,11           !klong
     *     ,  330,331,99,99,24           !etaprime
     *     ,  111,113,19,51,27           !rho0
     *     ,  121,213,99,52,25           !rho+
     *     , -121,-213,99,53,26          !rho-
     *     ,  221,223,99,50,32           !omega
     *     ,  131,323,99,63,28           !k*+
     *     , -131,-323,99,64,29          !k*-
     *     ,  231,313,99,62,30           !k*0
     *     , -231,-313,99,65,31          !k*0b
     *     ,  331,333,99,99,33           !phi
     *     , -140,421,8,116,99           !D0(1.864)
     *     ,  140,-421,8,119,99          !D0b(1.864)
     *     , -240,411,7,117,99           !D(1.869)+
     *     ,  240,-411,7,118,99          !Db(1.869)-
     *     , 1120,2212,2,14,13           !proton
     *     , 1220,2112,3,13,14           !neutron
     *     , 2130,3122,6,18,39           !lambda
     *     , 1130,3222,99,19,34          !sigma+
     *     , 1230,3212,99,20,35          !sigma0
     *     , 2230,3112,99,21,36          !sigma-
     *     , 1330,3322,99,22,37          !xi0
     *     , 2330,3312,99,23,38          !xi-
     *     , 1111,2224,99,54,40          !delta++
     *     , 1121,2214,99,55,41          !delta+
     *     , 1221,2114,99,56,42          !delta0
     *     , 2221,1114,99,57,43          !delta-
     *     , 1131,3224,99,99,44          !sigma*+
     *     , 1231,3214,99,99,45          !sigma*0
     *     , 2231,3114,99,99,46          !sigma*-
     *     , 1331, 3324,99,99,47         !xi*0
     *     , 2331, 3314,99,99,48         !xi*-
     *     , 3331, 3334,99,24,49         !omega-
     *     , 2140, 4122,9,137,99         !LambdaC(2.285)+
     *     ,17,1000010020,99,201,1002            !  Deuteron
     *     ,18,1000010030,99,301,1003            !  Triton
     *     ,19,1000020040,99,402,1004            !  Alpha
     *     ,0,0,99,0,0                  !  Air
     *     ,99,99,99,99,99 /             !  unknown
      data ((idt(i,j),i=1,ncode),j= 69,91)/
     $      -340,431,99,120,99           !  Ds+
     $     ,340,-431,99,121,99           !  Ds-
     $     ,-241,413,99,124,99           !  D*+
     $     ,241,-413,99,125,99           !  D*-
     $     ,-141,423,99,123,99           !  D*0
     $     ,141,-423,99,126,99           !  D*0b
     $     ,-341,433,99,127,99           !  Ds*+
     $     ,341,-433,99,128,99           !  Ds*-
     $     ,-440,441,99,122,99           !  etac
     $     ,440,-441,99,122,99           !  etacb
     $     ,-441,443,99,130,99           !  J/psi
     $     ,441,-443,99,130,99           !  J/psib
     $     ,2240,4112,99,142,99          !  sigmac0
     $     ,1240,4212,99,141,99          !  sigmac+
     $     ,1140,4222,99,140,99          !  sigmac++
     $     ,2241,4114,99,163,99          !  sigma*c0
     $     ,1241,4214,99,162,99          !  sigma*c+
     $     ,1141,4224,99,161,99          !  sigma*c++
     $     ,3240,4132,99,139,99          !  Xic0
     $     ,2340,4312,99,144,99          !  Xi'c0
     $     ,3140,4232,99,138,99          !  Xic+
     $     ,1340,4322,99,143,99          !  Xi'c+
     $     ,3340,4332,99,145,99 /        !  omegac0
      data ((idt(i,j),i=1,ncode),j= 92,nidt)/
     $       1112,32224,99,99,99         !  Delta(1600)++
     $     , 1112, 2222,99,99,99         !  Delta(1620)++
     $     , 1113,12224,99,99,99         !  Delta(1700)++
     $     , 1114,12222,99,99,99         !  Delta(1900)++
     $     , 1114, 2226,99,99,99         !  Delta(1905)++
     $     , 1114,22222,99,99,99         !  Delta(1910)++
     $     , 1114,22224,99,99,99         !  Delta(1920)++
     $     , 1114,12226,99,99,99         !  Delta(1930)++
     $     , 1114, 2228,99,99,99         !  Delta(1950)++
     $     , 2222,31114,99,99,99         !  Delta(1600)-
     $     , 2222, 1112,99,99,99         !  Delta(1620)-
     $     , 2223,11114,99,99,99         !  Delta(1700)-
     $     , 2224,11112,99,99,99         !  Delta(1900)-
     $     , 2224, 1116,99,99,99         !  Delta(1905)-
     $     , 2224,21112,99,99,99         !  Delta(1910)-
     $     ,2224,21114,99,99,99          !  Delta(1920)-
     $     ,2224,11116,99,99,99          !  Delta(1930)-
     $     ,2224, 1118,99,99,99          !  Delta(1950)-
     $     ,1122,12212,99,99,99          !  N(1440)+
     $     ,1123, 2124,99,99,99          !  N(1520)+
     $     ,1123,22212,99,99,99          !  N(1535)+
     $     ,1124,32214,99,99,99          !  Delta(1600)+
     $     ,1124, 2122,99,99,99          !  Delta(1620)+
     $     ,1125,32212,99,99,99          !  N(1650)+
     $     ,1125, 2216,99,99,99          !  N(1675)+
     $     ,1125,12216,99,99,99          !  N(1680)+
     $     ,1126,12214,99,99,99          !  Delta(1700)+
     $     ,1127,22124,99,99,99          !  N(1700)+
     $     ,1127,42212,99,99,99          !  N(1710)+
     $     ,1127,32124,99,99,99          !  N(1720)+
     $     ,1128,12122,99,99,99          !  Delta(1900)+
     $     ,1128, 2126,99,99,99          !  Delta(1905)+
     $     ,1128,22122,99,99,99          !  Delta(1910)+
     $     ,1128,22214,99,99,99          !  Delta(1920)+
     $     ,1128,12126,99,99,99          !  Delta(1930)+
     $     ,1128, 2218,99,99,99          !  Delta(1950)+
     $     ,1222,12112,99,99,99          !  N(1440)0
     $     ,1223, 1214,99,99,99          !  N(1520)0
     $     ,1223,22112,99,99,99          !  N(1535)0
     $     ,1224,32114,99,99,99          !  Delta(1600)0
     $     ,1224, 1212,99,99,99          !  Delta(1620)0
     $     ,1225,32112,99,99,99          !  N(1650)0
     $     ,1225, 2116,99,99,99          !  N(1675)0
     $     ,1225,12116,99,99,99          !  N(1680)0
     $     ,1226,12114,99,99,99          !  Delta(1700)0
     $     ,1227,21214,99,99,99          !  N(1700)0
     $     ,1227,42112,99,99,99          !  N(1710)0
     $     ,1227,31214,99,99,99          !  N(1720)0
     $     ,1228,11212,99,99,99          !  Delta(1900)0
     $     ,1228, 1216,99,99,99          !  Delta(1905)0
     $     ,1228,21212,99,99,99          !  Delta(1910)0
     $     ,1228,22114,99,99,99          !  Delta(1920)0
     $     ,1228,11216,99,99,99          !  Delta(1930)0
     $     ,1228, 2118,99,99,99          !  Delta(1950)0
     $     ,1233,13122,99,99,99          !  Lambda(1405)0
     $     ,1234, 3124,99,99,99          !  Lambda(1520)0
     $     ,1235,23122,99,99,99          !  Lambda(1600)0
     $     ,1235,33122,99,99,99          !  Lambda(1670)0
     $     ,1235,13124,99,99,99          !  Lambda(1690)0
     $     ,1236,13212,99,99,99          !  Sigma(1660)0
     $     ,1236,13214,99,99,99          !  Sigma(1670)0
     $     ,1237,23212,99,99,99          !  Sigma(1750)0
     $     ,1237, 3216,99,99,99          !  Sigma(1775)0
     $     ,1238,43122,99,99,99          !  Lambda(1800)0
     $     ,1238,53122,99,99,99          !  Lambda(1810)0
     $     ,1238, 3126,99,99,99          !  Lambda(1820)0
     $     ,1238,13126,99,99,99          !  Lambda(1830)0
     $     ,1238,23124,99,99,99          !  Lambda(1890)0
     $     ,1239,13216,99,99,99          !  Sigma(1915)0
     $     ,1239,23214,99,99,99          !  Sigma(1940)0
     $     ,1132,13222,99,99,99          !  Sigma(1660)+
     $     ,1132,13224,99,99,99          !  Sigma(1670)+
     $     ,1133,23222,99,99,99          !  Sigma(1750)+
     $     ,1133,3226,99,99,99           !  Sigma(1775)+
     $     ,1134,13226,99,99,99          !  Sigma(1915)+
     $     ,1134,23224,99,99,99          !  Sigma(1940)+
     $     ,2232,13112,99,99,99          !  Sigma(1660)-
     $     ,2232,13114,99,99,99          !  Sigma(1670)-
     $     ,2233,23112,99,99,99          !  Sigma(1750)-
     $     ,2233,3116,99,99,99           !  Sigma(1775)-
     $     ,2234,13116,99,99,99          !  Sigma(1915)-
     $     ,2234,23114,99,99,99          !  Sigma(1940)-
     $     ,5,7,99,99,99                 !  quark b'
     $     ,6,8,99,99,99                 !  quark t'
     $     ,16,17,99,99,99               !  lepton tau'
     $     ,15,18,99,99,99               !  lepton nu' tau
     $     ,90,23,99,99,99               !  Z0
     $     ,80,24,99,99,99               !  W+
     $     ,81,25,99,99,99               !  h0
     $     ,85,32,99,99,99               !  Z'0
     $     ,86,33,99,99,99               !  Z''0
     $     ,87,34,99,99,99               !  W'+
     $     ,82,35,99,99,99               !  H0
     $     ,83,36,99,99,99               !  A0
     $     ,84,37,99,99,99               !  H+
     $     ,1200,2101,99,99,99           !  diquark ud_0
     $     ,2300,3101,99,99,99           !  diquark sd_0
     $     ,1300,3201,99,99,99           !  diquark su_0
     $     ,2400,4101,99,99,99           !  diquark cd_0
     $     ,1400,4201,99,99,99           !  diquark cu_0
     $     ,3400,4301,99,99,99           !  diquark cs_0
     $     ,2500,5101,99,99,99           !  diquark bd_0
     $     ,1500,5201,99,99,99           !  diquark bu_0
     $     ,3500,5301,99,99,99           !  diquark bs_0
     $     ,4500,5401,99,99,99           !  diquark bc_0
     $     ,2200,1103,99,99,99           !  diquark dd_1
     $     ,1200,2103,99,99,99           !  diquark ud_1
     $     ,1100,2203,99,99,99           !  diquark uu_1
     $     ,2300,3103,99,99,99           !  diquark sd_1
     $     ,1300,3203,99,99,99           !  diquark su_1
     $     ,3300,3303,99,99,99           !  diquark ss_1
     $     ,2400,4103,99,99,99           !  diquark cd_1
     $     ,1400,4203,99,99,99           !  diquark cu_1
     $     ,3400,4303,99,99,99           !  diquark cs_1
     $     ,4400,4403,99,99,99           !  diquark cc_1
     $     ,2500,5103,99,99,99           !  diquark bd_1
     $     ,1500,5203,99,99,99           !  diquark bu_1
     $     ,3500,5303,99,99,99           !  diquark bs_1
     $     ,4500,5403,99,99,99           !  diquark bc_1
     $     ,5500,5503,99,99,99           !  diquark bb_1
     $     ,800000088,88,99,99,99        !  string junction  (pythia)
     $     ,800099999,99999,99,99,99     !  string (dmpjet)
     $     ,800000090,90,99,99,99        !  string (phojet)
     $     ,800000091,91,99,99,99        !  parton system in cluster fragmentation  (pythia)
     $     ,800000092,92,99,99,99        !  parton system in string fragmentation  (pythia)
     $     ,800000093,93,99,99,99        !  parton system in independent system  (pythia)
     $     ,800000094,94,99,99,99        !  CMshower (pythia)
     $     ,250,511,99,99,99             !  B0
     $     ,150,521,99,99,99             !  B+
     $     ,350,531,99,99,99             !  B0s+
     $     ,450,541,99,99,99             !  Bc+
     $     ,251,513,99,99,99             !  B*0
     $     ,151,523,99,99,99             !  B*+
     $     ,351,533,99,99,99             !  B*0s+
     $     ,451,543,99,99,99             !  B*c+
     $     ,550,551,99,99,99             !  etab
     $     ,551,553,99,99,99             !  Upsilon
     $     ,2341,4314,99,99,99           !  Xi*c0(2645)
     $     ,1341,4324,99,99,99           !  Xi*c+(2645)
     $     ,3341,4334,99,99,99           !  omega*c0
     $     ,2440,4412,99,99,99           !  dcc
     $     ,2441,4414,99,99,99           !  dcc*
     $     ,1440,4422,99,99,99           !  ucc
     $     ,1441,4424,99,99,99           !  ucc*
     $     ,3440,4432,99,99,99           !  scc
     $     ,3441,4434,99,99,99           !  scc*
     $     ,4441,4444,99,99,99           !  ccc*
     $     ,2250,5112,99,99,99           !  sigmab-
     $     ,2150,5122,99,99,99           !  lambdab0
     $     ,3250,5132,99,99,99           !  sdb
     $     ,4250,5142,99,99,99           !  cdb
     $     ,1250,5212,99,99,99           !  sigmab0
     $     ,1150,5222,99,99,99           !  sigmab+
     $     ,3150,5232,99,99,99           !  sub
     $     ,4150,5242,99,99,99           !  cub
     $     ,2350,5312,99,99,99           !  dsb
     $     ,1350,5322,99,99,99           !  usb
     $     ,3350,5332,99,99,99           !  ssb
     $     ,4350,5342,99,99,99           !  csb
     $     ,2450,5412,99,99,99           !  dcb
     $     ,1450,5422,99,99,99           !  ucb
     $     ,3450,5432,99,99,99           !  scb
     $     ,4450,5442,99,99,99           !  ccb
     $     ,2550,5512,99,99,99           !  dbb
     $     ,1550,5522,99,99,99           !  ubb
     $     ,3550,5532,99,99,99           !  sbb
     $     ,3550,5542,99,99,99           !  scb
     $     ,2251,5114,99,99,99           !  sigma*b-
     $     ,1251,5214,99,99,99           !  sigma*b0
     $     ,1151,5224,99,99,99           !  sigma*b+
     $     ,2351,5314,99,99,99           !  dsb*
     $     ,1351,5324,99,99,99           !  usb*
     $     ,3351,5334,99,99,99           !  ssb*
     $     ,2451,5414,99,99,99           !  dcb*
     $     ,1451,5424,99,99,99           !  ucb*
     $     ,3451,5434,99,99,99           !  scb*
     $     ,4451,5444,99,99,99           !  ccb*
     $     ,2551,5514,99,99,99           !  dbb*
     $     ,1551,5524,99,99,99           !  ubb*
     $     ,3551,5534,99,99,99           !  sbb*
     $     ,4551,5544,99,99,99           !  cbb*
     $     ,5551,5554,99,99,99           !  bbb*
     $     ,123,10213,99,99,99           !  b_1+
     $     ,122,10211,99,99,99           !  a_0+
     $     ,-143,10423,99,99,99          !  D_1(2420)0
     $     ,143,-10423,99,99,99          !  D_1(2420)0b
     $     ,-142,10421,99,99,99          !  D*_0(2400)0
     $     ,142,-10421,99,99,99          !  D*_0(2400)0b
     $     ,-243,10413,99,99,99          !  D_1(2420)+
     $     ,243,-10413,99,99,99          !  D_1(2420)-
     $     ,-242,10411,99,99,99          !  D*_0(2400)+
     $     ,242,-10411,99,99,99          !  D*_0(2400)-
     $     ,-343,20433,99,99,99          !  D_s1(2460)+
     $     ,343,-20433,99,99,99          !  D_s1(2460)-
     $     ,-342,10431,99,99,99          !  D_s0(2317)+
     $     ,342,-10431,99,99,99          !  D_s0(2317)-
     $     ,113,10113,99,99,99           !  b_10
     $     ,112,10111,99,99,99           !  a_00
     $     ,-443,20443,99,99,99          !  Xi_c1(1P)
     $     ,443,-20443,99,99,99          !  Xi_c1(1P)b
     $     ,-442,10441,99,99,99          !  Xi_c0(1P)
     $     ,442,-10441,99,99,99          !  Xi_c0(1P)b
     $     ,-444,10443,99,99,99          !  h_c(1P)
     $     ,444,-10443,99,99,99          !  h_c(1P)b
     $     ,253,10513,99,99,99           !  B_1(L)0 (db_10)
     $     ,252,10511,99,99,99           !  B*_00  (db*_00)
     $     ,153,10523,99,99,99           !  B_1(L)+ (ub_10)
     $     ,152,10521,99,99,99           !  B*_0+ (ub*_00)
     $     ,353,10533,99,99,99           !  B_s1(L)0 (sb_10)
     $     ,352,10531,99,99,99           !  B*_s0 (sb*_00)
     $     ,453,10543,99,99,99           !  B_c1(L)+ (cb_10)
     $     ,452,10541,99,99,99           !  B*_c0+ (cb*_00)
     $     ,553,20553,99,99,99           !  Xi_b1(1P) (Upsilon')
     $     ,552,10551,99,99,99           !  Xi_b0(1P) (Upsilon'*)
     $     ,124,20213,99,99,99           !  a_1+
     $     ,125,215,99,99,99             !  a_2+
     $     ,126,10215,99,99,99           !  pi_2+(1670)
     $     ,127,217,99,99,99             !  rho_3+(1690)
     $     ,232,10313,99,99,99           !  K0_1(1270)
     $     ,233,20313,99,99,99           !  K0_1(1400)
     $     ,234,315,99,99,99             !  K*0_2(1430)
     $     ,235,10311,99,99,99           !  K0(1460)
     $     ,132,10323,99,99,99           !  K+_1(1270)
     $     ,133,20323,99,99,99           !  K+_1(1400)
     $     ,134,325,99,99,99             !  K*+_2(1430)
     $     ,135,10321,99,99,99           !  K+(1460)
     $     ,-144,20423,99,99,99          !  D_1(2430)0
     $     ,144,-20423,99,99,99          !  D_1(2430)0b
     $     ,-145,425,99,99,99            !  D*_2(2460)0
     $     ,145,-425,99,99,99            !  D*_2(2460)0b
     $     ,-244,20413,99,99,99          !  D*_1(2430)+
     $     ,244,-20413,99,99,99          !  D*_1(2430)-
     $     ,-245,415,99,99,99            !  D*_2(2460)+
     $     ,245,-415,99,99,99            !  D*_2(2460)-
     $     ,-344,10433,99,99,99          !  D_s1(2536)+
     $     ,344,-10433,99,99,99          !  D_s1(2536)-
     $     ,-345,435,99,99,99            !  D*_s2(2573)+
     $     ,345,-435,99,99,99            !  D*_s2(2573)-
     $     ,114,20113,99,99,99           !  a_10
     $     ,115,115,99,99,99             !  a_20
     $     ,116,10115,99,99,99           !  pi_20(1670)
     $     ,117,117,99,99,99             !  rho_30(1690)
     $     ,222,9010221,99,99,99         !  f_00(980)
     $     ,223,10223,99,99,99           !  h_10(1170)
     $     ,224,225,99,99,99             !  f_20(1270)
     $     ,225,20223,99,99,99           !  f_10(1285)
     $     ,226,9030221,99,99,99         !  f_00(1500)
     $     ,332,20333,99,99,99           !  f_10(1420)
     $     ,333,335,99,99,99             !  f'_20(1525)
     $     ,-445,445,99,99,99            !  Xi_c2(1P)
     $     ,445,-445,99,99,99            !  Xi_c2(1P)b
     $     ,254,20513,99,99,99           !  B_1(H)0 (db*_10)
     $     ,255,515,99,99,99             !  B*_20 (db*_20)
     $     ,154,20523,99,99,99           !  B_1(H)+ (ub*_10)
     $     ,155,525,99,99,99             !  B*_2+ (ub*_20)
     $     ,354,20533,99,99,99           !  B_s1(H) (sb*_10)
     $     ,355,535,99,99,99             !  B*_s2 (sb*_20)
     $     ,454,20543,99,99,99           !  B*_c1 (cb*_10)
     $     ,455,545,99,99,99             !  B*_c2 (cb*_20)
     $     ,554,10553,99,99,99           !  h_b(1P) (bb*_10)
     $     ,555,555,99,99,99             !  Xi_b2(1P) (bb*_20)
     $     ,11099,9900110,99,99,99       !  diff pi0 state
     $     ,12099,9900210,99,99,99       !  diff pi+ state
     $     ,13099,9900320,99,99,99       !  diff K+ state
     $     ,22099,9900220,99,99,99       !  diff omega state
     $     ,2099,9900310,99,99,99       !  diff K0 state
     $     ,-2099,9900130,99,99,99       !  diff pi+ state
     $     ,33099,9900330,99,99,99       !  diff phi state
     $     ,44099,9900440,99,99,99       !  diff J/psi state
     $     ,112099,9902210,99,99,99      !  diff proton state
     $     ,122099,9902110,99,99,99      !  diff neutron state
     $     ,213099,9903120,99,99,99      !  diff lambda state
     $     ,800000110,110,99,99,99       !  Reggeon
     $     ,800000990,990,99,99,99 /     !  Pomeron


c      print *,'idtrafo',' ',code1,' ',code2,idi

      nidtmx=68
      id1=idi
      if(code1.eq.'nxs')then
        i=1
      elseif(code1.eq.'pdg')then
        i=2
      elseif(code1.eq.'qgs')then
        i=3
        if(id1.eq.-10)id1=19
      elseif(code1.eq.'cor')then
        i=4
      elseif(code1.eq.'sib')then
        i=5
      elseif(code1.eq.'ghe')then
        id1=ighenex(id1)
        i=1
      elseif(code1.eq.'flk')then
        id1=IFCTABL(id1)          !convert to corsika code
        i=4
      else
        stop "unknown code in idtrafo"
      endif
      if(code2.eq.'nxs')then
        j=1
        ji=j
        if(i.eq.2.and.id1.gt.1000000000)then   !nucleus from PDG
          idtrafo=id1 
          return
        elseif(i.eq.4.and.id1.gt.402)then               !nucleus from Corsika
          idtrafo=1000000000+mod(id1,100)*10000+(id1/100)*10   
          return
        elseif(i.eq.5.and.id1.gt.1004)then               !nucleus from Sibyll
          id1=(id1-1000)
          idtrafo=1000000000+id1/2*10000+id1*10   
          return
        elseif(id1.eq.130.and.i.eq.2)then
          idtrafo=-20
          return
        endif
        if(i.eq.2) nidtmx=nidt
        if(i.eq.4) nidtmx=89
      elseif(code2.eq.'pdg')then
        j=2
        ji=j
        if(i.eq.1.and.id1.gt.1000000000)then !nucleus from NEXUS
          idtrafo=id1 
          return
        elseif(i.eq.4.and.id1.gt.402)then               !nucleus from Corsika
          idtrafo=1000000000+mod(id1,100)*10000+(id1/100)*10   
          return
        elseif(i.eq.5.and.id1.gt.1004)then               !nucleus from Sibyll
          id1=(id1-1000)
          idtrafo=1000000000+id1/2*10000+id1*10   
          return
        elseif(id1.eq.-20.and.i.eq.1)then
          idtrafo=130
          return
        endif
        if(i.eq.1) nidtmx=nidt
        if(i.eq.4) nidtmx=89
       elseif(code2.eq.'qgs')then
        j=3
        ji=j
      elseif(code2.eq.'cor')then
        j=4
        ji=j
      elseif(code2.eq.'sib')then
        j=5
        ji=j
      elseif(code2.eq.'ghe')then
        j=4
        ji=6
      elseif(code2.eq.'flk')then
        j=4
        ji=7
        if(i.le.2) nidtmx=89
       else
        stop "unknown code in idtrafo"
      endif
      if(i.eq.4)then !corsika  id always >0 so convert antiparticles
        iadtr=id1
        if(iadtr.eq.25)then
          id1=-13
        elseif(iadtr.eq.15)then
          id1=-14
        elseif(iadtr.ge.26.and.iadtr.le.32)then
          id1=-iadtr+8
        elseif(iadtr.ge.58.and.iadtr.le.61)then
          id1=-iadtr+4
        elseif(iadtr.ge.149.and.iadtr.le.157)then
          id1=-iadtr+12
        elseif(iadtr.ge.171.and.iadtr.le.173)then
          id1=-iadtr+10
        endif
      endif
      iad1=abs(id1)
      isi=sign(1,id1)

      if(i.ne.j)then
      do n=1,nidtmx
        if(iad1.eq.abs(idt(i,n)))then
          m=1
          if(n+m.lt.nidt)then
            do while(abs(idt(i,n+m)).eq.iad1)
              m=m+1
            enddo
          endif
          mm=0
          if(m.gt.1)then
            if(m.eq.2.and.idt(i,n)*idt(i,n+1).lt.0)then
              if(id1.eq.idt(i,n+1))mm=1
              isi=1
            else
              mm=m
            endif
          else      !m=0 only one line, take care of sign
            if(idt(i,n).lt.0)isi=-isi
          endif
          idtrafo=idt(j,n+mm)*isi
          if(abs(idtrafo).eq.99)stop"New particle not allowed "
          if(idtrafo.lt.0.and.j.eq.4)then           !corsika  id always >0
            iadtr=abs(idtrafo)
            if(iadtr.eq.13)then
              idtrafo=25
            elseif(iadtr.eq.14)then
              idtrafo=15
            elseif(iadtr.ge.18.and.iadtr.le.24)then
              idtrafo=iadtr+8
            elseif(iadtr.ge.54.and.iadtr.le.57)then
              idtrafo=iadtr+4
            elseif(iadtr.ge.137.and.iadtr.le.145)then
              idtrafo=iadtr+12
            elseif(iadtr.ge.161.and.iadtr.le.163)then
              idtrafo=iadtr+10
            else
              idtrafo=iadtr
            endif
          elseif(idtrafo.eq.19.and.j.eq.3)then
            idtrafo=-10
          endif
          if(j.ne.ji)goto 100
          return
        endif
      enddo
      else
        idtrafo=id1
        if(j.ne.ji)goto 100
        return
      endif

c      print *, 'idtrafo: ',code1,' -> ', code2,id1,' not found.   '
c      stop
      idtrafo=0
      return

 100  if(j.eq.4)then            !corsika
        if(idtrafo.eq.201)then
          idtrafo=45
        elseif(idtrafo.eq.301)then
          idtrafo=46
        elseif(idtrafo.eq.402)then
          idtrafo=47
        elseif(idtrafo.eq.302)then
          idtrafo=48
        endif
        if(idtrafo.ne.0)then      !air
          if(ji.eq.6)then
            idtrafo=kipart(idtrafo)
          elseif(ji.eq.7)then
            idtrafo=ICFTABL(idtrafo)
          endif
        endif
        return
      else
        stop"Should not happen in idtrafo !"
      endif

      end

c-----------------------------------------------------------------------
      subroutine idmass(idi,amass)
c     returns the mass of the particle with ident code id.
c     (deuteron, triton and alpha mass come from Gheisha ???)
c-----------------------------------------------------------------------
      double precision amass
      dimension ammes(15,0:7),ambar0(30),ambar1(30)
      dimension amlep(52)
      parameter ( nqlep=41,nmes=2)
c-c   data amlep/.3,.3,.5,1.6,4.9,30.,-1.,-1.,0.,0.,
      data amlep/.002,.005,.100,1.50,4.20,171.,-1.,-1.,0.,0.,0.   !charm :1.30 , bottom: 4.75-> 4.20 
     *     ,.511003e-3
     *     ,0.,.105661,0.,1.807,1.87656,2.8167,3.755,.49767,.49767,
     *     100.3,100.3,100.5,101.6,104.9,130.,2*-1.,100.,0.,
     *     100.,100.005,100.,100.1,100.,101.8,2*-1.,100.,100.,
     *     11*0./
c     0- meson mass table
      data (ammes(m,0),m=1,15)
     *   /        .1349766,.13957018,.547853       !pi0,pi+-,eta
     *           ,.493677,.497614,.95778           !K+-, K0,eta'
     *    ,1.86483,1.86960,1.96847,2.9803          !D0,D+-,Ds,etac
     *    ,5.279,5.279,5.370,6.276,9.859/     !B+-,B0,Bs,Bc,etab
c     1- meson mass table
      data (ammes(m,1),m=1,15)
     *   /        .77549,.77549,.78265             !rho0,rho+-,omega
     *           ,.889166,.89594,1.019455          !K*+-,K0*,phi
     *     ,2.00693,2.01022,2.1123,3.096916        !D0*,D*+-,D*s,j/psi
     *     ,5.3251,5.3251,5.4154,6.610,9.46030/    !B*+-,B0*,B*s,B*c,upsilon
c     2- meson mass table
      data (ammes(m,2),m=1,15)
     *   /         .980,.980,.980             !a_00,a_0+-,f_0(980)
     *     ,1.270,1.270,1.4264                !K_1+-(1270),K_10(1270),f_1(1420)
     *     ,2.320,2.320,2.318,3.415           !D*_00,D*_0+,D_s0+,Xi_c0(1P)
     *     ,5.698,5.698,5.820,-1.,9.859/      !B*_0+,B*_00,B*_s0,B*_c0+,Xi_b0(1P)
c     3- meson mass table
      data (ammes(m,3),m=1,15)
     *   /        1.2295,1.2295,1.170          !b_10,b_1+-,h_1(1170)
     *     ,1.403,1.403,1.525                !K_1+-(1400),K_10(1400),f'_2(1525)
     *     ,2.422,2.422,2.460,3.511          !D_1(2420)0,D_1(2420)+,D_s1(2460)+,Xi_c1(1P)
     *     ,5.721,5.721,5.829,-1.,9.893/     !B_1(L)+,B_1(L)0,B_s1(L)0,B_c1(L)+,Xi_b1(1P)
c     4- meson mass table
      data (ammes(m,4),m=1,15)
     *   /        1.2300,1.2300,1.2751         !a_10,a_1+-,f_2(1270)
     *     ,1.4256,1.4256,-1.               !K*_2+-(1430),K*_20(1430))
     *     ,2.430,2.430,2.535,3.525         !D_1(2430)0,D_1(2430)+,D_s1(2536)+,h_c(1P)
     *     ,5.723,5.723,5.830,-1.,9.899/    !B_1(H)+,B_1(H)0,B_s1(H)0,B_c1(H)+,h_b(1P)
c     5- meson mass table
      data (ammes(m,5),m=1,15)
     *   /        1.3183,1.3183,1.2818         !a_20,a_2+-,f_1(1285)
     *     ,1.425,1.425,-1.                 !K*_0+-(1430),K*_00(1430)
     *     ,2.460,2.460,2.573,3.556         !D*_2(2460)0,D*_2(2460),D*_s2(2573)+,Xi_c2(1P)
     *     ,5.743,5.743,5.840,-1.,10.023/   !B*_2+,B*_20,B*_s2_B*_c2,Xi_b2(1P)
c     6- meson mass table
      data (ammes(m,6),m=1,15)
     *   /        1.6724,1.6724,1.5050         !pi_20(1670),pi_2+-,f_0(1500)
     *     ,3*-1.
     *     ,4*-1.
     *     ,5*-1./
c     7- meson mass table
      data (ammes(m,7),m=1,15)
     *   /        1.6888,1.6888,-1.            !rho_30(1690),rho_3+-
     *     ,3*-1.
     *     ,4*-1.
     *     ,5*-1./
c     1/2+ baryon mass table
      data ambar0/0.94,.93828,.93957,0.94,-1.,1.1894,1.1925,1.1974
     1     ,1.1156,1.3149,1.3213,3*-1.
     $     ,2.453               !15          sigma_c++!
     $     ,2.454               !            sigma_c+
     $     ,2.452               !            sigma_c0
     $     ,2.286               !            lambda_c+
     2     ,2.468               !19  1340   !Xi'_c+
     $     ,2.471               !20  2340   !Xi'_c0
     $     ,2.695               !21  3340   !omegac0
     $     ,2.471               !22  3240   !Xi_c0
     $     ,2.468               !23  3140   !Xi_c+
     $     ,3.55                !24  1440
     $     ,3.55                !25  2440
     $     ,3.70                !26  3440
     $     ,4*-1./
c     3/2+ baryon mass table
      data ambar1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     1     ,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     2     ,2.519               !15          sigma_c++
     $     ,2.52                !            sigma_c+
     $     ,2.517               !            sigma_c0
     $     ,-1.
     $     ,2.646               !Xi'_c+
     $     ,2.646               !Xi'_c0
     $     ,2.80
     $     ,2*-1.
     $     ,3.75
     $     ,3.75
     3     ,3.90
     $     ,4.80
     $     ,3*-1./
c     entry
      id=idi
      amass=0.
ctp060829      if(iabs(id).eq.30)then
ctp060829        amass=amhdibar
ctp060829        return
ctp060829      endif
      if(idi.eq.0)then
        id=1120                 !for air target
      elseif(abs(idi).ge.1000000000)then
        goto 500                !nucleus
      endif
      if(idi.gt.10000)return
      call idflav(id,ifl1,ifl2,ifl3,jspin,ind)
      if(id.ne.0.and.mod(id,100).eq.0) goto 400
      if(ifl2.eq.0) goto 200
      if(ifl1.eq.0) goto 100
      if(iabs(ifl1).ge.5.or.iabs(ifl2).ge.5.or.iabs(ifl3).ge.5)
     1     goto 300
c          baryons
      ind=ind-109*jspin-36*nmes-nqlep
      ind=ind-11
      amass=(1-jspin)*ambar0(ind)+jspin*ambar1(ind)
      return
c          mesons
100   continue
      if(jspin.gt.7)then     !meson with jspin>7 not defined
        amass=-1.
      else
      ind=ind-36*jspin-nqlep
      ind=ind-11
      amass=ammes(ind,jspin)
      endif
      return
c          quarks and leptons (+deuteron, triton, alpha, Ks and Kl)
200   continue
      amass=amlep(ind)
      return
c          b and t baryons
300   continue
      amass=amlep(iabs(ifl2))+amlep(iabs(ifl3))+1.07+.045*jspin
      amass=amass+amlep(iabs(ifl1))
      return
c          diquarks
400   amass=amlep(iabs(ifl1))+amlep(iabs(ifl2))
      return
c          nuclei
500   nbrpro=mod(abs(id/10000),1000)
      nbrneu=mod(abs(id/10),1000)-nbrpro
      amass=nbrpro*ambar0(2)+nbrneu*ambar0(3)
      return
      end
c-----------------------------------------------------------------------
      subroutine idflav(id,ifl1,ifl2,ifl3,jspin,indx)
c     unpacks the ident code id=+/-ijkl
c
c          mesons--
c          i=0, j<=k, +/- is sign for j
c          id=110 for pi0, id=220 for eta, etc.
c
c          baryons--
c          i<=j<=k in general
c          j<i<k for second state antisymmetric in (i,j), eg. l = 2130
c
c          other--
c          id=1,...,6 for quarks
c          id=9 for gluon
c          id=10 for photon
c          id=11,...,16 for leptons
c          i=17 for deuteron
c          i=18 for triton
c          i=19 for alpha
c          id=20 for ks, id=-20 for kl
c
c          i=21...26 for scalar quarks
c          i=29 for gluino
c
c          i=30 for h-dibaryon
c
c          i=31...36 for scalar leptons
c          i=39 for wino
c          i=40 for zino
c
c          id=80 for w+
c          id=81,...,83 for higgs mesons (h0, H0, A0, H+)
c          id=84,...,87 for excited bosons (Z'0, Z''0, W'+)
c          id=90 for z0
c
c          diquarks--
c          id=+/-ij00, i<j for diquark composed of i,j.
c
c
c          indx is a sequence number used internally
c          (indx=-1 if id doesn't exist)
c
c-----------------------------------------------------------------------
      parameter ( nqlep=41,nmes=2)
      ifl1=0
      ifl2=0
      ifl3=0
      jspin=0
      idabs=iabs(id)
      i=idabs/1000
      j=mod(idabs/100,10)
      k=mod(idabs/10,10)
      jspin=mod(idabs,10)
      if(id.ge.10000) goto 400
      if(id.ne.0.and.mod(id,100).eq.0) goto 300
      if(j.eq.0) goto 200
      if(i.eq.0) goto 100
c          baryons
c          only x,y baryons are qqx, qqy, q=u,d,s.
      ifl1=isign(i,id)
      ifl2=isign(j,id)
      ifl3=isign(k,id)
      if(k.le.6) then
        indx=max0(i-1,j-1)**2+i+max0(i-j,0)+(k-1)*k*(2*k-1)/6
     1  +109*jspin+36*nmes+nqlep+11
      else
        indx=max0(i-1,j-1)**2+i+max0(i-j,0)+9*(k-7)+91
     1  +109*jspin+36*nmes+nqlep+11
      endif
      return
c          mesons
100   continue
      ifl1=0
      ifl2=isign(j,id)
      ifl3=isign(k,-id)
      indx=j+k*(k-1)/2+36*jspin+nqlep
      indx=indx+11
      return
c          quarks, leptons, etc
200   continue
      ifl1=0
      ifl2=0
      ifl3=0
      jspin=0
      indx=idabs
      if(idabs.lt.20) return
c          define indx=20 for ks, indx=21 for kl
      indx=idabs+1
      if(id.eq.20) indx=20
c          indx=nqlep+1,...,nqlep+11 for w+, higgs, z0
      if(idabs.lt.80) return
      indx=nqlep+idabs-79
      return
300   ifl1=isign(i,id)
      ifl2=isign(j,id)
      ifl3=0
      jspin=0
      indx=0
      return
 400  indx=-1
      return
      end
c-----------------------------------------------------------------------
      subroutine Initialize
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer maxo,maxmm,musmm,iimom,i1mom,i2mom,i3mom
      parameter(maxo=2,maxmm=maxo)
      common/cxmom1/musmm,iimom(0:maxo,0:maxo,0:maxo)
     &           ,i1mom(0:maxmm),i2mom(0:maxmm),i3mom(0:maxmm)
      integer maxiz,numiz
      integer ngenmx
      parameter (ngenmx=100)
      double precision cntgen(0:ngenmx,2)
      common/countgen/cntgen
      double precision zamin,zamax,yieldz,yiex,spec
     &,tamin,tamax,ctime,eamin,eamax
      integer maxjz,maxin,maxie,numie
      common/cxtime/ tamin,tamax,ctime
      double precision ramin,ramax,xamin,xamax,yieldr,yieldx
     &,specr,spex,speca,yieldr2,yieldr1
      integer maxir,maxjr,maxiex,maxix,numix,numir,irfirst,modr
     &,iefirst,moden

      parameter(maxiz=10,maxjz=maxiz,maxir=101,maxjr=3
     &          ,maxie=151,maxix=101,maxin=12,maxiex=5 )
      common/cxlimits/zamin,zamax,ramin,ramax
     & ,eamin(2),eamax(2),xamin(maxiex),xamax(maxiex),numix(maxiex)
     & ,numie,numir,irfirst,modr,numiz,iefirst,moden
      common/cxyield/ yieldz(maxin,maxiz),yiex(maxin,maxiz,maxie)
     & ,yieldr(maxin,maxjz,maxir),yieldr1(maxin,maxjz,maxir)
     & ,yieldr2(maxin,maxjz,maxir)
     & ,yieldx(maxiex,maxin,maxjz,maxjr,maxix)
     & ,spec(0:maxmm,maxin,maxie,maxjz),spex(0:maxmm,maxin,maxie,maxjz)
     & ,specr(maxin,maxie,maxjz,maxjr),speca(maxin,maxjz,maxjr)

c Analyze for output

      numiz=0   !number of bins for height and time in histo (updated in corsread)
      musmm=-1 !maxmm
      numie=140          !number of bins for energy in MC analysis


c  analysis

      do i=0,ngenmx
        cntgen(i,1)=0d0
        cntgen(i,2)=0d0
      enddo

      zamin=1.d6        !min height for analysis
      zamax=0.d0        !max height for analysis udpated in corsread
      tamax=1.d5        !max time for analysis
      tamin=1.d3       !min time for analysis
      ctime=(tamax/tamin)**(1.d0/dble(numie))
      eamin(1)=0.0001d0      !minimum energy for leptons in MC analysis in GeV
      eamax(1)=1000.d0      !maximum energy for leptons in MC analysis in GeV (if 0, set to eprima)
      eamin(2)=0.1d0      !minimum energy for hadrons in MC analysis in GeV
      eamax(2)=1.d6      !maximum energy for hadrons in MC analysis in GeV (if 0, set to eprima)
      ramin=0.1d0    !moment analysis parameters
      ramax=10000d0
      numir=50
      irfirst=20
      modr=10
      xamin(1)=-1.d0      !cos(phi)
      xamax(1)=1.d0
      numix(1)=40
      xamin(2)=0d0        !sin2(theta)
      xamax(2)=1d0
      numix(2)=40
      xamin(3)=5000.      !time
      xamax(3)=15000.
      numix(3)=100
      xamin(4)=-3.14159        !delta phi (position/direction)
      xamax(4)=3.14159
      numix(4)=40
      xamin(5)=-10d0        !log10(sin2(theta))
      xamax(5)=0d0
      numix(5)=40
      iefirst=20
      moden=30

      ii=-1
      do n=0,maxo  !order
        do k=n,n  !block
          do i=n-k,0,-1
            ii=ii+1
            i1mom(ii)=k
            i2mom(ii)=i
            i3mom(ii)=n-k-i
            iimom(k,i,n-k-i)=ii
          enddo
        enddo
      enddo


      end


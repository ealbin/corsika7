c=======================================================================
c
c  s u m l i s t n k g i n f o . f
c  -------------------------------
c     sum all `NKG - AVERAGE` tabular quantities from all given 
c     `.lst`-files and calculate XMAX of the air shower 
c     (originally designed for parallel corsika simulations).
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check sumlistnkginfo.f -o sumlistnkginfo
c     f77 -fbounds-check sumlistnkginfo.f -o sumlistnkginfo
c     ifort -C -check bounds sumlistnkginfo.f -o sumlistnkginfo
c execution:
c     ls -1 DAT002345-*.lst > sumlistnkginfo.i002345
c     ./sumlistnkginfo < sumlistnkginfo.i002345
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program sumlistnkginfo

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (namax=50000) 

      character chlstname(namax)*100,chtable(210)*100,czeile*100

      dimension qthickn(210),qheight(210),qelects(210),lstlen(namax)

      parameter (ni=1000,nc=3,nn=12)       
      dimension zyg(nn),zx(nn),zy(nn),zv(ni),zw(ni),zh(nn)  

      double precision aatm(5),batm(5),catm(5),thickgr

      common /atmos/aatm,batm,catm

      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data zyg/nn*1./
      data qelects/210*0./

      cpi180 = atan(1.d0) / 45.d0
      c180pi = 45.d0 / atan(1.d0)

c - - - - - - read all Job files and keep names in a character array - -
      do  lst=1,namax
         read(*,'(a)',err=101,end=102) chlstname(lst)
         jl = 100 + 1 
  100    continue 
         jl = jl - 1
         if ( chlstname(lst)(jl:jl) .eq. ' ' ) goto 100
         lstlen(lst) = jl       
      enddo
      write(*,*) ' w a r n i n g    possibly >',namax,' file names'
      goto 102
  101 continue
      write(*,*) ' e r r o r   reading chlstname ',chlstname(lst)
  102 continue  
      lst = lst - 1
      read(chlstname(1)(4:9),'(i6)') lrunnr
      lflagthe = 0
      lflagphi = 0
      lflagegy = 0

c - - - - - - read first lst file - - - - - - - - - - - - - - - - - - -
      open(unit=1,file=chlstname(1)(1:lstlen(1)),
     +            form='formatted',status='old')
      do  is=1,9876
         read(1,'(a)',end=122,err=121) czeile
         if ( index(czeile,'PRIMARY ENERGY') .gt. 0 ) then
            if ( index(czeile,'STCKIN:') .gt. 0 ) then
               read(czeile(57:71),'(e15.7)') energy
               lflagegy = 1
            elseif ( index(czeile,'FIXED') .gt. 0 ) then
               read(czeile(39:48),'(f10.2)') energy
               lflagegy = 1
            elseif ( index(czeile,' ENERGY = ') .gt. 0 ) then
               read(czeile(18:30),'(f13.0)') energy
               lflagegy = 1
            endif
         endif
         if ( index(czeile,'THETA OF INCIDENCE') .gt. 0 ) then
            write(*,'(a)') czeile(1:74) 
            if ( index(czeile,'FIXED') .gt. 0 ) then
               read(czeile(35:42),'(f8.2)') thedeg
               lflagthe = 1
            elseif ( index(czeile,'CHOSEN') .gt. 0 ) then 
               read(czeile(35:42),'(f8.2)') thedeg
               read(czeile(48:55),'(f8.2)') thetb 
            endif
         endif
         if ( index(czeile,'PHI   OF INCIDENCE') .gt. 0 ) then
            write(*,'(a)') czeile(1:74) 
            if ( index(czeile,'FIXED') .gt. 0 ) then
               read(czeile(35:42),'(f8.2)') phideg
               lflagphi = 1
            elseif ( index(czeile,'CHOSEN') .gt. 0 ) then 
               read(czeile(35:42),'(f8.2)') phideg
               read(czeile(48:55),'(f8.2)') phib    
            endif
         endif
         if ( index(czeile,'PRIMARY ANGLES') .gt. 0 ) then
           if ( lflagthe .eq. 0 .and. lflagphi .eq. 0 ) then
             it = index(czeile,'THETA =')
             ip = index(czeile,'PHI =')
             read(czeile(it+7:it+13),'(f7.4)') therad
             read(czeile(ip+5:ip+12),'(f8.4)') phirad
             thedeg = c180pi * therad
             phideg = c180pi * phirad 
           endif 
         endif 
         if ( index(czeile,'NKG - AVERAGE') .gt. 0 ) then
            ! - - - - - - - nkg tabular reached:
            read(1,'(a25)',end=122,err=121) czeile
            read(1,'(a25)',end=122,err=121) czeile
            read(1,'(a25)',end=122,err=121) czeile
            read(1,'(a25)',end=122,err=121) czeile
            jl = 0
  106       continue
            jl = jl + 1
            read(1,'(a50)',end=122,err=121) chtable(jl)
            if ( chtable(jl)(1:5) .ne. '     ' ) goto 106
            jl = jl - 1
            ! - - - - - - - - - set first quantities:
            do  il=1,jl
               read(chtable(il),'(i5,f10.1,f12.1,f16.1)')
     +            ih,qthickn(il),qheight(il),qelects(il)
             ! write(*,'(3f20.2)')
     +       !    qthickn(il),qheight(il),qelects(il)
            enddo
            goto 108
         endif
      enddo
  108 continue
      close(unit=1)  
      write(*,'(5x,''run='',i7.6,8x,''theta='',f6.2,
     +          8x,''phi='',f8.2)') lrunnr,thedeg,phideg

c - - - - - - summ all NKG - AVERAGE quantities - - - - - - - - - - - -
      do  ll=2,lst
         open(unit=1,file=chlstname(ll)(1:lstlen(ll)),
     +               form='formatted',status='old')
         do  is=1,9876  
            read(1,'(a)',end=114,err=113) czeile
            if ( index(czeile,' NKG - AVERAGE ') .gt. 0 ) then
            ! - - - - - - - nkg tabular reached:
               read(1,'(a25)',end=114,err=113) czeile
               read(1,'(a25)',end=114,err=113) czeile
               read(1,'(a25)',end=114,err=113) czeile
               read(1,'(a25)',end=114,err=113) czeile
               jl = 0
  112          continue         
               jl = jl + 1    
               read(1,'(a50)',end=114,err=113) chtable(jl)
               if ( chtable(jl)(1:5) .ne. '     ' ) goto 112
               jl = jl - 1
               ! - - - - - - - - - sum up quantities: 
               do  il=1,jl
                  read(chtable(il),'(27x,f16.1)') elects 
                  qelects(il) = qelects(il) + elects
               enddo
               goto 115
            endif 
         enddo
  113    continue
         write(*,*) ' e r r o r   exit of',ll,' lst file ',chlstname(ll)
         goto 115 
  114    continue
         write(*,*) ' end-of-file exit of',ll,' lst file ',chlstname(ll)
  115    continue
         close(unit=1) 
      enddo
      write(*,'(5x,''lst='',i7,8x,''energy='',1p,e12.4,'' GeV'')')
     +   lst,energy
      do  il=1,jl 
      !  write(*,'(3f20.2)') qthickn(il),qheight(il),qelects(il)
      enddo
 
c - - - - - - search maximum between second and last quantities:
      np = jl - 1  
      do   ii=1,np
         zx(ii) = qthickn(ii+1)
         zh(ii) = qheight(ii+1)
         zy(ii) = qelects(ii+1)
      enddo 
      mp = np
      mv = 100  
      mq = mv  
      call asplin(mp,mq,zx,zy,zyg,zv,zw)            
      lflagmax = 0
      do   ii=1,mq-1
         if ( mod(ii,mv) .eq. 1 ) 
     +       write(*,'(3f20.2)') zv(ii),heightcm(zv(ii))*1.d-2,zw(ii)
         if ( lflagmax .eq. 0 .and. zw(ii+1) .lt. zw(ii) ) then
            write(*,'(3f20.2)') zv(ii),heightcm(zv(ii))*1.d-2,zw(ii) 
            lflagmax = 1
         endif 
      enddo
      goto 129

c - - - - - - end-of-file and error exits:
  121 continue
      write(*,*) ' e r r o r   exit of first lst file ',chlstname(1)   
      goto 129
  122 continue
      write(*,*) ' end-of-file exit of first lst file ',chlstname(1)
      goto 129
  123 continue
      write(*,*) ' e r r o r   exit of',ll,' lst file ',chlstname(ll)
      goto 129
  124 continue
      write(*,*) ' end-of-file exit of',ll,' lst file ',chlstname(ll)

c - - - - - - end of program sumlistnkginfo.
  129 continue
      stop
      end
c=======================================================================
c   es werden die koeffizienten a(i),b(i),c(i),d(i),i=0(1)n-1 einer
c   glaettenden natuerlichen kubischen splinefunktion s berechnet.
c   die splinefunktion s(x) wird in der form a(i) + b(i)*(x-xn(i)) +
c   + c(i)*(x-xn(i))**2 + d(i)*(x-xn(i))**3 fuer x element aus dem
c   intervall (xn(i),xn(i+1)),i=0,1,2,..,n-1 dargestellt.
c-voraussetzungen :   n > 4
c                     xn(i) < xn(i+1) , i = 0,1,2,..,n-1
c                     yw(i) > 0.0     , i = 0,1,2,..,n
c-eingabeparameter:
c     n  =  nummer des letzten knotens
c     xn(0:n) = knoten , i = 0,1,2,..,n
c     yn(0:n) = messwert an der stelle xn(i)
c     yw(0:n) = gewicht zum messwert yn(i)
c-hilfsfeld:
c     xhi :  1-dim feld (1:7*n-3)
c-ausgabeparameter:
c     xa(0;n)    die feldelemente 0 bis n-1 sind die
c     xb(0;n)    koeffizienten der splinefunktion s.
c     xc(0;n)    a(n),b(n),c(n),d(n) sind hilfsspeicher.
c     xd(0:n)
c     nfehl = 0, alles klar;
c           =-1, voraussetzungen verletzt;
c           = k, im k-ten eliminationsschritt division durch null.
c     ! nd = min(ni-2,nj-2)/np
c-----------------------------------------------------------------------
      subroutine asplin(np,nz,zx,zy,zyg,zv,zw)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (ni=1000,ipr=0)
      dimension zx(np),zy(np),zyg(np),zv(nz*np),zw(nz*np)
      dimension xn(0:ni),yn(0:ni),yw(0:ni),xhi(7*ni)
      dimension xa(0:ni),xb(0:ni),xc(0:ni),xd(0:ni)
      nd = 100
      nj = np*nd
      nz = nj
      nn = np-1
      nfehl = -1
      if (ni.lt.np-1) nfehl = -9
      do  43  n=1,np
         zv(n)=zx(n)
         zw(n)=zy(n)
         xn(n-1)=dble(zx(n))
         yn(n-1)=dble(zy(n))
         yw(n-1)=dble(zyg(n))
   43 continue
      if (nn.lt.5) then
         write(*,'(/11x,''zu wenig datenpunkte vorhanden.'')')
         goto 49
         endif
      do  44  n=0,nn-1
         if (xn(n).ge.xn(n+1)) then
            write(*,'(/11x,''datenpunkte falsch sortiert.'')')
            goto 49
            endif
   44 continue
      do  45  n=0,nn
         if (yw(n).le.0.0) then
            write(*,'(/11x,''gewicht eines datenpunktes <= 0.0'')')
            goto 49
            endif
   45 continue
c - - - - - berechnung der ausgleichs-spline-koeffizienten.
      call aspliz
     +   (nn,xn,yn,yw,xa,xb,xc,xd,xhi(1),xhi(nn),xhi(2*nn-1),
     +   xhi(3*nn-2),xhi(4*nn-3),xhi(5*nn-4),xhi(6*nn-4),nfehl)
      if (ipr.ge.2) then
         write(*,'(//9x,''i'',9x,''xa'',15x,''xb'',15x,''xc'',
     +   15x,''xd'')')
         do  46  i=0,nn
         write(*,'(6x,i5,4f17.12)') i,xa(i),xb(i),xc(i),xd(i)
   46    continue
         endif
      if (nd.le.0) nd = 1
      if (nd.ge.2) then
c - - - - - anzahl der punkte um den faktor nd>1 erhoehen.
         ! nz = 1 + nj
         if (ipr.ge.1) write(*,'(11x,''es gibt'',i4,
     +   '' punkte beim ausgleichsspline:''/)') 1+nj
         do  48  nb=0,nn-1
            na = 1 + nd * nb
            zxdn = (zx(nb+2)-zx(nb+1)) / nd
            if (nb.eq.nn-1) nd=nd+1
            do  47  n=0,nd-1
               zxd = zxdn * n
               zv(na+n) = zx(nb+1) + zxd
               zw(na+n) = ((xd(nb) * zxd + xc(nb)) *
     +               zxd + xb(nb)) * zxd + xa(nb)
               if (ipr.ge.1) write(*,*)
     +         na+n,zv(na+n),zw(na+n)
   47       continue
   48    continue
         endif
   49 continue
      if (nfehl.ne.0) write(*,'(6x,''nfehl ='',i3)') nfehl
      return
      end
c=======================================================================
      subroutine aspliz (nn,xn,yn,yw,xa,xb,xc,xd,xu1,xu2,xo1,xo2,xr,h1,
     +                   h2,nfehl)
      implicit doubleprecision(a-h,o-y)
      dimension xn(0:nn),yn(0:nn),yw(0:nn)
      dimension xa(0:nn),xb(0:nn),xc(0:nn),xd(0:nn)
      dimension xu1(1:nn-1),xu2(1:nn-1),xo1(1:nn-1),xo2(1:nn-1)
      dimension xr(1:nn-1),h1(0:nn-1),h2(0:nn-1)
c-berechnung der xc(i); die felder a,b,c,d werden zunaechst
c mit diesen hilfsgroessen besetzt.
      nk = nn-1
      do 10 i=0,nn-1
        h2(i) = xn(i+1)-xn(i)
        h1(i) = 1.d0/h2(i)
        xc(i) = h1(i)*h1(i)
        xb(i) = 6.d0/yw(i)
   10 continue
      xb(nn) = 6.d0/yw(nn)
      do 15 i=1,nn-1
        xd(i) = h1(i-1)+h1(i)
   15 continue
c-berechnung der matrixelemente und der rechten seite des fuenfdiago-
c nalen gleichungssystems zur bestimmung der xc(i) fuer i=1 bis nn-1
      do 20 i=3,nn-1
        xu2(i) = xb(i-1)*h1(i-2)*h1(i-1)
        xo2(i-2) = xu2(i)
   20 continue
      do 25 i=2,nn-1
        xu1(i) = h2(i-1)-xb(i-1)*h1(i-1)*xd(i-1)-xb(i)*h1(i-1)*xd(i)
        xo1(i-1) = xu1(i)
   25 continue
      do 30 i=1,nn-1
        xa(i) = 2.d0*(h2(i-1)+h2(i))+xb(i-1)*xc(i-1)+xb(i)*xd(i)*xd(i)+
     +          xb(i+1)*xc(i)
   30 continue
      xh1 = (yn(1)-yn(0))*h1(0)
      do 35 i=1,nn-1
        xh2 = (yn(i+1)-yn(i))*h1(i)
        xr(i) = (xh2-xh1)*3.d0
        xh1 = xh2
   35 continue
c-besetzen der koeffizienten xc(0),xc(nn)
      xc(0) = 0.d0
      xc(nn) = 0.d0
c-berechnung der loesungen xc(i),i=1 bis nk
c        des fuenfdiagonalen gleichungssystems
      nfehl = -1
      if (nk.lt.4) goto 70
c-zerlegung der matrix a
      xu2(1) = 0.d0
      xu2(2) = 0.d0
      xu1(1) = 0.d0
      xo1(nk) = 0.d0
      xo2(nk-1) = 0.d0
      xo2(nk) = 0.d0
      if (xa(1).eq.0.0) then
        nfehl = 1
        goto 45
      endif
      xh = 1.d0/xa(1)
      xo1(1) = xo1(1)*xh
      xo2(1) = xo2(1)*xh
      xa(2) = xa(2)-xu1(2)*xo1(1)
      xh = 1.d0/xa(2)
      xo1(2) = (xo1(2)-xu1(2)*xo2(1))*xh
      xo2(2) = xo2(2)*xh
      do 40 i=3,nk-2
        xu1(i) = xu1(i)-xu2(i)*xo1(i-2)
        xa(i) = xa(i)-xu2(i)*xo2(i-2)-xu1(i)*xo1(i-1)
        if (xa(i).eq.0.0) then
          nfehl = nn-1
          goto 45
        endif
        xh = 1.d0/xa(i)
        xo1(i) = (xo1(i)-xu1(i)*xo2(i-1))*xh
        xo2(i) = xo2(i)*xh
   40 continue
      xu1(nk-1) = xu1(nk-1)-xu2(nk-1)*xo1(nk-3)
      xa(nk-1) = xa(nk-1)-xu2(nk-1)*xo2(nk-3)-xu1(nk-1)*xo1(nk-2)
      if (xa(nk-1).eq.0.0) then
        nfehl = nk-1
        goto 45
      endif
      xh = 1.d0/xa(nk-1)
      xo1(nk-1) = (xo1(nk-1)-xu1(nk-1)*xo2(nk-2))*xh
      xu1(nk) = xu1(nk)-xu2(nk)*xo1(nk-2)
      xa(nk) = xa(nk)-xu2(nk)*xo2(nk-2)-xu1(nk)*xo1(nk-1)
      if (xa(nk).eq.0.0) then
        nfehl = nk
        goto 45
      endif
      nfehl = 0
   45 continue
c-falls nfehl=0: vorwaerts - und rueckwaerts-elimination
      if (nfehl.eq.0) then
        xr(1) = xr(1)/xa(1)
        xr(2) = (xr(2)-xu1(2)*xr(1))/xa(2)
        do 50 i=3,nk
          xr(i) = (xr(i)-xu2(i)*xr(i-2)-xu1(i)*xr(i-1))/xa(i)
   50   continue
        xc(nk) = xr(nk)
        xc(nk-1) = xr(nk-1)-xo1(nk-1)*xc(nk)
        do 55 i=nk-2,1,-1
          xc(i) = xr(i)-xo1(i)*xc(i+1)-xo2(i)*xc(i+2)
   55   continue
      endif
      if (nfehl.ne.0) goto 70
c-berechnung xa(i),i=0 bis nn bzw. xb(i) und xd(i),i=0 bis nn-1
      xa(0) = yn(0)+xb(0)/3.d0*h1(0)*(xc(0)-xc(1))
      do 60 i=1,nn-1
        xa(i) = yn(i)-xb(i)/3.d0*(xc(i-1)*h1(i-1)-xd(i)*xc(i)+xc(i+1)*
     +          h1(i))
   60 continue
      xa(nn) = yn(nn)-xb(nn)/3.d0*h1(nn-1)*(xc(nn-1)-xc(nn))
      do 65 i=0,nn-1
        xb(i) = h1(i)*(xa(i+1)-xa(i))-h2(i)/3.d0*(xc(i+1)+2.d0*xc(i))
        xd(i) = h1(i)/3.d0*(xc(i+1)-xc(i))
   65 continue
   70 continue
      if (nfehl.ne.0) write(*,'(/6x,''nfehl ='',i3)') nfehl
      return
      end
c=======================================================================

      double precision function heightcm( arg )

c-----------------------------------------------------------------------
c     height above sea level as function of gramms per cm^2
c-----------------------------------------------------------------------

      implicit double precision (a-h)
      double precision aatm(5),batm(5),catm(5),arg
      common /atmos/aatm,batm,catm

      if ( arg .gt. 631.1d0 ) then
        heightcm = catm(1) * log( batm(1) / (arg - aatm(1)) )
      elseif ( arg .gt. 271.7d0 ) then
        heightcm = catm(2) * log( batm(2) / (arg - aatm(2)) )
      elseif ( arg .gt. 3.0395d0 ) then
        heightcm = catm(3) * log( batm(3) / (arg - aatm(3)) )
      elseif ( arg .gt. 1.28292d-3 ) then
        heightcm = catm(4) * log( batm(4) / (arg - aatm(4)) )
      else
        heightcm = (aatm(5) - arg) / catm(5)
      endif

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
         thickgr = aatm(1) + batm(1) * exp( -arg / catm(1) )
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg / catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg / catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg / catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif

      return
      end
c=======================================================================
c    run= 002141        theta= 23.00        phi=   26.56
c    lst=         3014
c             200.00            11954.18         26039350.25
c             300.00             9347.05        240539017.23
c             400.00             7364.94        934583066.61
c             500.00             5748.86       2097964775.50
c             600.00             4384.47       3174452035.83
c             694.00             3263.19       3551052961.13
c             700.00             3195.68       3549236855.58
c             800.00             2133.13       3119941043.40
c             870.00             1451.62       2534534337.60
c=======================================================================

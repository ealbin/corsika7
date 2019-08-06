C=======================================================================

      program plottracks3c

c-----------------------------------------------------------------------
c  The CORSIKA PLOTSH Option writes track elements which are 
c  converted  into a plot with this program plottrack.
C  a figure in ppm format is produced which afterwards 
c  can be converted by the 'xv' procedure of the unix stations 
c  into a 'gif' or a 'jpg' picture.
c-----------------------------------------------------------------------
c  compilation:
C         gfortran -fbounds-check plottracks3c.f -o plottracks3c
C         f77 -fbounds-check -m32 plottracks3c.f -o plottracks3c
C         ifort -C plottracks3c.f -o plottracks3c
c-----------------------------------------------------------------------
c  execution:
c  provide the DAT<nnnnnn>.track_xx files in the directory of execution.
c         ./plottracks3c 
c
c-----------------------------------------------------------------------
c  Autor:        Johannes Knapp   24.11.1997
c  modifications: Dieter Heck       1.03.2004
c  modifications: Tanguy Pierog    15.02.2006
c-----------------------------------------------------------------------
      include 'work.inc'

      integer i,k,l
      real run
c-----------------------------------------------------------------------
      write(*,*) 
      write(*,*) 'produce a pixel graphics of a shower plot'
      write(*,*) '========================================='
      write(*,*) 
      write(*,*) 'This program will read :'
      write(*,*) 'DAT<run>.track_em,'
      write(*,*) 'DAT<run>.track_mu, and'
      write(*,*) 'DAT<run>.track_hd'
      write(*,*) '<run> being an integer between 0 and 999999 (input)'
      write(*,*) 'and being the same for all 3 files'
      write(*,*) 
      write(*,*) 'Outputs : track<run>.em.ppm, track<run>.mu.ppm'
      write(*,*) 'track<run>.had.ppm and the sum of all'
      write(*,*) 'track<run>.all.ppm '
      write(*,*) '(can be opened with xview (xv)) '
      write(*,*) 
      write(*,*) 'And now :'

      write(*,*) 'which projection ?  (1,2,3 = x-z,y-z,x-y)'
      read(*,'(a)') ich
      if ( ich .eq. ' ' .or. index(ich,'1') .ne. 0 ) then
        ipr = 1
      elseif ( index(ich,'2') .ne. 0 ) then
        ipr = 2
      else 
        ipr = 3
      endif
      write(*,*) 'ipr =',ipr

      write(*,*) 'radius range ?  (in km)'
      read(*,*) radi
      write(*,*) 'r < ',radi,' km'

      write(*,*) 'background colour ? (b/w)'
      read(*,'(a)') ich
      if ( ich .eq. 'b' .or. ich .eq. 'B' ) then
        write(*,*) 'black background'
        higr = 0
      else
        write(*,*) 'white background'
        higr = 255
      endif

      write(*,*) 'which energy cuts ? (em,mu,had) '
      read(*,'(a)') ich
      if ( ich .eq. ' ' ) then
        cutem = 0.
        cutmu = 0.
        cuthad = 0.
      else
        read(ich,*) cutem,cutmu,cuthad 
      endif
      write(*,*) 'cuts em,mu,had =',cutem,cutmu,cuthad

      write(*,*) 'which run ?'
      read(*,*) run
      if(run.lt.1e6)then
        irun=int(abs(run))
        write(*,*) 'irun =',irun
      else
        stop'Run number too big (>999999)'
      endif

c  fill the field arrays
      dsn = 'DAT000000.track'
      write(dsn(4:9),'(i6)') irun
      do  l = 4, 9
        if ( dsn(l:l) .eq. ' ' ) dsn(l:l) = '0'
      enddo
      dsn2 = 'track000000'
      write(dsn2(6:11),'(i6)') irun
      do  l = 6, 11
        if ( dsn2(l:l) .eq. ' ' ) dsn2(l:l) = '0'
      enddo

c  plot size
      if ( ipr .eq. 3 ) then
        nx = 500
        ny = 500
        xmin = -radi
        xmax =  radi
        ymin = -radi
        ymax =  radi
      else
        nx = 500
        ny = 1000
        xmin = -radi
        xmax =  radi
c  height : 28 km
        ymin =  0.110
        ymax = 28.
      endif       

      ncol = 255

c  pixels
      ixmin = 1
      ixmax = nx
      iymin = 1
      iymax = ny

      xconst = (ixmax-ixmin)/(xmax-xmin)
      yconst = (iymax-iymin)/(ymax-ymin)

c  electromagnetic
      do i=1,nx
      do k=1,ny
        ifeld(i,k) = 0
      enddo
      enddo
      cut = cutem
      icol = 2
      dsn(16:22) = '_em    '
      call getit
      dsn2(12:19) = '.em.raw '
      call rawenc
      dsn2(12:19) = '.em.ppm '
      call ppmenc

c  muonic
      do i=1,nx
      do k=1,ny
        ifeld(i,k) = 0
      enddo
      enddo
      cut = cutmu
      icol = 3
      dsn(16:22) = '_mu    '
      call getit
      dsn2(12:19) = '.mu.raw '
      call rawenc
      dsn2(12:19) = '.mu.ppm '
      call ppmenc

c  hadronic
      do i=1,nx
      do k=1,ny
        ifeld(i,k) = 0
      enddo
      enddo
      cut = cuthad
      icol = 4
      dsn(16:22) = '_hd    '
      call getit
      dsn2(12:20) = '.had.raw '
      call rawenc
      dsn2(12:20) = '.had.ppm '
      call ppmenc

c  overlay the plots
      call mergepl

      stop
      end
C=======================================================================

      subroutine getit

c-----------------------------------------------------------------------
c  reads the files produced by CORSIKA 
c-----------------------------------------------------------------------
      include 'work.inc'

cdh   real dd,ee,x1,x2,y1,y2,z1,z2
      real dd,ee,t1,t2,x1,x2,y1,y2,z1,z2
      integer nplem,ntot
c-----------------------------------------------------------------------
      open(unit=55,file=dsn,status='old',form='unformatted')
      nplem = 0
      ntot = 0

      write(*,*) ' '
      write(*,*) '.....  read tracks from ',dsn

 1    continue
      ntot = ntot + 1
      if ( ntot .gt. 51000000 ) goto 999 
cdh   read(55,end=999) dd,ee,x1,y1,z1,x2,y2,z2
      read(55,end=999) dd,ee,x1,y1,z1,t1,x2,y2,z2,t2
     
      if ( ee .ge. cut ) then
        nplem = nplem + 1
        if     ( ipr .eq. 1 ) then
          call linpl(x1*1.e-5,z1*1.e-5,x2*1.e-5,z2*1.e-5)
        elseif ( ipr .eq. 2 ) then
          call linpl(y1*1.e-5,z1*1.e-5,y2*1.e-5,z2*1.e-5)
        else
          call linpl(x1*1.e-5,y1*1.e-5,x2*1.e-5,y2*1.e-5)
        endif
      endif

      goto 1

 999  continue
      write(*,*) nplem,' tracks from ',dsn,'  above ',cut
      close(55)

      return
      end
C=======================================================================

      subroutine linpl(x1,y1,x2,y2)

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      include 'work.inc'

      real x1,x2,y1,y2,xx1,yy1
      integer ixx1,ixx2,iyy1,iyy2,ix1,ix2,iy1,iy2,iyy,ixx
      integer i
c-----------------------------------------------------------------------
      ixx1 = nint(ixmin + (x1-xmin) * xconst) 
      ixx2 = nint(ixmin + (x2-xmin) * xconst)
      iyy1 = nint(iymin + (y1-ymin) * yconst)
      iyy2 = nint(iymin + (y2-ymin) * yconst)
      ix1 = min(ixx1,ixx2)
      ix2 = max(ixx1,ixx2)
      iy1 = min(iyy1,iyy2)
      iy2 = max(iyy1,iyy2)
cc      write(*,*) 'linpl :  x1, y1, x2, y2', x1,y1,x2,y2
cc      write(*,*) 'linpl : ix1,iy1,ix2,iy2', ix1,iy1,ix2,iy2
      
      if ( ix1 .eq. ix2 ) then
        if ( iy1 .eq. iy2 ) then
c  only one point
            if ( ix1 .ge. ixmin ) then
              if ( ix1 .le. ixmax ) then
                if ( iy1 .ge. iymin ) then
                  if ( iy1 .le. iymax ) then
                    ifeld(ix1,iy1) = ifeld(ix1,iy1) + 1         
cc                    ifeld(ix1,iy1) = ifa
                  endif
                endif
              endif
            endif
cc            write(*,*) 'linpl point: ',ix1,iy1
        else
c  vertical line
          do i=iy1,iy2
            if ( ix1 .ge. ixmin ) then
              if ( ix1 .le. ixmax ) then
                if ( i .ge. iymin ) then
                  if ( i .le. iymax ) then
                    ifeld(ix1,i) = ifeld(ix1,i) + 1
cc                    ifeld(ix1,i) = ifa
                  endif
                endif
              endif
            endif
cc            write(*,*) 'linpl vert: ',ix1,i
          enddo
        endif
      elseif ( iy1 .eq. iy2 ) then
c  horizontal line
        do i=ix1,ix2
            if ( i .ge. ixmin ) then
              if ( i .le. ixmax ) then
                if ( iy1 .ge. iymin ) then
                  if ( iy1 .le. iymax ) then
                    ifeld(i,iy1) = ifeld(i,iy1) + 1
cc                    ifeld(i,iy1) = ifa
                  endif
                endif
              endif
            endif
cc            write(*,*) 'linpl hori: ',i,iy1
        enddo
      else
c  skew lines
c  along the x axis
        if ( abs(ix2-ix1) .gt. abs(iy2-iy1) ) then
          yyconst = (y2-y1)/(x2-x1)
          do i=ix1,ix2
            xx1 = (i-ixmin)/xconst + xmin
            yy1 = y1 + (xx1-x1) * yyconst
            iyy = nint(iymin + (yy1-ymin) * yconst)
cc            write(*,*) 'linpl entl x: ',xx1,yy1,i,iyy
            if ( i .ge. ixmin ) then
              if ( i .le. ixmax ) then
                if ( iyy .ge. iymin ) then
                  if ( iyy .le. iymax ) then
                    ifeld(i,iyy) = ifeld(i,iyy) + 1
cc                    ifeld(i,iyy) = ifa
                  endif
                endif
              endif
            endif
          enddo
        else
c  along th y axis
          xxconst = (x2-x1)/(y2-y1)
          do i=iy1,iy2
            yy1 = (i-iymin)/yconst + ymin
            xx1 = x1 + (yy1-y1) * xxconst
            ixx = nint(ixmin + (xx1-xmin) * xconst)
cc            write(*,*) 'linpl entl y: ',xx1,yy1,ixx,i
            if ( ixx .ge. ixmin ) then
              if ( ixx .le. ixmax ) then
                if ( i .ge. iymin ) then
                  if ( i .le. iymax ) then
                    ifeld(ixx,i) = ifeld(ixx,i) + 1
cc                    ifeld(ixx,i) = ifa
                  endif
                endif
              endif
            endif
          enddo
        endif
      endif

      return
      end
C=======================================================================

      subroutine ppmenc

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      include 'work.inc'

      integer r,g,b
      integer i,k
c-----------------------------------------------------------------------
      write(*,*) '.....  store picture in ',dsn2

      open(unit=11,file=dsn2,form='formatted',
     *     status='unknown')
      write(11,100) nx, ny, ncol
 100  format('P3',/,'# CREATOR: J. Knapp   24.11.1997',/,
     *       '# ',/,
     * 2i5,/,i3)

      do k=1,ny
        do i=1,nx
          call coldef(ifeld(i,k),r,g,b)
          write(11,101) r,g,b
 101      format(i3,2i4)
        enddo
      enddo

      close(11)

      return
      end
C=======================================================================

      subroutine rawenc

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      include 'work.inc'

c-----------------------------------------------------------------------
      write(*,*) '.....  store raw picture in ',dsn2

      open(unit=11,file=dsn2,form='unformatted',
     *     status='unknown')

      write(11) nx,ny,ncol
      write(11) ifeld

      close(11)

      return
      end
C=======================================================================

      subroutine coldef(iii,rr,gg,bb)

c-----------------------------------------------------------------------
c  defines the colours
c-----------------------------------------------------------------------
      include 'work.inc'

      integer iii,iiii,ict
      integer rr,gg,bb,f1(0:255),f2(0:255),f3(0:255)  
      logical first
      save
      data first /.true./
c-----------------------------------------------------------------------
      if ( first ) then
        first = .false.

c  define the background colour 
        f1(0) = higr
        f2(0) = higr
        f3(0) = higr
c  move to yellow
        do ict = 1,255
          f1(ict) = 255
          f2(ict) = ict
          f3(ict) = ict
        enddo
      endif

      iiii = min(iii,255)

      if ( icol .eq. 2 ) then
        rr = f1(iiii)
        gg = f2(iiii)
        bb = f3(iiii)
      elseif ( icol .eq. 3 ) then
        gg = f1(iiii)
        rr = f2(iiii)
        bb = f3(iiii)
      elseif ( icol .eq. 4 ) then
        bb = f1(iiii)
        rr = f2(iiii)
        gg = f3(iiii)
      else
        write(*,*) 'illegal colour'
        stop
      endif
      
      return
      end
C=======================================================================

      subroutine mergepl

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      include 'work.inc'

      integer i,k,rgb(3,3),ibl,iunit
      character*80 line
c-----------------------------------------------------------------------
      write(*,*) ' '
      write(*,*) '.....  merge pixel graphics of a shower plot'

c  open datasets
      dsn2(12:19) = '.em.ppm '
      open(unit=41,file=dsn2,form='formatted',status='old')
      dsn2(12:19) = '.mu.ppm '
      open(unit=42,file=dsn2,form='formatted',status='old')
      dsn2(12:20) = '.had.ppm '
      open(unit=43,file=dsn2,form='formatted',status='old')

      dsn2(12:20) = '.all.ppm '
      open(unit=44,file=dsn2,form='formatted',status='unknown')
      write(*,*) '       store it on ',dsn2

c  3 dummy lines
      read(41,109) line
 109  format(a)
      write(44,109) line(1:20)
      read(42,109) line
      read(43,109) line

      read(41,'(a)') line
      write(44,109)  line(1:20)
      read(42,'(a)') line
      read(43,'(a)') line

      read(41,'(a)') line
      write(44,109)  line(1:20)
      read(42,'(a)') line
      read(43,'(a)') line

c  nx ny
      do i=1,3
        iunit = 40+i
        read(iunit,*) rgb(i,1),rgb(i,2)
      enddo
      nx = rgb(1,1) 
      ny = rgb(1,2)
      if ( nx .ne. rgb(2,1) .or. nx .ne. rgb(3,1) .or.
     *     ny .ne. rgb(2,2) .or. ny .ne. rgb(3,2)      ) then
        write(*,*) 'incompatible file formats :'
        do i=1,3
          write(*,*) 'file ',i,': nx,ny =',rgb(i,1),rgb(i,2)
        enddo
        stop
      endif
      write(44,111) nx,ny
 111  format(2i5)

c  ncol
      do i=1,3
        iunit = 40+i
        read(iunit,*) rgb(i,1)
      enddo
      ncol = rgb(1,1) 
      if ( ncol .ne. rgb(2,1) .or. ncol .ne. rgb(3,1) ) then
        write(*,*) 'incompatible colour map :'
        do i=1,3
          write(*,*) 'file ',i,': ncol =',rgb(i,1)
        enddo
        stop
      endif
      write(44,112) ncol
 112  format(i3)

c  pixel data
 1    continue
c  read em mu and had data for this pixel
      do i=1,3
        iunit = 40+i
        read(iunit,113,end=999) rgb(i,1),rgb(i,2),rgb(i,3)
 113    format(i3,2i4)
      enddo
c------------------------------
c  muons and hadrons cover the electrons 
      if ( rgb(3,1) .eq. higr .and. 
     *     rgb(3,2) .eq. higr .and.
     *     rgb(3,3) .eq. higr        ) then
        if ( rgb(2,1) .eq. higr .and. 
     *       rgb(2,2) .eq. higr .and.
     *       rgb(2,3) .eq. higr      ) then
        else
          do k=1,3
            rgb(1,k) = rgb(2,k)
          enddo
        endif
      else
        do k=1,3
          rgb(1,k) = rgb(3,k)
        enddo
      endif
c------------------------------
c  mix the colours 
cc      do i=1,3
cc        rgb(1,i) = min(rgb(1,i) + rgb(2,i) + rgb(3,i),255)
cc      enddo
c------------------------------



      write(44,113) rgb(1,1),rgb(1,2),rgb(1,3)

      goto 1

 999  continue
      write(*,*) 'eof reached'
      close (41)
      close (42)
      close (43)
      close (44)
      
      return
      end

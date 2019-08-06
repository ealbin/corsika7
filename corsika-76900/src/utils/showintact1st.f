c=======================================================================
c
c  s h o w i n t a c t 1 s t . f
c  ----------------------------- 
c  read first particle data file of parallel corsika simulations 
c  and print out tabular of current heights of first interaction 
c  of each air shower simulation, i.e. additional tabular info to
c  showparallel.???-work-jobinfos by script showparallel.sh; 
c-----------------------------------------------------------------------
c compile+link:
c     gfortran -O0 -fbounds-check showintact1st.f -o showintact1st 
c     ifort -C -O0 -check bounds showintact1st.f -o showintact1st
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c running:
c     ./showintact1st.sh
c     ls -1 DAT000182-* | grep t -v | grep n -v > showintact1st.i000182
c     ./showintact1st < showintact1st.i000182 > showintact1st.out000182   
c-----------------------------------------------------------------------
c     input-files:
c          unit=*: name of particle data file.
c          unit=3: showintact1st.partfile 
c     output-file: 
c          unit=7: showintact1st.h1stme
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
  
      program showintact1st
 
      implicit double precision (a-h,o-z), integer (i-n) 
      parameter (lenthin=6552,lenstnd=5733) 
      character cdat*120, ch1stmeter*48, chusername*12
      double precision aatm(5),batm(5),catm(5)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/

c--read one name of a particle data file--------------------------------
      read(*,'(a)',end=479,err=479) cdat
      idat = index(cdat,'DAT')
      ilen = index(cdat,' ') - 1
      read(cdat(idat+3:idat+8),'(i6)') irun
      if ( ilen .le. 0 ) goto 479

c--found valid particle data file---------------------------------------
      call getenv('USER',chusername)
      iuse = index(chusername,' ') - 1
      if ( chusername .eq. 'jl5949' ) then
         ch1stmeter = 'showintact1st.h1stme'
      else
         ch1stmeter =
     +   '/cr/users/'//chusername(1:iuse)//'/showintact1st.h1stme'
      endif

c--read data record to get data and type of simulation------------------
      open(unit=3,file=cdat(1:ilen),status='old',
     +     form='unformatted',access='sequential')
      read(unit=3,err=479,end=479) (pdata(i),i=1,lenstnd)
      close(unit=3)
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         lenrec = lenthin
      endif
      lenblk = lenrec/21
      if ( lenblk .gt. 6543 ) then
         write(*,'(f14.4,'' m'',f13.6,'' gr'',i10.6,i7,4x,a)')
     +      h1stme*1.d-2,h1gram*1.d-2,irun,64,cdat(1:ilen)
         goto 479
      endif
      h1stme = 1.d-2 * pdata(lenblk+7)
      h1gram = thickgr(abs(dble(pdata(lenblk+7))))
      write(*,'(f14.4,'' m'',f13.6,'' gr'',i10.6,i7,4x,a)')
     +    h1stme,h1gram,irun,lenblk,cdat(1:ilen)
      open(unit=7,file=ch1stmeter,status='unknown',
     +     form='formatted',access='sequential')
      write(7,'(f14.4,'' m'',f13.6,'' gr'',i10.6,i7,4x,a)')
     +    h1stme,h1gram,irun,lenblk,cdat(1:ilen)
      close(unit=7)

c--end-of program showintact1st testing files.
  479 continue
      stop
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

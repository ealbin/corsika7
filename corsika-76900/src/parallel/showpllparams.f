c=======================================================================
c
c  s h o w p l l p a r a m s . f
c  -----------------------------
c  create a tabular of all available parallel steering files.
c-----------------------------------------------------------------------
c CompLink:
c     gfortran -O0 -fbounds-check showpllparams.f -o showpllparams
c     ifort -C -O0 -check bounds showpllparams.f -o showpllparams
c RunProg:
c     ./showpllparams.sh
c-----------------------------------------------------------------------
c HILOW       171.71
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program showpllparams

      implicit double precision (a-h,o-z), integer (i-n)
 
      character chzeile*240, chsteer*240, cpllfile*15, chdatum*12
      character cfmtflt*7, cfmtint*4, chthind*6

      dimension ecut(4), ecmnt(4), jcexp(4)

      data cfmtflt/'(f 6.2)'/,cfmtint/'(i1)'/

c - - - - - - print title line of tabular:
      write(*,'(/,''primary  lg(E)  theta   phi'',4x,''runnr'',5x,
     +      ''date'',7x,''ecutmax'',2x,''files'',3x,''Ratio'',2x,
     +      ''obslev'',3x,''Xmagn'',3x,''Zmagn'',2x,''ecutha'',2x,
     +      ''ecutmu'',2x,''ecutel'',2x,''ecutga'',2x,''thinning'',/)')

c - - - - - - read names of parallel steering files:
      jpll = 0
  101 continue
      read(*,'(a)',end=199,err=199) chsteer 
      jlen = 240 
  102 continue
      jlen = jlen - 1
      if ( chsteer(jlen:jlen) .eq. ' ' ) goto 102
      jpar = index(chsteer,'parallel')
      cpllfile = chsteer(jpar:jpar+14)
      chdatum = 'Feb 29  2000'
      if ( jpar .gt. 29 ) chdatum = chsteer(jpar-13:jpar-2)
      open(unit=1,file=cpllfile,form='formatted',
     +     access='sequential',status='old')
      jpll = jpll + 1
      if ( mod(jpll,60) .eq. 0 ) 
     +write(*,'(/,''primary  lg(E)  theta   phi'',4x,''runnr'',5x,
     +      ''date'',7x,''ecutmax'',2x,''files'',3x,''Ratio'',2x,
     +      ''obslev'',3x,''Xmagn'',3x,''Zmagn'',2x,''ecutha'',2x,
     +      ''ecutmu'',2x,''ecutel'',2x,''ecutga'',2x,''thinning'',/)')
      mthi = 0
      lthi = 0
      chthind = '      '

c - - - - - - test contents of parallel steering files:
  103 continue
      read(1,'(a)',end=195,err=195) chzeile
      if ( chzeile(1:4) .eq. 'EXIT' ) goto 195

      ! - - - parallel parameters - - - - - - - - - - - - - - - - - - -
      if ( chzeile(1:8) .eq. 'PARALLEL' ) then
            jbl1 = index( chzeile, ' ')
  110       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 110
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) ectmin
  111       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 111
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 -1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ectmax
      endif

      ! - - - primary particle code - - - - - - - - - - - - - - - - - - 
      if ( chzeile(1:6) .eq. 'PRMPAR' ) then
            jbl1 = index( chzeile(1:240), ' ')
  112       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 112
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtint) mcode
      endif  

      ! - - - primary particle energy - - - - - - - - - - - - - - - - - 
      if ( chzeile(1:6) .eq. 'ERANGE' ) then
            jbl1 = index( chzeile(1:240), ' ')
  113       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 113
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) energy
            energy = 9.+log10(energy)
      endif

      ! - - - theta angle in degrees - - - - - - - - - - - - - - - - - -
      if ( chzeile(1:6) .eq. 'THETAP' ) then
            jbl1 = index( chzeile, ' ')
  114       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 114
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) theta
      endif

      ! - - - azimuth angle in degrees - - - - - - - - - - - - - - - - -
      if ( chzeile(1:6) .eq. 'PHIP  ' ) then
            phia = 359.99d0
            jbl1 = index( chzeile, ' ')
  115       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 115
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) phia
  116       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 116
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) phib
      endif

      ! - - - observation level - - - - - - - - - - - - - - - - - - - -
      if ( chzeile(1:6) .eq. 'OBSLEV' ) then
            jbl1 = index( chzeile, ' ')
  119       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 119
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) obslev
            obslev = obslev / 100.
      endif

      ! - - - magnetic field components - - - - - - - - - - - - - - - -
      if ( chzeile(1:6) .eq. 'MAGNET' ) then
            jbl1 = index( chzeile, ' ')
  120       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 120
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) bxmag
  121       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 121
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) bzmag
      endif

      ! - - - energy cuts of four particle groups - - - - - - - - - - -
      if ( chzeile(1:6) .eq. 'ECUTS ' ) then
            jbl1 = index( chzeile, ' ')
  125       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 125
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) ecut(1)
  126       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 126
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ecut(2)
  127       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 127
            jbl4 = index( chzeile(jbl3:240), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) ecut(3)
  128       continue
            jbl4 = jbl4 + 1
            if ( chzeile(jbl4:jbl4) .eq. ' ' ) goto 128
            jbl5 = index( chzeile(jbl4:240), ' ') + jbl4 - 1
            write(cfmtflt(3:4),'(i2)') jbl5-jbl4
            read(chzeile(jbl4:jbl5-1),cfmtflt) ecut(4)
            do  j=1,4
               if ( ecut(j) .lt. 1. ) then
                  jcexp(j) = int(log10(ecut(j))*1.000003-0.999)
               else
                  jcexp(j) = int(log10(ecut(j))*1.000003)    
               endif  
               ecmnt(j) = ecut(j)/10.d0**jcexp(j)
            enddo  
         endif 

      ! - - - thinning specification - - - - - - - - - - - - - - - - - -
      if ( chzeile(1:7) .eq. 'MTHINH ' ) then
         chthind = chzeile(1:6)
      endif
      if ( chzeile(1:6) .eq. 'THIN ' ) then
         chthind = chzeile(1:5)//' ' 
            jbl1 = index( chzeile, ' ')
  140       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 140
            jbl2 = index( chzeile(jbl1:240), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) thinlev
            thinlev = log10(thinlev+4.67735141e-34)
  141       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 141
            jbl3 = index( chzeile(jbl2:240), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) thiwmax
  142       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 142
            jbl4 = index( chzeile(jbl3:240), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) thirad
            thirad = log10(thirad+4.67735141e-34) - 2. ! meters
            mthi = mthi + 1 
            lthi = 1 
      endif       

      goto 103

c - - - - - - end-of current parallel steering file reached.
  195 continue    
      close(unit=1)

c - - - - - - print tabular of current parallel steering file:
      deratio = 10.d0**(energy-9.)/ectmax
      ipfiles = int(10.d0**(energy-9.)*2.44/ectmax)
      if ( deratio .lt. 100000. ) then
         write(*,'(i6,f8.2,f6.1,f7.1,3x,a6,2x,a12,1p,e9.1,0p,i7,
     +      2f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2),2x,a6)')
     +      mcode, energy, theta, phia, cpllfile(10:15),
     +      chdatum, ectmax, ipfiles, deratio, obslev,
     +      bxmag,bzmag,(ecmnt(j),jcexp(j),j=1,4),chthind
      else
         if ( deratio .gt. 1000000. ) deratio = 999999.d0
         ipfiles = 1
         write(*,'(i6,f8.2,f6.1,f7.1,3x,a6,2x,a12,1p,e9.1,0p,i7,
     +      f8.0,f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2),2x,a6)')
     +      mcode, energy, theta, phia, cpllfile(10:15),
     +      chdatum, ectmax, ipfiles, deratio, obslev,
     +      bxmag,bzmag,(ecmnt(j),jcexp(j),j=1,4),chthind
      endif
      goto 101 

c - - - - - - print closing title line:
  199 continue
      write(*,'(/,''primary  lg(E)  theta   phi'',4x,''runnr'',5x,
     +      ''date'',7x,''ecutmax'',2x,''files'',3x,''Ratio'',2x,
     +      ''obslev'',3x,''Xmagn'',3x,''Zmagn'',2x,''ecutha'',2x,
     +      ''ecutmu'',2x,''ecutel'',2x,''ecutga'',2x,''thinning'',/)')
      stop
      end

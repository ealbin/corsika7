c=======================================================================
c
c  s u m m e . f
c  -------------
c     print sum of bytes of files of given names
c     gfortran -fbounds-check summe.f -o summe
c     ifort -C -check bounds summe.f -o summe
c
c-----------------------------------------------------------------------
      program summe
      implicit double precision (a-h,o-z), integer (i-n) 
      factor = 1.024d-3
      dsumm = 0
    1 continue
      read(*,*,end=9) dkbyt
      dsumm = dsumm + dkbyt
      goto 1
    9 continue
      dsumm = dsumm + 1.
      write(*,'(f35.6,'' GBy'')') factor*dsumm*1.d-3
      stop
      end

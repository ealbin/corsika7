c=======================================================================
c
c  c o r s p l i t e v t s . f
c  ===========================
c     split corsika particle data files - containing several events -
c     into corresponding single event files.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     CompLink: f77 -fbounds-check corsplitevts.f -o corsplitevts
c               gfortran -fbounds-check corsplitevts.f -o corsplitevts
c     RunProgr: ./corsplitevts 
c               `enter file name`
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     ndim=6552 for thinning, ndim=5733 for standard corsika simulation.
c     recl=26216 or recl=22940. 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c           runh=211285.281   evth=217433.078
c           long=52815.2969   evte=3397.39185   rune=3301.33252
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     input-files:
c     unit=14: /lxdata/d3lx67/joe/corsbsp/DAT106239
c     output-files:
c     unit=*: corsplitevts.out (or display)
c     unit=24: DAT106239_1, ......, DAT106239_3, ....... 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     J. Oehlschlaeger, KIT-CN, 01 Sep 2011
c=======================================================================

      program corsplitevts
      parameter (nmax=6552,nshx=99)
      character chdatin*120,chdatout(nshx)*120,chnew*4,cend*4
      dimension ddata(nmax*2),pdata(nmax),drunh(312)
      equivalence (cend,pend)
      data cend/'RUNE'/, jprint/0/

c = = = = 1 = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c - - - - - - read name of file to be splitted and number of showers:
      chdatin( 1:40) = '                                        '
      chdatin(41:80) = '                                        '
      chdatin(81:120)= '                                        '
      read(*,'(a)',end=499,err=499) chdatin

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - check and keep length of file name:
      is = 120 + 1
  406 continue
      is = is - 1
      if ( chdatin(is:is) .eq. ' ' ) goto 406
      lchd = is

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - create and keep names of the new single shower files:
      write(*,'(31('' -''))')
      do  408  is=1,nshx
         chnew = '    '
         chdatout(is)( 1:40) ='                                        '
         chdatout(is)(41:80) ='                                        '
         chdatout(is)(81:120)='                                        '
         lchn = 2 + log10(1.*nshx)
         lchn = lchn + 1  
         write(chnew(1:lchn),'(''_'',i2.2)') is 
         chdatout(is)(1:lchd+lchn) = chdatin(1:lchd)//chnew(1:lchn)
         if ( is .le. 1 ) write(*,*) chdatout(is) 
  408 continue
      read(chdatin(lchd-5:lchd),'(i6)') irun
      prun = irun 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - test first record on thinning or standard simulation:
      open(unit=14,file=chdatin,status='old',form='unformatted')
      read(14,err=494,end=494) (ddata(i),i=1,5733)
      close(14)
      lthi = 1 
      if (217433.0.lt.ddata(274) .and. ddata(274).lt.217433.2) lthi = 0 
      ndim = 5733 + lthi * 819  
      nblk = ndim / 21
      npar = nblk / 39 
      write(*,*) 'corsika particle data file splitting operates '
      if ( lthi .eq. 1 ) then
         write(*,*) '               on thinning simulation, nblk =',nblk
      else
         write(*,*) '               on standard simulation, nblk =',nblk
      endif 

c = = = = 2 = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c - - - - - - implicit loop for the first shower:
      open(unit=14,file=chdatin,status='old',form='unformatted')
      irec = 0
      isho = 0
      iend = 0
      ibshift = 0
      idatshi = ibshift * nblk
      mrecout = 0
      mbstart = 1
      lsho = 1

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - reading next record of the first shower:
  410 continue
      read(14,err=499,end=499) (ddata(i),i=1,ndim)
      irec = irec + 1
      if ( jprint .ge. 1 .and.
     +   mod(irec,1000).eq.0 ) write(*,*) '          irec =',irec
      jblkout = 0
      if ( irec .eq. 1 ) then
         write(*,'(31('' -''))')
         write(*,*) '  isho=',lsho
         do  ii=1,nblk*2,nblk
            write(*,*) '  sblk=',1+ii/nblk
            if ( jprint .ge. 2 ) 
     +         write(*,'(1p,8e14.6)') (ddata(i),i=ii,ii+npar-1) 
         enddo
         do  i=1,ndim 
            pdata(i) = 0.
         enddo
      endif

c - - - - - - - loop over the 21 subblocks of each record:
      do  412  ib=1,ndim,nblk 

c - - - - - - - - subbblock RUNH found as first subblock at all:
         if ( ddata(ib).ge.211285.2 .and. ddata(ib).lt.211285.4 ) then
            do  i=1,nblk
               drunh(i) = ddata(i)
            enddo
            do  i=1,nblk
               pdata(i) = drunh(i)
            enddo
            jblkout = 1

c - - - - - - - - subbblock RUNE found:
         elseif ( ddata(ib).ge.3301.3 .and. ddata(ib).lt.3301.5 ) then
            write(*,'(2x,''rune should not occur in the first'',
     +         1x,''record of thinned showers.'')')
            write(*,*) ' file contains only one shower!'
            goto 499 

c - - - - - - - - subbblock EVTH found:
         elseif (ddata(ib).ge.217433.0 .and. ddata(ib).lt.217433.2) then
            isho = isho + 1
            jblkout = jblkout + 1 
            ! - - - - - copy evth subblock to output array:
            do  i=1+nblk,nblk*2
               pdata(i) = ddata(i)
            enddo
            open(unit=24,file=chdatout(lsho),
     +           status='unknown',form='unformatted')

c - - - - - - - - subblock EVTE found:
         elseif ( ddata(ib).ge.3397.3 .and. ddata(ib).lt.3397.5 ) then
            jblkout = jblkout + 1 
            iend = iend + 1 
            do  i=1,ndim*2,nblk ! optional print out of full ddata:
               if ( i .lt. 0 ) write(*,'(i13,1p,3e14.6)')
     +            1+i/nblk,ddata(i),ddata(i+1),ddata(i+2)
            enddo
            if ( jprint .ge. 2 ) write(*,'(i13,1p,3e14.6)')
     +         1+ib/nblk,ddata(ib),ddata(ib+1),ddata(ib+2)
            ! - - - - - copy evte subblock to output array:
            do  i=ib,ib+nblk-1
               pdata(i) = ddata(i)
            enddo
            ipd = i-idatshi 
c - - - - - - - - - - new rune as subblock nr. 21:
            if ( jblkout .lt. 21 ) then
               jblkout = jblkout + 1
               mrecout = mrecout + 1
               do  n=ipd,ndim
                  pdata(n) = 0.
               enddo 
               pdata(ipd) = pend
               pdata(ipd+1) = prun
               pdata(ipd+2) = -lsho
               pdata(ipd+3) = 1.
               pdata(ipd+4) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
               goto 414
c - - - - - - - - - - new rune as subblock nr. 1:
            elseif ( jblkout .eq. 21 ) then
               mrecout = mrecout + 1
               write(24) (pdata(i),i=1,ndim)
               do  n=1,ndim 
                  pdata(n) = 0.
               enddo
               jblkout = 1
               mrecout = mrecout + 1
               pdata(1) = pend
               pdata(2) = prun
               pdata(3) = -lsho
               pdata(4) = 1.
               pdata(5) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
               goto 414
            else
               write(*,*) ' jblkout = ',jblkout,' > 21 is invalid!'
               goto 499
            endif

c - - - - - - - - subblock LONG found:
         elseif ( ddata(ib).ge.52815.2 .and. ddata(ib).lt.52815.4 ) then
            ! subblocks are ignored.

         else
c - - - - - - - - copy particle data subblock to output array:
            jblkout = jblkout + 1 
            do  i=ib,ib+nblk-1
               pdata(i) = ddata(i)
            enddo
            if ( jblkout .eq. 21 ) mrecout = mrecout + 1
         endif

c - - - - - - - end of loop over 21 subblocks of a corsika record.
  412 continue
      write(24) (pdata(i),i=1,ndim)
      goto 410

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - keep numbers of subblocks after the first shower:
  414 continue
      write(*,'(3x,''isho='',i2,2x,''finished at irec'','//
     +   'i7,2x,''and sblk'',i3,2x,''as rune'')') lsho,irec,jblkout
      write(*,'(42x,''(output records'',i7,'')'')') mrecout
      if ( jblkout .gt. 1 ) then
         do  ii=1+nblk*(jblkout-2),nblk*jblkout,nblk
            if ( jprint .ge. 2 )
     +         write(*,'(6x,''sblk='',i2,1p,5e14.6)')
     +            1+ii/nblk,(pdata(i),i=ii,ii+4) 
         enddo
      elseif ( jblkout .eq. 1 ) then
         do  ii=1,nblk,nblk
            if ( jprint .ge. 2 )
     +         write(*,'(6x,''sblk='',i2,1p,5e14.6)')
     +            1+ii/nblk,(pdata(i),i=ii,ii+4)
         enddo
         ! fill second part of ddata and dont overwrite:
         read(14,err=499,end=499) (ddata(i),i=ndim+1,ndim+ndim)
         irec = irec + 1
         jblkout = jblkout + 21
      endif
      if ( jprint .ge. 2 ) 
     +   write(*,'(10x,3(i13,''.''))') ii-nblk,ii-nblk+1,ii-nblk+2
      write(*,'(31('' -''))')
      do  i=1,ndim*2,nblk ! optional print out of full ddata:
         if ( i .lt. 0 ) write(*,'(i13,1p,3e14.6)')
     +      1+i/nblk,ddata(i),ddata(i+1),ddata(i+2)
      enddo
      if ( lsho .eq. nshx ) goto 499

c = = = = 3 = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c - - - - - - prepare new beginning of next (>=2.) shower:
  448 continue
      lsho = lsho + 1
      if ( lsho .gt. nshx ) goto 499
      write(*,*) chdatout(lsho)
      if ( jprint .ge. 2 ) 
     +   write(*,*) '     mbstart_old =',mbstart
      ibshift = jblkout - 2 
      mbstart = mbstart + ibshift
      if ( jprint .ge. 2 ) 
     +   write(*,*) '     mbstart_new =',mbstart 
      if ( mbstart .gt. 21 ) then
         mbstart = mbstart - 21
         if ( jprint .ge. 2 ) 
     +      write(*,*) '     mbstart_new =',mbstart,' reduced!' 
         do  i=1,ndim
            ddata(i) = ddata(i+ndim)
         enddo  
         read(14,err=499,end=499) (ddata(i),i=ndim+1,ndim+ndim)
         irec = irec + 1
      endif 
      ibshift = mbstart - 1
      idatshi = ibshift * nblk
      if ( ibshift .le. 0 ) ibshift = ibshift + 21
      write(*,*) '  isho=',lsho,'       ibl_shift=',ibshift
      ! - - - - - copy runh subblock to new position for next shower:
      mrecout = 0
      do  i=1,nblk
         ddata(i+idatshi) = drunh(i)
      enddo
      ! - - - - - second part of ddata alrerady read:
      if ( mbstart .eq. 21 ) goto 450

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - read next record as second part of ddata:
  449 continue
      read(14,err=499,end=499) (ddata(i),i=ndim+1,ndim+ndim)
      irec = irec + 1
  450 continue
      if ( jprint .ge. 2 .and. 
     +   mod(irec,1000).eq.0 ) write(*,*) '          irec =',irec
      jblkout = 0

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - loop over the subblocks of each shifted record:
      do  452  ib=idatshi+1,idatshi+ndim,nblk 

c - - - - - - - - subbblock RUNH found as subblock after ibshift:
         if ( ddata(ib).ge.211285.2 .and. ddata(ib).lt.211285.4 ) then
            do  ii=1,nblk*2,nblk
               write(*,*) '  sblk=',1+ii/nblk,
     +            '      sb_current=',mbstart+ii/nblk
               if ( jprint .ge. 2 )
     +         write(*,'(1p,8e14.6)') (ddata(i+idatshi),i=ii,ii+npar-1)
            enddo
            do  i=1,ndim 
               pdata(i) = 0.
            enddo
            do  i=1,nblk
               pdata(i) = drunh(i)
            enddo
            jblkout = 1

c - - - - - - - - subbblock RUNE found:
         elseif ( ddata(ib).ge.3301.2 .and. ddata(ib).lt.3301.4 ) then
            write(*,'(2x,''rune should not occur here.'')')
            goto 499 

c - - - - - - - - subbblock EVTH found:
         elseif (ddata(ib).ge.217433.0 .and. ddata(ib).lt.217433.2) then
            ! isho = isho + 1
            ! - - - - - copy evth subblock to output array:
            do  i=ib,ib+nblk-1
               pdata(i-idatshi) = ddata(i)
            enddo
            jblkout = jblkout + 1 
            open(unit=24,file=chdatout(lsho),
     +           status='unknown',form='unformatted')

c - - - - - - - - subblock EVTE found:
         elseif ( ddata(ib).ge.3397.3 .and. ddata(ib).lt.3397.5 ) then
            jblkout = jblkout + 1 
            iend = iend + 1
            do  i=1,ndim*2,nblk ! optional print out of full ddata:
               if ( i .lt. 0 .and. jprint .ge. 3 )
     +            write(*,'(i13,1p,3e14.6)')
     +            1+i/nblk,ddata(i),ddata(i+1),ddata(i+2)
            enddo
            if ( jprint .ge. 2 ) write(*,'(i13,1p,3e14.6)')
     +         1+ib/nblk,ddata(ib),ddata(ib+1),ddata(ib+2)
            jend = 1+ib/nblk
            ! - - - - - copy evte subblock to output array:
            do  i=ib,ib+nblk-1
               pdata(i-idatshi) = ddata(i)
            enddo
            ipd = i-idatshi 
c - - - - - - - - - - new rune as subblock nr. 21:
            if ( jblkout .lt. 21 ) then
               jblkout = jblkout + 1
               mrecout = mrecout + 1
               do  n=ipd,ndim
                  pdata(n) = 0.
               enddo 
               pdata(ipd) = pend
               pdata(ipd+1) = prun
               pdata(ipd+2) = -lsho
               pdata(ipd+3) = 1.
               pdata(ipd+4) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
               goto 454
c - - - - - - - - - - new rune as subblock nr. 1:
            elseif ( jblkout .eq. 21 ) then
               mrecout = mrecout + 1
               write(24) (pdata(i),i=1,ndim) 
               jblkout = 1
               do  n=1,ndim 
                  pdata(n) = 0.
               enddo
               mrecout = mrecout + 1
               ipd = 1
               pdata(ipd) = pend
               pdata(ipd+1) = prun
               pdata(ipd+2) = -lsho
               pdata(ipd+3) = 1.
               pdata(ipd+4) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
               goto 454
            else
               write(*,*) ' jblkout = ',jblkout,' > 21 is invalid!'
               goto 499
            endif

c - - - - - - - - subblock LONG found:
         elseif ( ddata(ib).ge.52815.2 .and. ddata(ib).lt.52815.4 ) then
            ! subblocks are ignored.

         else
c - - - - - - - - copy particle data subblock to output array:
            jblkout = jblkout + 1 
            do  i=ib,ib+nblk-1
               pdata(i-idatshi) = ddata(i)
            enddo
            if ( jblkout .eq. 21 ) then
               mrecout = mrecout + 1
               do  i=ib+nblk,ndim*2,nblk
               if ( ddata(i).ge.3301.2.and.ddata(i).lt.3301.4 ) goto 453
               enddo 
            endif
         endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - end of loop over all subblocks of shifted corsika record.
  452 continue
      write(24) (pdata(i),i=1,ndim)
      do  i=1,ndim
         ddata(i) = ddata(i+ndim)
      enddo 
      goto 449

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - RUNE subblock found, all possible showers copied:
  453 continue
      write(24) (pdata(i),i=1,ndim)
      do  i=1,ndim
         pdata(i) = 0.
      enddo
      idatshi = idatshi + ndim
      jblkout = 0
      do  n=ib+nblk,ndim*2,nblk
         if ( ddata(n).ge.3301.2 .and. ddata(n).lt.3301.4 ) then
            jblkout = jblkout + 1
            ipd = n-idatshi
            if ( jblkout .lt. 21 ) then
               mrecout = mrecout + 1
               pdata(ipd) = pend
               pdata(ipd+1) = prun
               pdata(ipd+2) = -lsho
               pdata(ipd+3) = 1.
               pdata(ipd+4) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
            elseif ( jblkout .eq. 21 ) then
               if ( jprint .ge. 2 )
     +            write(*,*) '   last rune as subblock nr. 1.'
               mrecout = mrecout + 1
               write(24) (pdata(i),i=1,ndim)
               jblkout = 1
               do  i=1,ndim
                  pdata(i) = 0.
               enddo
               mrecout = mrecout + 1
               ipd = 1
               pdata(ipd) = pend
               pdata(ipd+1) = prun
               pdata(ipd+2) = -lsho
               pdata(ipd+3) = 1.
               pdata(ipd+4) = -lsho ! test to check current shower.
               write(24) (pdata(i),i=1,ndim)
               close(24)
            endif 
         elseif ( ddata(n).ge.3397.3 .and. ddata(n).lt.3397.5 ) then
            jblkout = jblkout + 1
            iend = iend + 1
            if ( jprint .ge. 2 ) write(*,'(i13,1p,3e14.6)')
     +         1+n/nblk,ddata(n),ddata(n+1),ddata(n+2)
            do  i=n,n+nblk-1
               pdata(i-idatshi) = ddata(i)
            enddo
         elseif ( ddata(n) .ge. 1000. ) then 
            jblkout = jblkout + 1
            do  i=n,n+nblk-1
               pdata(i-idatshi) = ddata(i)
            enddo
         endif
      enddo 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - keep numbers of subblocks after >= second shower:
  454 continue
      write(*,'(3x,''isho='',i2,2x,''finished at irec'','//
     +   'i7,2x,''and sblk'',i3,2x,''as rune'')') lsho,irec,jblkout
      write(*,'(42x,''(output records'',i7,'')'')') mrecout
      if ( jblkout .gt. 1 .and. jblkout .le. 21 ) then
         do  ii=1+nblk*(jblkout-2),nblk*jblkout,nblk
            if ( jprint .ge. 2 )
     +         write(*,'(6x,''sblk='',i2,1p,5e14.6)')
     +            1+ii/nblk,(pdata(i),i=ii,ii+4)
         enddo
      elseif ( jblkout .eq. 1 ) then  
         do  ii=1,nblk,nblk
            if ( jprint .ge. 2 )
     +         write(*,'(6x,''sblk='',i2,1p,5e14.6)')
     +            1+ii/nblk,(pdata(i),i=ii,ii+4)
         enddo
      endif
      if ( jprint .ge. 2 )
     +   write(*,'(10x,3(i13,''.''))') ii-nblk,ii-nblk+1,ii-nblk+2
      write(*,'(31('' -''))')
      if ( jblkout .eq. 1 ) jblkout = jblkout + 21
      if ( lsho .lt. nshx ) goto 448
      goto 499

c = = = = 4 = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c - - - - - - end of program corsplitevts.
  494 continue
      write(*,'(31('' -''))')
      write(*,*) ' e r r o r: check kind of simulation, '
      write(*,*) '            thinning reclen = 6552, ',
     + 'standard reclen = 5733.  ithi =',ithi
  499 continue
      if ( iend .lt. lsho .and. iend .ge. 1 ) then
         write(*,*) chdatin(1:lchd),' contains only',iend,' showers!'
         write(*,'(31('' -''))')
      endif
      write(*,*) ' end of program corsplitevts.'
      stop
      end

C=======================================================================
C
C  c o r s i k a r e a d _ t h i n . f  with "thinning"
C           ====================================================
C                 READ  AND  PRINT  CORSIKA  SHOWER  DATA
C           ====================================================
C     Output format for particle output (blocklength = 26208+8 fixed)
C     each block consists of 21 subblocks of 312 words.
C----------------------------------------------------------------------
C     compilation:
C        gfortran -fbounds-check -frecord-marker=4 corsikaread_thin.f -o corsikaread_thin
C        f77 -fbounds-check -m32 corsikaread_thin.f -o corsikaread_thin
C        ifort -C corsikaread_thin.f -o corsikaread_thin
C----------------------------------------------------------------------
C     How to use this program:
C     1) Generate a file 'input' containing the path and name of the 
C        DATnnnnnn file to be analyzed by this program.
C        The name should not contain leading blanks but filled up
C        by trailing blanks to get a total length of >70 characters.
C     2) Execute this program with the file 'input' as standard input:
C              ./corsikaread_thin <input >output
C     3) The file 'output' will contain a short overview of the 
C        content of the DATnnnnnn file to be analyzed.
C     4) The file fort.8 will contain a detailed print out of the
C        content of DATnnnnnn.
C----------------------------------------------------------------------
C     J.Oehlschlaeger, D. Heck, 01 Sep 2011
C=======================================================================
      PROGRAM CORSIKAREAD_THIN
      CHARACTER CHV(6552)*4,CIDENT*4,CDAT*70,CBLK*70
      DIMENSION PDATA(6552)
      EQUIVALENCE (CHV(1),PDATA(1))
      COMMON /CHARS/CHV,CDAT,CBLK,CIDENT
      CBLK='                                                  '
      CDAT=CBLK
      IREC=0
 
C--READ FILE NAME-------------------------------------------------------
      READ(*,428,END=440,ERR=439) CDAT
  428 FORMAT(A)
  429 CONTINUE
      WRITE(*,430) CDAT
  430 FORMAT(1H ,'READ DATA FROM FILE = ',A)
      OPEN(UNIT=3,FILE=CDAT,STATUS='OLD',FORM='UNFORMATTED')
* - - - - - - read data records with 6552 words - - - -
  431 CONTINUE
      IREC = IREC + 1
      READ(UNIT=3,ERR=434,END=433) PDATA
      if ( mod(irec,100) .eq. 0 ) 
     +   WRITE(*,*)'         HAVE READ RECORD NR.',IREC
C-----------loop over subblocks-----------------------------------------
      DO    LIA=1,6552,312
        CIDENT(1:1) = CHV(LIA)(1:1)
        CIDENT(2:2) = CHV(LIA)(2:2)
        CIDENT(3:3) = CHV(LIA)(3:3)
        CIDENT(4:4) = CHV(LIA)(4:4)
        IF (PDATA(LIA).GE.211284.0 .AND.
     +      PDATA(LIA).LE.211286.0) THEN
          CIDENT = 'RUNH'
          WRITE(*,*)'RUNH'
        ENDIF
        IF (PDATA(LIA).GE.217432.0 .AND.
     +      PDATA(LIA).LE.217434.0) THEN
          CIDENT = 'EVTH'
          WRITE(*,*)'EVTH'
        ENDIF
        IF (PDATA(LIA).GE. 52814.0 .AND.
     +      PDATA(LIA).LE. 52816.0) THEN
          CIDENT = 'LONG'
          WRITE(*,*)'LONG'
        ENDIF
        IF (PDATA(LIA).GE.  3396.0 .AND.
     +      PDATA(LIA).LE.  3398.0) THEN
          CIDENT = 'EVTE'
          WRITE(*,*)'EVTE'
        ENDIF
        IF (PDATA(LIA).GE.  3300.0 .AND.
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
            DO    IL=LIA,LIA+311,8
              WRITE(8,'(1P,8E13.5)') (PDATA(II+IL),II=0,7)
            ENDDO
          ELSEIF ( CIDENT .EQ. 'EVTH' ) THEN
C----------------subblock event header----------------------------------
            PDATA(LIA) = 33333333.
            DO    IL=LIA,LIA+311,8
              WRITE(8,'(1P,8E13.5)') (PDATA(II+IL),II=0,7)
            ENDDO
C----------------subblock longitudinal data-----------------------------
          ELSEIF ( CIDENT .EQ. 'LONG' ) THEN
            PDATA(LIA) = 55555555.
            DO    IL=LIA,LIA+311,8
              WRITE(8,'(1P,8E13.5)') (PDATA(II+IL),II=0,7)
            ENDDO
C----------------subblock event end-------------------------------------
          ELSEIF ( CIDENT .EQ. 'EVTE' ) THEN
            PDATA(LIA) = 77777777.
            DO    IL=LIA,LIA+311,8
              WRITE(8,'(1P,8E13.5)') (PDATA(II+IL),II=0,7)
            ENDDO
C----------------subblock run end---------------------------------------
          ELSEIF ( CIDENT .EQ. 'RUNE' ) THEN
            PDATA(LIA) = 99999999.
            DO    IL=LIA,LIA+311,8
              WRITE(8,'(1P,8E13.5)') (PDATA(II+IL),II=0,7)
            ENDDO
            GOTO 929
          ENDIF
        ELSE
C-----------subblock with particle data---------------------------------
          DO    IL=LIA,LIA+311,8
            WRITE(8,'(1P,E16.8,7E13.5)') (PDATA(II+IL),II=0,7)
          ENDDO
        ENDIF
      ENDDO

  929 CONTINUE
      GOTO 431
 
C--END OF TEST----------------------------------------------------------
  433 CONTINUE
      WRITE(*,*)'         LAST RECORD ',irec-1
      CLOSE(UNIT=3)
      STOP
  434 CONTINUE
      WRITE(*,*)'         READ ERROR ON UNIT 3'
*     CLOSE(UNIT=3)
      GOTO 431
  439 CONTINUE
      WRITE(*,*)'         READ ERROR ON STANDARD INPUT'
      GOTO 429
  440 CONTINUE
      WRITE(*,*)'         READ END ON STANDARD INPUT'
      STOP
      END

C=======================================================================
C
C                           V E N U S  4.12
C                 (VERY ENERGETIC NUCLEAR SCATTERING)
C                      SUBROUTINE-TURBO-VERSION            MAY 04, 1994
C
C-----------------------------------------------------------------------
C
C   AUTHOR:
C   -------
C
C   KLAUS WERNER
C
C   UNIV. HEIDELBERG, INST. F. THEORETISCHE PHYSIK
C   PHILOSOPHENWEG 19, 6900 HEIDELBERG, GERMANY
C   WERNER@DHDMPI5.BITNET, WERNER@MINNIE.MPI-HD.MPG.DE, 28877::WERNER
C
C   NOW AT: 
C   UNIVERSITE DE NANTES, SUBATECH - ECOLE DES MINES
C   4, RUE ALFRED KASTLER, F44070 NANTES CEDEX 03, FRANCE
C   E-MAIL: WERNER@SUBATECH.IN2P3.FR
C   FAX:    (33)  51 85 84 79
C   TEL:    (33)  51 85 84 25
C
C
C   SUBROUTINE-TURBO-VERSION 4.12/5:
C   --------------------------------
C   DIETER HECK
C   FORSCHUNGSZENTRUM KARLSRUHE, INST. F. KERNPHYSIK 3
C   POSTFACH 3640, D-76021 KARLSRUHE, GERMANY
C   E-MAIL:   HECK@IK3.FZK.DE
C   FAX:     07247-82-4075
C   TEL:     07247-82-3777
C
C   LAST CHANGES:  OCT  05, 1995
C   rename of jdcay to jdecayv      oct. 18, 1999   by D.H.
C   include SAVE statements         Sept.26, 2000   by D.H.
C-----------------------------------------------------------------------
C
C   VENUS:
C   ------
C
C  VENUS IS A MONTE CARLO PROCEDURE TO SIMULATE HADRONIC INTERACTIONS AT
C  ULTRARELATIVISTIC ENERGIES (HADRON-HADRON, HADRON-NUCLEUS, NUCLEUS-
C  NUCLEUS SCATTERING), AND ALSO INTERACTIONS INVOLVING LEPTONS (E+E-
C  ANNIHILATION, LEPTON-NUCLEON, LEPTON-NUCLEUS SCATTERING).
C  VENUS IS BASED ON GRIBOV-REGGE THEORY (OF MULTIPLE POMERON EXCHANGE)
C  AND CLASSICAL RELATIVISTIC STRING DYNAMICS. A DETAILED DESCRIPTION CA
C  BE FOUND IN THE UNIV. HEIDELBERG PREPRINT HD-TVP-93-1 (270 PAGES),
C  WHICH IS PUBLISHED IN PHYSICS REPORTS 232 (1993) 87-299.
C
C   DISTRIBUTION:
C   -------------
C
C  THIS CODE SHOULD NOT BE DISTRIBUTED WITHOUT NOTIFYING THE AUTHOR, IN-
C  DICATING WHAT THE CODE IS GOING TO BE USED FOR. DEPENDING ON THE TYPE
C  OF REACTION AND THE KIND OF OBSERVABLES ANALYSED,  THE SYSTEMATIC UN-
C  CERTAINTIES OF THE VENUS SIMULATIONS VARY STRONGLY,  AND THIS SHOULD
C  BE CLARIFIED BEFORE USING VENUS.
C
C   IMPORTANT FEATURES:
C   -------------------
C
C  COVARIANT TREATMENT OF SECONDARY INTERACTIONS. EACH PRODUCED PARTICLE
C  IS ALLOWED TO REINTERACT WITH OTHER PRODUCED PARTICLES OR WITH
C  SPECTATORS. IMPORTANT FOR HADRON-NUCLEUS, NUCLEUS-NUCLEUS AND
C  LEPTON-NUCLEUS SCATTERING.
C  NO FINAL STATE INTERACTION, IF 'RADIAC' IS SET ZERO.
C
C  PARTICIPATION OF ANTIQUARKS (IN ADDITION TO QUARKS) IN THE
C  COLOUR EXCHANGE MECHANISM TO FORM STRINGS.
C
C  POSSIBILITY OF DIQUARK BREAKUP, LEADING TO MULTI-STRINGS, LIKE
C  A FORWARD QUARK LINKED VIA TWO (!) STRINGS TO TWO BACKWARD
C  QUARKS (DOUBLE-STRING). SUCH STRINGS FRAGMENT DIFFERENTLY THAN
C  QUARK-DIQUARK STRINGS. IN CASE OF THE DOUBLE-STRING, THE FORWARD
C  QUARK WILL FRAGMENT VIA TWO BREAKS INTO A LEADING BARYON.
C
C  SOPHISTICATED FRAGMENTATION PROCEDURE. SINCE SPACE-TIME
C  EVOLUTION IS AN IMPORTANT ISSUE CONCERNING FINAL STATE INTER-
C  ACTIONS, IT IS NOT ENOUGH TO HAVE A FRAGMENTATION MODEL WHICH
C  WORKS, IT SHOULD HAVE THE RIGHT SPACE-TIME DESCRIPTION! THERE-
C  FOR THE FIELD-FEYNMAN MODEL OF EARLIER VERSIONS (<3.00) HAS BEEN
C  ABANDONED AND REPLACED BY A VERY POWERFUL AND APPEALING PROCEDURE
C  SUGGESTED BY ARTRU/MENNESSIER.
C
C  VERY LARGE RESONANCE TABLE, INCLUDING FOR EXAMPLE ALL NUCLEON
C  RESONANCES UP TO 2 GEV. THE DECAY TABLE INCLUDES ALSO KSHORT
C  AND LAMBDA DECAYS. DECAY SUPPRESSION IS GOVERNED BY THE PARAMETERS
C  'NDECAY' AND 'NDECAX'. SETTING 'NDECAY' TO 1 SUPPRESSES ALL DECAYS.
C  FOR FURTHER DETAILS SEARCH FOR 'DECAY SUPPRESSION' IN SUBR. AINITL
C
C  FOR VERSION >= 4.01, WE USE GRIBOVS CUT-POMERON PROBABILITIES TO
C  DETERMINE THE WEIGHTS FOR MULTI-COLOUR-EXCHANGE.
C
C-----------------------------------------------------------------------
C
C   PARTICLE ID (SAME AS ISAJET - SEE ISAJET.DOC, F.E.PAIGE)
C   --------------------------------------------------------
C
C  QUARKS AND LEPTONS ARE NUMBERED IN ORDER OF MASS:
C        UP=1 DN=2 ST=3 CH=4 BT=5 TP=6
C        NUE=11 E-=12 NUM=13 MU-=14 NUT=15 TAU-=16
C  WITH A NEGATIVE SIGN FOR ANTIPARTICLES.
C  ARBITRARY CONVENTIONS ARE:
C        GL=9 GM=10 KS=20 KL=-20 W+=80 Z0=90 .
C  THE CODE FOR A MESON IS AN INTEGER +-JKL, WHERE J.LE.K ARE THE
C  QUARKS AND L IS THE SPIN. THE SIGN IS FOR THE J QUARK.
C  FLAVOR SINGLET MESONS ARE ORDERED BY MASS:
C        PI0=110 ETA=220 ETAP=330 ETAC=440 .
C  SIMILARLY, THE CODE FOR A BARYON IS A COMPOUND INTEGER +-IJKL
C  FORMED FROM THE THREE QUARKS I,J,K AND A SPIN LABEL L=0,1.
C  THE CODE FOR A DIQUARK IS +-IJ00.
C
C   LIST OF IDENT CODES:
C
C     IDENT     LABEL           MASS    CHARGE
C         1     UP        .30000E+00       .67
C        -1     UB        .30000E+00      -.67
C         2     DN        .30000E+00      -.33
C        -2     DB        .30000E+00       .33
C         3     ST        .50000E+00      -.33
C        -3     SB        .50000E+00       .33
C         4     CH        .16000E+01       .67
C        -4     CB        .16000E+01      -.67
C         5     BT        .49000E+01      -.33
C        -5     BB        .49000E+01       .33
C         6     TP        .30000E+02       .67
C        -6     TB        .30000E+02      -.67
C
C         9     GL       0.               0.00
C
C        10     GM       0.               0.00
C
C        11     NUE      0.               0.00
C       -11     ANUE     0.               0.00
C        12     E-        .51100E-03     -1.00
C       -12     E+        .51100E-03      1.00
C        13     NUM      0.               0.00
C       -13     ANUM     0.               0.00
C        14     MU-       .10566E+00     -1.00
C       -14     MU+       .10566E+00      1.00
C        15     NUT      0.               0.00
C       -15     ANUT     0.               0.00
C        16     TAU-      .18070E+01     -1.00
C       -16     TAU+      .18070E+01      1.00
C
C        20     KS        .49767E+00      0.00
C       -20     KL        .49767E+00      0.00
C
C        80     W+        SIN2W=.23       1.00
C       -80     W-        SIN2W=.23      -1.00
C        90     Z0        SIN2W=.23       0.00
C
C       110     PI0       .13496E+00      0.00
C       120     PI+       .13957E+00      1.00
C      -120     PI-       .13957E+00     -1.00
C       220     ETA       .54880E+00      0.00
C       130     K+        .49367E+00      1.00
C      -130     K-        .49367E+00     -1.00
C       230     K0        .49767E+00      0.00
C      -230     AK0       .49767E+00      0.00
C       330     ETAP      .95760E+00      0.00
C       140     AD0       .18633E+01      0.00
C      -140     D0        .18633E+01      0.00
C       240     D-        .18683E+01     -1.00
C      -240     D+        .18683E+01      1.00
C       340     F-        .20300E+01     -1.00
C      -340     F+        .20300E+01      1.00
C       440     ETAC      .29760E+01      0.00
C
C       111     RHO0      .77000E+00      0.00
C       121     RHO+      .77000E+00      1.00
C      -121     RHO-      .77000E+00     -1.00
C       221     OMEG      .78260E+00      0.00
C       131     K*+       .88810E+00      1.00
C      -131     K*-       .88810E+00     -1.00
C       231     K*0       .89220E+00      0.00
C      -231     AK*0      .89220E+00      0.00
C       331     PHI       .10196E+01      0.00
C       141     AD*0      .20060E+01      0.00
C      -141     D*0       .20060E+01      0.00
C       241     D*-       .20086E+01     -1.00
C      -241     D*+       .20086E+01      1.00
C       341     F*-       .21400E+01     -1.00
C      -341     F*+       .21400E+01      1.00
C       441     JPSI      .30970E+01      0.00
C
C      1120     P         .93828E+00      1.00
C     -1120     AP        .93828E+00     -1.00
C      1220     N         .93957E+00      0.00
C     -1220     AN        .93957E+00      0.00
C      1130     S+        .11894E+01      1.00
C     -1130     AS-       .11894E+01     -1.00
C      1230     S0        .11925E+01      0.00
C     -1230     AS0       .11925E+01      0.00
C      2130     L         .11156E+01      0.00
C     -2130     AL        .11156E+01      0.00
C      2230     S-        .11974E+01     -1.00
C     -2230     AS+       .11974E+01      1.00
C      1330     XI0       .13149E+01      0.00
C     -1330     AXI0      .13149E+01      0.00
C      2330     XI-       .13213E+01     -1.00
C     -2330     AXI+      .13213E+01      1.00
C      1140     SC++      .24300E+01      2.00
C     -1140     ASC--     .24300E+01     -2.00
C      1240     SC+       .24300E+01      1.00
C     -1240     ASC-      .24300E+01     -1.00
C      2140     LC+       .22600E+01      1.00
C     -2140     ALC-      .22600E+01     -1.00
C      2240     SC0       .24300E+01      0.00
C     -2240     ASC0      .24300E+01      0.00
C      1340     USC.      .25000E+01      1.00
C     -1340     AUSC.     .25000E+01     -1.00
C      3140     SUC.      .24000E+01      1.00
C     -3140     ASUC.     .24000E+01     -1.00
C      2340     DSC.      .25000E+01      0.00
C     -2340     ADSC.     .25000E+01      0.00
C      3240     SDC.      .24000E+01      0.00
C     -3240     ASDC.     .24000E+01      0.00
C      3340     SSC.      .26000E+01      0.00
C     -3340     ASSC.     .26000E+01      0.00
C      1440     UCC.      .35500E+01      2.00
C     -1440     AUCC.     .35500E+01     -2.00
C      2440     DCC.      .35500E+01      1.00
C     -2440     ADCC.     .35500E+01     -1.00
C      3440     SCC.      .37000E+01      1.00
C     -3440     ASCC.     .37000E+01     -1.00
C
C      1111     DL++      .12320E+01      2.00
C     -1111     ADL--     .12320E+01     -2.00
C      1121     DL+       .12320E+01      1.00
C     -1121     ADL-      .12320E+01     -1.00
C      1221     DL0       .12320E+01      0.00
C     -1221     ADL0      .12320E+01      0.00
C      2221     DL-       .12320E+01     -1.00
C     -2221     ADL+      .12320E+01      1.00
C      1131     S*+       .13823E+01      1.00
C     -1131     AS*-      .13823E+01     -1.00
C      1231     S*0       .13820E+01      0.00
C     -1231     AS*0      .13820E+01      0.00
C      2231     S*-       .13875E+01     -1.00
C     -2231     AS*+      .13875E+01      1.00
C      1331     XI*0      .15318E+01      0.00
C     -1331     AXI*0     .15318E+01      0.00
C      2331     XI*-      .15350E+01     -1.00
C     -2331     AXI*+     .15350E+01      1.00
C      3331     OM-       .16722E+01     -1.00
C     -3331     AOM+      .16722E+01      1.00
C      1141     UUC*      .26300E+01      2.00
C     -1141     AUUC*     .26300E+01     -2.00
C      1241     UDC*      .26300E+01      1.00
C     -1241     AUDC*     .26300E+01     -1.00
C      2241     DDC*      .26300E+01      0.00
C     -2241     ADDC*     .26300E+01      0.00
C      1341     USC*      .27000E+01      1.00
C     -1341     AUSC*     .27000E+01     -1.00
C      2341     DSC*      .27000E+01      0.00
C     -2341     ADSC*     .27000E+01      0.00
C      3341     SSC*      .28000E+01      0.00
C     -3341     ASSC*     .28000E+01      0.00
C      1441     UCC*      .37500E+01      2.00
C     -1441     AUCC*     .37500E+01     -2.00
C      2441     DCC*      .37500E+01      1.00
C     -2441     ADCC*     .37500E+01     -1.00
C      3441     SCC*      .39000E+01      1.00
C     -3441     ASCC*     .39000E+01     -1.00
C      4441     CCC*      .48000E+01      2.00
C     -4441     ACCC*     .48000E+01     -2.00
C-----------------------------------------------------------------------
C
C   LIST OF SUBROUTINES AND FUNCTIONS OF SUBROUTINE-TURBO-VERSION 4.12
C   ------------------------------------------------------------------
C
C     AVENUS  HAHABS  HAPAPA  HASI    HASTFC  HASTFL  HASTFR  HASTFS
C     HASTFW  HASTPR  HDECMP  HRESCL  IDCOMJ  IDCOMK  IDCOMP  IDDECO
C     IDENCO  IDFLAV  IDLABL  IDMASS  IDMIX   IDQUAC  IDRES   IDRESI
C     IDSPIN  IDTAU   IDTRA   IDTRAI  IDTRB   IDTRBI  IDTR4
C     JAMBR1  JAMBR2  JAMFRA  JCENTD  JCENTP  JCENTR  JCLUDE  JDECA
C     JDECAYV JDECIN  JESTPR  JETGEN  JFRADE  JINTA1  JINTA2
C     JINTCC  JINTCE  JINTCH  JINTCL  JINTEL  JINTFP  JINTFS
C     JINTFU  JINTPA  JRESCL  JSPLIT  LEPEXP  LEPSTR  LEPTAR
C     NUCINI  NUCLCO  NUCOGE  NUCOLL  NUCSTR  PVJPSF
C     RACPRO  RANSTC  RANXQ   SBET    SDENSI  SGAM    SGAU
C     SHOPAR  SJCENT  SJCEN4  SJCGAM  SMASS   SMASSI  SMASSP  SMASST
C     SPTF    SPTH    SPTJ    SPTQ    SSE0    SSE1    SSPLIT  SSPLIX
C     STAA    STXD    STXS    STXU    STXUS   STXZNE  STXZPR  SVA0
C     SVA1    UINTEG  UTACOS  UTAMNU  UTAMNX  UTAMNY  UTAMNZ  UTAMST
C     UTAXIS  UTCHM   UTCLEA  UTHIST  UTHSEA  UTINVT  UTKSIX
C     UTKSTR  UTLOBO  UTLOB2  UTLOC   UTMSG   UTMSGF  UTOVEL  UTPAGE
C     UTPART  UTPCM   UTQUAF  UTQZ    UTREMB  UTREPL  UTRESM  UTREST
C     UTROTA  UTROT2  UTSTOP  UTTAIN  UTTAIX  UTTAUS  UTTAUT  UTTUCL
C
C-----------------------------------------------------------------------
C
C   SUBROUTINE-TURBO-VERSION 4.12:
C   ------------------------------
C
C     IT IS TUNED FOR FASTER CALCULATION, ESPECIALLY IF ENERGY AND
C   TARGET ARE VARYING FROM COLLISION TO COLLISION, AS IS THE CASE IN
C   THE DEVELOPMENT OF AN AIR SHOWER CASCADE.
C     IT NEEDS LINKING ROUTINES FOR INPUT AND OUTPUT AND A RANDOM
C   GENERATOR. ALL PARAMETERS MUST BE SET IN THE LINKING ROUTINES.
C   FOR CONNECTION WITH THE AIR SHOWER PROGRAM 'CORSIKA', SUCH ROUTINES
C   ARE AVAILABLE:
C        RANGEN:    RANDOM GENERATOR
C        UTQSEA:    CALC. OF SEA     QUARK STRUCTURE FUNCTION INTEGRAL
C        UTQVAL:    CALC. OF VALENCE QUARK STRUCTURE FUNCTION INTEGRAL
C        VENDAT:    INITIALIZATION OF PARTICLE ID TABLE
C        VENINI:    FOR FIRST INITIALIZATION OF PARAMETERS
C        VENLNK:    FOR INITIALIZATION OF ENERGY DEPENDENT PARAMETERS
C        VSTORE:    TO STORE THE GENERATED SECONDARY PARTICLES
C
C   NOTE: THE COMMONS ARE IN GENERAL DIFFERENT FROM THE UNTUNED
C         VENUS VERSION 4.12
C-----------------------------------------------------------------------
C  COPYRIGHT AND ANY OTHER APPROPRIATE LEGAL PROTECTION OF THESE
C  COMPUTER PROGRAMS AND ASSOCIATED DOCUMENTATION RESERVED IN ALL
C  COUNTRIES OF THE WORLD.
C
C  FZK WELCOMES COMMENTS CONCERNING THE VENUS CODE BUT UNDERTAKES NO
C  OBLIGATION FOR MAINTENANCE OF THE PROGRAMS, NOR RESPONSIBILITY FOR
C  THEIR CORRECTNESS, AND ACCEPTS NO LIABILITY WHATSOEVER RESULTING
C  FROM THE USE OF ITS PROGRAMS.
C=======================================================================

      SUBROUTINE AVENUS

C-----------------------------------------------------------------------
C  GENERATES ONE VENUS EVENT
C-----------------------------------------------------------------------
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      DOUBLE PRECISION SEEDC,SEEDI
      COMMON /CSEED/   SEEDC,SEEDI
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      ISHS0=ISHSUB
      IF ( ISHEVT .NE. 0  .AND.  NREVT+1 .NE. ISHEVT ) THEN
        ISH=0
        ISHSUB=0
      ENDIF
CDH   IF ( NREVT.EQ.0 .AND. (ISH.EQ.13.OR.ISH.EQ.14) ) CALL UTTIMT
CDH   IF ( ISH .EQ. 14 ) CALL UTTIMA('*** AVENUS *** ')
      IF ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'AVENUS (ENTRY)'
      ENDIF
      IPAGE=0
      CALL UTPAGE
      IF ( ISH .GE. 91 )
     *  WRITE(IFCH,113)('-',L=1,79),IPAGE,NREVT+1,SEEDC,('-',L=1,79)
113   FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *  ,5X,'EVT:',I4,5X,'SEEDC:',D25.15,/,1X,79A1,/)
      IF ( ICHOIC. EQ. 4  .OR.  ICHOIC .EQ. 5 ) GOTO 1000
      BIMEVT=-1
      NTEVT0=NTEVT
3     NTEVT=NTEVT0
2     NTEVT=NTEVT+1
      IF     ( ICHOIC .EQ. 1 ) THEN
        CALL JETGEN(IER)
        IF ( IER .EQ. 1 ) GOTO 3
      ELSEIF ( ICHOIC .EQ. 2  .OR.  ICHOIC .EQ. 3 ) THEN
        CALL NUCOLL
      ENDIF
      IF ( ICHOIC .EQ. 1  .OR.  ICHOIC .EQ. 2 ) THEN
        CALL JFRADE(IER)
        IF ( IER. EQ. 1 ) GOTO 3
      ENDIF
      IF ( (ICHOIC.EQ.2 .OR. ICHOIC.EQ.3) .AND. NEVT.EQ.0 ) THEN
        BIMEVT=-1
        GOTO 2
      ENDIF
      IF ( JERR .GT. 0 ) THEN
        CALL UTSTOP('AVENUS: JERR > 0                        ')
      ENDIF
      NREVT=NREVT+1

1000  CONTINUE
CDH   IF ( ISH .EQ. 14 ) CALL UTTIMA('    AVENUS F   ')
      IF ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)'AVENUS (EXIT)'
      ENDIF
      ISH=ISH0
      ISHSUB=ISHS0
      RETURN
      END
C=======================================================================

      SUBROUTINE HAHABS(PROJ,TARG,IAP,IAT,ISKIP,IRETHH)

C-----------------------------------------------------------------------
C  PERFORMS A BASIC (ONE COLOR EXCHANGE) HADRON-HADRON COLLISION
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CKOL/    KOL
      COMMON /CLEAD/   COOAV3,COOAV4,LEAD
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNEW/    KOTRI,NEWCOL,NEWICO
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPZSTR/  ESTRL,PZSTRL,ISEA,ISTRL
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CSTSH/   NSTSH
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      PROJ(NSI,NHA),PSUM(5),PSUMX(5)
     *         ,SIL(NSI),SILP(NSI),SILT(NSI)
     *         ,SIX(NSI,NSIX),SIXP(NSI,NSIX),SIXT(NSI,NSIX)
     *         ,STRAP(NSI,NSIX+1),STRAT(NSI,NSIX+1)
     *         ,STRLP(NSI,NSIX+1),STRLT(NSI,NSIX+1)
     *         ,STR0(NSI,NSIX+1),TARG(NSI,NHA)
      INTEGER   ICVP(2),ICVT(2),JCVP(NFLAV,2),JCVT(NFLAV,2)
      CHARACTER XFLAP*3,XFLAT*3,XFLBP*3,XFLBT*3
      SAVE
C-----------------------------------------------------------------------
      R=RANGEN()
      NTRY=0
      IRETHH=0
      CALL UTREMB(PROJ,TARG,2)
      GOTO 9

9994  IRET=0
      CALL UTREST(PROJ,TARG,2)
 9    CONTINUE
      DO 8 NX=1,NSIX
        DO 8 N=1,NSI
          SIXP(N,NX)=0.
          SIXT(N,NX)=0.
 8    CONTINUE
      CALL HDECMP(PROJ,SIL,SIX)
      CALL UTKSIX(SIX,KMAXP)
      CALL HDECMP(TARG,SIL,SIX)
      CALL UTKSIX(SIX,KMAXT)
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,100) ( (PROJ(I,J),I=1,NSI), J=2,KMAXP+2 )
100     FORMAT ( ' PROJ:',4F13.5,2F8.0,/,20('      ',4F13.5,2F8.0,/) )
        WRITE(IFCH,102) ( (TARG(I,J),I=1,NSI), J=2,KMAXT+2 )
102     FORMAT ( ' TARG:',4F13.5,2F8.0,/,50('      ',4F13.5,2F8.0,/) )
      ENDIF
      IF ( ISKIP .EQ. 1 ) GOTO 9999

9997  CALL HDECMP(PROJ,SILP,SIXP)
      CALL HDECMP(TARG,SILT,SIXT)
 1    ICVP(1)=NINT(PROJ(5,1))
      ICVP(2)=NINT(PROJ(6,1))
      ICVT(1)=NINT(TARG(5,1))
      ICVT(2)=NINT(TARG(6,1))
      CALL HAPAPA(SILP,IFLAP,XFLAP,PTAP,PTAHP,IFLBP,XFLBP,PTBP,PTBHP
     *            ,NQAQP,ICVP)
      CALL HAPAPA(SILT,IFLAT,XFLAT,PTAT,PTAHT,IFLBT,XFLBT,PTBT,PTBHT
     *            ,NQAQT,ICVT)
      IF ( NQAQP*NQAQT .EQ. -1 ) GOTO 1

9998  NTRY=NTRY+1
      IF ( NTRY .GT. NTRYMX+1 ) THEN
        ISKIP=1
        GOTO 9999
      ENDIF

      CALL HDECMP(PROJ,SILP,SIXP)
      CALL HDECMP(TARG,SILT,SIXT)
      CALL HASTFS(SILP,SIXP,IFLAP,XFLAP,PTAP,PTAHP
     *            ,SILT,SIXT,IFLBT,XFLBT,PTBT,PTBHT,STRAP,IRET,1)
      IF ( IRET .EQ. 1 ) THEN
        ISKIP=1
        GOTO 9999
      ENDIF
      IF ( IRET .EQ. 2 ) GOTO 9997
      IF ( IRET .EQ. 3 ) THEN
        IF ( XFLBP(1:2) .EQ. 'VA' ) XFLBP(1:2)='SE'
        IF ( XFLAT(1:2) .EQ. 'VA' ) XFLAT(1:2)='SE'
        IFLAT=IFLBT
        IFLBP=IFLAP
        GOTO 9998
      ENDIF
      IF ( IRET .EQ. 5 ) THEN
        IF ( XFLAP(1:2) .EQ. 'SE'  .AND.  XFLBP(1:2) .EQ. 'SE' ) THEN
          IFLBP=-ABS(IFLBP)
          IFLAP=IFLBP
        ENDIF
        IF ( XFLAT(1:2) .EQ. 'SE'  .AND.  XFLBT(1:2) .EQ. 'SE' ) THEN
          IFLAT=-ABS(IFLAT)
          IFLBT=IFLAT
        ENDIF
        GOTO 9998
      ENDIF
      CALL HASTFS(SILT,SIXT,IFLAT,XFLAT,PTAT,PTAHT
     *,SILP,SIXP,IFLBP,XFLBP,PTBP,PTBHP,STRAT,IRET,2)
      IF ( IRET .EQ. 1 ) THEN
        ISKIP=1
        GOTO 9999
      ENDIF
      IF ( IRET .EQ. 2 ) GOTO 9997
      IF ( IRET .EQ. 3 ) THEN
        IF ( XFLAP(1:2) .EQ. 'VA' ) XFLAP(1:2)='SE'
        IF ( XFLBT(1:2) .EQ. 'VA' ) XFLBT(1:2)='SE'
        IFLBT=IFLAT
        IFLAP=IFLBP
        GOTO 9998
      ENDIF
      IF ( IRET .EQ. 5 ) THEN
        IF ( XFLAP(1:2) .EQ. 'SE'  .AND.  XFLBP(1:2) .EQ. 'SE' ) THEN
          IFLAP=-ABS(IFLAP)
          IFLBP=IFLAP
        ENDIF
        IF ( XFLAT(1:2) .EQ. 'SE'  .AND.  XFLBT(1:2) .EQ. 'SE' ) THEN
          IFLBT=-ABS(IFLBT)
          IFLAT=IFLBT
        ENDIF
        GOTO 9998
      ENDIF

      CALL HASTFC(SILP,SIXP,IRET)
      IF     ( IRET .EQ. 0 ) THEN
        CALL HASTFC(SILT,SIXT,IRET)
      ELSEIF ( IRET .EQ. 1 ) THEN
        IF ( ISH .GE. 91 ) WRITE(IFCH,*)'LIGHT STRING --> REDO HASTFS'
        GOTO 9998
      ELSEIF ( IRET .EQ. 2 ) THEN
        IF ( ISH .GE. 91 ) WRITE(IFCH,*)'JC>9 --> REDO HAPAPA'
        GOTO 9997
      ENDIF

9999  CONTINUE
      IF ( ISKIP .EQ. 1 ) THEN
        COLEVT=COLEVT-1./NCOLEX
        CALL HDECMP(PROJ,SILP,SIXP)
        CALL HDECMP(TARG,SILT,SIXT)
        IF ( ISH .GE. 91 ) THEN
          WRITE(IFCH,*)'SKIP'
          CALL HASTFW(SILP,SIXP)
          CALL HASTFW(SILT,SIXT)
        ENDIF
      ENDIF

      CALL HASTFL(SILP,SIXP,STRLP)
      CALL HASTFL(SILT,SIXT,STRLT)

      IF ( ISKIP .EQ. 1 ) GOTO 9995

      CALL IDDECO(ICVP,JCVP)
      CALL IDDECO(ICVT,JCVT)
      IF ( XFLAP(1:2).EQ.'VA') JCVP(IABS(IFLAP),1)=JCVP(IABS(IFLAP),1)-1
      IF ( XFLBP(1:2).EQ.'VA') JCVP(IABS(IFLBP),2)=JCVP(IABS(IFLBP),2)-1
      IF ( XFLAT(1:2).EQ.'VA') JCVT(IABS(IFLAT),1)=JCVT(IABS(IFLAT),1)-1
      IF ( XFLBT(1:2).EQ.'VA') JCVT(IABS(IFLBT),2)=JCVT(IABS(IFLBT),2)-1
      CALL IDENCO(JCVP,ICVP,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('HAPAPA: IDENCO RET CODE = 1             ')
      ENDIF
      CALL IDENCO(JCVT,ICVT,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('HAPAPA: IDENCO RET CODE = 1             ')
      ENDIF

      PROJ(5,1)=ICVP(1)
      PROJ(6,1)=ICVP(2)
      TARG(5,1)=ICVT(1)
      TARG(6,1)=ICVT(2)
      CALL UTKSIX(SIXP,KMAXPN)
      CALL UTKSIX(SIXT,KMAXTN)
      IF ( KMAXPN+1 .GT. NSIX ) THEN
        CALL UTSTOP('HAHABS: DIMENSION NSIX TOO SMALL        ')
      ENDIF
      IF ( KMAXTN+1 .GT. NSIX ) THEN
        CALL UTSTOP('HAHABS: DIMENSION NSIX TOO SMALL        ')
      ENDIF
      KMXP=KMAXPN+1
      KMXT=KMAXTN+1
      DO 5 N=1,NSI
        PROJ(N,2)=STRLP(N,1)
        TARG(N,2)=STRLT(N,1)
        DO 6 J=1,KMXP
          PROJ(N,2+J)=SIXP(N,J)
 6      CONTINUE
        DO 7 J=1,KMXT
          TARG(N,2+J)=SIXT(N,J)
 7      CONTINUE
        KPROJA(2,IPROJ)=KOL
        IF ( KMAXPN .GT. KMAXP ) THEN
          DO 2 K=KMAXP+1,KMAXPN
            KPROJA(2+K,IPROJ)=KOL
 2        CONTINUE
        ENDIF
        KTARGA(2,ITARG)=KOL
        IF ( KMAXTN .GT. KMAXT ) THEN
          DO 3 K=KMAXT+1,KMAXTN
            KTARGA(2+K,ITARG)=KOL
 3        CONTINUE
        ENDIF
 5    CONTINUE

      NSTSH=0
      LEAD=0
      ISPLT=0
ctp060302 11    CALL HASTPR(STRAP,ISPLT)
      CALL HASTPR(STRAP,ISPLT)
      IF ( ISPLT .NE. 0 ) THEN
        CALL UTSTOP('HAHABS: ISPLT /= 0                      ')
      ENDIF
      ISPLT=0
ctp060302 12    CALL HASTPR(STRAT,ISPLT)
      CALL HASTPR(STRAT,ISPLT)
      IF ( ISPLT .NE. 0 ) THEN
        CALL UTSTOP('HAHABS: ISPLT /= 0                      ')
      ENDIF

      IF ( KMAXPN .GT. KMAXP ) THEN
        PROJ(3,1)=PROJ(3,1)+COORD(3,KOL)
        PROJ(4,1)=PROJ(4,1)+COORD(4,KOL)
      ENDIF
      IF ( KMAXTN .GT. KMAXT ) THEN
        TARG(3,1)=TARG(3,1)+COORD(3,KOL)
        TARG(4,1)=TARG(4,1)+COORD(4,KOL)
      ENDIF
      KMAXP=KMAXPN
      KMAXT=KMAXTN

9995  LEAD=1

C  WRITE LEADING STRING (PROJ)
C  ---------------------------

      IF ( IAP .EQ. 1 ) THEN

        IF ( KMAXP .GT. 0 ) THEN
          COOAV3=PROJ(3,1)/KMAXP
          COOAV4=PROJ(4,1)/KMAXP
        ELSE
          COOAV3=COORD(3,KOL)
          COOAV4=COORD(4,KOL)
        ENDIF
        CALL UTKSTR(STRLP,KMAXOR)
        PZSTRL=0.
        ESTRL=0.
        ISEA=1
        PSUM(1)=0.
        PSUM(2)=0.
        PSUM(3)=0.
        PSUM(4)=0.
        DO 22 K=1,KMAXOR
          PSUM(1)=PSUM(1)+STRLP(1,K)
          PSUM(2)=PSUM(2)+STRLP(2,K)
          PSUM(3)=PSUM(3)+STRLP(3,K)
          PSUM(4)=PSUM(4)+STRLP(4,K)
          AMPR=0.
C-C       IF ( K .EQ. 1 ) AMPR=.94
          PZSTRL=PZSTRL+STRLP(3,K)
          ESTRL=ESTRL+SQRT(STRLP(4,K)**2+AMPR**2)
          IF (K.GE.2.AND.STRLP(5,K).GT.0..AND.STRLP(6,K).GT.0.) ISEA=0
22      CONTINUE
        PSUM(5)=PSUM(4)**2-PSUM(3)**2-PSUM(2)**2-PSUM(1)**2
        IF ( PSUM(5) .GT. 0. ) PSUM(5)=SQRT(PSUM(5))
        IF ( PZSTRL/PNLLX.GT.1.-0.850/ENGY**2.AND.ISKIP.NE.1 ) GOTO 1002
        ISTRL=0
        IF ( NSTSH .EQ. 0 .AND. ISEA .EQ. 1  .AND.
     *                     PZSTRL/PNLLX .GT. 1.-10.000/ENGY**2 ) ISTRL=1
        DO 23 K=1,KMAXOR
          DO 23 I=1,NSI
            STR0(I,K)=STRLP(I,K)
23      CONTINUE
        NSTR0=NSTR
17      NSTR=NSTR0
        ISPLT=0
15      CONTINUE
        DO 24 K=1,KMAXOR
          DO 24 I=1,NSI
            STRLP(I,K)=STR0(I,K)
24      CONTINUE
13      CALL HASTPR(STRLP,ISPLT)
        IF     ( ISPLT .EQ. -1 ) THEN
          IF ( ISKIP .EQ. 1 ) GOTO 9996
          GOTO 9994
        ELSEIF ( ISPLT .EQ. -3 ) THEN
          GOTO 15
        ELSEIF ( ISPLT .EQ. -4 ) THEN
          GOTO 1001
        ELSEIF ( ISPLT .EQ. -5 ) THEN
          GOTO 17
        ELSEIF ( ISPLT .GT. 0 ) THEN
          GOTO 13
        ENDIF
        IF ( NSTR .GT. NSTR0+1 ) THEN
          PSUMX(1)=0.
          PSUMX(2)=0.
          PSUMX(3)=0.
          PSUMX(4)=0.
          DO 25 J=NSTR0+1,NSTR
            PSUMX(1)=PSUMX(1)+PSTR(1,J)
            PSUMX(2)=PSUMX(2)+PSTR(2,J)
            PSUMX(3)=PSUMX(3)+PSTR(3,J)
            PSUMX(4)=PSUMX(4)+PSTR(4,J)
25        CONTINUE
          PSUMX(5)=PSUMX(4)**2-PSUMX(3)**2-PSUMX(2)**2-PSUMX(1)**2
          IF ( PSUMX(5) .GT. 0. ) PSUMX(5)=SQRT(PSUMX(5))
C-C       WRITE(6,*)' '
C-C       DO 26 J=NSTR0+1,NSTR
C-C       WRITE(6,109) J,(ICSTR(K,J)/100,K=1,4)
C-C  *      ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J)/PNLLX,PSTR(5,J)
C-C  *      ,IRLSTR(J)
C-C26     CONTINUE
C-C       WRITE(6,*)PSUM
C-C       WRITE(6,*)PSUMX
          DO 27 J=NSTR0+1,NSTR
            PSTR(1,J+NSTR-NSTR0)=PSTR(1,J)
            PSTR(2,J+NSTR-NSTR0)=PSTR(2,J)
            PSTR(3,J+NSTR-NSTR0)=PSTR(3,J)
            PSTR(4,J+NSTR-NSTR0)=PSTR(4,J)
27        CONTINUE
          CALL HRESCL(NSTR0+1,NSTR,PSUM,IFAIL)
          IF ( IFAIL .GT. 0 ) THEN
C-C         WRITE(6,*)'RESCALE FAILED'
            DO 28 J=NSTR0+1,NSTR
              PSTR(1,J)=PSTR(1,J+NSTR-NSTR0)
              PSTR(2,J)=PSTR(2,J+NSTR-NSTR0)
              PSTR(3,J)=PSTR(3,J+NSTR-NSTR0)
              PSTR(4,J)=PSTR(4,J+NSTR-NSTR0)
28          CONTINUE
          ENDIF
          PSUMX(1)=0.
          PSUMX(2)=0.
          PSUMX(3)=0.
          PSUMX(4)=0.
          DO 29 J=NSTR0+1,NSTR
            PSUMX(1)=PSUMX(1)+PSTR(1,J)
            PSUMX(2)=PSUMX(2)+PSTR(2,J)
            PSUMX(3)=PSUMX(3)+PSTR(3,J)
            PSUMX(4)=PSUMX(4)+PSTR(4,J)
29        CONTINUE
          PSUMX(5)=PSUMX(4)**2-PSUMX(3)**2-PSUMX(2)**2-PSUMX(1)**2
          IF ( PSUMX(5) .GT. 0. ) PSUMX(5)=SQRT(PSUMX(5))
C-C       DO 30 J=NSTR0+1,NSTR
C-C       WRITE(6,109) J,(ICSTR(K,J)/100,K=1,4)
C-C  *      ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J)/PNLLX,PSTR(5,J)
C-C  *      ,IRLSTR(J)
C-C30     CONTINUE
C-C       WRITE(6,*)PSUM
C-C       WRITE(6,*)PSUMX
C-C109       FORMAT(' /CSTR/',I4,3X,4I5,2X,3(E10.3),I4)
        ENDIF

      ENDIF

C  WRITE LEADING STRING (TARG)
C  ---------------------------

      IF ( IAT .EQ. 1 ) THEN

        IF ( KMAXT .GT. 0 ) THEN
          COOAV3=TARG(3,1)/KMAXT
          COOAV4=TARG(4,1)/KMAXT
        ELSE
          COOAV3=COORD(3,KOL)
          COOAV4=COORD(4,KOL)
        ENDIF
        CALL UTKSTR(STRLT,KMAXOR)
        PZSTRL=0.
        ESTRL=0.
        ISEA=1
        PSUM(1)=0.
        PSUM(2)=0.
        PSUM(3)=0.
        PSUM(4)=0.
        DO 33 K=1,KMAXOR
          PSUM(1)=PSUM(1)+STRLT(1,K)
          PSUM(2)=PSUM(2)+STRLT(2,K)
          PSUM(3)=PSUM(3)+STRLT(3,K)
          PSUM(4)=PSUM(4)+STRLT(4,K)
          AMPR=0.
C-C       IF ( K .EQ. 1 ) AMPR=.94
          PZSTRL=PZSTRL+STRLT(3,K)
          ESTRL=ESTRL+SQRT(STRLT(4,K)**2+AMPR**2)
          IF ( K.GE.2.AND.STRLT(5,K).GT.0..AND.STRLT(6,K).GT.0. ) ISEA=0
33      CONTINUE
        PSUM(5)=PSUM(4)**2-PSUM(3)**2-PSUM(2)**2-PSUM(1)**2
        IF ( PSUM(5) .GT. 0. ) PSUM(5)=SQRT(PSUM(5))
        IF ( -PZSTRL/PNLLX.GT.1.-0.850/ENGY**2.AND.ISKIP.NE.1 )GOTO 1002
        ISTRL=0
        IF ( NSTSH .EQ. 0  .AND.  ISEA .EQ. 1 .AND.
     *                 -PZSTRL/PNLLX .GT. 1.-10.000/ENGY**2 ) ISTRL=1
        DO 34 K=1,KMAXOR
          DO 34 I=1,NSI
            STR0(I,K)=STRLT(I,K)
34      CONTINUE
        NSTR0=NSTR
18      NSTR=NSTR0
        ISPLT=0
16      CONTINUE
        DO 35 K=1,KMAXOR
          DO 35 I=1,NSI
            STRLT(I,K)=STR0(I,K)
35      CONTINUE
14      CALL HASTPR(STRLT,ISPLT)
        IF     ( ISPLT .EQ. -1 ) THEN
           IF ( ISKIP .EQ. 1 ) GOTO 9996
           GOTO 9994
        ELSEIF ( ISPLT .EQ. -3 ) THEN
          GOTO 16
        ELSEIF ( ISPLT .EQ. -4 ) THEN
          GOTO 1001
        ELSEIF ( ISPLT .EQ. -5 ) THEN
          GOTO 18
        ELSEIF ( ISPLT .GT.  0 ) THEN
          GOTO 14
        ENDIF
        IF ( NSTR .GT. NSTR0+1 ) THEN
          PSUMX(1)=0.
          PSUMX(2)=0.
          PSUMX(3)=0.
          PSUMX(4)=0.
          DO 36 J=NSTR0+1,NSTR
            PSUMX(1)=PSUMX(1)+PSTR(1,J)
            PSUMX(2)=PSUMX(2)+PSTR(2,J)
            PSUMX(3)=PSUMX(3)+PSTR(3,J)
            PSUMX(4)=PSUMX(4)+PSTR(4,J)
36        CONTINUE
          PSUMX(5)=PSUMX(4)**2-PSUMX(3)**2-PSUMX(2)**2-PSUMX(1)**2
          IF ( PSUMX(5) .GT. 0. ) PSUMX(5)=SQRT(PSUMX(5))
C-C       WRITE(6,*)' '
C-C       DO 37 J=NSTR0+1,NSTR
C-C       WRITE(6,109) J,(ICSTR(K,J)/100,K=1,4)
C-C  *      ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J)/PNLLX,PSTR(5,J)
C-C  *      ,IRLSTR(J)
C-C37     CONTINUE
C-C       WRITE(6,*)PSUM
C-C       WRITE(6,*)PSUMX
          DO 38 J=NSTR0+1,NSTR
            PSTR(1,J+NSTR-NSTR0)=PSTR(1,J)
            PSTR(2,J+NSTR-NSTR0)=PSTR(2,J)
            PSTR(3,J+NSTR-NSTR0)=PSTR(3,J)
            PSTR(4,J+NSTR-NSTR0)=PSTR(4,J)
38        CONTINUE
          CALL HRESCL(NSTR0+1,NSTR,PSUM,IFAIL)
          IF ( IFAIL .GT. 0 ) THEN
C-C         WRITE(6,*)'RESCALE FAILED'
            DO 39 J=NSTR0+1,NSTR
              PSTR(1,J)=PSTR(1,J+NSTR-NSTR0)
              PSTR(2,J)=PSTR(2,J+NSTR-NSTR0)
              PSTR(3,J)=PSTR(3,J+NSTR-NSTR0)
              PSTR(4,J)=PSTR(4,J+NSTR-NSTR0)
39          CONTINUE
          ENDIF
          PSUMX(1)=0.
          PSUMX(2)=0.
          PSUMX(3)=0.
          PSUMX(4)=0.
          DO 40 J=NSTR0+1,NSTR
            PSUMX(1)=PSUMX(1)+PSTR(1,J)
            PSUMX(2)=PSUMX(2)+PSTR(2,J)
            PSUMX(3)=PSUMX(3)+PSTR(3,J)
            PSUMX(4)=PSUMX(4)+PSTR(4,J)
40        CONTINUE
          PSUMX(5)=PSUMX(4)**2-PSUMX(3)**2-PSUMX(2)**2-PSUMX(1)**2
          IF ( PSUMX(5) .GT. 0. ) PSUMX(5)=SQRT(PSUMX(5))
C-C       DO 41 J=NSTR0+1,NSTR
C-C       WRITE(6,109) J,(ICSTR(K,J)/100,K=1,4)
C-C  *      ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J)/PNLLX,PSTR(5,J)
C-C  *      ,IRLSTR(J)
C-C41     CONTINUE
C-C       WRITE(6,*)PSUM
C-C       WRITE(6,*)PSUMX
        ENDIF

      ENDIF

C  EXIT
C  ----

      IF ( ISH .LT. 91 ) RETURN

      WRITE(IFCH,100) ( (PROJ(I,J),I=1,NSI), J=2,KMAXP+2 )
      WRITE(IFCH,102) ( (TARG(I,J),I=1,NSI), J=2,KMAXT+2 )
ctp060203 101   FORMAT ( '  ',I3,'. TRIAL')
      RETURN

9996  CONTINUE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'SKIP IMPOSSIBLE, STRL CANT BE STORED'
      ENDIF
      ISKIP=2
      RETURN

1001  CONTINUE
        IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)
     *      'SR HAHABS: NOT JUSTIFIED COMPLETE ABSORPTION -> IRETHH=1'
      ENDIF
      IRETHH=1
      IF ( ISKIP .EQ. 1 ) GOTO 9996
      RETURN

1002  CONTINUE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'SR HAHABS: TOO FAST MULTI-STRING -> IRETHH=1'
      ENDIF
      IRETHH=1
      IF ( ISKIP .EQ. 1 ) THEN
        CALL UTSTOP('HAHABS: ISKIP=1 SHOULD NOT HAPPEN       ')
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE HAPAPA(SIL,IFLA,XFLA,PTA,PTAH,IFLB,XFLB,PTB,PTBH,NQAQ
     *,ICV)

C-----------------------------------------------------------------------
C  DETERMINES MOMENTA AND FLAVOR OF PARTICIPATING PARTONS IN A HADRON
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      PARAMETER (NPTQ=129)
      PARAMETER (NSI=6)
      PARAMETER (NSTRU=2049)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CNEW/    KOTRI,NEWCOL,NEWICO
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTQ/    QPTH(NPTQ),QPTQ(NPTQ),XPTQ(NPTQ),QPTQMX,QPTHMX
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      PROBAB(NFLAV),PROBSU(NFLAV+1),SIL(NSI)
      INTEGER   ICV(2),ID(2),JC(NFLAV,2),JCV(NFLAV,2)
      CHARACTER CHPT*1,XFLA*3,XFLB*3
      SAVE
C-----------------------------------------------------------------------
      SGN=SIGN(1.,SIL(3))
      IF     ( SGN .GT. 0. ) THEN
        QVALC=QVAPC
        QSEAC=QSEPC
        CHPT='P'
        W=WPROJ
      ELSE
        QVALC=QVATC
        QSEAC=QSETC
        CHPT='T'
        W=WTARG
      ENDIF
      CALL IDDECO(ICV,JCV)
      NVQ=0
      NVA=0
      DO 12 I=1,NFLAV
        IF ( ISH.GE.0  .AND.  (JCV(I,1).LT.0.OR.JCV(I,2).LT.0) ) THEN
          CALL UTMSG('HAHABS')
          WRITE(IFCH,*)'*****  NEGATIVE JCV'
          WRITE(IFCH,*)'JCV:'
          WRITE(IFCH,*)JCV
          CALL UTMSGF
        ENDIF
        NVQ=NVQ+JCV(I,1)
        NVA=NVA+JCV(I,2)
12    CONTINUE
      ID(1)=NINT(SIL(4+1))
      ID(2)=NINT(SIL(4+2))
      CALL IDDECO(ID,JC)
      NQU=0
      NAQ=0
      DO 11 I=1,NFLAV
        NQU=NQU+JC(I,1)
        NAQ=NAQ+JC(I,2)
11    CONTINUE
      NEQ=NQU-NAQ

C  DETERMINE NQAQ,IVAL1,IVAL2
C  --------------------------
      NQAQ=0
      IVAL=0
      IF ( RANGEN() .GT. W ) THEN
        IF ( NEQ .GT. NEQMN ) THEN
          PQ=1.
        ELSE
          PQ=0.
        ENDIF
        IF ( NEQ .LT. NEQMX ) THEN
          PA=1.
        ELSE
          PA=0.
        ENDIF
        SUM=(NVQ*PQ+NVA*PA*IAQU)*QVALC+QSEAC*.5*(PQ+PA*IAQU)
        IF ( SUM .EQ. 0. ) GOTO 14
        SUMI = 1./SUM
        PVQ=NVQ*PQ*QVALC*SUMI
        PVA=NVA*PA*IAQU*QVALC*SUMI
        PSQ=.5*PQ*QSEAC*SUMI
        PSA=.5*PA*IAQU*QSEAC*SUMI
        R=RANGEN()
        IF ( R .LT. PVQ+PVA ) IVAL=1
        NQAQ=-1
        IF ( IVAL .EQ. 1  .AND.  R .LT. PVQ ) NQAQ=1
        IF ( IVAL .EQ. 0  .AND.  R .LT. PVQ+PVA+PSQ ) NQAQ=1
        IF ( NQU-NAQ-NQAQ .LT. NEQMN ) THEN
          IF ( ISH .GE. 90 ) THEN
            CALL UTMSG('HAPAPA')
            WRITE(IFCH,*)'*****  NEQ < NEQMN'
            WRITE(IFCH,*)'PVQ/A PSQ/A:',PVQ,PVA,PSQ,PSA
            WRITE(IFCH,*)'NQAQ:',NQAQ
            CALL UTMSGF
          ENDIF
          NQAQ=0
        ENDIF
        IF ( NQU-NAQ-NQAQ .GT. NEQMX ) THEN
          IF ( ISH .GE. 90 ) THEN
            IF ( ISH .GE. 91 ) WRITE(IFCH,*)' '
            CALL UTMSG('HAPAPA')
            WRITE(IFCH,*)'*****  NEQ > NEQMX'
            WRITE(IFCH,*)'PVQ/A PSQ/A:',PVQ,PVA,PSQ,PSA
            WRITE(IFCH,*)'NQAQ:',NQAQ
            IF ( ISH .GE. 91 ) WRITE(IFCH,*)' '
          ENDIF
          NQAQ=0
        ENDIF
        IF ( NQU-NAQ-NQAQ .LT. NEQMN ) THEN
          CALL UTSTOP('HAPAPA: NEQ.LT.NEQMN                    ')
        ENDIF
        IF ( NQU-NAQ-NQAQ .GT. NEQMX ) THEN
          CALL UTSTOP('HAPAPA: NEQ.GT.NEQMX                    ')
        ENDIF
      ENDIF
14    CONTINUE
      IVAL1=IVAL
      IVAL2=IVAL
      IF ( NQAQ .EQ. 0 ) THEN
        IVAL1=0
        SUM=NVQ*QVALC+QSEAC*.5
        IF ( SUM .EQ. 0. ) GOTO 15
        SUMI = 1./SUM
        PVQ=NVQ*QVALC*SUMI
        PSQ=.5*QSEAC*SUMI
        R=RANGEN()
        IF ( R .LT. PVQ ) IVAL1=1
15      CONTINUE
C-C     IF ( NVQ.GT.0  .AND.  RANGEN().LT.PVALEN ) IVAL1=1
        IVAL2=0
        SUM=NVA*IAQU*QVALC+QSEAC*.5*IAQU
        IF ( SUM .EQ. 0. ) GOTO 16
        SUMI = 1./SUM
        PVA=NVA*IAQU*QVALC*SUMI
        PSA=.5*IAQU*QSEAC*SUMI
        R=RANGEN()
        IF ( R .LT. PVA ) IVAL2=1
16      CONTINUE
C-C     IF ( NVA.GT.0 .AND. RANGEN().LT.PVALEN ) IVAL2=1
        IF ( IVAL1 .EQ. 1  .AND. IVAL2 .EQ. 1 ) THEN
          R=RANGEN()
          IF ( R .LT. 0.5 ) THEN
            IVAL1=0
          ELSE
            IVAL2=0
          ENDIF
        ENDIF
      ENDIF

C  QUARK
C  -----
      XFLA='---'
      IF ( NQAQ .GE. 0 ) THEN
        IF ( IVAL1 .EQ. 1 ) THEN
          PROBAB(1)=JCV(1,1)
          PROBAB(2)=JCV(2,1)
          PROBAB(3)=JCV(3,1)
          PROBAB(4)=JCV(4,1)
          SU=PROBAB(1)+PROBAB(2)+PROBAB(3)+PROBAB(4)
          XFLA='VA'//CHPT
        ELSE
          PROBAB(1)=1.
          PROBAB(2)=1.
          PROBAB(3)=RSTRAS
          PROBAB(4)=0.
          SU=2.+PROBAB(3)
          XFLA='SE'//CHPT
        ENDIF
        PROBSU(1)=0.
        PROBSU(2)=1./SU*PROBAB(1)
        PROBSU(3)=1./SU*PROBAB(2)+PROBSU(2)
        PROBSU(4)=1./SU*PROBAB(3)+PROBSU(3)
        PROBSU(5)=1./SU*PROBAB(4)+PROBSU(4)
        R=RANGEN()
        IF     ( R .LE. PROBSU(2) ) THEN
          IFLA=1
        ELSEIF ( R .LE. PROBSU(3) ) THEN
          IFLA=2
        ELSEIF ( R .LE. PROBSU(4) ) THEN
          IFLA=3
        ELSE
          IFLA=4
        ENDIF
        R=RANGEN()
        IF     ( IOPTQ .EQ. 2 ) THEN
          PTA = SQRT( -4.*PTQ**2/PI * LOG(1.-R*QPTQMX) )
        ELSEIF ( IOPTQ .EQ. 3 ) THEN
          PTA = PTQ*SQRT( QPTQMX*R/(1.-QPTQMX*R) )
        ELSE
          PTA=UTINVT(NPTQ,XPTQ,QPTQ,R*QPTQ(NPTQ))
        ENDIF
        PTAH=PTH*SQRT( 1./SQRT(1.-2.*PTH**2*RANGEN()*QPTHMX) - 1. )
      ENDIF

C  ANTIQUARK
C  ---------
      XFLB='---'
      IF ( NQAQ .LE. 0 ) THEN
        IF ( IVAL2 .EQ. 1 ) THEN
          PROBAB(1)=JCV(1,2)
          PROBAB(2)=JCV(2,2)
          PROBAB(3)=JCV(3,2)
          PROBAB(4)=JCV(4,2)
          SU=PROBAB(1)+PROBAB(2)+PROBAB(3)+PROBAB(4)
          XFLB='VA'//CHPT
        ELSE
          PROBAB(1)=1.
          PROBAB(2)=1.
          PROBAB(3)=RSTRAS
          SU=2.+RSTRAS
          PROBAB(4)=0.
          XFLB='SE'//CHPT
        ENDIF
        PROBSU(1)=0.
        PROBSU(2)=1./SU*PROBAB(1)
        PROBSU(3)=1./SU*PROBAB(2)+PROBSU(2)
        PROBSU(4)=1./SU*PROBAB(3)+PROBSU(3)
        PROBSU(5)=1./SU*PROBAB(4)+PROBSU(4)
        R=RANGEN()
        IF     ( R .LE. PROBSU(2) ) THEN
          IFLB=1
        ELSEIF ( R .LE. PROBSU(3) ) THEN
          IFLB=2
        ELSEIF ( R .LE. PROBSU(4) ) THEN
          IFLB=3
        ELSE
          IFLB=4
        ENDIF
        IF ( NQAQ.EQ.0 .AND. IVAL1.EQ.0 .AND. IVAL2.EQ.0 ) IFLB=IFLA
C-C     IF ( NQAQ.EQ.0 .AND. IVAL2.EQ.0 ) IFLB=IFLA
C-C     IF ( NQAQ.EQ.0 .AND. IVAL1.EQ.0 ) IFLA=IFLB
        R=RANGEN()
        IF     ( IOPTQ .EQ. 2 ) THEN
          PTB = SQRT( -4.*PTQ**2/PI * LOG(1.-QPTQMX*R) )
        ELSEIF ( IOPTQ .EQ. 3 ) THEN
          PTB = PTQ*SQRT( QPTQMX*R/(1.-QPTQMX*R) )
        ELSE
          PTB=UTINVT(NPTQ,XPTQ,QPTQ,R*QPTQ(NPTQ))
        ENDIF
        PTBH=PTH*SQRT( 1./SQRT(1.-2.*PTH**2*RANGEN()*QPTHMX) - 1. )
      ENDIF

      IF ( NQAQ.EQ.0 .AND. IVAL1.EQ.0 .AND. IVAL2.EQ.0 ) THEN
        IFLA=-IFLA
        IFLB=-IFLB
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE HASI(QAQ,SIL,IFL,XFL,PT0,PT0H,SI,IRET,JORD,IXFLAB
     *               ,PTDIFF)

C-----------------------------------------------------------------------
C  DETERMINES STRING INGREDIENTS (=JETS)
C  IRET=0: OK
C  IRET=1: REMNANT CHANGES DIRECTION
C  IRET=2: JC(,)=10
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      PTDIFF(2),SI(NSI),SIL(NSI)
      INTEGER   IC(2),ID(2),JC(NFLAV,2)
      CHARACTER QAQ*5,XFL*3
      SAVE
C-----------------------------------------------------------------------
      IRET=0
      IRETX=0
      NTRY=0
      PTDIFF(1)=0.
      PTDIFF(2)=0.
      IFLA=ABS(IFL)
      IF ( XFL .EQ. '000' ) THEN
        DO 11 J=3,NSI
          SI(J)=0.
11      CONTINUE
      ELSE
        DO 10 J=1,NSI
          SI(J)=0.
10      CONTINUE
      ENDIF
      IF ( XFL .EQ. '---' ) GOTO 1000
      IF      ( JORD .EQ. 1 ) THEN
        IF     ( QAQ .EQ. 'QUARK' ) THEN
          AMS=AMPROJ
        ELSEIF ( QAQ .EQ. 'ANTIQ' ) THEN
          AMS=AMTARG
        ENDIF
      ELSEIF ( JORD .EQ. 2 ) THEN
        IF     ( QAQ .EQ. 'ANTIQ' ) THEN
          AMS=AMPROJ
        ELSEIF ( QAQ .EQ. 'QUARK' ) THEN
          AMS=AMTARG
        ENDIF
      ENDIF
C-C   ENLL=SQRT(AMS**2+PNLLX**2)
      ENLL=SIL(4)
      SGN=SIGN(1.,SIL(3))
      IF ( RANGEN() .LE. PHARD ) THEN
        IHARD=1
      ELSE
        IHARD=0
      ENDIF
      ID(1)=NINT(SIL(4+1))
      ID(2)=NINT(SIL(4+2))
      CALL IDDECO(ID,JC)
9999  NTRY=NTRY+1
      IF ( NTRY .GT. 20 ) THEN
        IRET=1
        GOTO 1000
      ENDIF
      IF     ( IHARD .EQ. 0 ) THEN
        PT=PT0
        IF ( PT .GT. ENLL ) PT=RANGEN()*ENLL
        ENMIN=0.
      ELSE
        PT=PT0H
        IF ( PT .GT. ENLL ) PT=RANGEN()*ENLL
        ENMIN=PT
      ENDIF
      PTFULL=PT
      LO=0
51    LO=LO+1
      IF ( XFL .EQ. '000' ) THEN
        EN=0.
        PT=0.
      ELSE
        EN=RANSTC(XFL,ENMIN/ENLL)*ENLL
      ENDIF
      IF ( PT .GT. EN ) THEN
        IF ( RANGEN() .LT. 0.5  .AND.  LO .LE. 10 ) GOTO 51
        PT=RANGEN()*EN
      ENDIF
      PHI=2.*PI*RANGEN()
      IF ( XFL .NE. '000' ) THEN
        SI(1)=PT*COS(PHI)
        SI(2)=PT*SIN(PHI)
        SI(3)=SGN*SQRT(EN**2-PT**2)
        IF ( SI(1).EQ.0.  .AND.  SI(2).EQ.0.  .AND.  SI(3).EQ.0. ) THEN
          IRETX=1
          IRET=1
          GOTO 1000
        ENDIF
C-C     PTDIFF(1)=(PTFULL-PT)*COS(PHI)
C-C     PTDIFF(2)=(PTFULL-PT)*SIN(PHI)
      ENDIF
      IF ( (SIL(3)-SI(3))*SIL(3) .LT. 0. ) THEN
        IF ( ISH .GE. 94 ) WRITE(IFCH,*)'SIL,SI,ENLL,EN',SIL,SI,ENLL,EN
        GOTO 9999
      ENDIF
      SI(4)=SQRT(SI(1)**2+SI(2)**2+SI(3)**2)
      IF     ( QAQ .EQ. 'QUARK' ) THEN
        SI(5)=10**(NFLAV-IFLA)
      ELSEIF ( QAQ .EQ. 'ANTIQ' ) THEN
        SI(6)=10**(NFLAV-IFLA)
      ENDIF

      IF ( IXFLAB .EQ. 0 ) THEN
        SIL(1)=SIL(1)-SI(1)
        SIL(2)=SIL(2)-SI(2)
        SIL(3)=SIL(3)-SI(3)
        SIL(4)=SQRT(SIL(1)**2+SIL(2)**2+SIL(3)**2)
      ENDIF

      IF ( JC(IFLA,1) .GT. 0 ) THEN
        L1=1
      ELSE
        L1=0
      ENDIF
      IF ( JC(IFLA,2) .GT. 0 ) THEN
        L2=1
      ELSE
        L2=0
      ENDIF
      IF     ( QAQ .EQ. 'QUARK' ) THEN
        IF     ( L1 .EQ. 0 ) THEN
          JC(IFLA,2)=JC(IFLA,2)+1
        ELSE
          JC(IFLA,1)=JC(IFLA,1)-1
        ENDIF
      ELSEIF ( QAQ .EQ. 'ANTIQ' ) THEN
        IF     ( L2 .EQ. 0 ) THEN
          JC(IFLA,1)=JC(IFLA,1)+1
        ELSE
          JC(IFLA,2)=JC(IFLA,2)-1
        ENDIF
      ENDIF
      IF ( JC(IFLA,1) .EQ. 10  .OR.  JC(IFLA,2) .EQ. 10 ) GOTO 9998
      CALL IDENCO(JC,IC,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('HASI  : IDENCO RET CODE = 1             ')
      ENDIF
      SIL(5)=IC(1)
      SIL(6)=IC(2)
      GOTO 1000

9998  IRET=2
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('HASI  ')
        WRITE(IFCH,*)'*****  JC(,)=10'
        WRITE(IFCH,*)JC
        CALL UTMSGF
      ENDIF

1000  CONTINUE
      IF ( ISH .GE. 90  .AND. IRETX .EQ. 1 ) THEN
        CALL UTMSG('HASI  ')
        WRITE(IFCH,*)'*****  SI(1/2/3)=0'
        WRITE(IFCH,*)SI
        CALL UTMSGF
      ENDIF
      IF ( ISH .GE. 93 ) THEN
        IF ( IRET .NE. 0  .OR.  NTRY .GT. 1 )
     *           WRITE(IFCH,*)'IRET=',IRET,'  NTRY=',NTRY
        IF     ( NTRY .GT. 0 ) THEN
          WRITE(IFCH,100)XFL,PT0,PT,EN,EN/ENLL
100       FORMAT(' HASI: XFL=',A3
     *          ,' PT0=',E10.3,' PT=',E10.3,' EN=',E10.3,' X=',E10.3)
        ELSEIF ( NTRY .EQ. 0 ) THEN
          WRITE(IFCH,101)XFL
101       FORMAT(1X,'HASI: XFL=',A3)
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTFC(SIL,SIX,IRET)

C-----------------------------------------------------------------------
C  CHECKS LONG STRINGS
C  IRET=0: OK
C  IRET=1: MASS .LT. MINIMAL MASS
C  IRET=2: JC(,) .GT. 9
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL    SI(NSI),SIL(NSI),SIX(NSI,NSIX)
      INTEGER IC(2),JC(NFLAV,2),JCP(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IRET=0

      CALL UTKSIX(SIX,KMAX)
      IF ( KMAX .LT. 1 ) GOTO 9000

      SI(1)=SIL(1)
      SI(2)=SIL(2)
      SI(3)=SIL(3)
      SI(4)=SIL(4)
      DO 110 K=1,KMAX
        SI(1)=SI(1)+SIX(1,K)
        SI(2)=SI(2)+SIX(2,K)
        SI(3)=SI(3)+SIX(3,K)
        SI(4)=SI(4)+SIX(4,K)
110   CONTINUE
      S=SI(4)**2-SI(1)**2-SI(2)**2-SI(3)**2

      IC(1)=NINT(SIL(5))
      IC(2)=NINT(SIL(6))
      CALL IDDECO(IC,JC)
      DO 130 K=1,KMAX
        IC(1)=NINT(ABS(SIX(5,K)))
        IC(2)=NINT(ABS(SIX(6,K)))
        CALL IDDECO(IC,JCP)
        DO 120 I=1,NFLAV
          JC(I,1)=JC(I,1)+JCP(I,1)
          JC(I,2)=JC(I,2)+JCP(I,2)
120     CONTINUE
130   CONTINUE

      DO 125 I=1,NFLAV
        IF ( JC(I,1)-JC(I,2) .GE.  10 ) IRET=2
        IF ( JC(I,1)-JC(I,2) .LE. -10 ) IRET=2
125   CONTINUE

      IF ( IRET .EQ. 2 ) GOTO 9000

      KEU=JC(1,1)-JC(1,2)
      KED=JC(2,1)-JC(2,2)
      KES=JC(3,1)-JC(3,2)
      KEC=JC(4,1)-JC(4,2)
      AMIN=UTAMNU(KEU,KED,KES,KEC,4)*0.5

C-C   IF ( S .LT. AMIN**2 ) IRET=1

9000  CONTINUE

      IF ( ISH .GE. 94 ) THEN
        WRITE(IFCH,*)' '
        IF ( IRET .EQ. 0 ) THEN
          WRITE(IFCH,*)('-',L=1,79)
        ELSE
          WRITE(IFCH,*)('#',L=1,79)
        ENDIF
        WRITE(IFCH,*)'IRET= ',IRET ,'  KMAX= ',KMAX,'  NREVT= ',NREVT
     *        ,'  NCES= ',NCES
        WRITE(IFCH,8004)SIL
8004    FORMAT(' SIL: ',4F13.5,2F8.0)
        IF ( KMAX .GT. 0 ) THEN
          WRITE(IFCH,8007)(SIX(I,1),I=1,NSI)
8007      FORMAT(' SIX: ',4F13.5,2F8.0)
          IF ( KMAX .GT. 1 ) THEN
            DO 1 J=2,KMAX
              WRITE(IFCH,8008)(SIX(I,J),I=1,NSI)
 1          CONTINUE
8008        FORMAT('      ',4F13.5,2F8.0)
          ENDIF
          WRITE(IFCH,8005)(SI(I),I=1,4)
8005      FORMAT('  SI: ',4F13.5)
          WRITE(IFCH,8006)(SI(I)**2,I=1,4)
8006      FORMAT(' SI^2:',4F13.5)
          WRITE(IFCH,*)'JC:'
          WRITE(IFCH,*)JC
          WRITE(IFCH,*)'KEU,KED,KES,KEC: ',KEU,KED,KES,KEC
          WRITE(IFCH,*)'S= ',S,'    AMIN**2= ',AMIN**2,'   AMIN= ',AMIN
        ENDIF
        IF ( IRET .NE. 0 ) THEN
          WRITE(IFCH,*)('#',L=1,79)
        ELSE
          WRITE(IFCH,*)('-',L=1,79)
        ENDIF
        WRITE(IFCH,*)' '
      ENDIF

      IF ( IRET .NE. 0 ) RETURN
      IF ( ISH .GE. 91 ) CALL HASTFW (SIL,SIX)
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTFL(SIL,SIX,STRL)

C-----------------------------------------------------------------------
C  FORMS A LEADING STRING
C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)

      REAL SIL(NSI),SIX(NSI,NSIX),STRL(NSI,NSIX+1)
      SAVE
C-----------------------------------------------------------------------
      DO 100 N=1,NSI
        STRL(N,1)=SIL(N)
        STRL(N,2)=0.
100   CONTINUE
      CALL UTKSIX(SIX,KMAX)
      IF ( KMAX .EQ. 0 ) RETURN
      DO 111 K=1,KMAX
        DO 110 N=1,NSI
          STRL(N,K+1)=SIX(N,K)
110     CONTINUE
111   CONTINUE
      DO 120 N=1,NSI
        STRL(N,KMAX+2)=0.
120   CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTFR(SILA,SILB,SIA,SIB,IRET)

C-----------------------------------------------------------------------
C  SUBTRACTS SIA-SIB FROM SILA-SILB (ONLY MOMENTUM COMPONENTS)
C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION A,B,D,DAUXIL,PAB(4),PAM,PAP,PA3,PA4
     *                ,PEM,PEP,PIM,PIP,PI3,PI4,PLAB(4)
     *                ,POT,POX,POY,PUT,PUX,PUY,PYM,PYP,SLA
     *                ,SSIA(NSI),SSIB(NSI),SSILA(NSI),SSILB(NSI)
      REAL             SIA(NSI),SIB(NSI),SILA(NSI),SILB(NSI)
      SAVE
C-----------------------------------------------------------------------
      IRET=0
      DO 12 I=1,4
        SSILB(I)=SILB(I)
        SSILA(I)=SILA(I)
        PLAB(I)=SSILA(I)+SSILB(I)
        SSIB(I)=SIB(I)
        SSIA(I)=SIA(I)
        PAB(I)=SSIA(I)+SSIB(I)
12    CONTINUE
      SLA=SIGN(1.D0, SSILA(3))

      IF ( ISH.GE.0 .AND. SSIA(3)*SSILA(3) .LT. 0.D0 ) THEN
        CALL UTMSG('HASTFR')
        WRITE(IFCH,*)'*****  SSIA(3)*SSILA(3)<0'
        WRITE(IFCH,*)SSIA(3),SSILA(3)
        CALL UTMSGF
      ENDIF

      A=0.D0
      D=0.D0
      PYP=0.D0
      PYM=0.D0
      PIP=0.D0
      PAM=0.D0

      POX=SSILA(1)-SSIA(1)
      POY=SSILA(2)-SSIA(2)
      POT=SQRT(POX**2+POY**2)
      PUX=SSILB(1)-SSIB(1)
      PUY=SSILB(2)-SSIB(2)
      PUT=SQRT(PUX**2+PUY**2)
      DAUXIL = SLA*(PLAB(3)-PAB(3))
      PEP= PLAB(4)-PAB(4) + DAUXIL
      PEM= PLAB(4)-PAB(4) - DAUXIL
      IF ( PEP .LT. 0.D0 ) GOTO 9001
      IF ( PEM .LT. 0.D0 ) GOTO 9001
      A=(PEM*PEP-PUT**2-POT**2)*0.5D0
      D=PUT*POT
      DAUXIL = A**2 - D**2
      IF ( DAUXIL .LT. 0.D0 ) GOTO 9001
      B=SQRT(DAUXIL)
      PYP=(A+PUT**2-B)/PEM
      PYM=(A+POT**2-B)/PEP
      IF ( PYP .LT. 0.D0 ) GOTO 9001
      IF ( PYM .LT. 0.D0 ) GOTO 9001
      PIP=PEP-PYP
      IF ( PIP .LT. 0.D0 ) GOTO 9001
      PIM=PYM
      PAP=PYP
      PAM=PEM-PYM
      IF ( PAM .LT. 0.D0 ) GOTO 9001
      PI3=(PIP-PIM)*0.5D0*SLA
      PI4=(PIP+PIM)*0.5D0
      PA3=(PAP-PAM)*0.5D0*SLA
      PA4=(PAP+PAM)*0.5D0
      IF ( PI3*SSILA(3) .LT. 0.D0 ) GOTO 9001
      IF ( PA3*SSILB(3) .LT. 0.D0 ) GOTO 9001
      SSILA(1)=POX
      SSILA(2)=POY
      SSILA(3)=PI3
      SSILA(4)=PI4
      SSILB(1)=PUX
      SSILB(2)=PUY
      SSILB(3)=PA3
      SSILB(4)=PA4

      DO 13 I=1,4
        SILA(I)=SSILA(I)
        SILB(I)=SSILB(I)
        SIA(I)=SSIA(I)
        SIB(I)=SSIB(I)
13    CONTINUE

      IF ( ISH .LT. 90 ) GOTO 9000

      IF ( ABS(PIP*PIM-POT**2) .GT. 1.D-4 ) THEN
        CALL UTMSG('HASTFR')
        WRITE(IFCH,*)'*****  PIP*PIM /= POT**2'
        WRITE(IFCH,*)'PIP*PIM=',PIP*PIM
        WRITE(IFCH,*)'POT**2=',POT**2
        WRITE(IFCH,*)'PIP=',PIP
        WRITE(IFCH,*)'PIM=',PIM
        WRITE(IFCH,*)'POT=',POT
        CALL UTMSGF
      ENDIF
      IF ( ABS(PAP*PAM-PUT**2) .GT. 1.D-4 ) THEN
        CALL UTMSG('HASTFR')
        WRITE(IFCH,*)'*****  PAP*PAM /= PUT**2'
        WRITE(IFCH,*)'PAP*PAM=',PAP*PAM
        WRITE(IFCH,*)'PUT**2=',PUT**2
        WRITE(IFCH,*)'PAP=',PAP
        WRITE(IFCH,*)'PAM=',PAM
        WRITE(IFCH,*)'PUT=',PUT
        CALL UTMSGF
      ENDIF
      IF ( ABS(SSILA(4)**2
     *        -SSILA(1)**2-SSILA(2)**2-SSILA(3)**2) .GT. 1.D-4 ) THEN
        CALL UTMSG('HASTFR')
        WRITE(IFCH,*)'*****  MASS**2 OF SSILA NONZERO'
        WRITE(IFCH,*)'MASS**2=',SSILA(4)**2
     *                         -SSILA(1)**2-SSILA(2)**2-SSILA(3)**2
        CALL UTMSGF
      ENDIF
      IF ( ABS(SSILB(4)**2
     *         -SSILB(1)**2-SSILB(2)**2-SSILB(3)**2) .GT. 1.D-4 ) THEN
        CALL UTMSG('HASTFR')
        WRITE(IFCH,*)'*****  MASS**2 OF SSILB NONZERO'
        WRITE(IFCH,*)'MASS**2=',SSILB(4)**2
     *                         -SSILB(1)**2-SSILB(2)**2-SSILB(3)**2
        CALL UTMSGF
      ENDIF
      DO 14 N=1,4
        IF ( ABS(PLAB(N)
     *         -SSILA(N)-SSILB(N)-SSIA(N)-SSIB(N)) .GT. 1.D-4 ) GOTO 15
14    CONTINUE
      GOTO 16
15    CONTINUE
      CALL UTMSG('HASTFR')
      WRITE(IFCH,*)'***** SSILA+SSILB /= SSILA_NEW+SSILB_NEW+SSIA+SSIB'
      WRITE(IFCH,*)'SSILA_OLD+SSILB_OLD:'
      WRITE(IFCH,*)PLAB
      WRITE(IFCH,*)'SSILA_NEW+SSILB_NEW+SSIA+SSIB:'
      WRITE(IFCH,*)((SSILA(K)+SSILB(K)+SSIA(K)+SSIB(K)),K=1,4)
      WRITE(IFCH,*)'SSILA_NEW:'
      WRITE(IFCH,*)(SSILA(N),N=1,4)
      WRITE(IFCH,*)'SSILB_NEW:'
      WRITE(IFCH,*)(SSILB(N),N=1,4)
      WRITE(IFCH,*)'SSIA:'
      WRITE(IFCH,*)(SSIA(N),N=1,4)
      WRITE(IFCH,*)'SSIB:'
      WRITE(IFCH,*)(SSIB(N),N=1,4)
      CALL UTMSGF
16    CONTINUE

9000  RETURN

9001  IRET=1
      IF ( ISH .LT. 90 ) RETURN
      CALL UTMSG('HASTFR')
      WRITE(IFCH,*)'STRING SUBTRACTION NOT POSSIBLE'
      WRITE(IFCH,*)'SNGL(SSILA):'
      WRITE(IFCH,*)(SNGL(SSILA(N)),N=1,4)
      WRITE(IFCH,*)'SNGL(SSILB):'
      WRITE(IFCH,*)(SNGL(SSILB(N)),N=1,4)
      WRITE(IFCH,*)'SNGL(SSIA):'
      WRITE(IFCH,*)(SNGL(SSIA(N)),N=1,4)
      WRITE(IFCH,*)'SNGL(SSIB):'
      WRITE(IFCH,*)(SNGL(SSIB(N)),N=1,4)
      WRITE(IFCH,*)'PEP,PEM:',PEP,PEM
      WRITE(IFCH,*)'POT,PUT:',POT,PUT
      WRITE(IFCH,*)'A,D:',A,D
      WRITE(IFCH,*)'A**2-D**2:',A**2-D**2
      WRITE(IFCH,*)'PYP,PYM:',PYP,PYM
      WRITE(IFCH,*)'PIP,PAM:',PIP,PAM
      WRITE(IFCH,*)'PI3,SILA(3):',PI3,SILA(3)
      WRITE(IFCH,*)'PA3,SILB(3):',PA3,SILB(3)
      CALL UTMSGF
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTFS(SILA,SIXA,IFLA,XFLA,PTA,PTAH
     *                 ,SILB,SIXB,IFLB,XFLB,PTB,PTBH  ,STR,IRET  ,JORD)

C-----------------------------------------------------------------------
C  FORMS A SHORT (=Q-QBAR) STRING
C  IRET=0: OK
C  IRET=1: IN SR HASI: REMNANT CHANGES DIRECTION OR ZERO SI().
C  IRET=2: IN SR HASI: JC(,)=10 .
C  IRET=3: STRING MASS TOO SMALL. NOT 4 , 5 .
C  IRET=4: STRING MASS TOO SMALL. EQUAL FLAVOUR,ZERO MOMENTUM,XFL='000'.
C  IRET=5: STRING MASS TOO SMALL. VALENCE QUARKS INVOLVED.
C-----------------------------------------------------------------------
      PARAMETER (NPTQ=129)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTQ/    QPTH(NPTQ),QPTQ(NPTQ),XPTQ(NPTQ),QPTQMX,QPTHMX
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      PTDIFF(2),SIA(NSI),SIB(NSI),SILA(NSI),SILB(NSI)
     *         ,SIXA(NSI,NSIX),SIXB(NSI,NSIX),STR(NSI,NSIX+1),STS(NSI,2)
      CHARACTER XFLA*3,XFLB*3
      SAVE
C-----------------------------------------------------------------------
      DO 50 N=1,NSI
        STR(N,1)=0.
        STR(N,2)=0.
        STR(N,3)=0.
50    CONTINUE

      CALL UTKSIX(SIXA,KMAXA)
      CALL UTKSIX(SIXB,KMAXB)

      R=RANGEN()
      IF     ( IOPTQ .EQ. 2 ) THEN
        PT = SQRT( -4.*PTQ**2/PI * LOG(1.-QPTQMX*R) )
      ELSEIF ( IOPTQ .EQ. 3 ) THEN
        PT = PTQ*SQRT( QPTQMX*R/(1.-QPTQMX*R) )
      ELSE
        PT=UTINVT(NPTQ,XPTQ,QPTQ,R*QPTQ(NPTQ))
      ENDIF
      PHI=2.*PI*RANGEN()
      SIA(1)=PT*COS(PHI)
      SIB(1)=-SIA(1)
      SIA(2)=PT*SIN(PHI)
      SIB(2)=-SIA(2)

      IXFLAB=0
      IF (  XFLA.NE.'---' .AND. XFLB.NE.'---' .AND.
     *       (XFLA .NE. '000'  .OR.  XFLB .NE. '000') ) IXFLAB=1
      CALL HASI('QUARK',SILA,IFLA,XFLA,PTA,PTAH,SIA,IRET,JORD,IXFLAB
     *          ,PTDIFF)
      IF ( IRET .EQ. 1 ) GOTO 9001
      IF ( IRET .EQ. 2 ) GOTO 9002
      SILA(1)=SILA(1)-PTDIFF(1)
      SILA(2)=SILA(2)-PTDIFF(2)
      SILB(1)=SILB(1)+PTDIFF(1)
      SILB(2)=SILB(2)+PTDIFF(2)
      CALL HASI('ANTIQ',SILB,IFLB,XFLB,PTB,PTBH,SIB,IRET,JORD,IXFLAB
     *          ,PTDIFF)
      IF ( IRET .EQ. 1 ) GOTO 9001
      IF ( IRET .EQ. 2 ) GOTO 9002
      SILB(1)=SILB(1)-PTDIFF(1)
      SILB(2)=SILB(2)-PTDIFF(2)
      SILA(1)=SILA(1)+PTDIFF(1)
      SILA(2)=SILA(2)+PTDIFF(2)
      IF ( SILA(5) .EQ. 0.  .AND.  SILA(6) .EQ. 0. ) THEN
        SILA(5)=100000.
        SILA(6)=100000.
      ENDIF
      IF ( SILB(5) .EQ. 0.  .AND.  SILB(6) .EQ. 0. ) THEN
        SILB(5)=100000.
        SILB(6)=100000.
      ENDIF
      IF     ( XFLA .EQ. '000'  .AND.  XFLB .EQ. '000' ) THEN
        GOTO 9003
      ELSEIF ( XFLA .NE. '---'  .AND.  XFLB .NE. '---' ) THEN
        CALL HASTFR(SILA,SILB,SIA,SIB,IRET)
        IF ( IRET .EQ. 1 ) GOTO 9003
        DO 65 N=1,NSI
          STS(N,1)=SIA(N)
          STS(N,2)=SIB(N)
65      CONTINUE
        CALL UTAMST(STS,AM,AMIN,IRET)
        IF ( IRET .NE. 0 ) GOTO 9003
        DO 70 N=1,NSI
          STR(N,1)=SIA(N)
          STR(N,2)=SIB(N)
70      CONTINUE
      ELSEIF ( XFLA .NE. '---'  .AND.  XFLB .EQ. '---' ) THEN
        IF ( KMAXB+1 .GT. NSIX ) THEN
          CALL UTSTOP('HASTFS: NSIX TOO SMALL                  ')
        ENDIF
        DO 60 N=1,NSI
          SIXB(N,KMAXB+1)=SIA(N)
          IF ( KMAXB+2 .LE. NSIX ) SIXB(N,KMAXB+2)=0.
60      CONTINUE
        IF ( IFLA .LT. 0 ) THEN
          SIXB(5,KMAXB+1)=-SIXB(5,KMAXB+1)
          SIXB(6,KMAXB+1)=-SIXB(6,KMAXB+1)
        ENDIF
        KMAXB=KMAXB+1
      ELSEIF ( XFLA .EQ. '---'  .AND.  XFLB .NE. '---' ) THEN
        IF ( KMAXA+1 .GT. NSIX ) THEN
          CALL UTSTOP('HASTFS: NSIX TOO SMALL                  ')
        ENDIF
        DO 80 N=1,NSI
          SIXA(N,KMAXA+1)=SIB(N)
          IF ( KMAXA+2 .LE. NSIX ) SIXA(N,KMAXA+2)=0.
80      CONTINUE
        IF ( IFLB .LT. 0 ) THEN
          SIXA(5,KMAXA+1)=-SIXA(5,KMAXA+1)
          SIXA(6,KMAXA+1)=-SIXA(6,KMAXA+1)
        ENDIF
        KMAXA=KMAXA+1
      ELSEIF ( XFLA .EQ. '---'  .AND.  XFLB .EQ. '---' ) THEN
C  NO ACTION
      ELSE
        CALL UTSTOP('HASTFS: IF/ELSE ERROR                   ')
      ENDIF

ctp060203 9000  IRET=0
      IRET=0
      IF ( ISH .GE. 91 ) WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
102   FORMAT(' SIA: ',4F13.5,2F8.0,2X,A3,
     *     /,' SIB: ',4F13.5,2F8.0,2X,A3)
      RETURN

9001  IRET=1
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
        WRITE(IFCH,*)'REMNANT CHANGES DIRECTION'
        WRITE(IFCH,*)' OR ZERO SI(1/2/3)'
        WRITE(IFCH,*)' '
      ENDIF
      RETURN

9002  IRET=2
C     JC(,)=10
      IF ( ISH .GE. 91 ) WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
      RETURN

9003  IRET=3
      IF ( XFLA(1:2) .EQ. 'VA'  .OR.  XFLB(1:2) .EQ. 'VA' ) GOTO 9005
      WT=ABS(SIA(1)+SIB(1))+ABS(SIA(2)+SIB(2))
     *                              +ABS(SIA(3)+SIB(3))
      IF ( ABS(IFLA).EQ.ABS(IFLB) .AND. WT.LT.1.E-3
     *            .AND. XFLA.EQ.'000' .AND. XFLB.EQ.'000' ) GOTO 9004
      IFLB=SIGN( IABS(IFLA), IFLB)
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
        WRITE(IFCH,*)'STRING MASS TOO SMALL (3)'
        WRITE(IFCH,*)' '
      ENDIF
      XFLA='000'
      XFLB='000'
      RETURN

9004  IRET=4
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
        WRITE(IFCH,*)'STRING MASS TOO SMALL (4)'
        WRITE(IFCH,*)' '
      ENDIF
      RETURN

9005  IRET=5
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,102)SIA,XFLA,SIB,XFLB
        WRITE(IFCH,*)'STRING MASS TOO SMALL (5)'
        WRITE(IFCH,*)' '
      ENDIF
      XFLA(1:2)='SE'
      XFLB(1:2)='SE'
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTFW(SIL,SIX)

C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL SIL(NSI),SIX(NSI,NSIX)
      SAVE
C-----------------------------------------------------------------------
      CALL UTKSIX(SIX,KMAX)
      WRITE(IFCH,103)SIL,((SIX(I,J),I=1,NSI),J=1,KMAX)
103   FORMAT(/,' SIL: ',4F13.5,2F8.0,/,' SIX: ',4F13.5,2F8.0,/
     *,50('      ',4F13.5,2F8.0,/))
      RETURN
      END
C=======================================================================

      SUBROUTINE HASTPR(STRO,ISPLT)

C-----------------------------------------------------------------------
C  PROCESSES A STRING
C  OUT: ISPLT=0 : PROCESSING OF STRO FINISHED
C            >0 : HASTPR TO BE REDONE TO PROCESS REDUCED STRING
C            -1 : ERROR
C            -3 : HASTPR TO BE REDONE WITH NEW EPART,APART
C            -4 : HAHA COLLISION TO BE REDONE BEC OF TOO FAST LD PTL
C                             (ONLY FOR ICHOIC=2)
C            -5 : HASTPR TO BE REDONE BECAUSE OF PRODUCED S=3/2 PARTICLE
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CJSPLI/  ALEAD,APART,ELEAD,EPART,SGNSIL,JPART,NSCC,NSCCX
      COMMON /CKOL/    KOL
      COMMON /CLEAD/   COOAV3,COOAV4,LEAD
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CNSUC/   NSUC
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CPZSTR/  ESTRL,PZSTRL,ISEA,ISTRL
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CSTSH/   NSTSH
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION DAUXIL,PPE,PPM,PPP,PPX,PPY,PPZ,PP1(4),PP2(4)
      REAL             P1(4),P2(4),STR(NSI,2),STRO(NSI,NSIX+1)
      INTEGER          IC(2),IC1(2),IC1X(2),IC2(2)
     *                ,JC(NFLAV,2),JC1(NFLAV,2),JC2(NFLAV,2)
     *                ,JC3(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IF     ( ISPLT .EQ. -3 ) THEN
        JPART=1
        ISPLT=0
      ELSEIF ( ISPLT .EQ.  0 ) THEN
        JPART=0
      ENDIF
      ISPLT0=ISPLT

      CALL UTKSTR(STRO,KMAX)

C  ZERO STRING
C  -----------
      IF ( KMAX .EQ. 0 ) THEN
        IF ( ISH .GE. 91 ) THEN
          WRITE(IFCH,*)'ZERO STRING'
          WRITE(IFCH,*)' '
        ENDIF
        ISPLT=0
        RETURN
      ENDIF

C  PRINT
C  -----
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,105)(STRO(I,1),I=1,4),(NINT(STRO(I,1)),I=5,6)
105     FORMAT(' STR: ',4F13.5,2I8)
        IF ( KMAX .GT. 1 ) THEN
          DO 8 K=2,KMAX
            WRITE(IFCH,104)(STRO(I,K),I=1,4),(NINT(STRO(I,K)),I=5,6)
104         FORMAT('      ',4F13.5,2I8)
 8        CONTINUE
        ENDIF
        WRITE(IFCH,*)' '
      ENDIF

C  CHECK LEADING OBJECT
C  --------------------
      IC1(1)=NINT(STRO(4+1,1))
      IC1(2)=NINT(STRO(4+2,1))
      CALL IDCOMP(IC1,IC1X,JC3,1)
      STRO(4+1,1)=IC1X(1)
      STRO(4+2,1)=IC1X(2)
      CALL IDDECO(IC1X,JC1)
      NPA=0
      DO 7 N=1,NFLAV
        NPA=NPA+JC1(N,1)+JC1(N,2)
 7    CONTINUE

C  SPLIT STRING
C  ------------
      IF ( KMAX .GT. 2  .OR. KMAX.EQ.2 .AND. ISPLT.GT.0
     *                  .OR. KMAX.EQ.2 .AND. NPA.GT.3
     *                  .OR. KMAX.EQ.2 .AND. LEAD.EQ.1 ) THEN
        IF ( ISPLT .EQ. 0 ) THEN
          NSUC=0
          KMAXOR=KMAX
          IF ( STRO(3,1) .LT. 0. ) THEN
            SGNSIL=-1
          ELSE
            SGNSIL=1
          ENDIF
          DO 17 N=1,NFLAV
            JC2(N,1)=0
            JC2(N,2)=0
17        CONTINUE
          DO 15 K=1,KMAX
            IC1(1)=NINT(ABS(STRO(4+1,K)))
            IC1(2)=NINT(ABS(STRO(4+2,K)))
            CALL IDDECO(IC1,JC1)
            DO 16 N=1,NFLAV
              JC2(N,1)=JC2(N,1)+JC1(N,1)
              JC2(N,2)=JC2(N,2)+JC1(N,2)
              IF ( N .GT. 4  .AND.  (JC2(N,1).NE.0 .OR.
     *                               JC2(N,2) .NE. 0) ) THEN
                 CALL UTSTOP('HASTPR: FLAVOUR > 4                     ')
               ENDIF
16          CONTINUE
15        CONTINUE
          KEU=JC2(1,1)-JC2(1,2)
          KED=JC2(2,1)-JC2(2,2)
          KES=JC2(3,1)-JC2(3,2)
          KEC=JC2(4,1)-JC2(4,2)
          ALEAD=UTAMNU(KEU,KED,KES,KEC,0)
C-C       ALEAD=0.
          ELEAD=STRO(4,1)
          IF ( JPART .EQ. 0 ) THEN
            NSCCX=KMAXOR-1
            NSTR0=NSTR
          ELSE
            NSCCX=MAX(1,NSCC)
            NSTR=NSTR0
          ENDIF
          APART=ALEAD/NSCCX
          EPART=ELEAD/NSCCX
          IF ( ALEAD .GT. ELEAD ) THEN
            IF ( ISH .GE. 90 ) THEN
              CALL UTMSG('HASTPR')
              WRITE(IFCH,*)'*****  ALEAD > ELEAD  ==> RET CODE = -1'
              WRITE(IFCH,*)'ALEAD=',ALEAD,'   ELEAD=',ELEAD
              CALL UTMSGF
            ENDIF
            ISPLT=-1
            RETURN
          ENDIF
        ENDIF
        CALL JSPLIT(STRO,STR,KOLSP,IER,KMAXOR)
        IF     ( IER .EQ. 1 ) THEN
          ISPLT=-1
          RETURN
        ELSEIF ( IER .EQ. 3 ) THEN
          ISPLT=-3
          RETURN
        ENDIF
        ISPLT=ISPLT+1
        IF ( IER .EQ. 2 ) THEN
          IF ( ISH .GE. 91 ) THEN
            WRITE(IFCH,*)'MULTISTRING: ABSORPTION OF ONE LEG'
            WRITE(IFCH,*)' '
          ENDIF
          RETURN
        ENDIF
        KOLZ=KOLSP
        IF ( KOLZ .LT. 1  .OR.  KOLZ .GT. KOLLMX ) THEN
          WRITE(IFCH,*)'KOLZ=',KOLZ
          CALL UTSTOP('HASTPR: KOLZ OUT OF RANGE (1)           ')
        ENDIF
        KIN=1
      ELSE
        IF ( ISPLT .GT. 0 ) THEN
          KIN=1
        ELSE
          KIN=0
        ENDIF
        DO 5 N=1,NSI
          IF ( N .LE. 4 ) THEN
            STR(N,1)=STRO(N,1)
            STR(N,2)=STRO(N,2)
          ELSE
            STR(N,1)=ABS(STRO(N,1))
            STR(N,2)=ABS(STRO(N,2))
          ENDIF
 5      CONTINUE
        IF     ( ISPLT .GT. 0 ) THEN
          IF ( SGNSIL .LT. 0. ) THEN
            KOLZ=KTARGA(2,ITARG)
          ELSE
            KOLZ=KPROJA(2,IPROJ)
          ENDIF
        ELSEIF ( LEAD .EQ. 1  .AND.  KMAX .EQ. 2 ) THEN
          IF ( STR(3,1) .LT. 0. ) THEN
            KOLZ=KTARGA(2,ITARG)
          ELSE
            KOLZ=KPROJA(2,IPROJ)
          ENDIF
        ELSE
          KOLZ=KOL
        ENDIF
        IF ( KOLZ .LT. 1  .OR.  KOLZ .GT. KOLLMX ) THEN
          IF(ISH.GE.90)THEN
            CALL UTMSG('HASTPR')
            WRITE(IFCH,*)'*****  KOLZ OUT OF RANGE;   KOLZ:',KOLZ
     *        ,'   SET TO:',KOL
            WRITE(IFCH,105)(STR(I,1),I=1,4),(NINT(STR(I,1)),I=5,6)
            WRITE(IFCH,104)(STR(I,2),I=1,4),(NINT(STR(I,2)),I=5,6)
            CALL UTMSGF
          ENDIF
          KOLZ=KOL
        ENDIF
        ISPLT=0
      ENDIF

C  ADD TWO JETS
C  ------------
      IC1(1)=NINT(STR(4+1,1))
      IC1(2)=NINT(STR(4+2,1))
      IC2(1)=NINT(STR(4+1,2))
      IC2(2)=NINT(STR(4+2,2))
      IA1=IC1(1)+IC1(2)
      IA2=IC2(1)+IC2(2)
      IF ( IA1 .EQ. 0  .AND.  IA2 .EQ. 0 ) THEN
        WRITE(IFCH,*)'STRO:'
        WRITE(IFCH,105)(STRO(I,1),I=1,4),(NINT(STRO(I,1)),I=5,6)
        IF ( KMAX .GT. 1 ) THEN
          DO 18 K=2,KMAX
            WRITE(IFCH,104)(STRO(I,K),I=1,4),(NINT(STRO(I,K)),I=5,6)
18        CONTINUE
        ENDIF
        WRITE(IFCH,*)'STR:'
        WRITE(IFCH,105)(STR(I,1),I=1,4),(NINT(STR(I,1)),I=5,6)
        WRITE(IFCH,104)(STR(I,2),I=1,4),(NINT(STR(I,2)),I=5,6)
        CALL UTSTOP('HASTPR: STR=0                           ')
      ENDIF
      CALL IDDECO(IC1,JC1)
      CALL IDDECO(IC2,JC2)
      N1=0
      N2=0
      DO 9 N=1,NFLAV
        N1=N1+JC1(N,1)+JC1(N,2)
        N2=N2+JC2(N,1)+JC2(N,2)
        JC(N,1)=JC1(N,1)+JC2(N,1)
        JC(N,2)=JC1(N,2)+JC2(N,2)
        IF ( N .GT. 4 .AND. (JC(N,1) .NE. 0 .OR.
     *                       JC(N,2) .NE. 0) ) THEN
          CALL UTSTOP ('HASTPR: FLAVOUR > 4                     ')
        ENDIF
 9    CONTINUE
      KEU=JC(1,1)-JC(1,2)
      KED=JC(2,1)-JC(2,2)
      KES=JC(3,1)-JC(3,2)
      KEC=JC(4,1)-JC(4,2)
      AMSTR=UTAMNU(KEU,KED,KES,KEC,6)
      CALL IDENCO(JC,IC,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('HASTPR: IDENCO RET CODE = 1             ')
      ENDIF
      IDSTR=IDTRA(IC,0,0,3)

C  DETERMINE ROT
C  -------------
      PP1(1)=STR(1,1)
      PP1(2)=STR(2,1)
      PP1(3)=STR(3,1)
      PP1(4)=STR(4,1)
      PP2(1)=STR(1,2)
      PP2(2)=STR(2,2)
      PP2(3)=STR(3,2)
      PP2(4)=STR(4,2)
      PPX=PP1(1)+PP2(1)
      PPY=PP1(2)+PP2(2)
      PPZ=PP1(3)+PP2(3)
      PPP=SQRT(PPX**2+PPY**2+PPZ**2)
      IF ( IA1 .NE. 0  .AND.  IA2 .NE. 0 ) THEN
        PP1(4)=SQRT(PP1(1)**2+PP1(2)**2+PP1(3)**2)
        PP2(4)=SQRT(PP2(1)**2+PP2(2)**2+PP2(3)**2)
        PPE=PP1(4)+PP2(4)
        PPM=SQRT((PPE-PPP)*(PPE+PPP))
        CALL UTLOB2(1,PPX,PPY,PPZ,PPE,PPM,PP1(1),PP1(2),PP1(3),PP1(4))
        CALL UTLOB2(1,PPX,PPY,PPZ,PPE,PPM,PP2(1),PP2(2),PP2(3),PP2(4))
      ELSE
        PPE=PP1(4)+PP2(4)
        DAUXIL=(PPE-PPP)*(PPE+PPP)
        IF ( DAUXIL .GT. 0.D0 ) THEN
          PPM=SQRT(DAUXIL)
        ELSE
          PPM=0.D0
          PPE=PPP
        ENDIF
        PP1(1)=0.D0
        PP1(2)=0.D0
        PP1(3)=0.D0
        PP1(4)=0.D0
        PP2(1)=0.D0
        PP2(2)=0.D0
        PP2(3)=0.D0
        PP2(4)=0.D0
      ENDIF
      PX=PPX
      PY=PPY
      PZ=PPZ
      E=PPE
      AM=PPM
      P=PPP
      P1(1)=PP1(1)
      P1(2)=PP1(2)
      P1(3)=PP1(3)
      P1(4)=PP1(4)
      P2(1)=PP2(1)
      P2(2)=PP2(2)
      P2(3)=PP2(3)
      P2(4)=PP2(4)
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,107)(P1(I),I=1,4),IC1
        WRITE(IFCH,107)(P2(I),I=1,4),IC2
107     FORMAT(' JET: ',4F13.5,2I8)
        WRITE(IFCH,*)' '
      ENDIF
CDH   IF ( P1(3) .NE. 0. ) THEN
CDH     DIRCN=SIGN(1.,P1(3))
CDH   ELSE
CDH     DIRCN=1.
CDH   ENDIF

C  MASS ADJUSTMENT
C  ---------------
      AM0=AM
      CALL IDRES(IDSTR,AM,IDSTRR,IADJ)

C  DIFFRACTIVE PARTICLE:
C-C   IF ( (IA1.EQ.0.OR.IA2.EQ.0) .AND. ISPLT0.EQ.0 ) THEN
C-C     IF ( MOD(IDSTRR,10).EQ.0 .AND. RANGEN().LT.0.1 ) THEN
C-C       AM=0.
C-C       CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
C-C       AM=AM+.300
C-C       CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
C-C     ENDIF
C-C   ENDIF

C  REMNANT:
      IF ( (IA1.EQ.0 .OR. IA2.EQ.0) .AND. ISPLT0.GT.0 ) THEN
        IF     ( NSUC .GT. 0 ) THEN
          AM=0.
          CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
          IF     ( N1+N2 .EQ. 2 ) THEN
            IF ( RANGEN() .LT. 0.33 ) AM=AM+0.400
          ELSEIF ( MOD(IDSTRR,10) .EQ. 0 ) THEN
            IF ( RANGEN() .LT. 0.33 ) AM=AM+0.300
          ELSE
            IF ( RANGEN() .GT. 0.33 ) THEN
              ISPLT=-5
              RETURN
            ENDIF
          ENDIF
          CALL IDRES(IDSTR,AM,IDSTRR,IADJ)

C  ALL LEGS ABSORBED
        ELSEIF ( NSUC .EQ. 0 ) THEN
          IF ( ISEA .EQ. 0 ) THEN
            ISPLT=-4
            RETURN
          ENDIF
          AM=0.
          IF ( ISTRL .EQ. 1 ) THEN
            AM=SQRT(ENGY**2*(1-ABS(PZSTRL)/PNLLX)+.94**2 )
          ELSE
            CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
            IF ( MOD(IDSTRR,10) .EQ. 0 ) THEN
              AM=AM+RANGEN()
            ELSE
              AM=AM+RANGEN()-0.30
            ENDIF
C-C         IF ( MOD(IDSTRR,10).EQ.0.AND.RANGEN().LT.0.33 ) AM=AM+0.30
          ENDIF
          CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
        ENDIF
      ENDIF

C  LEADING STRING = HADRON:
      IF ( LEAD.EQ.1 .AND. IA1.NE.0 .AND.IA2.NE.0
     *     .AND. ISPLT0.EQ.0 .AND. ISPLT.EQ.0 .AND. IDSTRR.NE.0 ) THEN
        AM=0.
        CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
        IF ( MOD(IDSTRR,10).EQ.0  .AND.  RANGEN().LT.0.33 ) AM=AM+0.300
        CALL IDRES(IDSTR,AM,IDSTRR,IADJ)
      ENDIF

      IDSTR=IDSTRR
      IF ( IDSTR .EQ. 0  .AND.  AM .LT. AMSTR ) AM=AMSTR
      PPM=AM
      E=SQRT(PPP**2+PPM**2)

C  WRITE /CSTR/
C  ------------
      NSTR=NSTR+1
      IF ( NSTR .GT. MXSTR ) THEN
        CALL UTSTOP('HASTPR: NSTR>MXSTR                      ')
      ENDIF
      IF ( LEAD .EQ. 0 ) NSTSH=1
      NSUC=NSUC+1
      IRLSTR(NSTR)=NSUC
      IF ( ISPLT0 .EQ. 0  .AND.  ISPLT .EQ. 0 ) IRLSTR(NSTR)=0
      ICSTR(1,NSTR)=IC1(1)
      ICSTR(2,NSTR)=IC1(2)
      ICSTR(3,NSTR)=IC2(1)
      ICSTR(4,NSTR)=IC2(2)
      PSTR(1,NSTR)=PX
      PSTR(2,NSTR)=PY
      PSTR(3,NSTR)=PZ
      PSTR(4,NSTR)=E
      PSTR(5,NSTR)=AM
      ROTSTR(1,NSTR)=P1(1)
      ROTSTR(2,NSTR)=P1(2)
      ROTSTR(3,NSTR)=P1(3)
      IF ( P1(1) .EQ. 0.  .AND.  P1(2) .EQ. 0.  .AND.  P1(3) .EQ. 0. )
     *                                                ROTSTR(3,NSTR)=1.
      XORSTR(1,NSTR)=COORD(1,KOLZ)
      XORSTR(2,NSTR)=COORD(2,KOLZ)
      XORSTR(3,NSTR)=COORD(3,KOLZ)
      XORSTR(4,NSTR)=COORD(4,KOLZ)
      IORSTR(NSTR)=-KOLZ
      AMSAC=AMSAC+AM
      IF ( ISH .GE. 91 ) THEN
        IF     ( ISPLT .GT. 0 ) THEN
          WRITE(IFCH,*)'SPLIT OFF STRING:'
        ELSEIF ( ISPLT0 .GT. 0 ) THEN
          WRITE(IFCH,*)'REMAINDER:'
        ELSEIF ( LEAD .EQ. 1 ) THEN
          WRITE(IFCH,*)'ORDINARY BARYONIC STRING:'
        ELSE
          WRITE(IFCH,*)'ORDINARY MESONIC STRING:'
        ENDIF
        WRITE(IFCH,106)NSTR,(ICSTR(K,NSTR)/100,K=1,4)
     *                      ,(PSTR(I,NSTR),I=3,5)
106     FORMAT(/,' /CSTR/',I4,3X,4I5,2X,3(E10.3),/)
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE HDECMP(BAR,SIL,SIX)

C-----------------------------------------------------------------------
C  DECOMPOSES BAR INTO SIL,SIX
C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      REAL BAR(NSI,NHA),SIL(NSI),SIX(NSI,NSIX)
      SAVE
C-----------------------------------------------------------------------
      DO 1 N=1,NSI
        SIL(N)=BAR(N,2)
 1    CONTINUE
      DO 3 M=1,NSIX
        SIXSQR=0.
        DO 2 N=1,NSI
          SIX(N,M)=BAR(N,2+M)
          SIXSQR=SIXSQR+SIX(N,M)**2
 2      CONTINUE
        IF ( SIXSQR .LE. 1.E-5 ) RETURN
 3    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE HRESCL(J1,J2,PSUM,IFAIL)

C-----------------------------------------------------------------------
C  RESCALES STRING MOMENTA OF STRINGS J1-J2 TO HAVE TOTAL MOM PSUM.
C-----------------------------------------------------------------------
      PARAMETER (MXSTR=3000)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION PADD(5),PP(5),PPSUM(5)
      REAL             PSUM(5)
      DATA             ERRLIM /.001/
      SAVE
C-----------------------------------------------------------------------
      IFAIL=1

      PPSUM(1)=PSUM(1)
      PPSUM(2)=PSUM(2)
      PPSUM(3)=PSUM(3)
      PPSUM(4)=PSUM(4)
      PPSUM(5)=PSUM(5)

      IF ( J1 .GE. J2 ) THEN
        CALL UTSTOP('HRESCL: J1 .GE. J2                      ')
      ENDIF

      PADD(1)=0.D0
      PADD(2)=0.D0
      PADD(3)=0.D0
      PADD(4)=0.D0
      PADD(5)=0.D0
      DO 110 J=J1,J2
        PADD(1)=PADD(1)+PSTR(1,J)
        PADD(2)=PADD(2)+PSTR(2,J)
        PADD(3)=PADD(3)+PSTR(3,J)
        PADD(4)=PADD(4)+PSTR(4,J)
        PADD(5)=PADD(5)+PSTR(5,J)
110   CONTINUE
      IF ( PADD(5) .GE. PPSUM(5) ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('HRESCL')
          WRITE(IFCH,*)'*****  SUM OF STR MASSES .GE. PSUM(5)'
          DO 1 J=J1,J2
            WRITE(IFCH,109)J,(ICSTR(K,J)/100,K=1,4)
     *                       ,(PSTR(L,J),L=3,5)
109         FORMAT(' /CSTR/',I4,3X,4I5,2X,3(E10.3))
 1        CONTINUE
          WRITE(IFCH,*)'PPSUM(345):',(SNGL(PPSUM(K)),K=3,5)
          CALL UTMSGF
        ENDIF
        RETURN
      ENDIF
      PADD(5)=PADD(4)**2-PADD(1)**2-PADD(2)**2-PADD(3)**2
      IF ( PADD(5) .LE. 0.D0 ) THEN
        DO 2 J=J1,J2
          WRITE(IFCH,108)J,(PSTR(L,J),L=1,5)
108       FORMAT(' /CSTR/',I4,2X,5(E10.3))
 2      CONTINUE
        CALL UTSTOP('HRESCL: MASS**2 OF STRING-SUM NEGATIVE  ')
      ENDIF
      PADD(5)=SQRT(PADD(5))

C  BOOST STRINGS TO REST
C  ---------------------
      DO 115 J=J1,J2
        PP(1)=PSTR(1,J)
        PP(2)=PSTR(2,J)
        PP(3)=PSTR(3,J)
        PP(4)=PSTR(4,J)
        CALL UTLOB2(1,PADD(1),PADD(2),PADD(3),PADD(4),PADD(5)
     *                ,PP(1),PP(2),PP(3),PP(4))
        PSTR(1,J)=PP(1)
        PSTR(2,J)=PP(2)
        PSTR(3,J)=PP(3)
        PSTR(4,J)=PP(4)
115   CONTINUE

C  RESCALE MOMENTA IN REST FRAME
C  -----------------------------
      SCAL=1.
      DO 301 IPASS=1,200
        SUM=0.
        DO 310 J=J1,J2
          PSTR(1,J)=SCAL*PSTR(1,J)
          PSTR(2,J)=SCAL*PSTR(2,J)
          PSTR(3,J)=SCAL*PSTR(3,J)
          PSTR(4,J)=SQRT(PSTR(1,J)**2+PSTR(2,J)**2+PSTR(3,J)**2
     *                  +PSTR(5,J)**2)
          SUM=SUM+PSTR(4,J)
310     CONTINUE
        SCAL=PSUM(5)/SUM
        IF ( ABS(SCAL-1.) .LE. ERRLIM ) GOTO 300
301   CONTINUE
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('HRESCL')
        WRITE(IFCH,*)'*****  SCAL=',SCAL
        CALL UTMSGF
      ENDIF
300   CONTINUE

C  BOOST BACK WITH PPSUM
C  ---------------------
      DO 315 J=J1,J2
        PP(1)=PSTR(1,J)
        PP(2)=PSTR(2,J)
        PP(3)=PSTR(3,J)
        PP(4)=PSTR(4,J)
        CALL UTLOB2(-1,PPSUM(1),PPSUM(2),PPSUM(3),PPSUM(4),PPSUM(5)
     *             ,PP(1),PP(2),PP(3),PP(4))
        PSTR(1,J)=PP(1)
        PSTR(2,J)=PP(2)
        PSTR(3,J)=PP(3)
        PSTR(4,J)=PP(4)
315   CONTINUE

      IFAIL=0
      RETURN
      END
C=======================================================================

      SUBROUTINE IDCOMJ(JC)

C-----------------------------------------------------------------------
C  COMPACTIFIES JC
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER IC(2),ICX(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      CALL IDCOMP(IC,ICX,JC,2)
      RETURN
      END
C=======================================================================

      SUBROUTINE IDCOMK(IC)

C-----------------------------------------------------------------------
C  COMPACTIFIES IC
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER IC(2),ICX(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      CALL IDCOMP(IC,ICX,JC,1)
      IC(1)=ICX(1)
      IC(2)=ICX(2)
      RETURN
      END
C=======================================================================

      SUBROUTINE IDCOMP(IC,ICX,JC,IM)

C-----------------------------------------------------------------------
C  COMPACTIFIES IC,JC
C  INPUT: IM (1 OR 2)
C         IC (IF IM=1)
C         JC (IF IM=2)
C  OUTPUT: ICX (IF IM=1)
C          JC
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER IC(2),ICX(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IF ( IM .EQ. 1 ) CALL IDDECO(IC,JC)
      ICX(1)=0
      ICX(2)=0
      DO 2 N=1,NFLAV
        IF ( JC(N,1) .NE. 0  .OR. JC(N,2) .NE. 0 ) GOTO 1
 2    CONTINUE
      RETURN
 1    L=0
      DO 3 N=1,NFLAV
        K=MIN(JC(N,1),JC(N,2))
        JC(N,1)=JC(N,1)-K
        JC(N,2)=JC(N,2)-K
        IF ( JC(N,1) .LT. 0  .OR.  JC(N,2) .LT. 0 ) THEN
          CALL UTSTOP('IDCOMP: JC NEGATIVE                     ')
        ENDIF
        L=L+JC(N,1)+JC(N,2)
 3    CONTINUE
      IF ( L .EQ. 0 ) THEN
        JC(1,1)=1
        JC(1,2)=1
      ENDIF
      IF ( IM .EQ. 1 ) THEN
        CALL IDENCO(JC,ICX,IRETEN)
        IF ( IRETEN .EQ. 1 ) THEN
          CALL UTSTOP('IDCOMP: IDENCO RET CODE = 1             ')
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE IDDECO(IC,JC)

C-----------------------------------------------------------------------
C  DECODE PARTICLE ID
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER IC(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      ICI=IC(1)
      JC(6,1)=MOD(ICI,10)
      JC(5,1)=MOD(ICI/10,10)
      JC(4,1)=MOD(ICI/100,10)
      JC(3,1)=MOD(ICI/1000,10)
      JC(2,1)=MOD(ICI/10000,10)
      JC(1,1)=MOD(ICI/100000,10)
      ICI=IC(2)
      JC(6,2)=MOD(ICI,10)
      JC(5,2)=MOD(ICI/10,10)
      JC(4,2)=MOD(ICI/100,10)
      JC(3,2)=MOD(ICI/1000,10)
      JC(2,2)=MOD(ICI/10000,10)
      JC(1,2)=MOD(ICI/100000,10)
      RETURN
      END
C=======================================================================

      SUBROUTINE IDENCO(JC,IC,IRETEN)

C-----------------------------------------------------------------------
C  ENCODE PARTICLE ID
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER IC(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IRETEN=0
      IC(1)=0
      DO 20 I=1,NFLAV
        IF ( JC(I,1) .GE. 10 ) GOTO 22
        IC(1)=IC(1)+JC(I,1)*10**(NFLAV-I)
20    CONTINUE
      IC(2)=0
      DO 21 I=1,NFLAV
        IF ( JC(I,2) .GE. 10 ) GOTO 22
        IC(2)=IC(2)+JC(I,2)*10**(NFLAV-I)
21    CONTINUE
      RETURN
22    IRETEN=1
      IC(1)=0
      IC(2)=0
      RETURN
      END
C=======================================================================

      SUBROUTINE IDFLAV(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)

C-----------------------------------------------------------------------
C  UNPACKS THE IDENT CODE ID=+/-IJKL
C
C          MESONS--
C          I=0, J<=K, +/- IS SIGN FOR J
C          ID=110 FOR PI0, ID=220 FOR ETA, ETC.
C
C          BARYONS--
C          I<=J<=K IN GENERAL
C          J<I<K FOR SECOND STATE ANTISYMMETRIC IN (I,J), EG. L = 2130
C
C          OTHER--
C          ID=1,...,6 FOR QUARKS
C          ID=9 FOR GLUON
C          ID=10 FOR PHOTON
C          ID=11,...,16 FOR LEPTONS
C          ID=20 FOR KS, ID=-20 FOR KL
C
C          I=21...26 FOR SCALAR QUARKS
C          I=29 FOR GLUINO
C          I=30 FOR PHOTINO
C          I=31...36 FOR SCALAR LEPTONS
C          I=39 FOR WINO
C          I=40 FOR ZINO
C
C          ID=80 FOR W+
C          ID=81,...,89 FOR HIGGS MESONS
C          ID=90 FOR Z0
C
C          DIQUARKS--
C          ID=+/-IJ00, I<J FOR DIQUARK COMPOSED OF I,J.
C
C          INDEX IS A SEQUENCE NUMBER USED INTERNALLY
C-----------------------------------------------------------------------
      PARAMETER (NMES=2)
      PARAMETER (NQLEP=41)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      SAVE
C-----------------------------------------------------------------------
      IDABS=ABS(ID)
      I=IDABS/1000
      J1=IDABS-I*1000
      J = J1/100
      ISI = SIGN(1,ID)
      IF ( ID .NE. 0  .AND.  MOD(ID,100) .EQ. 0 ) GOTO 300
      IF ( J .EQ. 0 ) GOTO 200
      K1= J1 - J*100
      K = K1/10
      JSPIN= K1 - K*10
      IF ( I .EQ. 0 ) GOTO 100
C  BARYONS
C  ONLY X,Y BARYONS ARE QQX, QQY, Q=U,D,S.
      IFL1=ISI*I
      IFL2=ISI*J
      IFL3=ISI*K
      IF ( K .LE. 6 ) THEN
        INDEX=MAX(I-1,J-1)**2+I+MAX(I-J,0)+(K-1)*K*(2*K-1)/6
     *                        +109*JSPIN+36*NMES+NQLEP+11
      ELSE
        INDEX=MAX(I-1,J-1)**2+I+MAX(I-J,0)+9*(K-7)+91
     *                        +109*JSPIN+36*NMES+NQLEP+11
      ENDIF
      RETURN
C  MESONS
100   CONTINUE
      IFL1=0
      IFL2=ISI*J
      IFL3=ISI*K
      INDEX=J+K*(K-1)/2+36*JSPIN+NQLEP
      INDEX=INDEX+11
      RETURN
C  QUARKS, LEPTONS, ETC
200   CONTINUE
      IFL1=0
      IFL2=0
      IFL3=0
      JSPIN=0
      INDEX=IDABS
      IF ( IDABS .LT. 20 ) RETURN
C  DEFINE INDEX=20 FOR KS, INDEX=21 FOR KL
      INDEX=IDABS+1
      IF ( ID .EQ. 20 ) INDEX=20
C  INDEX=NQLEP+1,...,NQLEP+11 FOR W+, HIGGS, Z0
      IF ( IDABS .LT. 80 ) RETURN
      INDEX=NQLEP+IDABS-79
      RETURN
300   CONTINUE
      IFL1=ISI*I
      IFL2=ISI*J
      IFL3=0
      JSPIN=0
      INDEX=0
      RETURN
      END
C=======================================================================

      FUNCTION IDLABL(ID)

C-----------------------------------------------------------------------
C  RETURNS THE CHARACTER*8 LABEL FOR THE PARTICLE ID
C-----------------------------------------------------------------------
      PARAMETER (NMES=2)
      PARAMETER (NQLEP=41)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP

      CHARACTER*8 IDLABL,LABAR0(109),LABAR1(109),LAQQ(21),LBAR0(109)
     *           ,LBAR1(109),LLEP(104),LMES0(64),LMES1(64),LQQ(21)

C          DIQUARK LABELS
      DATA LQQ/
     1'UU0. ','UD0. ','DD0. ','US0. ','DS0. ','SS0. ','UC0. ','DC0. ',
     2'SC0. ','CC0. ','UB0. ','DB0. ','SB0. ','CB0. ','BB0. ','UT0. ',
     3'DT0. ','ST0. ','CT0. ','BT0. ','TT0. '/
      DATA LAQQ/
     1'AUU0.','AUD0.','ADD0.','AUS0.','ADS0.','ASS0.','AUC0.','ADC0.',
     2'ASC0.','ACC0.','AUB0.','ADB0.','ASB0.','ACB0.','ABB0.','AUT0.',
     3'ADT0.','AST0.','ACT0.','ABT0.','ATT0.'/
C          QUARK AND LEPTON LABELS
      DATA LLEP/
     *'     ','UP   ','UB   ','DN   ','DB   ','ST   ','SB   ','CH   ',
     *'CB   ','BT   ','BB   ','TP   ','TB   ','Y    ','YB   ','X    ',
     *'XB   ','GL   ','ERR  ','GM   ','ERR  ','NUE  ','ANUE ','E-   ',
     *'E+   ','NUM  ','ANUM ','MU-  ','MU+  ','NUT  ','ANUT ','TAU- ',
     *'TAU+ ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','KS   ',
     *'ERR  ','ERR  ','KL   ',
     *'UPSS ','UBSS ','DNSS ','DBSS ','STSS ','SBSS ','CHSS ','CBSS ',
     *'BTSS ','BBSS ','TPSS ','TBSS ','ERR  ','ERR  ','ERR  ','ERR  ',
     *'GLSS ','ERR  ','GMSS ','ERR  ','NESS ','ANESS','E-SS ','E+SS ',
     *'NMSS ','ANMSS','MU-SS','MU+SS','NTSS ','ANTSS','T-SS ','T+SS ',
     *'ERR  ','ERR  ','ERR  ','ERR  ','W+SS ','W-SS ','Z0SS ','ERR  ',
     *'W+   ','W-   ','H10  ','AH10 ','H20  ','AH20 ','H30  ','AH30 ',
     *'H4+  ','H4-  ','H5+  ','H5-  ','H6+  ','H6-  ','H7++ ','H7-- ',
     *'H8++ ','H8-- ','H9++ ','H9-- ','Z0   '/
C          0- MESON LABELS
      DATA LMES0/
     1'PI0  ','PI+  ','ETA  ','PI-  ','K+   ','K0   ','ETAP ','AK0  ',
     2'K-   ','AD0  ','D-   ','F-   ','ETAC ','F+   ','D+   ','D0   ',
     2'UB.  ','DB.  ','SB.  ','CB.  ','BB.  ','BC.  ','BS.  ','BD.  ',
     3'BU.  ','UT.  ','DT.  ','ST.  ','CT.  ','BT.  ','TT.  ','TB.  ',
     4'TC.  ','TS.  ','TD.  ','TU.  ','UY.  ','DY.  ','SY.  ','CY.  ',
     5'BY.  ','TY.  ','YY.  ','YT.  ','YB.  ','YC.  ','YS.  ','YD.  ',
     6'YU.  ','UX.  ','DX.  ','SX.  ','CX.  ','BX.  ','TX.  ','YX.  ',
     7'XX.  ','XY.  ','XT.  ','XB.  ','XC.  ','XS.  ','XD.  ','XU.  '/
C          1- MESON LABELS
      DATA LMES1/
     1'RHO0 ','RHO+ ','OMEG ','RHO- ','K*+  ','K*0  ','PHI  ','AK*0 ',
     2'K*-  ','AD*0 ','D*-  ','F*-  ','JPSI ','F*+  ','D*+  ','D*0  ',
     3'UB*  ','DB*  ','SB*  ','CB*  ','UPSL ','BC*  ','BS*  ','BD*  ',
     4'BU*  ','UT*  ','DT*  ','ST*  ','CT*  ','BT*  ','TT*  ','TB*  ',
     5'TC*  ','TS*  ','TD*  ','TU*  ','UY*  ','DY*  ','SY*  ','CY*  ',
     6'BY*  ','TY*  ','YY*  ','YT*  ','YB*  ','YC*  ','YS*  ','YD*  ',
     7'YU*  ','UX*  ','DX*  ','SX*  ','CX*  ','BX*  ','TX*  ','YX*  ',
     8'XX*  ','XY*  ','XT*  ','XB*  ','XC*  ','XS*  ','XD*  ','XU*  '/
C          1/2+ BARYON LABELS
      DATA LBAR0/
     1'ERR  ','P    ','N    ','ERR  ','ERR  ','S+   ','S0   ','S-   ',
     2'L    ','XI0  ','XI-  ','ERR  ','ERR  ','ERR  ','SC++ ','SC+  ',
     3'SC0  ','LC+  ','USC. ','DSC. ','SSC. ','SDC. ','SUC. ','UCC. ',
     4'DCC. ','SCC. ','ERR  ','ERR  ','ERR  ','ERR  ','UUB. ','UDB. ',
     5'DDB. ','DUB. ','USB. ','DSB. ','SSB. ','SDB. ','SUB. ','UCB. ',
     6'DCB. ','SCB. ','CCB. ','CSB. ','CDB. ','CUB. ','UBB. ','DBB. ',
     7'SBB. ','CBB. ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','UTT. ',
     8'UDT. ','DDT. ','DUT. ','UST. ','DST. ','SST. ','SDT. ','SUT. ',
     9'UCT. ','DCT. ','SCT. ','CCT. ','CST. ','CDT. ','CUT. ','UBT. ',
     1'DBT. ','SBT. ','CBT. ','BBT. ','BCT. ','BST. ','BDT. ','BUT. ',
     2'UTT. ','DTT. ','STT. ','CTT. ','BTT. ','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY. ','UDY. ','DDY. ','DUY. ','USY. ',
     4'DSY. ','SSY. ','SDY. ','SUY. ','UUX. ','UDX. ','DDX. ','DUX. ',
     5'USX. ','DSX. ','SSX. ','SDX. ','SUX. '/
      DATA LABAR0/
     1'ERR  ','AP   ','AN   ','ERR  ','ERR  ','AS-  ','AS0  ','AS+  ',
     2'AL   ','AXI0 ','AXI+ ','ERR  ','ERR  ','ERR  ','ASC--','ASC- ',
     3'ASC0 ','ALC- ','AUSC.','ADSC.','ASSC.','ASDC.','ASUC.','AUCC.',
     4'ADCC.','ASCC.','ERR  ','ERR  ','ERR  ','ERR  ','AUUB.','AUDB.',
     5'ADDB.','ADUB.','AUSB.','ADSB.','ASSB.','ASDB.','ASUB.','AUCB.',
     6'ADCB.','ASCB.','ACCB.','ACSB.','ACDB.','ACUB.','AUBB.','ADBB.',
     7'ASBB.','ACBB.','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','AUTT.',
     8'AUDT.','ADDT.','ADUT.','AUST.','ADST.','ASST.','ASDT.','ASUT.',
     9'AUCT.','ADCT.','ASCT.','ACCT.','ACST.','ACDT.','ACUT.','AUBT.',
     1'ADBT.','ASBT.','ACBT.','ABBT.','ABCT.','ABST.','ABDT.','ABUT.',
     2'AUTT.','ADTT.','ASTT.','ACTT.','ABTT.','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY.','AUDY.','ADDY.','ADUY.','AUSY.',
     4'ADSY.','ASSY.','ASDY.','ASUY.','AUUX.','AUDX.','ADDX.','ADUX.',
     5'AUSX.','ADSX.','ASSX.','ASDX.','ASUX.'/
C          3/2+ BARYON LABELS
      DATA LBAR1/
     1'DL++ ','DL+  ','DL0  ','DL-  ','ERR  ','S*+  ','S*0  ','S*-  ',
     2'ERR  ','XI*0 ','XI*- ','OM-  ','ERR  ','ERR  ','UUC* ','UDC* ',
     3'DDC* ','ERR  ','USC* ','DSC* ','SSC* ','ERR  ','ERR  ','UCC* ',
     4'DCC* ','SCC* ','CCC* ','ERR  ','ERR  ','ERR  ','UUB* ','UDB* ',
     5'DDB* ','ERR  ','USB* ','DSB* ','SSB* ','ERR  ','ERR  ','UCB* ',
     6'DCB* ','SCB* ','CCB* ','ERR  ','ERR  ','ERR  ','UBB* ','DBB* ',
     7'SBB* ','CBB* ','BBB* ','ERR  ','ERR  ','ERR  ','ERR  ','UTT* ',
     8'UDT* ','DDT* ','ERR  ','UST* ','DST* ','SST* ','ERR  ','ERR  ',
     9'UCT* ','DCT* ','SCT* ','CCT* ','ERR  ','ERR  ','ERR  ','UBT* ',
     1'DBT* ','SBT* ','CBT* ','BBT* ','ERR  ','ERR  ','ERR  ','ERR  ',
     2'UTT* ','DTT* ','STT* ','CTT* ','BTT* ','TTT* ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY* ','UDY* ','DDY* ','ERR  ','USY* ',
     4'DSY* ','SSY* ','ERR  ','ERR  ','UUX* ','UDX* ','DDX* ','ERR  ',
     5'USX* ','DSX* ','SSX* ','ERR  ','ERR  '/
      DATA LABAR1/
     1'ADL--','ADL- ','ADL0 ','ADL+ ','ERR  ','AS*- ','AS*0 ','AS*+ ',
     2'ERR  ','AXI*0','AXI*+','AOM+ ','ERR  ','ERR  ','AUUC*','AUDC*',
     3'ADDC*','ERR  ','AUSC*','ADSC*','ASSC*','ERR  ','ERR  ','AUCC*',
     4'ADCC*','ASCC*','ACCC*','ERR  ','ERR  ','ERR  ','AUUB*','AUDB*',
     5'ADDB*','ERR  ','AUSB*','ADSB*','ASSB*','ERR  ','ERR  ','AUCB*',
     6'ADCB*','ASCB*','ACCB*','ERR  ','ERR  ','ERR  ','AUBB*','ADBB*',
     7'ASBB*','ACBB*','ABBB*','ERR  ','ERR  ','ERR  ','ERR  ','AUTT*',
     8'AUDT*','ADDT*','ERR  ','AUST*','ADST*','ASST*','ERR  ','ERR  ',
     9'AUCT*','ADCT*','ASCT*','ACCT*','ERR  ','ERR  ','ERR  ','AUBT*',
     1'ADBT*','ASBT*','ACBT*','ABBT*','ERR  ','ERR  ','ERR  ','ERR  ',
     2'AUTT*','ADTT*','ASTT*','ACTT*','ABTT*','ATTT*','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY*','AUDY*','ADDY*','ERR  ','AUSY*',
     4'ADSY*','ASSY*','ERR  ','ERR  ','AUUX*','AUDX*','ADDX*','ERR  ',
     5'AUSX*','ADSX*','ASSX*','ERR  ','ERR  '/
      SAVE
C-----------------------------------------------------------------------
      CALL IDFLAV(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF ( ABS(ID) .LT.   100 ) GOTO 200
      IF ( ABS(ID) .LT. 1000 ) GOTO 100
      IF ( ID .NE. 0  .AND.  MOD(ID,100) .EQ. 0 ) GOTO 300
C  BARYONS
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-11
      IF     ( JSPIN .EQ. 0 ) THEN
        IF     ( ID .GT. 0 ) THEN
          IDLABL=LBAR0(INDEX)
        ELSEIF ( ID .LT. 0 ) THEN
          IDLABL=LABAR0(INDEX)
        ENDIF
      ELSEIF ( JSPIN .EQ. 1 ) THEN
        IF     ( ID .GT. 0 ) THEN
          IDLABL=LBAR1(INDEX)
        ELSEIF ( ID .LT. 0 ) THEN
          IDLABL=LABAR1(INDEX)
        ENDIF
      ENDIF
      RETURN
C  MESONS
100   CONTINUE
      I=MAX(IFL2,IFL3)
      J=-MIN(IFL2,IFL3)
      INDEX=MAX(I-1,J-1)**2+I+MAX(I-J,0)
      IF     ( JSPIN .EQ. 0 ) THEN
        IDLABL=LMES0(INDEX)
      ELSEIF ( JSPIN .EQ. 1 ) THEN
        IDLABL=LMES1(INDEX)
      ENDIF
      RETURN
C  QUARKS, LEPTONS, ETC.
200   CONTINUE
      INDEX=2*INDEX
      IF ( ID .LE. 0 ) INDEX=INDEX+1
      IDLABL=LLEP(INDEX)
      RETURN
300   I=ABS(IFL1)
      J=ABS(IFL2)
      INDEX=I+J*(J-1)/2
      IF ( ID .GT. 0 ) THEN
        IDLABL=LQQ(INDEX)
      ELSE
        IDLABL=LAQQ(INDEX)
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE IDMASS(ID,AMASS)

C-----------------------------------------------------------------------
C  RETURNS THE MASS OF THE PARTICLE WITH IDENT CODE ID.
C-----------------------------------------------------------------------
      PARAMETER (NMES=2)
      PARAMETER (NQLEP=41)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      REAL    AMBAR0(30),AMBAR1(30),AMLEP(52),AMMES0(10),AMMES1(10)

      DATA AMLEP/.3,.3,.5,1.6,4.9,30.,-1.,-1.,0.,0.,
     *0.,.511003E-3,0.,.105661,0.,1.807,3*-1.,.49767,.49767,
     *100.3,100.3,100.5,101.6,104.9,130.,2*-1.,100.,0.,
     *100.,100.005,100.,100.1,100.,101.8,2*-1.,100.,100.,
     *11*0./
C          0- MESON MASS TABLE
      DATA AMMES0/.13496,.13957,.5488,.49367,.49767,.9576,1.8633
     1,1.8683,2.030,2.976/
C          1- MESON MASS TABLE
      DATA AMMES1/.770,.770,.7826,.8881,.8922,1.0196,2.006,2.0086
     1,2.140,3.097/
C          1/2+ BARYON MASS TABLE
      DATA AMBAR0/-1.,.93828,.93957,2*-1.,1.1894,1.1925,1.1974
     1,1.1156,1.3149,1.3213,3*-1.,2.43,2.43,2.43,2.26
     2,2.50,2.50,2.60,2.40,2.40,3.55,3.55,3.70,4*-1./
C          3/2+ BARYON MASS TABLE
      DATA AMBAR1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     1,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     2,2.63,2.63,2.63,-1.,2.70,2.70,2.80,2*-1.,3.75,3.75
     3,3.90,4.80,3*-1./
      SAVE
C-----------------------------------------------------------------------
      IDABS=ABS(ID)
      I=IDABS/1000
      J1=IDABS-I*1000
      J = J1/100
      IF ( ID .NE. 0  .AND.  MOD(ID,100) .EQ. 0 ) GOTO 400
      K1= J1 - J*100
      K = K1/10
      JSPIN= K1 - K*10
      IF ( I .GT. 4   .OR.  J .GT. 4   .OR. K .GT .4 ) GOTO 300
      IF ( J .EQ. 0 ) GOTO 200
      IF ( I .EQ. 0 ) GOTO 100
C  BARYONS
C  ONLY X,Y BARYONS ARE QQX, QQY, Q=U,D,S.
      IF ( K .LE. 6 ) THEN
        INDEX=MAX(I-1,J-1)**2+I+MAX(I-J,0)+(K-1)*K*(2*K-1)/6
      ELSE
        INDEX=MAX(I-1,J-1)**2+I+MAX(I-J,0)+9*(K-7)+91
      ENDIF
      AMASS=(1-JSPIN)*AMBAR0(INDEX)+JSPIN*AMBAR1(INDEX)
      RETURN
C  MESONS
100   CONTINUE
      INDEX = J + K*(K-1)/2
      AMASS = (1-JSPIN)*AMMES0(INDEX) + JSPIN*AMMES1(INDEX)
      RETURN
C  QUARKS, LEPTONS, ETC
200   CONTINUE
      IF     ( IDABS .LT. 20 ) THEN
        INDEX = IDABS
C  DEFINE INDEX=20 FOR KS, INDEX=21 FOR KL
      ELSEIF ( ID .EQ. 20 ) THEN
        INDEX = 20
C  INDEX=NQLEP+1,...,NQLEP+11 FOR W+, HIGGS, Z0
      ELSEIF ( IDABS .LT. 80 ) THEN
        INDEX = IDABS+1
      ELSE
        INDEX = NQLEP+IDABS-79
      ENDIF
      AMASS=AMLEP(INDEX)
      RETURN
C  B AND T PARTICLES
300   CONTINUE
      AMASS=AMLEP(J)+AMLEP(K)-.03+.04*JSPIN
      IF ( I .NE. 0 ) AMASS=AMASS+AMLEP(I)
      RETURN
C  DIQUARKS
400   CONTINUE
      AMASS=AMLEP(I)+AMLEP(J)
      RETURN
      END
C=======================================================================

      SUBROUTINE IDMIX(IC,JSPIN,ICM,IDM)

C-----------------------------------------------------------------------
C  ACCOUNTS FOR FLAVOUR MIXING
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      REAL    PMIX1(3,2),PMIX2(3,2)
      INTEGER IC(2),ICM(2)
      DATA PMIX1 /.25,.25,.5,0.,.5,1./, PMIX2 /.5,.5,1.,0.,0.,1./
      SAVE
C-----------------------------------------------------------------------
      ICM(1)=0
      ICM(2)=0
      IDM=0
      I=IC(1)
      IF ( I .NE. IC(2) ) RETURN
      ID=0
      IF ( I .EQ. 100000 ) ID=1
      IF ( I .EQ.  10000 ) ID=2
      IF ( I .EQ.   1000 ) ID=3
      IF ( ID .EQ. 0 ) RETURN
      RND=RANGEN()
      IDM=INT(PMIX1(ID,JSPIN+1)+RND)+INT(PMIX2(ID,JSPIN+1)+RND)+1
      ICM(1)=10**(NFLAV-IDM)
      ICM(2)=IC(1)
      IDM=IDM*100+IDM*10+JSPIN
      RETURN
      END
C=======================================================================

      SUBROUTINE IDQUAC(I,NQ,NS,NA,JC)

C-----------------------------------------------------------------------
C  RETURNS QUARK CONTENT OF PTL I FROM /CPTL/ .
C        NQ = # QUARKS - # ANTIQUARKS
C        NS = # STRANGE QUARKS - # STRANGE ANTIQUARKS
C        NA = # QUARKS + # ANTIQUARKS
C        JC(NFLAV,2) = JC-TYPE PARTICLE IDENTIFICATION CODE.
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (NFLAV=6)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      INTEGER IC(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IF ( ABS(IDPTL(I)) .EQ. 20 ) THEN
        IF ( RANGEN() .LT. .5 ) THEN
          IDPTL(I)=-230
        ELSE
          IDPTL(I)=230
        ENDIF
        GOTO 9999
      ENDIF

      IF ( ABS(IDPTL(I)) .LT. 100 ) THEN
        NQ=0
        NS=0
        DO 1 N=1,NFLAV
          JC(N,1)=0
          JC(N,2)=0
 1      CONTINUE
      RETURN
      ENDIF

9999  IF ( IDPTL(I)/100000000 .NE. 7 ) THEN
        CALL IDTR4(IDPTL(I),IC)
        CALL IDDECO(IC,JC)
      ELSE
        CALL IDTRB(IBPTL(1,I),IBPTL(2,I),IBPTL(3,I),IBPTL(4,I),JC)
      ENDIF
      NA=0
      NQ=0
      DO 53 N=1,NFLAV
        NA=NA+JC(N,1)+JC(N,2)
        NQ=NQ+JC(N,1)-JC(N,2)
53    CONTINUE
      NS=JC(3,1)-JC(3,2)
      RETURN
      END
C=======================================================================

      SUBROUTINE IDRES(ID,AM,IDR,IADJ)

C-----------------------------------------------------------------------
C  RETURNS RESONANCE ID IDR CORRESPONDING TO MASS AM.
C  PERFORMS MASS ADJUSTMENT, IF NECESSARY (IF SO IADJ=1, 0 ELSE)
C-----------------------------------------------------------------------
      PARAMETER (MXINDX=1000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXRE=100)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      IF ( AM .EQ. 0. ) AM=1.E-5

      IDI=ID
      AMI=AM
      IDR=0
      IADJ=0

      IF ( ID .EQ. 0 ) GOTO 9999
      IDABS = ABS(ID)
      DO 5 K=1,MXMX
        DO 3 M=2,MXMA
          IF ( IDABS .EQ. IDMX(M,K) ) THEN
            ID=SIGN(IDMX(1,K)*10,ID)
            GOTO 5
          ENDIF
3       CONTINUE
5     CONTINUE

      IX=IDABS/10
      IF ( IX .LT. 1  .OR.  IX .GT. MXINDX ) THEN
        CALL UTSTOP('IDRES: IX OUT OF RANGE.                 ')
      ENDIF
      I=INDX(IX)
      IF ( I .LT. 1  .OR.  I .GT. MXRE ) THEN
        CALL UTSTOP('IDRES: PARTICLE NOT IN TABLE            ')
      ENDIF
      DO 1 J=1,MXMA-1
        IF ( AM .GE. REMA(I,J)  .AND.  AM .LE. REMA(I,J+1) ) THEN
          IF ( J .GT. 10 ) THEN
            CALL UTSTOP('IDRES: SPIN > 9                         ')
          ENDIF
          IDR=ID/10*10+SIGN(J-1,ID)
          GOTO 2
        ENDIF
 1    CONTINUE
      GOTO 9999
 2    CONTINUE

      DO 4 K=1,MXMX
        IF ( IX .EQ. IDMX(1,K) ) THEN
          IF ( J .LT. 1  .OR.  J .GT. MXMA-1 ) THEN
            CALL UTSTOP('IDRES: INDEX J OUT OF RANGE             ')
          ENDIF
          IF ( IDMX(J+1,K) .NE. 0 ) IDR = SIGN(IDMX(J+1,K),ID)
        ENDIF
 4    CONTINUE

      IY=MOD(IABS(IDR),10)
      IF ( IY .GT. MAXRES ) THEN
        IADJ=0
        IDR=0
        GOTO 9999
      ENDIF

      IF ( IY .NE. 0  .AND.  IY .NE. 1 ) GOTO 9999

      CALL IDMASS(IDR,AM)
      IF ( AM .LT. 0. ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'*****  ERROR IN IDRES: '
     *      ,'NEG MASS RETURNED FROM IDMASS'
        WRITE(IFCH,*)'ID,AM(INPUT):',IDI,AMI
        WRITE(IFCH,*)'IDR,AM:',IDR,AM
        CALL UTSTOP('IDRES: NEG MASS RETURNED FROM IDMASS    ')
      ENDIF
      IADJ=1

9999  ID=IDI
      IF ( ISH .LT. 93 ) RETURN
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'RETURN FROM IDRES. ID,AMI,AM,IDR,IADJ:'
      WRITE(IFCH,*)ID,AMI,AM,IDR,IADJ
      WRITE(IFCH,*)' '
      RETURN
      END
C=======================================================================

      SUBROUTINE IDRESI

C-----------------------------------------------------------------------
C  INITIALIZES /CREMA/
C-----------------------------------------------------------------------
      PARAMETER (MXINDX=1000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXRE=100)
      PARAMETER (N=29)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      REAL    REMAI(N,MXMA),REWII(N,MXMA)
      INTEGER ICREI(N,2*MXMA),IDMXI(MXMA,MXMX)

      DATA (IDMXI(J,1),J=1,MXMA)/ 11, 110, 111,   0,   0,   0,   0, 4*0/
      DATA (IDMXI(J,2),J=1,MXMA)/ 22, 220, 330, 331,   0,   0,   0, 4*0/
      DATA (IDMXI(J,3),J=1,MXMA)/123,2130,1230,1231,   0,   0,   0, 4*0/
      DATA (IDMXI(J,4),J=1,MXMA)/124,2140,1240,1241,   0,   0,   0, 4*0/
      DATA (IDMXI(J,5),J=1,MXMA)/134,3140,1340,1341,   0,   0,   0, 4*0/
      DATA (IDMXI(J,6),J=1,MXMA)/234,3240,2340,2341,   0,   0,   0, 4*0/
      DATA ((ICREI(K,M),M=1,2*MXMA),K=1,10)/
     *                 111,000000, 9*300000,    11*0,
     *                 222,000000, 9*030000,    11*0,
     *                 112,       10*210000,    11*0,
     *                 122,       10*120000,    11*0,
     *                 113,       10*201000,    11*0,
     *                 223,       10*021000,    11*0,
     *                 123,       10*111000,    11*0,
     *                 133,       10*102000,    11*0,
     *                 233,       10*012000,    11*0,
     *                 333,000000, 9*003000,    11*0/
      DATA ((ICREI(K,M),M=1,2*MXMA),K=11,20)/
     *                 114,       10*200100,    11*0,
     *                 124,       10*110100,    11*0,
     *                 224,       10*020100,    11*0,
     *                 134,       10*101100,    11*0,
     *                 234,       10*011100,    11*0,
     *                 334,       10*002100,    11*0,
     *                 144,       10*100200,    11*0,
     *                 244,       10*010200,    11*0,
     *                 344,       10*001200,    11*0,
     *                 444,000000, 9*000300,    11*0/
      DATA ((ICREI(K,M),M=1,2*MXMA),K=21,29)/
     *                  11,  10*100000,    0,   10*100000,
     *                  22,  10*001000,    0,   10*001000,
     *                  12,  10*100000,    0,   10*010000,
     *                  13,  10*100000,    0,   10*001000,
     *                  23,  10*010000,    0,   10*001000,
     *                  14,  10*100000,    0,   10*000100,
     *                  24,  10*010000,    0,   10*000100,
     *                  34,  10*001000,    0,   10*000100,
     *                  44,  10*000100,    0,   10*000100/

      DATA ((REMAI(K,M),M=1,MXMA),K=1,10)/
     *111.,0.000,1.425,1.660,1.825,2.000,0.000,0.000,0.000,0.000,0.000,
     *222.,0.000,1.425,1.660,1.825,2.000,0.000,0.000,0.000,0.000,0.000,
     *112.,1.080,1.315,1.485,1.575,1.645,1.685,1.705,1.825,2.000,0.000,
     *122.,1.080,1.315,1.485,1.575,1.645,1.685,1.705,1.825,2.000,0.000,
     *113.,1.300,1.500,1.700,1.850,2.000,0.000,0.000,0.000,0.000,0.000,
     *223.,1.300,1.500,1.700,1.850,2.000,0.000,0.000,0.000,0.000,0.000,
     *123.,1.117,1.300,1.395,1.465,1.540,1.655,1.710,1.800,1.885,2.000,
     *133.,1.423,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *233.,1.428,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *333.,0.000,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
      DATA ((REMAI(K,M),M=1,MXMA),K=11,20)/
     *114.,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *124.,2.345,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *224.,2.530,2.730,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *134.,2.450,2.600,2.800,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *234.,2.450,2.600,2.800,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *334.,2.700,2.900,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *144.,3.650,3.850,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *244.,3.650,3.850,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *344.,3.800,4.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     *444.,0.000,5.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
      DATA ((REMAI(K,M),M=1,MXMA),K=21,29)/
     * 11.,0.450,0.950,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 22.,0.750,0.965,1.080,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 12.,0.450,0.950,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 13.,0.700,1.050,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 23.,0.700,1.050,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 14.,1.935,2.077,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 24.,1.938,2.079,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 34.,2.085,2.195,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     * 44.,3.037,3.158,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

      DATA ((REWII(K,M),M=1,MXMA),K=1,5)/
     *111.,0.000E+00,0.115E+00,0.140E+00,0.250E+00,0.250E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *222.,0.000E+00,0.115E+00,0.140E+00,0.250E+00,0.250E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *112.,0.000E+00,0.115E+00,0.200E+00,0.140E+00,0.140E+00,
     *     0.145E+00,0.250E+00,0.140E+00,0.250E+00,0.000E+00,
     *122.,0.000E+00,0.115E+00,0.200E+00,0.140E+00,0.140E+00,
     *     0.145E+00,0.250E+00,0.140E+00,0.250E+00,0.000E+00,
     *113.,0.824E-14,0.036E+00,0.080E+00,0.100E+00,0.170E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      DATA ((REWII(K,M),M=1,MXMA),K=6,10)/
     *223.,0.445E-14,0.039E+00,0.080E+00,0.100E+00,0.170E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *123.,0.250E-14,0.890E-05,0.036E+00,0.040E+00,0.016E+00,
     *     0.090E+00,0.080E+00,0.100E+00,0.145E+00,0.170E+00,
     *133.,0.227E-14,0.009E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *233.,0.400E-14,0.010E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *333.,0.000E+00,0.800E-14,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      DATA ((REWII(K,M),M=1,MXMA),K=11,15)/
     *114.,0.400E-11,0.010E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *124.,0.400E-11,0.400E-11,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *224.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *134.,0.150E-11,0.400E-11,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *234.,0.150E-11,0.400E-11,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      DATA ((REWII(K,M),M=1,MXMA),K=16,20)/
     *334.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *144.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *244.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *344.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     *444.,0.400E-11,0.010E+00,0.010E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      DATA ((REWII(K,M),M=1,MXMA),K=21,25)/
     * 11.,0.757E-08,0.153E+00,0.057E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 22.,0.105E-05,0.210E-03,0.034E+00,0.004E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 12.,0.000E+00,0.153E+00,0.057E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 13.,0.000E+00,0.051E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 23.,0.197E-02,0.051E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      DATA ((REWII(K,M),M=1,MXMA),K=26,29)/
     * 14.,0.154E-11,0.002E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 24.,0.615E-12,0.002E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 34.,0.150E-11,0.020E+00,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
     * 44.,0.010E+00,0.068E-03,0.000E+00,0.000E+00,0.000E+00,
     *     0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00/
      SAVE
C-----------------------------------------------------------------------
      NN=N
      DO 3 I=1,MXINDX
        INDX(I)=0
 3    CONTINUE
      DO 44 M=1,MXMA
        DO 4 K=1,MXRE
          REMA(K,M)=0.
 4      CONTINUE
44    CONTINUE
      DO 22 I=1,MXMX
        DO 2 J=1,MXMA
          IDMX(J,I)=IDMXI(J,I)
 2      CONTINUE
22    CONTINUE

      IF ( NN .GT. MXRE ) THEN
        CALL UTSTOP('IDRESI: DIMENSION MXRE TOO SMALL        ')
      ENDIF
      DO 1 K=1,N
        IX=NINT(REMAI(K,1))
        IX2=NINT(REWII(K,1))
        IX3=ICREI(K,1)
        IF ( IX .NE. IX2 ) THEN
          CALL UTSTOP('IDRESI: IX /= IX2                       ')
        ENDIF
        IF ( IX .NE. IX3 ) THEN
          CALL UTSTOP('IDRESI: IX /= IX3                       ')
        ENDIF
        IF ( IX .LT. 1  .OR.  IX .GT. MXINDX ) THEN
          CALL UTSTOP('IDRESI: IX OUT OF RANGE.                ')
        ENDIF
        INDX(IX)=K
        REMA(K,1)=0.
        REWI(K,1)=0.
        ICRE1(K,1)=0
        ICRE2(K,1)=0
        DO 5 M=2,MXMA
          REMA(K,M)=REMAI(K,M)
          REWI(K,M)=REWII(K,M)
          ICRE1(K,M)=ICREI(K,M)
          ICRE2(K,M)=ICREI(K,MXMA+M)
 5      CONTINUE
 1    CONTINUE

      INDX(33) =INDX(22)
      INDX(213)=INDX(123)
      INDX(214)=INDX(124)
      INDX(314)=INDX(134)
      INDX(324)=INDX(234)

      RETURN
      END
C=======================================================================

      SUBROUTINE IDSPIN(II,IC,ID,JSPIN)

C-----------------------------------------------------------------------
C  DETERMINES PARTICLE SPIN
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU

      REAL    PSPIN1(8)
      INTEGER IC(2),JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      PSPIN1(1)=PSPINL
      PSPIN1(2)=PSPINL
      PSPIN1(3)=PSPINL
      PSPIN1(4)=PSPINH
      PSPIN1(5)=PSPINH
      PSPIN1(6)=PSPINH
      PSPIN1(7)=PSPINH
      PSPIN1(8)=PSPINH
      CALL IDDECO(IC,JC)
      IHIGH=0
      IF ( II .EQ. 1 ) THEN
        DO 4 I=1,NFLAV
          IF ( JC(I,1) .NE. 0 ) IHIGH=I
          IF ( JC(I,2) .NE. 0 ) IHIGH=I
 4      CONTINUE
      ELSE
        CALL IDFLAV(ID,I1,I2,I3,IDU1,IDU2)
        DO 5 I=1,NFLAV
          IF ( I.EQ.ABS(I1)   .OR.  I .EQ. ABS(I2)
     *                        .OR.  I .EQ. ABS(I3) ) IHIGH=I
5       CONTINUE
      ENDIF
      JSPIN=INT(RANGEN()+PSPIN1(IHIGH))
      RETURN
      END
C=======================================================================

      SUBROUTINE IDTAU(ID,P4,P5,TAUGM)

C-----------------------------------------------------------------------
C  RETURNS LIFETIME*GAMMA FOR ID WITH ENERGY P4, MASS P5
C-----------------------------------------------------------------------
      PARAMETER (MXINDX=1000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXRE=100)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      SAVE
C-----------------------------------------------------------------------
      IDABS = ABS(ID)
      IF     ( IDABS .LT. 100  .AND.  ID .NE. 20 ) THEN
        WI=0.
      ELSEIF ( ID .EQ. 20 ) THEN
        WI=.197/2.675E13
      ELSEIF ( IDABS .LT. 100000000 ) THEN
        IX=IDABS/10
        IF ( IX .LT. 1  .OR.  IX .GT. MXINDX ) THEN
          CALL UTSTOP('IDTAU: IX OUT OF RANGE.                 ')
        ENDIF
        II=INDX(IX)
        JJ=MOD(IDABS,10)+2
        DO 75 IMX=1,MXMX
          DO 76 IMA=2,MXMA
            IF ( IDABS .EQ. IDMX(IMA,IMX) ) JJ=IMA
76        CONTINUE
75      CONTINUE
        IF ( II.LT.1 .OR. II.GT.MXRE .OR. JJ.LT.1 .OR. JJ.GT.MXMA ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'ID,II,JJ:',ID,'   ',II,JJ
          CALL UTSTOP('IDTAU: II OR JJ OUT OF RANGE            ')
        ENDIF
        WI=REWI(II,JJ)
      ELSE
        TAUZ=TAUNLL
C-C     TAUZ=MIN( 9./P5**2, TAUZ )
C-C     TAUZ=MAX( .2, TAUZ )
        WI=.197/TAUZ
      ENDIF
      IF ( WI .EQ. 0. ) THEN
        TAU=AINFIN
        TAUGM=AINFIN
        RETURN
      ELSE
        TAU=.197/WI
        IF ( TAU .GE. AINFIN ) THEN
          TAUGM = AINFIN
          RETURN
        ENDIF
      ENDIF
      IF ( P5 .EQ. 0. ) THEN
        GM=AINFIN
        TAUGM=AINFIN
      ELSE
        GM=P4/P5
        IF ( GM .GE. AINFIN ) THEN
          TAUGM = AINFIN
          RETURN
        ELSE
          TAUGM=TAU*GM
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================

      INTEGER FUNCTION IDTRA(IC,IER,IRES,IMIX)

C-----------------------------------------------------------------------
C  TRANFORMS FROM WERNER-ID TO PAIGE-ID
C-----------------------------------------------------------------------
      PARAMETER (NIDT=44)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU

      INTEGER IC(2),ICM(2),IDT(3,NIDT)
      DATA IDT/
     * 100000,000000,   1   ,010000,000000,   2   ,001000,000000,   3
     *,000100,000000,   4
     *,200000,000000,1100   ,110000,000000,1200   ,020000,000000,2200
     *,101000,000000,1300   ,011000,000000,2300   ,002000,000000,3300
     *,100100,000000,1400   ,010100,000000,2400   ,001100,000000,3400
     *,000200,000000,4400
     *,100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/
      SAVE
C-----------------------------------------------------------------------
      IF ( IC(1) .EQ. 0  .AND.  IC(2) .EQ. 0 ) THEN
        IDTRA=0
        RETURN
      ENDIF
      DO 1 I=1,NIDT
        IF ( IC(1).EQ.IDT(1,I) .AND. IC(2).EQ.IDT(2,I) ) THEN
          IDTRA=IDT(3,I)
          GOTO 2
        ENDIF
        IF ( IC(2).EQ.IDT(1,I) .AND. IC(1).EQ.IDT(2,I) ) THEN
          IDTRA=-IDT(3,I)
          GOTO 2
        ENDIF
 1    CONTINUE
      IDTRA=0
 2    CONTINUE
      IF ( IDTRA .NE. 0 ) THEN
        ISI=SIGN(1,IDTRA)
        JSPIN=0
        IF ( MOD(IDTRA,10).EQ.0 .AND. IRES.EQ.1 )
     *                                     CALL IDSPIN(1,IC,IDU,JSPIN)
      ELSE
        ISI=1
        JSPIN=0
      ENDIF
      IF     ( IMIX .EQ. 3 ) THEN
        IF     ( IDTRA .EQ. 220 ) THEN
          IDTRA=110
        ELSEIF ( IDTRA .EQ. 330 ) THEN
          IDTRA=220
        ENDIF
      ELSEIF ( IMIX .EQ. 2 ) THEN
        IF     ( IDTRA .EQ. 220 ) THEN
          IDTRA=110
        ELSEIF ( IDTRA .EQ. 330 ) THEN
          IDTRA=110
        ENDIF
      ELSEIF ( IMIX .EQ. 1 ) THEN
        CALL IDMIX(IC,JSPIN,ICM,IDTRAM)
        IF ( IDTRAM .NE. 0 ) IDTRA=IDTRAM
        IF ( JSPIN .EQ. 0 ) THEN
          IF ( RANGEN() .LT. PISPN ) THEN
            IF     ( ABS(IDTRA). EQ. 1230 ) THEN
              IDTRA=ISI*2130
            ELSEIF ( ABS(IDTRA) .EQ. 1240 ) THEN
              IDTRA=ISI*2140
            ELSEIF ( ABS(IDTRA) .EQ. 1340 ) THEN
              IDTRA=ISI*3140
            ELSEIF ( ABS(IDTRA) .EQ. 2340 ) THEN
              IDTRA=ISI*3240
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      IF ( IDTRA .NE. 0 ) IDTRA=IDTRA+JSPIN*ISI
      IF ( IDTRA .NE. 0 ) RETURN
      IF ( IER .NE. 1 ) RETURN
      JERR=JERR+1
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'***** ERROR IN IDTRA: UNKNOWN CODE'
      WRITE(IFCH,*)'IC = ',IC
      WRITE(IFCH,*)' '
      RETURN
      END
C=======================================================================

      INTEGER FUNCTION IDTRAI(NUM,ID,IER)

C-----------------------------------------------------------------------
      PARAMETER (NIDT=44)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      INTEGER IDT(3,NIDT)

      DATA IDT/
     * 100000,000000,   1   ,010000,000000,   2   ,001000,000000,   3
     *,000100,000000,   4
     *,200000,000000,1100   ,110000,000000,1200   ,020000,000000,2200
     *,101000,000000,1300   ,011000,000000,2300   ,002000,000000,3300
     *,100100,000000,1400   ,010100,000000,2400   ,001100,000000,3400
     *,000200,000000,4400
     *,100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/
      SAVE
C-----------------------------------------------------------------------
      IDABS = ABS(ID)
      DO 1 I=1,NIDT
        IF ( IDABS .EQ. IDT(3,I) ) THEN
          IF ( ID .LT. 0 ) THEN
            IDTRAI=IDT(3-NUM,I)
          ELSE
            IDTRAI=IDT(NUM,I)
          ENDIF
          RETURN
        ENDIF
1     CONTINUE
      IDTRAI=0
      IF ( IER .NE. 1 ) RETURN
      JERR=JERR+1
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'***** ERROR IN IDTRAI: UNKNOWN CODE'
      WRITE(IFCH,*)'ID = ',ID
      WRITE(IFCH,*)' '
      RETURN
      END
C=======================================================================

      SUBROUTINE IDTRB(IB1,IB2,IB3,IB4,JC)

C-----------------------------------------------------------------------
C  ID TRANSFORMATION IB -> JC
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      INTEGER JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      JC(1,1)=IB1/10000
      JC(2,1)=IB2/10000
      JC(3,1)=IB3/10000
      JC(4,1)=IB4/10000
      JC(5,1)=0
      JC(6,1)=0
      JC(1,2)=MOD(IB1,10000)
      JC(2,2)=MOD(IB2,10000)
      JC(3,2)=MOD(IB3,10000)
      JC(4,2)=MOD(IB4,10000)
      JC(5,2)=0
      JC(6,2)=0
      RETURN
      END
C=======================================================================

      SUBROUTINE IDTRBI(JC,IB1,IB2,IB3,IB4)

C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      INTEGER JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IB1=JC(1,1)*10000+JC(1,2)
      IB2=JC(2,1)*10000+JC(2,2)
      IB3=JC(3,1)*10000+JC(3,2)
      IB4=JC(4,1)*10000+JC(4,2)
      IB5=JC(5,1)*10000+JC(5,2)
      IB6=JC(6,1)*10000+JC(6,2)
      IF ( IB5 .NE. 0  .OR.  IB6 .NE. 0 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'***** ERROR IN IDTRBI: BOTTOM OR TOP QUARKS'
        WRITE(IFCH,*)'JC:'
        WRITE(IFCH,*)JC
        CALL UTSTOP('IDTRBI: BOTTOM OR TOP QUARKS            ')
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE IDTR4(ID,IC)

C-----------------------------------------------------------------------
C  TRANSFORMS GENERALIZED PAIGE_ID -> WERNER_ID  (FOR < 4 FLV)
C-----------------------------------------------------------------------
      PARAMETER (MXINDX=1000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXRE=100)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      INTEGER IC(2)
      SAVE
C-----------------------------------------------------------------------
      IDABS = ABS(ID)
      IF ( IDABS .LT. 100000000 ) THEN
        IX=IDABS/10
        IF ( IX .LT. 1  .OR.  IX .GT. MXINDX ) GOTO 9999
        II=INDX(IX)
        IF ( II .EQ. 0 ) GOTO 9998
        JJ=IDABS-10*IX + 2
        DO 28 IMX=1,MXMX
          DO 27 IMA=2,MXMA
            IF ( IDABS .EQ. IDMX(IMA,IMX) ) THEN
              JJ=IMA
              GOTO 29
            ENDIF
27        CONTINUE
28      CONTINUE
29      IF ( ID .GT. 0 ) THEN
          IC(1)=ICRE1(II,JJ)
          IC(2)=ICRE2(II,JJ)
        ELSE
          IC(2)=ICRE1(II,JJ)
          IC(1)=ICRE2(II,JJ)
        ENDIF
        IF ( IC(1) .EQ. 100000  .AND.  IC(2) .EQ. 100000
     *                          .AND.  RANGEN() .LT. 0.5 ) THEN
          IC(1)=010000
          IC(2)=010000
        ENDIF
      ELSEIF ( MOD(ID/100000000,10) .EQ. 8 ) THEN
        IC(1)=MOD(ID,100000000)/10000*100
        IC(2)=MOD(ID,10000)*100
      ELSE
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'***** ID: ',ID
        CALL UTSTOP('IDTR4: UNRECOGNIZED ID                  ')
      ENDIF
      RETURN

9998  WRITE(IFCH,*)' '
      WRITE(IFCH,*)'ID: ',ID
      CALL UTSTOP('IDTR4: INDX=0.                          ')
      RETURN

9999  WRITE(IFCH,*)' '
      WRITE(IFCH,*)'ID: ',ID
      CALL UTSTOP('IDTR4: IX OUT OF RANGE.                 ')
      RETURN
      END
C=======================================================================

      SUBROUTINE JAMBR1(NS,NSG,IRET)

C-----------------------------------------------------------------------
C  "AMOR" (KOCH+WERNER, 89)
C  BREAKS STRING NS ACCORDING TO A-M MODEL.
C  NS: POINTS TO THE CURRENT FRAGMENTING STRING.
C  NSG: POINTS TO THE LAST PRODUCED SUBSTRING (SUCC INCREASED)
C      INPUT:
C  PSG(,NS): MOMENTUM OF STRING NS IN PP-CMS
C  PJT(,2*NS-1),PJT(,2*NS): MOMENTUM OF END OF STRING NS IN PP-CMS
C      OUTPUT:
C  XORSG(,NSG): ORIGIN OF SUBSTRING NSG IN PP-CMS
C  PJT(,2*NSG-1),PJT(,2*NSG): MOMENTUM OF END OF SUBSTRING NSG IN PP-CMS
C  ICJT(,2*NSG-1),ICJT(,2*NSG): IC-CODE  OF END OF SUBSTRG NSG IN PP-CMS
C  PSG(,NSG): MOMENTUM OF SUBSTRING NSG IN PP-CMS
C  XBKPTL(,NPTLC-NPTL2): BREAKPOINT OF NS IN PP-CMS
C  ISPTL(,NPTLC-NPTL2): 1 IF STRING NS BREAKS
C-----------------------------------------------------------------------
      PARAMETER (MXPC=500)
      PARAMETER (MXSG=500)
      PARAMETER (MXJT=2*MXSG)
      PARAMETER (NFLAV=6)
      PARAMETER (NPTF=129)
      COMMON /CDELRE/  DELRER
      COMMON /CJAMBR/  NPTLC,NPTL2
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      DOUBLE PRECISION XBKPTL
      COMMON /CPC/     XBKPTL(2,MXPC),ISPTL(MXPC)
      COMMON /CPTF/    FPTFS,FPTFSS,FPTFU,FPTFUS,FPTFUU
     *                ,QPTFS(NPTF),QPTFSS(NPTF),QPTFU(NPTF),QPTFUS(NPTF)
     *                ,QPTFUU(NPTF),XPTF(NPTF)
      DOUBLE PRECISION PJT,PSG,ROTSG,XORSG
      COMMON /CSG/     PJT(5,MXJT),PSG(5,MXSG),ROTSG(3,MXSG)
     *                ,XORSG(4,MXSG)
     *                ,ICJT(2,MXJT),IORSG(MXSG),ISG(MXSG)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /QUARKM/  SMAS,SSMAS,USMAS,UUMAS

      DOUBLE PRECISION DAUXI1,DAUXI2,ETA,EETA,ETAM,ETAP
     *                ,PM,PP,PT,PX,PY,P0M,P0P,XBK(4),XOR1(4),XOR2(4)
      INTEGER          IC(2),ICM(2),ICMP(2),ICMS(2)
     *                ,ICP(2),ICPM(2),ICPS(2),ICX(2),ICY(2)
     *                ,JC(NFLAV,2)
     *                ,JCM(NFLAV,2),JCMP(NFLAV,2),JCMS(NFLAV,2)
     *                ,JCP(NFLAV,2),JCPM(NFLAV,2),JCPS(NFLAV,2)
     *                ,JCX(NFLAV,2),JCY(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      NCNT3=0
      NSG0=NSG
ctp060203 9993  NSG=NSG0
      NSG=NSG0

      IRET=0
      ICP(1)=ICJT(1,2*NS-1)
      ICP(2)=ICJT(2,2*NS-1)
      ICM(1)=ICJT(1,2*NS)
      ICM(2)=ICJT(2,2*NS)
      CALL IDDECO(ICP,JCP)
      CALL IDDECO(ICM,JCM)
      NP=0
      DO 7 NF=1,NFLAV
        JC(NF,1)=JCP(NF,1)+JCM(NF,1)
        JC(NF,2)=JCP(NF,2)+JCM(NF,2)
        NP=NP+JCP(NF,1)-JCP(NF,2)
7     CONTINUE
      CALL IDENCO(JC,IC,IRETEN)
      ID=IDTRA(IC,0,0,3)
      AMMS=UTAMNX(JCP,JCM)
      AM=PSG(5,NS)

C  SPLIT STRING
C  ------------
      J1 = 2*NS
      J2 = J1-1
      DAUXI1= PJT(4,J1)+PJT(4,J2)
      DAUXI2= PJT(3,J1)+PJT(3,J2)
      P0P = DAUXI1 + DAUXI2
      P0M = DAUXI1 - DAUXI2
9994  NCNT3=NCNT3+1
      IF ( NCNT3 .GT. 100 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JAMBR1')
          WRITE(IFCH,*)'*****  SPLIT KINEM NOT POSSIBLE.'
          WRITE(IFCH,112)
     *     (ICJT(J,2*NS-1),J=1,2),(ICJT(J,2*NS),J=1,2),PSG(5,NS)
112       FORMAT(1X,2I8,4X,2I8,4X,F7.2)
          CALL UTMSGF
        ENDIF
        IF ( NS .EQ. 1 ) THEN
          IRET=9999
          RETURN
        ENDIF
        IRET=9996
        RETURN
      ENDIF
      DO 17 NF=1,NFLAV
        JCPM(NF,1)=0
        JCPM(NF,2)=0
        JCMP(NF,1)=0
        JCMP(NF,2)=0
17    CONTINUE

C  DETERMINE FLAVOUR
C  -----------------
      NQU=0
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)'ORDINARY STRING FRAGMENTATION'
        WRITE(IFCH,*)' '
      ENDIF
      IF ( RANGEN() .LE. PDIQUA ) THEN
        NQU=2
      ELSE
        NQU=1
      ENDIF
      IF     ( MOD(NP+NQU,3) .EQ. 0 ) THEN
        II=1
      ELSEIF ( MOD(NP-NQU,3) .EQ. 0 ) THEN
        II=2
      ELSE
        CALL UTSTOP('JAMBR1: NO SINGLET CONSTRUCTION POSSIBLE')
      ENDIF
      IFLTT=0
      DO 8 N=1,NQU
        IFL=INT(RANGEN()/PUD)+1
        IFLTT=IFLTT*10+(IFL+1)/2
        JCPM(IFL,II)=JCPM(IFL,II)+1
        JCMP(IFL,3-II)=JCMP(IFL,3-II)+1
8     CONTINUE
      CALL IDENCO(JCPM,ICPM,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('JAMBR1: IDENCO(JCPM... RET.CODE=1       ')
      ENDIF
      CALL IDENCO(JCMP,ICMP,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('JAMBR1: IDENCO(JCMP... RET.CODE=1       ')
      ENDIF
      DO 25 NF=1,NFLAV
        JCPS(NF,1)=JCP(NF,1)+JCPM(NF,1)
        JCPS(NF,2)=JCP(NF,2)+JCPM(NF,2)
        JCMS(NF,1)=JCM(NF,1)+JCMP(NF,1)
        JCMS(NF,2)=JCM(NF,2)+JCMP(NF,2)
25    CONTINUE
      CALL IDENCO(JCPS,ICPS,IRETEN)
      CALL IDENCO(JCMS,ICMS,IRETEN)
      CALL IDCOMK(ICPS)
      CALL IDCOMK(ICMS)

C  CALCULATE P+,P-,PT OF STRING BREAKING
C  -------------------------------------
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,109)ICM(1),ICMP(1),ICPM(1),ICP(1)
     *     ,ICM(2),ICMP(2),ICPM(2),ICP(2)
109     FORMAT(1X,'FLAVORS:',2(I11,I7),/,9X,2(I11,I7),/)
        WRITE(IFCH,*)'IFLTT:',IFLTT
        WRITE(IFCH,*)' '
      ENDIF
      IDP=IDTRA(ICPS,0,0,3)
      IDM=IDTRA(ICMS,0,0,3)
      AMMP=UTAMNY(JCP,JCPM)
      AMMM=UTAMNY(JCM,JCMP)
      R = RANGEN()
      IF     ( IFLTT .EQ .1 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFU ,R*QPTFU(NPTF))
C##       WRITE(IFCH,*)'JAMBR1:PT(OLD)=',PT
        ELSE
          RPT = R*FPTFU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(1.+RPT*2./AUXIL))
        ENDIF
      ELSEIF ( IFLTT .EQ. 2 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFS ,R*QPTFS(NPTF))
        ELSE
          RPT = R*FPTFS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(SMAS**2/AUXIL)+RPT*2./AUXIL)-SMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 11 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUU,R*QPTFUU(NPTF))
        ELSE
          RPT = R*FPTFUU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(UUMAS**2/AUXIL)+RPT*2./AUXIL)-UUMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 12  .OR.  IFLTT .EQ. 21 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUS,R*QPTFUS(NPTF))
        ELSE
          RPT = R*FPTFUS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(USMAS**2/AUXIL)+RPT*2./AUXIL)-USMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 22 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFSS,R*QPTFSS(NPTF))
        ELSE
          RPT = R*FPTFSS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(SSMAS**2/AUXIL)+RPT*2./AUXIL)-SSMAS**2)
        ENDIF
      ENDIF
      LO=1
      GOTO 48
47    LO=LO+1
      PT=RANGEN()*PT
48    CONTINUE
      PHI=2.D0*PI*RANGEN()
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      TMMP=    (PT**2+AMMP**2)
      TMMM=    (PT**2+AMMM**2)
      AREA=-LOG(RANGEN())/PAREA
      AREART = SQRT(AREA)
      ETAM=LOG((TMMM   +AREA)/(AREART    *P0M))
      ETAP=LOG((AREART    *P0P)/(TMMP   +AREA))
      IF ( ETAM .GT. ETAP ) THEN
        IF ( LO .LT. 5 ) GOTO 47
        GOTO 9994
      ENDIF
      ETA=ETAM+RANGEN()*(ETAP-ETAM)
      AMP=SQRT(P0P*AREART    *EXP(-ETA)-AREA-PT**2)
      AMM=SQRT(P0M*AREART    *EXP( ETA)-AREA-PT**2)
      CALL IDRES(IDP,AMP,IDPR,IADJP)
      CALL IDRES(IDM,AMM,IDMR,IADJM)
      R=RANGEN()
      IF ( IDPR .EQ. 110  .AND. R .LT. 0.5 ) THEN
        IDP=220
        AMP=.549
        IF ( R .LT. 0.6666667 ) AMP=.958
        CALL IDRES(IDP,AMP,IDPR,IADJP)
        IADJP=1
      ENDIF
      R=RANGEN()
      IF ( IDMR .EQ. 110 .AND. R .LT. 0.5 ) THEN
        IDM=220
        AMM=.549
        IF ( R .LT. 0.6666667 ) AMM=.958
        CALL IDRES(IDM,AMM,IDMR,IADJM)
        IADJM=1
      ENDIF
      TMP2=(PT**2+AMP**2)
      TMM2=(PT**2+AMM**2)
      IF ( IADJP .EQ. 1  .AND.  IADJM .NE. 1 ) THEN
        ETA=LOG((AREART  *P0P)/(TMP2+AREA))
        IF ( ETA .LT. ETAM ) GOTO 9994
        AMM=SQRT(P0M*AREART  *EXP(ETA)-AREA-PT**2)
        CALL IDRES(IDM,AMM,IDMR,IADJM)
        TMM2=(PT**2+AMM**2)
      ENDIF
      IF ( IADJP .NE. 1  .AND.  IADJM .EQ. 1 ) THEN
        ETA=LOG((TMM2+AREA)/(AREART  *P0M))
        IF ( ETA .GT. ETAP ) GOTO 9994
        AMP=SQRT(P0P*AREART   *EXP(-ETA)-AREA-PT**2)
        CALL IDRES(IDP,AMP,IDPR,IADJP)
        TMP2=(PT**2+AMP**2)
      ENDIF
      IF ( IADJP .EQ. 1  .AND.  IADJM .EQ. 1 ) THEN
        TM=(P0P*P0M-TMM2-TMP2)*0.5
        IF ( TM .LT. 0. ) GOTO 9994
        IF ( TM**2-TMP2*TMM2 .LT. 0. ) GOTO 9994
        AREA=TM-SQRT(TM**2  -TMP2*TMM2)
        AREART = SQRT(AREA)
        EETA=P0P*AREART/(TMP2+AREA)
        IF ( EETA .LE. 0.D0 ) GOTO 9994
        ETA=LOG(EETA)
      ENDIF

      PP=AREART*EXP(ETA)
      PM=AREART*EXP(-ETA)
      IF ( P0P-PP-PT**2/PM .LT. 0.D0   .OR.
     *     P0M-PM-PT**2/PP .LT. 0.D0 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JAMBR1')
          WRITE(IFCH,*)'*****  NEGATIVE JET ENERGY. SPLIT AGAIN.'
          IF ( P0P-PP-PT**2/PM .LT. 0.D0 )
     *      WRITE(IFCH,*)'P0P: ',P0P,'   PP+PT**2/PM: ',PP+PT**2/PM
          IF ( P0M-PM-PT**2/PP .LT. 0.D0 )
     *      WRITE(IFCH,*)'P0M: ',P0M,'   PM+PT**2/PP: ',PM+PT**2/PP
          CALL UTMSGF
        ENDIF
        GOTO 9994
      ENDIF

C  DETERMINE BREAK PNT AND NEW ORIGINS
C  -----------------------------------
      XOR1(1)=0.D0
      XOR1(2)=0.D0
      XOR1(3)=PP*0.5D0/TENSN
      XOR1(4)=PP*0.5D0/TENSN
      IF ( ISH .GE. 93 ) WRITE(IFCH,118)XOR1
118   FORMAT(' NEW ORIGIN +: ',13X,4F7.2)
      CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *             ,XOR1(1),XOR1(2),XOR1(3))
      CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *             ,XOR1(1),XOR1(2),XOR1(3),XOR1(4))
      XOR1(1)=XOR1(1)+XORSG(1,NS)
      XOR1(2)=XOR1(2)+XORSG(2,NS)
      XOR1(3)=XOR1(3)+XORSG(3,NS)
      XOR1(4)=XOR1(4)+XORSG(4,NS)
      IF ( ISH .GE. 93 ) WRITE(IFCH,117)XOR1
      XBK(1)=0.D0
      XBK(2)=0.D0
      XBK(3)=0.5D0/TENSN*(PP-PM)
      XBK(4)=0.5D0/TENSN*(PP+PM)
      IF ( ISH .GE. 93 ) WRITE(IFCH,116)XBK
116   FORMAT(' BREAKING PNT: ',13X,4F7.2)
      CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *             ,XBK(1),XBK(2),XBK(3))
      CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *             ,XBK(1),XBK(2),XBK(3),XBK(4))
      XBK(1)=XBK(1)+XORSG(1,NS)
      XBK(2)=XBK(2)+XORSG(2,NS)
      XBK(3)=XBK(3)+XORSG(3,NS)
      XBK(4)=XBK(4)+XORSG(4,NS)
      IF ( ISH .GE. 93 ) WRITE(IFCH,117)XBK
117   FORMAT(28X,4F7.2,/)
      XOR2(1)=0.D0
      XOR2(2)=0.D0
      XOR2(3)=-PM*0.5D0/TENSN
      XOR2(4)= PM*0.5D0/TENSN
      IF ( ISH .GE. 93 ) WRITE(IFCH,119)XOR2
119   FORMAT(' NEW ORIGIN -: ',13X,4F7.2)
      CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *              ,XOR2(1),XOR2(2),XOR2(3))
      CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *              ,XOR2(1),XOR2(2),XOR2(3),XOR2(4))
      XOR2(1)=XOR2(1)+XORSG(1,NS)
      XOR2(2)=XOR2(2)+XORSG(2,NS)
      XOR2(3)=XOR2(3)+XORSG(3,NS)
      XOR2(4)=XOR2(4)+XORSG(4,NS)
      IF ( ISH .GE. 93 ) WRITE(IFCH,117)XOR2

C  STORE THE TWO SUBSTRINGS
C  ------------------------
      NSGB=NSG
      DO 9998 M=1,2
        NSG=NSG+1
        IF ( NSG .GT. MXSG ) THEN
          CALL UTSTOP('JAMBR1: NSG>MXSG                        ')
        ENDIF
        ISG(NSG)=NPTLC
        IF ( M .EQ. 1 ) THEN
          XORSG(1,NSG)=XOR1(1)
          XORSG(2,NSG)=XOR1(2)
          XORSG(3,NSG)=XOR1(3)
          XORSG(4,NSG)=XOR1(4)
          PJT(1,2*NSG-1)=0.D0
          PJT(2,2*NSG-1)=0.D0
          PJT(3,2*NSG-1)=(P0P-PP-PT**2/PM)*0.5D0
          PJT(4,2*NSG-1)=(P0P-PP-PT**2/PM)*0.5D0
          PJT(5,2*NSG-1)=0.D0
          ICJT(1,2*NSG-1)=ICJT(1,2*NS-1)
          ICJT(2,2*NSG-1)=ICJT(2,2*NS-1)
          PJT(1,2*NSG)=PX
          PJT(2,2*NSG)=PY
          PJT(3,2*NSG)=(PT**2/PM-PM)*0.5D0
          PJT(4,2*NSG)=(PT**2/PM+PM)*0.5D0
          PJT(5,2*NSG)=0.D0
          ICJT(1,2*NSG)=ICPM(1)
          ICJT(2,2*NSG)=ICPM(2)
          DO 11 NF=1,NFLAV
            JC(NF,1)=JCPS(NF,1)
            JC(NF,2)=JCPS(NF,2)
11        CONTINUE
        ELSE
          XORSG(1,NSG)=XOR2(1)
          XORSG(2,NSG)=XOR2(2)
          XORSG(3,NSG)=XOR2(3)
          XORSG(4,NSG)=XOR2(4)
          PJT(1,2*NSG-1)=-PX
          PJT(2,2*NSG-1)=-PY
          PJT(3,2*NSG-1)=(PP-PT**2/PP)*0.5D0
          PJT(4,2*NSG-1)=(PP+PT**2/PP)*0.5D0
          PJT(5,2*NSG-1)=0.D0
          ICJT(1,2*NSG-1)=ICMP(1)
          ICJT(2,2*NSG-1)=ICMP(2)
          DO 12 NF=1,NFLAV
            JC(NF,1)=JCMS(NF,1)
            JC(NF,2)=JCMS(NF,2)
12        CONTINUE
          PJT(1,2*NSG)=0.D0
          PJT(2,2*NSG)=0.D0
          PJT(3,2*NSG)=-(P0M-PM-PT**2/PP)*0.5D0
          PJT(4,2*NSG)= (P0M-PM-PT**2/PP)*0.5D0
          PJT(5,2*NSG)=0.D0
          ICJT(1,2*NSG)=ICJT(1,2*NS)
          ICJT(2,2*NSG)=ICJT(2,2*NS)
        ENDIF
        ICX(1)=ICJT(1,2*NSG-1)
        ICX(2)=ICJT(2,2*NSG-1)
        ICY(1)=ICJT(1,2*NSG)
        ICY(2)=ICJT(2,2*NSG)
        CALL IDDECO(ICX,JCX)
        CALL IDDECO(ICY,JCY)
        DO 28 N=1,NFLAV
          DO 29 I=1,2
            K=JCX(N,I)+JCY(N,I)-JC(N,I)
            IF ( K .LE. 0 ) GOTO 28
            DO 19 L=1,K
              JX=JCX(N,I)
              JY=JCY(N,I)
              IF     ( JX .EQ. 0  .AND.  JY .GT. 0 ) THEN
                JCY(N,I)=JCY(N,I)-1
              ELSEIF ( JX .GT. 0  .AND.  JY .EQ. 0 ) THEN
                JCX(N,I)=JCX(N,I)-1
              ELSEIF ( JX .GT. 0  .AND.  JY .GT. 0 ) THEN
                IF ( RANGEN() .LT. 0.5 ) THEN
                  JCX(N,I)=JCX(N,I)-1
                ELSE
                  JCY(N,I)=JCY(N,I)-1
                ENDIF
              ELSE
                WRITE(IFCH,*)('*',LP=1,71)
                WRITE(IFCH,*)'*****  IC;  ',IC
                WRITE(IFCH,*)'*****  ICX: ',ICX
                WRITE(IFCH,*)'*****  ICY: ',ICY
                WRITE(IFCH,*)'*****  N,I,K,L: ',N,I,K,L
                WRITE(IFCH,*)'*****  JX,JY: ',JX,JY
                WRITE(IFCH,*)('*',LP=1,71)
                CALL UTSTOP('JAMBR1: ERROR DURING JET COMPACTIFICATN ')
              ENDIF
19          CONTINUE
29        CONTINUE
28      CONTINUE
        CALL IDENCO(JCX,ICX,IRETEN)
        IF ( IRETEN .EQ. 1 ) THEN
          CALL UTSTOP('JAMBR1: IDENCO(JCX... RET.CODE=1        ')
        ENDIF
        CALL IDENCO(JCY,ICY,IRETEN)
        IF ( IRETEN .EQ. 1 ) THEN
          CALL UTSTOP('JAMBR1: IDENCO(JCY... RET.CODE=1        ')
        ENDIF
        ICJT(1,2*NSG-1)=ICX(1)
        ICJT(2,2*NSG-1)=ICX(2)
        ICJT(1,2*NSG)  =ICY(1)
        ICJT(2,2*NSG)  =ICY(2)
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,108)NSG,NS
     *       ,(ICJT(J,2*NSG-1),J=1,2),(SNGL(PJT(J,2*NSG-1)),J=1,5)
108       FORMAT(1X,I5,I4,I9,I7,1P,5E10.2)
          WRITE(IFCH,108)NSG,NS
     *       ,(ICJT(J,2*NSG  ),J=1,2),(SNGL(PJT(J,2*NSG  )),J=1,5)
        ENDIF
        PSG(1,NSG)=PJT(1,2*NSG-1)+PJT(1,2*NSG)
        PSG(2,NSG)=PJT(2,2*NSG-1)+PJT(2,2*NSG)
        PSG(3,NSG)=PJT(3,2*NSG-1)+PJT(3,2*NSG)
        PSG(4,NSG)=PJT(4,2*NSG-1)+PJT(4,2*NSG)
        PSG(5,NSG)=
     *     SQRT(PSG(4,NSG)**2-PSG(3,NSG)**2-PSG(2,NSG)**2-PSG(1,NSG)**2)
        IF ( ISH .GE. 93 ) WRITE(IFCH,110)NSG,NS,(PSG(J,NSG),J=1,5)
110     FORMAT(1X,I5,I4,16X,1P,5E10.2)
        CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *             ,PSG(1,NSG),PSG(2,NSG),PSG(3,NSG))
        CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *          ,PSG(1,NSG),PSG(2,NSG),PSG(3,NSG),PSG(4,NSG))
        IORSG(NSG)=NS
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,102)NSG,IORSG(NSG),(PSG(J,NSG),J=1,5)
102       FORMAT(1X,I5,I4,16X,1P,5E10.2)
          WRITE(IFCH,*)' '
        ENDIF

C-C   DO 56 I=1,2
C-C     ICPZ(I)=ICJT(I,2*NSG-1)
C-C56 ICMZ(I)=ICJT(I,2*NSG)
C-C   CALL IDDECO(ICPZ,JCPZ)
C-C   CALL IDDECO(ICMZ,JCMZ)
C-C   DO 57 NF=1,NFLAV
C-C     JCZ(NF,1)=JCPZ(NF,1)+JCMZ(NF,1)
C-C57 JCZ(NF,2)=JCPZ(NF,2)+JCMZ(NF,2)
C-C   CALL IDENCO(JCZ,ICZ,IRETEN)
C-C   IDZ=IDTRA(ICZ,0,0,3)
C-C   AMZ=PSG(5,NSG)
C-C   CALL IDRES(IDZ,AMZ,IDRZ,IADJ)
C-C   IF ( IDRZ.EQ.110 .AND. RANGEN().LT.0.5 ) GOTO 9993

9998  CONTINUE
      XBKPTL(1,NPTLC-NPTL2)=XBK(3)
      XBKPTL(2,NPTLC-NPTL2)=XBK(4)
      ISPTL(NPTLC-NPTL2)=1

      RETURN
      END
C=======================================================================

      SUBROUTINE JAMBR2(NS,NSG,IRET)

C-----------------------------------------------------------------------
C  "SAMBA" (SCHOLTEN+WERNER, MAR 92)
C  BREAKS STRING NS ACCORDING TO A-M MODEL.
C  NS: POINTS TO THE CURRENT FRAGMENTING STRING.
C  NSG: POINTS TO THE LAST PRODUCED SUBSTRING (SUCC INCREASED)
C     INPUT:
C  PSG(,NS): MOMENTUM OF STRING NS IN PP-CMS
C  PJT(,2*NS-1),PJT(,2*NS): MOMENTUM OF END OF STRING NS IN PP-CMS
C     OUTPUT:
C  XORSG(,NSG): ORIGIN OF SUBSTRING NSG IN PP-CMS
C  PJT(,2*NSG-1),PJT(,2*NSG): MOMENTUM OF END OF SUBSTRING NSG IN PP-CMS
C  ICJT(,2*NSG-1),ICJT(,2*NSG): IC-CODE  OF END OF SUBSTRG NSG IN PP-CMS
C  PSG(,NSG): MOMENTUM OF SUBSTRING NSG IN PP-CMS
C  XBKPTL(,NPTLC-NPTL2): BREAKPOINT OF NS IN PP-CMS
C  ISPTL(,NPTLC-NPTL2): 1 IF STRING NS BREAKS
C-----------------------------------------------------------------------
      PARAMETER (MXPC=500)
      PARAMETER (MXSG=500)
      PARAMETER (MXJT=2*MXSG)
      PARAMETER (NBRM=99)
      PARAMETER (NFLAV=6)
      PARAMETER (NPTF=129)
      DOUBLE PRECISION XBKPTL
      COMMON /CDELRE/  DELRER
      COMMON /CJAMBR/  NPTLC,NPTL2
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPC/     XBKPTL(2,MXPC),ISPTL(MXPC)
      COMMON /CPTF/    FPTFS,FPTFSS,FPTFU,FPTFUS,FPTFUU
     *                ,QPTFS(NPTF),QPTFSS(NPTF),QPTFU(NPTF),QPTFUS(NPTF)
     *                ,QPTFUU(NPTF),XPTF(NPTF)
      DOUBLE PRECISION PJT,PSG,ROTSG,XORSG
      COMMON /CSG/     PJT(5,MXJT),PSG(5,MXSG),ROTSG(3,MXSG)
     *                ,XORSG(4,MXSG)
     *                ,ICJT(2,MXJT),IORSG(MXSG),ISG(MXSG)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /QUARKM/  SMAS,SSMAS,USMAS,UUMAS

      DOUBLE PRECISION A,AMAM,AMAX,ANEW,DA,DALFA,DATWID
     *                ,DAUXIL,DR,DY,DYST,DYT,DY1,DY2,PHI,PM,PP,PT
     *                ,PW(4),PX,PY,P0M,P0P,RA,RMX,RMY,XBK(4)
     *                ,XBR,XOR2(4),XP,X2,YBR,YP,YTD,YNEW
      REAL     PXBRAK(100),PYBRAK(100),XBREAK(100),YBREAK(100)
      INTEGER  ICM(2),ICMP(2),ICMPBR(2,100),ICMZ(2)
     *        ,ICP(2),ICPM(2),ICPMBR(2,100),ICPZ(2)
     *        ,ICZ(2),ITRD(100),JC(NFLAV,2)
     *        ,JCM(NFLAV,2),JCMP(NFLAV,2),JCMS(NFLAV,2),JCMZ(NFLAV,2)
     *        ,JCP(NFLAV,2),JCPM(NFLAV,2),JCPZ(NFLAV,2),JCZ(NFLAV,2)
      LOGICAL  LAST
      SAVE
C-----------------------------------------------------------------------
      NCNT3=0
      NSG0=NSG
 9993 NSG=NSG0

C  INITIALIZATION
C  --------------

      IRET=0
      ICP(1)=ICJT(1,2*NS-1)
      ICP(2)=ICJT(2,2*NS-1)
      ICM(1)=ICJT(1,2*NS)
      ICM(2)=ICJT(2,2*NS)
      ICPMBR(1,1)=ICM(1)
      ICPMBR(2,1)=ICM(2)

      CALL IDDECO(ICP,JCP)
      CALL IDDECO(ICM,JCM)
      NP=0
      DO 7 NF=1,NFLAV
        NP=NP+JCP(NF,1)-JCP(NF,2)
 7    CONTINUE
      AMMS=UTAMNX(JCP,JCM)
      AM=PSG(5,NS)
C
      J = 2*NS
      DAUXIL= PJT(3,J-1)+PJT(3,J)
      P0P=PJT(4,J-1)+PJT(4,J)+DAUXIL
      P0M=PJT(4,J-1)+PJT(4,J)-DAUXIL
      XT=P0M
      YT=P0P
C
      DALFA=DBLE(PAREA)
C  INITIALIZE
      XBREAK(1)=XT
      YBREAK(1)=0.
      PXBRAK(1)=0.
      PYBRAK(1)=0.

C  REDO
C  ----

      NCNT3=0
9994  NCNT3=NCNT3+1
      IBR=1
      XP=XT
      YP=0.D0
      YTD=YT
      IF ( NCNT3 .GT. 100 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JAMBR2')
          WRITE(IFCH,*)'*****  SPLIT KINEM NOT POSSIBLE.'
          WRITE(IFCH,112)
     *     (ICJT(J,2*NS-1),J=1,2),(ICJT(J,2*NS),J=1,2),PSG(5,NS)
112       FORMAT(1X,2I8,4X,2I8,4X,F7.2)
          CALL UTMSGF
        ENDIF
        IF ( NS .EQ. 1 ) THEN
          IRET=9999
          RETURN
        ENDIF
        IRET=9996
        RETURN
      ENDIF

C  SEARCH FOR BREAKPOINTS
C  ----------------------

 9    CONTINUE
      DO 17 NF=1,NFLAV
        JCPM(NF,1)=0
        JCPM(NF,2)=0
        JCMP(NF,1)=0
        JCMP(NF,2)=0
17    CONTINUE

C  ORDINARY STRINGS
C  ----------------
C     NQU=0
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'ORDINARY STRING FRAGMENTATION'
        WRITE(IFCH,*)' '
      ENDIF
      IF ( RANGEN() .LE. PDIQUA ) THEN
        NQU=2
      ELSE
        NQU=1
      ENDIF
      IF     ( MOD(NP+NQU,3) .EQ. 0 ) THEN
        II=1
      ELSEIF ( MOD(NP-NQU,3) .EQ. 0 ) THEN
        II=2
      ELSE
        CALL UTSTOP('JAMBR2: NO SINGLET CONSTRUCTION POSSIBLE')
      ENDIF
      IFLTT=0
      DO 18 N=1,NQU
        IFL=INT(RANGEN()/PUD)+1
        IFLTT=IFLTT*10+(IFL+1)/2
        JCPM(IFL,II)=JCPM(IFL,II)+1
        JCMP(IFL,3-II)=JCMP(IFL,3-II)+1
18    CONTINUE
      CALL IDENCO(JCPM,ICPM,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('JAMBR2: IDENCO(JCPM... RET.CODE=1       ')
      ENDIF
      CALL IDENCO(JCMP,ICMP,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('JAMBR2: IDENCO(JCMP... RET.CODE=1       ')
      ENDIF

C  CALCULATE PT OF STRING BREAKING
C  -------------------------------------
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,109)ICM(1),ICMP(1),ICPM(1),ICP(1)
     *                ,ICM(2),ICMP(2),ICPM(2),ICP(2)
109     FORMAT(1X,'FLAVORS:',2(I11,I7),/,9X,2(I11,I7),/)
        WRITE(IFCH,*)'IFLTT:',IFLTT
        WRITE(IFCH,*)' '
      ENDIF
      AMMP=UTAMNY(JCP,JCPM)
      AMMM=UTAMNY(JCMP,JCM)
      R = RANGEN()
      IF     ( IFLTT .EQ. 1 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFU ,R*QPTFU(NPTF))
C##       WRITE(IFCH,*)'JAMBR2:PT(OLD)=',PT
        ELSE
          RPT = R*FPTFU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(1.+RPT*2./AUXIL))
        ENDIF
      ELSEIF ( IFLTT .EQ. 2 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFS ,R*QPTFS(NPTF))
        ELSE
          RPT = R*FPTFS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(SMAS **2/AUXIL)+RPT*2./AUXIL)-SMAS **2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 11 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUU,R*QPTFUU(NPTF))
        ELSE
          RPT = R*FPTFUU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(UUMAS**2/AUXIL)+RPT*2./AUXIL)-UUMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 12  .OR.  IFLTT .EQ. 21 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUS,R*QPTFUS(NPTF))
        ELSE
          RPT = R*FPTFUS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(USMAS**2/AUXIL)+RPT*2./AUXIL)-USMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 22 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFSS,R*QPTFSS(NPTF))
        ELSE
          RPT = R*FPTFSS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(SSMAS**2/AUXIL)+RPT*2./AUXIL)-SSMAS**2)
        ENDIF
      ENDIF

      IF ( IBR .EQ. 1 ) THEN
        RMX=AMMM**2+PT**2
      ELSE
        RMX=0.001D0
      ENDIF
      RMY=AMMP**2+PT**2
      IF ( RMY .GT. XP*(YTD-YP) ) GOTO 8

C  SINGLE ARTRU-MENNESSIER BREAK (SAMB)
C  ------------------------------------
C  IN: XP,YP,YTD,RMX,RMY
C  IN: DALFA
C  OUT: XBR,YBR

      DYT=YTD-YP
C  Y-STEPPING RANGE
      AMAM=1.D0+(RMX-RMY)/(XP*DYT)
      A=(1.D0-4.D0*RMX/(AMAM*AMAM*XP*DYT))
      IF ( A .LT. 0.D0 ) GOTO 8
      A=SQRT(A)
      DY1=DYT*AMAM*(1.D0-A)*0.5D0
      DY2=DYT*AMAM*(1.D0+A)*0.5D0
      DYST=0.0001D0*(DYT-DY2)
ctp060203  3    CONTINUE
      DR=DBLE(1.-RANGEN())
      AMAX=XP*DYT*A*AMAM+RMX*LOG(DY1/DY2)+RMY*LOG((DYT-DY2)/(DYT-DY1))
      RA=-LOG(DR)/DALFA
      IF ( RA .GT. AMAX  .AND. IBR .GT. 1 ) GOTO 8
      RA=MOD(RA,AMAX)
      DY=DYT*SQRT(RMX)/(SQRT(RMX)+SQRT(RMY))
      DATWID=XP*(DY-DY1)+RMX*LOG(DY1/DY)+RMY*LOG((DYT-DY)/(DYT-DY1))
      ANEW=RA-DATWID
      ICOUNT=0
      IF ( ANEW .GT. 0.D0 ) GOTO 1
 2    CONTINUE
      DA=XP - RMX/DY - RMY/(DYT-DY)
      YNEW=DY+ANEW/DA
      DY=INT(YNEW/DYST)*DYST
      IF ( DY .LT. DY1 ) THEN
        WRITE(IFCH,*) 'DY,DY1',DY,DY1
        GOTO 4
      ENDIF
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))
      IF ( ANEW .LT. 0.D0 ) GOTO 2
      GOTO 4
 1    CONTINUE
      DA=XP - RMX/DY - RMY/(DYT-DY)
      YNEW=DY+ANEW/DA
      DY=(INT(YNEW/DYST)+1)*DYST
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))
      IF ( ANEW .GT. 0.D0 ) GOTO 1
      DY=DY-DYST
      ANEW=RA-XP*(DY-DY1)-RMX*LOG(DY1/DY)-RMY*LOG((DYT-DY)/(DYT-DY1))
 4    CONTINUE
      YBR=MIN( DY+RANGEN()*DYST, DY2 )
      X2=(XP-RMX/YBR)
      XBR=MIN( ANEW/DYST+RMY/(DYT-YBR), X2 )
C  BETTER: SOLVE FOR YBR FROM ANEW=0
C  FIND XBR FROM HOMOGENEOUS (NOT EXP) DISTR, X1<X<X2
C        X1=RMY/(DYT-YBR)
      YBR=YBR+YP
C
C  END SAMB
C
      IF ( IBR .GE. NBRM ) THEN
        CALL UTSTOP('JAMBR2: IBR>NBRM                        ')
      ENDIF
      IBR=IBR+1
      XBREAK(IBR)=XBR
      YBREAK(IBR)=YBR
      PHI=2.D0*PI*RANGEN()
      PXBRAK(IBR)=PT*COS(PHI)
      PYBRAK(IBR)=PT*SIN(PHI)
      ICPMBR(1,IBR)=ICPM(1)
      ICPMBR(2,IBR)=ICPM(2)
      ICMPBR(1,IBR)=ICMP(1)
      ICMPBR(2,IBR)=ICMP(2)
      XP=XBR
      YP=YBR
      AMLEFT=SQRT(YBR*(XT-XBR))
      AMRIGT=SQRT(XBR*(YT-YBR))
      GOTO 9
 8    CONTINUE

      IF ( IBR .EQ. 1 ) GOTO 9994

C  INITIALIZE TAIL END
C  -------------------
      XBREAK(IBR+1)=0.
      YBREAK(IBR+1)=YT
      PXBRAK(IBR+1)=0.
      PYBRAK(IBR+1)=0.
      ICMPBR(1,IBR+1)=ICP(1)
      ICMPBR(2,IBR+1)=ICP(2)

C  PRINT
C  -----
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'ICMPBR(1/2,)  X/YBREAK  PX/YBREAK:'
        DO 52 IB=1,IBR+1
          WRITE(IFCH,100)ICMPBR(1,IB),ICMPBR(2,IB)
     *            ,XBREAK(IB),YBREAK(IB),PXBRAK(IB),PYBRAK(IB)
100       FORMAT(1X,I10,I10,2E11.3,7X,2E11.3)
          WRITE(IFCH,100)ICPMBR(1,IB),ICPMBR(2,IB)
52      CONTINUE
        WRITE(IFCH,*)' '
      ENDIF

C  TIME ORDER BREAKPOINTS
C  ----------------------

      T1=0.
      DO 203 J=2,IBR
        T2=2.
        DO 20 I=2,IBR
C         T=XBREAK(I)/XT+YBREAK(I)/YT
          T=XBREAK(I)*YBREAK(I)/(XT*YT)
          IF ( T .LE. T1  .OR.  T .GT. T2 ) GOTO 20
          T2=T
          NT=I
20      CONTINUE
        T1=T2
        ITRD(J)=NT
203   CONTINUE

C  PRINT
C  -----
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'ITRD:'
        WRITE(IFCH,*)(ITRD(IB),IB=2,IBR)
        WRITE(IFCH,*)' '
      ENDIF

C  CHECK ACCEPTANCE CRITERIUM
C  --------------------------

      NBR=0
      DO 12 J=2,IBR
        I=ITRD(J)
        XBRI=XBREAK(I)
        YBRI=YBREAK(I)
C  FIND NEIGHBORING EARLIER BREAK POINTS
        IR=IBR+1
        IL=1
        DO 204 JN=2,J
          IN=ITRD(JN)
          IF ( IN .GT. I  .AND.  IN .LT. IR ) IR=IN
          IF ( IN .LT. I  .AND.  IN .GT. IL ) IL=IN
204     CONTINUE
        PML=(XBREAK(IL)-XBRI)
        PPL=(YBRI-YBREAK(IL))
        PMR=(XBRI-XBREAK(IR))
        PPR=(YBREAK(IR)-YBRI)
        PX=PXBRAK(I)
        PY=PYBRAK(I)
        PXL=PX-PXBRAK(IL)
        PYL=PY-PYBRAK(IL)
        AMMM=PML*PPL-(PXL*PXL+PYL*PYL)
        PXR=PXBRAK(IR)-PX
        PYR=PYBRAK(IR)-PY
        AMMP=PMR*PPR-(PXR*PXR+PYR*PYR)
        IF ( AMMP .LE. 0. ) GOTO 13
        IF ( AMMM .LE. 0. ) GOTO 13
        AMMP=SQRT(AMMP)
        AMMM=SQRT(AMMM)
        AMP=AMMP
        AMM=AMMM
        CALL UTRESM(ICMPBR(1,IR),ICMPBR(2,IR)
     *             ,ICPMBR(1,I),ICPMBR(2,I),AMP,IDPR,IADJP,IRETEN)
        CALL UTRESM(ICPMBR(1,IL),ICPMBR(2,IL)
     *             ,ICMPBR(1,I),ICMPBR(2,I),AMM,IDMR,IADJM,IRETEN)
        IF ( AMP .GT. AMMP ) GOTO 13
        IF ( AMM .GT. AMMM ) GOTO 13
        PT2=PX*PX+PY*PY
        PTL2=PXBRAK(IL)**2+PYBRAK(IL)**2
        PTR2=PXBRAK(IR)**2+PYBRAK(IR)**2
        D=(PML*PPL-PTL2-PT2)*0.5
        D=D*D-PT2*PTL2
        IF ( D .LE. 0. ) GOTO 13
        D=(PMR*PPR-PTR2-PT2)*0.5
        D=D*D-PT2*PTR2
        IF ( D .LE. 0. ) GOTO 13
        NBR=NBR+1
        GOTO 12
13      CONTINUE
        ITRD(J)=-1
        XBREAK(I)=-1.
        YBREAK(I)=-1.
12    CONTINUE
      IF ( NBR .EQ. 0 ) GOTO 9994

C  PRINT
C  -----
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'ICMPBR(1/2,)  X/YBREAK  PX/YBREAK:'
        DO 53 IB=1,IBR+1
          WRITE(IFCH,100)ICMPBR(1,IB),ICMPBR(2,IB)
     *         ,XBREAK(IB),YBREAK(IB),PXBRAK(IB),PYBRAK(IB)
          WRITE(IFCH,100)ICPMBR(1,IB),ICPMBR(2,IB)
53      CONTINUE
        WRITE(IFCH,*)' '
      ENDIF

C  BUILD NEW STRINGS
C  -----------------
      NBREAK=0
      TAUAVE=0.
      IL=1
      IN=IL
      IR=1
      LAST=.FALSE.
11    CONTINUE
      IR=IR+1
      IF ( IR .EQ. IBR+1 ) LAST=.TRUE.
      IF ( XBREAK(IR) .LT. 0. ) GOTO 11
      IF ( IN .EQ. 1 ) GOTO 14
ctp060203 15    CONTINUE
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'IL IN IR:  ',IL,IN,IR
        WRITE(IFCH,*)' '
      ENDIF
      PML=XBREAK(IL)-XBREAK(IN)
      PPL=YBREAK(IN)-YBREAK(IL)
      PMR=(XBREAK(IN)-XBREAK(IR))
      PPR=(YBREAK(IR)-YBREAK(IN))
      PX=PXBRAK(IN)
      PY=PYBRAK(IN)
      PXL=PX-PXBRAK(IL)
      PYL=PY-PYBRAK(IL)
      PTL2=(PXL*PXL+PYL*PYL)
      AMM=SQRT(PML*PPL-PTL2)
      PXR=PXBRAK(IR)-PX
      PYR=PYBRAK(IR)-PY
      PTR2=(PXR*PXR+PYR*PYR)
      AMP=SQRT(PMR*PPR-PTR2)
      CALL UTRESM(ICMPBR(1,IR),ICMPBR(2,IR)
     *           ,ICPMBR(1,IN),ICPMBR(2,IN),AMP,IDPR,IADJP,IRETEN)
      CALL UTRESM(ICPMBR(1,IL),ICPMBR(2,IL)
     *           ,ICMPBR(1,IN),ICMPBR(2,IN),AMM,IDMR,IADJM,IRETEN)
      AREA=PPL*PMR
      P0P=PPL+PPR
      P0M=PML+PMR
      TMM2=(PTL2+AMM**2)
      TMP2=(PTR2+AMP**2)
      IF(ISH.GE.90)THEN
        IF ( TMP2 .GT. PPR*PMR+1.E-4 ) THEN
          CALL UTMSG('JAMBR2')
          WRITE(IFCH,*)'*****  TMP*TMP.GT.PPR*PMR'
          WRITE(IFCH,*)'TMP*TMP PPR*PMR: ',TMP2,PPR*PMR
          CALL UTMSGF
        ENDIF
        IF ( TMM2. GT. PPL*PML+1.E-4 ) THEN
          CALL UTMSG('JAMBR2')
          WRITE(IFCH,*)'*****  TMM*TMM.GT.PPL*PML'
          WRITE(IFCH,*)'TMM*TMM PPL*PML: ',TMM2,PPL*PML
          CALL UTMSGF
        ENDIF
      ENDIF
      ETA=0.5*LOG(PPL/PMR)
C ------ ADJUST P NOT M
      IF     ( IADJP .EQ. 1  .AND.  IADJM .NE. 1 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)'ADJUST P NOT M'
        ETANEW=LOG((SQRT(AREA)*P0P)/(TMP2+AREA))
        IF ( ETANEW .LT. ETA-1.E-4 ) THEN
          CALL UTMSG('JAMBR2')
          WRITE(IFCH,*)'*****  ETANEW.LT.ETA'
          WRITE(IFCH,*)'ETANEW ETA: ',ETANEW,ETA
          CALL UTMSGF
          ETANEW=ETA
        ENDIF
        XBREAK(IN)=SQRT(AREA)*EXP(-ETANEW)+XBREAK(IR)
        YBREAK(IN)=SQRT(AREA)*EXP(ETANEW)+YBREAK(IL)
C ------ ADJUST M NOT P
      ELSEIF ( IADJM .EQ. 1  .AND.  IADJP .NE. 1 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)'ADJUST M NOT P'
        ETANEW=-LOG((SQRT(AREA)*P0M)/(TMM2+AREA))
        IF ( ETANEW .GT. ETA+1.E-4 ) THEN
          IF(ISH.GE.90)THEN
            CALL UTMSG('JAMBR2')
            WRITE(IFCH,*)'*****  ETANEW.GT.ETA'
            WRITE(IFCH,*)'ETANEW ETA: ',ETANEW,ETA
            CALL UTMSGF
          ENDIF
          ETANEW=ETA
        ENDIF
        XBREAK(IN)=SQRT(AREA)*EXP(-ETANEW)+XBREAK(IR)
        YBREAK(IN)=SQRT(AREA)*EXP(ETANEW)+YBREAK(IL)
C ------ ADJUST BOTH
      ELSEIF ( IADJP .EQ. 1  .AND.  IADJM .EQ. 1 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)'ADJUST BOTH'
        D=(P0P*P0M-TMP2-TMM2)**2-4*TMP2*TMM2
        IF ( D .LT. 0. ) THEN
          IF(ISH.GE.90)THEN
            CALL UTMSG('JAMBR2')
            WRITE(IFCH,*)'***** NEGATIVE D'
            WRITE(IFCH,*)'D: ',D
            CALL UTMSGF
          ENDIF
          D=0.
        ENDIF
        D=SQRT(D)
        T=P0P*P0M-TMP2+TMM2
        P1=(T+D)/(2.*P0M)
        P2=(T-D)/(2.*P0M)
        IF ( (P1+P2)*0.5 .GT. PPL ) THEN
          P=P2
        ELSE
          P=P1
        ENDIF
        XBREAK(IN)=P0M-TMM2/P+XBREAK(IR)
        YBREAK(IN)=P+YBREAK(IL)
      ENDIF
C  ------
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'X/YBREAK:',XBREAK(IN),YBREAK(IN)
        WRITE(IFCH,*)' '
      ENDIF

C  WRITE SG
C  --------
      NBREAK=NBREAK+1
      TAUAVE=TAUAVE+XBREAK(IN)*YBREAK(IN)
      MM=1
      IF ( LAST ) MM=2
      DO 9998 M=1,MM
        NSG=NSG+1
        IF ( NSG .GT. MXSG ) THEN
          CALL UTSTOP('JAMBR2: NSG>MXSG                        ')
        ENDIF
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)'NSG LAST:  ',NSG,LAST
        ISG(NSG)=NPTLC
        IF ( M .EQ. 2 ) THEN
          IL=IN
          IN=IR
        ENDIF
        XOR2(1)=0.D0
        XOR2(2)=0.D0
        XOR2(3)=(YBREAK(IL)-XBREAK(IN))*0.5D0/TENSN
        XOR2(4)=(YBREAK(IL)+XBREAK(IN))*0.5D0/TENSN
        IF ( ISH .GE. 92 ) WRITE(IFCH,119)XOR2
119     FORMAT(/,' ORIGIN: ',13X,4F7.2,/)
        CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *              ,XOR2(1),XOR2(2),XOR2(3))
        CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *              ,XOR2(1),XOR2(2),XOR2(3),XOR2(4))
        XOR2(1)=XOR2(1)+XORSG(1,NS)
        XORSG(1,NSG)=XOR2(1)
        XOR2(2)=XOR2(2)+XORSG(2,NS)
        XORSG(2,NSG)=XOR2(2)
        XOR2(3)=XOR2(3)+XORSG(3,NS)
        XORSG(3,NSG)=XOR2(3)
        XOR2(4)=XOR2(4)+XORSG(4,NS)
        XORSG(4,NSG)=XOR2(4)
C       IF ( ISH .GE. 92 ) WRITE(IFCH,117)XOR2
        PP=YBREAK(IN)-YBREAK(IL)
        PM=XBREAK(IL)-XBREAK(IN)
        PSG(1,NSG)=PXBRAK(IN)-PXBRAK(IL)
        PSG(2,NSG)=PYBRAK(IN)-PYBRAK(IL)
        PSG(3,NSG)=(PP-PM)*0.5D0
        PSG(4,NSG)=(PP+PM)*0.5D0
        SS=PSG(4,NSG)**2-PSG(3,NSG)**2-PSG(2,NSG)**2-PSG(1,NSG)**2
        PSG(5,NSG)=SQRT(SS)
        IF ( ISH .GE. 92 ) WRITE(IFCH,110)NSG,NS,(PSG(J,NSG),J=1,5)
110     FORMAT(1X,I5,I4,16X,1P,5E10.2)
        PW(1)=PSG(1,NSG)
        PW(2)=PSG(2,NSG)
        PW(3)=PSG(3,NSG)
        PW(4)=PSG(4,NSG)
        CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *              ,PSG(1,NSG),PSG(2,NSG),PSG(3,NSG))
        CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *              ,PSG(1,NSG),PSG(2,NSG),PSG(3,NSG),PSG(4,NSG))
        IORSG(NSG)=NS
        IF ( ISH .GE. 92 ) THEN
          WRITE(IFCH,102)NSG,IORSG(NSG),(PSG(J,NSG),J=1,5)
102       FORMAT(1X,I5,I4,16X,1P,5E10.2)
          WRITE(IFCH,*)' '
        ENDIF

C  WRITE JT
C  --------
        ICJT(1,2*NSG-1)=ICMPBR(1,IN)
        ICJT(2,2*NSG-1)=ICMPBR(2,IN)
        DO 41 NF=1,NFLAV
CDH ########  JCMS IST NICHT GESETZT!!!!
          JC(NF,1)=JCMS(NF,1)
          JC(NF,2)=JCMS(NF,2)
41      CONTINUE
        ICJT(1,2*NSG)=ICPMBR(1,IL)
        ICJT(2,2*NSG)=ICPMBR(2,IL)
        IF ( M .EQ. 2 ) THEN
          IDX=IADJP
        ELSE
          IDX=IADJM
        ENDIF
        IF ( IDX .EQ. 0 ) THEN
          PJT(1,2*NSG-1)=PXBRAK(IN)
          PJT(2,2*NSG-1)=PYBRAK(IN)
          PTJR2=PXBRAK(IN)**2+PYBRAK(IN)**2
          PTJL2=PXBRAK(IL)**2+PYBRAK(IL)**2
          AA=PM*PP-PTJL2-PTJR2
          AUXIL=SQRT(AA*AA*0.25-PTJL2*PTJR2)
          PRM=(AA*0.5+PTJR2-AUXIL)/PP
          PLP=(AA*0.5+PTJL2-AUXIL)/PM
          PJT(3,2*NSG-1)=(PP-PLP-PRM)*0.5D0
          PJT(4,2*NSG-1)=(PP-PLP+PRM)*0.5D0
          PJT(5,2*NSG-1)=0.D0
          PJT(1,2*NSG)=-PXBRAK(IL)
          PJT(2,2*NSG)=-PYBRAK(IL)
          PJT(3,2*NSG)=(PLP-PM+PRM)*0.5D0
          PJT(4,2*NSG)=(PLP+PM-PRM)*0.5D0
          PJT(5,2*NSG)=0.D0
          IF ( ISH .GE. 92 ) THEN
            WRITE(IFCH,108)NSG,NS
     *          ,(ICJT(J,2*NSG-1),J=1,2),(SNGL(PJT(J,2*NSG-1)),J=1,5)
108         FORMAT(2X,I3,I4,2X,2I7,5(E10.2))
            WRITE(IFCH,108)NSG,NS
     *          ,(ICJT(J,2*NSG  ),J=1,2),(SNGL(PJT(J,2*NSG  )),J=1,5)
          ENDIF
          ERR=    (PW(1)-PJT(1,2*NSG-1)-PJT(1,2*NSG))**2
          ERR=ERR+(PW(2)-PJT(2,2*NSG-1)-PJT(2,2*NSG))**2
          ERR=ERR+(PW(3)-PJT(3,2*NSG-1)-PJT(3,2*NSG))**2
          ERR=ERR+(PW(4)-PJT(4,2*NSG-1)-PJT(4,2*NSG))**2
        ENDIF
        ICPZ(1)=ICJT(1,2*NSG-1)
        ICPZ(2)=ICJT(2,2*NSG-1)
        ICMZ(1)=ICJT(1,2*NSG)
        ICMZ(2)=ICJT(2,2*NSG)
        CALL IDDECO(ICPZ,JCPZ)
        CALL IDDECO(ICMZ,JCMZ)
        DO 57 NF=1,NFLAV
          JCZ(NF,1)=JCPZ(NF,1)+JCMZ(NF,1)
          JCZ(NF,2)=JCPZ(NF,2)+JCMZ(NF,2)
57      CONTINUE
        CALL IDENCO(JCZ,ICZ,IRETEN)
        IDZ=IDTRA(ICZ,0,0,3)
        AMZ=PSG(5,NSG)
        CALL IDRES(IDZ,AMZ,IDRZ,IADJ)
        IF ( IDRZ .EQ. 110  .AND.  RANGEN() .LT. 0.5 ) GOTO 9993
9998  CONTINUE

14    CONTINUE
      IF ( LAST ) GOTO 10
      IL=IN
      IN=IR
      GOTO 11
10    CONTINUE
C
C  DETERMINE BREAK PNT
C  -----------------------------------
      TAUAVE=TAUAVE/NBREAK
      XBK(1)=0.D0
      XBK(2)=0.D0
      AUXIL1=SQRT(TAUAVE*YT/XT)
      AUXIL2=SQRT(TAUAVE*XT/YT)
      XBK(3)=(AUXIL1-AUXIL2)*0.5D0/TENSN
      XBK(4)=(AUXIL1+AUXIL2)*0.5D0/TENSN
      IF ( ISH .GE. 92 ) WRITE(IFCH,116) XBK
116   FORMAT(' BREAKING PNT: ',13X,4F7.2)
      CALL UTROT2(-1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *             ,XBK(1),XBK(2),XBK(3))
      CALL UTLOB2(-1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *             ,XBK(1),XBK(2),XBK(3),XBK(4))
      XBK(1)=XBK(1)+XORSG(1,NS)
      XBK(2)=XBK(2)+XORSG(2,NS)
      XBK(3)=XBK(3)+XORSG(3,NS)
      XBK(4)=XBK(4)+XORSG(4,NS)
      IF ( ISH .GE. 92 ) WRITE(IFCH,117)XBK
117   FORMAT(28X,4F7.2,/)
C
      XBKPTL(1,NPTLC-NPTL2)=XBK(3)
      XBKPTL(2,NPTLC-NPTL2)=XBK(4)
      ISPTL(NPTLC-NPTL2)=1
      RETURN
      END
C=======================================================================

      SUBROUTINE JAMFRA(JS,NEWEVT)

C-----------------------------------------------------------------------
C  FRAGMENTS STRING JS ACCORDING TO A-M MODEL.
C  VERSION MAR-92 (CALLS JAMBR1 OR JAMBR2)
C-----------------------------------------------------------------------
      PARAMETER (MXPC=500)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSG=500)
      PARAMETER (MXJT=2*MXSG)
      PARAMETER (MXSTR=3000)
      PARAMETER (NFLAV=6)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CDELRE/  DELRER
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CJAMBR/  NPTLC,NPTL2
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      DOUBLE PRECISION XBKPTL
      COMMON /CPC/     XBKPTL(2,MXPC),ISPTL(MXPC)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      DOUBLE PRECISION PJT,PSG,ROTSG,XORSG
      COMMON /CSG/     PJT(5,MXJT),PSG(5,MXSG),ROTSG(3,MXSG)
     *                ,XORSG(4,MXSG)
     *                ,ICJT(2,MXJT),IORSG(MXSG),ISG(MXSG)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION ARM(4),ARP(4),TOR,ZOR
      INTEGER          IC(2),ICJ1(2),ICJ2(2)
     *                ,ICM(2),ICP(2),ICUM(2),ICUP(2)
     *                ,JC(NFLAV,2) ,JCJ1(NFLAV,2),JCJ2(NFLAV,2)
     *                ,JCM(NFLAV,2),JCP(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      DELREC=0.600
      DELREX=0.050
      NEWEVT=0

      PSG(1,1)=PSTR(1,JS)
      PSG(2,1)=PSTR(2,JS)
      PSG(3,1)=PSTR(3,JS)
      PSG(4,1)=PSTR(4,JS)
      PSG(5,1)=PSTR(5,JS)
      ROTSG(1,1)=ROTSTR(1,JS)
      ROTSG(2,1)=ROTSTR(2,JS)
      ROTSG(3,1)=ROTSTR(3,JS)
      XORSG(1,1)=XORSTR(1,JS)
      XORSG(2,1)=XORSTR(2,JS)
      XORSG(3,1)=XORSTR(3,JS)
      XORSG(4,1)=XORSTR(4,JS)
      SQRTS=ABS(PSTR(5,JS))
      S=PSTR(5,JS)**2
      ISG(1)=IORSTR(JS)
      IORSG(1)=0
      PJT(1,1)=0.D0
      PJT(1,2)=0.D0
      PJT(2,1)=0.D0
      PJT(2,2)=0.D0
      PJT(3,1)=SQRTS*0.5D0
      PJT(3,2)=-SQRTS*0.5D0
      PJT(4,1)=SQRTS*0.5D0
      PJT(4,2)=SQRTS*0.5D0
      PJT(5,1)=0.D0
      PJT(5,2)=0.D0
      IF ( ROTSTR(3,JS) .LT. 0. ) THEN
        ICJT(1,1)=ICSTR(3,JS)
        ICJT(2,1)=ICSTR(4,JS)
        ICJT(1,2)=ICSTR(1,JS)
        ICJT(2,2)=ICSTR(2,JS)
      ELSE
        ICJT(1,1)=ICSTR(1,JS)
        ICJT(2,1)=ICSTR(2,JS)
        ICJT(1,2)=ICSTR(3,JS)
        ICJT(2,2)=ICSTR(4,JS)
      ENDIF
      ICJ1(1)=ICJT(1,1)
      ICJ1(2)=ICJT(2,1)
      ICJ2(1)=ICJT(1,2)
      ICJ2(2)=ICJT(2,2)
      CALL IDDECO(ICJ1,JCJ1)
      CALL IDDECO(ICJ2,JCJ2)
      NQJ1=0
      NQJ2=0
      DO 4 NF=1,NFLAV
        NQJ1=NQJ1+JCJ1(NF,1)-JCJ1(NF,2)
        NQJ2=NQJ2+JCJ2(NF,1)-JCJ2(NF,2)
4     CONTINUE
      IF ( NQJ1 .GE. 0 ) THEN
        IF ( NQJ2 .GE. 0 ) THEN
          ISI = 0
        ELSE
          ISI = 1
        ENDIF
      ELSE
        IF ( NQJ2 .GE. 0 ) THEN
          ISI = 2
        ELSE
          ISI = 3
        ENDIF
      ENDIF
      NQJSTR=ISI*1000000+ABS(NQJ1)*1000+ABS(NQJ2)

C  ENTRY STRING FRAGMENTATION
C  --------------------------
      NPTL2=NPTL
      NCNT2=0
9996  NCNT2=NCNT2+1
      NPTL=NPTL2
      NSG=1
      IF ( NCNT2 .GT. 1000 ) GOTO 1001

C  ENTRY SUBSTRING PROCESSING
C  --------------------------
      NS=0
9999  NS=NS+1
      IF ( NS .GT. NSG ) GOTO 9997
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'ENTRY SUBSTRING PROCESSING'
        WRITE(IFCH,*)' '
        WRITE(IFCH,107)NS
     *        ,(ICJT(J,2*NS-1 ),J=1,2),(SNGL(PJT(J,2*NS-1 )),J=1,5)
107     FORMAT(2X,I3,3X,3X,2I7,5(E10.2))
        WRITE(IFCH,107)NS
     *        ,(ICJT(J,2*NS   ),J=1,2),(SNGL(PJT(J,2*NS   )),J=1,5)
        WRITE(IFCH,101)NS,(PSG(J,NS),J=1,5)
101     FORMAT(2X,I3,23X,5(E10.2),/)
ctp060203 114     FORMAT(' LEFT BREAKING PNT: ',8X,4F7.2,/)
      ENDIF

C  DETERMINE ID AND MIN.MASS
C  -------------------------
      DELRER=DELREX+RANGEN()*(DELREM-DELREX)
      ICP(1)=ICJT(1,2*NS-1)
      ICP(2)=ICJT(2,2*NS-1)
      ICM(1)=ICJT(1,2*NS)
      ICM(2)=ICJT(2,2*NS)
      CALL IDCOMK(ICP)
      CALL IDCOMK(ICM)
      ICJT(1,2*NS-1)=ICP(1)
      ICJT(2,2*NS-1)=ICP(2)
      ICJT(1,2*NS)=ICM(1)
      ICJT(2,2*NS)=ICM(2)
      CALL IDDECO(ICP,JCP)
      CALL IDDECO(ICM,JCM)
      AMMS=UTAMNX(JCM,JCP)
      NUBAR=0
      DO 7 NF=1,NFLAV
        JC(NF,1)=JCP(NF,1)+JCM(NF,1)
        JC(NF,2)=JCP(NF,2)+JCM(NF,2)
        IF ( NF .GT. 4  .AND.  (JC(NF,1).NE.0 .OR. JC(NF,2).NE.0) ) THEN
          CALL UTSTOP('JAMFRA: FLAVOUR > 4                     ')
        ENDIF
        NUBAR=NUBAR+JC(NF,1)-JC(NF,2)
7     CONTINUE
      IF ( NS .EQ. 1 ) NUMBAR=NUBAR
      CALL IDENCO(JC,IC,IRETEN)
      CALL IDCOMK(IC)
      ID=IDTRA(IC,0,0,3)
      IDK=ID
      AM=PSG(5,NS)
      IREMN=0
      IF ( ICP(1)+ICP(2) .EQ. 0  .OR.  ICM(1)+ICM(2) .EQ. 0 ) IREMN=1
      IF ( AM .GT. AMMS+DELRER  .AND.  IREMN .EQ. 0 ) ID=0

C  MASS ADJUSTMENT => FRAGMENT AGAIN
C  ---------------------------------
      AMC=PSG(5,NS)
      CALL IDRES(ID,AMC,IDR,IADJ)
      IF ( ABS(AMC-PSG(5,NS)) .GT. 1.E-3 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JAMFRA')
          WRITE(IFCH,*)'*****  MASS CHANGED. FRAGMENT AGAIN.'
          WRITE(IFCH,*)'MASS BEFORE, AFTER: ',PSG(5,NS),AMC
          WRITE(IFCH,*)'IC,ID,IDR: ',IC,ID,IDR
          CALL UTMSGF
        ENDIF
        GOTO 9996
      ENDIF

C  IF MASS < MIN.MASS ==> AGAIN
C  ----------------------------
C-C   IF ( IDR.EQ.0.AND.AM.LT.AMMS-1.E-3 ) THEN
C-C     IF ( ISH .GE. 90 ) THEN
C-C       CALL UTMSG('JAMFRA')
C-C       WRITE(IFCH,*)'*****  MASS < MIN MASS. ',AM,AMMS
C-C       WRITE(IFCH,*)'IDK,ID,IDR,PSG(5,NS):'
C-C       WRITE(IFCH,*)IDK,ID,IDR,PSG(5,NS)
C-C       WRITE(IFCH,*)'P_JET:'
C-C       WRITE(IFCH,107)NS
C-C  *       ,(ICJT(J,2*NS-1 ),J=1,2),(SNGL(PJT(J,2*NS-1 )),J=1,5)
C-C       WRITE(IFCH,107)NS
C-C  *       ,(ICJT(J,2*NS   ),J=1,2),(SNGL(PJT(J,2*NS   )),J=1,5)
C-C       WRITE(IFCH,*)'P_STR:'
C-C       WRITE(IFCH,101)NS,(PSG(J,NS),J=1,5)
C-C       CALL UTMSGF
C-C     ENDIF
C-C     IF ( NS .EQ. 1 ) GOTO 1001
C-C     GOTO 9996
C-C   ENDIF

C  WRITE ON /CPTL/ (1)
C  -------------------
      IORI=ISG(NS)
      NPTLB=NPTL
ctp060203 9995  NPTL=NPTL+1
      NPTL=NPTL+1
      NPTLC=NPTL
      IF ( NPTL .GT. MXPTL ) THEN
        CALL UTSTOP('JAMFRA: NPTL>MXPTL                      ')
      ENDIF
      IF ( NPTL-NPTL2 .GT. MXPC ) THEN
        CALL UTSTOP('JAMFRA: NPTL-NPTL2>MXPC                 ')
      ENDIF

      PPTL(1,NPTL)=PSG(1,NS)
      PPTL(2,NPTL)=PSG(2,NS)
      PPTL(3,NPTL)=PSG(3,NS)
      PPTL(4,NPTL)=PSG(4,NS)
      PPTL(5,NPTL)=PSG(5,NS)
      XORPTL(1,NPTL)=XORSG(1,NS)
      XORPTL(2,NPTL)=XORSG(2,NS)
      XORPTL(3,NPTL)=XORSG(3,NS)
      XORPTL(4,NPTL)=XORSG(4,NS)
      XBKPTL(1,NPTL-NPTL2)=AINFIN
      XBKPTL(2,NPTL-NPTL2)=AINFIN
      TIVPTL(1,NPTL)=-AINFIN
      TIVPTL(2,NPTL)=AINFIN
      IFRPTL(1,NPTL)=0
      IFRPTL(2,NPTL)=0
      ICLPTL(NPTL)=0
      ISPTL(NPTL-NPTL2)=0
      IORPTL(NPTL)=IORI
      JORPTL(NPTL)=0
      IF     ( NPTL .GT. NPTL2+1 ) THEN
        NQJPTL(NPTL)=-NQJSTR
      ELSEIF ( NPTL .EQ. NPTL2+1 ) THEN
        NQJPTL(NPTL)= NQJSTR
      ENDIF
      IDPTL(NPTL)=IDR
      IF ( IDR .EQ. 0 ) THEN
        IF ( IC(1) .EQ. 0  .AND.  IC(2) .EQ. 0 ) THEN
          IDPTL(NPTL)=700000000
          CALL IDTRBI(JC,IBPTL(1,NPTL),IBPTL(2,NPTL)
     *                ,IBPTL(3,NPTL),IBPTL(4,NPTL))
        ELSE
          IB5=JC(5,1)*10000+JC(5,2)
          IB6=JC(6,1)*10000+JC(6,2)
          IF ( IB5 .NE. 0  .OR.  IB6 .NE. 0 ) THEN
            WRITE(IFCH,*)' '
            WRITE(IFCH,*)'***** ERROR IN JAMFRA: BOTTOM OR TOP QUARKS'
            WRITE(IFCH,*)'JC:'
            WRITE(IFCH,*)JC
            CALL UTSTOP('JAMFRA: BOTTOM OR TOP QUARKS            ')
          ENDIF
          IDPTL(NPTL)=800000000+IC(1)*100+IC(2)/100
        ENDIF
      ENDIF
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,115)NPTL,IORPTL(NPTL),IDPTL(NPTL)
     *                 ,(PPTL(J,NPTL),J=1,5)
115     FORMAT(1X,'/CPTL/',I6,I7,I10,5(E10.2))
        IF ( IDPTL(NPTL) .EQ. 700000000 )
     *                               WRITE(IFCH,*)(IBPTL(I,NPTL),I=1,4)
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)' '
      ENDIF
      IF ( IDR .NE. 0 ) GOTO 9999

      NQP=0
      NQM=0
      NAP=0
      NAM=0
      DO 23 NF=1,NFLAV
        NQP=NQP+JCP(NF,1)
        NQM=NQM+JCM(NF,1)
        NAP=NAP+JCP(NF,2)
        NAM=NAM+JCM(NF,2)
23    CONTINUE
      NP=NQP-NAP
      NM=NQM-NAM
      JP=NQP+NAP
      JM=NQM+NAM

C  QUARK-CLUSTER
C  -------------
      IF ( ICP(1)+ICP(2) .EQ. 0  .OR.  ICM(1)+ICM(2) .EQ. 0 ) GOTO 78
      IF ( AM .GT. AMMS+DELRER ) GOTO 77
      IF ( IDK .EQ. 0  .AND.  AM .GT. AMMS+DELREC ) GOTO 77
      IF ( IDK .NE. 0  .AND.  AM .GT. AMMS+DELREX ) GOTO 77
78    CONTINUE
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'QUARK-CLUSTER'
        WRITE(IFCH,*)' '
      ENDIF
      IF ( ICP(1)+ICP(2) .EQ. 0  .AND. ICM(1)+ICM(2) .EQ. 0 ) THEN
        CALL UTSTOP('JAMFRA: ZERO STRING.                    ')
      ENDIF
      GOTO 9999
77    CONTINUE

C  JET TRAFOS
C  ----------
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'JET TRAFOS'
        WRITE(IFCH,*)' '
      ENDIF
      ARP(1)=PJT(1,2*NS-1)
      ARP(2)=PJT(2,2*NS-1)
      ARP(3)=PJT(3,2*NS-1)
      ARP(4)=PJT(4,2*NS-1)
      ARM(1)=PJT(1,2*NS)
      ARM(2)=PJT(2,2*NS)
      ARM(3)=PJT(3,2*NS)
      ARM(4)=PJT(4,2*NS)
      IF ( ISH .GE. 90 ) CALL UTCHM(ARP,ARM,1)
      NSO=IORSG(NS)
      IF ( NSO .GT. 0 ) THEN
        CALL UTROT2(-1,ROTSG(1,NSO),ROTSG(2,NSO),ROTSG(3,NSO)
     *             ,ARP(1),ARP(2),ARP(3))
        CALL UTLOB2(-1,PSG(1,NSO),PSG(2,NSO),PSG(3,NSO),PSG(4,NSO)
     *              ,PSG(5,NSO),ARP(1),ARP(2),ARP(3),ARP(4))
        CALL UTROT2(-1,ROTSG(1,NSO),ROTSG(2,NSO),ROTSG(3,NSO)
     *               ,ARM(1),ARM(2),ARM(3))
        CALL UTLOB2(-1,PSG(1,NSO),PSG(2,NSO),PSG(3,NSO),PSG(4,NSO)
     *               ,PSG(5,NSO),ARM(1),ARM(2),ARM(3),ARM(4))
        IF ( ISH .GE. 90 ) CALL UTCHM(ARP,ARM,2)
        CALL UTLOB2(1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *              ,ARP(1),ARP(2),ARP(3),ARP(4))
        CALL UTLOB2(1,PSG(1,NS),PSG(2,NS),PSG(3,NS),PSG(4,NS),PSG(5,NS)
     *               ,ARM(1),ARM(2),ARM(3),ARM(4))
        IF ( ISH .GE. 90 ) CALL UTCHM(ARP,ARM,3)
        ROTSG(1,NS)=(ARP(1)-ARM(1))*0.5D0
        ROTSG(2,NS)=(ARP(2)-ARM(2))*0.5D0
        ROTSG(3,NS)=(ARP(3)-ARM(3))*0.5D0
        CALL UTROT2(1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *              ,ARP(1),ARP(2),ARP(3))
        CALL UTROT2(1,ROTSG(1,NS),ROTSG(2,NS),ROTSG(3,NS)
     *              ,ARM(1),ARM(2),ARM(3))
      ENDIF
      ICUP(1)=ICJT(1,2*NS-1)
      ICUP(2)=ICJT(2,2*NS-1)
      ICUM(1)=ICJT(1,2*NS)
      ICUM(2)=ICJT(2,2*NS)
      IF ( ARP(3) .LT. 0.D0 ) THEN
        ICJT(1,2*NS-1)=ICUM(1)
        ICJT(2,2*NS-1)=ICUM(2)
        PJT(1,2*NS-1)=ARM(1)
        PJT(2,2*NS-1)=ARM(2)
        PJT(3,2*NS-1)=ARM(3)
        PJT(4,2*NS-1)=ARM(4)
        ICJT(1,2*NS)=ICUP(1)
        ICJT(2,2*NS)=ICUP(2)
        PJT(1,2*NS)=ARP(1)
        PJT(2,2*NS)=ARP(2)
        PJT(3,2*NS)=ARP(3)
        PJT(4,2*NS)=ARP(4)
      ELSE
        ICJT(1,2*NS-1)=ICUP(1)
        ICJT(2,2*NS-1)=ICUP(2)
        PJT(1,2*NS-1)=ARP(1)
        PJT(2,2*NS-1)=ARP(2)
        PJT(3,2*NS-1)=ARP(3)
        PJT(4,2*NS-1)=ARP(4)
        ICJT(1,2*NS)=ICUM(1)
        ICJT(2,2*NS)=ICUM(2)
        PJT(1,2*NS)=ARM(1)
        PJT(2,2*NS)=ARM(2)
        PJT(3,2*NS)=ARM(3)
        PJT(4,2*NS)=ARM(4)
      ENDIF
      S=PSG(5,NS)**2
      IF ( ISH .GE. 90 ) THEN

        AUXIL=100.*ABS( SNGL(PJT(3,2*NS-1)-PJT(3,2*NS))-ABS(PSG(5,NS)) )
        IF ( AUXIL .GT. 1.  .AND.  AUXIL .GT. ABS(PSG(5,NS)) ) THEN
          CALL UTMSG('JAMFRA')
          WRITE(IFCH,*)'*****  PZ1-PZ2-SQRT(S) NONZERO'
          WRITE(IFCH,*)'VALUE:   '
     *                 ,SNGL(PJT(3,2*NS-1)-PJT(3,2*NS))-SQRT(S)
          WRITE(IFCH,*)'SQRT(S): ',SQRT(S)
          WRITE(IFCH,*)(SNGL(PJT(I,2*NS-1)),I=1,4)
          WRITE(IFCH,*)(SNGL(PJT(I,2*NS)),I=1,4)
          CALL UTMSGF
        ENDIF
        IF ( ISH .GE. 92 ) THEN
          WRITE(IFCH,103)NS,NSO
     *      ,(ICJT(J,2*NS-1 ),J=1,2),(SNGL(PJT(J,2*NS-1 )),J=1,5)
103       FORMAT(2X,I3,I3,3X,2I7,5(E10.2))
          WRITE(IFCH,104)
     *       (ICJT(J,2*NS   ),J=1,2),(SNGL(PJT(J,2*NS   )),J=1,5)
104       FORMAT(2X,9X,2I7,5(E10.2),/)
        ENDIF
      ENDIF

      IF     ( IOPBRK .EQ. 1 ) THEN
        CALL JAMBR1(NS,NSG,IRET)
      ELSEIF ( IOPBRK .EQ. 2 ) THEN
        CALL JAMBR2(NS,NSG,IRET)
      ENDIF
      IF ( IRET .EQ. 9996 ) GOTO 9996
      GOTO 9999

C  WRITE ON /CPTL/ (2)
C  -------------------
9997  CONTINUE
      IF ( NPTL .LE. NPTL2 ) THEN
        CALL UTSTOP('JAMFRA: NPTL<=NPTL2                     ')
      ENDIF
      IF ( ISH .GE. 91 ) WRITE(IFCH,*)' '
      DO 67 I=NPTL2+1,NPTL
        ISTPTL(I)=ISPTL(I-NPTL2)
        IO=IORPTL(I)
        IF ( IO .GT. 0 ) THEN
          IF ( IFRPTL(1,IO) .EQ. 0 ) IFRPTL(1,IO)=I
          IFRPTL(2,IO)=I
        ELSE
          IORPTL(I)=0
        ENDIF
        ZOR=XORSG(3,1)
        TOR=XORSG(4,1)
        R=RANGEN()
        TAURAN=-TAUREA*LOG(R)
        CALL UTTAIX(I,TAURAN,ZOR,TOR,ZIS,TIS)
        TIVPTL(1,I)=MAX(TIS,XORPTL(4,I))
        IF     ( ISTPTL(I) .NE. 0 ) THEN
          TAUBR=SQRT
     *          ((XBKPTL(2,I-NPTL2)-TOR)**2-(XBKPTL(1,I-NPTL2)-ZOR)**2)
          CALL UTTAIX(I,TAUBR,ZOR,TOR,ZIS,TIS)
          TIVPTL(2,I)=TIS
        ELSEIF ( IO .GT. 0 ) THEN
          TAUBR=SQRT
     *         ((XBKPTL(2,IO-NPTL2)-TOR)**2-(XBKPTL(1,IO-NPTL2)-ZOR)**2)
          CALL UTTAIX(I,TAUBR,XORSG(3,1),XORSG(4,1),ZISBR,TISBR)
          CALL IDTAU(IDPTL(I),PPTL(4,I),PPTL(5,I),TAUGM)
          TIVPTL(2,I)=TISBR+TAUGM
        ELSE
          CALL IDTAU(IDPTL(I),PPTL(4,I),PPTL(5,I),TAUGM)
          TIVPTL(2,I)=XORPTL(4,I)+TAUGM
        ENDIF
        IF ( ISH .GE. 92 ) WRITE(IFCH,120)I,IORPTL(I),IDPTL(I)
     *               ,(SNGL(XBKPTL(J,I-NPTL2)),J=1,2)
     *               ,(XORPTL(J,I),J=3,4),(TIVPTL(J,I),J=1,2)
120     FORMAT(1X,'/CPTL/',I6,I7,I10
     *         ,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2)
67    CONTINUE

1000  RETURN

1001  NEWEVT=1
      GOTO 1000

      END
C=======================================================================

      SUBROUTINE JCENTD

C----------------------------------------------------------------------
      PARAMETER (KPARX=15)
      PARAMETER (NQUAX=12)
      COMMON /CENTRO/  ENTRO(1+KPARX,1+NQUAX)
      REAL ENTROX(1+KPARX,1+NQUAX)

      DATA (ENTROX( 1,1+N),N=0,12)/
     *.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,
     *.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,.00000E+00,
     *.00000E+00/
      DATA (ENTROX( 4,1+N),N=0,12)/
     *.40254E+01,.59349E+01,.74855E+01,.87464E+01,.98267E+01,.10770E+02,
     *.11611E+02,.12368E+02,.13058E+02,.13691E+02,.14277E+02,.14822E+02,
     *.15332E+02/
      DATA (ENTROX( 7,1+N),N=0,12)/
     *.61944E+01,.89306E+01,.11114E+02,.12969E+02,.14595E+02,.16054E+02,
     *.17380E+02,.18599E+02,.19728E+02,.20780E+02,.21767E+02,.22696E+02,
     *.23574E+02/
      DATA (ENTROX(10,1+N),N=0,12)/
     *.68876E+01,.10183E+02,.12855E+02,.15159E+02,.17201E+02,.19046E+02,
     *.20734E+02,.22296E+02,.23751E+02,.25115E+02,.26400E+02,.27616E+02,
     *.28771E+02/
      DATA (ENTROX(13,1+N),N=0,12)/
     *.61944E+01,.99602E+01,.13068E+02,.15784E+02,.18208E+02,.20406E+02,
     *.22425E+02,.24295E+02,.26041E+02,.27681E+02,.29228E+02,.30696E+02,
     *.32092E+02/
      DATA (ENTROX(16,1+N),N=0,12)/
     *.40254E+01,.82375E+01,.11781E+02,.14923E+02,.17745E+02,.20311E+02,
     *.22669E+02,.24853E+02,.26890E+02,.28803E+02,.30609E+02,.32320E+02,
     *.33948E+02/
      SAVE
C----------------------------------------------------------------------
      KPH=KPARX
      NQH=NQUAX
      DO 2 N=1,1+NQUAX
        DO 1 K=1,1+KPARX
          ENTRO(K,N)=ENTROX(K,N)
1       CONTINUE
2     CONTINUE

      IF ( KPH .NE. 15 ) THEN
        CALL UTSTOP('ICENTD: INSUFFICIENT INITIALIZATION;   K')
      ENDIF
      IF ( NQH .NE. 12 ) THEN
        CALL UTSTOP('ICENTD: INSUFFICIENT INITIALIZATION;   N')
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE JCENTP

C----------------------------------------------------------------------
C  PLOTS ENTRO(,) AND FTN SJCENT
C----------------------------------------------------------------------
      PARAMETER (KPARX=15)
      PARAMETER (NQUAX=12)
      COMMON /CENTRO/  ENTRO(1+KPARX,1+NQUAX)
      COMMON /CJCENT/  IGX,NSYMX
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      XPLOT(101),YPLOT(101)
      INTEGER   IPLOT(5)
      CHARACTER TEXT*50
      DATA      IPLOT/0,0,0,1,1/
      SAVE
C----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 13 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .LT. 95 ) GOTO 1000

      IF ( IPLOT(1) .EQ. 1 ) THEN
        TEXT='TITLE IG=     NSYM=     K= 3,6,9,12,15$       '
        WRITE(TEXT(10:11),122)IGX
        WRITE(TEXT(20:21),122)NSYMX
122     FORMAT(I2)
        DO 54 KX=1,5
          K=3*KX
          DO 55 N=1,13
            XPLOT(N)=N-1
            YPLOT(N)=ENTRO(1+K,N)
55        CONTINUE
          CALL UTHIST(0.,12.,0.,40.,1+12,XPLOT,YPLOT,'LIN','LINLIN'
     *      ,'XAXIS ENERGY / OMEGA            $                 '
     *      ,'YAXIS ENTROPY                   $                 ',TEXT)
54      CONTINUE
      ENDIF

      IF ( IPLOT(2) .EQ. 1 ) THEN
        TEXT='TITLE IG=     NSYM=     K= 3,6,9,12,15$       '
        WRITE(TEXT(10:11),122)IGX
        WRITE(TEXT(20:21),122)NSYMX
        DO 64 KX=1,5
          K=3*KX
          DO 65 N=1,13
            XPLOT(N)=(N-1.)/FLOAT(K)
            YPLOT(N)=ENTRO(1+K,N)/K
65        CONTINUE
          CALL UTHIST(0.,5.,0.,6.,1+12,XPLOT,YPLOT,'LIN','LINLIN'
     *      ,'XAXIS ENERGY / OMEGA / PARTICLE $                 '
     *      ,'YAXIS ENTROPY / PARTICLE        $                 ',TEXT)
64      CONTINUE
      ENDIF

      IF ( IPLOT(3) .EQ. 1 ) THEN
        DO 62 KX=1,10
          K=3*KX
          DO 63 N=0,100
            X=N*0.05
            XPLOT(1+N)=X
            YPLOT(1+N)=SJCENT(K,K,K*X)/K
63        CONTINUE
          CALL UTHIST(0.,5.,0.,6.,101,XPLOT,YPLOT,'LIN','LINLIN'
     *        ,'XAXIS ENERGY / OMEGA / PARTICLE $                 '
     *        ,'YAXIS ENTROPY / PARTICLE        $                 '
     *        ,'TITLE INTER(EXTRA)POLATED       $                 ')
62      CONTINUE
      ENDIF

      IF ( IPLOT(4) .EQ. 1 ) THEN
        DO 58 KX=1,10
          K=3*KX
          DO 59 N=0,100
            X=N*0.5
            XPLOT(1+N)=X
            YPLOT(1+N)=SJCENT(K,K,X)
59        CONTINUE
          CALL UTHIST(0.,50.,0.,80.,101,XPLOT,YPLOT,'LIN','LINLIN'
     *      ,'XAXIS ENERGY / OMEGA            $                 '
     *      ,'YAXIS ENTROPY                   $                 '
     *      ,'TITLE INTER(EXTRA)POLATED  KU=K $                 ')
58      CONTINUE
      ENDIF

      IF ( IPLOT(5) .EQ. 1 ) THEN
        DO 48 KX=1,10
          K=3*KX
          DO 49 N=0,100
            X=N*0.5
            XPLOT(1+N)=X
            YPLOT(1+N)=SJCENT(K,0,X)
49        CONTINUE
          CALL UTHIST(0.,50.,0.,80.,101,XPLOT,YPLOT,'LIN','LINLIN'
     *      ,'XAXIS ENERGY / OMEGA            $                 '
     *      ,'YAXIS ENTROPY                   $                 '
     *      ,'TITLE INTER(EXTRA)POLATED  KU=0 $                 ')
48      CONTINUE
      ENDIF

1000  CONTINUE
      ISH=ISH0
      RETURN
      END
C=======================================================================

      SUBROUTINE JCENTR(NSYM,IG,IDI,INIT)

C----------------------------------------------------------------------
C  FILLS ARRAY DEGEN(1+K,1+N) CONTAINING
C    THE NUMBER OF MIXED SYMMETRIC (BOX OF WIDTH NSYM)
C    K-PARTICLE STATES WITH ENERGY N (UNITS OF OMEGA), BASED ON
C    OSCILLATOR (DIM: IDI) WAVEFUNCTIONS (ADDIT. DEGENERACY: IG) .
C  FILLS ARRAY ENTRO(1+K,1+N) = LOG( DEGEN(1+K,1+N) )  .
C  ENTRO(,) WRITTEN IN FORM OF DATA FOR JCENTD IF ISH=95
C    (ISHSUB=12.. SELECTS THIS).
C  INIT MUST BE SET 1 FOR THE FIRST CALL OF JCENTR
C    (INIT=1 CALLS UTTUCL, UTPART, UTOVEL)
C----------------------------------------------------------------------
      PARAMETER (IOVMAX=100)
      PARAMETER (ITAMAX=1000)
      PARAMETER (JOVMAX=100)
      PARAMETER (KPARX=15)
      PARAMETER (KTUMAX=100)
      PARAMETER (NQUAX=12)
      PARAMETER (JPAMAX=NQUAX*NQUAX*NQUAX)
      PARAMETER (KKPMAX=NQUAX*NQUAX*2)
      PARAMETER (NSYMAX=20)
      PARAMETER (NTUMAX=100)
      PARAMETER (NYMAX=1000)
      COMMON /CDEGEN/  DEGEN(1+KPARX,1+NQUAX)
      COMMON /CENTRO/  ENTRO(1+KPARX,1+NQUAX)
      COMMON /CJCENT/  IGX,NSYMX
      COMMON /COVEL/   OVEL(1+IOVMAX,1+JOVMAX)
      COMMON /CPARTA/  PARTA(NQUAX),IPART(NQUAX,JPAMAX)
      DOUBLE PRECISION TUCL
      COMMON /CTUCL/   TUCL(1+KTUMAX,1+NTUMAX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      YOFA(ITAMAX)
      INTEGER   IYO(NYMAX,KPARX),IYOL(ITAMAX),IYOM(ITAMAX),IYOO(ITAMAX)
     *         ,IYOR(NSYMAX),IYOX(KPARX),IYOZ(ITAMAX,NSYMAX,2)
     *         ,JYO(NYMAX),KKK(KPARX)
      CHARACTER CIGA*7
      CHARACTER*1 CYOX(1+NQUAX,KPARX),DELI(1+KPARX,1+NQUAX)
      SAVE
C----------------------------------------------------------------------
      IF ( INIT .EQ. 1 ) THEN
        CALL UTTUCL
        CALL UTPART
        CALL UTOVEL
      ENDIF

      WRITE(IFMT,*)'EXECUTE SR JCENTR ...'
      IF ( NSYM .GE. NSYMAX ) THEN
        CALL UTSTOP('JCENTR: DIMENSION NSYMAX TOO SMALL      ')
      ENDIF
      IF ( MOD(KPARX,NSYM) .NE. 0 ) THEN
        CALL UTSTOP('JCENTR: KPARX SHOULD BE MULTIPLE OF NSYM')
      ENDIF
      MMAX=1+(KPARX-1)/NSYM
      IZERO=0
      NSYMX=NSYM
      IGX=IG

C  LOOP OVER ENERGY
C  ----------------
      DO 6000 N=0,NQUAX
        WRITE(IFMT,*)'SR JCENTR: N=',N
        IF ( N .EQ. 0 ) THEN
          KKP=1
        ELSE
          KKP=-IPART(N,1)
          IF ( KKP .GT. KKPMAX ) THEN
            CALL UTSTOP('JCENTR: DIMENSION KKPMAX TOO SMALL      ')
          ENDIF
          IF ( -IPART(N,1) .NE. NINT(PARTA(N)) ) THEN
            CALL UTSTOP('JCENTR: # OF PARTITIONS DISAGREE        ')
          ENDIF
        ENDIF
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,103)('-',IC=1,79),N,KKP,('-',IC=1,79)
103       FORMAT(/,1X,79A1,/,7X,'N = ',I2,
     *                   '   --->   ',I4,' PARTITION(S)',/,1X,79A1,/)
          JJ=2
          DO 8 KK=1,KKP
            IF     ( N .EQ. 0 ) THEN
              WRITE(IFCH,102)KK,IZERO
            ELSEIF ( N .GT. 0 ) THEN
              LL=-IPART(N,JJ)
              WRITE(IFCH,102)KK,(IPART(N,I),I=JJ+1,JJ+MIN(20,LL))
102           FORMAT(2X,I2,'. PARTITION:',3X,20I3)
            ENDIF
            JJ=JJ+LL+1
 8        CONTINUE
        ENDIF

C  ZERO PARTICLES
C  --------------
        DEGEN(1,1+N)=1.
        ENTRO(1,1+N)=0.
        IF ( N .EQ. NQUAX ) THEN
          DELI(1,1+N)='/'
        ELSE
          DELI(1,1+N)=','
        ENDIF

C  LOOP OVER PARTICLE NUMBER
C  -------------------------
        DO 6001 KX=1,KPARX/3
          K=3*KX
          IF ( K .GT. KPARX ) THEN
            CALL UTSTOP('JCENTR: DIMENSION KPARX TOO SMALL       ')
          ENDIF
          DEGEN(1+K,1+N)=0.
          IF ( N .EQ. NQUAX ) THEN
            DELI(1+K,1+N)='/'
          ELSE
            DELI(1+K,1+N)=','
          ENDIF
          KKK(K)=KKP
          MMAXK=1+(K-1)/NSYM
          IF ( K .LT. KPARX ) THEN
            DO 5 L=K+1,KPARX
              DO 5 NY=1,NYMAX
                IYO(NY,L)=-1
 5          CONTINUE
          ENDIF
          IF ( ISH.GE.93 ) WRITE(IFCH,100)('-',IC=1,11),K,('-',IC=1,11)
100       FORMAT(/,1X,11A1,/,3X,'K = ',I2,/,1X,11A1,/)

C  LOOP OVER PARTITIONS
C  --------------------
          JJ=2
          DO 6002 KK=1,KKP
            IF     ( N .GT. 0 ) THEN
              LL=-IPART(N,JJ)
            ELSEIF ( N .EQ. 0 ) THEN
              LL=1
            ENDIF
            IF ( LL .GT. K ) THEN
              KKK(K)=KK-1
              GOTO 6003
            ENDIF
            IF ( ISH .GE. 93 ) THEN
              IF     ( N .EQ. 0 ) THEN
                WRITE(IFCH,102)KK,IZERO
              ELSEIF ( N .GT. 0 ) THEN
                WRITE(IFCH,102)KK,(IPART(N,I),I=JJ+1,JJ+MIN(20,LL))
              ENDIF
              WRITE(IFCH,*)' '
            ENDIF

C  CONSTRUCT YOUNG TABLEAUS
C  ------------------------
            NY1=1
            NY2=1
            JYO(NY1)=1
            DO 46 L=1,K
              IYO(NY1,L)=0
46          CONTINUE
            IF ( N .GT. 0 ) THEN
              IYO(NY1,K)=IPART(N,JJ+1)
              NY1=NY1-1
45            NY1=NY1+1
              IF ( JYO(NY1) .EQ. LL ) GOTO 51
              DO 43 LX=1,NSYM
                L=NSYM+1-LX
                DO 44 MX=1,MMAXK
                  M=MMAXK+1-MX
                  IF ( IYO(NY1,(M-1)*NSYM+L) .EQ. 0 ) THEN
                    IF ( L.EQ.NSYM .OR. (L.NE.NSYM.AND.
     *                           IYO(NY1,(M-1)*NSYM+L+1).NE.0) ) THEN
                      NY2=NY2+1
                      IF ( NY2 .GT. NYMAX ) THEN
                        CALL UTSTOP
     *                    ('JCENTR: DIMENSION NYMAX TOO SMALL       ')
                      ENDIF
                      JY=JYO(NY1)+1
                      JYO(NY2)=JY
                      DO 50 LP=1,K
                        IYO(NY2,LP)=IYO(NY1,LP)
50                    CONTINUE
                      IYO(NY2,(M-1)*NSYM+L)=IPART(N,JJ+JY)
                      IF ( NY2 .GT. 1 ) THEN
                        DO 47 NY3=1,NY2-1
                          IF ( JYO(NY3) .NE. JYO(NY2) ) GOTO 47
                          DO 48 LP=1,K
                            IF ( IYO(NY3,LP) .NE. IYO(NY2,LP) ) GOTO 47
48                        CONTINUE
                          NY2=NY2-1
                          GOTO 49
47                      CONTINUE
49                      CONTINUE
                      ENDIF
                    ENDIF
                    GOTO 43
                  ENDIF
44              CONTINUE
43            CONTINUE
              GOTO 45
51            CONTINUE
            ENDIF

C  LOOP OVER YOUNG TABLEAUS
C  ------------------------
            DO 6005 NY=NY1,NY2

              IF ( ISH .GE. 93  .AND.  NSYM .EQ. 3 ) THEN
                WRITE(IFCH,117)NY-NY1+1,((IYO(NY,(M-1)*NSYM+I),
     *                                               I=1,NSYM),M=1,1)
117             FORMAT(2X,I2,'. TABLEAU:',5X,3I2)
                IF ( MMAXK .GT. 1 ) WRITE(IFCH,110)
     *               ((IYO(NY,(M-1)*NSYM+I),I=1,NSYM),M=2,MMAXK)
110             FORMAT(19X,3I2)
                WRITE(IFCH,*)' '
              ENDIF

C  SELECT NN-SUBTABLEAUS
C  ---------------------
              DITAB=1
              DO 6004 NN=0,N
                IGA=IG*NINT(EXP(TUCL(1+IDI,1+NN)))

                DO 15 L=1,KPARX
                  IYOX(L)=-1
15              CONTINUE
                DO 11 M=1,MMAXK
                  MEMP=M-1
                  MS=(M-1)*NSYM
                  DO 11 I=1,NSYM
                    IF ( IYO(NY,MS+I) .EQ. NN ) GOTO 12
11              CONTINUE
12              CONTINUE
                MMAXKS=MMAXK-MEMP
                DO 13 L=MS+1,K
                  IF ( IYO(NY,L) .EQ. NN ) IYOX(L-MS)=IYO(NY,L)
13              CONTINUE
                IF ( NSYM .GT. 1 ) THEN
                  DO 17 I=1,NSYM-1
                    DO 16 M=1,MMAXK
                      MS=(M-1)*NSYM
                      IF ( IYOX(MS+1) .NE. -1 ) GOTO 18
16                  CONTINUE
                    DO 19 M=1,MMAXK
                      MS=(M-1)*NSYM
                      DO 20 L=1,NSYM-1
                        IYOX(MS+L)=IYOX(MS+L+1)
20                    CONTINUE
                      IYOX(MS+NSYM)=-1
19                  CONTINUE
17                CONTINUE
18                CONTINUE
                ENDIF
                DO 38 L=1,KPARX
                  CYOX(1+NN,L)=' '
                  IF ( IYOX(L) .GE. 0 ) CYOX(1+NN,L)='X'
38              CONTINUE

C  SKIP FOR EMPTY TABLEAUS
C  -----------------------
                DO 21 L=1,K
                  IF ( IYOX(L) .NE. -1 ) GOTO 22
21              CONTINUE
                GOTO 6004
22              CONTINUE

C  PRINT
C  -----
                IF ( ISH .GE. 93  .AND.  NSYM .EQ .3 ) THEN
                  IF ( IGA.GE.  1.AND.IGA.LT.  10 ) WRITE(CIGA,107)IGA
107               FORMAT('GL(',I1,')  ')
                  IF ( IGA.GT. 10.AND.IGA.LT. 100 ) WRITE(CIGA,108)IGA
108               FORMAT('GL(',I2,') ')
                  IF ( IGA.GT.100.AND.IGA.LT.1000 ) WRITE(CIGA,109)IGA
109               FORMAT('GL(',I3,')')
                  WRITE(IFCH,116)NN,((CYOX(1+NN,(M-1)*NSYM+I)
     *                                        ,I=1,NSYM),M=1,1),CIGA
116               FORMAT(2X,I2,'-SUBTABLEAU:',3X,3(1X,A1),3X,A7)
                  IF ( MMAXKS .GT. 1 ) THEN
                    DO 39 M=2,MMAXKS
                      WRITE(IFCH,106)(CYOX(1+NN,(M-1)*NSYM+I),I=1,NSYM)
106                   FORMAT(19X,3(1X,A1))
39                  CONTINUE
                  ENDIF
                  IF ( ISH .GE. 94 ) WRITE(IFCH,*)' '
                ENDIF

C  REDUCTIONS OF TABLEAUS 1+IYOX()
C  -------------------------------
                DISUTA=0.
                N1=1
                N3=1
                N300=1
                DO 24 I=1,NSYM
                  MA=0
                  ME=0
                  DO 25 M=1,MMAXK
                    IF ( 1+IYOX((M-1)*NSYM+I) .LE. 0 ) GOTO 25
                    IF ( MA .EQ. 0 ) MA=M
                    ME=M
25                CONTINUE
                  IYOZ(N1,I,1)=MA
                  IYOZ(N1,I,2)=ME
24              CONTINUE
                IYOL(N1)=0
                YOFA(N1)=1.
                IYOO(N1)=0
                IYOM(N1)=1
                N1=N1-1
5000            N1=N1+1
                N30=N3
                IF ( N1.GT.1 .AND. IYOL(N1).GT.IYOL(N1-1) ) N300=N3
                N2=N1
                DO 26 I=1,NSYM
                  IYOR(I)=MAX(0,IYOZ(N1,I,2)-1)
26              CONTINUE
                LEV=IYOL(N1)
                FAC=YOFA(N1)
                IF ( IGA-1.LE.IOVMAX .AND. LEV.LE.JOVMAX
     *                                    .AND. LEV.LE.IGA-1 ) THEN
                  FACX=EXP(OVEL(IGA,1+LEV))
                  IF ( ABS(FACX-FAC) .GT. 1.E-5*FAC ) THEN
                    WRITE(IFCH,*)' '
                    WRITE(IFCH,*)'N=',N,'    K=',K,'    KK=',KK
     *                                                ,'    NY=',NY
                    WRITE(IFCH,*)'FACX=',FACX,'     FAC=',FAC
                    CALL UTSTOP
     *                     ('JCENTR: BINOMIALS DIFFER                ')
                  ENDIF
                ENDIF
                IF ( LEV .GT. IGA-1 ) GOTO 5003
                IHEIM=0
                DO 32 I=1,NSYM
                  IHEI=IYOZ(N1,I,2)-IYOZ(N1,I,1)
                  IF ( IHEI .GT. IHEIM ) IHEIM=IHEI
32              CONTINUE
                IF ( IHEIM .EQ. 0 ) DISUTA=DISUTA+FAC*IYOM(N1)
                IF ( ISH.GE.94 .AND. NSYM.EQ.3 ) THEN
                  WRITE(IFCH,112)N1,(IYOZ(N1,I,1),I=1,NSYM),IYOO(N1)
     *                ,LEV,FAC,IYOM(N1),IHEIM,DISUTA
     *                ,(IYOZ(N1,I,2),I=1,NSYM)
112               FORMAT(3X,I3,2X,3I2,2X,'ORI:',I3,3X
     *                ,'LEV:',I2,3X,'FAC:',F8.1
     *                ,3X,'MUL:',I2,3X,'HEI:',I2,3X,'SUM:',F8.1,/,8X
     *                ,3I2,/)
                ENDIF
                N2=N2-1
5001            N2=N2+1
                DO 27 I=1,NSYM
                  IF ( IYOZ(N2,I,2)-1 .LT. IYOR(I) ) GOTO 27
                  IF ( I.LT.NSYM  .AND.
     *                         IYOZ(N2,I,2)-1.LT.IYOZ(N2,I+1,2) ) GOTO27
                  N3=N3+1
                  IF ( N3 .GT. ITAMAX ) THEN
                    CALL UTSTOP
     *                      ('JCENTR: DIMENSION ITAMAX TOO SMALL      ')
                  ENDIF
                  IYOL(N3)=IYOL(N1)+1
                  YOFA(N3)=( YOFA(N1)*(IGA-IYOL(N3)) )/IYOL(N3)
                  IYOO(N3)=N1
                  IYOM(N3)=IYOM(N1)
                  DO 28 J=1,NSYM
                    IYOZ(N3,J,1)=IYOZ(N2,J,1)
                    IYOZ(N3,J,2)=IYOZ(N2,J,2)
28                CONTINUE
                  IYOZ(N3,I,2)=IYOZ(N3,I,2)-1
                  IF ( IYOZ(N3,I,2) .LT. IYOZ(N3,I,1) ) THEN
                    IYOZ(N3,I,1)=0
                    IYOZ(N3,I,2)=0
                  ENDIF
                  IF ( N30+1 .LT. N3 ) THEN
                    DO 30 NCH=N30+1,N3-1
                      DO 31 IJ=1,NSYM
                        IF ( IYOZ(NCH,IJ,1).NE.IYOZ(N3,IJ,1) .OR.
     *                      IYOZ(NCH,IJ,2).NE.IYOZ(N3,IJ,2) ) GOTO 30
31                    CONTINUE
                      N3=N3-1
                      GOTO 27
30                  CONTINUE
                  ENDIF
                  IF ( ISH .GE. 95  .AND.  NSYM .EQ. 3 ) THEN
                    WRITE(IFCH,113)N1,(IYOZ(N1,IJ,1),IJ=1,NSYM),N2
     *           ,(IYOZ(N2,IJ,1),IJ=1,NSYM),N3,(IYOZ(N3,IJ,1),IJ=1,NSYM)
     *           ,(IYOZ(N1,IJ,2),IJ=1,NSYM)
     *           ,(IYOZ(N2,IJ,2),IJ=1,NSYM),(IYOZ(N3,IJ,2),IJ=1,NSYM)
113                 FORMAT(3X,'N1: ',I2,3X,3I2,4X,'N2: ',I2,3X,3I2,4X
     *                 ,'N3: ',I2,3X,3I2,/,12X,3I2,13X,3I2,13X,3I2,/)
                  ENDIF
27              CONTINUE
                IF ( N2 .EQ. N1 ) N2=N30
                IF ( N2 .LT. N3 ) GOTO 5001
                IF ( N30 .LT. N3  .AND.  N300 .LT. N30 ) THEN
                  N3S=N3
                  N3M=N3
                  N3=N30
                  DO 33 NCH3=N30+1,N3S
                    N3=N3+1
                    DO 34 NCH=N300+1,N30
                      DO 35 IJ=1,NSYM
                        IF ( IYOZ(NCH,IJ,1).NE.IYOZ(N3,IJ,1) .OR.
     *                       IYOZ(NCH,IJ,2).NE.IYOZ(N3,IJ,2) ) GOTO 34
35                    CONTINUE
                      IYOM(NCH)=IYOM(NCH)+IYOM(N3)
                      N3=N3-1
                      N3M=N3M-1
                      IF ( N3+1 .LE. N3M ) THEN
                        DO 36 NM=N3+1,N3M
                          IYOL(NM)=IYOL(NM+1)
                          YOFA(NM)=YOFA(NM+1)
                          IYOO(NM)=IYOO(NM+1)
                          IYOM(NM)=IYOM(NM+1)
                          DO 37 J=1,NSYM
                            IYOZ(NM,J,1)=IYOZ(NM+1,J,1)
                            IYOZ(NM,J,2)=IYOZ(NM+1,J,2)
37                        CONTINUE
36                      CONTINUE
                      ENDIF
                      GOTO 33
34                  CONTINUE
33                CONTINUE
                ENDIF
                IF ( N1 .LT. N3 ) GOTO 5000
5003            CONTINUE

                DIALT=-99999.
                IF ( IYOZ(1,1,1) .EQ. 1 ) THEN
                  DIALT=1.
                  DO 40 I=1,NSYM
                    IF     ( IYOZ(1,I,1) .GT. 1 ) THEN
                      CALL UTSTOP
     *                     ('JCENTR: IYOZ(,,1).GT.1                  ')
                    ELSEIF ( IYOZ(1,I,1) .EQ. 0 ) THEN
                      GOTO 40
                    ENDIF
                    DO 41 M=1,IYOZ(1,I,2)
                      HAK=IYOZ(1,I,2)-M+1
                      IF ( I .LT. NSYM ) THEN
                        DO 42 J=I+1,NSYM
                          IF ( IYOZ(1,J,2) .GE. M ) HAK=HAK+1
42                      CONTINUE
                      ENDIF
                      DIALT=DIALT*(IGA+I-M)/HAK
41                  CONTINUE
40                CONTINUE
                ENDIF

                IF ( DIALT.NE.-99999.  .AND.
     *                   ABS(DIALT-DISUTA).GT.1.E-5*DIALT ) THEN
                  WRITE(IFCH,*)' '
                  WRITE(IFCH,*)'N=',N,'    K=',K,'    KK=',KK
     *                                                   ,'    NY=',NY
                  WRITE(IFCH,*)'DISUTA=',DISUTA,'     DIALT=',DIALT
                  CALL UTSTOP
     *                   ('JCENTR: DIMENSIONS DIFFER               ')
                ENDIF

                IF ( ISH .GE. 93 ) WRITE(IFCH,118)DISUTA
118             FORMAT(28X,'SUBTAB.-DIMENSION:',F8.1,/)

                DITAB=DITAB*DISUTA

6004          CONTINUE

              IF ( ISH .GE. 93 ) WRITE(IFCH,119)DITAB
119           FORMAT(35X,'TABLEAU-DIMENSION:',F11.1,/)

              DEGEN(1+K,1+N)=DEGEN(1+K,1+N)+DITAB

6005        CONTINUE

            JJ=JJ+LL+1
6002      CONTINUE
6003      CONTINUE
          IF ( KKK(K) .EQ. 0 ) THEN
            CALL UTSTOP('JCENTR: NO ALLOWED PARTITION            ')
          ENDIF
          ENTRO(1+K,1+N)=-9999999.
          IF ( DEGEN(1+K,1+N) .GT. 0. )
     *                         ENTRO(1+K,1+N)=LOG(DEGEN(1+K,1+N))
          IF ( ISH .GE. 93 ) WRITE(IFCH,120)DEGEN(1+K,1+N)
120       FORMAT(49X,'DEGENERACY:',F11.1,/)

6001    CONTINUE
6000  CONTINUE

      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,121)('-',IC=1,79)
121     FORMAT(1X,79A1)
        WRITE(IFCH,*)'   DEGENERACY(K,N)         IG=',IG,'          '
     *                ,'NSYM=',NSYM
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   K:','   3','   6','   9','   12'
        WRITE(IFCH,121)('-',IC=1,79)
        DO 52 N=0,NQUAX
          WRITE(IFCH,*)N,(DEGEN(1+3*K,1+N),K=1,4)
52      CONTINUE
        WRITE(IFCH,*)' '
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   ENTROPY(K,N)            IG=',IG,'          '
     *                ,'NSYM=',NSYM
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   K:','   3','   6','   9','   12'
        WRITE(IFCH,121)('-',IC=1,79)
        DO 53 N=0,NQUAX
          WRITE(IFCH,*)N,(ENTRO(1+3*K,1+N),K=1,4)
53      CONTINUE
      ENDIF

      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 12 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 95 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   ENTROPY(K,N)            IG=',IG
     *                ,'          NSYM=',NSYM
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,125)KPARX
125     FORMAT(6X,'IF(KPARX.NE.',I2,')'/5X,
     *      56H*CALL UTSTOP('JCENTD: INSUFFICIENT INITIALIZATION;   K'))
        WRITE(IFCH,126)NQUAX
126     FORMAT(6X,'IF(NQUAX.NE.',I2,')'/5X,
     *      56H*CALL UTSTOP('JCENTD: INSUFFICIENT INITIALIZATION;   N'))
        DO 58 KX=0,KPARX/3
          K=3*KX
          WRITE(IFCH,123)K+1,NQUAX
123       FORMAT(6X,'DATA (ENTRO(',I2,',1+N),N=0,',I2,')/')
          WRITE(IFCH,124)(ENTRO(1+K,1+N),DELI(1+K,1+N), N=0,NQUAX)
124       FORMAT(10(5X,'*',6(E10.5,A1),/))
58      CONTINUE
      ENDIF
      ISH=ISH0

      RETURN
      END
C=======================================================================

      SUBROUTINE JCLUDE(IP,IRET)

C-----------------------------------------------------------------------
C  DECAYS CLUSTER IP FROM /CPTL/ .
C  REQIRES JCENTR OR JCENTD TO BE CALLED BEFORE.
C-----------------------------------------------------------------------
      PARAMETER (IOMAX=54)
      PARAMETER (IOMAXM=25)
      PARAMETER (MOXMAX=30)
      PARAMETER (MOMAX=MOXMAX*IOMAX)
      PARAMETER (MXPTL=70000)
      PARAMETER (NFLAV=6)
      COMMON /CENTEX/  ENTEXP
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CSCAL/   SCAL
      COMMON /CSJCGA/  AMEGAM,AMNULL,ASUHA(7),ENTRPY,NOPHA,NSUHA(7)
      COMMON /CTIMEL/  NTC
      COMMON /CUTINV/  LUTINV
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION PA(5),PE(5),PPT(5)
      REAL             OO(1+MOMAX),PA0(5),POL(IOMAX)
     *                ,QO(1+MOMAX),U(3),XO(1+MOMAX),YO(1+MOMAX)
      INTEGER          ICA(2),IDOL(IOMAX),IDOLIS(IOMAX),IFOLIS(IOMAX)
     *                ,JCA(NFLAV,2),JCA0(NFLAV,2),JCE(NFLAV,2)
     *                ,JCO(NFLAV,2),KO(1+MOMAX)
      DATA IDOLIS/
     *             110,  120,-120,  220, 130, -130, 230, -230, 330
     *           ,1120,-1120,1220,-1220,1130,-1130,1230,-1230,2130,-2130
     *           ,2230,-2230,1330,-1330,2330,-2330
     *           ,1111,-1111,2221,-2221,3331,-3331
     *            ,111,  121,-121,  221, 131, -131, 231, -231, 331
     *           ,1121,-1121,1221,-1221,1131,-1131,1231,-1231
     *           ,2231,-2231,1331,-1331,2331,-2331/
      DATA IFOLIS/9*1,16*2,6*4,9*3,14*4/
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 1 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)'   CLUSTER DECAY OF',IP,IDPTL(IP),PPTL(5,IP)
        WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)' '
      ENDIF
      DELPOX=.01
      NPTLB=NPTL
      LOOP=0
      IRET=0
      IPOXRA=0
      ISTPFL=0
      EBAMIN=1.5
      LOOPMX=20
      NOPHAX=20

C  ORIGINAL CLUSTER --> PA,JCA
C  ---------------------------

6010  CONTINUE
      PA(1)=PPTL(1,IP)
      PA(2)=PPTL(2,IP)
      PA(3)=PPTL(3,IP)
      PA(4)=PPTL(4,IP)
      PA(5)=PPTL(5,IP)
      CALL IDQUAC(IP,NDU,NDU,NDU,JCA)
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'INITIAL JCA:'
        WRITE(IFCH,*)JCA
      ENDIF
      CALL IDCOMJ(JCA)
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'INITIAL JCA AFTER COMPACTIFICATION:'
        WRITE(IFCH,*)JCA
      ENDIF

      PA0(1)=PA(1)
      PA0(2)=PA(2)
      PA0(3)=PA(3)
      PA0(4)=PA(4)
      PA0(5)=PA(5)
      NQN=0
      DO 2 NF=1,NFLAV
        NQN=NQN+JCA(NF,1)-JCA(NF,2)
        JCA0(NF,1)=JCA(NF,1)
        JCA0(NF,2)=JCA(NF,2)
 2    CONTINUE
      NBA=ABS(NQN)/3
      IF ( NBA .GT. 0 ) THEN
        EBA=PA(5)/NBA
      ELSE
        EBA=AINFIN
      ENDIF
      ISH00=ISH

C  INITIALIZATION FOR REDO
C  -----------------------

6001  LOOP=LOOP+1
      IF ( ISH00 .EQ. 90 ) THEN
        IF ( LOOP .EQ. LOOPMX ) THEN
          ISH=92
          WRITE(IFCH,117)('-',L=1,75),LOOP,ISH
117       FORMAT(/,/,1X,75A1,/,1X,I2,'. ATTEMPT TO DECAY THIS CLUSTER!'
     *          ,/,1X,'ISH SET TO: ',I2)
        ENDIF
      ENDIF
      IF ( LOOP .GT. LOOPMX ) THEN
        ISH=ISH00
        GOTO 1001
      ENDIF
      NOPHA=0
      ISMALL=0
      NPTL=NPTLB
      PA(1)=PA0(1)
      PA(2)=PA0(2)
      PA(3)=PA0(3)
      PA(4)=PA0(4)
      PA(5)=PA0(5)
      DO 3 NF=1,NFLAV
        JCA(NF,1)=JCA0(NF,1)
        JCA(NF,2)=JCA0(NF,2)
 3    CONTINUE

C  LOOP OVER SEQUENTIAL DECAYS
C  ---------------------------

6002  NPTL=NPTL+1
      IF ( NPTL .GT. MXPTL ) THEN
        CALL UTSTOP('JCLUDE: NPTL>MXPTL                      ')
      ENDIF
      ISJCA=0
      DO 210 NF=1,NFLAV
        ISJCA=ISJCA+ABS(JCA(NF,1))+ABS(JCA(NF,2))
210   CONTINUE
      IF ( ISJCA .EQ. 0 ) THEN
        JCA(1,1)=1
        JCA(1,2)=1
      ENDIF
      CALL IDCOMJ(JCA)
      CALL IDENCO(JCA,ICA,IRETEN)
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,100)(PA(K),K=1,5)
100     FORMAT(1X,'PA:',20X,5(E10.3,1X))
        WRITE(IFCH,116)JCA
116     FORMAT(1X,'JCA: ',6I5,/,6X,6I5)
      ENDIF

      KAU=JCA(1,1)-JCA(1,2)
      KAD=JCA(2,1)-JCA(2,2)
      KAS=JCA(3,1)-JCA(3,2)
      KAC=JCA(4,1)-JCA(4,2)
      AMICL=UTAMNU(KAU,KAD,KAS,KAC,5)

      IF ( NOPHA.EQ.NOPHAX .OR. NOPHA.GT.1.AND.PA(5).GT.2.*AMICL ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JCLUDE')
          IF     ( ISMALL .EQ. 0 ) THEN
            WRITE(IFCH,*)'*****  NO PHASE SPACE --> REDO DECAY'
          ELSEIF ( ISMALL .EQ. 1 ) THEN
            WRITE(IFCH,*)'*****  SMALL PHASE SP --> REDO DECAY'
          ELSE
            WRITE(IFCH,*)'***** SMALL/NO PHASE SPACE --> REDO DECAY'
            WRITE(IFCH,*)'ISMALL=',ISMALL
          ENDIF
          WRITE(IFCH,*)'NT=',NTC,'   LOOP=',LOOP,'   NOPHA=',NOPHA
          IF ( ISMALL .EQ. 1 )
     *            WRITE(IFCH,*)'YO_HIT=',YOHIT,'   YO_MAX=',YOMAX
          WRITE(IFCH,*)'M_MIN=',AMICL,'   M=',PA(5)
          WRITE(IFCH,100)(PA(K),K=1,5)
          WRITE(IFCH,116)JCA
          CALL UTMSGF
        ENDIF
        GOTO 6001
      ENDIF

      IF ( NOPHA .GT. 1 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JCLUDE')
          IF     ( ISMALL .EQ. 0 ) THEN
            WRITE(IFCH,*)'*****  NO PHASE SPACE --> INCRS MASS'
          ELSEIF ( ISMALL .EQ. 1 ) THEN
            WRITE(IFCH,*)'*****  SMALL PHASE SP --> INCRS MASS'
          ELSE
            WRITE(IFCH,*)'***** SMALL/NO PHASE SPACE --> INCRS MASS'
            WRITE(IFCH,*)'ISMALL=',ISMALL
          ENDIF
          WRITE(IFCH,*)'NT=',NTC,'   LOOP=',LOOP,'   NOPHA=',NOPHA
          IF ( ISMALL .EQ. 1 )
     *             WRITE(IFCH,*)'YO_HIT=',YOHIT,'   YO_MAX=',YOMAX
          WRITE(IFCH,*)'M_MIN=',AMICL,'   M=',PA(5)
          WRITE(IFCH,100)(PA(K),K=1,5)
          WRITE(IFCH,116)JCA
        ENDIF
        PA(5)=PA(5)*1.05
        PA(4)=SQRT(PA(1)**2+PA(2)**2+PA(3)**2+PA(5)**2)
        IF ( ISH .GE. 90 ) THEN
          WRITE(IFCH,100)(PA(K),K=1,5)
          CALL UTMSGF
        ENDIF
      ENDIF

      IDAR=0
      IF ( NOPHA .EQ. 0 ) THEN
        ICH=0
        IF ( IRETEN .EQ. 0 ) THEN
          IDA=IDTRA(ICA,0,0,3)
        ELSE
          IDA=0
        ENDIF
        AMA=PA(5)
        CALL IDRES(IDA,AMA,IDAR,IADJ)
        IF ( IDAR .NE. 0  .AND.  EBA .LT. EBAMIN ) THEN
          AMA=0.
          CALL IDRES(IDA,AMA,IDAR,IADJ)
          IF ( MOD(IDAR,10) .NE. 0  .AND. LOOP .LT. LOOPMX/2 ) GOTO 6001
        ENDIF
        IF ( IDAR .NE. IDPTL(IP) ) ICH=1
      ENDIF

      IF ( IDAR .NE. 0 ) THEN
        IF ( NPTL .GT. NPTLB+1  .OR.  ICH .EQ. 1 ) THEN
          IDPTL(NPTL)=IDAR
          PA(5)=AMA
          PA(4)=SQRT(AMA**2+PA(1)**2+PA(2)**2+PA(3)**2)
          PPTL(1,NPTL)=PA(1)
          PPTL(2,NPTL)=PA(2)
          PPTL(3,NPTL)=PA(3)
          PPTL(4,NPTL)=PA(4)
          PPTL(5,NPTL)=PA(5)
          IF ( ISH .GE. 92 )
     *            WRITE(IFCH,101)NPTL,IDPTL(NPTL),(PPTL(K,NPTL),K=1,5)
101       FORMAT(1X,'/CPTL/',I6,I11,2X,5(E10.3,1X))
        ELSE
          NPTL=NPTL-1
        ENDIF
        GOTO 7001
      ENDIF

      AMA=PA(5)

C  LOOP OVER HADRONS
C  -----------------

      MO=1
      PO=0.
      XO(1)=0.
      QO(1)=0.
      YO(1)=0.
      OO(1)=0.
      KO(1)=0
      IOM=0
      QOL=0.
      IF ( ISH .GE. 92 ) WRITE(IFCH,106)
106   FORMAT(' PARTIAL DECAY WIDTHS:')

      DO 6003 IO=1,IOMAX
        IF ( EBA .LT. EBAMIN  .AND.  IO .GT. IOMAXM ) GOTO 6003

        IDO=IDOLIS(IO)
        IOM=IOM+1
        POL(IOM)=PO
        IDOL(IOM)=IDO
        CALL IDMASS(IDO,AMO)
        IDPTL(NPTL)=IDO
        CALL IDQUAC(NPTL,NDU,NDU,NDU,JCO)
        DO 22 NF=1,NFLAV
          JCE(NF,1)=JCA(NF,1)-JCO(NF,1)
          JCE(NF,2)=JCA(NF,2)-JCO(NF,2)
          IF ( JCE(NF,1) .LT. 0 ) THEN
            JCE(NF,2)=JCE(NF,2)-JCE(NF,1)
            JCE(NF,1)=0
          ENDIF
          IF ( JCE(NF,2) .LT. 0 ) THEN
            JCE(NF,1)=JCE(NF,1)-JCE(NF,2)
            JCE(NF,2)=0
          ENDIF
22      CONTINUE
        DO 13 I=5,NFLAV
          IF ( JCE(I,1).NE.0 .OR. JCE(I,2).NE.0 .OR. ISTPFL.GT.0 ) THEN
            IF ( ISTPFL .EQ. 2 ) THEN
              CALL UTSTOP('JCLUDE: MORE THAN 4 FLAVOURS            ')
            ENDIF
            ISTPFL=ISTPFL+1
            ISH=93
            IF ( ISTPFL .EQ. 2 ) GOTO 6010
            GOTO 6001
          ENDIF
13      CONTINUE
        KEU=JCE(1,1)-JCE(1,2)
        KED=JCE(2,1)-JCE(2,2)
        KES=JCE(3,1)-JCE(3,2)
        KEC=JCE(4,1)-JCE(4,2)

C  LOOP OVER HADRON MOMENTA
C  ------------------------

        IF ( MO+MOXMAX .GT. 1+MOMAX ) THEN
          CALL UTSTOP('JCLUDE: DIMENSION MOMAX TOO SMALL       ')
        ENDIF
        POX=0.
        MOM=MO
        POM=PO
        QOM=QO(MO)
        GAMXM=0.
        POXM=0.
        GAMY=0.
        DELPO=DELPOX
        YOINT=0.

        DO 6004 MOX=1,MOXMAX

          MO=MO+1
          DELPO=DELPO*1.2
          POX=POX+DELPO
          PO=PO+DELPO
          XO(MO)=PO
          IF ( MOX .GT. 1 ) THEN
            KO(MO)=0
          ELSE
            KO(MO)=IDO
          ENDIF
          GAMX=SJCGAM(KEU,KED,KES,KEC,AMA,AMO,POX,MOX)
          YO(MO)=IFOLIS(IO)*DELPO*GAMX
          OO(MO)=ENTEXP
          YOINT=YOINT+IFOLIS(IO)*DELPO*(GAMX+GAMY)*0.5
          IF ( GAMX .EQ. 0. ) THEN
            I6005=1
          ELSE
            IF ( GAMX .LT. 1.E-2*GAMXM  .AND.  GAMXM .GT. 0.
     *        .AND.  MOX .GT. 1 ) THEN
              I6005=1
            ELSE
              I6005=0
            ENDIF
          ENDIF
          IF ( ISH.GE.93 .AND. (I6005.EQ.0 .OR. MOX.GT.1) ) THEN
            IF ( MOX .EQ. 1 )
     *         WRITE(IFCH,109)KEU,KED,KES,NSUHA,AMA,ASUHA,AMNULL,IDO,AMO
109         FORMAT(/,' U_D_S:',3I3,'  N:',7I5,11X,'  A:',F10.2,
     *             /,18X,'M:',7F5.2,'  M0:',F6.2,'  O:',I5,F5.2,
     *             /,'   MO     XO   POX      EO AMA-EO     ENTRPY     '
     *             ,'  GAMX         YO         OO ')
            WRITE(IFCH,105)MO,XO(MO),POX,SQRT(POX**2+AMO**2)
     *                    ,AMA-SQRT(POX**2+AMO**2),ENTRPY,GAMX
     *                    ,YO(MO),OO(MO)
105         FORMAT(1X,I5,1X,F6.2,F6.2,2X,F6.2,F6.2
     *             ,E12.3,E11.3,E11.3,E11.3)
          ENDIF
          IF ( I6005 .EQ. 1 ) GOTO 6005
          IF ( GAMX .GT. GAMXM ) THEN
            GAMXM=GAMX
            POXM=POX
          ENDIF
          GAMY=GAMX

6004    CONTINUE

        IF ( ISH .GE. 90  .AND.  IPOXRA .EQ. 0 ) THEN
          IPOXRA=1
          CALL UTMSG('JCLUDE')
          WRITE(IFCH,*)'*****  POX-RANGE TOO SMALL'
          WRITE(IFCH,100)(PA(K),K=1,5)
          WRITE(IFCH,*)'JCA:'
          WRITE(IFCH,*)JCA
          WRITE(IFCH,*)'POX= ',POX, '   GAMX= ',GAMX
          WRITE(IFCH,*)'POXM=',POXM,'   GAMXM=',GAMXM
          CALL UTMSGF
        ENDIF

6005    CONTINUE

        IF ( YOINT .EQ. 0. ) THEN
          IF ( ISH .GE. 94 ) THEN
            WRITE(IFCH,*)' '
            WRITE(IFCH,*)'NO PHASE SPACE FOR',IDO
          ENDIF
          MO=MOM
          PO=POM
          IOM=IOM-1
        ELSE
          IF ( ISH .GE. 92 ) WRITE(IFCH,107)IDO,POXM,YOINT,OO(MO)
107       FORMAT(' IDO,POXM,YOINT,OO:',I6,3X,F5.2,E12.3,E12.3)
        ENDIF

6003  CONTINUE
      IF ( ISH .GE. 92 ) WRITE(IFCH,*)' '

C  NO PHASE SPACE
C  --------------
      IF ( IOM .EQ. 0 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)'NO PHASE SPACE'
        NPTL=NPTL-1
        NOPHA=NOPHA+1
        GOTO 6002
      ENDIF

C  DETERMINE QO(M)
C  ---------------
      OOMAX=0.
      DO 211 M=1,MO
        IF ( OO(M) .GT. OOMAX ) OOMAX=OO(M)
  211 CONTINUE
      OOX=OOMAX-10.
      DO 212 M=1,MO
        IF ( OO(M) .LT. OOX ) THEN
          YO(M)=0.
        ELSE
          YO(M)=YO(M)*EXP(OO(M)-OOX)
        ENDIF
  212 CONTINUE
      QO(1)=0.
      DO 213 M=2,MO
        IF ( KO(M) .NE. 0 ) YOM=0.
        QO(M)=QO(M-1)+(YOM+YO(M))*0.5
        YOM=YO(M)
  213 CONTINUE
      YOMAX=0.
      DO 214 M=1,MO
        IF ( YO(M) .GT. YOMAX ) YOMAX=YO(M)
  214 CONTINUE

      IF ( ISH .GE. 93 ) THEN
        IDO=0
        XOX=0.
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'ACCUMULATED DECAY PROBABILITIES:'
        WRITE(IFCH,*)' '
        DO 215 M=1,MO
          IF ( KO(M) .NE. 0 ) THEN
            WRITE(IFCH,*)' '
            IDO=KO(M)
          ENDIF
          WRITE(IFCH,110)M,IDO,XO(M)-XOX,QO(M)
110       FORMAT(1X,'M,IDO,XO,QO: ',I6,I6,3X,F6.2,E12.3)
          IF ( M .LT. MO  .AND.  KO(M+1) .NE. 0 ) XOX=XO(M)
  215   CONTINUE
      ENDIF

C  SELECT RANDOMLY HADRON+MOMENTUM
C  -------------------------------

      LO=0
24    LO=LO+1
      IF ( LO .GT. 10 ) THEN
        IF ( ISH.GE.92 ) WRITE(IFCH,*)'SMALL PHASE SPACE:',YOHIT,YOMAX
        NPTL=NPTL-1
        NOPHA=NOPHA+1
        ISMALL=1
        GOTO 6002
      ENDIF

      XOS=UTINVT(MO,XO,QO,RANGEN()*QO(MO))
      MHIT=LUTINV
      DO 28 I=2,IOM
        IX=I-1
        IF ( XOS .LT. POL(I) ) GOTO 29
28    CONTINUE
      IX=IOM
29    CONTINUE
      POS=XOS-POL(IX)
      IDS=IDOL(IX)
      CALL IDMASS(IDS,AMS)
      IDPTL(NPTL)=IDS

      CALL IDQUAC(NPTL,NDU,NDU,NDU,JCO)
      DO 26 NF=1,NFLAV
        JCE(NF,1)=JCA(NF,1)-JCO(NF,1)
        JCE(NF,2)=JCA(NF,2)-JCO(NF,2)
        IF ( JCE(NF,1) .LT. 0 ) THEN
          JCE(NF,2)=JCE(NF,2)-JCE(NF,1)
          JCE(NF,1)=0
        ENDIF
        IF ( JCE(NF,2) .LT. 0 ) THEN
          JCE(NF,1)=JCE(NF,1)-JCE(NF,2)
          JCE(NF,2)=0
        ENDIF
26    CONTINUE
      DO 27 I=5,NFLAV
        IF ( JCE(I,1) .NE. 0  .OR. JCE(I,2) .NE. 0 ) THEN
          CALL UTSTOP('JCLUDE: FLAVOUR > 4                     ')
        ENDIF
27    CONTINUE

      KEU=JCE(1,1)-JCE(1,2)
      KED=JCE(2,1)-JCE(2,2)
      KES=JCE(3,1)-JCE(3,2)
      KEC=JCE(4,1)-JCE(4,2)
      GAMXHT=SJCGAM(KEU,KED,KES,KEC,AMA,AMS,POS,1)
      YOHIT=GAMXHT*(XO(MHIT+1)-XO(MHIT))*4.
      IF ( YOHIT .LT. 1.E-5  .AND.  YOHIT .LT. 1.E-5*YOMAX ) GOTO 24

      PPTL(5,NPTL)=AMS
      U(3)=2.*RANGEN()-1.
      PHI=2.*PI*RANGEN()
      U(1)=SQRT(1.-U(3)**2)*COS(PHI)
      U(2)=SQRT(1.-U(3)**2)*SIN(PHI)
      PPT(1)=POS*U(1)
      PPTL(1,NPTL)=PPT(1)
      PPT(2)=POS*U(2)
      PPTL(2,NPTL)=PPT(2)
      PPT(3)=POS*U(3)
      PPTL(3,NPTL)=PPT(3)
      PPT(4)=SQRT(PPTL(5,NPTL)**2+POS**2)
      PPTL(4,NPTL)=PPT(4)
      IF ( ISH .GE. 93 )
     *    WRITE(IFCH,101)NPTL,IDPTL(NPTL),(PPTL(K,NPTL),K=1,5)
      CALL UTLOB2(-1,PA(1),PA(2),PA(3),PA(4),PA(5)
     *        ,PPT(1),PPT(2),PPT(3),PPT(4))
      PPTL(1,NPTL)=PPT(1)
      PPTL(2,NPTL)=PPT(2)
      PPTL(3,NPTL)=PPT(3)
      PPTL(4,NPTL)=PPT(4)
      IF ( ISH .GE. 92 )
     *        WRITE(IFCH,101)NPTL,IDPTL(NPTL),(PPTL(K,NPTL),K=1,5)

      PE(5)=AMEGAM
      PE(1)=-POS*U(1)
      PE(2)=-POS*U(2)
      PE(3)=-POS*U(3)
      PE(4)=SQRT(PE(5)**2+POS**2)
      CALL UTLOB2(-1,PA(1),PA(2),PA(3),PA(4),PA(5)
     *            ,PE(1),PE(2),PE(3),PE(4))

      PA(1)=PE(1)
      PA(2)=PE(2)
      PA(3)=PE(3)
      PA(4)=PE(4)
      PA(5)=PE(5)
      DO 23 NF=1,NFLAV
        JCA(NF,1)=JCE(NF,1)
        JCA(NF,2)=JCE(NF,2)
23    CONTINUE
      NOPHA=0

      GOTO 6002

7001  CONTINUE

C  CHECK ENERGY CONSERVATION
C  -------------------------
      IF ( NPTL .LE. NPTLB ) GOTO 1000

      IFAIL=1
      IF ( NPTL .GT. NPTLB+1 ) THEN
        ISHRSC=ISH
        ISH=0
        CALL JRESCL(NPTLB+1,NPTL,PA0,IFAIL)
        ISH=ISHRSC
C-C     IF ( IFAIL.NE.0 .AND. ISH.GE.90 ) THEN
C-C       CALL UTMSG('JCLUDE')
C-C       WRITE(IFCH,*)'*****  IFAIL_JRESCL=',IFAIL
C-C       CALL UTMSGF
C-C     ENDIF
      ENDIF

C-C   DO 114 N=NPTLB+1,NPTL
C-C     P=PPTL(3,N)
C-C     E=PPTL(4,N)
C-C     Y=100.
C-C     DY=3.
C-C     IF ( E-P.NE.0. .AND. E+P.NE.0. ) Y=.5*LOG((E+P)/(E-P))+DY
C-C     IDA=ABS(IDPTL(N))
C-C     IF ( IDA.GT.1000 .AND. MOD(IDA,10).NE.0 .AND. Y.LT.10. ) THEN
C-C       WRITE(6,*)('-',K=1,69)
C-C       P=PPTL(3,IP)
C-C       E=PPTL(4,IP)
C-C       Y=100.
C-C       DY=3.
C-C       IF ( E-P.NE.0. .AND. E+P.NE.0. ) Y=.5*LOG((E+P)/(E-P))+DY
C-C       WRITE(6,115)IP,IDPTL(IP)
C-C  *         ,(PPTL(K,IP),K=3,5),Y
C-C       WRITE(6,*)'-------> '
C-C       DO 113 M=NPTLB+1,NPTL
C-C         P=PPTL(3,M)
C-C         E=PPTL(4,M)
C-C         Y=100.
C-C         DY=3.
C-C         IF ( E-P.NE.0. .AND. E+P.NE.0. ) Y=.5*LOG((E+P)/(E-P))+DY
C-C         WRITE(6,115)M,IDPTL(M)
C-C  *            ,(PPTL(K,M),K=3,5),Y
C113      CONTINUE
C-C     ENDIF
C114  CONTINUE
C-C 115  FORMAT(1X,'/CPTL/',I6,I10
C-C     *,3(E10.2),2(E10.2))

      IF ( ISH .GE. 0 ) THEN
        PX=0.
        PY=0.
        PZ=0.
        E=0.
        DO 10 N=NPTLB+1,NPTL
          PX=PX+PPTL(1,N)
          PY=PY+PPTL(2,N)
          PZ=PZ+PPTL(3,N)
          E=E+  PPTL(4,N)
10      CONTINUE
      ENDIF

      IF ( ISH .GE. 90 ) THEN
        IF ( ISH .GE. 92 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'CHECK ENERGY CONSERVATION'
          DO 11 N=NPTLB+1,NPTL
            WRITE(IFCH,101)N,IDPTL(N),(PPTL(K,N),K=1,5)
11        CONTINUE
          AM=SQRT(E**2-PX**2-PY**2-PZ**2)
          WRITE(IFCH,103)PX,PY,PZ,E,AM
103       FORMAT(1X,'P_SUM:  ',15X,5(E10.3,1X))
          WRITE(IFCH,104)(PA0(K),K=1,5)
104       FORMAT(1X,'P_CLU:  ',15X,5(E10.3,1X))
        ENDIF

        IF ( IFAIL.EQ.0 .AND.
     *   (ABS(PX-PA0(1)).GT.1.E-2*ABS(PX).AND.ABS(PX-PA0(1)).GT.1.E-2
     *.OR.ABS(PY-PA0(2)).GT.1.E-2*ABS(PY).AND.ABS(PY-PA0(2)).GT.1.E-2
     *.OR.ABS(PZ-PA0(3)).GT.1.E-2*ABS(PZ).AND.ABS(PZ-PA0(3)).GT.1.E-2
     *.OR.ABS(E -PA0(4)).GT.1.E-2*ABS(E ).AND.ABS(E -PA0(4)).GT.1.E-2)
     *.OR.
     *       IFAIL.NE.0 .AND.
     *   (ABS(PX-PA0(1)).GT.1.E-2*ABS(PX).AND.ABS(PX-PA0(1)).GT.1.E-2
     *.OR.ABS(PY-PA0(2)).GT.1.E-2*ABS(PY).AND.ABS(PY-PA0(2)).GT.1.E-2
     *.OR.ABS(PZ-PA0(3)).GT.1.E-2*ABS(PZ).AND.ABS(PZ-PA0(3)).GT.1.E-2
     *.OR.ABS(E -PA0(4)).GT.35.E-1*ABS(E ).AND.ABS(E -PA0(4)).GT.35.E-1)
     *   ) THEN
          CALL UTMSG('JCLUDE')
          WRITE(IFCH,*)'*****  P_SUM /= P_CLU'
          WRITE(IFCH,*)'IFAIL_JRESCL:',IFAIL,'   SCAL:',SCAL
          DO 30 N=NPTLB+1,NPTL
            WRITE(IFCH,101)N,IDPTL(N),(PPTL(K,N),K=1,5)
30        CONTINUE
          WRITE(IFCH,103)PX,PY,PZ,E
          WRITE(IFCH,104)(PA0(K),K=1,4)
          CALL UTMSGF
        ENDIF
      ENDIF

1000  CONTINUE
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)('-',L=1,25)
        WRITE(IFCH,*)'   RETURN FROM JCLUDE   '
        WRITE(IFCH,*)('-',L=1,25)
        WRITE(IFCH,*)' '
      ENDIF
      ISH=ISH0
      RETURN

1001  IRET=1
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('JCLUDE')
        WRITE(IFCH,*)'*****  CLUSTER DECAY NOT POSSIBLE  --> IRET=1'
        WRITE(IFCH,100)(PA0(K),K=1,5)
        WRITE(IFCH,*)'JCA:'
        WRITE(IFCH,*)JCA0
        CALL UTMSGF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE JDECA(I,IRET)

C-----------------------------------------------------------------------
C  DECAYS I (CALLS JDECAYV)
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (NFLAV=6)
      COMMON /CCLUDE/  KCLUDE
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION TOR,ZOR
      INTEGER          JCDU(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IRET=0
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 5 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)'ENTRY JDECA. DECAY OF',I,IDPTL(I),PPTL(5,I)
      ENDIF
CDH   CALL IDMASS(111,AMRHO0)
CDH   CALL IDMASS(221,AMOMEG)
      IOI=IORPTL(I)
      IF ( .NOT.( IOI.GT.0 .AND. ABS(IDPTL(IOI)).LT.10000
     *                    .AND. JORPTL(I).EQ.0 ) ) THEN
        IF ( IDPTL(I) .EQ. 111 ) IDPTL(I)=221
        IF ( IDPTL(I) .EQ. 221 .AND. RANGEN() .GT. 0.5 ) IDPTL(I)=111
      ENDIF
      IF ( MOD(NDECAY        ,10) .EQ.1 ) GOTO 1000
      IDA=ABS(IDPTL(I))
      IF ( MOD(NDECAY/10     ,10) .EQ.1 .AND. IDA .EQ.  20 ) GOTO 1000
      IF ( MOD(NDECAY/100    ,10) .EQ.1 .AND. IDA .EQ.2130 ) GOTO 1000
      IF ( MOD(NDECAY/1000   ,10) .EQ.1 ) THEN
        IF ( IDA .EQ.1130 ) GOTO 1000
        IF ( IDA .EQ.2230 ) GOTO 1000
      ENDIF
      IF ( MOD(NDECAY/10000  ,10) .EQ.1 ) THEN
        IF ( IDA .EQ.2330 ) GOTO 1000
        IF ( IDA .EQ.1330 ) GOTO 1000
      ENDIF
      IF ( MOD(NDECAY/100000 ,10) .EQ.1 .AND. IDA .EQ.3331 ) GOTO 1000
      IF ( MOD(NDECAY/1000000,10) .EQ.1 .AND. IDA .EQ. 110 ) GOTO 1000
      IF ( MOD(NDECAX        ,10) .EQ.1 .AND. IDA .EQ. 441 ) GOTO 1000
      IF ( MOD(NDECAX/10     ,10) .EQ.1 .AND. IDA .EQ. 230 ) GOTO 1000
      IF ( MOD(NDECAX/100    ,10) .EQ.1 ) THEN
        IF ( IDA .EQ.1111 ) GOTO 1000
        IF ( IDA .EQ.1121 ) GOTO 1000
        IF ( IDA .EQ.1221 ) GOTO 1000
        IF ( IDA .EQ.2221 ) GOTO 1000
      ENDIF

      IF ( MOD(NDECAX/1000   ,10) .EQ.1 ) THEN
        IF ( IDA .EQ. 111 ) GOTO 1000
        IF ( IDA .EQ. 121 ) GOTO 1000
        IF ( IDA .EQ. 221 ) GOTO 1000
        IF ( IDA .EQ. 331 ) GOTO 1000
      ENDIF
      IF ( MOD(NDECAX/10000  ,10) .EQ.1 .AND. IDA .EQ. 220 ) GOTO 1000
      IF ( MOD(NDECAX/100000 ,10) .EQ.1 .AND. IDA .EQ. 330 ) GOTO 1000
      IF ( MOD(NDECAX/1000000,10) .EQ.1 ) THEN
        IF ( IDA .EQ. 112 ) GOTO 1000
        IF ( IDA .EQ. 122 ) GOTO 1000
      ENDIF
      IF ( MOD(NDECAW        ,10) .EQ.1 .AND. IDA .EQ. 332 ) GOTO 1000
      IF ( MOD(NDECAW/10     ,10) .EQ.1 ) THEN
        IF ( IDA .EQ. 131 ) GOTO 1000
        IF ( IDA .EQ.-131 ) GOTO 1000
        IF ( IDA .EQ. 231 ) GOTO 1000
        IF ( IDA .EQ.-231 ) GOTO 1000
      ENDIF
      T=TIVPTL(2,I)
      NPTLB=NPTL
      IF ( NPTL .GT. MXPTL-10 ) THEN
        CALL UTSTOP('JDECA: MXPTL TOO SMALL                  ')
      ENDIF
      ISH=ISH0
      CALL JDECAYV(I,IRET)
      IF ( ISHSUB/100 .EQ. 5 ) ISH=MOD(ISHSUB,100)
      IF ( IRET .EQ. 1 ) GOTO 1000
      IF ( NPTL .LE. NPTLB ) GOTO 1000
      ISH00=ISH
      IF ( ISHSUB/100.EQ.14 .AND. KCLUDE.EQ.1 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,101)SNGL(TTAUS)
101     FORMAT(1X,'DECAY  AT TAU =',E10.3)
        WRITE(IFCH,115)I,IDPTL(I)
     *                 ,(PPTL(K,I),K=3,5),(TIVPTL(K,I),K=1,2)
115     FORMAT(1X,'/CPTL/',I6,I10
     *           ,1X,3(E10.2),1X,2(E10.2))
      ENDIF
      ISTPTL(I)=1
      IFRPTL(1,I)=NPTLB+1
      IFRPTL(2,I)=NPTL
      X=XORPTL(1,I)+(T-XORPTL(4,I))*PPTL(1,I)/PPTL(4,I)
      Y=XORPTL(2,I)+(T-XORPTL(4,I))*PPTL(2,I)/PPTL(4,I)
      Z=XORPTL(3,I)+(T-XORPTL(4,I))*PPTL(3,I)/PPTL(4,I)
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)
     *            'LOOP OVER DECAY PRODUCTS ',NPTLB+1,' - ',NPTL,' :'
      DO 20 N=NPTLB+1,NPTL
        IF ( ISH .GE. 93 ) WRITE(IFCH,*)'PARTICLE: ',N,IDPTL(N)
        IORPTL(N)=I
        JORPTL(N)=0
        ISTPTL(N)=0
        IFRPTL(1,N)=0
        IFRPTL(2,N)=0
        XORPTL(1,N)=X
        XORPTL(2,N)=Y
        XORPTL(3,N)=Z
        XORPTL(4,N)=T
        NQJPTL(N)=NQJPTL(I)
        IO=N
 1      IO=IORPTL(IO)
        IF ( ISH .GE. 93 ) WRITE(IFCH,*)'IO = ',IO
        IF ( IORPTL(IO) .GT. 0 ) GOTO 1
        IF ( ISH.GE. 93 ) WRITE(IFCH,*)'ORIGIN: ',IO,IDPTL(IO)
        ZOR=XORPTL(3,IO)
        TOR=XORPTL(4,IO)
        CALL IDQUAC(IO,NQ,NDU,NDU,JCDU)
        R=RANGEN()
        TAURAN=-TAUREA*LOG(R)
        CALL UTTAIX(N,TAURAN,ZOR,TOR,ZIS,TIS)
        TIVPTL(1,N)=MAX(T,TIS)
        CALL IDTAU(IDPTL(N),PPTL(4,N),PPTL(5,N),TAUGM)
        TIVPTL(2,N)=T+TAUGM
        ICLPTL(N)=1
        IF ( ISH .GE. 91 ) WRITE(IFCH,115)N,IDPTL(N)
     *                    ,(PPTL(K,N),K=3,5),(TIVPTL(K,N),K=1,2)
20    CONTINUE
      ISH=ISH00

1000  CONTINUE
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)('-',L=1,79)
      ISH=ISH0
      RETURN
      END
C=======================================================================

      SUBROUTINE JDECAYV(IP,IRET)

C-----------------------------------------------------------------------
C  DECAYS PARTICLE IP FROM /CPTL/
C-----------------------------------------------------------------------
      PARAMETER (MXDKY=2000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXPTL=70000)
      COMMON /CCLUDE/  KCLUDE
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /WCO/     WGAM2,WMASS2

      REAL    BETA(3),PGEN(5,5),PREST(4,5)
     *       ,REDUCE(5),RND(5),U(3)
      DATA REDUCE/1.,1.,2.,5.,15./,TWOME/1.022006E-3/
      SAVE
C-----------------------------------------------------------------------
C  FCTN DEFINITIONS
      DOT(I1,I2)=PREST(4,I1)*PREST(4,I2)-PREST(1,I1)*PREST(1,I2)
     *          -PREST(2,I1)*PREST(2,I2)-PREST(3,I1)*PREST(3,I2)
C  CHARGED W PROPAGATOR.
      WPROP(Z)=(Z-WMASS2**2)**2+(WMASS2*WGAM2)**2
C-----------------------------------------------------------------------

      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 4 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)'DECAY OF',IP,IDPTL(IP),PPTL(5,IP)
      ENDIF
      ISH=ISH0

      IRET=0
      KCLUDE=0

C  NO K_LONG DECAY
C  ---------------
      IF ( IDPTL(IP) .EQ. -20 ) GOTO 1000

C  CLUSTER DECAY
C  -------------
      IF ( ABS(IDPTL(IP)) .GT. 100000000 ) THEN
        KCLUDE=1
        CALL JCLUDE(IP,IRET)
        GOTO 1000
      ENDIF

      IF ( ISHSUB/100 .EQ. 4 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)'ORDINARY DECAY'

C  SELECT DECAY MODE
C  -----------------
      NTRY=0
 2    NTRY=NTRY+1
      IF ( NTRY .GT. 100 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JDECAYV')
          WRITE(IFCH,*)'*****  DECAY NOT POSSIBLE. IRET = 1.'
          WRITE(IFCH,*)'ID,MASS: ',IDPTL(IP),PPTL(5,IP)
          CALL UTMSGF
        ENDIF
        IRET=1
        GOTO 1000
      ENDIF
      IDLV1=IDPTL(IP)
      AMSS=PPTL(5,IP)
ctp060203  1    CONTINUE
      IPOINT=LOOK(IABS(IDLV1))-1
      IF ( IPOINT .LT. 0 ) GOTO 1000
      TRY=RANGEN()
100   IPOINT=IPOINT+1
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)'IPOINT,CBR,TRY'
     *                                ,IPOINT,CBR(IPOINT),TRY
      IF ( TRY .GT. CBR(IPOINT) ) GOTO 100
      NADD=0
      SUM=0.
      NSTART=NPTL+1
      DO 110 I=1,5
        IF ( MODE(I,IPOINT) .EQ. 0 ) GOTO 110
        IF ( NPTL+NADD+1 .GT. MXPTL ) GOTO 9999
        NADD=NADD+1
        NEW=NPTL+NADD
        IDPTL(NEW)=MODE(I,IPOINT)
        IDLV1=IDPTL(NEW)
        CALL IDMASS(IDLV1,PPTL(5,NEW))
        SUM=SUM+PPTL(5,NEW)
110   CONTINUE
      IF ( NADD .NE. 1  .AND.  SUM+1.E-2 .GE. AMSS ) GOTO 2
      NADD1=NADD-1
      DO 120 J=1,5
        PGEN(J,1)=PPTL(J,IP)
120   CONTINUE
      PGEN(5,NADD)=PPTL(5,NPTL+NADD)
      IF ( NADD .EQ. 1 ) GOTO 700
      IF ( NADD .EQ. 2 ) GOTO 400

C  USE KROLL-WADA DISTRIBUTION FOR DALITZ DECAYS.
C  ----------------------------------------------
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)'>= 3 BODY DECAY'
      IF ( .NOT. ( (IDPTL(IP).EQ.110 .OR. IDPTL(IP).EQ.220) .AND.
     *             ABS(IDPTL(NPTL+2)).EQ.12 ) ) GOTO 130
      NTRY=0
125   NTRY=NTRY+1
      IF ( NTRY .GT. 10 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JDECAYV')
          WRITE(IFCH,*)'*****  NTRY > 10. IRET = 1.'
          WRITE(IFCH,*)'AMEE,REE,WTEE',AMEE,REE,WTEE
          CALL UTMSGF
        ENDIF
        IRET=1
        GOTO 1000
      ENDIF
      AMEE=TWOME*(PPTL(5,IP)/TWOME)**RANGEN()
      REE=(TWOME/AMEE)**2
      WTEE=(1.-(AMEE/PPTL(5,IP))**2)**3*SQRT(1.-REE)*(1.+.5*REE)
      IF ( WTEE .LT. RANGEN() ) GOTO 125
      PGEN(5,2)=AMEE
      GOTO 400
130   CONTINUE

C  CALCULATE MAXIMUM PHASE-SPACE WEIGHT
C  ------------------------------------
      WTMAX=1./REDUCE(NADD)
      SUM1=PGEN(5,1)
      SUM2=SUM-PPTL(5,NPTL+1)
      DO 200 I=1,NADD1
        WTMAX=WTMAX*UTPCM(SUM1,SUM2,PPTL(5,NPTL+I))
        SUM1=SUM1-PPTL(5,NPTL+I)
        SUM2=SUM2-PPTL(5,NPTL+I+1)
200   CONTINUE

C  GENERATE UNIFORM NADD-BODY PHASE SPACE
C  --------------------------------------
      NTRY=0
300   NTRY=NTRY+1
      IF ( NTRY .GT. 10000 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JDECAYV')
          WRITE(IFCH,*)'*****  INFINITE LOOP (2). IRET = 1.'
          WRITE(IFCH,*)'IP,IDPTL(IP),PPTL(5,IP):'
     *                 ,IP,IDPTL(IP),PPTL(5,IP)
          WRITE(IFCH,*)'WT,WTMAX:',WT,WTMAX
          WRITE(IFCH,*)'I,PGEN(5,I),PPTL(5,NPTL+I),IDPTL(NPTL+I):'
          DO 305 I=1,NADD
            WRITE(IFCH,*)I,PGEN(5,I),PPTL(5,NPTL+I),IDPTL(NPTL+I)
305       CONTINUE
          CALL UTMSGF
        ENDIF
        IRET=1
        GOTO 1000
      ENDIF
      RND(1)=1.
      DO 310 I=2,NADD1
        RNEW=RANGEN()
        I1=I-1
        DO 320 JJ1=1,I1
          J=I-JJ1
          JSAVE=J+1
          IF ( RNEW .LE. RND(J) ) GOTO 315
          RND(JSAVE)=RND(J)
320     CONTINUE
315     RND(JSAVE)=RNEW
310   CONTINUE
      RND(NADD)=0.
      WT=1.
      SUM1=SUM
      DO 330 I=2,NADD
        SUM1=SUM1-PPTL(5,NPTL+I-1)
        PGEN(5,I)=SUM1+RND(I)*(PGEN(5,1)-SUM)
        A=PGEN(5,I-1)
        B=PGEN(5,I)
        C=PPTL(5,NPTL+I-1)
        WT=WT*UTPCM(A,B,C)
330   CONTINUE
      IF ( WT .LT. RANGEN()*WTMAX ) GOTO 300

C  CARRY OUT TWO-BODY DECAYS IN PGEN FRAMES
C  ----------------------------------------
400   CONTINUE
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)'2 BODY DECAY'
      DO 410 I=1,NADD1
        QCM=UTPCM(PGEN(5,I),PGEN(5,I+1),PPTL(5,NPTL+I))
        U(3)=2.*RANGEN()-1.
        PHI=2.*PI*RANGEN()
        U(1)=SQRT(1.-U(3)**2)*COS(PHI)
        U(2)=SQRT(1.-U(3)**2)*SIN(PHI)
        PPTL(1,NPTL+I)=QCM*U(1)
        PGEN(1,I+1)=-PPTL(1,NPTL+I)
        PPTL(2,NPTL+I)=QCM*U(2)
        PGEN(2,I+1)=-PPTL(2,NPTL+I)
        PPTL(3,NPTL+I)=QCM*U(3)
        PGEN(3,I+1)=-PPTL(3,NPTL+I)
        PPTL(4,NPTL+I)=SQRT(QCM**2+PPTL(5,NPTL+I)**2)
        PGEN(4,I+1)=SQRT(QCM**2+PGEN(5,I+1)**2)
410   CONTINUE
      PPTL(1,NPTL+NADD)=PGEN(1,NADD)
      PPTL(2,NPTL+NADD)=PGEN(2,NADD)
      PPTL(3,NPTL+NADD)=PGEN(3,NADD)
      PPTL(4,NPTL+NADD)=PGEN(4,NADD)

C  BOOST PGEN FRAMES TO LAB FRAME
C       ALSO SAVE MOMENTA IN REST FRAME (LAST FRAME)
C  -------------------------------------------------
      DO 500 II=1,NADD1
        I=NADD-II
        BETA(1)=1./PGEN(4,I)*PGEN(1,I)
        BETA(2)=1./PGEN(4,I)*PGEN(2,I)
        BETA(3)=1./PGEN(4,I)*PGEN(3,I)
        GAMMA=PGEN(4,I)/PGEN(5,I)
        DO 520 K=I,NADD
          K1=NPTL+K
          BP=BETA(1)*PPTL(1,K1)+BETA(2)*PPTL(2,K1)+BETA(3)*PPTL(3,K1)
          AUXIL=GAMMA*(PPTL(4,K1)+BP*GAMMA/(GAMMA+1.))
          PREST(1,K)=PPTL(1,K1)
          PPTL(1,K1)=PPTL(1,K1)+BETA(1)*AUXIL
          PREST(2,K)=PPTL(2,K1)
          PPTL(2,K1)=PPTL(2,K1)+BETA(2)*AUXIL
          PREST(3,K)=PPTL(3,K1)
          PPTL(3,K1)=PPTL(3,K1)+BETA(3)*AUXIL
          PREST(4,K)=PPTL(4,K1)
          PPTL(4,K1)=GAMMA*(PPTL(4,K1)+BP)
520     CONTINUE
500   CONTINUE

C  MATRIX ELEMENTS
C  ---------------
      IF ( NADD .EQ. 3 ) THEN
        IF ( IDPTL(IP) .EQ. 221  .OR.  IDPTL(IP) .EQ. 331 ) GOTO 610
        IF ( ABS(IDPTL(NPTL+1)) .LT. 20   .AND.
     *             IDPTL(NPTL+1) .NE. 10 ) GOTO 620
      ENDIF
      GOTO 800

C  OMEG AND PHI DECAY
C       USE VECTORS IN REST FRAME
C  ------------------------------
610   WT=(PPTL(5,NPTL+1)*PPTL(5,NPTL+2)*PPTL(5,NPTL+3))**2
     *                 -(PPTL(5,NPTL+1)*DOT(2,3))**2
     *                 -(PPTL(5,NPTL+2)*DOT(1,3))**2
     *                 -(PPTL(5,NPTL+3)*DOT(1,2))**2
     *                +2.*DOT(1,2)*DOT(2,3)*DOT(1,3)
      IF ( WT .LT. RANGEN()*PPTL(5,IP)**6/108. ) GOTO 300
      GOTO 800

C  SEMILEPTONIC AND QUARK DECAYS
C       USE VECTORS IN REST FRAME, WHERE IP HAS (M,0,0,0)
C       INCLUDE W PROPAGATOR
C  ------------------------------------------------------
620   WT=(PPTL(5,IP)*PREST(4,2))*DOT(1,3)
      S12=PPTL(5,NPTL+1)**2+PPTL(5,NPTL+2)**2+2.*DOT(1,2)
      S12MAX=PPTL(5,IP)**2
      WT=WT*WPROP(S12MAX)/WPROP(S12)
      IF ( WT .LT. RANGEN()*PPTL(5,IP)**4/16. ) GOTO 300
      GOTO 800

C  ONE-PARTICLE DECAYS
C  -------------------
700   CONTINUE
      DO 710 J=1,5
        PPTL(J,NPTL+1)=PPTL(J,IP)
710   CONTINUE

C  SWAP PARTICLES AND ANTIPARTICLES IF IDPTL(IP)<0
C  -----------------------------------------------
800   CONTINUE
      IF ( IDPTL(IP).GE.0 .OR. ABS(IDPTL(IP)).EQ.20 ) GOTO 900
      DO 810 I=1,NADD
        IDABS=ABS(IDPTL(NPTL+I))
        IFL1=IDABS/1000
        IFL2=MOD(IDABS/100,10)
        IFL3=MOD(IDABS/10,10)
        IF ( IFL1.EQ.0 .AND. IFL2.NE.0 .AND. IFL2.EQ.IFL3 ) GOTO 810
        IF ( IDABS.EQ.9  .OR. IDABS.EQ.10 .OR. IDABS.EQ.20 ) GOTO 810
        IF ( IDABS.EQ.29 .OR. IDABS.EQ.30 .OR. IDABS.EQ.40 ) GOTO 810
        IDPTL(NPTL+I)=-IDPTL(NPTL+I)
810   CONTINUE

900   CONTINUE
      NPTL=NPTL+NADD
      IF ( NPTL .GT. MXPTL ) THEN
        CALL UTSTOP('JDECAYV: NPTL>MXPTL                      ')
      ENDIF
      NQK=0
      IF ( ABS(IDPTL(NPTL)).LT.10 .OR. MOD(IDPTL(NPTL),100).EQ.0 ) THEN
        CALL UTSTOP('JDECAYV: DECAY PTCL IS PARTON            ')
      ENDIF

1000  CONTINUE
      IF ( ISH .GE. 93 ) WRITE(IFCH,*)('-',L=1,79)
      ISH=ISH0
      RETURN

9999  CALL UTSTOP('JDECAYV: MXPTL TOO SMALL                 ')
      RETURN
      END
C=======================================================================

      SUBROUTINE JDECIN(LPRINT)

C-----------------------------------------------------------------------
C  SETS UP /DKYTAB/
C-----------------------------------------------------------------------
      PARAMETER (MXDKY=2000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (NDECTB=1171)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      LOGICAL          NODCAY,NOETA,NOEVOL,NOHADR,NONUNU,NOPI0
      COMMON /NODCAY/  NODCAY,NOETA,NOEVOL,NOHADR,NONUNU,NOPI0
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /WCO/     WGAM2,WMASS2

      REAL        DECTAB(7,NDECTB)
      INTEGER     IMODE(6)
      CHARACTER*8 IBLANK,IDLABL,LMODE(6),LRES
      LOGICAL     LPRINT

      DATA IBLANK/' '/
      DATA ((DECTAB(I,J),I=1,7),J=  1, 18)/
     *  110., .98850,  10.,  10.,   0.,   0.,   0.
     *, 110.,1.00000,  10.,  12., -12.,   0.,   0.
     *, 220., .38000,  10.,  10.,   0.,   0.,   0.
     *, 220., .71000, 110., 110., 110.,   0.,   0.
     *, 220., .94600, 120.,-120., 110.,   0.,   0.
     *, 220., .99500, 120.,-120.,  10.,   0.,   0.
     *, 220.,1.00000,  10.,  12., -12.,   0.,   0.
     *, 330., .44100, 220., 120.,-120.,   0.,   0.
     *, 330., .66100, 220., 110., 110.,   0.,   0.
     *, 330., .95900, 111.,  10.,   0.,   0.,   0.
     *, 330., .98000, 221.,  10.,   0.,   0.,   0.
     *, 330.,1.00000,  10.,  10.,   0.,   0.,   0.
     *, 121.,1.00000, 120., 110.,   0.,   0.,   0.
     *, 111., .99989, 120.,-120.,   0.,   0.,   0.
     *, 111., .99993,  12., -12.,   0.,   0.,   0.
     *, 111.,1.00000,  14., -14.,   0.,   0.,   0.
     *, 221., .89900, 120.,-120., 110.,   0.,   0.
     *, 221., .91200, 120.,-120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J= 19, 36)/
     *  221., .99992, 110.,  10.,   0.,   0.,   0.
     *, 221.,1.00000,  12., -12.,   0.,   0.,   0.
     *, 331., .48600, 130.,-130.,   0.,   0.,   0.
     *, 331., .83700, 230.,-230.,   0.,   0.,   0.
     *, 331., .98400, 120.,-120., 110.,   0.,   0.
     *, 331., .99944, 220.,  10.,   0.,   0.,   0.
     *, 331., .99975,  12., -12.,   0.,   0.,   0.
     *, 331.,1.00000,  14., -14.,   0.,   0.,   0.
     *, 230., .50000,  20.,   0.,   0.,   0.,   0.
     *, 230.,1.00000, -20.,   0.,   0.,   0.,   0.
     *, 131., .66670, 230., 120.,   0.,   0.,   0.
     *, 131.,1.00000, 130., 110.,   0.,   0.,   0.
     *, 231., .66670, 130.,-120.,   0.,   0.,   0.
     *, 231.,1.00000, 230., 110.,   0.,   0.,   0.
     *, 240., .11000,  12., -11., 230.,   0.,   0.
     *, 240., .17000,  12., -11., 231.,   0.,   0.
     *, 240., .28000,  14., -13., 230.,   0.,   0.
     *, 240., .34000,  14., -13., 231.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J= 37, 54)/
     *  240., .37800, 230.,-120.,   0.,   0.,   0.
     *, 240., .56300, 230.,-121.,   0.,   0.,   0.
     *, 240., .60800, 231.,-120.,   0.,   0.,   0.
     *, 240., .62100, 230.,-120., 110.,   0.,   0.
     *, 240., .71000, 130.,-120.,-120.,   0.,   0.
     *, 240., .80100, 230.,-120.,-120., 120.,   0.
     *, 240., .87900, 130.,-120.,-120., 110.,   0.
     *, 240., .95400, 230.,-120., 110., 110.,   0.
     *, 240., .96600, 230.,-130.,   0.,   0.,   0.
     *, 240., .97600, 331.,-120.,   0.,   0.,   0.
     *, 240., .98800,-130., 231.,   0.,   0.,   0.
     *, 240.,1.00000,-131., 230.,   0.,   0.,   0.
     *, 140., .04500,  12., -11., 130.,   0.,   0.
     *, 140., .07500,  12., -11., 131.,   0.,   0.
     *, 140., .12000,  14., -13., 130.,   0.,   0.
     *, 140., .15000,  14., -13., 131.,   0.,   0.
     *, 140., .20300, 130.,-120.,   0.,   0.,   0.
     *, 140., .22700, 230., 110.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J= 55, 72)/
     *  140., .24700, 230., 220.,   0.,   0.,   0.
     *, 140., .28900, 230., 221.,   0.,   0.,   0.
     *, 140., .45100, 130.,-121.,   0.,   0.,   0.
     *, 140., .53600, 131.,-120.,   0.,   0.,   0.
     *, 140., .56200, 231., 110.,   0.,   0.,   0.
     *, 140., .57600, 230., 111.,   0.,   0.,   0.
     *, 140., .58700, 130.,-120., 110.,   0.,   0.
     *, 140., .60300, 230.,-120., 120.,   0.,   0.
     *, 140., .72700, 130.,-120.,-120., 120.,   0.
     *, 140., .87600, 230.,-120., 120., 110.,   0.
     *, 140., .96900, 130.,-120., 110., 110.,   0.
     *, 140.,1.00000, 230., 110., 110., 110.,   0.
     *, 340., .03250,  12., -11., 220.,   0.,   0.
     *, 340., .06500,  12., -11., 331.,   0.,   0.
     *, 340., .09750,  14., -13., 220.,   0.,   0.
     *, 340., .13000,  14., -13., 331.,   0.,   0.
     *, 340., .17900,-130., 230.,   0.,   0.,   0.
     *, 340., .22800,-120., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J= 73, 90)/
     *  340., .33800,-131., 230.,   0.,   0.,   0.
     *, 340., .44800,-130., 231.,   0.,   0.,   0.
     *, 340., .55800,-120., 331.,   0.,   0.,   0.
     *, 340., .57500,-130., 230., 110.,   0.,   0.
     *, 340., .59200,-230., 230.,-120.,   0.,   0.
     *, 340., .69400,-130., 230.,-120., 120.,   0.
     *, 340., .79600,-130., 230., 110., 110.,   0.
     *, 340., .89800,-130., 130.,-120., 110.,   0.
     *, 340.,1.00000,-230., 230.,-120., 110.,   0.
     *, 241., .64000, 140.,-120.,   0.,   0.,   0.
     *, 241., .92000, 240., 110.,   0.,   0.,   0.
     *, 241.,1.00000, 240.,  10.,   0.,   0.,   0.
     *, 141., .55000, 140., 110.,   0.,   0.,   0.
     *, 141.,1.00000, 140.,  10.,   0.,   0.,   0.
     *, 341.,1.00000, 340.,  10.,   0.,   0.,   0.
     *, 441., .07400,  12., -12.,   0.,   0.,   0.
     *, 441., .14800,  14., -14.,   0.,   0.,   0.
     *, 441., .15210,-121., 120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J= 91,108)/
     *  441., .15620, 111., 110.,   0.,   0.,   0.
     *, 441., .16020, 121.,-120.,   0.,   0.,   0.
     *, 441., .16300,-121., 111., 120.,   0.,   0.
     *, 441., .16580, 121.,-121., 110.,   0.,   0.
     *, 441., .16860, 121., 111.,-120.,   0.,   0.
     *, 441., .28740, 120.,-120., 130.,-130.,   0.
     *, 441., .40620, 110., 110., 130.,-130.,   0.
     *, 441., .52500, 120.,-120., 120.,-120.,   0.
     *, 441., .64380, 120.,-120., 110., 110.,   0.
     *, 441., .76260, 110., 110., 110., 110.,   0.
     *, 441., .88130, 120.,-120., 230.,-230.,   0.
     *, 441.,1.00000, 110., 110., 230., 230.,   0.
     *, 150., .06000, -12.,  11., 140.,   0.,   0.
     *, 150., .12000, -12.,  11., 141.,   0.,   0.
     *, 150., .18000, -14.,  13., 140.,   0.,   0.
     *, 150., .24000, -14.,  13., 141.,   0.,   0.
     *, 150., .25500, -16.,  15., 140.,   0.,   0.
     *, 150., .27000, -16.,  15., 141.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=109,126)/
     *  150., .28050, 140., 120.,   0.,   0.,   0.
     *, 150., .29100, 140., 121.,   0.,   0.,   0.
     *, 150., .30150, 141., 120.,   0.,   0.,   0.
     *, 150., .31200, 141., 121.,   0.,   0.,   0.
     *, 150., .32650, 140.,-340.,   0.,   0.,   0.
     *, 150., .34100, 140.,-341.,   0.,   0.,   0.
     *, 150., .35550, 141.,-340.,   0.,   0.,   0.
     *, 150., .37000, 141.,-341.,   0.,   0.,   0.
     *, 150., .39800, 140., 120., 110.,   0.,   0.
     *, 150., .42600, 140., 120., 220.,   0.,   0.
     *, 150., .45400, 140., 120., 111.,   0.,   0.
     *, 150., .48200, 140., 120., 221.,   0.,   0.
     *, 150., .51000, 140., 121., 110.,   0.,   0.
     *, 150., .53800, 140., 121., 220.,   0.,   0.
     *, 150., .56600, 140., 121., 111.,   0.,   0.
     *, 150., .59400, 140., 121., 221.,   0.,   0.
     *, 150., .62200, 141., 120., 110.,   0.,   0.
     *, 150., .65000, 141., 120., 220.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=127,144)/
     *  150., .67800, 141., 120., 111.,   0.,   0.
     *, 150., .70600, 141., 120., 221.,   0.,   0.
     *, 150., .73400, 141., 121., 110.,   0.,   0.
     *, 150., .76200, 141., 121., 220.,   0.,   0.
     *, 150., .79000, 141., 121., 111.,   0.,   0.
     *, 150., .81800, 141., 121., 221.,   0.,   0.
     *, 150., .83200, 140., 130.,-230.,   0.,   0.
     *, 150., .84600, 140., 130.,-231.,   0.,   0.
     *, 150., .86000, 140., 131.,-230.,   0.,   0.
     *, 150., .87400, 140., 131.,-231.,   0.,   0.
     *, 150., .88800, 141., 130.,-230.,   0.,   0.
     *, 150., .90200, 141., 130.,-231.,   0.,   0.
     *, 150., .91600, 141., 131.,-230.,   0.,   0.
     *, 150., .93000, 141., 131.,-231.,   0.,   0.
     *, 150., .93300, 140.,-140., 130.,   0.,   0.
     *, 150., .93600, 140.,-140., 131.,   0.,   0.
     *, 150., .93900, 140.,-141., 130.,   0.,   0.
     *, 150., .94200, 140.,-141., 131.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=145,162)/
     *  150., .94500, 141.,-140., 130.,   0.,   0.
     *, 150., .94800, 141.,-140., 131.,   0.,   0.
     *, 150., .95100, 141.,-141., 130.,   0.,   0.
     *, 150., .95400, 141.,-141., 131.,   0.,   0.
     *, 150., .95700, 140.,-240., 230.,   0.,   0.
     *, 150., .96000, 140.,-240., 231.,   0.,   0.
     *, 150., .96300, 140.,-241., 230.,   0.,   0.
     *, 150., .96600, 140.,-241., 231.,   0.,   0.
     *, 150., .96900, 141.,-240., 230.,   0.,   0.
     *, 150., .97200, 141.,-240., 231.,   0.,   0.
     *, 150., .97500, 141.,-241., 230.,   0.,   0.
     *, 150., .97800, 141.,-241., 231.,   0.,   0.
     *, 150., .97950, 140.,-340., 330.,   0.,   0.
     *, 150., .98100, 140.,-340., 331.,   0.,   0.
     *, 150., .98250, 140.,-341., 331.,   0.,   0.
     *, 150., .98400, 140.,-341., 331.,   0.,   0.
     *, 150., .98550, 141.,-340., 330.,   0.,   0.
     *, 150., .98700, 141.,-340., 331.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=163,180)/
     *  150., .98850, 141.,-341., 331.,   0.,   0.
     *, 150., .99000, 141.,-341., 331.,   0.,   0.
     *, 150., .99200, 441., 130., 110.,   0.,   0.
     *, 150., .99400, 441., 131., 110.,   0.,   0.
     *, 150., .99600, 441., 230., 120.,   0.,   0.
     *, 150., .99800, 441., 231., 120.,   0.,   0.
     *, 150., .99900, 441., 330., 130.,   0.,   0.
     *, 150.,1.00000, 441., 331., 130.,   0.,   0.
     *, 250., .06000, -12.,  11., 240.,   0.,   0.
     *, 250., .12000, -12.,  11., 241.,   0.,   0.
     *, 250., .18000, -14.,  13., 240.,   0.,   0.
     *, 250., .24000, -14.,  13., 241.,   0.,   0.
     *, 250., .25500, -16.,  15., 240.,   0.,   0.
     *, 250., .27000, -16.,  15., 241.,   0.,   0.
     *, 250., .28050, 240., 120.,   0.,   0.,   0.
     *, 250., .29100, 240., 121.,   0.,   0.,   0.
     *, 250., .30150, 241., 120.,   0.,   0.,   0.
     *, 250., .31200, 241., 121.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=181,198)/
     *  250., .32650, 240.,-340.,   0.,   0.,   0.
     *, 250., .34100, 240.,-341.,   0.,   0.,   0.
     *, 250., .35550, 241.,-340.,   0.,   0.,   0.
     *, 250., .37000, 241.,-341.,   0.,   0.,   0.
     *, 250., .39800, 240., 120., 110.,   0.,   0.
     *, 250., .42600, 240., 120., 220.,   0.,   0.
     *, 250., .45400, 240., 120., 111.,   0.,   0.
     *, 250., .48200, 240., 120., 221.,   0.,   0.
     *, 250., .51000, 240., 121., 110.,   0.,   0.
     *, 250., .53800, 240., 121., 220.,   0.,   0.
     *, 250., .56600, 240., 121., 111.,   0.,   0.
     *, 250., .59400, 240., 121., 221.,   0.,   0.
     *, 250., .62200, 241., 120., 110.,   0.,   0.
     *, 250., .65000, 241., 120., 220.,   0.,   0.
     *, 250., .67800, 241., 120., 111.,   0.,   0.
     *, 250., .70600, 241., 120., 221.,   0.,   0.
     *, 250., .73400, 241., 121., 110.,   0.,   0.
     *, 250., .76200, 241., 121., 220.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=199,216)/
     *  250., .79000, 241., 121., 111.,   0.,   0.
     *, 250., .81800, 241., 121., 221.,   0.,   0.
     *, 250., .83200, 240., 130.,-230.,   0.,   0.
     *, 250., .84600, 240., 130.,-231.,   0.,   0.
     *, 250., .86000, 240., 131.,-230.,   0.,   0.
     *, 250., .87400, 240., 131.,-231.,   0.,   0.
     *, 250., .88800, 241., 130.,-230.,   0.,   0.
     *, 250., .90200, 241., 130.,-231.,   0.,   0.
     *, 250., .91600, 241., 131.,-230.,   0.,   0.
     *, 250., .93000, 241., 131.,-231.,   0.,   0.
     *, 250., .93300, 240.,-140., 130.,   0.,   0.
     *, 250., .93600, 240.,-140., 131.,   0.,   0.
     *, 250., .93900, 240.,-141., 130.,   0.,   0.
     *, 250., .94200, 240.,-141., 131.,   0.,   0.
     *, 250., .94500, 241.,-140., 130.,   0.,   0.
     *, 250., .94800, 241.,-140., 131.,   0.,   0.
     *, 250., .95100, 241.,-141., 130.,   0.,   0.
     *, 250., .95400, 241.,-141., 131.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=217,234)/
     *  250., .95700, 240.,-240., 230.,   0.,   0.
     *, 250., .96000, 240.,-240., 231.,   0.,   0.
     *, 250., .96300, 240.,-241., 230.,   0.,   0.
     *, 250., .96600, 240.,-241., 231.,   0.,   0.
     *, 250., .96900, 241.,-240., 230.,   0.,   0.
     *, 250., .97200, 241.,-240., 231.,   0.,   0.
     *, 250., .97500, 241.,-241., 230.,   0.,   0.
     *, 250., .97800, 241.,-241., 231.,   0.,   0.
     *, 250., .97950, 240.,-340., 330.,   0.,   0.
     *, 250., .98100, 240.,-340., 331.,   0.,   0.
     *, 250., .98250, 240.,-341., 331.,   0.,   0.
     *, 250., .98400, 240.,-341., 331.,   0.,   0.
     *, 250., .98550, 241.,-340., 330.,   0.,   0.
     *, 250., .98700, 241.,-340., 331.,   0.,   0.
     *, 250., .98850, 241.,-341., 331.,   0.,   0.
     *, 250., .99000, 241.,-341., 331.,   0.,   0.
     *, 250., .99200, 441., 130.,-120.,   0.,   0.
     *, 250., .99400, 441., 131.,-120.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=235,252)/
     *  250., .99600, 441., 230., 220.,   0.,   0.
     *, 250., .99800, 441., 231., 221.,   0.,   0.
     *, 250., .99900, 441., 330., 230.,   0.,   0.
     *, 250.,1.00000, 441., 331., 230.,   0.,   0.
     *, 350., .06000, -12.,  11., 340.,   0.,   0.
     *, 350., .12000, -12.,  11., 341.,   0.,   0.
     *, 350., .18000, -14.,  13., 340.,   0.,   0.
     *, 350., .24000, -14.,  13., 341.,   0.,   0.
     *, 350., .25500, -16.,  15., 340.,   0.,   0.
     *, 350., .27000, -16.,  15., 341.,   0.,   0.
     *, 350., .28050, 340., 120.,   0.,   0.,   0.
     *, 350., .29100, 340., 121.,   0.,   0.,   0.
     *, 350., .30150, 341., 120.,   0.,   0.,   0.
     *, 350., .31200, 341., 121.,   0.,   0.,   0.
     *, 350., .32650, 340.,-340.,   0.,   0.,   0.
     *, 350., .34100, 340.,-341.,   0.,   0.,   0.
     *, 350., .35550, 341.,-340.,   0.,   0.,   0.
     *, 350., .37000, 341.,-341.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=253,270)/
     *  350., .39800, 340., 120., 110.,   0.,   0.
     *, 350., .42600, 340., 120., 220.,   0.,   0.
     *, 350., .45400, 340., 120., 111.,   0.,   0.
     *, 350., .48200, 340., 120., 221.,   0.,   0.
     *, 350., .51000, 340., 121., 110.,   0.,   0.
     *, 350., .53800, 340., 121., 220.,   0.,   0.
     *, 350., .56600, 340., 121., 111.,   0.,   0.
     *, 350., .59400, 340., 121., 221.,   0.,   0.
     *, 350., .62200, 341., 120., 110.,   0.,   0.
     *, 350., .65000, 341., 120., 220.,   0.,   0.
     *, 350., .67800, 341., 120., 111.,   0.,   0.
     *, 350., .70600, 341., 120., 221.,   0.,   0.
     *, 350., .73400, 341., 121., 110.,   0.,   0.
     *, 350., .76200, 341., 121., 220.,   0.,   0.
     *, 350., .79000, 341., 121., 111.,   0.,   0.
     *, 350., .81800, 341., 121., 221.,   0.,   0.
     *, 350., .83200, 340., 130.,-230.,   0.,   0.
     *, 350., .84600, 340., 130.,-231.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=271,288)/
     *  350., .86000, 340., 131.,-230.,   0.,   0.
     *, 350., .87400, 340., 131.,-231.,   0.,   0.
     *, 350., .88800, 341., 130.,-230.,   0.,   0.
     *, 350., .90200, 341., 130.,-231.,   0.,   0.
     *, 350., .91600, 341., 131.,-230.,   0.,   0.
     *, 350., .93000, 341., 131.,-231.,   0.,   0.
     *, 350., .93300, 340.,-140., 130.,   0.,   0.
     *, 350., .93600, 340.,-140., 131.,   0.,   0.
     *, 350., .93900, 340.,-141., 130.,   0.,   0.
     *, 350., .94200, 340.,-141., 131.,   0.,   0.
     *, 350., .94500, 341.,-140., 130.,   0.,   0.
     *, 350., .94800, 341.,-140., 131.,   0.,   0.
     *, 350., .95100, 341.,-141., 130.,   0.,   0.
     *, 350., .95400, 341.,-141., 131.,   0.,   0.
     *, 350., .95700, 340.,-240., 230.,   0.,   0.
     *, 350., .96000, 340.,-240., 231.,   0.,   0.
     *, 350., .96300, 340.,-241., 230.,   0.,   0.
     *, 350., .96600, 340.,-241., 231.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=289,306)/
     *  350., .96900, 341.,-240., 230.,   0.,   0.
     *, 350., .97200, 341.,-240., 231.,   0.,   0.
     *, 350., .97500, 341.,-241., 230.,   0.,   0.
     *, 350., .97800, 341.,-241., 231.,   0.,   0.
     *, 350., .97950, 340.,-340., 330.,   0.,   0.
     *, 350., .98100, 340.,-340., 331.,   0.,   0.
     *, 350., .98250, 340.,-341., 331.,   0.,   0.
     *, 350., .98400, 340.,-341., 331.,   0.,   0.
     *, 350., .98550, 341.,-340., 330.,   0.,   0.
     *, 350., .98700, 341.,-340., 331.,   0.,   0.
     *, 350., .98850, 341.,-341., 331.,   0.,   0.
     *, 350., .99000, 341.,-341., 331.,   0.,   0.
     *, 350., .99200, 441., 130.,-130.,   0.,   0.
     *, 350., .99400, 441., 131.,-130.,   0.,   0.
     *, 350., .99600, 441., 230.,-230.,   0.,   0.
     *, 350., .99800, 441., 231.,-230.,   0.,   0.
     *, 350., .99900, 441., 330., 330.,   0.,   0.
     *, 350.,1.00000, 441., 331., 331.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=307,324)/
     *  160., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 160., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 160., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 160., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 160.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 260., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 260., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 260., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 260., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 260.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 360., .33330,  -1.,   2.,  -5.,   0.,   0.
     *, 360., .66660,  -4.,   3.,  -5.,   0.,   0.
     *, 360., .77770, -11.,  12.,  -5.,   0.,   0.
     *, 360., .88880, -13.,  14.,  -5.,   0.,   0.
     *, 360.,1.00000, -15.,  16.,  -5.,   0.,   0.
     *, 151.,1.00000, 150.,  10.,   0.,   0.,   0.
     *, 251.,1.00000, 250.,  10.,   0.,   0.,   0.
     *, 351.,1.00000, 350.,  10.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=325,342)/
     *  161.,1.00000, 160.,  10.,   0.,   0.,   0.
     *, 261.,1.00000, 260.,  10.,   0.,   0.,   0.
     *, 361.,1.00000, 360.,  10.,   0.,   0.,   0.
     *,1230.,1.00000,2130.,  10.,   0.,   0.,   0.
     *,1111.,1.00000,1120., 120.,   0.,   0.,   0.
     *,1121., .66670,1120., 110.,   0.,   0.,   0.
     *,1121.,1.00000,1220., 120.,   0.,   0.,   0.
     *,1221., .66670,1220., 110.,   0.,   0.,   0.
     *,1221.,1.00000,1120.,-120.,   0.,   0.,   0.
     *,2221.,1.00000,1220.,-120.,   0.,   0.,   0.
     *,1131., .88000,2130., 120.,   0.,   0.,   0.
     *,1131., .94000,1130., 110.,   0.,   0.,   0.
     *,1131.,1.00000,1230., 120.,   0.,   0.,   0.
     *,1231., .88000,2130., 110.,   0.,   0.,   0.
     *,1231., .94000,1130.,-120.,   0.,   0.,   0.
     *,1231.,1.00000,2230., 120.,   0.,   0.,   0.
     *,2231., .88000,2130.,-120.,   0.,   0.,   0.
     *,2231., .94000,1230.,-120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=343,360)/
     * 2231.,1.00000,2230., 110.,   0.,   0.,   0.
     *,1331., .66670,2330., 120.,   0.,   0.,   0.
     *,1331.,1.00000,1330., 110.,   0.,   0.,   0.
     *,2331., .66670,1330.,-120.,   0.,   0.,   0.
     *,2331.,1.00000,2330., 110.,   0.,   0.,   0.
     *,  16., .18000,  12., -11.,  15.,   0.,   0.
     *,  16., .36000,  14., -13.,  15.,   0.,   0.
     *,  16., .45100,-120.,  15.,   0.,   0.,   0.
     *,  16., .66000,-121.,  15.,   0.,   0.,   0.
     *,  16., .78000, 110., 110.,-120.,  15.,   0.
     *,  16., .83600, 120.,-120.,-120.,  15.,   0.
     *,  16.,1.00000, 120., 110.,-120.,-120.,  15.
     *,2140., .03750, -12.,  11.,2130.,   0.,   0.
     *,2140., .07500, -12.,  11.,1231.,   0.,   0.
     *,2140., .11250, -14.,  13.,2130.,   0.,   0.
     *,2140., .15000, -14.,  13.,1231.,   0.,   0.
     *,2140., .18200,2130., 120.,   0.,   0.,   0.
     *,2140., .21300,1230., 110.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=361,378)/
     * 2140., .24400,1120.,-230.,   0.,   0.,   0.
     *,2140., .29500,1131., 110.,   0.,   0.,   0.
     *,2140., .34600,1231., 120.,   0.,   0.,   0.
     *,2140., .39700,1121.,-230.,   0.,   0.,   0.
     *,2140., .44800,1111.,-130.,   0.,   0.,   0.
     *,2140., .49900,1130., 111.,   0.,   0.,   0.
     *,2140., .55000,1230., 121.,   0.,   0.,   0.
     *,2140., .60100,1120.,-231.,   0.,   0.,   0.
     *,2140., .65800,1120.,-230., 120.,-120.,   0.
     *,2140., .71500,1120.,-230., 110., 110.,   0.
     *,2140., .77200,1120.,-130., 120., 110.,   0.
     *,2140., .82900,1220.,-230., 120., 110.,   0.
     *,2140., .88600,1220.,-130., 120., 120.,   0.
     *,2140., .94300,2130., 120., 120.,-120.,   0.
     *,2140.,1.00000,2130., 120., 110., 110.,   0.
     *,1140.,1.00000,2140., 120.,   0.,   0.,   0.
     *,1240.,1.00000,2140., 110.,   0.,   0.,   0.
     *,2240.,1.00000,2140.,-120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=379,396)/
     * 1340., .03750, -12.,  11.,1330.,   0.,   0.
     *,1340., .07500, -12.,  11.,1331.,   0.,   0.
     *,1340., .11250, -14.,  13.,1330.,   0.,   0.
     *,1340., .15000, -14.,  13.,1331.,   0.,   0.
     *,1340., .19900,1330., 120.,   0.,   0.,   0.
     *,1340., .24800,1231.,-230.,   0.,   0.,   0.
     *,1340., .28800,1330., 120.,   0.,   0.,   0.
     *,1340., .32800,1131.,-230.,   0.,   0.,   0.
     *,1340., .36800,1330., 121.,   0.,   0.,   0.
     *,1340., .40800,1130.,-230.,   0.,   0.,   0.
     *,1340., .44800,1330., 120., 110.,   0.,   0.
     *,1340., .48800,2330., 120., 120.,   0.,   0.
     *,1340., .52800,1130.,-130., 120.,   0.,   0.
     *,1340., .56800,1130.,-230., 110.,   0.,   0.
     *,1340., .60800,1230.,-230., 120.,   0.,   0.
     *,1340., .66400,2130.,-230., 120., 110.,   0.
     *,1340., .72000,2130.,-130., 120., 120.,   0.
     *,1340., .77600,1130.,-230., 120., 120.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=397,414)/
     * 1340., .83200,1130.,-230., 110., 110.,   0.
     *,1340., .88800,1330., 120., 120.,-120.,   0.
     *,1340., .94400,1330., 120., 110., 110.,   0.
     *,1340.,1.00000,2330., 120., 120., 110.,   0.
     *,3140., .03750, -12.,  11.,1330.,   0.,   0.
     *,3140., .07500, -12.,  11.,1331.,   0.,   0.
     *,3140., .11250, -14.,  13.,1330.,   0.,   0.
     *,3140., .15000, -14.,  13.,1331.,   0.,   0.
     *,3140., .19900,1330., 120.,   0.,   0.,   0.
     *,3140., .24800,1231.,-230.,   0.,   0.,   0.
     *,3140., .28800,1330., 120.,   0.,   0.,   0.
     *,3140., .32800,1131.,-230.,   0.,   0.,   0.
     *,3140., .36800,1330., 121.,   0.,   0.,   0.
     *,3140., .40800,1130.,-230.,   0.,   0.,   0.
     *,3140., .44800,1330., 120., 110.,   0.,   0.
     *,3140., .48800,2330., 120., 120.,   0.,   0.
     *,3140., .52800,1130.,-130., 120.,   0.,   0.
     *,3140., .56800,1130.,-230., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=415,432)/
     * 3140., .60800,1230.,-230., 120.,   0.,   0.
     *,3140., .66400,2130.,-230., 120., 110.,   0.
     *,3140., .72000,2130.,-130., 120., 120.,   0.
     *,3140., .77600,1130.,-230., 120., 120.,   0.
     *,3140., .83200,1130.,-230., 110., 110.,   0.
     *,3140., .88800,1330., 120., 120.,-120.,   0.
     *,3140., .94400,1330., 120., 110., 110.,   0.
     *,3140.,1.00000,2330., 120., 120., 110.,   0.
     *,2340., .03750, -12.,  11.,2330.,   0.,   0.
     *,2340., .07500, -12.,  11.,2331.,   0.,   0.
     *,2340., .11250, -14.,  13.,2330.,   0.,   0.
     *,2340., .15000, -14.,  13.,2331.,   0.,   0.
     *,2340., .17500,2330., 120.,   0.,   0.,   0.
     *,2340., .20000,1330., 110.,   0.,   0.,   0.
     *,2340., .22500,1130.,-130.,   0.,   0.,   0.
     *,2340., .25000,1230.,-230.,   0.,   0.,   0.
     *,2340., .29500,2331., 120.,   0.,   0.,   0.
     *,2340., .34000,1331., 110.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=433,450)/
     * 2340., .38500,1131.,-130.,   0.,   0.,   0.
     *,2340., .43000,1231.,-230.,   0.,   0.,   0.
     *,2340., .47500,2330., 121.,   0.,   0.,   0.
     *,2340., .52000,1330., 111.,   0.,   0.,   0.
     *,2340., .56500,1130.,-131.,   0.,   0.,   0.
     *,2340., .61000,1230.,-231.,   0.,   0.,   0.
     *,2340., .64900,2130.,-230., 120.,-120.,   0.
     *,2340., .68800,2130.,-230., 110., 110.,   0.
     *,2340., .72700,2130.,-130., 120., 110.,   0.
     *,2340., .76600,1130.,-230.,-120., 110.,   0.
     *,2340., .80500,1130.,-130., 120.,-120.,   0.
     *,2340., .84400,1130.,-130., 110., 110.,   0.
     *,2340., .88300,1330., 120.,-120., 110.,   0.
     *,2340., .92200,1330., 110., 110., 110.,   0.
     *,2340., .96100,2330., 120., 120.,-120.,   0.
     *,2340.,1.00000,2330., 120., 110., 110.,   0.
     *,3240., .03750, -12.,  11.,2330.,   0.,   0.
     *,3240., .07500, -12.,  11.,2331.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=451,468)/
     * 3240., .11250, -14.,  13.,2330.,   0.,   0.
     *,3240., .15000, -14.,  13.,2331.,   0.,   0.
     *,3240., .17500,2330., 120.,   0.,   0.,   0.
     *,3240., .20000,1330., 110.,   0.,   0.,   0.
     *,3240., .22500,1130.,-130.,   0.,   0.,   0.
     *,3240., .25000,1230.,-230.,   0.,   0.,   0.
     *,3240., .29500,2331., 120.,   0.,   0.,   0.
     *,3240., .34000,1331., 110.,   0.,   0.,   0.
     *,3240., .38500,1131.,-130.,   0.,   0.,   0.
     *,3240., .43000,1231.,-230.,   0.,   0.,   0.
     *,3240., .47500,2330., 121.,   0.,   0.,   0.
     *,3240., .52000,1330., 111.,   0.,   0.,   0.
     *,3240., .56500,1130.,-131.,   0.,   0.,   0.
     *,3240., .61000,1230.,-231.,   0.,   0.,   0.
     *,3240., .64900,2130.,-230., 120.,-120.,   0.
     *,3240., .68800,2130.,-230., 110., 110.,   0.
     *,3240., .72700,2130.,-130., 120., 110.,   0.
     *,3240., .76600,1130.,-230.,-120., 110.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=469,486)/
     * 3240., .80500,1130.,-130., 120.,-120.,   0.
     *,3240., .84400,1130.,-130., 110., 110.,   0.
     *,3240., .88300,1330., 120.,-120., 110.,   0.
     *,3240., .92200,1330., 110., 110., 110.,   0.
     *,3240., .96100,2330., 120., 120.,-120.,   0.
     *,3240.,1.00000,2330., 120., 110., 110.,   0.
     *,3340., .07500, -12.,  11.,3331.,   0.,   0.
     *,3340., .15000, -12.,  11.,3331.,   0.,   0.
     *,3340., .25000,1330.,-230.,   0.,   0.,   0.
     *,3340., .31000,3331., 120.,   0.,   0.,   0.
     *,3340., .37000,1331.,-230.,   0.,   0.,   0.
     *,3340., .43000,1330.,-231.,   0.,   0.,   0.
     *,3340., .49000,2330.,-230., 120.,   0.,   0.
     *,3340., .55000,1330.,-230., 110.,   0.,   0.
     *,3340., .61000,1330.,-130., 120.,   0.,   0.
     *,3340., .67500,3331., 120., 120.,-120.,   0.
     *,3340., .74000,3331., 120., 110., 110.,   0.
     *,3340., .80500,1330.,-230., 120.,-120.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=487,504)/
     * 3340., .87000,1330.,-230., 110., 110.,   0.
     *,3340., .93500,2330.,-230., 120., 110.,   0.
     *,3340.,1.00000,2330.,-130., 120., 120.,   0.
     *,1141.,1.00000,2140., 120.,   0.,   0.,   0.
     *,1241.,1.00000,2140., 110.,   0.,   0.,   0.
     *,2241.,1.00000,2140.,-120.,   0.,   0.,   0.
     *,1341., .66670,2340., 120.,   0.,   0.,   0.
     *,1341.,1.00000,1340., 110.,   0.,   0.,   0.
     *,2341., .66670,1340.,-120.,   0.,   0.,   0.
     *,2341.,1.00000,2340., 110.,   0.,   0.,   0.
     *,3341.,1.00000,3340., 110.,   0.,   0.,   0.
     *,1150., .06000,  12., -11.,1140.,   0.,   0.
     *,1150., .12000,  12., -11.,1141.,   0.,   0.
     *,1150., .18000,  14., -13.,1140.,   0.,   0.
     *,1150., .24000,  14., -13.,1141.,   0.,   0.
     *,1150., .25500,  16., -15.,1140.,   0.,   0.
     *,1150., .27000,  16., -15.,1141.,   0.,   0.
     *,1150., .28925,1140.,-120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=505,522)/
     * 1150., .30850,1140.,-121.,   0.,   0.,   0.
     *,1150., .32775,1141.,-120.,   0.,   0.,   0.
     *,1150., .34700,1141.,-121.,   0.,   0.,   0.
     *,1150., .35775,1140., 340.,   0.,   0.,   0.
     *,1150., .36850,1140., 341.,   0.,   0.,   0.
     *,1150., .37925,1141., 340.,   0.,   0.,   0.
     *,1150., .39000,1141., 341.,   0.,   0.,   0.
     *,1150., .42050,1140.,-120., 110.,   0.,   0.
     *,1150., .45100,1140.,-120., 220.,   0.,   0.
     *,1150., .48150,1140.,-120., 111.,   0.,   0.
     *,1150., .51200,1140.,-120., 221.,   0.,   0.
     *,1150., .54250,1140.,-121., 110.,   0.,   0.
     *,1150., .57300,1140.,-121., 220.,   0.,   0.
     *,1150., .60350,1140.,-121., 111.,   0.,   0.
     *,1150., .63400,1140.,-121., 221.,   0.,   0.
     *,1150., .66450,1141.,-120., 110.,   0.,   0.
     *,1150., .69500,1141.,-120., 220.,   0.,   0.
     *,1150., .72550,1141.,-120., 111.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=523,540)/
     * 1150., .75600,1141.,-120., 221.,   0.,   0.
     *,1150., .78650,1141.,-121., 110.,   0.,   0.
     *,1150., .81700,1141.,-121., 220.,   0.,   0.
     *,1150., .84750,1141.,-121., 111.,   0.,   0.
     *,1150., .87800,1141.,-121., 221.,   0.,   0.
     *,1150., .89325,1140.,-130., 230.,   0.,   0.
     *,1150., .90850,1140.,-130., 231.,   0.,   0.
     *,1150., .92375,1140.,-131., 230.,   0.,   0.
     *,1150., .93900,1140.,-131., 231.,   0.,   0.
     *,1150., .95425,1141.,-130., 230.,   0.,   0.
     *,1150., .96950,1141.,-130., 231.,   0.,   0.
     *,1150., .98475,1141.,-131., 230.,   0.,   0.
     *,1150.,1.00000,1141.,-131., 231.,   0.,   0.
     *,1250., .06000,  12., -11.,1240.,   0.,   0.
     *,1250., .12000,  12., -11.,1241.,   0.,   0.
     *,1250., .18000,  14., -13.,1240.,   0.,   0.
     *,1250., .24000,  14., -13.,1241.,   0.,   0.
     *,1250., .25500,  16., -15.,1240.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=541,558)/
     * 1250., .27000,  16., -15.,1241.,   0.,   0.
     *,1250., .28925,1240.,-120.,   0.,   0.,   0.
     *,1250., .30850,1240.,-121.,   0.,   0.,   0.
     *,1250., .32775,1241.,-120.,   0.,   0.,   0.
     *,1250., .34700,1241.,-121.,   0.,   0.,   0.
     *,1250., .35775,1240., 340.,   0.,   0.,   0.
     *,1250., .36850,1240., 341.,   0.,   0.,   0.
     *,1250., .37925,1241., 340.,   0.,   0.,   0.
     *,1250., .39000,1241., 341.,   0.,   0.,   0.
     *,1250., .42050,1240.,-120., 110.,   0.,   0.
     *,1250., .45100,1240.,-120., 220.,   0.,   0.
     *,1250., .48150,1240.,-120., 111.,   0.,   0.
     *,1250., .51200,1240.,-120., 221.,   0.,   0.
     *,1250., .54250,1240.,-121., 110.,   0.,   0.
     *,1250., .57300,1240.,-121., 220.,   0.,   0.
     *,1250., .60350,1240.,-121., 111.,   0.,   0.
     *,1250., .63400,1240.,-121., 221.,   0.,   0.
     *,1250., .66450,1241.,-120., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=559,576)/
     * 1250., .69500,1241.,-120., 220.,   0.,   0.
     *,1250., .72550,1241.,-120., 111.,   0.,   0.
     *,1250., .75600,1241.,-120., 221.,   0.,   0.
     *,1250., .78650,1241.,-121., 110.,   0.,   0.
     *,1250., .81700,1241.,-121., 220.,   0.,   0.
     *,1250., .84750,1241.,-121., 111.,   0.,   0.
     *,1250., .87800,1241.,-121., 221.,   0.,   0.
     *,1250., .89325,1240.,-130., 230.,   0.,   0.
     *,1250., .90850,1240.,-130., 231.,   0.,   0.
     *,1250., .92375,1240.,-131., 230.,   0.,   0.
     *,1250., .93900,1240.,-131., 231.,   0.,   0.
     *,1250., .95425,1241.,-130., 230.,   0.,   0.
     *,1250., .96950,1241.,-130., 231.,   0.,   0.
     *,1250., .98475,1241.,-131., 230.,   0.,   0.
     *,1250.,1.00000,1241.,-131., 231.,   0.,   0.
     *,1350., .06000,  12., -11.,1340.,   0.,   0.
     *,1350., .12000,  12., -11.,1341.,   0.,   0.
     *,1350., .18000,  14., -13.,1340.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=577,594)/
     * 1350., .24000,  14., -13.,1341.,   0.,   0.
     *,1350., .25500,  16., -15.,1340.,   0.,   0.
     *,1350., .27000,  16., -15.,1341.,   0.,   0.
     *,1350., .28925,1340.,-120.,   0.,   0.,   0.
     *,1350., .30850,1340.,-121.,   0.,   0.,   0.
     *,1350., .32775,1341.,-120.,   0.,   0.,   0.
     *,1350., .34700,1341.,-121.,   0.,   0.,   0.
     *,1350., .35775,1340., 340.,   0.,   0.,   0.
     *,1350., .36850,1340., 341.,   0.,   0.,   0.
     *,1350., .37925,1341., 340.,   0.,   0.,   0.
     *,1350., .39000,1341., 341.,   0.,   0.,   0.
     *,1350., .42050,1340.,-120., 110.,   0.,   0.
     *,1350., .45100,1340.,-120., 220.,   0.,   0.
     *,1350., .48150,1340.,-120., 111.,   0.,   0.
     *,1350., .51200,1340.,-120., 221.,   0.,   0.
     *,1350., .54250,1340.,-121., 110.,   0.,   0.
     *,1350., .57300,1340.,-121., 220.,   0.,   0.
     *,1350., .60350,1340.,-121., 111.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=595,612)/
     * 1350., .63400,1340.,-121., 221.,   0.,   0.
     *,1350., .66450,1341.,-120., 110.,   0.,   0.
     *,1350., .69500,1341.,-120., 220.,   0.,   0.
     *,1350., .72550,1341.,-120., 111.,   0.,   0.
     *,1350., .75600,1341.,-120., 221.,   0.,   0.
     *,1350., .78650,1341.,-121., 110.,   0.,   0.
     *,1350., .81700,1341.,-121., 220.,   0.,   0.
     *,1350., .84750,1341.,-121., 111.,   0.,   0.
     *,1350., .87800,1341.,-121., 221.,   0.,   0.
     *,1350., .89325,1340.,-130., 230.,   0.,   0.
     *,1350., .90850,1340.,-130., 231.,   0.,   0.
     *,1350., .92375,1340.,-131., 230.,   0.,   0.
     *,1350., .93900,1340.,-131., 231.,   0.,   0.
     *,1350., .95425,1341.,-130., 230.,   0.,   0.
     *,1350., .96950,1341.,-130., 231.,   0.,   0.
     *,1350., .98475,1341.,-131., 230.,   0.,   0.
     *,1350.,1.00000,1341.,-131., 231.,   0.,   0.
     *,2150., .06000,  12., -11.,2140.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=613,630)/
     * 2150., .12000,  12., -11.,1241.,   0.,   0.
     *,2150., .18000,  14., -13.,2140.,   0.,   0.
     *,2150., .24000,  14., -13.,1241.,   0.,   0.
     *,2150., .25500,  16., -15.,2140.,   0.,   0.
     *,2150., .27000,  16., -15.,1241.,   0.,   0.
     *,2150., .28925,2140.,-120.,   0.,   0.,   0.
     *,2150., .30850,2140.,-121.,   0.,   0.,   0.
     *,2150., .32775,1241.,-120.,   0.,   0.,   0.
     *,2150., .34700,1241.,-121.,   0.,   0.,   0.
     *,2150., .35775,2140., 340.,   0.,   0.,   0.
     *,2150., .36850,2140., 341.,   0.,   0.,   0.
     *,2150., .37925,1241., 340.,   0.,   0.,   0.
     *,2150., .39000,1241., 341.,   0.,   0.,   0.
     *,2150., .42050,2140.,-120., 110.,   0.,   0.
     *,2150., .45100,2140.,-120., 220.,   0.,   0.
     *,2150., .48150,2140.,-120., 111.,   0.,   0.
     *,2150., .51200,2140.,-120., 221.,   0.,   0.
     *,2150., .54250,2140.,-121., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=631,648)/
     * 2150., .57300,2140.,-121., 220.,   0.,   0.
     *,2150., .60350,2140.,-121., 111.,   0.,   0.
     *,2150., .63400,2140.,-121., 221.,   0.,   0.
     *,2150., .66450,1241.,-120., 110.,   0.,   0.
     *,2150., .69500,1241.,-120., 220.,   0.,   0.
     *,2150., .72550,1241.,-120., 111.,   0.,   0.
     *,2150., .75600,1241.,-120., 221.,   0.,   0.
     *,2150., .78650,1241.,-121., 110.,   0.,   0.
     *,2150., .81700,1241.,-121., 220.,   0.,   0.
     *,2150., .84750,1241.,-121., 111.,   0.,   0.
     *,2150., .87800,1241.,-121., 221.,   0.,   0.
     *,2150., .89325,2140.,-130., 230.,   0.,   0.
     *,2150., .90850,2140.,-130., 231.,   0.,   0.
     *,2150., .92375,2140.,-131., 230.,   0.,   0.
     *,2150., .93900,2140.,-131., 231.,   0.,   0.
     *,2150., .95425,1241.,-130., 230.,   0.,   0.
     *,2150., .96950,1241.,-130., 231.,   0.,   0.
     *,2150., .98475,1241.,-131., 230.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=649,666)/
     * 2150.,1.00000,1241.,-131., 231.,   0.,   0.
     *,2250., .06000,  12., -11.,2240.,   0.,   0.
     *,2250., .12000,  12., -11.,2241.,   0.,   0.
     *,2250., .18000,  14., -13.,2240.,   0.,   0.
     *,2250., .24000,  14., -13.,2241.,   0.,   0.
     *,2250., .25500,  16., -15.,2240.,   0.,   0.
     *,2250., .27000,  16., -15.,2241.,   0.,   0.
     *,2250., .28925,2240.,-120.,   0.,   0.,   0.
     *,2250., .30850,2240.,-121.,   0.,   0.,   0.
     *,2250., .32775,2241.,-120.,   0.,   0.,   0.
     *,2250., .34700,2241.,-121.,   0.,   0.,   0.
     *,2250., .35775,2240., 340.,   0.,   0.,   0.
     *,2250., .36850,2240., 341.,   0.,   0.,   0.
     *,2250., .37925,2241., 340.,   0.,   0.,   0.
     *,2250., .39000,2241., 341.,   0.,   0.,   0.
     *,2250., .42050,2240.,-120., 110.,   0.,   0.
     *,2250., .45100,2240.,-120., 220.,   0.,   0.
     *,2250., .48150,2240.,-120., 111.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=667,684)/
     * 2250., .51200,2240.,-120., 221.,   0.,   0.
     *,2250., .54250,2240.,-121., 110.,   0.,   0.
     *,2250., .57300,2240.,-121., 220.,   0.,   0.
     *,2250., .60350,2240.,-121., 111.,   0.,   0.
     *,2250., .63400,2240.,-121., 221.,   0.,   0.
     *,2250., .66450,2241.,-120., 110.,   0.,   0.
     *,2250., .69500,2241.,-120., 220.,   0.,   0.
     *,2250., .72550,2241.,-120., 111.,   0.,   0.
     *,2250., .75600,2241.,-120., 221.,   0.,   0.
     *,2250., .78650,2241.,-121., 110.,   0.,   0.
     *,2250., .81700,2241.,-121., 220.,   0.,   0.
     *,2250., .84750,2241.,-121., 111.,   0.,   0.
     *,2250., .87800,2241.,-121., 221.,   0.,   0.
     *,2250., .89325,2240.,-130., 230.,   0.,   0.
     *,2250., .90850,2240.,-130., 231.,   0.,   0.
     *,2250., .92375,2240.,-131., 230.,   0.,   0.
     *,2250., .93900,2240.,-131., 231.,   0.,   0.
     *,2250., .95425,2241.,-130., 230.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=685,702)/
     * 2250., .96950,2241.,-130., 231.,   0.,   0.
     *,2250., .98475,2241.,-131., 230.,   0.,   0.
     *,2250.,1.00000,2241.,-131., 231.,   0.,   0.
     *,2350., .06000,  12., -11.,2340.,   0.,   0.
     *,2350., .12000,  12., -11.,2341.,   0.,   0.
     *,2350., .18000,  14., -13.,2340.,   0.,   0.
     *,2350., .24000,  14., -13.,2341.,   0.,   0.
     *,2350., .25500,  16., -15.,2340.,   0.,   0.
     *,2350., .27000,  16., -15.,2341.,   0.,   0.
     *,2350., .28925,2340.,-120.,   0.,   0.,   0.
     *,2350., .30850,2340.,-121.,   0.,   0.,   0.
     *,2350., .32775,2341.,-120.,   0.,   0.,   0.
     *,2350., .34700,2341.,-121.,   0.,   0.,   0.
     *,2350., .35775,2340., 340.,   0.,   0.,   0.
     *,2350., .36850,2340., 341.,   0.,   0.,   0.
     *,2350., .37925,2341., 340.,   0.,   0.,   0.
     *,2350., .39000,2341., 341.,   0.,   0.,   0.
     *,2350., .42050,2340.,-120., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=703,720)/
     * 2350., .45100,2340.,-120., 220.,   0.,   0.
     *,2350., .48150,2340.,-120., 111.,   0.,   0.
     *,2350., .51200,2340.,-120., 221.,   0.,   0.
     *,2350., .54250,2340.,-121., 110.,   0.,   0.
     *,2350., .57300,2340.,-121., 220.,   0.,   0.
     *,2350., .60350,2340.,-121., 111.,   0.,   0.
     *,2350., .63400,2340.,-121., 221.,   0.,   0.
     *,2350., .66450,2341.,-120., 110.,   0.,   0.
     *,2350., .69500,2341.,-120., 220.,   0.,   0.
     *,2350., .72550,2341.,-120., 111.,   0.,   0.
     *,2350., .75600,2341.,-120., 221.,   0.,   0.
     *,2350., .78650,2341.,-121., 110.,   0.,   0.
     *,2350., .81700,2341.,-121., 220.,   0.,   0.
     *,2350., .84750,2341.,-121., 111.,   0.,   0.
     *,2350., .87800,2341.,-121., 221.,   0.,   0.
     *,2350., .89325,2340.,-130., 230.,   0.,   0.
     *,2350., .90850,2340.,-130., 231.,   0.,   0.
     *,2350., .92375,2340.,-131., 230.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=721,738)/
     * 2350., .93900,2340.,-131., 231.,   0.,   0.
     *,2350., .95425,2341.,-130., 230.,   0.,   0.
     *,2350., .96950,2341.,-130., 231.,   0.,   0.
     *,2350., .98475,2341.,-131., 230.,   0.,   0.
     *,2350.,1.00000,2341.,-131., 231.,   0.,   0.
     *,3150., .06000,  12., -11.,3140.,   0.,   0.
     *,3150., .12000,  12., -11.,1341.,   0.,   0.
     *,3150., .18000,  14., -13.,3140.,   0.,   0.
     *,3150., .24000,  14., -13.,1341.,   0.,   0.
     *,3150., .25500,  16., -15.,3140.,   0.,   0.
     *,3150., .27000,  16., -15.,1341.,   0.,   0.
     *,3150., .28925,3140.,-120.,   0.,   0.,   0.
     *,3150., .30850,3140.,-121.,   0.,   0.,   0.
     *,3150., .32775,1341.,-120.,   0.,   0.,   0.
     *,3150., .34700,1341.,-121.,   0.,   0.,   0.
     *,3150., .35775,3140., 340.,   0.,   0.,   0.
     *,3150., .36850,3140., 341.,   0.,   0.,   0.
     *,3150., .37925,1341., 340.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=739,756)/
     * 3150., .39000,1341., 341.,   0.,   0.,   0.
     *,3150., .42050,3140.,-120., 110.,   0.,   0.
     *,3150., .45100,3140.,-120., 220.,   0.,   0.
     *,3150., .48150,3140.,-120., 111.,   0.,   0.
     *,3150., .51200,3140.,-120., 221.,   0.,   0.
     *,3150., .54250,3140.,-121., 110.,   0.,   0.
     *,3150., .57300,3140.,-121., 220.,   0.,   0.
     *,3150., .60350,3140.,-121., 111.,   0.,   0.
     *,3150., .63400,3140.,-121., 221.,   0.,   0.
     *,3150., .66450,1341.,-120., 110.,   0.,   0.
     *,3150., .69500,1341.,-120., 220.,   0.,   0.
     *,3150., .72550,1341.,-120., 111.,   0.,   0.
     *,3150., .75600,1341.,-120., 221.,   0.,   0.
     *,3150., .78650,1341.,-121., 110.,   0.,   0.
     *,3150., .81700,1341.,-121., 220.,   0.,   0.
     *,3150., .84750,1341.,-121., 111.,   0.,   0.
     *,3150., .87800,1341.,-121., 221.,   0.,   0.
     *,3150., .89325,3140.,-130., 230.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=757,774)/
     * 3150., .90850,3140.,-130., 231.,   0.,   0.
     *,3150., .92375,3140.,-131., 230.,   0.,   0.
     *,3150., .93900,3140.,-131., 231.,   0.,   0.
     *,3150., .95425,1341.,-130., 230.,   0.,   0.
     *,3150., .96950,1341.,-130., 231.,   0.,   0.
     *,3150., .98475,1341.,-131., 230.,   0.,   0.
     *,3150.,1.00000,1341.,-131., 231.,   0.,   0.
     *,3250., .06000,  12., -11.,3240.,   0.,   0.
     *,3250., .12000,  12., -11.,2341.,   0.,   0.
     *,3250., .18000,  14., -13.,3240.,   0.,   0.
     *,3250., .24000,  14., -13.,2341.,   0.,   0.
     *,3250., .25500,  16., -15.,3240.,   0.,   0.
     *,3250., .27000,  16., -15.,2341.,   0.,   0.
     *,3250., .28925,3240.,-120.,   0.,   0.,   0.
     *,3250., .30850,3240.,-121.,   0.,   0.,   0.
     *,3250., .32775,2341.,-120.,   0.,   0.,   0.
     *,3250., .34700,2341.,-121.,   0.,   0.,   0.
     *,3250., .35775,3240., 340.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=775,792)/
     * 3250., .36850,3240., 341.,   0.,   0.,   0.
     *,3250., .37925,2341., 340.,   0.,   0.,   0.
     *,3250., .39000,2341., 341.,   0.,   0.,   0.
     *,3250., .42050,3240.,-120., 110.,   0.,   0.
     *,3250., .45100,3240.,-120., 220.,   0.,   0.
     *,3250., .48150,3240.,-120., 111.,   0.,   0.
     *,3250., .51200,3240.,-120., 221.,   0.,   0.
     *,3250., .54250,3240.,-121., 110.,   0.,   0.
     *,3250., .57300,3240.,-121., 220.,   0.,   0.
     *,3250., .60350,3240.,-121., 111.,   0.,   0.
     *,3250., .63400,3240.,-121., 221.,   0.,   0.
     *,3250., .66450,2341.,-120., 110.,   0.,   0.
     *,3250., .69500,2341.,-120., 220.,   0.,   0.
     *,3250., .72550,2341.,-120., 111.,   0.,   0.
     *,3250., .75600,2341.,-120., 221.,   0.,   0.
     *,3250., .78650,2341.,-121., 110.,   0.,   0.
     *,3250., .81700,2341.,-121., 220.,   0.,   0.
     *,3250., .84750,2341.,-121., 111.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=793,810)/
     * 3250., .87800,2341.,-121., 221.,   0.,   0.
     *,3250., .89325,3240.,-130., 230.,   0.,   0.
     *,3250., .90850,3240.,-130., 231.,   0.,   0.
     *,3250., .92375,3240.,-131., 230.,   0.,   0.
     *,3250., .93900,3240.,-131., 231.,   0.,   0.
     *,3250., .95425,2341.,-130., 230.,   0.,   0.
     *,3250., .96950,2341.,-130., 231.,   0.,   0.
     *,3250., .98475,2341.,-131., 230.,   0.,   0.
     *,3250.,1.00000,2341.,-131., 231.,   0.,   0.
     *,3350., .06000,  12., -11.,3340.,   0.,   0.
     *,3350., .12000,  12., -11.,3341.,   0.,   0.
     *,3350., .18000,  14., -13.,3340.,   0.,   0.
     *,3350., .24000,  14., -13.,3341.,   0.,   0.
     *,3350., .25500,  16., -15.,3340.,   0.,   0.
     *,3350., .27000,  16., -15.,3341.,   0.,   0.
     *,3350., .28925,3340.,-120.,   0.,   0.,   0.
     *,3350., .30850,3340.,-121.,   0.,   0.,   0.
     *,3350., .32775,3341.,-120.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=811,828)/
     * 3350., .34700,3341.,-121.,   0.,   0.,   0.
     *,3350., .35775,3340., 340.,   0.,   0.,   0.
     *,3350., .36850,3340., 341.,   0.,   0.,   0.
     *,3350., .37925,3341., 340.,   0.,   0.,   0.
     *,3350., .39000,3341., 341.,   0.,   0.,   0.
     *,3350., .42050,3340.,-120., 110.,   0.,   0.
     *,3350., .45100,3340.,-120., 220.,   0.,   0.
     *,3350., .48150,3340.,-120., 111.,   0.,   0.
     *,3350., .51200,3340.,-120., 221.,   0.,   0.
     *,3350., .54250,3340.,-121., 110.,   0.,   0.
     *,3350., .57300,3340.,-121., 220.,   0.,   0.
     *,3350., .60350,3340.,-121., 111.,   0.,   0.
     *,3350., .63400,3340.,-121., 221.,   0.,   0.
     *,3350., .66450,3341.,-120., 110.,   0.,   0.
     *,3350., .69500,3341.,-120., 220.,   0.,   0.
     *,3350., .72550,3341.,-120., 111.,   0.,   0.
     *,3350., .75600,3341.,-120., 221.,   0.,   0.
     *,3350., .78650,3341.,-121., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=829,846)/
     * 3350., .81700,3341.,-121., 220.,   0.,   0.
     *,3350., .84750,3341.,-121., 111.,   0.,   0.
     *,3350., .87800,3341.,-121., 221.,   0.,   0.
     *,3350., .89325,3340.,-130., 230.,   0.,   0.
     *,3350., .90850,3340.,-130., 231.,   0.,   0.
     *,3350., .92375,3340.,-131., 230.,   0.,   0.
     *,3350., .93900,3340.,-131., 231.,   0.,   0.
     *,3350., .95425,3341.,-130., 230.,   0.,   0.
     *,3350., .96950,3341.,-130., 231.,   0.,   0.
     *,3350., .98475,3341.,-131., 230.,   0.,   0.
     *,3350.,1.00000,3341.,-131., 231.,   0.,   0.
     *,1160., .33300,   1.,  -2.,1500.,   0.,   0.
     *,1160., .66700,   4.,  -3.,1500.,   0.,   0.
     *,1160., .77800, -12.,  11.,1500.,   0.,   0.
     *,1160., .88900, -14.,  13.,1500.,   0.,   0.
     *,1160.,1.00000, -16.,  15.,1500.,   0.,   0.
     *,1260., .33300,   1.,  -2.,2500.,   0.,   0.
     *,1260., .66700,   4.,  -3.,2500.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=847,864)/
     * 1260., .77800, -12.,  11.,2500.,   0.,   0.
     *,1260., .88900, -14.,  13.,2500.,   0.,   0.
     *,1260.,1.00000, -16.,  15.,2500.,   0.,   0.
     *,2260., .33300,   1.,  -2.,2500.,   0.,   0.
     *,2260., .66700,   4.,  -3.,2500.,   0.,   0.
     *,2260., .77800, -12.,  11.,2500.,   0.,   0.
     *,2260., .88900, -14.,  13.,2500.,   0.,   0.
     *,2260.,1.00000, -16.,  15.,2500.,   0.,   0.
     *,2160., .33300,   1.,  -2.,1500.,   0.,   0.
     *,2160., .66700,   4.,  -3.,1500.,   0.,   0.
     *,2160., .77800, -12.,  11.,1500.,   0.,   0.
     *,2160., .88900, -14.,  13.,1500.,   0.,   0.
     *,2160.,1.00000, -16.,  15.,1500.,   0.,   0.
     *,1360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,1360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,1360., .77800, -12.,  11.,3500.,   0.,   0.
     *,1360., .88900, -14.,  13.,3500.,   0.,   0.
     *,1360.,1.00000, -16.,  15.,3500.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=865,882)/
     * 2360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,2360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,2360., .77800, -12.,  11.,3500.,   0.,   0.
     *,2360., .88900, -14.,  13.,3500.,   0.,   0.
     *,2360.,1.00000, -16.,  15.,3500.,   0.,   0.
     *,3360., .33300,   1.,  -2.,3500.,   0.,   0.
     *,3360., .66700,   4.,  -3.,3500.,   0.,   0.
     *,3360., .77800, -12.,  11.,3500.,   0.,   0.
     *,3360., .88900, -14.,  13.,3500.,   0.,   0.
     *,3360.,1.00000, -16.,  15.,3500.,   0.,   0.
     *,1151.,1.00000,1150.,  10.,   0.,   0.,   0.
     *,1251.,1.00000,1250.,  10.,   0.,   0.,   0.
     *,2251.,1.00000,2250.,  10.,   0.,   0.,   0.
     *,1351.,1.00000,1350.,  10.,   0.,   0.,   0.
     *,2351.,1.00000,2350.,  10.,   0.,   0.,   0.
     *,3351.,1.00000,3350.,  10.,   0.,   0.,   0.
     *,1161.,1.00000,1160.,  10.,   0.,   0.,   0.
     *,1261.,1.00000,1260.,  10.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=883,886)/
     * 2261.,1.00000,2260.,  10.,   0.,   0.,   0.
     *,1361.,1.00000,1360.,  10.,   0.,   0.,   0.
     *,2361.,1.00000,2360.,  10.,   0.,   0.,   0.
     *,3361.,1.00000,3360.,  10.,   0.,   0.,   0./
C    *---------------------------------------------
C    *    DELTA++ RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=887,900)/
C    *--DL++(1620)---------------------------------
     * 1112., .30000,1120., 120.,   0.,   0.,   0.
     *,1112., .66000,1111., 110.,   0.,   0.,   0.
     *,1112., .90000,1121., 120.,   0.,   0.,   0.
     *,1112.,1.00000,1120., 120., 110.,   0.,   0.
C    *--DL++(1700)---------------------------------
     *,1113., .15000,1120., 120.,   0.,   0.,   0.
     *,1113., .51000,1111., 110.,   0.,   0.,   0.
     *,1113., .75000,1121., 120.,   0.,   0.,   0.
     *,1113.,1.00000,1120., 120., 110.,   0.,   0.
C    *--DL++(1925)---------------------------------
     *,1114., .28000,1120., 120.,   0.,   0.,   0.
     *,1114., .40600,1111., 110.,   0.,   0.,   0.
     *,1114., .49000,1121., 120.,   0.,   0.,   0.
     *,1114., .69000,1120., 121.,   0.,   0.,   0.
     *,1114., .70000,1130., 130.,   0.,   0.,   0.
     *,1114.,1.00000,1122., 120.,   0.,   0.,   0./
C    *---------------------------------------------
C    *    DELTA- RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=901,914)/
C    *--DL-(1620)----------------------------------
     * 2222., .30000,1220.,-120.,   0.,   0.,   0.
     *,2222., .66000,2221., 110.,   0.,   0.,   0.
     *,2222., .90000,1221.,-120.,   0.,   0.,   0.
     *,2222.,1.00000,1220., 110.,-120.,   0.,   0.
C    *--DL-(1700)----------------------------------
     *,2223., .15000,1220.,-120.,   0.,   0.,   0.
     *,2223., .51000,2221., 110.,   0.,   0.,   0.
     *,2223., .75000,1221.,-120.,   0.,   0.,   0.
     *,2223.,1.00000,1220., 110.,-120.,   0.,   0.
C    *--DL-(1925)----------------------------------
     *,2224., .28000,1220.,-120.,   0.,   0.,   0.
     *,2224., .40600,2221., 110.,   0.,   0.,   0.
     *,2224., .49000,1221.,-120.,   0.,   0.,   0.
     *,2224., .69000,1220.,-121.,   0.,   0.,   0.
     *,2224., .70000,2230., 230.,   0.,   0.,   0.
     *,2224.,1.00000,1222.,-120.,   0.,   0.,   0./
C    *---------------------------------------------
C    *    N*+ RESONANCES + DELTA+ RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=915,931)/
C    *--N*+(1440)----------------------------------
     * 1122., .20000,1120., 110.,   0.,   0.,   0.
     *,1122., .60000,1220., 120.,   0.,   0.,   0.
     *,1122., .68000,1111.,-120.,   0.,   0.,   0.
     *,1122., .73000,1121., 110.,   0.,   0.,   0.
     *,1122., .76000,1221., 120.,   0.,   0.,   0.
     *,1122., .84000,1120., 120.,-120.,   0.,   0.
     *,1122., .87000,1120., 110., 110.,   0.,   0.
     *,1122.,1.00000,1220., 120., 110.,   0.,   0.
C    *--N*+(1530)----------------------------------
     *,1123., .17000,1120., 110.,   0.,   0.,   0.
     *,1123., .51000,1220., 120.,   0.,   0.,   0.
     *,1123., .57000,1111.,-120.,   0.,   0.,   0.
     *,1123., .61000,1121., 110.,   0.,   0.,   0.
     *,1123., .63000,1221., 120.,   0.,   0.,   0.
     *,1123., .67000,1120., 120.,-120.,   0.,   0.
     *,1123., .68000,1120., 110., 110.,   0.,   0.
     *,1123., .75000,1220., 120., 110.,   0.,   0.
     *,1123.,1.00000,1120., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=932,948)/
C    *--DL+(1620)----------------------------------
     * 1124., .20000,1120., 110.,   0.,   0.,   0.
     *,1124., .30000,1220., 120.,   0.,   0.,   0.
     *,1124., .54000,1111.,-120.,   0.,   0.,   0.
     *,1124., .58000,1121., 110.,   0.,   0.,   0.
     *,1124., .90000,1221., 120.,   0.,   0.,   0.
     *,1124., .96000,1120., 120.,-120.,   0.,   0.
     *,1124.,1.00000,1220., 120., 110.,   0.,   0.
C    *--N*+(1665)----------------------------------
     *,1125., .16700,1120., 110.,   0.,   0.,   0.
     *,1125., .49970,1220., 120.,   0.,   0.,   0.
     *,1125., .62470,1111.,-120.,   0.,   0.,   0.
     *,1125., .70800,1121., 110.,   0.,   0.,   0.
     *,1125., .74970,1221., 120.,   0.,   0.,   0.
     *,1125., .82080,1120., 120.,-120.,   0.,   0.
     *,1125., .85190,1120., 110., 110.,   0.,   0.
     *,1125., .96300,1220., 120., 110.,   0.,   0.
     *,1125., .97300,1120., 220.,   0.,   0.,   0.
     *,1125.,1.00000,2130., 130.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=949,955)/
C    *--DL+(1700)----------------------------------
     * 1126., .10000,1120., 110.,   0.,   0.,   0.
     *,1126., .15000,1220., 120.,   0.,   0.,   0.
     *,1126., .39000,1111.,-120.,   0.,   0.,   0.
     *,1126., .43000,1121., 110.,   0.,   0.,   0.
     *,1126., .75000,1221., 120.,   0.,   0.,   0.
     *,1126., .91500,1120., 120.,-120.,   0.,   0.
     *,1126.,1.00000,1220., 120., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=956,969)/
C    *--N*+(1710)----------------------------------
     * 1127., .04430,1120., 110.,   0.,   0.,   0.
     *,1127., .13290,1220., 120.,   0.,   0.,   0.
     *,1127., .23790,1111.,-120.,   0.,   0.,   0.
     *,1127., .30790,1121., 110.,   0.,   0.,   0.
     *,1127., .34290,1221., 120.,   0.,   0.,   0.
     *,1127., .41190,1120., 120.,-120.,   0.,   0.
     *,1127., .48090,1120., 110., 110.,   0.,   0.
     *,1127., .54990,1220., 120., 110.,   0.,   0.
     *,1127., .66070,1120., 220.,   0.,   0.,   0.
     *,1127., .72800,2130., 130.,   0.,   0.,   0.
     *,1127., .74930,1230., 130.,   0.,   0.,   0.
     *,1127., .76000,1130., 230.,   0.,   0.,   0.
     *,1127., .84000,1120., 111.,   0.,   0.,   0.
     *,1127.,1.00000,1220., 121.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=970,980)/
C    *--DL+(1925)----------------------------------
     * 1128., .18700,1120., 110.,   0.,   0.,   0.
     *,1128., .28000,1220., 120.,   0.,   0.,   0.
     *,1128., .36400,1111.,-120.,   0.,   0.,   0.
     *,1128., .37800,1121., 110.,   0.,   0.,   0.
     *,1128., .49000,1221., 120.,   0.,   0.,   0.
     *,1128., .62300,1120., 111.,   0.,   0.,   0.
     *,1128., .69000,1220., 121.,   0.,   0.,   0.
     *,1128., .69350,1130., 230.,   0.,   0.,   0.
     *,1128., .69900,1230., 130.,   0.,   0.,   0.
     *,1128., .89900,1122., 110.,   0.,   0.,   0.
     *,1128.,1.00000,1222., 120.,   0.,   0.,   0./
C    *---------------------------------------------
C    *    N*0  RESONANCES + DELTA0 RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=981,997)/
C    *----------N*0(1440)--------------------------
     * 1222., .20000,1220., 110.,   0.,   0.,   0.
     *,1222., .60000,1120.,-120.,   0.,   0.,   0.
     *,1222., .68000,2221., 120.,   0.,   0.,   0.
     *,1222., .73000,1221., 110.,   0.,   0.,   0.
     *,1222., .76000,1121.,-120.,   0.,   0.,   0.
     *,1222., .84000,1220., 120.,-120.,   0.,   0.
     *,1222., .87000,1220., 110., 110.,   0.,   0.
     *,1222.,1.00000,1120.,-120., 110.,   0.,   0.
C    *----------N*0(1530)--------------------------
     *,1223., .17000,1220., 110.,   0.,   0.,   0.
     *,1223., .51000,1120.,-120.,   0.,   0.,   0.
     *,1223., .57000,2221., 120.,   0.,   0.,   0.
     *,1223., .61000,1221., 110.,   0.,   0.,   0.
     *,1223., .63000,1121.,-120.,   0.,   0.,   0.
     *,1223., .67000,1220., 120.,-120.,   0.,   0.
     *,1223., .68000,1220., 110., 110.,   0.,   0.
     *,1223., .75000,1120.,-120., 110.,   0.,   0.
     *,1223.,1.00000,1220., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=998,1014)/
C    *----------DL0(1620)--------------------------
     * 1224., .20000,1220., 110.,   0.,   0.,   0.
     *,1224., .30000,1120.,-120.,   0.,   0.,   0.
     *,1224., .54000,2221., 120.,   0.,   0.,   0.
     *,1224., .58000,1221., 110.,   0.,   0.,   0.
     *,1224., .90000,1121.,-120.,   0.,   0.,   0.
     *,1224., .96500,1220., 120.,-120.,   0.,   0.
     *,1224.,1.00000,1120.,-120., 110.,   0.,   0.
C    *----------N*0(1665)--------------------------
     *,1225., .16700,1220., 110.,   0.,   0.,   0.
     *,1225., .49970,1120.,-120.,   0.,   0.,   0.
     *,1225., .62470,2221., 120.,   0.,   0.,   0.
     *,1225., .70800,1221., 110.,   0.,   0.,   0.
     *,1225., .74970,1121.,-120.,   0.,   0.,   0.
     *,1225., .82080,1220., 120.,-120.,   0.,   0.
     *,1225., .85190,1220., 110., 110.,   0.,   0.
     *,1225., .96300,1120.,-120., 110.,   0.,   0.
     *,1225., .97300,1220., 220.,   0.,   0.,   0.
     *,1225.,1.00000,2130., 230.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1015,1021)/
C    *----------DL0(1700)--------------------------
     * 1226., .10000,1220., 110.,   0.,   0.,   0.
     *,1226., .15000,1120.,-120.,   0.,   0.,   0.
     *,1226., .39000,2221., 120.,   0.,   0.,   0.
     *,1226., .43000,1221., 110.,   0.,   0.,   0.
     *,1226., .75000,1121.,-120.,   0.,   0.,   0.
     *,1226., .91500,1220., 120.,-120.,   0.,   0.
     *,1226.,1.00000,1120.,-120., 110.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1022,1035)/
C    *----------N*0(1710)--------------------------
     * 1227., .04430,1220., 110.,   0.,   0.,   0.
     *,1227., .13290,1120.,-120.,   0.,   0.,   0.
     *,1227., .23790,2221., 120.,   0.,   0.,   0.
     *,1227., .30790,1221., 110.,   0.,   0.,   0.
     *,1227., .34290,1121.,-120.,   0.,   0.,   0.
     *,1227., .41190,1220., 120.,-120.,   0.,   0.
     *,1227., .48090,1220., 110., 110.,   0.,   0.
     *,1227., .54990,1120.,-120., 110.,   0.,   0.
     *,1227., .66070,1220., 220.,   0.,   0.,   0.
     *,1227., .72800,2130., 230.,   0.,   0.,   0.
     *,1227., .73870,1230., 230.,   0.,   0.,   0.
     *,1227., .76000,2230., 130.,   0.,   0.,   0.
     *,1227., .92000,1120.,-121.,   0.,   0.,   0.
     *,1227.,1.00000,1220., 111.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1036,1046)/
C    *----------DL0(1925)--------------------------
     * 1228., .18700,1220., 110.,   0.,   0.,   0.
     *,1228., .28000,1120.,-120.,   0.,   0.,   0.
     *,1228., .36400,2221., 120.,   0.,   0.,   0.
     *,1228., .37800,1221., 110.,   0.,   0.,   0.
     *,1228., .49000,1121.,-120.,   0.,   0.,   0.
     *,1228., .55700,1220., 111.,   0.,   0.,   0.
     *,1228., .69000,1120.,-121.,   0.,   0.,   0.
     *,1228., .69350,2230., 130.,   0.,   0.,   0.
     *,1228., .70000,1230., 230.,   0.,   0.,   0.
     *,1228., .80000,1122.,-120.,   0.,   0.,   0.
     *,1228.,1.00000,1222., 110.,   0.,   0.,   0./
C    *---------------------------------------------
C    *   LAMBDA RESONANCES + SIGMA0 RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1047,1059)/
C    *----------LAMBDA(1405)-----------------------
     * 1233., .33000,1230., 110.,   0.,   0.,   0.
     *,1233., .66000,2230., 120.,   0.,   0.,   0.
     *,1233.,1.00000,1130.,-120.,   0.,   0.,   0.
C    *----------LAMBDA(1520)-----------------------
     *,1234., .22500,1120.,-130.,   0.,   0.,   0.
     *,1234., .48000,1220.,-230.,   0.,   0.,   0.
     *,1234., .62000,1230., 110.,   0.,   0.,   0.
     *,1234., .76000,2230., 120.,   0.,   0.,   0.
     *,1234., .90000,1130.,-120.,   0.,   0.,   0.
     *,1234., .96000,2130., 120.,-120.,   0.,   0.
     *,1234., .99000,2130., 110., 110.,   0.,   0.
     *,1234., .99330,1130.,-120., 110.,   0.,   0.
     *,1234., .99660,2230., 120., 110.,   0.,   0.
     *,1234.,1.00000,1230., 120.,-120.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1060,1075)/
C    *----------LAMBDA(1645)-----------------------
     * 1235., .10000,1120.,-130.,   0.,   0.,   0.
     *,1235., .20000,1220.,-230.,   0.,   0.,   0.
     *,1235., .35000,1230., 110.,   0.,   0.,   0.
     *,1235., .50000,2230., 120.,   0.,   0.,   0.
     *,1235., .65000,1130.,-120.,   0.,   0.,   0.
     *,1235., .75000,2130., 120.,-120.,   0.,   0.
     *,1235., .80000,2130., 110., 110.,   0.,   0.
     *,1235., .84500,1130.,-120., 110.,   0.,   0.
     *,1235., .89000,2230., 120., 110.,   0.,   0.
     *,1235., .93500,1230., 120.,-120.,   0.,   0.
     *,1235.,1.00000,2130., 220.,   0.,   0.,   0.
C    *----------SIGMA0(1665)-----------------------
     *,1236., .10000,1120.,-130.,   0.,   0.,   0.
     *,1236., .20000,1220.,-230.,   0.,   0.,   0.
     *,1236., .40000,2230., 120.,   0.,   0.,   0.
     *,1236., .60000,1130.,-120.,   0.,   0.,   0.
     *,1236.,1.00000,2130., 110.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1076,1084)/
C    *----------SIGMA0(1776)-----------------------
     * 1237., .17500,1120.,-130.,   0.,   0.,   0.
     *,1237., .35000,1220.,-230.,   0.,   0.,   0.
     *,1237., .38750,2230., 120.,   0.,   0.,   0.
     *,1237., .42500,1130.,-120.,   0.,   0.,   0.
     *,1237., .57500,2130., 110.,   0.,   0.,   0.
     *,1237., .60000,2231., 120.,   0.,   0.,   0.
     *,1237., .62500,1131.,-120.,   0.,   0.,   0.
     *,1237., .75000,1234., 110.,   0.,   0.,   0.
     *,1237.,1.00000,1230., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1085,1094)/
C    *----------LAMBDA(1845)-----------------------
     * 1238., .17000,1120.,-130.,   0.,   0.,   0.
     *,1238., .34000,1220.,-230.,   0.,   0.,   0.
     *,1238., .44000,1230., 110.,   0.,   0.,   0.
     *,1238., .54000,2230., 120.,   0.,   0.,   0.
     *,1238., .64000,1130.,-120.,   0.,   0.,   0.
     *,1238., .70000,1231., 110.,   0.,   0.,   0.
     *,1238., .76000,2231., 120.,   0.,   0.,   0.
     *,1238., .82000,1131.,-120.,   0.,   0.,   0.
     *,1238., .91000,1120.,-131.,   0.,   0.,   0.
     *,1238.,1.00000,1220.,-231.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1095,1106)/
C    *----------SIGMA0(1930)-----------------------
     * 1239., .07500,1120.,-130.,   0.,   0.,   0.
     *,1239., .15000,1220.,-230.,   0.,   0.,   0.
     *,1239., .20000,1121.,-130.,   0.,   0.,   0.
     *,1239., .25000,1221.,-230.,   0.,   0.,   0.
     *,1239., .32500,1120.,-131.,   0.,   0.,   0.
     *,1239., .40000,1220.,-231.,   0.,   0.,   0.
     *,1239., .47500,2230., 120.,   0.,   0.,   0.
     *,1239., .55000,1130.,-120.,   0.,   0.,   0.
     *,1239., .70000,2130., 110.,   0.,   0.,   0.
     *,1239., .77500,2231., 120.,   0.,   0.,   0.
     *,1239., .85000,1131.,-120.,   0.,   0.,   0.
     *,1239.,1.00000,1234., 110.,   0.,   0.,   0./
C    *---------------------------------------------
C    *            SIGMA+ RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1107,1118)/
C    *----------SIGMA+(1665)-----------------------
     * 1132., .20000,1120.,-230.,   0.,   0.,   0.
     *,1132., .40000,1130., 110.,   0.,   0.,   0.
     *,1132., .60000,1230., 120.,   0.,   0.,   0.
     *,1132.,1.00000,2130., 120.,   0.,   0.,   0.
C    *----------SIGMA+(1776)-----------------------
     *,1133., .35000,1120.,-230.,   0.,   0.,   0.
     *,1133., .38750,1130., 110.,   0.,   0.,   0.
     *,1133., .42500,1230., 120.,   0.,   0.,   0.
     *,1133., .57500,2130., 120.,   0.,   0.,   0.
     *,1133., .60000,1131., 110.,   0.,   0.,   0.
     *,1133., .62500,1231., 120.,   0.,   0.,   0.
     *,1133., .75000,1234., 120.,   0.,   0.,   0.
     *,1133.,1.00000,1130., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1119,1128)/
C    *----------SIGMA+(1930)-----------------------
     * 1134., .15000,1120.,-230.,   0.,   0.,   0.
     *,1134., .22500,1111.,-130.,   0.,   0.,   0.
     *,1134., .25000,1121.,-230.,   0.,   0.,   0.
     *,1134., .40000,1120.,-231.,   0.,   0.,   0.
     *,1134., .47500,1130., 110.,   0.,   0.,   0.
     *,1134., .55000,1230., 120.,   0.,   0.,   0.
     *,1134., .70000,2130., 120.,   0.,   0.,   0.
     *,1134., .77500,1131., 110.,   0.,   0.,   0.
     *,1134., .85000,1231., 120.,   0.,   0.,   0.
     *,1134.,1.00000,1234., 120.,   0.,   0.,   0./
C    *---------------------------------------------
C    *            SIGMA- RESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1129,1140)/
C    *----------SIGMA-(1665)-----------------------
     * 2232., .20000,1220.,-130.,   0.,   0.,   0.
     *,2232., .40000,2230., 110.,   0.,   0.,   0.
     *,2232., .60000,1230.,-120.,   0.,   0.,   0.
     *,2232.,1.00000,2130.,-120.,   0.,   0.,   0.
C    *----------SIGMA-(1776)-----------------------
     *,2233., .35000,1220.,-130.,   0.,   0.,   0.
     *,2233., .38750,2230., 110.,   0.,   0.,   0.
     *,2233., .42500,1230.,-120.,   0.,   0.,   0.
     *,2233., .57500,2130.,-120.,   0.,   0.,   0.
     *,2233., .60000,2231., 110.,   0.,   0.,   0.
     *,2233., .62500,1231.,-120.,   0.,   0.,   0.
     *,2233., .75000,1234.,-120.,   0.,   0.,   0.
     *,2233.,1.00000,2230., 220.,   0.,   0.,   0./
      DATA ((DECTAB(I,J),I=1,7),J=1141,1150)/
C    *----------SIGMA-(1930)-----------------------
     * 2234., .15000,1220.,-130.,   0.,   0.,   0.
     *,2234., .17500,1221.,-130.,   0.,   0.,   0.
     *,2234., .25000,2221.,-230.,   0.,   0.,   0.
     *,2234., .40000,1220.,-131.,   0.,   0.,   0.
     *,2234., .47500,2230., 110.,   0.,   0.,   0.
     *,2234., .55000,1230.,-120.,   0.,   0.,   0.
     *,2234., .70000,2130.,-120.,   0.,   0.,   0.
     *,2234., .77500,2231., 110.,   0.,   0.,   0.
     *,2234., .85000,1231.,-120.,   0.,   0.,   0.
     *,2234.,1.00000,1234.,-120.,   0.,   0.,   0./
C    *---------------------------------------------
C    *      ADDITIONAL MESONRESONANCES
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1151,1159)/
C    *-----------F0(975)---------------------------
     *  332., .50000, 120.,-120.,   0.,   0.,   0.
     *, 332., .75000, 110., 110.,   0.,   0.,   0.
     *, 332., .87500, 130.,-130.,   0.,   0.,   0.
     *, 332.,1.00000, 230.,-230.,   0.,   0.,   0.
C    *-----------A0(980)---------------------------
     *, 112., .56000, 110., 220.,   0.,   0.,   0.
     *, 112., .78000, 130.,-130.,   0.,   0.,   0.
     *, 112.,1.00000, 230.,-230.,   0.,   0.,   0.
C    *-----------A+(980)---------------------------
     *, 122., .60000, 120., 220.,   0.,   0.,   0.
     *, 122.,1.00000, 130.,-230.,   0.,   0.,   0./
C    *---------------------------------------------
C    *      WEAK BARYON DECAYS
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1160,1169)/
C    *-----------LAMBDA(1116)----------------------
     * 2130.,0.64200,1120.,-120.,   0.,   0.,   0.
     *,2130.,1.00000,1220., 110.,   0.,   0.,   0.
C    *-----------SIGMA+(1180)----------------------
     *,1130.,0.51580,1120., 110.,   0.,   0.,   0.
     *,1130.,1.00000,1220., 120.,   0.,   0.,   0.
C    *-----------SIGMA-(1180)----------------------
     *,2230.,1.00000,1220.,-120.,   0.,   0.,   0.
C    *---------KASKADE-(1360)----------------------
     *,2330.,1.00000,2130.,-120.,   0.,   0.,   0.
C    *---------KASKADE0(1360)----------------------
     *,1330.,1.00000,2130., 110.,   0.,   0.,   0.
C    *---------OMEGA-(1680)------------------------
     *,3331.,0.68000,2130.,-130.,   0.,   0.,   0.
     *,3331.,0.82000,1330.,-120.,   0.,   0.,   0.
     *,3331.,1.00000,2330., 110.,   0.,   0.,   0./
C    *---------------------------------------------
C    *      WEAK MESON DECAYS
C    *---------------------------------------------
      DATA ((DECTAB(I,J),I=1,7),J=1170,1171)/
C    *-----------K0S(975)--------------------------
     *   20., .68610, 120.,-120.,   0.,   0.,   0.
     *,  20.,1.00000, 110., 110.,   0.,   0.,   0./
C    *---------------------------------------------
      DATA ALFA /0.00729735/, GF /1.16570E-5/, SIN2W /.215/
C  SINW = SQRT(SIN2W), COSW=SQRT(1.-SIN2W)
      DATA SINW /.463681/,COSW /.886002/
      SAVE
C-----------------------------------------------------------------------
      IF     ( IENTRO .EQ. 1 ) THEN
        CALL JCENTR(3,6,3,1)
      ELSEIF ( IENTRO .EQ. 2 ) THEN
        CALL JCENTD
      ENDIF
      CALL JCENTP

      CALL IDRESI

C  DETERMINE WMASS2,WGAM2
C  ----------------------
      AMW=SQRT(PI*ALFA/(.9304*1.41421356*GF))/SINW
      WMASS2=AMW
      CALL IDMASS(5,AMLEP5)
      CALL IDMASS(6,AMLEP6)
      IF ( AMLEP5+AMLEP6 .GT. AMW ) THEN
        NGAM=9
      ELSE
        NGAM=12
      ENDIF
      WGAM2=GF*AMW**3/(6.*PI*1.41421356)*NGAM

      IRD=0
      DO 1 I=1,MXLOOK
        LOOK(I)=0
 1    CONTINUE
      DO 2 I=1,MXDKY
        MODE(1,I)=0
        MODE(2,I)=0
        MODE(3,I)=0
        MODE(4,I)=0
        MODE(5,I)=0
        CBR(I)=0.
 2    CONTINUE
      NODCAY=.FALSE.
      NOETA=.FALSE.
      NOPI0=.FALSE.
      NONUNU=.FALSE.
      NOEVOL=.FALSE.
      NOHADR=.FALSE.
      IF ( LPRINT ) WRITE(IFCH,10)
10    FORMAT('1',30('*'),/,' *',28X,'*',/,
     *        ' *',5X,'ISAJET DECAY TABLE',5X,'*',/,
     *        ' *',28X,'*',/,' ',30('*'),/,/,
     *        6X,'PART',18X,'DECAY MODE',19X,'CUM BR',15X,'IDENT',17X,
     *        'DECAY IDENT',/)
      LOOP=0
      IOLD=0
      IF ( NODCAY ) RETURN

200   LOOP=LOOP+1
      IF ( LOOP .GT. MXDKY ) GOTO 9999
220   CONTINUE
      IMODE(1)=0
      IMODE(2)=0
      IMODE(3)=0
      IMODE(4)=0
      IMODE(5)=0
      LMODE(1)=IBLANK
      LMODE(2)=IBLANK
      LMODE(3)=IBLANK
      LMODE(4)=IBLANK
      LMODE(5)=IBLANK
      IRD=IRD+1
      IF ( IRD .GT. NDECTB ) RETURN
      IRES=NINT(DECTAB(1,IRD))
      BR=DECTAB(2,IRD)
      IMODE(1)=NINT(DECTAB(2+1,IRD))
      IMODE(2)=NINT(DECTAB(2+2,IRD))
      IMODE(3)=NINT(DECTAB(2+3,IRD))
      IMODE(4)=NINT(DECTAB(2+4,IRD))
      IMODE(5)=NINT(DECTAB(2+5,IRD))
      IF ( NOPI0  .AND.  IRES .EQ. 110 ) GOTO 220
      IF ( NOETA  .AND.  IRES .EQ. 220 ) GOTO 220
      IF ( IRES .EQ. IOLD ) GOTO 230
      IF ( IRES .LT. 0  .OR.  IRES .GT. MXLOOK ) THEN
        CALL UTSTOP('JDECIN: IRES OUT OF RANGE               ')
      ENDIF
      LOOK(IRES)=LOOP
230   IOLD=IRES
      CBR(LOOP)=BR
      MODE(1,LOOP)=IMODE(1)
      MODE(2,LOOP)=IMODE(2)
      MODE(3,LOOP)=IMODE(3)
      MODE(4,LOOP)=IMODE(4)
      MODE(5,LOOP)=IMODE(5)
      IF ( LPRINT ) THEN
        IF ( IMODE(1) .NE. 0 ) LMODE(1)=IDLABL(IMODE(1))
        IF ( IMODE(2) .NE. 0 ) LMODE(2)=IDLABL(IMODE(2))
        IF ( IMODE(3) .NE. 0 ) LMODE(3)=IDLABL(IMODE(3))
        IF ( IMODE(4) .NE. 0 ) LMODE(4)=IDLABL(IMODE(4))
        IF ( IMODE(5) .NE. 0 ) LMODE(5)=IDLABL(IMODE(5))
        LRES=IDLABL(IRES)
        WRITE(IFCH,20) LRES,(LMODE(K),K=1,5), BR,IRES,(IMODE(K),K=1,5)
20      FORMAT(6X,A5,6X,5(A5,2X),3X,F8.5,15X,I5,4X,5(I5,2X))
      ENDIF
      GOTO 200

9999  WRITE(IFCH,*)'LOOP=', LOOP
      CALL UTSTOP('JDECIN: LOOP > MXDKY                    ')
      RETURN
      END
C=======================================================================

      SUBROUTINE JESTPR(IC1,IC2,AM,IER)

C-----------------------------------------------------------------------
C  PROCESSES STRINGS
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CKOL/    KOL
      COMMON /CLEAD/   COOAV3,COOAV4,LEAD
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CPZSTR/  ESTRL,PZSTRL,ISEA,ISTRL
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL    STRO(NSI,NSIX+1)
      INTEGER IC1(2),IC2(2),JC(NFLAV,2),JC1(NFLAV,2),JC2(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IER=0

C  PMAX
C  ----
      CALL IDDECO(IC1,JC)
      NQ=0
      DO 7 NF=1,NFLAV
        NQ=NQ+JC(NF,1)-JC(NF,2)
 7    CONTINUE
      IF ( ABS(NQ) .GE. 2 ) THEN
        AM1MIN=PROM
      ELSE
        AM1MIN=PIOM
      ENDIF
      CALL IDDECO(IC2,JC)
      NQ=0
      DO 8 NF=1,NFLAV
        NQ=NQ+JC(NF,1)-JC(NF,2)
 8    CONTINUE
      IF ( ABS(NQ) .GE. 2 ) THEN
        AM2MIN=PROM
      ELSE
        AM2MIN=PIOM
      ENDIF
      IF ( AM .LE. AM1MIN+AM2MIN ) THEN
        PMAX=AM*0.5
      ELSE
        PMAX=UTPCM(AM,AM1MIN,AM2MIN)
      ENDIF

C  HASTPR
C  ------
      IPROJ=1
      ITARG=1
      KPROJA(2,1)=1
      KTARGA(2,1)=1
      KPROJA(3,1)=1
      KTARGA(3,1)=1
      KOL=1
      COORD(1,1)=0.
      COORD(2,1)=0.
      COORD(3,1)=0.
      COORD(4,1)=0.
      CALL IDDECO(IC1,JC1)
      NPA1=0
      DO 2 N=1,NFLAV
        NPA1=NPA1+JC1(N,1)+JC1(N,2)
 2    CONTINUE
      CALL IDDECO(IC2,JC2)
      NPA2=0
      DO 3 N=1,NFLAV
        NPA2=NPA2+JC2(N,1)+JC2(N,2)
 3    CONTINUE
      IF ( NPA2 .GT. NPA1 ) THEN
        II=2
      ELSE
        II=1
      ENDIF
      IF ( NPA1 .GT. 1  .OR.  NPA2 .GT. 1 ) THEN
        LEAD=1
      ELSE
        LEAD=0
      ENDIF

      NSTR0=NSTR
17    NSTR=NSTR0

      DO 15 N=1,NSI
        STRO(N,1)=0.
        STRO(N,2)=0.
        STRO(N,3)=0.
15    CONTINUE
      STRO(3,II)=AM*0.5
      STRO(4,II)=AM*0.5
      STRO(5,II)=IC1(1)
      STRO(6,II)=IC1(2)
      STRO(3,3-II)=-AM*0.5
      STRO(4,3-II)=AM*0.5
      STRO(5,3-II)=IC2(1)
      STRO(6,3-II)=IC2(2)

      PZSTRL=STRO(3,1)
      ESTRL=STRO(4,1)
      PZSTRL=PZSTRL+STRO(3,2)
      ESTRL=ESTRL+STRO(4,2)
      ISTRL=0

      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,110)('-',L=1,79),IPAGE,('-',L=1,79)
110     FORMAT(1X,79A1,/,1X,I5,'.PAGE            '
     *           ,'STRING GENERATION',/,1X,79A1,/)
        WRITE(IFCH,105)(STRO(I,1),I=1,4),(NINT(STRO(I,1)),I=5,6)
105     FORMAT(' STR: ',4F13.5,2I8)
        WRITE(IFCH,104)(STRO(I,2),I=1,4),(NINT(STRO(I,2)),I=5,6)
104     FORMAT('      ',4F13.5,2I8,/)
      ENDIF

      ISEA=1
      ISPLT=0
14    CALL HASTPR(STRO,ISPLT)
      IF ( ISPLT .EQ. -1 ) GOTO 9001
      IF ( ISPLT .EQ. -3 ) GOTO 9001
      IF ( ISPLT .EQ. -4 ) THEN
        CALL UTSTOP('JESTPR: ISPLT=-4                        ')
      ENDIF
      IF ( ISPLT .EQ. -5 ) GOTO 17
      IF ( ISPLT .GT.  0 ) GOTO 14

      PMXEVT=PMAX
      EGYEVT=AM

      RETURN

9001  IER=1
      RETURN
      END
C=======================================================================

      SUBROUTINE JETGEN(IER)

C-----------------------------------------------------------------------
C  GENERATES STRINGS
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      PARAMETER (MAMX=56)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CLEP/    ICINPU,IDSCAT
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CNTEVM/  NTEVM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      REAL    TARGA(NSI,NHA,MAMX)
      INTEGER ICSTRI(4,18),IC1(2),IC2(2),IC4(2)

      DATA ((ICSTRI(I,J),I=1,4),J=1,8)/
     *100000,000000,110000,000000,
     *010000,000000,200000,000000,
     *100000,000000,210000,100000,
     *010000,000000,210000,010000,
     *001000,000000,210000,001000,
     *000000,100000,310000,000000,
     *000000,010000,220000,000000,
     *000000,001000,211000,000000/
      DATA ((ICSTRI(I,J),I=1,4),J=11,18)/
     *100000,000000,020000,000000,
     *010000,000000,110000,000000,
     *100000,000000,120000,100000,
     *010000,000000,120000,010000,
     *001000,000000,120000,001000,
     *000000,100000,220000,000000,
     *000000,010000,130000,000000,
     *000000,001000,121000,000000/
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 2 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)'STRING GENERATION. SR JETGEN.'
      ENDIF

      IER=0

      NEVT=1
      NSTR=0
      NPTL=0

      IF ( ICINPU .GE. 1 ) THEN
3       CALL LEPEXP(XBJ,QSQ)
        RNU=QSQ/(2.* PROM * XBJ)
        WSQ=PROM**2 + 2.* PROM * RNU - QSQ
        AMST=SQRT(WSQ)
        ELEPTO=ELEPTI-RNU
        COSANG=1.-QSQ/(2.*ELEPTI*ELEPTO)
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,*)'ELEPTI,ELEPTO,RNU: ',ELEPTI,ELEPTO,RNU
          WRITE(IFCH,*)'XBJ,QSQ,COSANG=1-QSQ/2/ELEPTI/ELEPTO: '
     *                 ,XBJ,QSQ,COSANG
        ENDIF
        IF ( RNU .GT. ELEPTI ) THEN
          IF ( ISH .GE. 93 ) WRITE(IFCH,*)'*****  Q0 TOO LARGE: ',RNU
          GOTO 3
        ENDIF
        IF     ( COSANG.GT.1.  .AND.  COSANG.LT.1.001 ) THEN
          COSANG=1.
        ELSEIF ( COSANG.LT.-1.  .OR.  COSANG.GT.1.    ) THEN
          COSANO=COSANG
          IF ( COSANG .GE. 0. ) COSANG= 1.
          IF ( COSANG .LT. 0. ) COSANG=-1.
          IF(ISH.GE.90)THEN
            CALL UTMSG('JETGEN')
            WRITE(IFCH,*)'*****  COSANG OUT OF RANGE'
            WRITE(IFCH,*)'ELEPTI,ELEPTO,RNU: ',ELEPTI,ELEPTO,RNU
            WRITE(IFCH,*)'XBJ,QSQ,COSANG=1-QSQ/2/ELEPTI/ELEPTO: '
     *                   ,XBJ,QSQ,COSANO
            WRITE(IFCH,*)'COSANG_NEW: ',COSANG
            CALL UTMSGF
          ENDIF
        ENDIF
        ANGMUE=ACOS(COSANG)
        CALL LEPTAR(XBJ,QSQ,MATARG,LATARG,IDSCAT)
      ENDIF

      IF ( MATARG .GT. 0 ) THEN
        CALL NUCOGE
        NPTL=0
        CALL NUCINI('STR',TARGA,LATARG,MATARG,-1)
        CALL NUCSTR(IER)
        IF ( IER .EQ. 1 ) GOTO 99999
      ENDIF

      IF     ( ICINPU .EQ. 0 ) THEN
        R=RANGEN()
        PS=0.
        DO 1 K=1,99
          PS=PS+PROB(K)
          IF ( R .LE. PS ) GOTO 2
 1      CONTINUE
        CALL UTSTOP('JETGEN: NO K FOUND                      ')
 2      CONTINUE
        IC1(1)=ICFOR(K,1)
        IC1(2)=ICFOR(K,2)
        IC2(1)=ICBAC(K,1)
        IC2(2)=ICBAC(K,2)
        AMST=ENGY
      ELSEIF ( ICINPU .GT. 0 ) THEN
        CALL LEPSTR(IDSCAT,XBJ,QSQ,IDS)
        IC1(1)=ICSTRI(1,IDS)
        IC1(2)=ICSTRI(2,IDS)
        IC2(1)=ICSTRI(3,IDS)
        IC2(2)=ICSTRI(4,IDS)
      ENDIF
      CALL JESTPR(IC1,IC2,AMST,IER)

99999 ISH=ISH0
      RETURN
      END
C=======================================================================

      SUBROUTINE JFRADE(IER)

C-----------------------------------------------------------------------
C  PERFORMS STRING FRAGMENTATION/DECAY AND FIN. STATE INTERACTIONS
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXDKY=2000)
      PARAMETER (MXINDX=1000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXRE=100)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NPRBMS=20)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CKOL/    KOL
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNFR/    NRFRA
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPRBMS / PRBMS(NPRBMS)
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      COMMON /DIDIB/   NDIDIB
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      DOUBLE PRECISION DIFF,ZFI
      REAL             PSUM(5)
      INTEGER          IC4(2)
      SAVE
C-----------------------------------------------------------------------
C  INITIALIZATION
C  --------------
CDH   IF ( ISH .EQ. 13 .OR. ISH .EQ. 14 ) CALL UTTIMA('*** JFRADE *** ')
      IER=0
      IRET=0
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 16 ) ISH=MOD(ISHSUB,100)
      IF ( NEVT .NE. 1  .OR.  IFRADE .EQ. 0 ) GOTO 1000
      NPTLPT=ABS(MAPROJ)+ABS(MATARG)
      IF ( ICHOIC .EQ. 1  .OR.  ICHOIC .EQ. 4 ) THEN
        YCMMAX=LOG(EGYEVT*2.5)
        ETAPRO=YCMMAX*.6667
        ETATAR=-ETAPRO
      ELSE
        ETAPRO=(YPJTL-YHAHA)*.6667
        ETATAR=-YHAHA*.6667
      ENDIF
      DETAP=ETAPRO
      DETAT=ETATAR
      TPRO=COSH(DETAP)
      ZPRO=SINH(DETAP)
      TTAR=COSH(DETAT)
      ZTAR=SINH(DETAT)

C  WRITE
C  -----
      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,110)('-',L=1,79),IPAGE,('-',L=1,79)
110     FORMAT(1X,79A1,/,1X,I5,'.PAGE            '
     *           ,'STRINGS BEFORE RESCALING',/,1X,79A1,/)
        DO 9 J=1,NSTR
          WRITE(IFCH,109)J,(ICSTR(K,J)/100,K=1,4)
     *          ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J),PSTR(5,J)
     *          ,IRLSTR(J)
109       FORMAT(' /CSTR/',I4,3X,4I5,3(E11.3),I4)
 9      CONTINUE
      ENDIF

C  STRING RESCALING
C  ----------------
      IF ( ICHOIC.LE.2 .AND. NSTR.GT.1 .AND. IRESCL.EQ.1 ) THEN
        PSUM(1)=0.
        PSUM(2)=0.
        IF ( ICHOIC .EQ. 1 ) THEN
          PSUM(3)=0.
          PSUM(4)=EGYEVT
        ELSE
          PSUM(3)=(NPJEVT-NTGEVT)*PNLLX
          PSUM(4)=NPJEVT*SQRT(AMPROJ**2+PNLLX**2)
     *           +NTGEVT*SQRT(AMTARG**2+PNLLX**2)
        ENDIF
        PSUM(5)=SQRT(PSUM(4)**2-PSUM(3)**2)
        CALL HRESCL(1,NSTR,PSUM,IFAIL)
        IF ( IFAIL .NE. 0 ) GOTO 1001
      ENDIF

C  WRITE
C  -----
      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,111)('-',L=1,79),IPAGE,('-',L=1,79)
111     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                       ,'STRINGS AFTER RESCALING',/,1X,79A1,/)
        DO 10 J=1,NSTR
          WRITE(IFCH,109)J,(ICSTR(K,J)/100,K=1,4)
     *        ,SQRT(PSTR(1,J)**2+PSTR(2,J)**2),PSTR(3,J),PSTR(5,J)
     *        ,IRLSTR(J)
10      CONTINUE
      ENDIF

C  FRAGMENTATION
C  -------------
      IF ( ISHSUB/100 .EQ. 3 ) ISH=MOD(ISHSUB,100)
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('               ')
      DO 3 J=1,NSTR
        CALL UTPAGE
        IF ( ISH .GE. 91 ) WRITE(IFCH,102)('-',L=1,79),IPAGE,J
     *       ,(ICSTR(K,J),K=1,4),SQRT(PSTR(1,J)**2+PSTR(2,J)**2
     *            +PSTR(3,J)**2),PSTR(4,J),PSTR(5,J),('-',L=1,79)
102     FORMAT(/,1X,79A1,
     *         /,1X,I5,'.PAGE  STR:',I3,4I7,3(E10.2),/,1X,79A1,/)
        CALL JAMFRA(J,NEWEVT)
        IF ( NEWEVT .EQ. 1 ) GOTO 10011
 3    CONTINUE
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('FRAGMENTATION  ')

C  PRINT /CPTL/
C  ------------
      IF ( ISHSUB/100 .EQ. 9 ) ISH=MOD(ISHSUB,100)
      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,104)('-',L=1,79),IPAGE,('-',L=1,79)
104     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                       ,'PTLS AFTER FRAGMENTATION',/,1X,79A1,/)
        DO 18 N=1,NPTL
          IF ( N.GT.NPTLPT .OR. ISTPTL(N).NE.0 )
     *      WRITE(IFCH,116)IORPTL(N),JORPTL(N),N,IFRPTL(1,N),IFRPTL(2,N)
     *                 ,IDPTL(N),PPTL(3,N),PPTL(4,N),PPTL(5,N),ISTPTL(N)
18      CONTINUE
      ENDIF

C  INITIAL DECAY
C  -------------
      IF ( ISHSUB/100 .EQ. 10 ) ISH=MOD(ISHSUB,100)
      CALL UTPAGE
CDH   IF ( ISH.EQ.13 ) CALL UTTIMA('               ')
      IF ( RADIAC .GT. 0. ) THEN
        TTAUS=TAUMIN
        IACN=1
      ELSE
        TTAUS=AINFIN
        IACN=0
      ENDIF
      TTP=TTAUS*TPRO
      TTT=TTAUS*TTAR
      ZZP=TTAUS*ZPRO
      ZZT=TTAUS*ZTAR
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,119)('-',L=1,79),IPAGE,SNGL(TTAUS),('-',L=1,79)
119     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                       ,'DECAY BEFORE TAU =',E10.3,/,1X,79A1,/)
      ENDIF
      NP1=1
21    NP2=NPTL
      DO 5 I=NP1,NP2
        IF ( TTAUS .LE. 0.D0   ) THEN
          DIFF = TIVPTL(2,I) - TTAUS
        ELSE
          ZFI=XORPTL(3,I)+(TIVPTL(2,I)-XORPTL(4,I))*PPTL(3,I)/PPTL(4,I)
          CALL UTTAUT(SNGL(ZFI),TZFI)
          DIFF = TIVPTL(2,I) - TZFI
        ENDIF
        IF ( DIFF .LE. 0. .OR. IACN.EQ.0 ) THEN
          IF ( ISTPTL(I) .EQ. 0 ) THEN
            CALL JDECA(I,IRET)
            IF ( IRET .EQ. 1 ) GOTO 1001
          ENDIF
        ENDIF
5     CONTINUE
      NP1=NP2+1
      IF ( NP1 .LE. NPTL ) GOTO 21
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('DECAY INI      ')

C  INTERACTION AND DECAY
C  ---------------------
      IF ( IACN .EQ. 0 ) GOTO 5000
      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,118)('-',L=1,79),IPAGE,('-',L=1,79)
118     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                       ,'INTERACTIONS AND DECAY',/,1X,79A1,/)
      ENDIF
      DTAUS=1./(NUMTAU-1.)*(TAUMAX-TAUMIN)
      DO 23 NT=1,NUMTAU
        NTC=NT
        TTAUS=TAUMIN+(NT-1.)*DTAUS
        TTP=TTAUS*TPRO
        TTT=TTAUS*TTAR
        ZZP=TTAUS*ZPRO
        ZZT=TTAUS*ZTAR
        IF     ( IOJINT .EQ. 1 ) THEN
          CALL JINTA1
        ELSEIF ( IOJINT .EQ. 2 ) THEN
          CALL JINTA2
        ENDIF
        TTAUS=TTAUS+DTAUS
        TTP=TTAUS*TPRO
        TTT=TTAUS*TTAR
        ZZP=TTAUS*ZPRO
        ZZT=TTAUS*ZTAR
        NP1=1
36      NP2=NPTL
        DO 37 IP=NP1,NP2
          IF ( TTAUS .LE. 0.D0   ) THEN
            DIFF = TIVPTL(2,IP) - TTAUS
          ELSE
            ZFI=XORPTL(3,IP)+(TIVPTL(2,IP)-XORPTL(4,IP))
     *                                            *PPTL(3,IP)/PPTL(4,IP)
            IF     ( ZFI .LE. ZZT ) THEN
              DIFF = TIVPTL(2,IP) - TTT - (ZFI-ZZT)*ZZT/TTT
            ELSEIF ( ZFI .GE. ZZP ) THEN
              DIFF = TIVPTL(2,IP) - TTP - (ZFI-ZZP)*ZZP/TTP
            ELSE
              IF ( TTAUS .LT. AINFIN ) THEN
                IF ( TIVPTL(2,IP) .GE. 0. ) THEN
                  DIFF = TIVPTL(2,IP)
                  DIFF = DIFF**2  - (TTAUS**2+ZFI**2)
                ELSE
                  DIFF = TIVPTL(2,IP) - SQRT(TTAUS**2+ZFI**2)
                ENDIF
              ELSE
                DIFF = TIVPTL(2,IP) - TTAUS
                IF ( ISH .GE. 90 ) THEN
                  CALL UTMSG('JFRADE')
                  WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
                  WRITE(IFCH,*)'TTAUS=',TTAUS,'ZFI=',ZFI
                  CALL UTMSGF
                ENDIF
              ENDIF
            ENDIF
          ENDIF

          IF ( DIFF .LE. 0. ) THEN
            IF ( ISTPTL(IP) .EQ. 0 ) THEN
              CALL JDECA(IP,IRET)
              IF ( IRET .EQ. 1 ) GOTO 1001
            ENDIF
          ENDIF
37      CONTINUE
        NP1=NP2+1
        IF ( NP1 .LE. NPTL ) THEN
          DO 4 IP=NP1,NPTL
            IAAPTL(IP)=1
 4        CONTINUE
          GOTO 36
        ENDIF
CDH     IF ( ISH .EQ. 13 ) CALL UTTIMA('DECAY          ')
23    CONTINUE

C  FINAL DECAY
C  -----------
      CALL UTPAGE
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,120)('-',L=1,79),IPAGE,('-',L=1,79)
120     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                        ,'FINAL DECAY',/,1X,79A1,/)
      ENDIF
      NBEF=NPTL
CDH   N0BEF=NPTL0
      NAFT=NPTL
CDH   N0AFT=NPTL0
      IF ( NCLEAN .GT. 0 ) THEN
        CALL UTCLEA(NPTL0)
        NAFT=NPTL
CDH     N0AFT=NPTL0
      ENDIF
      ISHNPT=ISH
      IF ( ISHSUB/100 .EQ. 19 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .EQ. 22 ) THEN
        WRITE(IFMT,131)NBEF,NAFT
131     FORMAT(1X,'BEF FIN DE: NBEF=',I8,4X,'NAFT=',I8)
      ENDIF
      ISH=ISHNPT
      NP1=1
41    NP2=NPTL
      DO 42 IP=NP1,NP2
        IF ( ISTPTL(IP) .EQ. 0 ) THEN
          CALL JDECA(IP,IRET)
          IF ( IRET .EQ. 1 ) GOTO 1001
        ENDIF
42    CONTINUE
      NP1=NP2+1
      IF ( NP1 .LE. NPTL ) THEN
        DO 6 IP=NP1,NPTL
          IAAPTL(IP)=1
 6      CONTINUE
        GOTO 41
      ENDIF
      ISHNPT=ISH
      IF ( ISHSUB/100 .EQ. 19 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .EQ. 22 ) THEN
        WRITE(IFMT,132)NPTL
132     FORMAT(1X,'AFT FIN DE: NPTL=',I8)
      ENDIF
      ISH=ISHNPT
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('DECAY FIN      ')

C  PRINT /CPTL/
C  ------------
      CALL UTPAGE
      IF ( ISHSUB/100 .EQ. 11 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,117)('-',L=1,79),IPAGE,('-',L=1,79)
117     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *         ,'PTLS AFTER PERFORMING INTERACTIONS',/,1X,79A1,/
     *         1X,' IOR',' JOR',4X,'   N',4X,'IFR1IFR2',10X,'ID',3X,
     *         5X,'PTR',7X,'PZ',4X,'MASS ','IST',/)
        DO 34 N=1,NPTL
          IF ( N .GT. NPTLPT  .OR.  ISTPTL(N) .NE. 0 )
     *      WRITE(IFCH,116)IORPTL(N),JORPTL(N),N,IFRPTL(1,N),IFRPTL(2,N)
     *     ,IDPTL(N),SQRT(PPTL(1,N)**2+PPTL(2,N)**2),PPTL(3,N),PPTL(5,N)
     *     ,ISTPTL(N),IAAPTL(N)
116       FORMAT(1X,I7,I7,4X,I7,4X,I4,I4,I12,3(E10.2),1X,I2,I2)
34      CONTINUE
      ENDIF

C  TRAFO -> LAB CM
C  ---------------
5000  CONTINUE
      IF ( LABSYS .EQ. 1 ) THEN
        DO 7 I=1,NPTL
          AMT=SQRT(PPTL(5,I)**2+PPTL(1,I)**2+PPTL(2,I)**2)
          PZ=PPTL(3,I)
          E=PPTL(4,I)
          YI=SIGN( LOG((E+ABS(PZ))/AMT), PZ )
          Y=YI+YHAHA
          PPTL(3,I)=AMT*SINH(Y)
          PPTL(4,I)=AMT*COSH(Y)
 7      CONTINUE
      ENDIF

C  FINISH
C  ------

      CALL UTPAGE
      IF     ( ISH .GE. 91 ) THEN
        WRITE(IFCH,113)('-',L=1,79),IPAGE,('-',L=1,79)
113     FORMAT(/,1X,79A1,/,1X,I5,'.PAGE            '
     *                 ,'PARTICLE PRODUCTION FINISHED',/,1X,79A1,/)
CDH   ELSEIF ( ISH .EQ. 13 ) THEN
CDH     CALL UTTIMA('               ')
      ENDIF
      GOTO 1000

10011 CONTINUE
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('JFRADE')
        WRITE(IFCH,*)'*****  FRAGMENTATION NOT POSSIBLE'
        WRITE(IFCH,112)J,(ICSTR(K,J)/100,K=1,4)
     *                  ,(PSTR(L,J),L=1,5)
112     FORMAT( ' STR:',3X,I4,3X,3X,4I5,5F7.2,F7.2)
        CALL UTMSGF
      ENDIF
      IER=1
      GOTO 1000

1001  CONTINUE
      IER=1
      GOTO 1000

1000  CONTINUE
      IF ( IER .EQ. 1 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JFRADE')
          WRITE(IFCH,*)'*****  REDO EVENT',NREVT+1
          CALL UTMSGF
        ENDIF
      ENDIF
CDH   IF ( ISH .EQ. 14 ) CALL UTTIMA('    JFRADE F   ')
      ISH=ISH0

      RETURN
      END
C=======================================================================

      SUBROUTINE JINTA1

C----------------------------------------------------------------------
C  SECONDARY INTERACTIONS
C----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXDKY=2000)
      PARAMETER (MXINDX=1000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXIFR=MXPTL)
      PARAMETER (MXRE=100)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NFLAV=6)
      PARAMETER (NPRBMS=20)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CIFRIJ/  IFRIJ(MXIFR)
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CKOL/    KOL
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNFR/    NRFRA
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPRBMS/  PRBMS(NPRBMS)
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      COMMON /C4PTL/   OPTL(MXPTL),TPTL(MXPTL),UPTL(MXPTL)
     *                ,XPTL(MXPTL),YPTL(MXPTL),ZPTL(MXPTL)
      COMMON /DIDIB/   NDIDIB
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      DOUBLE PRECISION DD,DERR,TI1,TI2,TT,VV,VVP,VVT,XO3,XO4,ZZ,ZZA
      INTEGER          IC4(2),JC(NFLAV,2),JCDU(NFLAV,2),JCI(NFLAV,2)
      LOGICAL          IACPTL(MXPTL+10)
      DATA DERR/1.D-2/
      SAVE
C----------------------------------------------------------------------
C  INITIALIZATION FOR NTC=1
C  ------------------------
      IF ( NTC .EQ. 1 ) THEN
        IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
          IF ( MAPROJ .EQ. 1 ) THEN
            RNUP=0.
          ELSE
            RNUP=1.19*MAPROJ**(0.3333333)-1.61*MAPROJ**(-0.3333333)
          ENDIF
          IF ( MATARG .EQ. 1 ) THEN
            RNUT=0.
          ELSE
            RNUT=1.19*MATARG**(0.3333333)-1.61*MATARG**(-0.3333333)
          ENDIF
          RNUP=RNUP+OVERLP
          RNUT=RNUT+OVERLP
          VELP=PPTL(3,1)/PPTL(4,1)
          VELT=PPTL(3,IABS(MAPROJ)+1)/PPTL(4,IABS(MAPROJ)+1)
          BX=COS(PHIEVT)*BIMEVT
          BY=SIN(PHIEVT)*BIMEVT
        ENDIF
        NPTL0=NPTL
        DO 51 I=1,NPTL0
          IAAPTL(I)=1
          CALL IDQUAC(I,NQI,NDUMMY,NDUMMY,JCI)
          IF ( NQI .NE. 0 ) THEN
            RADPTL(I)=RADIAC
            AMIPTL(I)=PROM+AMSIAC
          ELSE
            RADPTL(I)=RADIAS
            AMIPTL(I)=PIOM+AMSIAC
          ENDIF
51      CONTINUE
      ENDIF

C  INITIALIZATION FOR EACH NTC
C  ---------------------------
      NT=NTC
      CALL UTPAGE
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('               ')
      NBEF=NPTL
CDH   N0BEF=NPTL0
      NAFT=NPTL
CDH   N0AFT=NPTL0
      IF ( NCLEAN .GT. 0 ) THEN
        IF ( MOD(NT-1,NCLEAN) .EQ. 0 ) THEN
          CALL UTCLEA(NPTL0)
          NAFT=NPTL
CDH       N0AFT=NPTL0
        ENDIF
      ENDIF
      TAUS=TTAUS
      NPTLPT=ABS(MAPROJ)+ABS(MATARG)
      IF ( NPTL .GT. NPTL0 ) THEN
        DO 52 I=NPTL0+1,NPTL
          CALL IDQUAC(I,NQI,NDUMMY,NDUMMY,JCI)
          IF ( NQI .NE. 0 ) THEN
            RADPTL(I)=RADIAC
            AMIPTL(I)=PROM+AMSIAC
          ELSE
            RADPTL(I)=RADIAS
            AMIPTL(I)=PIOM+AMSIAC
          ENDIF
52      CONTINUE
      ENDIF
      DO 44 I=1,NPTL
        IF ( ICLPTL(I).EQ.0 .AND. PPTL(5,I).GT.AMIPTL(I) ) THEN
          IACPTL(I)=.TRUE.
        ELSE
          IACPTL(I)=.FALSE.
        ENDIF
        IF ( IORPTL(I).EQ.-1 ) IACPTL(I)=.TRUE.
C  CALL UTTAIN IS REPLACED HERE TO AVOID OVERHEAD
C       CALL UTTAIN(I,X,Y,Z,T,N,0)
        XO4=XORPTL(4,I)
C*      IF     ( IOPT .EQ. 0 ) THEN
          TI1=TIVPTL(1,I)
C*      ELSEIF ( IOPT .EQ. 1 ) THEN
C*        TI1=XO4
C*      ENDIF
        TI2=TIVPTL(2,I)
        IF ( TI1 .GT. TI2 ) GOTO 1009
        PPT4I = 1./PPTL(4,I)
        VV=PPTL(3,I)*PPT4I
        XO3=XORPTL(3,I)
        ZZ=XO3+(TI2-XO4)*VV
        IF ( TTAUS .LE. 0.D0 ) THEN
          TZ=TTAUS
        ELSE
          IF     ( ZZ .LE. ZZT ) THEN
            TZ=TTT+(ZZ-ZZT)*ZZT/TTT
          ELSEIF ( ZZ .GE. ZZP ) THEN
            TZ=TTP+(ZZ-ZZP)*ZZP/TTP
          ELSE
            IF ( TTAUS .GE. AINFIN ) THEN
              TZ=TTAUS
              IF ( ISH .GE. 90 ) THEN
                CALL UTMSG('JINTA1')
                WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
                WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
                CALL UTMSGF
              ENDIF
            ELSE
C*DH          TZ=SQRT(TTAUS**2+ZZ**2)
              IF ( TI2 .LT. 0.D0 ) GOTO 1002
              IF ( TTAUS**2+ZZ**2 .GE. TI2**2 ) GOTO 1002
              GOTO 1006
            ENDIF
          ENDIF
        ENDIF
        IF ( TZ .GE. TI2 ) GOTO 1002
 1006   ZZ=XO3+(TI1-XO4)*VV
        IF ( TTAUS .GT. 0.D0 ) THEN
          IF     ( ZZ .LE. ZZT ) THEN
            TZ=TTT+(ZZ-ZZT)*ZZT/TTT
          ELSEIF ( ZZ .GE. ZZP ) THEN
            TZ=TTP+(ZZ-ZZP)*ZZP/TTP
          ELSE
            IF ( TTAUS .GE. AINFIN ) THEN
              TZ=TTAUS
              IF ( ISH .GE. 90 ) THEN
                CALL UTMSG('JINTA1')
                WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
                WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
                CALL UTMSGF
              ENDIF
            ELSE
C*DH          TZ=SQRT(TTAUS**2+ZZ**2)
              IF ( TI1 .LT. 0.D0 ) GOTO 1007
              IF ( TTAUS**2+ZZ**2 .LE. TI1**2 ) GOTO 1001
              GOTO 1007
            ENDIF
          ENDIF
        ENDIF
        IF ( TZ .LE. TI1 ) GOTO 1001
 1007   IF ( TTAUS .LE. 0.D0 ) THEN
          TT=TTAUS
          ZZ=XO3+(TT-XO4)*VV
          IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1031
        ELSE
          ZZA=XO3-XO4*VV
          VVT=ZZT/TTT
          TT=(TTT+(ZZA-ZZT)*VVT)/(1.D0-VV*VVT)
          ZZ=XO3+(TT-XO4)*VV
          IF ( ZZ .LE. ZZT ) THEN
            IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1032
            GOTO 1000
          ENDIF
          VVP=ZZP/TTP
          TT=(TTP+(ZZA-ZZP)*VVP)/(1.D0-VV*VVP)
          ZZ=XO3+(TT-XO4)*VV
          IF ( ZZ .GE. ZZP ) THEN
            IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1033
            GOTO 1000
          ENDIF
          DD=1.D0-VV**2
          IF ( DD .EQ. 0.D0 ) THEN
            TT=-VV*(TTAUS**2+ZZA**2)*0.5D0/ZZA
          ELSE
            TT=(ZZA*VV+SQRT(ZZA**2+TTAUS**2*DD))/DD
          ENDIF
          ZZ=XO3+(TT-XO4)*VV
          IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1034
          IF ( TT .LT. 0.D0 ) GOTO 1035
          IF ( ZZ.LE.ZZT .OR. ZZ.GE.ZZP ) GOTO 1004
          IF ( ABS(TTAUS**2-(TT+ZZ)*(TT-ZZ))
     *                         .GT. DERR*TTAUS**2 ) GOTO 1005
        ENDIF
 1000   N=0
 1011   T=TT
        Z=ZZ
        X=XORPTL(1,I)+(T-XO4)*PPTL(1,I)*PPT4I
        Y=XORPTL(2,I)+(T-XO4)*PPTL(2,I)*PPT4I
        GOTO 53
 1001   N=1
        GOTO 53
 1002   N=2
        GOTO 53
 1031   N=31
        GOTO 1003

 1032   N=32
        GOTO 1003
 1033   N=33
        GOTO 1003
 1034   N=34
 1003   IF ( ABS(TT-TI1) .LE. DERR*ABS(TT) ) GOTO 1000
        IF ( ABS(TT-TI2) .LE. DERR*ABS(TT) ) GOTO 1000
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTA1')
          WRITE(IFCH,*)'*****  TI1 < TT < TI2   NOT FULFILLED - ',N
          WRITE(IFCH,*)SNGL(TI1),SNGL(TT),SNGL(TI2)
          CALL UTMSGF
        ENDIF
        GOTO 1011
 1035   CONTINUE
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTA1')
          WRITE(IFCH,*)'*****  TT < 0     ( ',TT,' )'
          WRITE(IFCH,*)'VV,DD:',VV,DD
          WRITE(IFCH,*)'ZZA,TTAUS:',ZZA,TTAUS
          CALL UTMSGF
        ENDIF
        GOTO 1011
 1004   N=4
        IF ( ABS(ZZ-ZZT) .LE. DERR*ABS(ZZ) ) GOTO 1000
        IF ( ABS(ZZ-ZZP) .LE. DERR*ABS(ZZ) ) GOTO 1000
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTA1')
          WRITE(IFCH,*)'*****  ZZT < ZZ < ZZP   NOT FULFILLED'
          WRITE(IFCH,*)SNGL(ZZT),SNGL(ZZ),SNGL(ZZP)
          CALL UTMSGF
        ENDIF
        GOTO 1011
 1005   N=5
        IF ( ABS(TTAUS**2-(TT+ZZ)*(TT-ZZ)) .LE. DERR ) GOTO 1000
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTA1')
          WRITE(IFCH,*)'*****  TTAUS**2 .NE. (TT+ZZ)*(TT-ZZ)'
          WRITE(IFCH,*)SNGL(TTAUS**2),SNGL((TT+ZZ)*(TT-ZZ))
          CALL UTMSGF
        ENDIF
        GOTO 1011
 1009   N=9
 53     CONTINUE

        IF ( N.EQ.2 ) THEN
          IF ( NCLEAN .GT. 0 ) THEN
            IF ( MOD(NT,NCLEAN).EQ.0 ) THEN
              ISTPTL(I)=2
            ENDIF
          ENDIF
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( N.EQ.9 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( N.EQ.1 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( IAAPTL(I) .EQ. 0 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( ISTPTL(I) .EQ. 2 ) IACPTL(I)=.TRUE.
        IF ( ABS(IDPTL(I)) .LT. 100
     *       .AND.  ABS(IDPTL(I)) .NE. 20 ) IACPTL(I)=.TRUE.
        IF ( IDPTL(I) .EQ. 441  .AND.  JPSIFI .EQ. 0 ) IACPTL(I)=.TRUE.
 54     XPTL(I)=X
        YPTL(I)=Y
        ZPTL(I)=Z
        TPTL(I)=T
        IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
          IF ( (X-BX*0.5)**2+(Y-BY*0.5)**2+(Z-VELP*TAUS)**2
     *                                    .LT. RNUP**2  .AND.
     *         (X+BX*0.5)**2+(Y+BY*0.5)**2+(Z-VELT*TAUS)**2
     *                                .LT. RNUT**2  )  IACPTL(I)=.TRUE.
        ENDIF
44    CONTINUE
      NPTL0=NPTL
      I=0

C  I LOOP --> 24
C  -------------
9999  I=I+1
      IF ( IACPTL(I) ) GOTO 24
      J0=NPTLPT+1
      IF ( I .GT. NPTLPT ) J0=I+1
      IF ( I .GT. NPTL0 ) J0=1

C  J LOOP --> 250
C  -------------
      J = J0-1
25    J = J+1
      IF ( J .GT. NPTL ) GOTO 24
        IF ( IACPTL(J) ) THEN
          IF ( IACPTL(J+1) ) THEN
            IF ( IACPTL(J+2) ) THEN
              IF ( IACPTL(J+3) ) THEN
                IF ( IACPTL(J+4) ) THEN
                  IF ( IACPTL(J+5) ) THEN
                    IF ( IACPTL(J+6) ) THEN
                      IF ( IACPTL(J+7) ) THEN
                        IF ( IACPTL(J+8) ) THEN
                          IF ( IACPTL(J+9) ) THEN
                            IF ( IACPTL(J+10) ) THEN
                              J=J+10
                              GOTO 25
                            ELSE
                              J=J+10
                              IF ( J .GT. NPTL ) GOTO 24
                              GOTO 26
                            ENDIF
                          ELSE
                            J=J+9
                            IF ( J .GT. NPTL ) GOTO 24
                            GOTO 26
                          ENDIF
                        ELSE
                          J=J+8
                          IF ( J .GT. NPTL ) GOTO 24
                          GOTO 26
                        ENDIF
                      ELSE
                        J=J+7
                        IF ( J .GT. NPTL ) GOTO 24
                        GOTO 26
                      ENDIF
                    ELSE
                      J=J+6
                      IF ( J .GT. NPTL ) GOTO 24
                      GOTO 26
                    ENDIF
                  ELSE
                    J=J+5
                    IF ( J .GT. NPTL ) GOTO 24
                    GOTO 26
                  ENDIF
                ELSE
                  J=J+4
                  IF ( J .GT. NPTL ) GOTO 24
                  GOTO 26
                ENDIF
              ELSE
                J=J+3
                IF ( J .GT. NPTL ) GOTO 24
                GOTO 26
              ENDIF
            ELSE
              J=J+2
              IF ( J .GT. NPTL ) GOTO 24
              GOTO 26
            ENDIF
          ELSE
            J=J+1
            IF ( J .GT. NPTL ) GOTO 24
            GOTO 26
          ENDIF
        ENDIF
 26     CONTINUE
        RADSQR=(RADPTL(I)+RADPTL(J))**2
        IF ( (ZPTL(I)-ZPTL(J))**2 .GT. RADSQR ) GOTO 25
        IF ( (YPTL(I)-YPTL(J))**2 .GT. RADSQR ) GOTO 25
        IF ( (XPTL(I)-XPTL(J))**2 .GT. RADSQR ) GOTO 25
        IF ( I .EQ. J ) GOTO 25
        IF ( IORPTL(I).GT.0 .AND. IORPTL(J).EQ.IORPTL(I) ) GOTO 25
        IF ( IORPTL(I) .EQ. J ) GOTO 25
        IF ( IORPTL(J) .EQ. I ) GOTO 25
        PDE=(PPTL(3,I)+PPTL(3,J))/(PPTL(4,I)+PPTL(4,J))
        GAM2I=1.-PDE**2
        IF ( GAM2I .EQ. 0. ) GOTO 25
        IF ( (ZPTL(I)-ZPTL(J)-(TPTL(I)-TPTL(J))*PDE)**2
     *                                 .GT. RADSQR*GAM2I ) GOTO 25
        IF ( (XPTL(I)-XPTL(J))**2+(YPTL(I)-YPTL(J))**2+
     *    1./GAM2I*(ZPTL(I)-ZPTL(J)-(TPTL(I)-TPTL(J))*PDE)**2
     *                                       .GT. RADSQR ) GOTO 25
        CALL JINTCC(I,J,IRET)
        IF ( IRET .EQ. 1 ) GOTO 25

        IACTN=0
        NPTL00=NPTL
        NSTR00=NSTR
        XAVER(1)=(XPTL(I)+XPTL(J))*0.5
        XAVER(2)=(YPTL(I)+YPTL(J))*0.5
        XAVER(3)=(ZPTL(I)+ZPTL(J))*0.5
        XAVER(4)=(TPTL(I)+TPTL(J))*0.5

        CALL JINTFS(I,J,NQIFUS,JC,AMIM,IRET)
        IF ( IRET .EQ. 1 ) GOTO 25

        CALL JINTCE(I,J,AMIM,IACTN,IRET)
        IF ( IRET .EQ. 25 ) GOTO25

        TIVPTL(2,I)=TPTL(I)
        TIVPTL(2,J)=TPTL(J)
        ISTPTL(I)=1
        ISTPTL(J)=1
        IACPTL(I)=.TRUE.
        IACPTL(J)=.TRUE.

        CALL JINTCH(I,J,KMAX)
        DO 30 K=1,KMAX
          N=IFRIJ(K)
          ISTPTL(N)=2
          IACPTL(N)=.TRUE.
30      CONTINUE
        CALL JINTPA(I,J,KMAX)
        DO 31 K=1,KMAX
          N=IFRIJ(K)
          IACPTL(N)=.TRUE.
31      CONTINUE

        IF ( IACTN .EQ. 1 ) THEN
          DO 32 N=NPTL00+1,NPTL
            IAAPTL(N)=1
            CALL IDQUAC(N,NQI,NDUMMY,NDUMMY,JCDU)
            IF ( NQI .EQ. 0 ) THEN
              RADPTL(N)=RADIAS
              AMIPTL(N)=PIOM+AMSIAC
            ELSE
              RADPTL(N)=RADIAC
              AMIPTL(N)=PROM+AMSIAC
            ENDIF
            IACPTL(N)=.FALSE.
            IF ( PPTL(5,N) .GT. AMIPTL(N) ) IACPTL(N)=.TRUE.
            CALL UTTAIN(N,X,Y,Z,T,K,0)
            IF ( K.EQ.1 .OR. K.EQ.2 .OR. K.EQ.9 ) IACPTL(N)=.TRUE.
            IF ( ABS(IDPTL(N)) .LT. 100
     *           .AND. ABS(IDPTL(N)) .NE. 20 ) IACPTL(N)=.TRUE.
            XPTL(N)=X
            YPTL(N)=Y
            ZPTL(N)=Z
            TPTL(N)=T
            IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
              IF ( (X-BX*.5)**2+(Y-BY*.5)**2+(Z-VELP*TAUS)**2
     *                                 .LT. RNUP**2  .AND.
     *             (X+BX*.5)**2+(Y+BY*.5)**2+(Z-VELT*TAUS)**2
     *                                 .LT. RNUT**2 )  IACPTL(N)=.TRUE.
            ENDIF
32        CONTINUE
          GOTO 24
        ENDIF

        CALL JINTEL(I,J,AMIM,IACTN)
        IF ( IACTN .EQ. 2 ) THEN
          DO 33 N=NPTL00+1,NPTL
            IF     ( N .EQ. NPTL00+1 ) THEN
              IJ=I
            ELSEIF ( N .EQ. NPTL00+2 ) THEN
              IJ=J
            ENDIF
            RADPTL(N)=RADPTL(IJ)
            AMIPTL(N)=AMIPTL(IJ)
            XPTL(N)=XAVER(1)
            YPTL(N)=XAVER(2)
            ZPTL(N)=XAVER(3)
            TPTL(N)=XAVER(4)
            IACPTL(N)=.FALSE.
            IAAPTL(N)=1
33        CONTINUE
          GOTO 24
        ENDIF

        CALL JINTFU(I,J,JC,IACTN)
        DO 34 N=NPTL00+1,NPTL
          IF ( NQIFUS .EQ. 0 ) THEN
            RADPTL(N)=RADIAS
            AMIPTL(N)=PIOM+AMSIAC
          ELSE
            RADPTL(N)=RADIAC
            AMIPTL(N)=PROM+AMSIAC
          ENDIF
          XPTL(N)=XAVER(1)
          YPTL(N)=XAVER(2)
          ZPTL(N)=XAVER(3)
          TPTL(N)=XAVER(4)
          IACPTL(N)=.FALSE.
          IAAPTL(N)=1
34      CONTINUE
*       GOTO 24

*250   GOTO 25

24    CONTINUE
      IF ( I .LT. NPTL-1 ) GOTO 9999

      ISHNPT=ISH
      IF ( ISHSUB/100 .EQ. 19 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .EQ. 22 ) THEN
        WRITE(IFMT,131)NT,NBEF,NAFT,NPTL
131     FORMAT(1X,'NT=',I5,4X,'NBEF=',I8,4X,'NAFT=',I8,4X,'NPTL=',I8)
      ENDIF
      ISH=ISHNPT
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('INTERACTIONS   ')
      RETURN
      END
C=======================================================================

      SUBROUTINE JINTA2

C----------------------------------------------------------------------
C  SECONDARY INTERACTIONS: PERCOLATION--CLUSTER MODEL
C----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXDKY=2000)
      PARAMETER (MXINDX=1000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXIFR=MXPTL)
      PARAMETER (MXRE=100)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NFLAV=6)
      PARAMETER (NPRBMS=20)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CIFRIJ/  IFRIJ(MXIFR)
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CKOL/    KOL
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNFR/    NRFRA
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPRBMS/  PRBMS(NPRBMS)
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      COMMON /C3PTL/   DESPTL(MXPTL),DEZPTL(MXPTL)
      COMMON /C4PTL/   OPTL(MXPTL),TPTL(MXPTL),UPTL(MXPTL)
     *                ,XPTL(MXPTL),YPTL(MXPTL),ZPTL(MXPTL)
      COMMON /DIDIB/   NDIDIB
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      INTEGER          IC4(2),JC(NFLAV,2),JCDU(NFLAV,2),JCI(NFLAV,2)
      LOGICAL          IACPTL(MXPTL+10)
      SAVE
C----------------------------------------------------------------------
C  INITIALIZATION FOR NTC=1
C  ------------------------
      IF ( NTC .EQ. 1 ) THEN
        IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
          IF ( MAPROJ .EQ. 1 ) THEN
            RNUP=0.
          ELSE
            RNUP=1.19*MAPROJ**(.3333333)-1.61*MAPROJ**(-.3333333)
          ENDIF
          RNUP=RNUP+OVERLP
          IF ( MATARG .EQ. 1 ) THEN
            RNUT=0.
          ELSE
            RNUT=1.19*MATARG**(.3333333)-1.61*MATARG**(-.3333333)
          ENDIF
          RNUT=RNUT+OVERLP
          VELP=PPTL(3,1)/PPTL(4,1)
          VELT=PPTL(3,IABS(MAPROJ)+1)/PPTL(4,IABS(MAPROJ)+1)
          BX=COS(PHIEVT)*BIMEVT
          BY=SIN(PHIEVT)*BIMEVT
        ENDIF
        NPTL0=NPTL
        DO 51 I=1,NPTL0
          IAAPTL(I)=1
          CALL IDQUAC(I,NQI,NDUMMY,NDUMMY,JCI)
          IF ( NQI .EQ. 0 ) THEN
            RADPTL(I)=RADIAS
            DESPTL(I)=RADIAS
            DEZPTL(I)=0.
            AMIPTL(I)=PIOM+AMSIAC
          ELSE
            RADPTL(I)=RADIAC
            DESPTL(I)=RADIAC
            DEZPTL(I)=0.
            AMIPTL(I)=PROM+AMSIAC
          ENDIF
51      CONTINUE
      ENDIF

C  INITIALIZATION FOR EACH NTC
C  ---------------------------
      NT=NTC
      CALL UTPAGE
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('               ')
      NBEF=NPTL
CDH   N0BEF=NPTL0
      NAFT=NPTL
CDH   N0AFT=NPTL0
      IF ( NCLEAN .GT. 0 ) THEN
        IF ( MOD(NT-1,NCLEAN) .EQ. 0 ) THEN
          CALL UTCLEA(NPTL0)
          NAFT=NPTL
CDH       N0AFT=NPTL0
        ENDIF
      ENDIF
      TAUS=TTAUS
      NPTLPT=ABS(MAPROJ)+ABS(MATARG)
      IF ( NPTL .GT. NPTL0 ) THEN
        DO 52 I=NPTL0+1,NPTL
          CALL IDQUAC(I,NQI,NDUMMY,NDUMMY,JCI)
          IF ( NQI .EQ. 0 ) THEN
            RADPTL(I)=RADIAS
            DESPTL(I)=RADIAS
            DEZPTL(I)=0.
            AMIPTL(I)=PIOM+AMSIAC
          ELSE
            RADPTL(I)=RADIAC
            DESPTL(I)=RADIAC
            DEZPTL(I)=0.
            AMIPTL(I)=PROM+AMSIAC
          ENDIF
52      CONTINUE
      ENDIF
      DO 44 I=1,NPTL
        IF ( ICLPTL(I).EQ.0 .AND. PPTL(5,I).GT.AMIPTL(I) ) THEN
          IACPTL(I)=.TRUE.
        ELSE
          IACPTL(I)=.FALSE.
        ENDIF
        IF ( IORPTL(I).EQ.-1 ) IACPTL(I)=.TRUE.
        CALL UTTAIN(I,X,Y,Z,T,N,0)
        IF ( N.EQ.2 ) THEN
          IF ( NCLEAN .GT. 0 ) THEN
            IF ( MOD(NT,NCLEAN).EQ.0 ) THEN
              ISTPTL(I)=2
            ENDIF
          ENDIF
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( N.EQ.9 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( N.EQ.1 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( IAAPTL(I) .EQ. 0 ) THEN
          IACPTL(I)=.TRUE.
          GOTO 54
        ENDIF
        IF ( ISTPTL(I) .EQ. 2 ) IACPTL(I)=.TRUE.
        IF ( ABS(IDPTL(I)) .LT. 100
     *       .AND.  ABS(IDPTL(I)) .NE. 20 ) IACPTL(I)=.TRUE.
        IF ( IDPTL(I) .EQ. 441  .AND.  JPSIFI .EQ. 0 ) IACPTL(I)=.TRUE.
 54     XPTL(I)=X
        YPTL(I)=Y
        ZPTL(I)=Z
        TPTL(I)=T
        CALL UTTAUS(Z,SZ)
        OPTL(I)=SZ+DESPTL(I)
        UPTL(I)=SZ-DESPTL(I)
        IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
          IF ( (X-BX*.5)**2+(Y-BY*.5)**2+(Z-VELP*TAUS)**2 .LT. RNUP**2
     *         .AND. (X+BX*.5)**2+(Y+BY*.5)**2+(Z-VELT*TAUS)**2
     *                                  .LT. RNUT**2 ) IACPTL(I)=.TRUE.
        ENDIF
44    CONTINUE
      NPTL0=NPTL
      I=0

C  I LOOP --> 24
C  -------------
9999  I=I+1
      IF ( IACPTL(I) ) GOTO 24
      J0=NPTLPT+1
      IF ( I .GT. NPTLPT ) J0=I+1
      IF ( I .GT. NPTL0 ) J0=1
      UPTLI = UPTL(I)
      OPTLI = OPTL(I)

C  J LOOP --> 250
C  -------------
      J = J0-1
25    J = J+1
      IF ( J .GT. NPTL ) GOTO 24
        IF ( IACPTL(J) ) THEN
          IF ( IACPTL(J+1) ) THEN
            IF ( IACPTL(J+2) ) THEN
              IF ( IACPTL(J+3) ) THEN
                IF ( IACPTL(J+4) ) THEN
                  IF ( IACPTL(J+5) ) THEN
                    IF ( IACPTL(J+6) ) THEN
                      IF ( IACPTL(J+7) ) THEN
                        IF ( IACPTL(J+8) ) THEN
                          IF ( IACPTL(J+9) ) THEN
                            IF ( IACPTL(J+10) ) THEN
                              J=J+10
                              GOTO 25
                            ELSE
                              J=J+10
                              IF ( J .GT. NPTL ) GOTO 24
                              GOTO 26
                            ENDIF
                          ELSE
                            J=J+9
                            IF ( J .GT. NPTL ) GOTO 24
                            GOTO 26
                          ENDIF
                        ELSE
                          J=J+8
                          IF ( J .GT. NPTL ) GOTO 24
                          GOTO 26
                        ENDIF
                      ELSE
                        J=J+7
                        IF ( J .GT. NPTL ) GOTO 24
                        GOTO 26
                      ENDIF
                    ELSE
                      J=J+6
                      IF ( J .GT. NPTL ) GOTO 24
                      GOTO 26
                    ENDIF
                  ELSE
                    J=J+5
                    IF ( J .GT. NPTL ) GOTO 24
                    GOTO 26
                  ENDIF
                ELSE
                  J=J+4
                  IF ( J .GT. NPTL ) GOTO 24
                  GOTO 26
                ENDIF
              ELSE
                J=J+3
                IF ( J .GT. NPTL ) GOTO 24
                GOTO 26
              ENDIF
            ELSE
              J=J+2
              IF ( J .GT. NPTL ) GOTO 24
              GOTO 26
            ENDIF
          ELSE
            J=J+1
            IF ( J .GT. NPTL ) GOTO 24
            GOTO 26
          ENDIF
        ENDIF
 26     IF ( OPTLI .LT. UPTL(J) ) GOTO 25
        IF ( UPTLI .GE. OPTL(J) ) GOTO 25
        RADSQR=(RADPTL(I)+RADPTL(J))**2
        IF ( (YPTL(I)-YPTL(J))**2 .GT. RADSQR ) GOTO 25
        IF ( (XPTL(I)-XPTL(J))**2 .GT. RADSQR ) GOTO 25
        IF ( I .EQ. J ) GOTO 25
        IF ( IORPTL(I).GT.0 .AND. IORPTL(J).EQ.IORPTL(I) ) GOTO 25
        IF ( IORPTL(I) .EQ. J ) GOTO 25
        IF ( IORPTL(J) .EQ. I ) GOTO 25
        IF ( (XPTL(I)-XPTL(J))**2+(YPTL(I)-YPTL(J))**2 .GT. RADSQR )
     *                                                       GOTO 25
        CALL JINTCC(I,J,IRET)
        IF ( IRET .EQ. 1 ) GOTO 25

        IACTN=0
        NPTL00=NPTL
        NSTR00=NSTR
        XAVER(1)=(XPTL(I)+XPTL(J))*0.5
        XAVER(2)=(YPTL(I)+YPTL(J))*0.5
        XAVER(3)=(ZPTL(I)+ZPTL(J))*0.5
        XAVER(4)=(TPTL(I)+TPTL(J))*0.5

        CALL JINTFS(I,J,NQIFUS,JC,AMIM,IRET)
        IF ( IRET .EQ. 1 ) GOTO 25

        CALL JINTCE(I,J,AMIM,IACTN,IRET)
        IF ( IRET .EQ. 25 ) GOTO 25

        TIVPTL(2,I)=TPTL(I)
        TIVPTL(2,J)=TPTL(J)
        ISTPTL(I)=1
        ISTPTL(J)=1
        IACPTL(I)=.TRUE.
        IACPTL(J)=.TRUE.

        CALL JINTCH(I,J,KMAX)
        DO 30 K=1,KMAX
          N=IFRIJ(K)
          ISTPTL(N)=2
          IACPTL(N)=.TRUE.
30      CONTINUE
        CALL JINTPA(I,J,KMAX)
        DO 31 K=1,KMAX
          N=IFRIJ(K)
          IACPTL(N)=.TRUE.
31      CONTINUE

        IF ( IACTN .EQ. 1 ) THEN
          DO 32 N=NPTL00+1,NPTL
            IAAPTL(N)=1
            CALL IDQUAC(N,NQI,NDUMMY,NDUMMY,JCDU)
            IF ( NQI .EQ. 0 ) THEN
              RADPTL(N)=RADIAS
              DESPTL(N)=RADIAS
              DEZPTL(N)=0.
              AMIPTL(N)=PIOM+AMSIAC
            ELSE
              RADPTL(N)=RADIAC
              DESPTL(N)=RADIAC
              DEZPTL(N)=0.
              AMIPTL(N)=PROM+AMSIAC
            ENDIF
            IACPTL(N)=.FALSE.
            IF ( PPTL(5,N) .GT. AMIPTL(N) ) IACPTL(N)=.TRUE.
            CALL UTTAIN(N,X,Y,Z,T,K,0)
            IF ( K.EQ.1 .OR. K.EQ.2 .OR. K.EQ.9 ) IACPTL(N)=.TRUE.
            IF ( ABS(IDPTL(N)).LT.100 .AND. ABS(IDPTL(N)).NE.20 )
     *                                              IACPTL(N)=.TRUE.
            XPTL(N)=X
            YPTL(N)=Y
            ZPTL(N)=Z
            TPTL(N)=T
            CALL UTTAUS(Z,SZ)
            OPTL(N)=SZ+DESPTL(N)
            UPTL(N)=SZ-DESPTL(N)
            IF ( OVERLP.GE.0. .AND. MAPROJ.NE.0 .AND. MATARG.NE.0 ) THEN
              IF ( (X-BX*.5)**2+(Y-BY*.5)**2+(Z-VELP*TAUS)**2.LT.RNUP**2
     *           .AND. (X+BX*.5)**2+(Y+BY*.5)**2+(Z-VELT*TAUS)**2
     *                                    .LT.RNUT**2 ) IACPTL(N)=.TRUE.
            ENDIF
32        CONTINUE
          GOTO 24
        ENDIF

        CALL JINTEL(I,J,AMIM,IACTN)
        IF ( IACTN .EQ. 2 ) THEN
          DO 33 N=NPTL00+1,NPTL
            IF     ( N .EQ. NPTL00+1 ) THEN
              IJ=I
            ELSEIF ( N .EQ. NPTL00+2 ) THEN
              IJ=J
            ENDIF
            RADPTL(N)=RADPTL(IJ)
            DESPTL(N)=DESPTL(IJ)
            DEZPTL(N)=DEZPTL(IJ)
            AMIPTL(N)=AMIPTL(IJ)
            XPTL(N)=XAVER(1)
            YPTL(N)=XAVER(2)
            ZPTL(N)=XAVER(3)
            TPTL(N)=XAVER(4)
            Z=ZPTL(N)
            CALL UTTAUS(Z,SZ)
            OPTL(N)=SZ+DESPTL(N)
            UPTL(N)=SZ-DESPTL(N)
            IACPTL(N)=.FALSE.
            IAAPTL(N)=1
33        CONTINUE
          GOTO 24
        ENDIF

        CALL JINTFU(I,J,JC,IACTN)
        DO 34 N=NPTL00+1,NPTL
          IF ( NQIFUS .EQ. 0 ) THEN
            RADPTL(N)=RADIAS
            DESPTL(N)=RADIAS
            DEZPTL(N)=0.
            AMIPTL(N)=PIOM+AMSIAC
          ELSE
            RADPTL(N)=RADIAC
            DESPTL(N)=RADIAC
            DEZPTL(N)=0.
            AMIPTL(N)=PROM+AMSIAC
          ENDIF
          XPTL(N)=XAVER(1)
          YPTL(N)=XAVER(2)
          ZPTL(N)=XAVER(3)
          TPTL(N)=XAVER(4)
          Z=ZPTL(N)
          CALL UTTAUS(Z,SZ)
          OPTL(N)=SZ+DESPTL(N)
          UPTL(N)=SZ-DESPTL(N)
          IACPTL(N)=.FALSE.
          IAAPTL(N)=1
34      CONTINUE
        IF ( ISH .EQ. 23  .AND.  PPTL(5,NPTL) .GE. AMPRIF ) THEN
          N=NPTL
          CALL JINTFP(I,J,N,
     *           XPTL(I),YPTL(I),RADPTL(I),OPTL(I),UPTL(I),
     *           XPTL(J),YPTL(J),RADPTL(J),OPTL(J),UPTL(J),
     *           XPTL(N),YPTL(N),RADPTL(N),OPTL(N),UPTL(N))
        ENDIF
*       GOTO 24

*250  GOTO 25

24    CONTINUE
      IF ( I .LT. NPTL-1 ) GOTO 9999

      IF ( ISH .EQ. 24 ) THEN
        IF ( TAUS.EQ.1.  .OR.  TAUS.EQ.2.  .OR.
     *       TAUS.EQ.4.  .OR.  TAUS.EQ.8.      ) THEN
          DO 45 I=1,NPTL
            IF ( IACPTL(I)  .AND.  I .NE. NPTL ) GOTO 45
            CALL JINTCL(I,XPTL(I),YPTL(I),RADPTL(I)
     *                                       ,OPTL(I),UPTL(I),IACPTL(I))
45        CONTINUE
        ENDIF
      ENDIF

      ISHNPT=ISH
      IF ( ISHSUB/100 .EQ. 19 ) ISH=MOD(ISHSUB,100)
      IF ( ISH .EQ. 22 ) THEN
        WRITE(IFMT,131)NT,NBEF,NAFT,NPTL
131     FORMAT(1X,'NT=',I5,4X,'NBEF=',I8,4X,'NAFT=',I8,4X,'NPTL=',I8)
      ENDIF
      ISH=ISHNPT
CDH   IF ( ISH .EQ. 13 ) CALL UTTIMA('INTERACTIONS   ')
      RETURN
      END
C=======================================================================

      SUBROUTINE JINTCC(I,J,IRET)

C---------------------------------------------------------------------
C  IRET=1 IF I = CHILD OF J
C  IRET=1 IF J = CHILD OF I
C  IRET=0 ELSE
C---------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (MXIFR=MXPTL)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)

      INTEGER IFR(MXIFR)
      SAVE
C---------------------------------------------------------------------
      IRET=0
      DO 57 KK=1,2
        IF ( KK .EQ. 1 ) THEN
          N=I
          M=J
        ELSE
          N=J
          M=I
        ENDIF
        K1=0
        K2=0
        LOOP=0
55      LOOP=LOOP+1
        IF ( IFRPTL(1,N) .GT. 0 ) THEN
          DO 56 K=IFRPTL(1,N),IFRPTL(2,N)
            IF ( K .EQ. M ) GOTO 1001
            K2=K2+1
            IF ( K2 .GT. MXIFR ) THEN
              CALL UTSTOP('JINTCC: K2 > MXIFR                      ')
            ENDIF
            IFR(K2)=K
56        CONTINUE
        ENDIF
        K1=K1+1
        IF ( K1 .LE. K2 ) THEN
          N=IFR(K1)
          GOTO 55
        ENDIF
57    CONTINUE
      GOTO 1000

1001  IRET=1

1000  RETURN
      END
C=======================================================================

      SUBROUTINE JINTCE(I,J,AMIM,IACTN,IRET)

C----------------------------------------------------------------------
C  COLOUR EXCHANGE INTERACTION OF PTLS I,J
C  INPUT:
C  I,J: PTL NUMBERS; AMF: MASS, AMIM: MIN MASS, OF FUSED OBJ
C  IACTN=1: CE DONE
C  IRET=0: CE DONE  OR  CRITERIA FOR CE NOT FULFILLED
C  IRET=25: SKIP
C----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXDKY=2000)
      PARAMETER (MXINDX=1000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXMA=11)
      PARAMETER (MXMX=6)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXRE=100)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NPRBMS=20)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CKOL/    KOL
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNFR/    NRFRA
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPRBMS/  PRBMS(NPRBMS)
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CREMA/   REMA(MXRE,MXMA),REWI(MXRE,MXMA)
     *                ,ICRE1(MXRE,MXMA),ICRE2(MXRE,MXMA)
     *                ,IDMX(MXMA,MXMX),INDX(MXINDX)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      COMMON /DIDIB/   NDIDIB
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      DOUBLE PRECISION ARM(5),ARP(5),ARQ,BOO(5),ROT(3)
      REAL             PROJ(NSI,NHA),PSUM(5),TARG(NSI,NHA)
      INTEGER          IC(2),IC4(2)
      SAVE
C----------------------------------------------------------------------
C  INITIALIZATION
C  --------------
      IRET=0
      DELRAP=1.5
      DELAMF=1.0
      ISH00=ISH
      PNLLX0=PNLLX
      TAUS=TTAUS
      AMF=PPTL(5,NPTL+1)

C  CHECK WHETHER CE CRITERIA FULFILLED
C  -----------------------------------
      VEI=PPTL(3,I)/PPTL(4,I)
      VEJ=PPTL(3,J)/PPTL(4,J)
      IF ( ABS(VEI) .LT. 1.  .AND.  ABS(VEJ) .LT. 1. ) THEN
        RAI=0.5*LOG((1.+VEI)/(1.-VEI))
        RAJ=0.5*LOG((1.+VEJ)/(1.-VEJ))
      ELSE
        RAI=0.
        RAJ=0.
      ENDIF
      IF ( .NOT. (AMF.GT.AMIM+DELAMF .AND. ABS(RAI-RAJ).GT.DELRAP
     *            .AND. ABS(IDPTL(I)).LT.10000
     *            .AND. ABS(IDPTL(J)).LT.10000) ) GOTO 1000
      IACTN=1

C  PRINT
C  -----
      IF ( ISH .GE. 91 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)' '
        WRITE(IFCH,101)NTC,TAUS
101     FORMAT(1X,'CO.EX. INTERACTION --- NT=',I3,' --- TAUS=',F6.2)
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)' '
        WRITE(IFCH,115)I,IDPTL(I)
     *                  ,(PPTL(K,I),K=3,5),(XORPTL(K,I),K=3,4)
115     FORMAT(1X,'/CPTL/',I8,I10
     *            ,1X,2(E10.2),E10.2,1X,2(E10.2),2X,I4)
        WRITE(IFCH,115)J,IDPTL(J)
     *                  ,(PPTL(K,J),K=3,5),(XORPTL(K,J),K=3,4)
        IF ( ISH .GE. 92 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'AMF,AMIM: ',AMF,AMIM
          WRITE(IFCH,*)'RAI,RAJ:  ',RAI,RAJ
          WRITE(IFCH,*)' '
        ENDIF
      ENDIF

C  INITIALIZATION
C  --------------
      NPTL00=NPTL
      NSTR00=NSTR

      IPROJ=1
      ITARG=1
      DO 11 L=1,NHA
        KPROJA(L,1)=1
        KTARGA(L,1)=1
11    CONTINUE
      KOL=1
      COORD(1,1)=XAVER(1)
      COORD(2,1)=XAVER(2)
      COORD(3,1)=XAVER(3)
      COORD(4,1)=XAVER(4)

      IDP=IDPTL(I)
      IDM=IDPTL(J)
      DO 12 L=1,5
        ARP(L)=PPTL(L,I)
        ARM(L)=PPTL(L,J)
        BOO(L)=PPTL(L,NPTL+1)
12    CONTINUE
      SROOT = ABS(PPTL(5,NPTL+1))
      S=PPTL(5,NPTL+1)**2

C  BOOSTS INTO I-J CM
C  ------------------
      CALL UTLOB2(1,BOO(1),BOO(2),BOO(3),BOO(4),BOO(5)
     *             ,ARP(1),ARP(2),ARP(3),ARP(4))
      CALL UTLOB2(1,BOO(1),BOO(2),BOO(3),BOO(4),BOO(5)
     *             ,ARM(1),ARM(2),ARM(3),ARM(4))
      IF ( ARP(3) .LT. 0.D0 ) THEN
        IDQ=IDM
        IDM=IDP
        IDP=IDQ
        DO 14 L=1,5
          ARQ=ARM(L)
          ARM(L)=ARP(L)
          ARP(L)=ARQ
14      CONTINUE
      ENDIF
      PNLLX=ARP(3)
      ROT(1)=(ARP(1)-ARM(1))*0.5D0
      ROT(2)=(ARP(2)-ARM(2))*0.5D0
      ROT(3)=(ARP(3)-ARM(3))*0.5D0
      CALL UTROT2(1,ROT(1),ROT(2),ROT(3)
     *             ,ARP(1),ARP(2),ARP(3))
      CALL UTROT2(1,ROT(1),ROT(2),ROT(3)
     *             ,ARM(1),ARM(2),ARM(3))

C  CHECKS
C  ------
      IF ( ARP(3) .LT. 0.D0 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTCE')
          WRITE(IFCH,*)'*****  Z-COMP OF +JET < 0.'
          WRITE(IFCH,*)(SNGL(ARP(L)),L=1,4)
          WRITE(IFCH,*)(SNGL(ARM(L)),L=1,4)
          CALL UTMSGF
        ENDIF
        GOTO 10025
      ENDIF
      IF ( ABS(SNGL(ARP(4)+ARM(4))-SROOT) .GT.
     *                                 2.E-2*SROOT ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTCE')
          WRITE(IFCH,*)'*****  ARP(4)+ARM(4)-SQRT(S) NONZERO'
          WRITE(IFCH,*)'VALUE:   ',SNGL(ARP(4)+ARM(4))-SROOT
          WRITE(IFCH,*)'SQRT(S): ',SROOT
          WRITE(IFCH,*)(SNGL(ARP(I)),I=1,4)
          WRITE(IFCH,*)(SNGL(ARM(I)),I=1,4)
          CALL UTMSGF
        ENDIF
      ENDIF

C  FILL PROJ, TARG
C  ---------------
      CALL IDTR4(IDP,IC)
      DO 25 M=1,NHA
        DO 24 N=1,NSI
          PROJ(N,M)=0.
24      CONTINUE
25    CONTINUE
      PROJ(5,1)=IC(1)
      PROJ(6,1)=IC(2)
      PROJ(5,2)=IC(1)
      PROJ(6,2)=IC(2)
      PROJ(1,2)=ARP(1)
      PROJ(2,2)=ARP(2)
      PROJ(3,2)=ARP(3)
      PROJ(4,2)=ARP(4)
      CALL IDTR4(IDM,IC)
      DO 28 M=1,NHA
        DO 27 N=1,NSI
          TARG(N,M)=0.
27      CONTINUE
28    CONTINUE
      TARG(5,1)=IC(1)
      TARG(6,1)=IC(2)
      TARG(5,2)=IC(1)
      TARG(6,2)=IC(2)
      TARG(1,2)=ARM(1)
      TARG(2,2)=ARM(2)
      TARG(3,2)=ARM(3)
      TARG(4,2)=ARM(4)

C  REDO
C  ----
      CALL UTREMB(PROJ,TARG,1)
      LOO=1
      GOTO 58
57    CONTINUE
      LOO=LOO+1
      IF ( LOO .GT. 20 ) GOTO 10025
      IF ( ISH .GE. 91 ) WRITE(IFCH,*)'REDO HH COLLISION'
      CALL UTREST(PROJ,TARG,1)
58    CONTINUE

C  NUMBER OF CE'S
C  --------------
      NCOLEX=1
C-C   LO=0
C-C16 LO=LO+1
C-C   IF ( LO .EQ. 3 ) THEN
C-C     IF ( ISH .GE. 90 ) THEN
C-C       CALL UTMSG('JINTCE')
C-C       WRITE(IFCH,*)'*****  LO=3'
C-C       CALL UTMSGF
C-C     ENDIF
C-C   ENDIF
C-C   R=RANGEN()
C-C   NCOLEX=0
C-C15 NCOLEX=NCOLEX+1
C-C   IF ( NCOLEX .GT. NPRBMS ) GOTO 16
C-C   IF ( R .GT. PRBMS(NCOLEX) ) GOTO 15

C  DO NCOLEX CE'S
C  --------------
      ISKIP=0
      DO 30 NCE=1,NCOLEX
        IF ( ISH .GE. 92 .AND. NCE .GT. 1 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)NCE,'. COLOUR EXCHANGE'
          WRITE(IFCH,*)' '
        ENDIF
        NCES=NCE
        ISH=ISH-2
        CALL HAHABS(PROJ,TARG,NCE/NCOLEX,NCE/NCOLEX,ISKIP,IRETHH)
        ISH=ISH+2
        IF ( ISKIP .NE. 0 ) GOTO 10025
        IF ( IRETHH .EQ. 1 ) GOTO 57
30    CONTINUE

C  FRAGMENTATION
C  -------------
      IF ( NSTR .LE. NSTR00 ) GOTO 10050

      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,123)
123     FORMAT(/,1X,'STRINGS BEFORE RESCALING',/)
        DO 31 L=NSTR00+1,NSTR
          WRITE(IFCH,109)L,(ICSTR(K,L)/100,K=1,4)
     *          ,SQRT(PSTR(1,L)**2+PSTR(2,L)**2),PSTR(3,L),PSTR(5,L)
     *          ,IRLSTR(L)
109       FORMAT(' /CSTR/',I5,3X,4I5,3(E11.3),I5)
31      CONTINUE
      ENDIF

      IF ( IRESCL .EQ. 1 ) THEN
        PSUM(1)=0.
        PSUM(2)=0.
        PSUM(3)=0.
        PSUM(4)=SROOT
        PSUM(5)=SROOT
        ISH=ISH-2
        CALL HRESCL(NSTR00+1,NSTR,PSUM,IFAIL)
        ISH=ISH+2
        IF ( IFAIL .NE. 0 ) GOTO 57
      ENDIF

      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,124)
124     FORMAT(/,1X,'STRINGS AFTER RESCALING',/)
        DO 32 L=NSTR00+1,NSTR
          WRITE(IFCH,109)L,(ICSTR(K,L)/100,K=1,4)
     *          ,SQRT(PSTR(1,L)**2+PSTR(2,L)**2),PSTR(3,L),PSTR(5,L)
     *          ,IRLSTR(L)
32      CONTINUE
      ENDIF

      DO 33 L=NSTR00+1,NSTR
        IF ( ISH .GE. 92 ) WRITE(IFCH,122)L
     *             ,(ICSTR(K,L),K=1,4)
     *             ,SQRT(PSTR(1,L)**2+PSTR(2,L)**2),PSTR(3,L),PSTR(5,L)
122     FORMAT(/,' STR:',I7,4I7,3(E10.2),/)
        ISH=ISH-2
        CALL JAMFRA(L,NEWEVT)
        ISH=ISH+2
        IF ( NEWEVT .EQ. 1 ) GOTO 10025
33    CONTINUE

C  BOOST PARTICLES FROM STRING FRAGMENTATION
C  -----------------------------------------
      IF ( NPTL .LE. NPTL00 ) GOTO 10051

      IORPTL(NPTL00+1)=I
      JORPTL(NPTL00+1)=J
      DO 34 L=1,5
        PSUM(L)=0.
34    CONTINUE

      DO 38 N=NPTL00+1,NPTL

        DO 35 L=1,5
          ARP(L)=PPTL(L,N)
35      CONTINUE
        CALL UTROT2(-1,ROT(1),ROT(2),ROT(3)
     *                ,ARP(1),ARP(2),ARP(3))
        CALL UTLOB2(-1,BOO(1),BOO(2),BOO(3),BOO(4),BOO(5)
     *                ,ARP(1),ARP(2),ARP(3),ARP(4))
        DO 36 L=1,5
          PPTL(L,N)=ARP(L)
36      CONTINUE
        NQJPTL(N)=0

        IF ( ISTPTL(N) .EQ. 0 ) THEN
          PSUM(1)=PSUM(1)+PPTL(1,N)
          PSUM(2)=PSUM(2)+PPTL(2,N)
          PSUM(3)=PSUM(3)+PPTL(3,N)
          PSUM(4)=PSUM(4)+PPTL(4,N)
        ENDIF

        IF ( ISH .GE. 91 ) THEN
          WRITE(IFCH,115)N,IDPTL(N)
     *              ,(PPTL(K,N),K=3,5),(XORPTL(K,N),K=3,4),ISTPTL(N)
        ENDIF

38    CONTINUE

      IF ( ISH .GE. 90 ) THEN
        IF ( ISH .GE. 92 ) THEN
          WRITE(IFCH,*)'P_FIN:',(PSUM(K),K=1,4)
          WRITE(IFCH,*)'P_INI:',(BOO(K),K=1,4)
        ENDIF

        IF ( (ABS(BOO(1)-PSUM(1)) .GT. 5.E-3*ABS(BOO(1))
     *         .OR. ABS(BOO(2)-PSUM(2)) .GT. 5.E-3*ABS(BOO(2))
     *         .OR. ABS(BOO(3)-PSUM(3)) .GT. 5.E-3*ABS(BOO(3))
     *         .OR. ABS(BOO(4)-PSUM(4)) .GT. 1.E-1*ABS(BOO(4)))
     *    .AND.
     *       (ABS(BOO(1)-PSUM(1)) .GT. 5.E-3
     *         .OR. ABS(BOO(2)-PSUM(2)) .GT. 5.E-3
     *         .OR. ABS(BOO(3)-PSUM(3)) .GT. 5.E-3
     *         .OR. ABS(BOO(4)-PSUM(4)) .GT. 1.E-1) ) THEN
          CALL UTMSG('JINTCE')
          WRITE(IFCH,*)'*****  P_INI /= P_FIN'
          WRITE(IFCH,*)'FINAL PARTICLES:'
          DO 39 N=NPTL00+1,NPTL
            WRITE(IFCH,125)N,IDPTL(N),(PPTL(K,N),K=1,5),ISTPTL(N)
125         FORMAT(1X,'/CPTL/',I6,I10,5(E11.3),I2)
39        CONTINUE
          WRITE(IFCH,*)'P_FIN:',(PSUM(K),K=1,4)
          WRITE(IFCH,*)'P_INI:',(SNGL(BOO(K)),K=1,4)
          CALL UTMSGF
        ENDIF
      ENDIF
      GOTO 1000

C  FINISH
C  ------

10025 CONTINUE
      IRET=25
      IF ( ISH .GE. 91 ) WRITE(IFCH,*)'SKIP'
      NPTL=NPTL00
      NSTR=NSTR00
      GOTO 1000

10050 CONTINUE
      CALL UTMSG('JINTCE')
      WRITE(IFCH,*)'INCIDENT PARTICLES:'
      WRITE(IFCH,115)I,IDPTL(I),(PPTL(K,I),K=3,5)
      WRITE(IFCH,115)J,IDPTL(J),(PPTL(K,J),K=3,5)
      WRITE(IFCH,*)'NSTR=',NSTR,'    NSTR00=',NSTR00
      CALL UTMSGF
      CALL UTSTOP('JINTCE: NSTR .LE. NSTR00                ')

10051 CONTINUE
      CALL UTMSG('JINTCE')
      WRITE(IFCH,*)'INCIDENT PARTICLES:'
      WRITE(IFCH,115)I,IDPTL(I),(PPTL(K,I),K=3,5)
      WRITE(IFCH,115)J,IDPTL(J),(PPTL(K,J),K=3,5)
      WRITE(IFCH,*)'NSTR=',NSTR,'    NSTR00=',NSTR00
      WRITE(IFCH,*)'NPTL=',NPTL,'    NPTL00=',NPTL00
      CALL UTMSGF
      CALL UTSTOP('JINTCE: NPTL .LE. NPTL00                ')

1000  CONTINUE
      PNLLX=PNLLX0
      ISH=ISH00
      RETURN
      END
C=======================================================================

      SUBROUTINE JINTCH(I,J,KMAX)

C----------------------------------------------------------------------
C  WRITES CHILDREN OF I,J TO IFRIJ(1-KMAX)
C----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (MXIFR=MXPTL)
      COMMON /CIFRIJ/  IFRIJ(MXIFR)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C----------------------------------------------------------------------
      KCUR=0
      KMAX=0
      LOOP=0
31    LOOP=LOOP+1
      IF     ( LOOP .EQ. 1 ) THEN
        N=I
      ELSEIF ( LOOP .EQ. 2 ) THEN
        N=J
      ENDIF
      IF ( IFRPTL(1,N) .GT. 0 ) THEN
        IF ( ISH.GE.92 ) WRITE(IFCH,*)N,' ---> ',IFRPTL(1,N),IFRPTL(2,N)
        DO 30 K=IFRPTL(1,N),IFRPTL(2,N)
          KMAX=KMAX+1
          IF ( KMAX .GT. MXIFR ) THEN
            CALL UTSTOP('JINTCH: KMAX > MXIFR                    ')
          ENDIF
          IFRIJ(KMAX)=K
30      CONTINUE
      ENDIF
      IF ( LOOP .EQ. 1 ) GOTO 31
      KCUR=KCUR+1
      IF ( KCUR .LE. KMAX ) THEN
        N=IFRIJ(KCUR)
        GOTO 31
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE JINTCL(I,X,Y,RAD,O,U,IAC)

C---------------------------------------------------------------------
C  FILLS HISTOGRAM CONCERNING CLUSTER CHARACTERISTICS
C---------------------------------------------------------------------
      PARAMETER (MXEPS=10)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXTAU=4)
      PARAMETER (MXVOL=10)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CJINTC/  CLUST(MXTAU,MXVOL,MXEPS)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CVSN/    IVERSN
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /NEVNT/   NEVNT
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARO5/   DELEPS,DELVOL

      CHARACTER AX*2
      LOGICAL   IAC
      SAVE
CTP060202 TO AVOID WARNINGS WITH GFORTRAN COMPILATION
      LOGICAL CTP060202
      CTP060202=.FALSE.
      IF(CTP060202)WRITE(*,*)X,Y
CTP060202 END WARNING
C---------------------------------------------------------------------
      TAU=SNGL(TTAUS)
CDH   ITAU=1+LOG(TAU)/LOG(2.)
      ITAU=1+LOG(TAU)*1.442695

      IF ( IAC ) GOTO 1

      VOL=(O-U)*RAD**2*PI
      IVOL=1+VOL/DELVOL
      EPS=PPTL(5,I)/VOL
      IEPS=1+EPS/DELEPS
      IF ( ITAU.GE.1 .AND. ITAU.LE.MXTAU .AND.
     *     IVOL.GE.1 .AND. IVOL.LE.MXVOL .AND.
     *     IEPS.GE.1 .AND. IEPS.LE.MXEPS )
     *                CLUST(ITAU,IVOL,IEPS)=CLUST(ITAU,IVOL,IEPS)+1

 1    CONTINUE
      IF ( I .LT. NPTL  .OR.  NREVT+1 .LT. NEVNT ) GOTO 1000

      IF ( ITAU .EQ. 1 ) WRITE(IFCH,105)MAPROJ,MATARG,ENGY,IVERSN/100.
105   FORMAT(/,1X,'PROJ=',I3,3X,'TARG=',I3,3X,'ENGY=',F7.2
     *          ,3X,'VENUS ',F4.2,' (TURBOVENUS)')
      WRITE(IFCH,100)TAU,ITAU,NEVNT,DELVOL,DELEPS
100   FORMAT(/,1X,'TAU=',F5.2,3X,'ITAU=',I1
     *          ,3X,'NEVNT=',I4,3X,'DELVOL=',F5.1,3X,'DELEPS=',F5.1,/)
      WRITE(IFCH,101)
101   FORMAT(9X,'IVOL=1 IVOL=2 IVOL=3 IVOL=4 IVOL=5 '
     *         ,'IVOL=6 IVOL=7 IVOL=8 IVOL=9 IVOL=10 ')
      DO 106 IE=1,MXEPS
        IF ( IE .LE. 9 ) THEN
          WRITE(AX,102)IE
102       FORMAT('0',I1)
        ELSE
          WRITE(AX,103)IE
103       FORMAT(I2)
        ENDIF
        WRITE(IFCH,104)AX,(NINT(CLUST(ITAU,IV,IE)),IV=1,MXVOL)
104     FORMAT(' IEPS=',A2,10I7)
106   CONTINUE

1000  CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE JINTEL(I,J,AMIM,IACTN)

C----------------------------------------------------------------------
C  ELASTIC SCATTERING OF PTLS I,J
C  EL SCATT REQUIR NOT FULF: RETURN WITHOUT ACTION, IACTN UNCHANGED
C    ELSE: ELASTIC SCATTERING DONE, IACTN=2
C----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (NPTQ=129)
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CPTQ/    QPTH(NPTQ),QPTQ(NPTQ),XPTQ(NPTQ),QPTQMX,QPTHMX
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL             P(5),PEI(5),PEJ(5),U(5)
      SAVE
C----------------------------------------------------------------------
C  CHECK
C  -----
      AMF=PPTL(5,NPTL+1)
      IF ( AMF .GE. AMIM ) GOTO 1000

C  PRINT
C  -----
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,101)NTC,SNGL(TTAUS)
101     FORMAT(1X,'ELASTIC SCATTERING --- NT=',I3,' --- TAUS=',F6.2)
        WRITE(IFCH,115)I,IDPTL(I)
     *                  ,(PPTL(K,I),K=3,5),(XORPTL(K,I),K=3,4)
        WRITE(IFCH,115)J,IDPTL(J)
     *                  ,(PPTL(K,J),K=3,5),(XORPTL(K,I),K=3,4)
115     FORMAT(1X,'/CPTL/',I6,I10
     *           ,1X,2(E10.2),E11.2,1X,2(E10.2),2X,I4)
      ENDIF

C  INITIALIZATION
C  --------------
      IACTN=2
      DO 125 K=1,5
        P(K)=PPTL(K,NPTL+1)
125   CONTINUE

C  DETERMINE MOMENTA OF OUTGOING PARTICLES (PEI,PEJ)
C  -------------------------------------------------
      IF     ( P(5) .LE. (PPTL(5,I)+PPTL(5,J))*.99 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTEL')
          WRITE(IFCH,132)P(5),PPTL(5,I)+PPTL(5,J)
132       FORMAT(1X,'*****  M_FUS < M_I+M_J ---> QCM SET ZERO    ( '
     *             ,2F10.3,' )')
          WRITE(IFCH,133)'P_I:  ',(PPTL(K,I),K=1,5)
          WRITE(IFCH,133)'P_J:  ',(PPTL(K,J),K=1,5)
          WRITE(IFCH,133)'P_FUS:',(P(K),K=1,5)
133       FORMAT(1X,A6,3X,5F10.3)
          CALL UTMSGF
        ENDIF
        QCM=0.
      ELSEIF ( P(5) .LE. PPTL(5,I)+PPTL(5,J) ) THEN
        QCM=0.
      ELSE
        QCM=UTPCM(P(5),PPTL(5,I),PPTL(5,J))
      ENDIF

C  ISOTROPIC
      U(3)=2.*RANGEN()-1.
      PHI=2.*PI*RANGEN()
      AUXIL=SQRT(1.-U(3)**2)
      U(1)=AUXIL*COS(PHI)
      U(2)=AUXIL*SIN(PHI)
      PEI(1)= QCM*U(1)
      PEJ(1)=-QCM*U(1)
      PEI(2)= QCM*U(2)
      PEJ(2)=-QCM*U(2)
      PEI(3)= QCM*U(3)
      PEJ(3)=-QCM*U(3)

C  NONISOTROPIC
C-C   R=RANGEN()
C-C   IF     ( IOPTQ .EQ. 2 ) THEN
C-C     PT = SQRT( -4.*PTQ**2/PI * LOG(1.-QPTQMX*R) )
C-C   ELSEIF ( IOPTQ .EQ. 3 ) THEN
C-C     PT = PTQ*SQRT( QPTQMX*R/(1.-QPTQMX*R) )
C-C   ELSE
C-C     PT=UTINVT(NPTQ,XPTQ,QPTQ,R*QPTQ(NPTQ))
C-C   ENDIF
C-C   IF ( PT .GE. QCM ) PT=RANGEN()*QCM
C-C   QPL=SQRT(QCM**2-PT**2)
C-C   U(3)=QPL
C-C   PHI=2.*PI*RANGEN()
C-C   U(1)=PT*COS(PHI)
C-C   U(2)=PT*SIN(PHI)
C-C   CALL UTAXIS(I,J,A1,A2,A3)
C-C   IVT=1
C-C   IF ( A3 .LT. 0. ) THEN
C-C     A1=-A1
C-C     A2=-A2
C-C     A3=-A3
C-C     IVT=-1
C-C   ENDIF
C-C   CALL UTROTA(-1,A1,A2,A3,U(1),U(2),U(3))
C-C   DO 47 K=1,3
C-C     PEI(K)= U(K)*IVT
C-C     PEJ(K)=-U(K)*IVT
C-C47 CONTINUE
      PEI(4)=SQRT(QCM**2+PPTL(5,I)**2)
      PEJ(4)=SQRT(QCM**2+PPTL(5,J)**2)
      PEI(5)=PPTL(5,I)
      PEJ(5)=PPTL(5,J)
      CALL UTLOBO(-1,P(1),P(2),P(3),P(4),P(5)
     *            ,PEI(1),PEI(2),PEI(3),PEI(4))
      CALL UTLOBO(-1,P(1),P(2),P(3),P(4),P(5)
     *            ,PEJ(1),PEJ(2),PEJ(3),PEJ(4))

C  FILL /CPTL/
C  -----------
      DO 49 LO=1,2
        NPTL=NPTL+1
        IF ( LO .EQ. 1 ) THEN
          IJ=I
          PPTL(1,NPTL)=PEI(1)
          PPTL(2,NPTL)=PEI(2)
          PPTL(3,NPTL)=PEI(3)
          PPTL(4,NPTL)=PEI(4)
          PPTL(5,NPTL)=PEI(5)
        ELSE
          IJ=J
          PPTL(1,NPTL)=PEJ(1)
          PPTL(2,NPTL)=PEJ(2)
          PPTL(3,NPTL)=PEJ(3)
          PPTL(4,NPTL)=PEJ(4)
          PPTL(5,NPTL)=PEJ(5)
        ENDIF
        ISTPTL(NPTL) =0
        IDPTL(NPTL)  =IDPTL(IJ)
        IBPTL(1,NPTL)=IBPTL(1,IJ)
        IBPTL(2,NPTL)=IBPTL(2,IJ)
        IBPTL(3,NPTL)=IBPTL(3,IJ)
        IBPTL(4,NPTL)=IBPTL(4,IJ)
        XORPTL(1,NPTL)=XAVER(1)
        XORPTL(2,NPTL)=XAVER(2)
        XORPTL(3,NPTL)=XAVER(3)
        XORPTL(4,NPTL)=XAVER(4)
        IORPTL(NPTL)=I
        JORPTL(NPTL)=J
        TIVPTL(1,NPTL)=XAVER(4)
        CALL IDTAU(IDPTL(NPTL),PPTL(4,NPTL),PPTL(5,NPTL),TAUGM)
        TIVPTL(2,NPTL)=TIVPTL(1,NPTL)+TAUGM
        IFRPTL(1,NPTL)=0
        IFRPTL(2,NPTL)=0
        ICLPTL(NPTL)=ICLPTL(IJ)
        NQJPTL(NPTL)=NQJPTL(IJ)
        IF ( ISH .GE. 91 ) WRITE(IFCH,115)NPTL,IDPTL(NPTL)
     *                    ,(PPTL(K,NPTL),K=3,5),(XORPTL(K,NPTL),K=3,4)
49    CONTINUE

1000  RETURN
      END
C=======================================================================

      SUBROUTINE JINTFP(I,J,N,     XPL1,YPL1,RADPL1,OPL1,UPL1,
     *  XPL2,YPL2,RADPL2,OPL2,UPL2,XPL3,YPL3,RADPL3,OPL3,UPL3)

C----------------------------------------------------------------------
C  PRINTOUT
C----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP

      REAL       OPL(3),RADPL(3),UPL(3),XPL(3),YPL(3)
      CHARACTER  AX*10,LIN*59,MARK*1,TXT*8
      SAVE
C----------------------------------------------------------------------
      XPL(1)  =XPL1
      YPL(1)  =YPL1
      RADPL(1)=RADPL1
      OPL(1)  =OPL1
      UPL(1)  =UPL1
      XPL(2)  =XPL2
      YPL(2)  =YPL2
      RADPL(2)=RADPL2
      OPL(2)  =OPL2
      UPL(2)  =UPL2
      XPL(3)  =XPL3
      YPL(3)  =YPL3
      RADPL(3)=RADPL3
      OPL(3)  =OPL3
      UPL(3)  =UPL3

      V1=(OPL1-UPL1)*RADPL1**2*PI
      V2=(OPL2-UPL2)*RADPL2**2*PI
      V3=(OPL3-UPL3)*RADPL3**2*PI

      WRITE(IFCH,103)SNGL(TTAUS)
103   FORMAT(' MONITORING FUSION INTERACTION AT TTAUS=',F5.2)
      TX=TTAUS
      TX=MAX(TX,1.)
      WRITE(TXT(3:7),104)TX
104   FORMAT(F5.2)
      WRITE(IFCH,105)'INCOMING:',I,IDPTL(I),PPTL(5,I),V1,PPTL(5,I)/V1
      WRITE(IFCH,105)'INCOMING:',J,IDPTL(J),PPTL(5,J),V2,PPTL(5,J)/V2
      WRITE(IFCH,105)'FUSED:   ',N,IDPTL(N),PPTL(5,N),V3,PPTL(5,N)/V3
105   FORMAT(1X,A9,2X,'NR=',I6,2X,'ID=',I9,2X,'MASS=',F5.1,2X
     *         ,'VOL=',F5.1,2X,'EPS=',F4.1)

      DO 110 MM=1,3
        IF     ( MM .EQ. 1 ) THEN
          TXT(1:2)='S/'
          TXT(8:8)=':'
        ELSEIF ( MM .EQ. 2 ) THEN
          TXT='X:      '
        ELSEIF ( MM .EQ. 3 ) THEN
          TXT='Y:      '
        ENDIF
        AX='---------!'
        WRITE(IFCH,100)TXT,(AX,L=1,6)
100     FORMAT(2X,A8,3X,'-3',8X,'-2',8X,'-1',9X,'0',9X,'1',9X,'2',9X,'3'
     *         ,/,13X,' !',6A10)
        DO 109 K=1,3
          DO 108 L=1,59
            WO=-3.+L*0.1 + 0.05
            WU=-3.+L*0.1 - 0.05
            IF     ( MM .EQ. 1 ) THEN
              O=OPL(K)/TX
              U=UPL(K)/TX
            ELSEIF ( MM .EQ. 2 ) THEN
              O=XPL(K)+RADPL(K)
              U=XPL(K)-RADPL(K)
            ELSEIF ( MM .EQ. 3 ) THEN
              O=YPL(K)+RADPL(K)
              U=YPL(K)-RADPL(K)
            ENDIF
            IF ( K .EQ. 3 ) THEN
              MARK='X'
            ELSE
              MARK='O'
            ENDIF
            IF ( WU .LE. O  .AND.  WO .GE. U ) THEN
              LIN(L:L)=MARK
            ELSE
              LIN(L:L)=' '
            ENDIF
108       CONTINUE
          WRITE(IFCH,101)U,O,LIN
101       FORMAT(1X,2F6.2,' !',A59,'!')
109     CONTINUE
        WRITE(IFCH,102)('-',L=1,59)
102     FORMAT(13X,' !',59A1,'!')
110   CONTINUE

      RETURN
      END
C=======================================================================

      SUBROUTINE JINTFS(I,J,NQI,JC,AMIM,IRET)

C----------------------------------------------------------------------
C  INPUT: PTL NUMBERS I,J
C  OUTPUT: PPFUS(5): MOMENTUM , NQI: NET QUARK NUMBER,
C              JC: JC-CODE, AMIM: MINIMUM MASS,    OF FUSED PTL
C          IRET=0 IF OK, 1 ELSE
C          PPFUS() WRITTEN TO PPTL(,NPTL+1)
C----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (NFLAV=6)
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION PPFUS(5),PP52
      REAL             P(5)
      INTEGER          JC(NFLAV,2),JCI(NFLAV,2),JCJ(NFLAV,2)
      SAVE
C----------------------------------------------------------------------
      IRET=0

      P(1)=PPTL(1,I)+PPTL(1,J)
      PPFUS(1)=P(1)
      P(2)=PPTL(2,I)+PPTL(2,J)
      PPFUS(2)=P(2)
      P(3)=PPTL(3,I)+PPTL(3,J)
      PPFUS(3)=P(3)
      P(4)=PPTL(4,I)+PPTL(4,J)
      PPFUS(4)=P(4)
      PP52=PPFUS(4)**2-PPFUS(3)**2-PPFUS(2)**2-PPFUS(1)**2
      IF ( PP52 .LE. 0.D0 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JINTFS')
          WRITE(IFCH,*)'*****  MFUS**2 < 0    (',PP52,' )'
          WRITE(IFCH,*)(PPFUS(M),M=1,4)
          CALL UTMSGF
        ENDIF
        GOTO 1001
      ENDIF
      PPFUS(5)=SQRT(PP52)
      IF ( NPTL+1 .GT. MXPTL ) THEN
        CALL UTSTOP('JINTFS: NPTL>MXPTL                      ')
      ENDIF
      DO 36 K=1,5
        PPTL(K,NPTL+1)=PPFUS(K)
36    CONTINUE
      AMF=PPFUS(5)

      ISTPTL(NPTL+1)=0

      CALL IDQUAC(I,NDUMMY,NDUMMY,NDUMMY,JCI)
      CALL IDQUAC(J,NDUMMY,NDUMMY,NDUMMY,JCJ)
      NQI=0
      DO 29 N=1,NFLAV
        JC(N,1)=JCI(N,1)+JCJ(N,1)
        JC(N,2)=JCI(N,2)+JCJ(N,2)
        NQI=NQI+JC(N,1)-JC(N,2)
29    CONTINUE

      CALL IDCOMJ(JC)
      AMIM=UTAMNZ(JC,5)+.200
      IF ( AMF.LT.AMIM .AND. I.GT.NPTL0 .AND. J.GT.NPTL0 ) GOTO 1001
      GOTO 1000

1001  IRET=1

1000  RETURN
      END
C=======================================================================

      SUBROUTINE JINTFU(I,J,JC,IACTN)

C----------------------------------------------------------------------
C  FUSION OF PTLS I,J: DETERMINE CLUSTER
C----------------------------------------------------------------------
      PARAMETER (MXDKY=2000)
      PARAMETER (MXLOOK=10000)
      PARAMETER (MXPTL=70000)
      PARAMETER (NFLAV=6)
      COMMON /CJINT/   BX,BY,RNUP,RNUT,VELP,VELT,XAVER(4),NPTL0
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CTIMEL/  NTC
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /DKYTAB/  CBR(MXDKY),LOOK(MXLOOK),MODE(5,MXDKY)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      INTEGER IB(4),IC(2),JC(NFLAV,2)
      SAVE
C----------------------------------------------------------------------
C  PRINT
C  -----
      IF ( ISH .GE. 91 ) THEN
        WRITE(IFCH,101)NTC,SNGL(TTAUS)
101     FORMAT(1X,'FUSION INTERACTION --- NT=',I3,' --- TAUS=',F6.2)
        WRITE(IFCH,115)I,IDPTL(I)
     *                  ,(PPTL(K,I),K=3,5),(XORPTL(K,I),K=3,4)
        WRITE(IFCH,115)J,IDPTL(J)
     *                  ,(PPTL(K,J),K=3,5),(XORPTL(K,I),K=3,4)
115     FORMAT(1X,'/CPTL/',I6,I10
     *           ,1X,2(E10.2),E11.3,1X,2(E10.2),2X,I4)
      ENDIF

C  INITIALIZATION
C  --------------
      IACTN=3
      AMF=PPTL(5,NPTL+1)
      NPTL=NPTL+1

C  DETERMINE IDR, IB(1-4)
C  ----------------------
      IDR=0
      DO 40 NF=1,NFLAV
        IF ( JC(NF,1) .GE. 10 ) IDR=700000000
        IF ( JC(NF,2) .GE. 10 ) IDR=700000000
40    CONTINUE
      IF ( IDR/100000000 .NE. 7 ) THEN
        CALL IDENCO(JC,IC,IRETEN)
        IF ( IRETEN .EQ. 1 ) THEN
          CALL UTSTOP('JINTFU: IDENCO RET CODE = 1             ')
        ENDIF
        ID=IDTRA(IC,0,0,3)
43      AMC=AMF
        CALL IDRES(ID,AMC,IDR,IADJ)
        IF ( IDR .NE. 0 ) THEN
          LID=LOOK(IABS(IDR))
          IF ( LID.LE.0  .OR.  LID.GT.0 .AND. MODE(2,LID).EQ.0 ) THEN
            IF ( PPTL(5,NPTL) .GT. AMC+1.E-3 ) THEN
              AMF=AMF+0.010
              GOTO 43
            ENDIF
            IF ( ABS(AMC-PPTL(5,NPTL)) .GT. 1.E-3 ) THEN
              IF ( ISH .GE. 90 ) THEN
                CALL UTMSG('JINTFU')
                WRITE(IFCH,*)'*****  NOT ON MASS SHELL AFTER FUSION: '
     *                       ,PPTL(5,NPTL),AMC
                CALL UTMSGF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF ( IDR .EQ. 0 ) THEN
          IF ( MOD(IC(1),100).NE.0 .OR. MOD(IC(2),100).NE.0 ) THEN
            IDR=900000000
          ELSE
            IDR=800000000+IC(1)*100+IC(2)/100
          ENDIF
        ENDIF
      ELSE
        CALL IDTRBI(JC,IB(1),IB(2),IB(3),IB(4))
        IDR=IDR
     *      +MOD(JC(1,1)+JC(2,1)+JC(3,1)+JC(4,1),10000)*10000
     *      +MOD(JC(1,2)+JC(2,2)+JC(3,2)+JC(4,2),10000)
        IF ( ISH .GE. 93 ) WRITE(IFCH,*) 'IB:',(IB(KK),KK=1,4)
        IBPTL(1,NPTL)=IB(1)
        IBPTL(2,NPTL)=IB(2)
        IBPTL(3,NPTL)=IB(3)
        IBPTL(4,NPTL)=IB(4)
      ENDIF

C  FILL /CPTL/
C  -----------
      IDPTL(NPTL)=IDR
      XORPTL(1,NPTL)=XAVER(1)
      XORPTL(2,NPTL)=XAVER(2)
      XORPTL(3,NPTL)=XAVER(3)
      XORPTL(4,NPTL)=XAVER(4)
      IORPTL(NPTL)=I
      JORPTL(NPTL)=J
      TIVPTL(1,NPTL)=XAVER(4)
      CALL IDTAU(IDPTL(NPTL),PPTL(4,NPTL),PPTL(5,NPTL),TAUGM)
      TIVPTL(2,NPTL)=TIVPTL(1,NPTL)+TAUGM
      IFRPTL(1,NPTL)=0
      IFRPTL(2,NPTL)=0
      ICLPTL(NPTL)=1
      NQJPTL(NPTL)=0

C  PRINT + RETURN
C  --------------
      IF ( ISH .GE. 91 ) THEN
        N=NPTL
        WRITE(IFCH,115)N,IDPTL(N)
     *                  ,(PPTL(K,N),K=3,5),(XORPTL(K,N),K=3,4)
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE JINTPA(I,J,KMAX)

C----------------------------------------------------------------------
C  WRITES PARENTS OF I,J TO IFRIJ(1-KMAX)
C  SETS IAAPTL()=0 FOR PARENTS
C----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (MXIFR=MXPTL)
      COMMON /CIFRIJ/  IFRIJ(MXIFR)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C----------------------------------------------------------------------
      KCUR=0
      KMAX=0
      LOOP=0
12    LOOP=LOOP+1
      IF     ( LOOP .EQ. 1 ) THEN
        N=I
      ELSEIF ( LOOP .EQ. 2 ) THEN
        N=J
      ENDIF
      IF ( IORPTL(N) .GT. 0 ) THEN
        IF ( ISH .GE. 92 ) WRITE(IFCH,*)N,' <--- ',IORPTL(N),JORPTL(N)
        IF ( KMAX+2 .GT. MXIFR ) THEN
          CALL UTSTOP('JINTPA: KMAX+2 > MXIFR                  ')
        ENDIF
        IF ( IAAPTL(IORPTL(N)) .NE. 0 ) THEN
          KMAX=KMAX+1
          IFRIJ(KMAX)=IORPTL(N)
          IAAPTL(IORPTL(N))=0
        ENDIF
        IF ( JORPTL(N) .GT. 0 ) THEN
          IF ( IAAPTL(JORPTL(N)) .NE. 0 ) THEN
            KMAX=KMAX+1
            IFRIJ(KMAX)=JORPTL(N)
            IAAPTL(JORPTL(N))=0
          ENDIF
        ENDIF
      ENDIF
      IF ( LOOP .EQ. 1 ) GOTO 12
ctp060203 8     KCUR=KCUR+1
      KCUR=KCUR+1
      IF ( KCUR .LE. KMAX ) THEN
        N=IFRIJ(KCUR)
        IF ( IAAPTL(N) .NE. 0 ) THEN
          CALL UTSTOP('JINTPA: SHOULD NOT HAPPEN               ')
        ENDIF
        GOTO 12
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE JRESCL(J1,J2,PSUM,IFAIL)

C-----------------------------------------------------------------------
C  RESCALES PTL MOMENTA OF PTLS J1-J2 TO HAVE TOTAL MOM PSUM.
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CSCAL/   SCAL
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION ENE,PADD(5),PP(5),PPSUM(5)
      REAL             PSUM(5)
      DATA ERRLIM/.001/
      SAVE
C-----------------------------------------------------------------------
      IFAIL=1

      IF ( J1 .GE. J2 ) THEN
        CALL UTSTOP('JRESCL: J1 .GE. J2                      ')
      ENDIF

      DO 100 K=1,5
        PPSUM(K)=PSUM(K)
        PADD(K)=0.D0
100   CONTINUE
      DO 110 J=J1,J2
        DO 110 K=1,5
          PADD(K)=PADD(K)+PPTL(K,J)
110   CONTINUE
      IF ( PADD (5) .GE. PPSUM(5) ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JRESCL')
          WRITE(IFCH,*)'*****  SUM OF PTL MASSES .GE. PPSUM(5)'
          DO 1 J=J1,J2
            WRITE(IFCH,109)J,IDPTL(J),(PPTL(L,J),L=3,5)
109         FORMAT(' /CPTL/',I6,I10,3(E11.3))
 1        CONTINUE
          WRITE(IFCH,*)'PPSUM(345):',(SNGL(PPSUM(K)),K=3,5)
          CALL UTMSGF
        ENDIF
        RETURN
      ENDIF
      PADD(5)=PADD(4)**2-PADD(1)**2-PADD(2)**2-PADD(3)**2
      IF ( PADD(5) .LE. 0.D0 ) THEN
        ENE = 0.D0
        DO 111 J=J1,J2
          ENE = ENE + DSQRT( PPTL(1,J)**2 + PPTL(2,J)**2
     *                + DBLE(PPTL(3,J))**2 + PPTL(5,J)**2 )
111     CONTINUE
        PADD(5) = ENE**2 - PADD(1)**2 - PADD(2)**2 - PADD(3)**2
        IF ( PADD(5) .LE. 0.D0 ) THEN
          DO 2 J=J1,J2
            WRITE(IFCH,108)J,(PPTL(L,J),L=1,5)
108         FORMAT(' /CPTL/',I4,5(E11.3))
 2        CONTINUE
          CALL UTSTOP('JRESCL: MASS**2 OF STRING-SUM NEGATIVE  ')
        ENDIF
      ENDIF
      PADD(5)=SQRT(PADD(5))

C  BOOST PTLS TO REST
C  ------------------
      DO 200 J=J1,J2
        PP(1)=PPTL(1,J)
        PP(2)=PPTL(2,J)
        PP(3)=PPTL(3,J)
        PP(4)=PPTL(4,J)
        CALL UTLOB2(1,PADD(1),PADD(2),PADD(3),PADD(4),PADD(5)
     *                 ,PP(1),PP(2),PP(3),PP(4))
        PPTL(1,J)=PP(1)
        PPTL(2,J)=PP(2)
        PPTL(3,J)=PP(3)
        PPTL(4,J)=PP(4)
200   CONTINUE

C  RESCALE MOMENTA IN REST FRAME
C  -----------------------------
      SCAL=1.
      DO 301 IPASS=1,200
        SUM=0.
        DO 310 J=J1,J2
          PPTL(1,J)=SCAL*PPTL(1,J)
          PPTL(2,J)=SCAL*PPTL(2,J)
          PPTL(3,J)=SCAL*PPTL(3,J)
          PPTL(4,J)=SQRT(PPTL(1,J)**2+PPTL(2,J)**2+PPTL(3,J)**2
     *                  +PPTL(5,J)**2)
          SUM=SUM+PPTL(4,J)
310     CONTINUE
        SCAL=PSUM(5)/SUM
        IF ( ABS(SCAL-1.) .LE. ERRLIM ) GOTO 300
301   CONTINUE
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('JRESCL')
        WRITE(IFCH,*)'*****  SCAL=',SCAL
        CALL UTMSGF
      ENDIF
300   CONTINUE

C  BOOST BACK WITH PPSUM
C  ---------------------
      DO 330 J=J1,J2
        PP(1)=PPTL(1,J)
        PP(2)=PPTL(2,J)
        PP(3)=PPTL(3,J)
        PP(4)=PPTL(4,J)
        CALL UTLOB2(-1,PPSUM(1),PPSUM(2),PPSUM(3),PPSUM(4),PPSUM(5)
     *                   ,PP(1),PP(2),PP(3),PP(4))
        PPTL(1,J)=PP(1)
        PPTL(2,J)=PP(2)
        PPTL(3,J)=PP(3)
        PPTL(4,J)=PP(4)
330   CONTINUE

      IFAIL=0
      RETURN
      END
C=======================================================================

      SUBROUTINE JSPLIT(STRO,STR,KOLSP,IER,KMAXOR)

C-----------------------------------------------------------------------
C  SPLITS STRING STRO INTO Q-QBAR STRING STR AND REMAINDER (->STRO)
C  DIMENSIONS: STRO(NSI,NSIX+1),STR(NSI,2)
C  IER=0: OK ; IER=1: ERROR ; IER=2: ABSORPTION ;
C  IER=3: AGAIN WITH NEW APART, EPART
C-----------------------------------------------------------------------
      PARAMETER (MAMX=56)
      PARAMETER (NFLAV=6)
      PARAMETER (NPTF=129)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      PARAMETER (NSPLIT=129)
      COMMON /CJSPLI/  ALEAD,APART,ELEAD,EPART,SGNSIL,JPART,NSCC,NSCCX
      COMMON /CNFUSN/  NFUSN(NSIX+1)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTF/    FPTFS,FPTFSS,FPTFU,FPTFUS,FPTFUU
     *                ,QPTFS(NPTF),QPTFSS(NPTF),QPTFU(NPTF),QPTFUS(NPTF)
     *                ,QPTFUU(NPTF),XPTF(NPTF)
      COMMON /CPZSTR/  ESTRL,PZSTRL,ISEA,ISTRL
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /QUARKM/  SMAS,SSMAS,USMAS,UUMAS

      DOUBLE PRECISION A,D,DAUXIL,PAM,PAP,PEM,PEP,PIM,PIP,PNL3,PNL4
     *                ,PM,PO(5),POM,POP,POT,POX,POY,PO123,PP
     *                ,PUM,PUP,PUT,PUX,PUY,PYM,PYP,SSTR(NSI,2)
      REAL             STR(NSI,2),STRO(NSI,NSIX+1),STRO0(NSI,NSIX+1)
      INTEGER          IC(2),ICX(2),JC(NFLAV,2),JCX(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      PUDX=PUD*.945
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 15 ) ISH=MOD(ISHSUB,100)

      CALL UTKSTR(STRO,KMAX)
      IF ( KMAX .EQ. KMAXOR ) THEN
        NSCC=0
        IF ( JPART .EQ. 0 ) THEN
          DO  5 I=1,KMAXOR
            NFUSN(I)=0
 5        CONTINUE
        ENDIF
      ENDIF
      KMAX0=KMAX
      DO 10 K=1,KMAX
        DO 10 I=1,NSI
          STRO0(I,K)=STRO(I,K)
10    CONTINUE
C-C   PDIQ=PDIQUA*0.5
C-C   IF ( KMAXOR .GT. KUTDIQ ) PDIQ=0.
C-C   PDIQ=PDIQUA*0.5*NSCCX/(KMAXOR-1.)
      PDIQ=0.
      LOOP=0
      NCORN=0
      XUNTER=0.
      XOBER=1.
      MESSCO=0
      IF ( SGNSIL .LT. 0. ) THEN
        KOLSP=KTARGA(KMAX+1,ITARG)
      ELSE
        KOLSP=KPROJA(KMAX+1,IPROJ)
      ENDIF
      PAP=2.D0*EPART
      PAM=0.D0
      PUX=STRO(1,KMAX)
      PUY=STRO(2,KMAX)
      PUT=SQRT(PUX**2+PUY**2)
      PUP=STRO(4,KMAX)-ABS(STRO(3,KMAX))
      PUM=STRO(4,KMAX)+ABS(STRO(3,KMAX))

C  PRINT
C  -----
      IF ( ISH .GE. 92 ) THEN
        IF ( ISH .GE. 93 ) WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)'MULTI-STRING DETECTED. SPLIT OFF Q-QBAR STRING:'
        IF ( ISH .GE. 93 ) WRITE(IFCH,*)('-',L=1,79)
        WRITE(IFCH,*)' '
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,*)'INPUT STRING STRO:'
          WRITE(IFCH,*)' '
          WRITE(IFCH,105)(STRO(I,1),I=1,4),(NINT(STRO(I,1)),I=5,6)
          DO 8 K=2,KMAX
            WRITE(IFCH,104)(STRO(I,K),I=1,4),(NINT(STRO(I,K)),I=5,6)
 8        CONTINUE
          WRITE(IFCH,*)' '
        ENDIF
      ENDIF

C  ABSORPTION
C  ----------
      IF ( ISTRL .EQ. 1 ) GOTO 1002

C  RESET
C  -----
ctp060203 5001  LOOP=0
      LOOP=0
5000  LOOP=LOOP+1
      IF ( LOOP .GE. 5 ) GOTO 1002
      KMAX=KMAX0
      DO 11 I=1,NSI
        DO 11 K=1,KMAX
        STRO(I,K)=STRO0(I,K)
11    CONTINUE

C  SPLIT OFF HADRON + DETERMINE REMAINDER STRING (MOMENTA)
C  -------------------------------------------------------
      IF ( RANGEN() .LT. PDIQ ) THEN
        NQU=2
      ELSE
        NQU=1
      ENDIF
      IFLTT=0
      IFLTO=0
      DO 24 I=1,NQU
        IFL=INT(RANGEN()/PUDX)+1
        IFLTO=IFLTO*10+IFL
        IFLTT=IFLTT*10+(IFL+1)/2
24    CONTINUE
      R=RANGEN()
      IF     ( IFLTT .EQ. 1 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFU ,R*QPTFU(NPTF))
C##       WRITE(IFCH,*)'JSPLIT:PT(OLD)=',PT
        ELSE
          RPT = R*FPTFU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(1.+RPT*2./AUXIL))
        ENDIF
      ELSEIF ( IFLTT .EQ. 2 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFS ,R*QPTFS(NPTF))
        ELSE
          RPT = R*FPTFS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP( SMAS**2/AUXIL)+RPT*2./AUXIL)- SMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 11 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUU,R*QPTFUU(NPTF))
        ELSE
          RPT = R*FPTFUU
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(UUMAS**2/AUXIL)+RPT*2./AUXIL)-UUMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 12  .OR.  IFLTT .EQ. 21 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFUS,R*QPTFUS(NPTF))
        ELSE
          RPT = R*FPTFUS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(USMAS**2/AUXIL)+RPT*2./AUXIL)-USMAS**2)
        ENDIF
      ELSEIF ( IFLTT .EQ. 22 ) THEN
        IF ( IOPTF .EQ. 1 ) THEN
          PT=UTINVT(NPTF,XPTF,QPTFSS,R*QPTFSS(NPTF))
        ELSE
          RPT = R*FPTFSS
          AUXIL=-4.*PTF**2/PI
          PT=SQRT(AUXIL*LOG(EXP(SSMAS**2/AUXIL)+RPT*2./AUXIL)-SSMAS**2)
        ENDIF
      ENDIF
      AT=SQRT(APART**2+PT**2)
      R=RANGEN()
      AUXIL=2.*R-1.
      IF ( AUXIL .LT. 0. ) THEN
        X = SQRT( 0.5-COS( (ACOS(-AUXIL)+PI)*.33333333 ) )
      ELSE
        X = SQRT( 0.5+COS( (ACOS(AUXIL)+PI)*.33333333 ) )
      ENDIF
      IF ( ESTRL .LE. ABS(PZSTRL) ) THEN
        XUNTER=1.
      ELSE
        YSTRL=.5*LOG((ESTRL+PZSTRL)/(ESTRL-PZSTRL))
        XUNTER=ABS(0.94*SINH(YSTRL))/EPART
        IF ( XUNTER .GT. 1. ) XUNTER=1.
      ENDIF
      XUNTER=XUNTER-(XOBER-XUNTER)
      X=XUNTER+X*(XOBER-XUNTER)
C-C   EN=APART+X*(EPART-APART)
      PZ=X*EPART
      IF ( PZ .LT. 0. ) THEN
        SGNPO=-1.
      ELSE
        SGNPO=1.
      ENDIF
      EN=SQRT(PZ**2+AT**2)
C-C   IF ( AT .GT. EN ) AT=APART+RANGEN()*(EN-APART)
C-C   PT=SQRT((AT-APART)*(AT+APART))
      PHI=2.*PI*RANGEN()
      PO(1)=PT*COS(PHI)
      PO(2)=PT*SIN(PHI)
C-C   PO(3)=SGNSIL*SQRT(EN**2-AT**2)
      PO(3)=SGNSIL*PZ
      PO(4)=EN
      PO123=SQRT(PO(3)**2+PO(2)**2+PO(1)**2)
      IF ( PO(4)-PO123 .GT. 0.D0 ) THEN
        PO(5)=SQRT((PO(4)-PO123)*(PO(4)+PO123))
      ELSE
        PO(5)=0.D0
        IF ( PO(4)-PO123 .LT. -1.D-4*PO(4) ) THEN
          IF(ISH.GE.90)THEN
            CALL UTMSG('JSPLIT')
            WRITE(IFCH,*)'*****  !PO123! > PO4'
            WRITE(IFCH,*)'PO123,PO4:',PO123,PO(4)
            WRITE(IFCH,*)'PO1:',PO(1)
            WRITE(IFCH,*)'PO2:',PO(2)
            WRITE(IFCH,*)'PO3:',PO(3)
            WRITE(IFCH,*)'PO5:',PO(5)
            CALL UTMSGF
          ENDIF
        ENDIF
      ENDIF
      POX=PO(1)
      POY=PO(2)
      POT=SQRT(POX**2+POY**2)
      POP=PO(4)+ABS(PO(3))*SGNPO
      POM=PO(4)-ABS(PO(3))*SGNPO
      STRO(1,1)=STRO(1,1)+PO(1)
      STRO(2,1)=STRO(2,1)+PO(2)
      STRO(3,1)=STRO(3,1)-EPART*SGNSIL +PO(3)
      STRO(4,1)=STRO(4,1)-EPART +PO(4)
      STRO(1,KMAX)=0.
      STRO(2,KMAX)=0.
      STRO(3,KMAX)=0.
      STRO(4,KMAX)=0.
      IF ( PO(4)-ABS(PO(3)) .GT. 0.D0 ) THEN
        YLD=0.5*LOG((PO(4)+ABS(PO(3)))/(PO(4)-ABS(PO(3))))
      ELSE
        YLD=100.
      ENDIF

C  DETERMINE SPLIT STRING (MOMENTA)
C  --------------------------------
      PP=PAP+PUP
      PM=PAM+PUM
      PEP=PP-POP
      PEM=PM-POM
      IF ( PEP .LE. 0.D0 ) GOTO 5000
      IF ( PEM .LE. 0.D0 ) GOTO 5000
      A=(PEM*PEP-PUT**2-POT**2)*0.5D0
      D=PUT*POT
      DAUXIL = A**2-D**2
      IF ( DAUXIL .LT. 0.D0 ) GOTO 5000
      DAUXIL=SQRT(DAUXIL)
      PYP=A+PUT**2-DAUXIL
      IF ( PYP .LT. 0.D0  .AND.  PYP .GT. -1.D-6 ) PYP=0.D0
      PYP=PYP/PEM
      PYM=A+POT**2-DAUXIL
      IF ( PYM .LT. 0.D0  .AND.  PYM .GT. -1.D-6 ) PYM=0.D0
      PYM=PYM/PEP
      IF ( PYP .LT. 0.D0 ) GOTO 5000
      IF ( PYM .LT. 0.D0 ) GOTO 5000
      PIP=PEP-PYP
      IF ( PIP .LT. 0.D0 ) GOTO 5000
      PIM=PYM
      PAP=PYP
      PAM=PEM-PYM
      IF ( PAM .LT. 0.D0 ) GOTO 5000
      SSTR(1,1)=-POX
      SSTR(2,1)=-POY
      SSTR(3,1)=SGNSIL*(PIP-PIM)*0.5D0
      SSTR(4,1)=(PIP+PIM)*0.5D0
      SSTR(1,2)=PUX
      SSTR(2,2)=PUY
      SSTR(3,2)=SGNSIL*(PAP-PAM)*0.5D0
      SSTR(4,2)=(PAP+PAM)*0.5D0
      PNL3=SSTR(3,1)+SSTR(3,2)
      PNL4=SSTR(4,1)+SSTR(4,2)
      IF ( PNL4-ABS(PNL3) .NE. 0. ) THEN
        YNL=0.5*LOG((PNL4+ABS(PNL3))/(PNL4-ABS(PNL3)))
      ELSE
        YNL=100.
      ENDIF
      STR(1,1)=SSTR(1,1)
      STR(2,1)=SSTR(2,1)
      STR(3,1)=SSTR(3,1)
      STR(4,1)=SSTR(4,1)
      STR(1,2)=SSTR(1,2)
      STR(2,2)=SSTR(2,2)
      STR(3,2)=SSTR(3,2)
      STR(4,2)=SSTR(4,2)

C  CHECKS
C  ------
      IF ( ISH .GE. 90 ) THEN
        IF ( ABS(PIP*PIM-POT**2) .GT. 1.D-4 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  PIP*PIM /= POT**2'
          WRITE(IFCH,*)'PIP*PIM=',PIP*PIM
          WRITE(IFCH,*)'POT**2=',POT**2
          WRITE(IFCH,*)'PIP=',PIP
          WRITE(IFCH,*)'PIM=',PIM
          WRITE(IFCH,*)'POT=',POT
          CALL UTMSGF
        ENDIF
        IF ( ABS(PAP*PAM-PUT**2) .GT. 1.D-4 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  PAP*PAM /= PUT**2'
          WRITE(IFCH,*)'PAP*PAM=',PAP*PAM
          WRITE(IFCH,*)'PUT**2=',PUT**2
          WRITE(IFCH,*)'PAP=',PAP
          WRITE(IFCH,*)'PAM=',PAM
          WRITE(IFCH,*)'PUT=',PUT
          CALL UTMSGF
        ENDIF
        IF ( ABS(SSTR(4,1)**2
     *     -SSTR(1,1)**2-SSTR(2,1)**2-SSTR(3,1)**2) .GT. 1.D-4 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  MASS**2 OF STRING END 1 NONZERO'
          WRITE(IFCH,*)'MASS**2=',SSTR(4,1)**2
     *                         -SSTR(1,1)**2-SSTR(2,1)**2-SSTR(3,1)**2
          CALL UTMSGF
        ENDIF
        IF ( ABS(SSTR(4,2)**2
     *     -SSTR(1,2)**2-SSTR(2,2)**2-SSTR(3,2)**2) .GT. 1.D-4 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  MASS**2 OF STRING END 2 NONZERO'
          WRITE(IFCH,*)'MASS**2=',SSTR(4,2)**2
     *                         -SSTR(1,2)**2-SSTR(2,2)**2-SSTR(3,2)**2
          CALL UTMSGF
        ENDIF
        DO 14 N=1,4

          IF ( ABS(STR(N,1)+STRO(N,1)-STRO0(N,1)
     *            +STR(N,2)+STRO(N,KMAX)-STRO0(N,KMAX)) .GT. 1.E-4
     *      .AND.  ABS(STR(N,1)+STRO(N,1)-STRO0(N,1)
     *                +STR(N,2)+STRO(N,KMAX)-STRO0(N,KMAX))
     *           .GT. 1.E-4*ABS(STRO0(N,1)+STRO0(N,KMAX)) ) GOTO 15
14      CONTINUE
        GOTO 16
15      CONTINUE
        CALL UTMSG('JSPLIT')
        WRITE(IFCH,*)'*****  P_STR + P_STRO /= P_STRO0'
        WRITE(IFCH,*)'P_STR + P_STRO:'
        WRITE(IFCH,104)((STR(N,1)+STR(N,2)+STRO(N,1)
     *                                    +STRO(N,KMAX)),N=1,4)
        WRITE(IFCH,*)'P_STRO0:'
        WRITE(IFCH,104)((STRO0(N,1)+STRO0(N,KMAX)),N=1,4)
        WRITE(IFCH,*)'STR:'
        WRITE(IFCH,104)(STR(N,1),N=1,4)
        WRITE(IFCH,104)(STR(N,2),N=1,4)
        WRITE(IFCH,*)'STRO:'
        WRITE(IFCH,104)(STRO(N,1),N=1,4)
        WRITE(IFCH,104)(STRO(N,KMAX),N=1,4)
        WRITE(IFCH,*)'STRO0:'
        WRITE(IFCH,104)(STRO0(N,1),N=1,4)
        WRITE(IFCH,104)(STRO0(N,KMAX),N=1,4)
        CALL UTMSGF
16      CONTINUE
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,100)APART,EPART,EN,PT
100       FORMAT(3X,'APART,EPART,EN,PT:',4F13.5)
          WRITE(IFCH,101)(SNGL(PO(K)),K=1,5)
101       FORMAT(3X,'PO:',5F11.5,/)
          IF ( NQU .EQ. 2 ) THEN
            WRITE(IFCH,*)'DIQUARK-ANTIDIQUARK BREAK'
            WRITE(IFCH,*)' '
          ENDIF
        ENDIF
      ENDIF

C  FLAVOUR
C  -------
      IC(1)=NINT(ABS(STRO(4+1,KMAX)))
      IC(2)=NINT(ABS(STRO(4+2,KMAX)))
      STR(4+1,2)=IC(1)
      STR(4+2,2)=IC(2)
      DO 4 N=1,NFLAV
        JCX(N,1)=0
        JCX(N,2)=0
4     CONTINUE
      M=0
      IF     ( NQU .EQ. 1 ) THEN
        IF ( IC(1) .GT. 0 ) M=2
        IF ( IC(2) .GT. 0 ) M=1
        NFL=MOD(IFLTO,10)
        JCX(NFL,M)=JCX(NFL,M)+1
      ELSEIF ( NQU .EQ. 2 ) THEN
        IF ( IC(1) .GT. 0 ) M=1
        IF ( IC(2) .GT. 0 ) M=2
        NFL=MOD(IFLTO,10)
        JCX(NFL,M)=JCX(NFL,M)+1
        NFL=IFLTO/10
        JCX(NFL,M)=JCX(NFL,M)+1
      ENDIF
      IF ( M .EQ. 0 ) THEN
        CALL UTSTOP('JSPLIT: M = 0                           ')
      ENDIF
      CALL IDENCO(JCX,ICX,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        CALL UTSTOP('JSPLIT: IDENCO RET CODE = 1             ')
      ENDIF
      STR(4+1,1)=ICX(1)
      STR(4+2,1)=ICX(2)
      CALL UTAMST(STR,AM,AMIN,IRET)
      IF ( IRET .NE. 0 ) GOTO 1002
      IF ( NFUSN(KMAX) .EQ. 1 ) GOTO 1002
      IC(1)=NINT(STRO(4+1,1))
      IC(2)=NINT(STRO(4+2,1))
      CALL IDDECO(IC,JC)
      DO 26 N=1,NQU
        IF ( N .EQ. 1 ) THEN
          NFL=MOD(IFLTO,10)
        ELSE
          NFL=IFLTO/10
        ENDIF
        IF ( JC(NFL,M) .GT. 0 ) THEN
          JC(NFL,M)=JC(NFL,M)-1
        ELSE
          JC(NFL,3-M)=JC(NFL,3-M)+1
        ENDIF
26    CONTINUE
      NN=0
      DO 27 N=1,NFLAV
        NN=NN+JC(N,1)+JC(N,2)
27    CONTINUE
      IF ( NN .EQ. 0 ) THEN
        NFL=INT(RANGEN()/PUDX)+1
        JC(NFL,1)=1
        JC(NFL,2)=1
      ENDIF
      CALL IDENCO(JC,IC,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        IF ( ISH .GE. 90  .AND.  MESSCO .EQ. 0 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  IDENCO RET CODE = 1.   REDO JSPLIT'
          WRITE(IFCH,*)'JC:'
          WRITE(IFCH,*)JC
          CALL UTMSGF
          MESSCO=1
        ENDIF
        GOTO 5000
      ENDIF
      STRO(4+1,1)=IC(1)
      STRO(4+2,1)=IC(2)
      STRO(4+1,KMAX)=0.
      STRO(4+2,KMAX)=0.

C  OK
C  --
      IER=0
ctp060203 1000  IER=0
      NSCC=NSCC+1
      GOTO 10002

C  ERROR
C  -----
1001  IER=1
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('JSPLIT')
        WRITE(IFCH,*)'*****  SPLIT NOT POSSIBLE'
        CALL UTMSGF
      ENDIF
      GOTO 10001

C  ABSORPTION
C  ----------
1002  IER=2
      NFUSN(KMAX)=1
      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)'ABSORPTION OF STRO(,KMAX)'
        WRITE(IFCH,*)' '
      ENDIF
      DO 18 I=1,NSI
        STR(I,1)=0.
        STR(I,2)=0.
18    CONTINUE
      DO 17 K=1,KMAX
        DO 17 I=1,NSI
          STRO(I,K)=STRO0(I,K)
17    CONTINUE
      STRO(1,1)=STRO(1,1)+STRO(1,KMAX)
      STRO(2,1)=STRO(2,1)+STRO(2,KMAX)
      STRO(3,1)=STRO(3,1)+STRO(3,KMAX)
      STRO(4,1)=STRO(4,1)+STRO(4,KMAX)
      STRO(1,KMAX)=0.
      STRO(2,KMAX)=0.
      STRO(3,KMAX)=0.
      STRO(4,KMAX)=0.
      IC(1)=NINT(STRO(4+1,1))
      IC(2)=NINT(STRO(4+2,1))
      CALL IDDECO(IC,JC)
      ICX(1)=NINT(ABS(STRO(4+1,KMAX)))
      ICX(2)=NINT(ABS(STRO(4+2,KMAX)))
      CALL IDDECO(ICX,JCX)
      DO 22 NF=1,NFLAV
        JC(NF,1)=JC(NF,1)+JCX(NF,1)
        JC(NF,2)=JC(NF,2)+JCX(NF,2)
22    CONTINUE
      CALL IDENCO(JC,IC,IRETEN)
      IF ( IRETEN .EQ. 1 ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('JSPLIT')
          WRITE(IFCH,*)'*****  IDENCO RET CODE = 1'
     *                ,'   (AFTER ABSORPTION)'
          WRITE(IFCH,*)'JC:'
          WRITE(IFCH,*)JC
          CALL UTMSGF
        ENDIF
        GOTO 1001
      ENDIF
      STRO(5,1)=IC(1)
      STRO(6,1)=IC(2)
      STRO(5,KMAX)=0.
      STRO(6,KMAX)=0.
      GOTO 10002

C  CHECK NSCC
C  ----------
10002 CONTINUE
      IF ( KMAXOR.GT.2 .AND. KMAX.EQ.2 .AND. MAX(1,NSCC).NE.NSCCX ) THEN
        IER=3
        IF ( ISH .GE. 91 ) THEN
          WRITE(IFCH,*)'REDO STRING PROCESSING WITH NEW APART, EPART'
          WRITE(IFCH,*)' '
        ENDIF
        GOTO 10001
      ENDIF

C  PRINT
C  -----
      IF ( ISH .GE. 92 ) THEN
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,*)'REMAINDER STRING:'
          WRITE(IFCH,*)' '
          WRITE(IFCH,105)(STRO(I,1),I=1,4),(NINT(STRO(I,1)),I=5,6)
          DO 9 K=2,KMAX
            WRITE(IFCH,104)(STRO(I,K),I=1,4),(NINT(STRO(I,K)),I=5,6)
 9        CONTINUE
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'SPLIT STRING:'
          WRITE(IFCH,*)' '
        ENDIF
        WRITE(IFCH,105)(STR(I,1),I=1,4),(NINT(STR(I,1)),I=5,6)
105     FORMAT(' STR: ',4F13.5,2I8)
        WRITE(IFCH,104)(STR(I,2),I=1,4),(NINT(STR(I,2)),I=5,6)
104     FORMAT('      ',4F13.5,2I8)
        WRITE(IFCH,*)' '
      ENDIF
10001 CONTINUE
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)('-',L=1,25)
        WRITE(IFCH,*)'   RETURN FROM JSPLIT   '
        WRITE(IFCH,*)('-',L=1,25)
        WRITE(IFCH,*)' '
      ENDIF
      ISH=ISH0
      RETURN
      END
C=======================================================================

      SUBROUTINE LEPEXP(RXBJ,RQSQ)

C-----------------------------------------------------------------------
C  GENERATES X_BJORKEN AND Q**2 ACCORDING TO AN EXPERIMENTAL
C  DISTRIBUTION ( GIVEN IN ARRAY XQ(NXBJ,NQSQ) ).
C-----------------------------------------------------------------------
      PARAMETER (NQSQ=10)
      PARAMETER (NXBJ=10)
      REAL  XQ(NXBJ,NQSQ),VXQ(NXBJ*NQSQ)
      EQUIVALENCE (XQ(1,1),VXQ(1))

      DATA VXQ/1304.02,   366.40,    19.84,    10.79,     6.42,
     *            4.54,     4.15,     3.38,     2.03,     1.56,
     *          241.63,  1637.26,   427.36,   164.51,    73.72,
     *           43.07,    20.73,    12.78,     9.34,     5.83,
     *            0.01,   724.66,   563.79,   275.08,   176.13,
     *          106.44,    85.82,    54.52,    37.12,    28.65,
     *            0.01,   202.40,   491.10,   245.13,   157.07,
     *          104.43,    61.05,    49.42,    37.84,    26.79,
     *            0.01,     3.77,   316.38,   226.92,   133.45,
     *           90.30,    63.67,    48.42,    35.73,    28.04,
     *            0.01,     0.01,   153.74,   213.09,   114.14,
     *           76.26,    60.02,    43.15,    43.47,    25.60,
     *            0.01,     0.01,    39.31,   185.74,   108.56,
     *           88.40,    47.29,    39.35,    31.80,    22.91,
     *            0.01,     0.01,     0.01,   104.61,   107.01,
     *           66.24,    45.34,    37.45,    33.44,    23.78,
     *            0.01,     0.01,     0.01,    56.58,    99.39,
     *           67.78,    43.28,    35.98,    34.63,    18.31,
     *            0.01,     0.01,     0.01,    13.56,    76.25,
     *           64.30,    42.80,    28.56,    21.19,    20.75 /
      DATA QSQMIN/4./,QSQWID/4./,XBJMIN/0./,XBJWID/.025/,INIT/0/
      SAVE
C-----------------------------------------------------------------------
      INIT=INIT+1
      IF ( INIT .EQ. 1 ) THEN
        N=NXBJ*NQSQ
        SUM=VXQ(1)
        DO 1 I=2,N
          SUM=SUM+VXQ(I)
          VXQ(I)=VXQ(I)+VXQ(I-1)
 1      CONTINUE
        DO 3 I=1,N
          VXQ(I)=VXQ(I)/SUM
 3      CONTINUE
      ENDIF

      N=NXBJ*NQSQ
      R=RANGEN()
      CALL UTLOC(VXQ,N,R,ILOC)
      IF ( ILOC .GE. N ) ILOC=ILOC-1
      I=MOD(ILOC,NXBJ)+1
      IF ( I .EQ. 0 ) I=NXBJ
      J=ILOC/NXBJ + 1
      IF ( ILOC .GT. 0 ) THEN
        DXINT=VXQ(ILOC+1)-VXQ(ILOC)
      ELSE
        DXINT=VXQ(1)
      ENDIF
      DXBJ=XBJWID*ABS(R-VXQ(ILOC+1))/DXINT
      DY  =QSQWID*RANGEN()
      RXBJ=XBJMIN+XBJWID*FLOAT(I-1)+DXBJ
      RQSQ=QSQMIN+QSQWID*FLOAT(J-1)+DY
      RETURN
      END
C=======================================================================

      SUBROUTINE LEPSTR(NUCLON,XBJ,QSQ,NSTRNG)

C-----------------------------------------------------------------------
C  RETURNS STRING CODE NSTRNG FOR GIVEN NUCLON, XBJ, QSQ.
C  NUCLON  : THE STRUCK NUCLEON (1120/1220 = PROTON/NEUTRON)
C  XBJ,QSQ : X-BJORKEN AND Q**2
C  NSTRNG  : STRING CODE:
C        PROTON STRINGS:          NEUTRON STRINGS:
C        1 : U  ---  UD           11 : U  ---  DD
C        2 : D  ---  UU           12 : D  ---  UD
C        3 : U  ---  UB(UUD)      13 : U  ---  UB(UDD)
C        4 : D  ---  DB(UUD)      14 : D  ---  DB(UDD)
C        5 : S  ---  SB(UUD)      15 : S  ---  SB(UDD)
C        6 : UB ---  U(UUD)       16 : UB ---  U(UDD)
C        7 : DB ---  D(UUD)       17 : DB ---  D(UDD)
C        8 : SB ---  S(UUD)       18 : SB ---  S(UDD)
C-----------------------------------------------------------------------
      REAL QUARKS(9)
      SAVE
C-----------------------------------------------------------------------
      NSTRNG=0

C  PROTON-STRING (VALENCE PART)
      IF     ( NUCLON .EQ. 1120 ) THEN
        QUARKS(1) = 4.* STXU(XBJ,QSQ)
        QUARKS(2) =     STXD(XBJ,QSQ)

C  NEUTRON-STRING (VALENCE PART)
      ELSEIF ( NUCLON .EQ. 1220 ) THEN
        QUARKS(1) = 4.* STXD(XBJ,QSQ)
        QUARKS(2) =     STXU(XBJ,QSQ)

      ELSE
        RETURN
      ENDIF

C  THE SEA CONTRIBUTIONS (PROTON/NEUTRON)
      UDSEA     =     STXUS(XBJ,QSQ)
      SSEA      =     STXS(XBJ,QSQ)
      QUARKS(3) = 4.* UDSEA
      QUARKS(4) =     UDSEA
      QUARKS(5) =     SSEA
      QUARKS(6) = 4.* UDSEA
      QUARKS(7) =     UDSEA
      QUARKS(8) =     SSEA

      QUARKS(9) = 0.
      DO 11 I=1,8
        QUARKS(9) = QUARKS(9)+QUARKS(I)
11    CONTINUE

      R = RANGEN() * QUARKS(9)
      SUMQ = 0.
      DO 12 I=1,8
        NSTRNG = I
        SUMQ = SUMQ + QUARKS(I)
        IF ( R .LE. SUMQ ) GO TO 13
12    CONTINUE
13    CONTINUE

      IF ( NUCLON .EQ. 1220 ) NSTRNG=NSTRNG+10

      RETURN
      END
C=======================================================================

      SUBROUTINE LEPTAR(XBJ,QSQ,MATARG,LATARG,NUCLON)

C-----------------------------------------------------------------------
C  RETURNS NUCLON = ID OF HIT TARGET NUCLEON.
C  XBJ   : X BJORKEN
C  QSQ   : Q SQUARED
C  MATARG : A  OF TARGET
C  LATARG : Z  OF TARGET
C  NUCLON: ID OF TARGET NUCLEON (1120/1220 FOR PROTON/NEUTRON)
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
C  PROTON PART OF THE STRUCTURE FUNCTION:
      F2P = LATARG * STXZPR(XBJ,QSQ)
C  NEUTRON PART OF THE STRUCTURE FUNCTION:
      F2N = (MATARG-LATARG) * STXZNE(XBJ,QSQ)
C  STRUCTURE FUNCTION OF NUCLEUS:
      F2A = F2P + F2N
C  SELECT  THE TARGET-NUCLEON  ( PROTON OR NEUTRON ):
      RN = RANGEN()
      F2RNDM = F2A * RN
      IF ( F2RNDM .LT. F2N ) THEN
        NUCLON = 1220
      ELSE
        NUCLON = 1120
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE NUCINI(OPT,ANUC,LA,MA,ISI)

C-----------------------------------------------------------------------
C  INITIALIZES NUCLEON-MOMENTA.
C  WRITES NUCLEONS ON /CPTL/ (P,IFR,ICL).
C-----------------------------------------------------------------------
      PARAMETER (MAMX=56)
      PARAMETER (MAMX2=MAMX*2)
      PARAMETER (MXPTL=70000)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /CNNN/    NNNPTL(MAMX2)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL      ANUC(NSI,NHA,MAMX)
      CHARACTER OPT*3
      SAVE
C-----------------------------------------------------------------------
      NHAX=3
      IF ( NHA .LT. NHAX ) THEN
        CALL UTSTOP('NUCINI: NHA TOO SMALL                   ')
      ENDIF
      LAS=0
      MAS=0

      DO 1 L=1,MA
        DO 2 M=1,NHA
          DO 3 N=1,NSI
            ANUC(N,M,L)=0.
 3        CONTINUE
 2      CONTINUE

        IF     ( OPT .EQ. 'NUC' ) THEN
          IF     ( ISI .GT. 0  .AND.  LAPROJ .LT. 0 ) THEN
            ID=IDPROJ
          ELSEIF ( ISI .LT. 0  .AND.  LATARG .LT. 0 ) THEN
            ID=IDTARG
          ELSE
            IF ( RANGEN() .LE. (LA-LAS)/FLOAT(MA-MAS)) THEN
              ID=1120
              LAS=LAS+1
            ELSE
              ID=1220
            ENDIF
            MAS=MAS+1
          ENDIF
          CALL IDMASS(ID,AMS)
          IC1=IDTRAI(1,ID,1)
          IC2=IDTRAI(2,ID,1)
          P1=0.
          P2=0.
          P3=ISI*PNLLX
          P0=SQRT(PNLLX**2+AMS**2)
          P0X=PNLLX
          P5=AMS

        ELSEIF ( OPT .EQ. 'STR' ) THEN
          IF ( RANGEN() .LE. (LA-LAS)/FLOAT(MA-MAS) ) THEN
            ID=1120
            LAS=LAS+1
          ELSE
            ID=1220
          ENDIF
          MAS=MAS+1
          IC1=IDTRAI(1,ID,1)
          IC2=IDTRAI(2,ID,1)
          P1=0.
          P2=0.
          P3=0.
          P0=PROM
          P5=PROM
        ENDIF

        ANUC(5,1,L)=IC1
        ANUC(6,1,L)=IC2
        ANUC(3,2,L)=P3
        ANUC(4,2,L)=P0X
        ANUC(5,2,L)=IC1
        ANUC(6,2,L)=IC2
        NPTL=NPTL+1
        IF ( NPTL .NE. NNNPTL(NPTL) ) THEN
          CALL UTSTOP('NUCINI: NPTL AND NNNPTL DONT MATCH      ')
        ENDIF
        IDPTL(NPTL)=ID
        PPTL(1,NPTL)=P1
        PPTL(2,NPTL)=P2
        PPTL(3,NPTL)=P3
        PPTL(4,NPTL)=P0
        PPTL(5,NPTL)=P5
        IFRPTL(1,NPTL)=0
        IFRPTL(2,NPTL)=0
        ICLPTL(NPTL)=1
        NQJPTL(NPTL)=0
 1    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE NUCLCO(MASSNR,N,X,Y,Z,YNUC)

C-----------------------------------------------------------------------
C  CALCULATES COORDINATES OF THE NUCLEONS IN A NUCLEUS.
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MXPTL=70000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL X(N),Y(N),Z(N)
      SAVE
C-----------------------------------------------------------------------
      IF ( MASSNR .EQ. 0 ) RETURN
      IF ( MASSNR .GT. N ) THEN
        CALL UTSTOP('NUCLCO: MASSNR.GT.N                     ')
      ENDIF
      IF ( MASSNR .EQ. 1 ) THEN
        X(1)=0.
        Y(1)=0.
        Z(1)=0.
        RETURN
      ENDIF
      DO 3 I=1,MASSNR
10      CONTINUE
        IF ( MASSNR .EQ. MAPROJ )
     *                      R=UTINVT(NDEP,XDEP,QDEP,RANGEN()*QDEP(NDEP))
        IF ( MASSNR .EQ. MATARG ) THEN
          IF     ( LTARG .EQ. 1 ) THEN
            R=UTINVT(NDET,XDET14,QDET14,RANGEN()*QDET14(NDET))
          ELSEIF ( LTARG .EQ. 2 ) THEN
            R=UTINVT(NDET,XDET16,QDET16,RANGEN()*QDET16(NDET))
          ELSEIF ( LTARG .EQ. 3 ) THEN
            R=UTINVT(NDET,XDET40,QDET40,RANGEN()*QDET40(NDET))
          ELSE
            R=UTINVT(NDET,XDET99,QDET99,RANGEN()*QDET99(NDET))
          ENDIF
        ENDIF
        IF ( MASSNR.NE.MAPROJ .AND. MASSNR.NE.MATARG ) THEN
          CALL UTSTOP('NUCLCO: NUCLEUS NEITHER PROJ NOR TARG   ')
        ENDIF
        COSTHE=1.-2.*RANGEN()
        SINTHE= SQRT(1. - COSTHE**2)
        PHI=2.*PI*RANGEN()
        X(I)=R*SINTHE*COS(PHI)
        Y(I)=R*SINTHE*SIN(PHI)
        Z(I)=R*COSTHE
        IF ( I .EQ. 1 ) GOTO 3
        IF ( CORE .EQ. 0. ) GOTO 3
        DO 2 J=1,I-1
          IF ( (X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2
     *                               .LT. CORE**2 ) GOTO 10
2       CONTINUE
3     CONTINUE
      IF ( ISH .GE. 93 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'NUCLEON COORDINATES:'
      ENDIF
      AUXIL = 1./COSH(YNUC)
      DO 4 I=1,MASSNR
        Z(I)=Z(I)*AUXIL
4     CONTINUE
      IF ( ISH .GE. 93 ) THEN
        DO 5 I=1,MASSNR
          WRITE(IFCH,*)'I X Y Z: ',I,X(I),Y(I),Z(I)
5       CONTINUE
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE NUCOGE

C-----------------------------------------------------------------------
C  WRITES NUCLEONS ON /CPTL/ (XOR,TIV,IST,IOR,JOR,NST).
C  DETERMINES SEQUENCE OF COLLISIONS ACCORDING TO NUCLEAR GEOMETRY
C    IF MAPROJ>0.
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MAMX2=MAMX*2)
      PARAMETER (MXPTL=70000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CNCL/    XPROJ(MAMX),XTARG(MAMX),YPROJ(MAMX)
     *                ,YTARG(MAMX),ZPROJ(MAMX),ZTARG(MAMX)
      COMMON /CNNN/    NNNPTL(MAMX2)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      IF ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'NUCOGE (ENTRY)'
      ENDIF
      VELI=1./TANH(YPJTL-YHAHA)+TANH(YHAHA)
      NPROJ=1
      NTARG=1
      DO 15 N=1,KOLLMX
        NORD(N)=N
        COORD(1,N)=0.
        COORD(2,N)=0.
        COORD(3,N)=0.
        COORD(4,N)=0.
15    CONTINUE

      IF     ( MATARG .LT. 0 ) THEN
        KOLL=-MATARG
        NTARG=KOLL
        BIMP=0.
        NPTL=NPTL+1
        XORPTL(1,NPTL)=0.
        XORPTL(2,NPTL)=0.
        XORPTL(3,NPTL)=0.
        XORPTL(4,NPTL)=0.
        TIVPTL(1,NPTL)=-AINFIN
        TIVPTL(2,NPTL)=0.
        ISTPTL(NPTL)=1
        IORPTL(NPTL)=-1
        JORPTL(NPTL)=0
        NNNPTL(NPTL)=NPTL
        DO 1 K=1,KOLL
          NRPROJ(K)=1
          NRTARG(K)=K
          NPTL=NPTL+1
          XORPTL(1,NPTL)=0.
          XORPTL(2,NPTL)=0.
          XORPTL(3,NPTL)=0.
          XORPTL(4,NPTL)=0.
          TIVPTL(1,NPTL)=-AINFIN
          TIVPTL(2,NPTL)=0.
          ISTPTL(NPTL)=1
          IORPTL(NPTL)=-1
          JORPTL(NPTL)=0
          NNNPTL(NPTL)=NPTL
1       CONTINUE
        GOTO 9999

      ELSEIF ( MAPROJ .EQ. 1  .AND.  MATARG .EQ. 1 ) THEN
        KOLL=1
        BIMP=0.
        NRPROJ(1)=1
        NRTARG(1)=1
        DO 5 II=1,2
          NPTL=NPTL+1
          XORPTL(1,NPTL)=0.
          XORPTL(2,NPTL)=0.
          XORPTL(3,NPTL)=0.
          XORPTL(4,NPTL)=0.
          TIVPTL(1,NPTL)=-AINFIN
          TIVPTL(2,NPTL)= 0.
          ISTPTL(NPTL)  = 1
          IORPTL(NPTL)  =-1
          JORPTL(NPTL)  = 0
          NNNPTL(NPTL)  = NPTL
 5      CONTINUE
        GOTO 9999
      ENDIF

      CALL NUCLCO(MAPROJ,MAMX,XPROJ,YPROJ,ZPROJ,YPJTL-YHAHA)
      CALL NUCLCO(MATARG,MAMX,XTARG,YTARG,ZTARG,YHAHA)
      BX=0.
      BY=0.
      IF ( MAPROJ .GT. 0 ) THEN
        IF ( BIMEVT .LT. 0. ) THEN
          B1=BMINIM
          B2=MIN(RMPROJ+RMTARG(LTARG),BMAXIM)
          IF ( B1 .GT. B2 ) THEN
            CALL UTSTOP('NUCOGE: BMIN > BMAX                     ')
          ENDIF
          BIMP=SQRT(B1**2+(B2**2-B1**2)*RANGEN())
          PHI=2.*PI*RANGEN()
        ELSE
          PHI=PHIEVT
          BIMP=BIMEVT
        ENDIF
        BX=COS(PHI)*BIMP
        BY=SIN(PHI)*BIMP
        DO 6 I=1,MAPROJ
          NPTL=NPTL+1
          XORPTL(1,NPTL)=XPROJ(I)+BX*0.5
          XORPTL(2,NPTL)=YPROJ(I)+BY*0.5
          XORPTL(3,NPTL)=ZPROJ(I)
          XORPTL(4,NPTL)=0.
          TIVPTL(1,NPTL)=-AINFIN
          TIVPTL(2,NPTL)= AINFIN
          ISTPTL(NPTL)=0
          IORPTL(NPTL)=0
          JORPTL(NPTL)=0
          NNNPTL(NPTL)=NPTL
 6      CONTINUE
      ENDIF
      DO 7 I=1,MATARG
        NPTL=NPTL+1
        XORPTL(1,NPTL)=XTARG(I)-BX*0.5
        XORPTL(2,NPTL)=YTARG(I)-BY*0.5
        XORPTL(3,NPTL)=ZTARG(I)
        XORPTL(4,NPTL)=0.
        TIVPTL(1,NPTL)=-AINFIN
        TIVPTL(2,NPTL)= AINFIN
        ISTPTL(NPTL)=0
        IORPTL(NPTL)=0
        JORPTL(NPTL)=0
        NNNPTL(NPTL)=NPTL
 7    CONTINUE
      IF ( MAPROJ .EQ. 0 ) GOTO 1000
      KOLL=0
      DO 12 I=1,MAPROJ
        DO 12 J=1,MATARG
          IF ( PI*( (XPROJ(I)+BX-XTARG(J))**2
     *        +(YPROJ(I)+BY-YTARG(J))**2 ) .GT. SIGPPI ) GOTO 12
          IF ( KOLL .GT. 0 ) THEN
            KP=0
            KT=0
            DO 30 KL=1,KOLL
              IF ( NRPROJ(KL) .EQ. I ) KP=1
              IF ( NRTARG(KL) .EQ. J ) KT=1
30          CONTINUE
CDH         FTRSIG=1.0
            FTR=1.0
CDH         IF ( KP .EQ. 1 ) FTR=FTR*FTRSIG
CDH         IF ( KT .EQ. 1 ) FTR=FTR*FTRSIG
            IF ( KP.EQ.1 .OR. KT.EQ.1 ) THEN
              IF ( PI*( (XPROJ(I)+BX-XTARG(J))**2
     *           +(YPROJ(I)+BY-YTARG(J))**2 ) .GT. FTR*SIGPPI ) GOTO 12
            ENDIF
          ENDIF
          KOLL=KOLL+1
          IF ( KOLL .GT. KOLLMX ) THEN
            CALL UTSTOP('NUCOGE: KOLLMX TOO SMALL                ')
          ENDIF
          NRPROJ(KOLL)=I
          NRTARG(KOLL)=J
          DISTCE(KOLL)=ZTARG(J)-ZPROJ(I)
          COORD(1,KOLL)=(XPROJ(I)+XTARG(J))*0.5
          COORD(2,KOLL)=(YPROJ(I)+YTARG(J))*0.5
          COORD(3,KOLL)=(ZPROJ(I)+ZTARG(J))*0.5
          COORD(4,KOLL)=DISTCE(KOLL)*VELI
          ISTPTL(I)=1
          IORPTL(I)=-1
          TIVPTL(2,I)=COORD(4,KOLL)
          ISTPTL(MAPROJ+J)=1
          IORPTL(MAPROJ+J)=-1
          TIVPTL(2,MAPROJ+J)=COORD(4,KOLL)
12    CONTINUE
      IF ( KOLL .LE. 1 ) GOTO 9999

      DO 21 N=2,KOLL
        DO 22 M=1,N-1
          IF ( NRPROJ(M) .EQ. NRPROJ(N) ) GOTO 21
22      CONTINUE
        NPROJ=NPROJ+1
21    CONTINUE
      DO 23 N=2,KOLL
        DO 24 M=1,N-1
          IF ( NRTARG(M) .EQ. NRTARG(N) ) GOTO 23
24      CONTINUE
        NTARG=NTARG+1
23    CONTINUE

      DO 20 N=1,KOLL-1
        DO 20 M=N+1,KOLL
          IF ( DISTCE(NORD(M)) .LT. DISTCE(NORD(N)) ) THEN
            NORDM=NORD(M)
            NORD(M)=NORD(N)
            NORD(N)=NORDM
          ENDIF
20    CONTINUE

9999  CONTINUE
      IF ( KOLL .LE. 0 ) GOTO 1000
      IF ( KOLL .LT. KO1KO2/10000  .OR.  KOLL .GT. MOD(KO1KO2,10000) )
     *                                                     GOTO 1000
      NEVT=1
      BIMEVT=BIMP
      PHIEVT=PHI
      KOLEVT=KOLL
      NPJEVT=NPROJ
      NTGEVT=NTARG
      PMXEVT=PNLL
      EGYEVT=ENGY

1000  CONTINUE
      IF ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)'NUCOGE (EXIT)'
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE NUCOLL

C-----------------------------------------------------------------------
C  PERFORMS A  NUCLEUS-NUCLEUS COLLISION (INCL. NUCLEON-NUCLEON)
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NPRBMS=20)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CIV/     IVP,IVP0,IVT,IVT0
      COMMON /CKOL/    KOL
      COMMON /CNCE/    NCES,NCOLEX
      COMMON /CNEW/    KOTRI,NEWCOL,NEWICO
      COMMON /CNFR/    NRFRA
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CNTEVM/  NTEVM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPRBMS/  PRBMS(NPRBMS)
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARTNR/  PEX,PEY,PEZET,PE0,PX4,PY4,SUMMAS
     *                ,IC4,IPTNR,JS4,NPS

      REAL      PROJ(NSI,NHA),PROJA(NSI,NHA,MAMX)
     *         ,TARG(NSI,NHA),TARGA(NSI,NHA,MAMX)
      INTEGER   IC4(2)
      CHARACTER DASH*1
      DATA DASH/'-'/
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 6 ) ISH=MOD(ISHSUB,100)
      IF     ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'NUCOLL (ENTRY)'
CDH   ELSEIF ( ISH .EQ. 14 ) THEN
CDH     CALL UTTIMA('*** NUCOLL *** ')
      ENDIF

      IF ( ISHSUB/100 .EQ. 7 ) ISH=MOD(ISHSUB,100)
      NEVT=0
      NPTL=0
      CALL NUCOGE
      IF ( ICHOIC .EQ. 3  .OR.  KOLL .EQ. 0 ) GOTO 1000
      NAEVT=NAEVT+1
      IF ( KOLL .LT. KO1KO2/10000  .OR.
     *     KOLL .GT. MOD(KO1KO2,10000) ) GOTO 1000
      COLEVT=KOLL
      AMSEV=SQRT((NPJEVT*SQRT(AMPROJ**2+PNLLX**2)
     *           +NTGEVT*SQRT(AMTARG**2+PNLLX**2))**2
     *                -((NPJEVT-NTGEVT)*PNLLX)**2)
      IF ( ISHSUB/100 .EQ. 7 ) ISH=ISH0

      ITRY=0
38    CONTINUE
      ITRY=ITRY+1
      NSTR=0
      NPTL=0
      SUMPT2=0.
      AMSAC=0.
      CALL NUCINI('NUC',PROJA,LAPROJ,MAPROJ,1)
      CALL NUCINI('NUC',TARGA,LATARG,IABS(MATARG),-1)

      KOLRAN=RANGEN()*KOLL+1
      KOLRAN=MIN(KOLRAN,KOLL)
      KOLRAN=MAX(KOLRAN,1)

      DO 2 KOLS=1,KOLL
        KOL=KOLS
        ISKIP=0
        I=NRPROJ(NORD(KOL))
        J=NRTARG(NORD(KOL))
        IPROJ=I
        ITARG=J
        NRECOP=0
        NRECOT=0
        IF ( KOL .LT. KOLL ) THEN
          DO 33 K=KOL+1,KOLL
            IF ( NRPROJ(NORD(K)) .EQ. I ) NRECOP=NRECOP+1
            IF ( NRTARG(NORD(K)) .EQ. J ) NRECOT=NRECOT+1
33        CONTINUE
        ENDIF
        NCOP=0
        NCOT=0
        IF ( KOL .GT. 1 ) THEN
          DO 34 K=1,KOL-1
            IF ( NRPROJ(NORD(K)) .EQ. I ) NCOP=NCOP+1
            IF ( NRTARG(NORD(K)) .EQ. J ) NCOT=NCOT+1
34        CONTINUE
        ENDIF
        IF ( NCOP.GE.NCOLMX  .OR.  NCOT.GE.NCOLMX ) ISKIP=1

        DO 56 M=1,NHA
          SM=0.
          DO 57 N=1,NSI
            PROJ(N,M)=PROJA(N,M,I)
            SM=SM+PROJ(N,M)**2
57        CONTINUE
          IF ( M .GE. 3  .AND.  SM .LT. 1.E-5 ) GOTO 58
56      CONTINUE
58      CONTINUE
        DO 53 M=1,NHA
          SM=0.
          DO 54 N=1,NSI
            TARG(N,M)=TARGA(N,M,J)
            SM=SM+TARG(N,M)**2
54        CONTINUE
          IF ( M .GE. 3  .AND.  SM .LT. 1.E-5 ) GOTO 55
53      CONTINUE
55      CONTINUE

        IF ( KOL .EQ. KOLRAN  .AND.  JPSI .EQ. 1 ) THEN
          CALL PVJPSF(PROJ,TARG)
          ISKIP=1
        ENDIF

C  NR OF COLOUR EXCHANGES
C  ----------------------

        AMSAC0=AMSAC
        ISKIP0=ISKIP
        SMPT2=SUMPT2
        CALL UTREMB(PROJ,TARG,1)
        GOTO 4
 3      CONTINUE
        IF ( ISH .GE. 91 ) WRITE(IFCH,*)'REDO HH COLLISION'
        CALL UTREST(PROJ,TARG,1)
        SUMPT2=SMPT2
        ISKIP=ISKIP0
        AMSAC=AMSAC0
 4      CONTINUE

        NCOLEX=1
        IF ( ISKIP .NE. 1 ) THEN
          LO=0
16        LO=LO+1
          IF ( LO .EQ. 3 ) THEN
            IF ( ISH .GE. 90 ) THEN
              CALL UTMSG('NUCOLL')
              WRITE(IFCH,*)'*****  LO=3'
              CALL UTMSGF
            ENDIF
          ENDIF
          R=RANGEN()
          NCOLEX=0
15        NCOLEX=NCOLEX+1
          IF ( NCOLEX .GT. NPRBMS ) GOTO 16
          IF ( R .GT. PRBMS(NCOLEX) ) GOTO 15
        ENDIF

C  HADRON-HADRON COLLISION
C  -----------------------

        DO 31 NCE=1,NCOLEX
          NCES=NCE
          CALL UTPAGE
          IF ( ISH .GE. 91 ) THEN
            WRITE(IFCH,101)(DASH,L=1,79),IPAGE,KOL,NCE
     *        ,NRPROJ(NORD(KOL)),NRTARG(NORD(KOL)),(DASH,L=1,79)
101         FORMAT(/,1X,79A1,/,1X,I5,'.PAGE             COL: ',I2
     *           ,'   CEX: ',I2,'   PRJ: ',I3,'   TRG: ',I3,/,1X,79A1,/)
          ENDIF
          IF ( NRECOP .EQ. 0 ) THEN
            ISTORP=1
          ELSE
            ISTORP=0
          ENDIF
          IF ( NRECOT .EQ. 0 ) THEN
            ISTORT=1
          ELSE
            ISTORT=0
          ENDIF
          CALL HAHABS(PROJ,TARG
     *           ,ISTORP*(NCE/NCOLEX),ISTORT*(NCE/NCOLEX),ISKIP,IRETHH)
          IF ( ISKIP .GE. 2 ) GOTO 9997
          IF ( IRETHH .EQ. 1 ) GOTO 3
          IF ( IRESCL .EQ. 1  .AND.  AMSAC .GT. AMSEV ) GOTO 9998
          SUMPT2=SUMPT2+
     *            PROJ(1,2)**2+PROJ(2,2)**2+TARG(1,2)**2+TARG(2,2)**2
          IF ( ISH .EQ. 11 ) WRITE(IFCH,*)'SUMPT2:',NREVT,KOL,NCE,SUMPT2
          ISKIP=0
31      CONTINUE

        DO 40 M=1,NHA
          SM=0.
          DO 41 N=1,NSI
            PROJA(N,M,I)=PROJ(N,M)
            SM=SM+PROJ(N,M)**2
41        CONTINUE
          IF ( M .GE. 3  .AND.  SM .LT. 1.E-5 ) GOTO 42
40      CONTINUE
42      CONTINUE
        DO 43 M=1,NHA
          SM=0.
          DO 44 N=1,NSI
            TARGA(N,M,J)=TARG(N,M)
            SM=SM+TARG(N,M)**2
44        CONTINUE
          IF ( M .GE. 3  .AND.  SM .LT. 1.E-5 ) GOTO 45
43      CONTINUE
45      CONTINUE
  2   CONTINUE

      IF ( SUMPT2 .LT. 1.E-5 ) GOTO 9999

1000  CONTINUE
CDH   IF ( ISH .EQ. 14 ) CALL UTTIMA('    NUCOLL F   ')
      IF ( ISH .EQ. 17  .OR.  ISH .GT. 92 ) THEN
        WRITE(IFCH,*)'NUCOLL (EXIT)'
      ENDIF
      ISH=ISH0
      RETURN

9999  INOIAC=INOIAC+1
      IF ( ISH .GE. 91  .OR.  ISH .EQ. 11 ) THEN
        CALL UTMSG('NUCOLL')
        WRITE(IFCH,*)'*****  NO INTERACTION. REDO NUCOLL'
        CALL UTMSGF
      ENDIF
      GOTO 38

9998  ILAMAS=ILAMAS+1
      IF ( ISH .GE. 91  .OR. ISH .EQ. 12 ) THEN
        CALL UTMSG('NUCOLL')
        WRITE(IFCH,*)'*****  AMSAC>AMSEV: ',AMSAC,AMSEV
     *    ,' . REDO NUCOLL'
        DO 1 NS=1,NSTR
          WRITE(IFCH,109)NS,(ICSTR(L,NS)/100,L=1,4)
     *              ,(PSTR(L,NS),L=3,5)
109       FORMAT(1X,I3,3X,4I5,3(E11.3))
1       CONTINUE
        CALL UTMSGF
      ENDIF
      GOTO 38

9997  CONTINUE
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('NUCOLL')
        WRITE(IFCH,*)'*****  ISKIP>=2. REDO NUCOLL'
        CALL UTMSGF
      ENDIF
      GOTO 38

      END
C=======================================================================

      SUBROUTINE NUCSTR(IER)

C-----------------------------------------------------------------------
C  PERFORMES X AND P TRAFOS FOR NUCLEONS FOR STRING DECAY IN NUCLEUS
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CLEP/    ICINPU,IDSCAT
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      IF ( MATARG .NE. NPTL ) THEN
        CALL UTSTOP('NUCSTR: MATARG /= NPTL                  ')
      ENDIF
      NCNT=0
3     NCNT=NCNT+1
      IF ( NCNT .GE. 10 ) GOTO 1001
      N0=MIN(1+INT(RANGEN()*NPTL),NPTL)
      IF ( ICINPU .GT. 0  .AND.  IDPTL(N0) .NE. IDSCAT ) GOTO 3
      ISTPTL(N0)=1
      IORPTL(N0)=-1
      TIVPTL(2,N0)=0.
      DO 2 N=1,NPTL
        XORPTL(1,N)=XORPTL(1,N)-XORPTL(1,N0)
        XORPTL(2,N)=XORPTL(2,N)-XORPTL(2,N0)
        XORPTL(3,N)=XORPTL(3,N)-XORPTL(3,N0)
        PHI=2.*PI*RANGEN()
        P1=      -ELEPTO*SIN(ANGMUE)*SIN(PHI)
        P2=      -ELEPTO*SIN(ANGMUE)*COS(PHI)
        P3=ELEPTI-ELEPTO*COS(ANGMUE)
        CALL UTROTA(1,P1,P2,P3,XORPTL(1,N),XORPTL(2,N),XORPTL(3,N))
        P3=SQRT(P1**2+P2**2+P3**2)
        P1=0.
        P2=0.
        P4=PROM+ELEPTI-ELEPTO
        CALL UTLOBO(1,P1,P2,P3,P4,SQRT(P4**2-P3**2-P2**2-P1**2)
     *        ,XORPTL(1,N),XORPTL(2,N),XORPTL(3,N),XORPTL(4,N))
        CALL UTLOBO(1,P1,P2,P3,P4,SQRT(P4**2-P3**2-P2**2-P1**2)
     *         ,PPTL(1,N),PPTL(2,N),PPTL(3,N),PPTL(4,N))
2     CONTINUE
      IER=0
      RETURN

1001  IER=1
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('NUCSTR')
        WRITE(IFCH,*)'*****  IDSCAT NOT POSSIBLE ==> REDO EVENT.'
        CALL UTMSGF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE PVJPSF(PROJ,TARG)

C-----------------------------------------------------------------------
C  FORMS A JPSI
C-----------------------------------------------------------------------
      PARAMETER (KOLLMX=2500)
      PARAMETER (MAMX=56)
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NDEP=129)
      PARAMETER (NDET=129)
      PARAMETER (NGAU=129)
      PARAMETER (NPTJ=129)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CGAU/    QGAU(NGAU),XGAU(NGAU)
      COMMON /CKOL/    KOL
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /COL/     BIMP,BMAX,COORD(4,KOLLMX),DISTCE(KOLLMX)
     *                ,QDEP(NDEP),QDET14(NDET),QDET16(NDET),QDET40(NDET)
     *                ,QDET99(NDET),RMPROJ,RMTARG(4),XDEP(NDEP)
     *                ,XDET14(NDET),XDET16(NDET),XDET40(NDET)
     *                ,XDET99(NDET)
     *                ,KOLL,LTARG,NORD(KOLLMX),NPROJ,NRPROJ(KOLLMX)
     *                ,NRTARG(KOLLMX),NTARG
      COMMON /CPROJA/  IPROJ,ITARG,KPROJA(NHA,MAMX),KTARGA(NHA,MAMX)
      COMMON /CPTJ/    QPTJ(NPTJ),XPTJ(NPTJ)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL PROJ(NSI,NHA),PR(4),TARG(NSI,NHA),TG(4)
      SAVE
C-----------------------------------------------------------------------
      PAX=0.19
      PROX=PROJ(3,2)
      TARX=TARG(3,2)
      LOOP=0
      PR(1)=PROJ(1,2)
      PR(2)=PROJ(2,2)
      PR(3)=PROJ(3,2)
      PR(4)=PROJ(4,2)
      TG(1)=TARG(1,2)
      TG(2)=TARG(2,2)
      TG(3)=TARG(3,2)
      TG(4)=TARG(4,2)

      CALL UTPAGE
      IF ( ISH .GE. 91 ) WRITE(IFCH,110)('-',L=1,79),IPAGE,('-',L=1,79)
110   FORMAT(1X,79A1,/,1X,I5,'.PAGE            '
     *         ,'JPSI FORMATION',/,1X,79A1,/)

5000  LOOP=LOOP+1
      IF ( LOOP .GT. 100 ) THEN
        IF(ISH.GE.90)THEN
          CALL UTMSG('PVJPSF')
          WRITE(IFCH,*)'*****  JPSI FORMATION NOT POSSIBLE'
          CALL UTMSGF
        ENDIF
        GOTO 1000
      ENDIF
      PROJ(1,2)=PR(1)
      PROJ(2,2)=PR(2)
      PROJ(3,2)=PR(3)
      PROJ(4,2)=PR(4)
      TARG(1,2)=TG(1)
      TARG(2,2)=TG(2)
      TARG(3,2)=TG(3)
      TARG(4,2)=TG(4)

C  JPSI MOMENTA
C  ------------
      ID=441
      CALL IDMASS(ID,AM)
      S=AM**2
      PT=UTINVT(NPTJ,XPTJ,QPTJ,RANGEN()*QPTJ(NPTJ))
      PHI=2.*PI*RANGEN()
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      LO=0
 1    LO=LO+1
      IF ( LO .GT. 10 ) THEN
        CALL UTSTOP('PVJPSF: LO > 10                         ')
      ENDIF
      Z=PAX*UTINVT(NGAU,XGAU,QGAU,RANGEN()*QGAU(NGAU))
      IF ( Z .GT. 1. ) GOTO 1
      IF ( RANGEN() .LE. 0.5 ) THEN
        PZ=Z*PNLLX
      ELSE
        PZ=-Z*PNLLX
      ENDIF
      E=SQRT(S+PX**2+PY**2+PZ**2)
      PP=E+PZ
      PM=E-PZ

C  PROJ AND TARG MOMENTA
C  ---------------------
      R=RANGEN()
      POX=PROJ(1,2)-R*PX
      POY=PROJ(2,2)-R*PY
      POT2=(POX**2+POY**2)
      R=1.-R
      PUX=TARG(1,2)-R*PX
      PUY=TARG(2,2)-R*PY
      PUT2=(PUX**2+PUY**2)
      POP=PROJ(4,2)+PROJ(3,2)
      POM=PROJ(4,2)-PROJ(3,2)
      PUP=TARG(4,2)+TARG(3,2)
      PUM=TARG(4,2)-TARG(3,2)
      PEP=POP+PUP-PP
      PEM=POM+PUM-PM
      A=(PEM*PEP-PUT2-POT2)*0.5
      D2=PUT2*POT2
      AUXIL = A**2 - D2
      IF ( AUXIL .LT. 0. ) GOTO 5000
      AUXIL = SQRT(AUXIL)
      PYP=(A+PUT2-AUXIL)/PEM
      PYM=(A+POT2-AUXIL)/PEP
      IF ( PYP .LT. 0. ) GOTO 5000
      IF ( PYM. LT. 0. ) GOTO 5000
      PIP=PEP-PYP
      IF ( PIP .LT. 0. ) GOTO 5000
      PIM=PYM
      PAP=PYP
      PAM=PEM-PYM
      IF ( PAM .LT. 0. ) GOTO 5000
      PROJ(1,2)=POX
      PROJ(2,2)=POY
      PROJ(3,2)=(PIP-PIM)*0.5
      PROJ(4,2)=(PIP+PIM)*0.5
      TARG(1,2)=PUX
      TARG(2,2)=PUY
      TARG(3,2)=(PAP-PAM)*0.5
      TARG(4,2)=(PAP+PAM)*0.5
      KPROJA(2,IPROJ)=KOL
      KTARGA(2,ITARG)=KOL
      IF ( PROJ(3,2)*PROX .LT. 0. ) GOTO 5000
      IF ( TARG(3,2)*TARX .LT. 0. ) GOTO 5000

C  JPSI STRING
C  -----------
      NSTR=NSTR+1
      IF ( NSTR .GT. MXSTR ) THEN
        CALL UTSTOP('PVJPSF: NSTR>MXSTR                      ')
      ENDIF
      IORSTR(NSTR)=-KOL
      ICSTR(1,NSTR)=000100
      ICSTR(2,NSTR)=0
      ICSTR(3,NSTR)=0
      ICSTR(4,NSTR)=000100
      PSTR(1,NSTR)=PX
      PSTR(2,NSTR)=PY
      PSTR(3,NSTR)=PZ
      PSTR(4,NSTR)=E
      PSTR(5,NSTR)=AM
      ROTSTR(1,NSTR)=0.
      ROTSTR(2,NSTR)=0.
      ROTSTR(3,NSTR)=1.
      XORSTR(1,NSTR)=COORD(1,KOL)
      XORSTR(2,NSTR)=COORD(2,KOL)
      XORSTR(3,NSTR)=COORD(3,KOL)
      XORSTR(4,NSTR)=COORD(4,KOL)
      IF ( ISH .GE. 90 ) THEN
        IF ( ISH .GE. 91 ) THEN
          J=NSTR
          WRITE(IFCH,100)J,(ICSTR(K,J)/100,K=1,4)
     *           ,PSTR(3,J),PSTR(4,J),PSTR(5,J)
100       FORMAT(' /CSTR/',I4,3X,4I5,3(E11.3))
        ENDIF

C  CHECKS
C  ------
        IF ( ABS(PIP*PIM-POT2) .GT. 1.E-4 ) THEN
          CALL UTMSG('PVJPSF')
          WRITE(IFCH,*)'*****  PIP*PIM /= POT**2'
          WRITE(IFCH,*)'PIP*PIM=',PIP*PIM,'POT**2=',POT2
          WRITE(IFCH,*)'PIP=',PIP,'   PIM=',PIM
          CALL UTMSGF
        ENDIF
        IF ( ABS(PAP*PAM-PUT2) .GT. 1.E-4 ) THEN
          CALL UTMSG('PVJPSF')
          WRITE(IFCH,*)'*****  PAP*PAM /= PUT**2'
          WRITE(IFCH,*)'PAP*PAM=',PAP*PAM,'PUT**2=',PUT2
          WRITE(IFCH,*)'PAP=',PAP,'   PAM=',PAM
          CALL UTMSGF
        ENDIF
        IF ( ABS(PROJ(4,2)**2
     *      -PROJ(1,2)**2-PROJ(2,2)**2-PROJ(3,2)**2) .GT. 1.E-4 ) THEN
          CALL UTMSG('PVJPSF')
          WRITE(IFCH,*)'*****  MASS**2 OF PROJ NONZERO'
          WRITE(IFCH,*)'MASS**2=',PROJ(4,2)**2
     *              -PROJ(1,2)**2-PROJ(2,2)**2-PROJ(3,2)**2
          CALL UTMSGF
        ENDIF
        IF ( ABS(TARG(4,2)**2
     *       -TARG(1,2)**2-TARG(2,2)**2-TARG(3,2)**2) .GT. 1.E-4 ) THEN
          CALL UTMSG('PVJPSF')
          WRITE(IFCH,*)'*****  MASS**2 OF TARG NONZERO'
          WRITE(IFCH,*)'MASS**2=',TARG(4,2)**2
     *             -TARG(1,2)**2-TARG(2,2)**2-TARG(3,2)**2
          CALL UTMSGF
        ENDIF
        DO 14 N=1,4
          IF ( ABS(PR(N)+TG(N)
     *       -PROJ(N,2)-TARG(N,2)-PSTR(N,NSTR)) .GT. 1.E-4 ) GOTO 15
14      CONTINUE
        GOTO 16
15      CONTINUE
        CALL UTMSG('PVJPSF')
        WRITE(IFCH,*)'*****  PROJ + TARG /= PROJ_NEW + TARG_NEW +JPSI'
        WRITE(IFCH,*)'PROJ,TARG:'
        WRITE(IFCH,*)PR
        WRITE(IFCH,*)TG
        WRITE(IFCH,*)'PROJ_NEW,TARG_NEW,JPSI:'
        WRITE(IFCH,*)(PROJ(N,2),N=1,4)
        WRITE(IFCH,*)(TARG(N,2),N=1,4)
        WRITE(IFCH,*)(PSTR(N,NSTR),N=1,4)
        CALL UTMSGF
16      CONTINUE
      ENDIF

      AMSAC=AMSAC+AM

1000  RETURN
      END
C=======================================================================

      SUBROUTINE RACPRO(TYP,QMU,N,ACPROB)

C-----------------------------------------------------------------------
C  RETURNS THE ARRAY ACPROB CONTAINING ACCUMULATED PROB FOR:
C    EXPONENTIAL OR POISSON DISTRIBUTION (FOR TYP = EXP OR POI)
C    MULTI POMERON CUTS ACC TO GRIBOV (FOR TYP = GRI).
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARO3/   ASUHAX(7),ASUHAY(7),OMEGA,SIGPPD,SIGPPE,UENTRO
     *                ,IWZZZZ
      COMMON /PARO4/   GRICEL,GRIDEL,GRIGAM,GRIRSQ,GRISLO

      DOUBLE PRECISION DENOM,FZ,FZH,Q,QZ,QZH,Z
      REAL             ACPROB(N)
      CHARACTER        TYP*3
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 17 ) ISH=MOD(ISHSUB,100)

C  GRIBOV-REGGE
C  ------------
      IF     ( TYP .EQ. 'GRI' ) THEN
        XI=LOG(ENGY**2)
        AUXIL = EXP(XI*GRIDEL)*GRIGAM
        SIG=8.*PI*AUXIL
        Z=2.D0*GRICEL/(GRIRSQ+GRISLO*XI)*AUXIL

C  0-CUT-POMERON PROBABILITY
C  -------------------------
        QZ=-1.D0/Z
        QZH=-2.D0/Z
        FZ=0.D0
        FZH=0.D0
        DO 10 I=1,20
          DENOM=1.D0/I
          QZ=(-Z)*DENOM*QZ
          QZH=(-Z)*DENOM*0.5*QZH
          FZ=DENOM*QZ+FZ
          FZH=DENOM*QZH+FZH
10      CONTINUE
        RELER=ABS( QZ/(I*FZ) )
        RELERH=ABS( QZH/(I*FZH) )
        SIGPPZ=SIG*(FZH-FZ)
        SIGPPE=SIGPPZ/GRICEL
        SIGPPD=SIGPPE*(GRICEL-1)

C  N-CUT-POMERON PROBABILITIES, N>0
C  --------------------------------
        AUXIL = EXP(-Z)
        Q=1.
        R=1.
        ACPROB(1)=SIG/Z*(1.-AUXIL)
        DO 21 I=2,N
          Q=Q*Z/(I-1)
          R=R+Q
          SI=SIG/(I*Z)*(1.-R*AUXIL)
          ACPROB(I)=ACPROB(I-1)+SI
21      CONTINUE
        IF ( SIGPPI .LT. 0. ) SIGPPI=ACPROB(N)
        AUXIL2 = 1./ACPROB(N)
        DO 22 I=1,N
          ACPROB(I)=AUXIL2*ACPROB(I)
22      CONTINUE

C  WARNINGS
C  --------
        IF ( ISH.GE.90 ) THEN
          IF ( RELER .GT. 1.E-3 ) THEN
            CALL UTMSG('RACPRO')
            WRITE(IFCH,*)'*****  RELER TOO LARGE'
            WRITE(IFCH,*)'RELER,QZ/I,FZ:',RELER,QZ/I,FZ
            CALL UTMSGF
          ENDIF
          IF ( RELERH .GT. 1.E-3 ) THEN
            CALL UTMSG('RACPRO')
            WRITE(IFCH,*)'*****  RELERH TOO LARGE'
            WRITE(IFCH,*)'RELERH,QZH/I,FZH:',RELERH,QZH/I,FZH
            CALL UTMSGF
          ENDIF
          IF ( SIGPPZ .LT. 0. ) THEN
            CALL UTMSG('RACPRO')
            WRITE(IFCH,*)'*****  NEGATIVE SIGPPZ'
            WRITE(IFCH,*)'SIGPPZ,SIG,FZH,FZ:',SIGPPZ,SIG,FZH,FZ
            CALL UTMSGF
          ENDIF
          DO 24 I=1,N
            IF ( ACPROB(I) .LT. 0. ) THEN
              CALL UTMSG('RACPRO')
              WRITE(IFCH,*)'*****  NEGATIVE ACPROB(I)'
              WRITE(IFCH,*)'I,ACPROB(I):',I,ACPROB(I)
              CALL UTMSGF
            ENDIF
24        CONTINUE
          RELERN=ACPROB(N)-ACPROB(N-1)
          IF ( RELERN .GT. 1.E-3 ) THEN
            CALL UTMSG('RACPRO')
            WRITE(IFCH,*)'*****  RELERN TOO LARGE'
            WRITE(IFCH,*)'RELERN:',RELERN
            CALL UTMSGF
          ENDIF
        ENDIF

C  PRINT
C  -----
        IF ( ISH .GE. 92 ) THEN
100       FORMAT(1X,79A1)
          WRITE(IFCH,*)' '
          WRITE(IFCH,100)('-',IC=1,79)
          WRITE(IFCH,*)'   CROSS SECTIONS AND',
     *      ' CUT-POMERON WEIGHTS ACC TO GRIBOV-REGGE-THEORY'
          WRITE(IFCH,100)('-',IC=1,79)
          WRITE(IFCH,*)'   CMS-ENERGY (GEV):',ENGY
          WRITE(IFCH,100)('-',IC=1,79)
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'   PARAMETERS: '
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'     GAMMA: ',GRIGAM
          WRITE(IFCH,*)'     R**2:  ',GRIRSQ
          WRITE(IFCH,*)'     DELTA: ',GRIDEL
          WRITE(IFCH,*)'     SLOPE: ',GRISLO
          WRITE(IFCH,*)'     C:     ',GRICEL
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'   CROSS SECTIONS:'
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'     ELASTIC:     ',SIGPPE
          WRITE(IFCH,*)'     DIFFRACTIVE: ',SIGPPD
          WRITE(IFCH,*)'     INELASTIC:   ',SIGPPI
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'   WEIGHTS W(N) OF N CUT POMERONS:'
          WRITE(IFCH,*)'      ( N - W(N) - W(N)_EXP )'
          WRITE(IFCH,*)' '
          A=QMU/(QMU+1.)
          I=1
          WRITE(IFCH,*)I,ACPROB(I),(1-A)*A**(I-1)
          DO 25 I=2,N
            WRITE(IFCH,*)I,ACPROB(I)-ACPROB(I-1),(1-A)*A**(I-1)
25        CONTINUE
          WRITE(IFCH,*)' '
        ENDIF

C  POISSON
C  -------
      ELSEIF ( TYP .EQ. 'POI' ) THEN
        Z=QMU
        IF ( Z .GE. N-1. ) THEN
          CALL UTSTOP('RACPRO: Z >= N-1.                       ')
        ENDIF
        K=MAX( 0.D0, Z-(N-1-Z) )+1.
        IF ( K .GE. N ) THEN
          CALL UTSTOP('RACPRO: K >= N.                         ')
        ENDIF
        IF     ( K .EQ. 1 ) THEN
          PRBAB=EXP(-Z)
          ACPROB(1)=PRBAB
          DO 1 I=2,N
            PRBAB=PRBAB*Z/(I-1)
            ACPROB(I)=ACPROB(I-1)+PRBAB
1         CONTINUE
        ELSEIF ( K .GT. 1 ) THEN
          X=Z*EXP(-Z/(K-1.))
          PRBAB=1.
          DO 3 I=1,K-1
            PRBAB=PRBAB*X/I
            ACPROB(I)=0.
3         CONTINUE
          ACPROB(K)=PRBAB
          DO 4 I=K+1,N
            PRBAB=PRBAB*Z/(I-1)
            ACPROB(I)=ACPROB(I-1)+PRBAB
4         CONTINUE
        ELSE
          CALL UTSTOP('RACPRO: K <= 0.                         ')
        ENDIF

C  EXPONENTIAL
C  -----------
      ELSEIF ( TYP .EQ. 'EXP' ) THEN
        A=QMU/(QMU+1.)
        PRBAB=1.-A
        ACPROB(1)=PRBAB
        DO 2 I=2,N
          PRBAB=PRBAB*A
          ACPROB(I)=ACPROB(I-1)+PRBAB
2       CONTINUE
      ENDIF

      ISH=ISH0
      RETURN
      END
C=======================================================================

      FUNCTION RANSTC(XFL,XMIN)

C-----------------------------------------------------------------------
C  RETURNS RANDOM NUMBER ACCORDING TO A QUARK STRUCTURE FCTN
C  WITH X>=XMIN
C
C  CHANGES  : D. HECK    IK3  KFK KARLSRUHE
C  DATE     : MAR  22, 1994
C-----------------------------------------------------------------------
      PARAMETER (NSTRU=2049)
      COMMON /CIPIO/   IPIO
      COMMON /CUTINV/  LUTINV
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /STRU/    QSEP(NSTRU),QSET(NSTRU),QVAP(NSTRU)
     *                ,QVAT(NSTRU),XCUTAR,XSTRU(NSTRU)
     *                ,IDTG
      COMMON /STRU2/   DELTA0,DELTA1,QSEH(NSTRU),QSEPI(NSTRU)
     *                ,QVAH(NSTRU),QVAPI(NSTRU),XSE(NSTRU),XVA(NSTRU)
      CHARACTER XFL*3
      SAVE
C-----------------------------------------------------------------------
      NSTRUC = NSTRU
      IF     ( XFL .EQ. 'SEP' ) THEN
        RANSTC=RANXQ(NSTRUC,XSE  ,QSEP,XMIN)
      ELSEIF ( XFL .EQ. 'SET' ) THEN
        RANSTC=RANXQ(NSTRUC,XSE  ,QSET,XMIN)
      ELSEIF ( XFL .EQ. 'VAP' ) THEN
        RANSTC=RANXQ(NSTRUC,XVA  ,QVAP,XMIN)
      ELSEIF ( XFL .EQ. 'VAT' ) THEN
        RANSTC=RANXQ(NSTRUC,XVA  ,QVAT,XMIN)
      ENDIF
      RETURN
      END
C=======================================================================

      FUNCTION RANXQ(N,X,Q,XMIN)

C-----------------------------------------------------------------------
C  RETURNS RANDOM NUMBER ACCORDING TO X(I) Q(I) WITH X>=XMIN
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      REAL X(N),Q(N)
      REAL XMIN,QRAN
      SAVE
C-----------------------------------------------------------------------
      IMIN=1
      IF ( XMIN .EQ. 0. ) GOTO 3
      I1=1
      I2=N
 1    I=I1+(I2-I1)/2
      IF     ( X(I) .LT. XMIN ) THEN
        I1=I
      ELSEIF ( X(I) .GT. XMIN ) THEN
        I2=I
      ELSE
        IMIN=I
        GOTO 3
      ENDIF
      IF ( I2-I1 .GT. 1 ) GOTO 1
      IMIN=I2
 3    CONTINUE
      IF ( Q(IMIN) .GT. Q(N)*.9999 ) THEN
        RANXQ=XMIN
        GOTO 4
      ENDIF
      QRAN=Q(IMIN)+RANGEN()*(Q(N)-Q(IMIN))
      RANXQ=UTINVT(N,X,Q,QRAN)
 4    CONTINUE
      IF ( RANXQ .LT. XMIN ) THEN
        IF(ISH.GE.90)THEN
          CALL UTMSG('RANXQ ')
          WRITE(IFCH,*)'*****  RANXQ=',RANXQ,' <       XMIN=',XMIN
          WRITE(IFCH,*)'Q(IMIN) Q Q(N):',Q(IMIN),QRAN,Q(N)
          WRITE(IFCH,*)'X(IMIN) X X(N):',X(IMIN),RANXQ,X(N)
          CALL UTMSGF
        ENDIF
        RANXQ=XMIN
      ENDIF
      IF ( ISH .GT. 91 ) THEN
        WRITE(IFCH,*)'RANXQ:'
        WRITE(IFCH,*)'   Q(IMIN) Q Q(N):',Q(IMIN),QRAN,Q(N)
        WRITE(IFCH,*)'   X(IMIN) X X(N):',X(IMIN),RANXQ,X(N)
      ENDIF

      RETURN
      END
C=======================================================================

      FUNCTION SBET(Z,W)

C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      SBET=SGAM(Z)*SGAM(W)/SGAM(Z+W)
      RETURN
      END
C=======================================================================

      FUNCTION SDENSI(R)

C-----------------------------------------------------------------------
C  NUCLEAR DENSITY
C-----------------------------------------------------------------------
      COMMON /CDEN/   MASSNR,RMX,R0
C  AI IS INVERSE OF A=0.54
      REAL AI
      DATA AI/1.85185185/
      SAVE
C-----------------------------------------------------------------------
      SDENSI=R**2 / ( 1. + EXP((R-R0)*AI) )
      RETURN
      END
C=======================================================================

      FUNCTION SGAM(X)

C-----------------------------------------------------------------------
C  GAMMA FUNCTION
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      DOUBLE PRECISION GC(13),GF,GP,GX,GY,GZ
      DATA GC/
     *  0.000539698958808D0, 0.002619307282746D0, 0.020449630823590D0,
     *  0.073094836414370D0, 0.279643691578538D0, 0.553387692385769D0,
     *  0.999999999999998D0,-0.000832724708684D0, 0.004698658079622D0,
     *  0.022523834747260D0,-0.170447932874746D0,-0.056810335086194D0,
     *  1.130603357286556D0/
      DATA GP / 3.141592653589793D0 /
      SAVE
C-----------------------------------------------------------------------
      GF = 0.D0
      GX = DBLE(X)
      GZ = GX
      IF ( GX .GT. 0.D0 ) GOTO 1
      IF ( GX .EQ. DINT(GX) ) THEN
        WRITE(IFCH,'('' ARGUMENT OF GAMMA FUNCTION ='',E20.5)') X
        CALL UTSTOP('SGAM : NEGATIVE INTEGER ARGUMENT        ')
         GOTO 5
      ENDIF
      GZ = 1.D0 - GZ
 1    CONTINUE
      GY = 1.D0 / GZ
      IF ( GZ .LE. 1.D0 ) GOTO 4
      GY = 1.D0
 2    CONTINUE
      IF ( GZ .LT. 2.D0 ) GOTO 3
      GZ = GZ - 1.D0
      GY = GY * GZ
      GOTO 2
 3    CONTINUE
      GZ = GZ - 1.D0
 4    CONTINUE
      GF = GY * ((((((GC(1)*GZ+GC(2))*GZ+GC(3))*GZ+GC(4))*GZ+
     *       GC(5))*GZ+GC(6))*GZ+GC(7))/((((((GC(8)*GZ+GC(9))*GZ+
     *       GC(10))*GZ+GC(11))*GZ+GC(12))*GZ+GC(13))*GZ+1.D0)
      IF ( GX .LE. 0.D0 ) GF = GP / ( SIN(GP*GX) * GF )
 5    CONTINUE
      SGAM = GF
      RETURN
      END
C=======================================================================

      FUNCTION SGAU(X)

C-----------------------------------------------------------------------
C  RETURNS GAUSSIAN DISTRIBUTION (NOT NORMALIZED)
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      SGAU=EXP(-0.5*X**2)
      RETURN
      END
C=======================================================================

      SUBROUTINE SHOPAR

C-----------------------------------------------------------------------
      DOUBLE PRECISION SEEDC,SEEDI
      COMMON /CSEED/   SEEDC,SEEDI
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARO3/   ASUHAX(7),ASUHAY(7),OMEGA,SIGPPD,SIGPPE,UENTRO
     *                ,IWZZZZ
      COMMON /PARO4/   GRICEL,GRIDEL,GRIGAM,GRIRSQ,GRISLO
      COMMON /PARO5/   DELEPS,DELVOL
      COMMON /QUARKM/  SMAS,SSMAS,USMAS,UUMAS
      SAVE
C-----------------------------------------------------------------------
      IF ( ISUP .NE. 1  .AND.  ISH .GE. 90 ) THEN
        WRITE(IFMT,102)('*',I=1,69)
     *    ,TAUNLL,MAXRES,PTF,PTQ,IOPTQ,PTMX,NEQMN,IAQU,WTARG
     *    ,WPROJ,QMUST,SIGPPI,CORE,FCTRMX,NCOLMX,LABSYS,IRESCL,OVERLP
     *    ,NTRYMX,DELMSS,SEEDI,GAUMX,BMAXIM,PUD,PSPINL,PSPINH,PISPN
     *    ,NCLEAN,JPSI,CUTMSS,RSTRAS,ISTMAX,TENSN,NEQMX,IPAGI,NDECAY
     *    ,PDIQUA,PAREA,DELREM,TAUMX,NSTTAU,SIGJ,JPSIFI,BMINIM
     *    ,RADIAC,TAUMIN,TAUMAX,NUMTAU,THEMAS,AMSIAC,ELEPTI
102     FORMAT(' ',69A1
     *,/,' *  TAUNLL=',F7.3,'  MAXRES=',I7,'  PTF   =',F7.3
     *,'  PTQ   =',F7.3,'   *',/,' *  IOPTQ =',I7,'  PTMX  =',F7.3
     *,'  NEQMN =',I7,'  IAQU  =',I7,'   *',/,' *  WTARG =',F7.3
     *,'  WPROJ =',F7.3,'  QMUST =',F7.3,'  SIGPPI=',F7.3,'   *'
     *,/,' *  CORE  =',F7.3,'  FCTRMX=',F7.3,'  NCOLMX=',I7
     *,'  LABSYS=',  I7,'   *',/,' *  IRESCL=',I7,'  OVERLP=',F7.3
     *,'  NTRYMX=',I7,'  DELMSS=',F7.3,'   *',/,' *  SEEDI=',D24.17
     *,'  GAUMX =',F7.3,'  BMAXIM=',F7.1,'   *',/,' *  PUD   =',F7.3
     *,'  PSPINL=',F7.3,'  PSPINH=',F7.3,'  PISPN =',F7.3,'   *'
     *,/,' *  NCLEAN=',I7,'  JPSI  =',I7,'  CUTMSS=',F7.3
     *,'  RSTRAS=',F7.3,'   *',/,' *  ISTMAX=',I7,'  TENSN =',F7.3
     *,'  NEQMX =',I7,'  IPAGI =',I7,'   *',/,' *  NDECAY=',I7
     *,'  PDIQUA=',F7.3,'  PAREA =',F7.3,'  DELREM=',F7.3,'   *'
     *,/,' *  TAUMX =',F7.3,'  NSTTAU=',I7,'  SIGJ  =',F7.3
     *,'  JPSIFI=',I7,'   *',/,' *  BMINIM=',F7.1,'  RADIAC=',F7.3
     *,'  TAUMIN=',F7.3,'  TAUMAX=',F7.3,'   *',/,' *  NUMTAU=',I7
     *,'  THEMAS=',F7.3,'  AMSIAC=',F7.3,'  ELEPTI=',F7.3,'   *')
        WRITE(IFMT,108)
     *    ELEPTO,ANGMUE,SMAS,UUMAS,USMAS,SSMAS,IOPBRK,NDECAW
     *   ,IMIHIS,KENTRO,RADIAS,ISPHIS,RHOPHI,ICLHIS,ISHSUB,IOPENU
     *   ,IOPENT,CUTMSQ,KUTDIQ,IDPM,TAUREA,ISPALL,YMXIMI,PTH,PHARD
     *   ,IOPTF,PROSEA,NDECAX,WTMINI,WTSTEP,IWCENT,ISHEVT,PVALEN
     *   ,IFRADE,IENTRO,GRIGAM,GRIRSQ,GRIDEL,GRISLO,GRICEL,IJPHIS
     *   ,UENTRO,IWZZZZ,IOJINT,AMPRIF,DELVOL,DELEPS
     *   ,('*',I=1,69)
108     FORMAT(
     * ' *  ELEPTO=',F7.3,'  ANGMUE=',F7.3,'  SMAS  =',F7.3
     *,'  UUMAS =',F7.3,'   *',/,' *  USMAS =',F7.3,'  SSMAS =',F7.3
     *,'  IOPBRK=',I7,'  NDECAW=',I7  ,'   *',/,' *  IMIHIS=',I7
     *,'  KENTRO=',I7,'  RADIAS=',F7.3,'  ISPHIS=',I7  ,'   *'
     *,/,' *  RHOPHI=',F7.3,'  ICLHIS=',I7  ,'  ISHSUB=',I7
     *,'  IOPENU=',I7  ,'   *',/,' *  IOPENT=',I7,'  CUTMSQ=',F7.3
     *,'  KUTDIQ=',I7  ,'  IDPM  =',I7  ,'   *',/,' *  TAUREA=',F7.3
     *,'  ISPALL=',I7  ,'  YMXIMI=',F7.3,'  PTH   =',F7.3,'   *'
     *,/,' *  PHARD =',F7.3,'  IOPTF =',I7  ,'  PROSEA=',F7.3
     *,'  NDECAX=',I7  ,'   *',/,' *  WTMINI=',F7.3,'  WTSTEP=',F7.3
     *,'  IWCENT=',I7  ,'  ISHEVT=',I7  ,'   *',/,' *  PVALEN=',F7.3
     *,'  IFRADE=',I7  ,'  IENTRO=',I7  ,'  GRIGAM=',F7.3,'   *'
     *,/,' *  GRIRSQ=',F7.3,'  GRIDEL=',F7.3,'  GRISLO=',F7.3
     *,'  GRICEL=',F7.3,'   *',/,' *  IJPHIS=',I7  ,'  UENTRO=',F7.3
     *,'  IWZZZZ=',I7  ,'  IOJINT=',I7  ,'   *',/,' *  AMPRIF=',F7.3
     *,'  DELVOL=',F7.3,'  DELEPS=',F7.3,'         ',7X  ,'   *'
     *,/,' ',69A1)
      ENDIF
      RETURN
      END
C=======================================================================

      FUNCTION SJCENT(K,KU,U)

C----------------------------------------------------------------------
C  RETURNS ENTROPY.
C  INPUT: QUARK NUMBER K; ENERGY U (GEV).
C  IOPENT=1: OSCILLATOR MODEL.
C    INTERPOLATES AND EXTRAPOLATES ENTRO(1+K,1+N)
C    FROM SR JCENTR (JCENTD).
C    SR JCENTR (JCENTD) HAS TO BE CALLED BEFORE!!
C  IOPENT=2,3: FERMI GAS MODEL; JOERG  AICHELIN.
C    IOPENT=2: CONST VOLUME, IOPENT=3: CONST DENSITY
C  IOPENT=4: FERMI GAS (NEW)
C  IOPENT=5: HAGEDORN
C----------------------------------------------------------------------
      PARAMETER (KPARX=15)
      PARAMETER (NQUAX=12)
      COMMON /CENTRO/  ENTRO(1+KPARX,1+NQUAX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /PARO3/   ASUHAX(7),ASUHAY(7),OMEGA,SIGPPD,SIGPPE,UENTRO
     *                ,IWZZZZ
      SAVE
C----------------------------------------------------------------------
      SJCENT=0.
      IF     ( IOPENT .EQ. 0 ) THEN
        RETURN

      ELSEIF ( IOPENT .EQ. 1 ) THEN
        IF ( K .LE. KENTRO ) RETURN
        IF ( MOD(K,3) .NE. 0 ) THEN
          CALL UTSTOP('SJCENT: K MUST BE MULTIPLE OF 3         ')
        ENDIF
        E=U/OMEGA
        IF ( K .GT. KPARX ) GOTO 5001
        N=INT(E)
        IF ( N .LT. 0 ) THEN
          CALL UTSTOP('SJCENT: NEGATIVE ENERGY                 ')
        ENDIF
        IF ( N .GE. NQUAX ) N=NQUAX-1
        SJCENT=ENTRO(1+K,1+N)+(E-N)*(ENTRO(1+K,1+N+1)-ENTRO(1+K,1+N))
        IF ( K .LE. 3 ) RETURN
5002    E3=E*3./K
        N3=INT(E3)
        IF ( N3 .GE. NQUAX ) N3=NQUAX-1
        SENTR3=ENTRO(1+3,1+N3)+(E3-N3)
     *                              *(ENTRO(1+3,1+N3+1)-ENTRO(1+3,1+N3))
        SJCENT=MIN(SJCENT,K/3.*SENTR3)
        RETURN
5001    CONTINUE
        L=KPARX
        EL=(E*L)/K
        NL=INT(EL)
        IF ( NL .GE. NQUAX ) NL=NQUAX-1
        SENTRL=ENTRO(1+L,1+NL)+(EL-NL)
     *                              *(ENTRO(1+L,1+NL+1)-ENTRO(1+L,1+NL))
        SJCENT=FLOAT(K)/L*SENTRL
        GOTO 5002

      ELSEIF ( IOPENT .EQ. 2 ) THEN
C  CONSTANT VOLUME 5 FM**3
        IF ( K .LE. KENTRO ) RETURN
        X1=12.96*K**(-.315)
        IF ( K .LT. 19 ) THEN
          X2=.785 + .005*K
        ELSE
          X2=.88
        ENDIF
        SJCENT=K*X1*(U/K)**X2
        RETURN

      ELSEIF ( IOPENT .EQ. 3 ) THEN
C  CONSTANT DENSITY
        IF ( K .LE. KENTRO ) RETURN
        X1=9.785
        X2=.7926
        SJCENT=K*X1*(U/K)**X2
        RETURN

      ELSEIF ( IOPENT .EQ. 4 ) THEN
C  CONSTANT VOLUME 15 FM**3
        IF ( K .LE. KENTRO ) RETURN
        SJCENT=SJCEN4(K,KU,U)
        RETURN

      ELSEIF ( IOPENT .EQ. 5 ) THEN
C  HAGEDORN
        IF ( U .LE. UENTRO  .AND.  K .LE. KENTRO ) THEN
          RETURN
        ENDIF
CDH     THAGED=.250
CDH     SJCENT = U/THAGED
        SJCENT = U*4.

      ELSE
        CALL UTSTOP('SJCENT: INVALID OPTION IOPENT           ')
      ENDIF

      RETURN
      END
C=======================================================================

      FUNCTION SJCEN4(K,KU,U)

C----------------------------------------------------------------------
C  RETURNS TOTAL ENTROPY.
C  CONSTANT VOLUME 15 FM**3.
C  INPUT: TOTAL QUARK NUMBER K; UP AND DOWN QUARKS KU;
C  TOTAL EXCITATION ENERGY U(GEV).
C  INITIAL CALL OF SJCEN0 REQUIRED!!!
C----------------------------------------------------------------------
      COMMON /CSJCEN/  ENT(16000)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C----------------------------------------------------------------------
      SJCEN4 = 0.

      IF ( K .LT. KU ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('SJCEN4')
          WRITE(IFCH,*)'*****  K < KU'
          WRITE(IFCH,*)'*****  K: ',K,'   KU: ',KU
          CALL UTMSGF
        ENDIF
        RETURN
      ENDIF

      IF ( K .GT. 45 ) THEN
        KO=K
        K=45
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('SJCEN4')
          WRITE(IFCH,*)'*****  K > 45'
          WRITE(IFCH,*)'*****  K: ',KO,'   K_NEW: ',K
          CALL UTMSGF
        ENDIF
      ENDIF

      IF ( MOD(K,3) .NE. 0 ) THEN
        CALL UTSTOP('SJCEN4: NONINTEGER BARYON NUMBER        ')
      ENDIF

      IK=0
      DO 1 I=3,K-3,3
        IK=IK+I+1
 1    CONTINUE
      IF ( U .LT. 10. ) THEN
        DU=MOD(U,.25)
        IU=INT(U/.25)+1
        IENTRY=(IK+KU)*41+IU
        SJCEN4=ENT(IENTRY)+(ENT(IENTRY+1)-ENT(IENTRY))*DU*4.
      ELSE
        IU=40
        IENTRY=(IK+KU)*41+IU
        SJCEN4=ENT(IENTRY)+(ENT(IENTRY+1)-ENT(IENTRY))*4.*(U-9.75)
      ENDIF

      RETURN
      END
C=======================================================================

      FUNCTION SJCGAM(KEUX,KEDX,KESX,KECX,AMA,AMO,PO,MOX)

C----------------------------------------------------------------------
C  RETURNS PARTIAL DECAY WIDTH DGAMMA = PHASE SPACE * DENSITY
C    FOR DECAY OF CLUSTER INTO CLUSTER AND HADRON.
C  KE*X: NET QUARK NUMBER
C  AMA: CLUSTER MASS;  AMO: HADRON MASS;  PO: HADRON MOMENTUM
C----------------------------------------------------------------------
      COMMON /CENTEX/  ENTEXP
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CSJCGA/  AMEGAM,AMNULL,ASUHA(7),ENTRPY,NOPHA,NSUHA(7)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO3/   ASUHAX(7),ASUHAY(7),OMEGA,SIGPPD,SIGPPE,UENTRO
     *                ,IWZZZZ
      SAVE
C----------------------------------------------------------------------
      ENTRPY=-99.9
      IF ( MOX .EQ. 1 ) THEN
        ENTEXP=0.
        AMNULL=UTAMNU(KEUX,KEDX,KESX,KECX,IOPENU)
      ELSE
        ENTRPY=ENTRPY-ENTEXP
      ENDIF
      EO=SQRT(AMO**2+PO**2)
      IF ( EO .GT. AMA ) GOTO 999
      AME2=(AMA-EO)**2-PO**2
      IF ( AME2 .LT. 0. ) GOTO 999
      AMEGAM=SQRT(AME2)
      E=AMEGAM-AMNULL
      IF ( E .LT. 0. ) GOTO 998
      KE=ABS(KEUX+KEDX+KESX+KECX)

      IF ( IOPENT .EQ. 5 ) THEN
        IF ( E .LE. UENTRO  .AND.  KE .LE. KENTRO ) THEN
          ENTRPY=0.
        ELSE
C  ENTROPY AFTER HAGEDORN
          ENTRPY = E*4.
        ENDIF
      ELSE
        KU=ABS(KEUX+KEDX)
        ENTRPY=SJCENT(KE,KU,E)
      ENDIF

      IF ( MOX .EQ. 1 ) THEN
        ENTEXP=ENTRPY
        ENTRPY=0.
        SJCGAM=.125* PO**2/( PI**2*AMA*EO )
      ELSE
        ENTRPY=ENTRPY-ENTEXP
        SJCGAM=.125*EXP(ENTRPY)* PO**2/( PI**2*AMA*EO )
      ENDIF
      RETURN
 999  AMEGAM = 0.
 998  SJCGAM= 0.
      RETURN
      END
C=======================================================================

      FUNCTION SMASS(A,Y,Z)

C-----------------------------------------------------------------------
C  RETURNS CLUSTER MASS (IN GEV) (PER CLUSTER, NOT (!) PER NUCLEON)
C  ACCORDING TO BERGER/JAFFE MASS FORMULA, PRC35(1987)213 EQ.2.31,
C  SEE ALSO C. DOVER, BNL-46322, INTERSECTIONS-MEETING, TUCSON, 91.
C  A: MASSNR, Y: HYPERCHARGE, Z: CHARGE,
C-----------------------------------------------------------------------
      COMMON /CMASS/   AC,AS,CZ,DY,DZ,EPSI,RZERO,SIGMA,THET,YM,ZM
      SAVE
C-----------------------------------------------------------------------
      YMIN=YM*A
      ZMIN=CZ/(DZ/A+ZM/A**.3333333)
      SMASS=EPSI*A +AS*A**.6666667
     *             +(AC/A**.3333333 +DZ/A*0.5)*(Z-ZMIN)**2
     *              +DY/A*0.5*(Y-YMIN)**2
      RETURN
      END
C=======================================================================

      SUBROUTINE SMASSI(THETA)

C-----------------------------------------------------------------------
C  INITIALIZATION FOR SMASS.
C  CALCULATES PARAMETERS FOR BERGER/JAFFE MASS FORMULA
C  (PRC35(1987)213 EQ.2.31, SEE ALSO C. DOVER, BNL-46322).
C  THETA: PARAMETER THAT DETERMINES ALL PARAMETERS IN MASS FORMULA.
C-----------------------------------------------------------------------
      COMMON /CMASS/   AC,AS,CZ,DY,DZ,EPSI,RZERO,SIGMA,THET,YM,ZM
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      DATA ALP /0.007297145/
      SAVE
C-----------------------------------------------------------------------
      THET=THETA

      ASTR=.150

      CO=COS(THETA)
      SI=SIN(THETA)
      BET=(1.+CO**3)*0.5
      RZERO=SI/ASTR/(0.6666667/PI*(1.+CO**3)  )**0.3333333
      CS=ASTR/SI
      CZ=CS*(1.-BET**0.3333333 )
      SIGMA=0.75/PI*(ASTR/SI)**3*0.3333333*(CO**2*0.5 -SI**2*(1-SI)
     *    -1./PI*(PI*0.5-THETA-SIN(2*THETA)+SI**3*LOG((1+CO)/SI) ) )

      EPSI=ASTR*(BET**0.3333333+2.)/SI
      AS=4.*PI*SIGMA*RZERO**2
      AC=0.6*ALP/RZERO
      DZ=ASTR/SI*BET**0.3333333 *CO**2 *
     *  (CO**4*(1.+BET**0.3333333)+(1+BET)**2)/
     *( (2.*CO**2+BET**.3333333)*(CO**4*(1+BET**.6666667)+(1.+BET)**2)-
     *        (CO**4+BET**0.3333333*(1+BET))
     *                     *( (2.*BET**0.6666667-1.)*CO**2+1.+BET) )
      DY=ASTR/6.*(1.+CO**3)**3/SI*
     *       (  1.+(1.+CO)/(4.*(1.+CO**3))**0.6666667  )/
     *       ( CO**6 + CO + CO*(.5*(1+CO**3))**1.333333 )
      ZM=6.*ALP/(5.*RZERO)
      YM=(1.-CO**3)/(1.+CO**3)

      RETURN
      END
C=======================================================================

      SUBROUTINE SMASSP

C-----------------------------------------------------------------------
C  PRINTS SMASS.
C-----------------------------------------------------------------------
      COMMON /CMASS/   AC,AS,CZ,DY,DZ,EPSI,RZERO,SIGMA,THET,YM,ZM
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      REAL ENG(14),YMI(14),ZMI(14)
      SAVE
C-----------------------------------------------------------------------
      WRITE(IFCH,*)'PARAMETERS OF MASS FORMULA:'
      WRITE(IFCH,*)'---------------------------'
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'THETA=',THET,'   EPSI=',EPSI
      WRITE(IFCH,*)'AS=',AS,'   AC=',AC
      WRITE(IFCH,*)'DY=',DY,'   DZ=',DZ
      WRITE(IFCH,*)'YM=',YM
      WRITE(IFCH,*)'CZ DZ ZM=',CZ,DZ,ZM
      WRITE(IFCH,*)'SIGMA**1/3=',SIGMA**(1./3.),'   RZERO=',RZERO
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'MASS:'
      WRITE(IFCH,*)'-----'
      WRITE(IFCH,5000)(J,J=1,14)
5000  FORMAT(/,5X,'A:',14I5,/)
      DO 4 J=1,14
        A=J
        YMI(J)=YM*A
        ZMI(J)=CZ/(DZ/A+ZM/A**0.3333333)
 4    CONTINUE
      WRITE(IFCH,5002)(YMI(J),J=1,14)
5002  FORMAT(1X,'YMIN: ',14F5.2,/)
      WRITE(IFCH,5003)(ZMI(J),J=1,14)
5003  FORMAT(1X,'ZMIN: ',14F5.2,/)
      DO 2 I=1,15
        NS=11-I
        DO 3 J=1,14
          A=J
          Y=A-NS
          Z=0.
          ENG(J)=SMASS(A,Y,Z)/A
 3      CONTINUE
        WRITE(IFCH,5001)NS,(ENG(J),J=1,14)
5001    FORMAT(1X,'S=',I2,2X,14F5.2)
 2    CONTINUE
      WRITE(IFCH,*)' '
      WRITE(IFCH,*)'MASS-MASS(FREE):'
      WRITE(IFCH,*)'----------------'
      WRITE(IFCH,5000)(J,J=1,14)
      DO 5 I=1,15
        NS=11-I
        DO 6 J=1,14
          A=J
          Y=A-NS
          Z=0.
          SG=SIGN(1.,A)
          AX=SG*A
          YX=SG*Y
          ZX=SG*Z
          KU=NINT(AX+ZX)
          KD=NINT(AX-ZX+YX)
          KS=NINT(AX-YX)
          KC=0
          ENG(J)=(SMASS(A,Y,Z)-UTAMNU(KU,KD,KS,KC,3))/A
 6      CONTINUE
        WRITE(IFCH,5001)NS,(ENG(J),J=1,14)
5     CONTINUE

      RETURN
      END
C=======================================================================

      SUBROUTINE SMASST(KUX,KDX,KSX,KCX,A,Y,Z)

C-----------------------------------------------------------------------
C  INPUT: KUX,KDX,KSX,KCX = NET QUARK NUMBERS (FOR U,D,S,C QUARKS).
C  OUTPUT: MASSNR A, HYPERCHARGE Y AND CHARGE Z.
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      SG=1.
      IF ( KUX+KDX+KSX+KCX .LT. 0 ) SG=-1.
      KU=SG*KUX
      KD=SG*KDX
      KS=SG*KSX
      KC=SG*KCX
      K=KU+KD+KS+KC
      IF ( MOD(K,3) .NE. 0 ) THEN
        CALL UTSTOP('SMASST: NONINTEGER BARYON NUMBER        ')
      ENDIF
      A=K/3
      Y=A-KS
      NZ=2*KU-KD-KS+2*KC
      IF ( MOD(NZ,3) .NE. 0 ) THEN
        CALL UTSTOP('SMASST: NONINTEGER CHARGE               ')
      ENDIF
      Z=NZ/3
      RETURN
      END
C=======================================================================

      FUNCTION SPTF(X)

C-----------------------------------------------------------------------
C  RETURNS PT-DISTRIBUTION FOR FRAGMENTATION
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CQUAMA/  QUAMA
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      SAVE
C-----------------------------------------------------------------------
      IF     ( IOPTF .EQ. 1 ) THEN
        SPTF = X * EXP(-2./PTF*SQRT(X**2+QUAMA**2))
      ELSEIF ( IOPTF .EQ. 2 ) THEN
        SPTF = X * EXP(-PI*0.25*(X**2+QUAMA**2)/PTF**2)
      ENDIF
      RETURN
      END
C=======================================================================

      FUNCTION SPTH(X)

C-----------------------------------------------------------------------
C  RETURNS PT-DISTRIBUTION FOR HARD SCATTERING
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      SAVE
C-----------------------------------------------------------------------
      SPTH=PTH**2*2.*X/(X**2+PTH**2)**3
      RETURN
      END
C=======================================================================

      FUNCTION SPTJ(X)

C-----------------------------------------------------------------------
C  JPSI PT-DISTRIBUTION IN 200 GEV PP
C-----------------------------------------------------------------------
      DATA AI/1.052631579/, C/2.75482/, C1/1.75482/, CC/16.30716419/
C  AI IS THE INVERSE OF A
      SAVE
C-----------------------------------------------------------------------
CDH   A=0.95
CDH   C=1./0.363
CDH   Z=X/A
      Z=X*AI
CDH   SPTJ=1./A*C**C/SGAM(C)*Z**(C-1.)*EXP(-C*Z)
      SPTJ = AI * CC * 0.619299158 * Z**C1 * EXP(-C*Z)
      RETURN
      END
C=======================================================================

      FUNCTION SPTQ(X)

C-----------------------------------------------------------------------
C  RETURNS PT-DISTRIBUTION OF QUARKS IN NUCLEONS
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      SAVE
C-----------------------------------------------------------------------
      IF     ( IOPTQ .EQ. 2 ) THEN
        AUXIL2 = PI/PTQ**2
        SPTQ=(0.5*X)* AUXIL2 * EXP(-(0.5*X)**2*AUXIL2)
      ELSEIF ( IOPTQ .EQ. 1 ) THEN
        AUXIL1 = 2./PTQ
        SPTQ=AUXIL1**2 * X * EXP(-X*AUXIL1)
      ELSEIF ( IOPTQ .EQ. 3 ) THEN
        SPTQ=PTQ**2 * 2. * X / (X**2+PTQ**2)**2
      ENDIF
      RETURN
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION SSE0(Z)

C-----------------------------------------------------------------------
C  SEA QUARK STRUCTURE FUNCTION FOR HADRONS
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      DOUBLE PRECISION Z
      SAVE
C-----------------------------------------------------------------------
      SSE0 = (1.D0-Z)**8.05D0 * 1.265D0 / SQRT(Z**2+XCUT**2)
      RETURN
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION SSE1(Z)

C-----------------------------------------------------------------------
C  SEA QUARK STRUCTURE FUNCTION FOR PIONS
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      DOUBLE PRECISION Z
      SAVE
C-----------------------------------------------------------------------
      SSE1 = (1.D0-Z)**5.D0 * 0.9D0 / SQRT(Z**2+XCUT**2)
      RETURN
      END
C=======================================================================

      FUNCTION SSPLIT(X)

C-----------------------------------------------------------------------
C  RETURNS SPLITTING FUNCTION
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
C-C   SSPLIT= ( 0.3 + 10.0*X**7 - 10.3*X**15 ) /SQRT(X**2+.2**2)
C-C   SSPLIT=   0.2 +  2.4*X**2 -  2.6*X**5
C-C   SSPLIT=   0.4 +  2.0*X    -  2.4*X**3
C-C   SSPLIT= ( 0.8 +  2.0*X    -  2.8*X**3  ) /SQRT(X**2+.2**2)
      SSPLIT= X**3 - X**5
      RETURN
      END
C=======================================================================

      FUNCTION SSPLIX(X)

C-----------------------------------------------------------------------
C  RETURNS SPLITTING FUNCTION
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
C-C   SSPLIX= (1-X) - (1-X)**7
C-C   SSPLIX= 1
C-C   Z=2*X-.5
C-C   SSPLIX= Z**3 - Z**5
C-C   IF ( Z .LT. 0. ) SSPLIX=0.
C-C   IF ( Z .GT. 1. ) SSPLIX=0.
      SSPLIX= X**3 - X**5
      RETURN
      END
C=======================================================================

      SUBROUTINE STAA(X,Q2I,Z,S)

C-----------------------------------------------------------------------
C  STRUCTURE FUNCTIONS.
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      SAVE
C-----------------------------------------------------------------------
      Z=SQRT(X**2+XCUT**2)
      Q2=MAX(4.00001,Q2I)
CDH   S=LOG(LOG(Q2/.2**2)/LOG(4/.2**2))
      S=LOG(LOG(Q2*25.  )*0.21714724)
      RETURN
      END
C=======================================================================

      FUNCTION STXD(X,Q2)

C-----------------------------------------------------------------------
C  VALENCE D-QUARK DISTRIBUTION.
C  FROM GLUECK, HOFFMANN, REYA, Z. PHYS. C13 (1982) 119.
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      CALL STAA(X,Q2,Z,S)
      A=.364-.0368*S
      C=2.-.5414*S**.8
      D=5.09+.3463*S
      STXD=C*X**A*(1.-X**C)**D/SBET(D+1.,A/C)
      RETURN
      END
C=======================================================================

      FUNCTION STXS(X,Q2)

C-----------------------------------------------------------------------
C  S-SEA DISTRIBUTION.
C  FROM GLUECK, HOFFMANN, REYA, Z. PHYS. C13 (1982) 119.
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      CALL STAA(X,Q2,Z,S)
      A=.0625+.1132*S**1.3
      B=12.64*S-51.70*S**1.8+38.02*S**2
      C=4.448*S
      D=7.+1.562*S
      E=.3081*S**2.5
      F=47.24+67.91*S
      STXS=A*(1.+B*X+C*X**2)*(1.-X)**D + E*EXP(-F*X)
      RETURN
      END
C=======================================================================

      FUNCTION STXU(X,Q2)

C-----------------------------------------------------------------------
C  VALENCE U-QUARK DISTRIBUTION.
C  FROM GLUECK, HOFFMANN, REYA, Z. PHYS. C13 (1982) 119.
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      CALL STAA(X,Q2,Z,S)
      A=.421-.0412*S
      C=2.-.6223*S**.8
      D=3.37+.4319*S
      STXU=2.*C*X**A*(1.-X**C)**D/SBET(D+1.,A/C)
      RETURN
      END
C=======================================================================

      FUNCTION STXUS(X,Q2)

C-----------------------------------------------------------------------
C  U-SEA DISTRIBUTION.
C  FROM GLUECK, HOFFMANN, REYA, Z. PHYS. C13 (1982) 119.
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      CALL STAA(X,Q2,Z,S)
      A=.25+.088*S**1.3
      B=.8128*S-2.003*S**1.8+.0831*S**2
      C=3.97*S
      D=7.+1.666*S
      E=.2487*S**2.5
      F=27.8+59.68*S
      STXUS=A*(1.+B*X+C*X**2)*(1.-X)**D + E*EXP(-F*X)
      RETURN
      END
C=======================================================================

      FUNCTION STXZNE(X,Q2)

C-----------------------------------------------------------------------
C  STRUCTURE FUNCTION OF NEUTRON
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      STXZNE=(4.*STXD(X,Q2) + STXU(X,Q2) +
     *        10.*STXUS(X,Q2) + 2.*STXS(X,Q2))/9.
      RETURN
      END
C=======================================================================

      FUNCTION STXZPR(X,Q2)

C-----------------------------------------------------------------------
C  STRUCTURE FUNCTION OF PROTON
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      STXZPR=(4.*STXU(X,Q2) + STXD(X,Q2) +
     *        10.*STXUS(X,Q2) + 2.*STXS(X,Q2))/9.
      RETURN
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION SVA0(Z)

C-----------------------------------------------------------------------
C  VALENCE QUARK STRUCTURE FUNCTION FOR HADRONS
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      DOUBLE PRECISION Z
      SAVE
C-----------------------------------------------------------------------
      IF ( Z .NE. 0.D0 ) THEN
        SVA0=(1.D0-Z)**3.46 * Z**0.419 * (2.74793064D0*Z + 0.62452969D0)
     *                   / SQRT(Z**2+XCUT**2)
      ELSE
        SVA0=0.D0
      ENDIF
      RETURN
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION SVA1(Z)

C-----------------------------------------------------------------------
C  VALENCE QUARK STRUCTURE FUNCTION FOR PIONS
C-----------------------------------------------------------------------
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      DOUBLE PRECISION Z
      SAVE
C-----------------------------------------------------------------------
      IF ( Z .NE. 0.D0 ) THEN
        SVA1 = (1.D0-Z)**0.7D0 * Z**0.4D0 * 0.1730725D0
     *                                   / SQRT(Z**2+XCUT**2)
      ELSE
        SVA1=0.D0
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UINTEG(VAL,FUNC,A,B,AERR,RERR,LEVEL,ERROR,IFLAG)

C-----------------------------------------------------------------------
C  CACLULATION OF DEFINITE INTEGRAL OF FUNC(X) FROM A TO B
C-----------------------------------------------------------------------
C        RELERR=1.D-12
C        ABSERR=0.D0
C        LEVEL=1
C        CALL UINTEG(VALUE,FUNCTN,A,B,ABSERR,RELERR,LEVEL,ERROR,IFLAG)
C        IF (IFLAG.GT.3) WRITE(*,'('' IFLAG ='',I7)') IFLAG
C-----------------------------------------------------------------------
      IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(10,10),R(10),AIT(10),DIF(10),RN(4)
      DIMENSION TS(2049),IBEGS(30),BEGIN(30),FINIS(30),EST(30)
      LOGICAL   H2CONV,AITKEN,RIGHT,REGLAR,REGLSV(30)
      DATA      TOLSAV,AITLOW,H2TOL,AITTOL,VJUMP,MAXTS,MAXTBL,MXSTGE
     *          / 1.D-16, 1.1D0, .15D0, .1D0, .01D0, 2049, 10, 30 /
      DATA      RN /.71420053D0,.34662815D0,.843751D0,.12633046D0 /
      DATA      ALG402 /.3010299956639795D0 /
      SAVE
C-----------------------------------------------------------------------
      TOLMCH=TOLSAV
      VAL=0.D0
      ERROR=0.D0
      IFLAG=1
      VLONG=DABS(B-A)
      IF ( VLONG .EQ. 0.D0 ) RETURN
      ERRR=DMIN1( .1D0, DMAX1(DABS(RERR),1.D1*TOLMCH) )
      ERRA=DABS(AERR)
      STEPMN=DMAX1(VLONG/FLOAT(2**MXSTGE),
     *   DMAX1(VLONG,DABS(A),DABS(B))*TOLMCH)
      STAGE=.5D0
      ISTAGE=1
      CUREST=0.D0
      FNSIZE=0.D0
      PREVER=0.D0
      REGLAR=.FALSE.
      BEGI=A
      FBEG=FUNC(A)/2.D0
      TS(1)=FBEG
      IBEG=1
      ENDE=B
      FEND=FUNC(B)/2.D0
      TS(2)=FEND
      IEND=2
 60   CONTINUE
      RIGHT=.FALSE.
 61   CONTINUE
      STEP=ENDE-BEGI
      ASTEP=DABS(STEP)
      IF ( ASTEP .LT. STEPMN ) GOTO 97
      IF ( LEVEL .GE. 3 ) WRITE(*,101) BEGI,STEP,ISTAGE
101   FORMAT(' BEGI,STEP',1P,2E16.7,I5)
      T(1,1)=FBEG+FEND
      TABS=DABS(FBEG)+DABS(FEND)
      L=1
      N=1
      H2CONV=.FALSE.
      AITKEN=.FALSE.
      GOTO 63
 62   IF ( LEVEL .GE. 4 ) WRITE(*,102) L,T(1,LM1)
102   FORMAT(1X,I5,7E16.8,/,1X,3E16.8)
 63   LM1=L
      L=L+1
      N2=N*2
      FN=N2
      ISTEP=(IEND-IBEG)/N
      IF ( ISTEP .LE. 1 ) THEN
        II=IEND
        IEND=IEND+N
        IF ( IEND .GT. MAXTS ) GOTO 96
        HOVN=STEP/FN
        III=IEND
        DO  64  I=1,N2,2
          TS(III)=TS(II)
          TS(III-1)=FUNC(ENDE-FLOAT(I)*HOVN)
          III=III-2
          II=II-1
 64     CONTINUE
        ISTEP=2
      ENDIF
      ISTEP2=IBEG+ISTEP/2
      SUM=0.D0
      SUMABS=0.D0
      DO  65  I=ISTEP2,IEND,ISTEP
        SUM=SUM+TS(I)
        SUMABS=SUMABS+DABS(TS(I))
 65   CONTINUE
      T(L,1)=T(L-1,1)/2.D0 + SUM/FN
      TABS=TABS/2.D0+SUMABS/FN
      ABSI=ASTEP*TABS
      N=N2
      IT=1
      VINT=STEP*T(L,1)
      TABTLM=TABS*TOLMCH
      FNSIZE=DMAX1(FNSIZE,DABS(T(L,1)))
      ERGOAL=DMAX1(ASTEP*TOLMCH*FNSIZE,
     *       STAGE*DMAX1(ERRA,ERRR*DABS((CUREST)+VINT)))
      FEXTRP=1.D0
      DO  66  I=1,LM1
        FEXTRP=FEXTRP*4.D0
        T(I,L)=T(L,I)-T(L-1,I)
        T(L,I+1)=T(L,I)+T(I,L)/(FEXTRP-1.D0)
 66   CONTINUE
      ERRER=ASTEP*DABS(T(1,L))
      IF ( L .LE. 2 ) THEN
        IF ( DABS(T(1,2)) .LE. TABTLM)  GOTO 81
        GOTO 63
      ENDIF
      DO  67  I=2,LM1
        DIFF=0.D0
        IF ( DABS(T(I-1,L)) .GT. TABTLM ) DIFF=T(I-1,LM1)/T(I-1,L)
        T(I-1,LM1)=DIFF
 67   CONTINUE
      IF ( DABS(4.D0-T(1,LM1)) .LE. H2TOL ) GOTO 69
      IF ( T(1,LM1) .EQ. 0.D0 ) GOTO 68
      IF ( DABS(2.D0-DABS(T(1,LM1))) .LT. VJUMP ) GOTO 80
      IF ( L .EQ. 3 ) GOTO 62
      H2CONV=.FALSE.
      IF ( DABS((T(1,LM1)-T(1,L-2))/T(1,LM1)) .LE. AITTOL ) GOTO 72
      IF ( REGLAR ) GOTO 68
      IF ( L .EQ. 4 ) GOTO 62
 68   IF ( ERRER .LE. ERGOAL ) GOTO 83
      IF ( LEVEL .GE. 4 ) WRITE(*,103) L,T(1,LM1)
103   FORMAT(1X,I5,7E16.8,/,1X,3E16.8)
      GOTO 90
 69   CONTINUE
      IF ( LEVEL .GE. 4 ) WRITE(*,104) L,T(1,LM1)
104   FORMAT(1X,I5,E16.8,5X,'H2CONV')
      IF ( H2CONV ) GOTO 70
      AITKEN=.FALSE.
      H2CONV=.TRUE.
      IF ( LEVEL .GE. 3 ) WRITE(*,105) L
105   FORMAT(' H2 CONVERGENCE AT ROW',I4)
 70   FEXTRP=4.D0
 71   IT=IT+1
      VINT=STEP*T(L,IT)
      ERRER=DABS(STEP/(FEXTRP-1.D0)*T(IT-1,L))
      IF ( ERRER .LE. ERGOAL ) GOTO 86
      IF ( IT .EQ. LM1 ) GOTO 79
      IF ( T(IT,LM1) .EQ. 0.D0 ) GOTO 71
      IF ( T(IT,LM1) .LE. FEXTRP ) GOTO 79
      IF ( DABS(T(IT,LM1)/4.D0-FEXTRP)/FEXTRP .LT. AITTOL )
     *                                       FEXTRP=FEXTRP*4.D0
      GOTO 71
 72   IF ( LEVEL .GE. 4 ) WRITE(*,106) L,T(1,LM1)
106   FORMAT(1X,I5,E16.8,5X,'AITKEN')
      IF ( T(1,LM1) .LT. AITLOW ) GOTO 90
      IF ( AITKEN ) GOTO 73
      H2CONV=.FALSE.
      AITKEN=.TRUE.
      IF ( LEVEL .GE. 3 ) WRITE(*,107) L
107   FORMAT(' AITKEN AT ROW',I4)
 73   FEXTRP=T(L-2,LM1)
      IF ( FEXTRP .GT. 4.5 ) GOTO 70
      IF ( FEXTRP .LT. AITLOW ) GOTO 90
      IF ( DABS(FEXTRP-T(L-3,LM1))/T(1,LM1) .GT. H2TOL ) GOTO 90
      IF ( LEVEL .GE. 3 ) WRITE(*,108) FEXTRP
108   FORMAT(' RATIO',F13.8)
      SING=FEXTRP
      FEXTM1=FEXTRP-1
      AIT(1)=0.
      DO  74  I=2,L
        AIT(I)=T(I,1)+(T(I,1)-T(I-1,1))/FEXTM1
        R(I)=T(1,I-1)
        DIF(I)=AIT(I)-AIT(I-1)
 74   CONTINUE
      IT=2
 75   VINT=STEP*AIT(L)
      IF ( LEVEL .GE. 5 ) THEN
        WRITE(*,109) (R(I+1),I=IT,LM1)
109     FORMAT(1X,8E15.8)
        WRITE(*,109) (AIT(I),I=IT,L)
        WRITE(*,109) (DIF(I+1),I=IT,LM1)
      ENDIF
      ERRER=ERRER/FEXTM1
      IF ( ERRER .GT. ERGOAL ) GOTO 76
      ALPHA=DLOG10(SING)/ALG402-1.D0
      IF ( LEVEL .GE. 2 ) WRITE(*,110) ALPHA,BEGI,ENDE
110   FORMAT(11X,'INTEGRAND SHOWS SINGULAR ',
     *   'BEHAVIOUR OF TYPE X**(',F5.2,') BETWEEN',1P,E15.7,
     *   ' AND',1P,E15.7)
      IFLAG=MAX0(IFLAG,2)
      GOTO 86
 76   IT=IT+1
      IF ( IT .EQ. LM1 ) GOTO 79
      IF ( IT .LE. 3 ) THEN
        H2NEXT=4.D0
        SINGNX=2.D0*SING
      ENDIF
      IF ( H2NEXT .GE. SINGNX)  THEN
        FEXTRP=SINGNX
        SINGNX=2.D0*SINGNX
      ELSE
        FEXTRP=H2NEXT
        H2NEXT=4.D0*H2NEXT
      ENDIF
      DO 77 I=IT,LM1
        R(I+1)=0.D0
        IF ( DABS(DIF(I+1)) .GT. TABTLM ) R(I+1)=DIF(I)/DIF(I+1)
 77   CONTINUE
      IF ( LEVEL .GE. 4 ) WRITE(*,111) FEXTRP,R(L-1),R(L)
111   FORMAT(' FEXTRP + RATIOS',1P,3E15.7)
      H2TFEX=-H2TOL*FEXTRP
      IF ( R(L)-FEXTRP .LT. H2TFEX ) GOTO 79
      IF ( R(L-1)-FEXTRP .LT. H2TFEX ) GOTO 79
      ERRER=ASTEP*DABS(DIF(L))
      FEXTM1=FEXTRP-1.D0
      DO 78 I=IT,L
        AIT(I)=AIT(I)+DIF(I)/FEXTM1
        DIF(I)=AIT(I)-AIT(I-1)
 78   CONTINUE
      GOTO 75
 79   FEXTRP=DMAX1(PREVER/ERRER,AITLOW)
      PREVER=ERRER
      IF ( L .LT. 5 ) GOTO 63
      IF ( LEVEL .GE. 3 ) WRITE(*,112)  ERRER,ERGOAL,FEXTRP,IT
112   FORMAT(' ERRER,ERGOAL,FEXTRP,IT',1P,2E15.7,1P,E14.5,0P,I3)
      IF ( L-IT .GT. 2  .AND.  ISTAGE .LT. MXSTGE ) GOTO 89
      IF ( ERRER/FEXTRP**(MAXTBL-L) .LT. ERGOAL ) GOTO 63
      GOTO 89
 80   IF ( LEVEL .GE. 4 ) WRITE(*,113) L,T(1,LM1)
113   FORMAT(1X,I5,E16.8,5X,'JUMP')
      IF ( ERRER .GT. ERGOAL ) GOTO 89
      DIFF=DABS(T(1,L))*2.D0*FN
      IF ( LEVEL .GE. 2 ) WRITE(*,114) DIFF,BEGI,ENDE
114   FORMAT(13X,'INTEGRAND SEEMS TO HAVE JUMP OF SIZE',
     *   1P,E15.7,' BETWEEN',1P,E15.7,' AND',1P,E15.7)
      GOTO 86
 81   IF ( LEVEL .GE. 4 ) WRITE(*,115) L
115   FORMAT(1X,I5,21X,'STRAIGHT LINE')
      SLOPE=(FEND-FBEG)*2.D0
      FBEG2=FBEG*2.D0
      DO  82  I=1,4
        DIFF=DABS(FUNC(BEGI+RN(I)*STEP)-FBEG2-RN(I)*SLOPE)
        IF ( DIFF .GT. TABTLM)  GOTO 85
 82   CONTINUE
      IF ( LEVEL .GE. 3 ) WRITE(*,116) BEGI,ENDE
116   FORMAT(27X,'INTEGRAND SEEMS TO BE STRAIGHT LINE BETWEEN',
     *       1P,E15.7,' AND',1P,E15.7)
      GOTO 86
 83   IF ( LEVEL .GE. 4 ) WRITE(*,117) L,T(1,LM1)
117   FORMAT(1X,I5,1P,E15.7,5X,'NOISE')
      SLOPE=(FEND-FBEG)*2.D0
      FBEG2=FBEG*2.D0
      I=1
 84   DIFF=DABS(FUNC(BEGI+RN(I)*STEP)-FBEG2-RN(I)*SLOPE)
 85   ERRER=DMAX1(ERRER,ASTEP*DIFF)
      IF ( ERRER .GT. ERGOAL ) GOTO 90
      I=I+1
      IF ( I .LE. 4 ) GOTO 84
      IF ( LEVEL .GE. 3 ) WRITE(*,118) BEGI,ENDE
118   FORMAT(' NOISE BETWEEN',1P,E15.7,' AND',1P,E15.7)
      IFLAG=3
 86   VAL=VAL+VINT
      ERROR=ERROR+ERRER
      IF ( LEVEL .GE. 3 ) THEN
        IF ( LEVEL .GE. 5 ) THEN
          DO 87 I=1,L
            IF ( LEVEL .GE. 4 ) WRITE(*,119) I,(T(I,J),J=1,L)
119         FORMAT(1X,I5,7E16.8,/,1X,3E16.8)
 87       CONTINUE
        ENDIF
        WRITE(*,120) VINT,ERRER,L,IT
120     FORMAT(' INTEGRAL IS',1P,E16.8,', ERROR',1P,E16.8,
     *         '  FROM T(',I1,',',I1,'1H)')
      ENDIF
      IF ( RIGHT ) GOTO 88
      ISTAGE=ISTAGE-1
      IF ( ISTAGE .EQ. 0 ) RETURN
      REGLAR=REGLSV(ISTAGE)
      BEGI=BEGIN(ISTAGE)
      ENDE=FINIS(ISTAGE)
      CUREST=CUREST-EST(ISTAGE+1)+VINT
      IEND=IBEG-1
      FEND=TS(IEND)
      IBEG=IBEGS(ISTAGE)
      GOTO 92
 88   CUREST=CUREST+VINT
      STAGE=STAGE*2.D0
      IEND=IBEG
      IBEG=IBEGS(ISTAGE)
      ENDE=BEGI
      BEGI=BEGIN(ISTAGE)
      FEND=FBEG
      FBEG=TS(IBEG)
      GOTO 60
 89   REGLAR=.TRUE.
 90   IF ( ISTAGE .EQ. MXSTGE ) GOTO 97
      IF ( LEVEL .GE. 5 ) THEN
        DO 91 I=1,L
          IF ( LEVEL .GE. 4 ) WRITE(*,121) I,(T(I,J),J=1,L)
121       FORMAT(1X,I5,7E16.8,/,1X,3E16.8)
 91     CONTINUE
      ENDIF
      IF ( RIGHT ) GOTO 93
      REGLSV(ISTAGE+1)=REGLAR
      BEGIN(ISTAGE)=BEGI
      IBEGS(ISTAGE)=IBEG
      STAGE=STAGE/2.D0
 92   RIGHT=.TRUE.
      BEGI=(BEGI+ENDE)/2.D0
      IBEG=(IBEG+IEND)/2
      TS(IBEG)=TS(IBEG)/2.D0
      FBEG=TS(IBEG)
      GOTO 61
 93   NNLEFT=IBEG-IBEGS(ISTAGE)
      IF ( IEND+NNLEFT .GE. MAXTS ) GOTO 96
      III=IBEGS(ISTAGE)
      II=IEND
      DO 94 I=III,IBEG
        II=II+1
        TS(II)=TS(I)
 94   CONTINUE
      DO 95 I=IBEG,II
        TS(III)=TS(I)
        III=III+1
 95   CONTINUE
      IEND=IEND+1
      IBEG=IEND-NNLEFT
      FEND=FBEG
      FBEG=TS(IBEG)
      FINIS(ISTAGE)=ENDE
      ENDE=BEGI
      BEGI=BEGIN(ISTAGE)
      BEGIN(ISTAGE)=ENDE
      REGLSV(ISTAGE)=REGLAR
      ISTAGE=ISTAGE+1
      REGLAR=REGLSV(ISTAGE)
      EST(ISTAGE)=VINT
      CUREST=CUREST+EST(ISTAGE)
      GOTO 60
 96   CONTINUE
      IF ( LEVEL .GE. 2 ) WRITE(*,122) BEGI,ENDE
122   FORMAT(' TOO MANY FUNCTION EVALUATIONS AROUND',/,
     *       11X,1P,E15.7,' AND',1P,E15.7)
      IFLAG=4
      GOTO 99
 97   CONTINUE
      IFLAG=5
      IF ( LEVEL .GE. 2 ) THEN
        IF ( LEVEL .GE. 5 ) THEN
          DO  98  I =1,L
            IF ( LEVEL .GE. 4 ) WRITE(*,123) I,(T(I,J),J=1,L)
123         FORMAT(1X,I5,7E16.8,/,1X,3E16.8)
 98       CONTINUE
        ENDIF
        WRITE(*,124) BEGI,ENDE
124     FORMAT(11X,'INTEGRAND SHOWS SINGULAR BEHAVIOUR OF ',
     *      'UNKNOWN TYPE BETWEEN',1P,E15.7,' AND',1P,E15.7)
      ENDIF
 99   CONTINUE
      VAL=CUREST+VINT
      RETURN
      END
C=======================================================================

      FUNCTION UTACOS(X)

C-----------------------------------------------------------------------
C  RETURNS ACOS(X) FOR -1 <= X <= 1 , ACOS(+-1) ELSE
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      ARGUM=X
      IF     ( X .LT. -1. ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('UTACOS')
          WRITE(IFCH,*)'*****  ARGUM = ',ARGUM,' SET -1'
          CALL UTMSGF
        ENDIF
        ARGUM=-1.
      ELSEIF ( X .GT.  1. ) THEN
        IF ( ISH .GE .90 ) THEN
          CALL UTMSG('UTACOS')
          WRITE(IFCH,*)'*****  ARGUM = ',ARGUM,' SET 1'
          CALL UTMSGF
        ENDIF
        ARGUM=1.
      ENDIF
      UTACOS=ACOS(ARGUM)
      RETURN
      END
C=======================================================================

      FUNCTION UTAMNU(KEUX,KEDX,KESX,KECX,MODUS)

C----------------------------------------------------------------------
C  RETURNS MIN MASS OF CLUSTER WITH GIVEN U,D,S,C CONTENT
C  KEUX: NET U QUARK NUMBER
C  KEDX: NET D QUARK NUMBER
C  KESX: NET S QUARK NUMBER
C  KECX: NET C QUARK NUMBER
C  MODUS:0,1,2,3,4,5,6
C----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CSJCGA/  AMEGAM,AMNULL,ASUHA(7),ENTRPY,NOPHA,NSUHA(7)
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO3/   ASUHAX(7),ASUHAY(7),OMEGA,SIGPPD,SIGPPE,UENTRO
     *                ,IWZZZZ
      SAVE
C----------------------------------------------------------------------
      AMNULL=0.

      IF     ( MODUS .EQ. 0 ) THEN
        DO 21 I=1,7
          ASUHA(I)=ASUHAY(I)
          NSUHA(I)=0
21      CONTINUE
      ELSEIF ( MODUS .EQ. 1 ) THEN
        IF     ( NOPHA .EQ. 0 ) THEN
          DO 22 I=1,7
            ASUHA(I)=ASUHAX(I)
            NSUHA(I)=0
22        CONTINUE
        ELSEIF ( NOPHA .GE. 1 ) THEN
          DO 23 I=1,7
            ASUHA(I)=ASUHAY(I)
            NSUHA(I)=0
23        CONTINUE
        ENDIF
      ELSEIF ( MODUS .EQ. 2 ) THEN
        DO 24 I=1,7
          ASUHA(I)=ASUHAY(I)
          NSUHA(I)=0
24      CONTINUE
      ELSEIF ( MODUS .EQ. 3 ) THEN
        DO 25 I=1,7
          ASUHA(I)=ASUHAY(I)
          NSUHA(I)=0
25      CONTINUE
      ELSEIF ( MODUS .EQ. 4 ) THEN
        DO 26 I=1,7
          ASUHA(I)=ASUHAX(I)
          NSUHA(I)=0
26      CONTINUE
      ELSEIF ( MODUS .EQ. 5 ) THEN
        DO 27 I=1,7
          ASUHA(I)=ASUHAY(I)
          NSUHA(I)=0
27      CONTINUE
      ELSEIF ( MODUS .EQ. 6 ) THEN
        DO 28 I=1,7
          ASUHA(I)=SQRT(ASUHAY(I)**2+DELMSS**2)
          NSUHA(I)=0
28      CONTINUE
      ENDIF
      IAUX=KEUX+KEDX+KESX+KECX
      KE=ABS(IAUX)

      IF ( IAUX .GE. 0 ) THEN
        KEU=KEUX
        KED=KEDX
        KES=KESX
        KEC=KECX
      ELSE
        KEU=-KEUX
        KED=-KEDX
        KES=-KESX
        KEC=-KECX
      ENDIF

      IF ( KEC .NE. 0 ) THEN
10      CONTINUE
        IF   ( KEC .LT. 0 ) THEN
          KEC=KEC+1
          IF ( KEU .GT. KED ) THEN
            KEU=KEU-1
          ELSE
            KED=KED-1
          ENDIF
          AMNULL=AMNULL+1.87
          GOTO 10
        ENDIF
11      CONTINUE
        IF ( KEC .GT. 0 ) THEN
          KEC=KEC-1
          IF ( KEU .LT. KED ) THEN
            KEU=KEU+1
          ELSE
            KED=KED+1
          ENDIF
          AMNULL=AMNULL+1.87
          GOTO 11
        ENDIF
      ENDIF

 5    CONTINUE
      IF ( KES .LT. 0 ) THEN
        AMNULL=AMNULL+ASUHA(6)
        IF ( KEU .GE. KED ) THEN
          KEU=KEU-1
        ELSE
          KED=KED-1
        ENDIF
        KES=KES+1
        GOTO 5
      ENDIF

 6    CONTINUE
      IF ( KED .LT. 0 ) THEN
        IF ( KEU .GE. KES ) THEN
          AMNULL=AMNULL+ASUHA(5)
          KEU=KEU-1
        ELSE
          AMNULL=AMNULL+ASUHA(6)
          KES=KES-1
        ENDIF
        KED=KED+1
        GOTO 6
      ENDIF

 7    CONTINUE
      IF ( KEU .LT. 0 ) THEN
        IF ( KED .GE. KES ) THEN
          AMNULL=AMNULL+ASUHA(5)
          KED=KED-1
        ELSE
          AMNULL=AMNULL+ASUHA(6)
          KES=KES-1
        ENDIF
        KEU=KEU+1
        GOTO 7
      ENDIF

      IF ( KEU+KED+KES+KEC .NE. KE ) THEN
        CALL UTSTOP('UTAMNU: SUM_KEI /= KE                   ')
      ENDIF
      KEQ=KEU+KED
      KEQX=KEQ

      IF ( MODUS .EQ. 2  .AND.  KE .GT. 3 ) THEN
        CALL SMASST(KEU,KED,KES,KEC,A,Y,Z)
        AMNUZ=SMASS(A,Y,Z)
      ENDIF

      AMNUX=0.

      I=4
 2    I=I-1
 3    CONTINUE
      IF ( (4-I)*KES .GT. (I-1)*KEQ ) THEN
        AMNUX=AMNUX+ASUHA(1+I)
        KEQ=KEQ-3+I
        KES=KES-I
        IF ( KES .LT. 0 ) THEN
          CALL UTSTOP('UTAMNU: NEGATIVE KES                    ')
        ENDIF
        IF ( KEQ .LT. 0 ) THEN
          CALL UTSTOP('UTAMNU: NEGATIVE KEQ                    ')
        ENDIF
        GOTO 3
      ENDIF
      IF ( I .GT. 1 ) GOTO 2

      IF ( KEQX .GT. KEQ ) THEN
        DO 8 K=1,KEQX-KEQ
          IF ( KEU .GE. KED ) THEN
            KEU=KEU-1
          ELSE
            KED=KED-1
          ENDIF
 8      CONTINUE
      ENDIF
      IF ( KEU+KED .NE. KEQ ) THEN
        CALL UTSTOP('UTAMNU: KEU+KED /= KEQ                  ')
      ENDIF

 9    CONTINUE
      IF ( KEU .GT. 2*KED ) THEN
        AMNUX=AMNUX+ASUHA(7)
        KEU=KEU-3
        IF ( KEU .LT. 0 ) THEN
          CALL UTSTOP('UTAMNU: NEGATIVE KEU                    ')
        ENDIF
        GOTO 9
      ENDIF
      IF ( KED .GT. 2*KEU ) THEN
        AMNUX=AMNUX+ASUHA(7)
        KED=KED-3
        IF ( KED .LT. 0 ) THEN
          CALL UTSTOP('UTAMNU: NEGATIVE KED                    ')
        ENDIF
        GOTO 9
      ENDIF

      KEQ=KEU+KED
      IF ( MOD(KEQ,3) .NE. 0 ) THEN
        CALL UTSTOP('UTAMNU: MOD(KEQ,3) /= 0                 ')
      ENDIF
      AMNUX=AMNUX+ASUHA(1)*KEQ/3

      AMNU=AMNUX
      IF ( MODUS.EQ.2 .AND. KE.GT.3 ) AMNU=MIN(AMNUX,AMNUZ)
      AMNULL=AMNULL+AMNU

      IF ( AMNULL.EQ.0. .AND. MODUS.GT.0 ) AMNULL=ASUHA(5)

ctp060203 1000  UTAMNU=AMNULL
      UTAMNU=AMNULL
      RETURN
      END
C=======================================================================

      FUNCTION UTAMNX(JCP,JCM)

C-----------------------------------------------------------------------
C  RETURNS MINIMUM MASS FOR THE DECAY OF JCP---JCM (BY CALLING UTAMNU).
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER JCM(NFLAV,2),JCP(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      DO 3 I=1,NFLAV
        IF ( JCP(I,1) .NE. 0 ) GOTO 1
        IF ( JCP(I,2) .NE. 0 ) GOTO 1
 3    CONTINUE
      KEU=JCM(1,1)-JCM(1,2)
      KED=JCM(2,1)-JCM(2,2)
      KES=JCM(3,1)-JCM(3,2)
      KEC=JCM(4,1)-JCM(4,2)
      UTAMNX=UTAMNU(KEU,KED,KES,KEC,5)
      RETURN
 1    CONTINUE

      DO  4 I=1,NFLAV
        IF ( JCM(I,1) .NE. 0 ) GOTO 2
        IF ( JCM(I,2) .NE. 0 ) GOTO 2
 4    CONTINUE
      KEU=JCP(1,1)-JCP(1,2)
      KED=JCP(2,1)-JCP(2,2)
      KES=JCP(3,1)-JCP(3,2)
      KEC=JCP(4,1)-JCP(4,2)
      UTAMNX=UTAMNU(KEU,KED,KES,KEC,5)
      RETURN
 2    CONTINUE

      KEU=JCP(1,1)-JCP(1,2)
      KED=JCP(2,1)-JCP(2,2)
      KES=JCP(3,1)-JCP(3,2)
      KEC=JCP(4,1)-JCP(4,2)
      KE=KEU+KED+KES+KEC
      IF     ( MOD(KE+1,3) .EQ. 0 ) THEN
        KEU=KEU+1
      ELSEIF ( MOD(KE-1,3) .EQ. 0 ) THEN
        KEU=KEU-1
      ELSE
        CALL UTSTOP('UTAMNX: NO SINGLET POSSIBLE (1)         ')
      ENDIF
      AMMS=UTAMNU(KEU,KED,KES,KEC,5)
      KEU=JCM(1,1)-JCM(1,2)
      KED=JCM(2,1)-JCM(2,2)
      KES=JCM(3,1)-JCM(3,2)
      KEC=JCM(4,1)-JCM(4,2)
      KE=KEU+KED+KES+KEC
      IF     ( MOD(KE+1,3) .EQ. 0 ) THEN
        KEU=KEU+1
      ELSEIF ( MOD(KE-1,3) .EQ. 0 ) THEN
        KEU=KEU-1
      ELSE
        CALL UTSTOP('UTAMNX: NO SINGLET POSSIBLE (2)         ')
      ENDIF
      UTAMNX=AMMS+UTAMNU(KEU,KED,KES,KEC,5)
      RETURN
      END
C=======================================================================

      FUNCTION UTAMNY(JCP,JCM)

C-----------------------------------------------------------------------
C  RETURNS MINIMUM MASS OF JCP+JCM (BY CALLING UTAMNU).
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER JC(NFLAV,2),JCM(NFLAV,2),JCP(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      DO 7 NF=1,NFLAV
        JC(NF,1)=JCP(NF,1)+JCM(NF,1)
        JC(NF,2)=JCP(NF,2)+JCM(NF,2)
 7    CONTINUE
      KEU=JC(1,1)-JC(1,2)
      KED=JC(2,1)-JC(2,2)
      KES=JC(3,1)-JC(3,2)
      KEC=JC(4,1)-JC(4,2)
      UTAMNY=UTAMNU(KEU,KED,KES,KEC,5)
      RETURN
      END
C=======================================================================

      FUNCTION UTAMNZ(JC,MODUS)

C-----------------------------------------------------------------------
C  RETURNS MINIMUM MASS OF JC (BY CALLING UTAMNU).
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER JC(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      KEU=JC(1,1)-JC(1,2)
      KED=JC(2,1)-JC(2,2)
      KES=JC(3,1)-JC(3,2)
      KEC=JC(4,1)-JC(4,2)
      UTAMNZ=UTAMNU(KEU,KED,KES,KEC,MODUS)
      RETURN
      END
C=======================================================================

      SUBROUTINE UTAMST(STS,AM,AMIN,IRET)

C-----------------------------------------------------------------------
C  INPUT:  STS   = STRING  (SINGLE)
C  OUTPUT: AM    = MASS
C          AMIN  = MINIMUM MASS
C          IRET  = RETURN CODE (=3 IF AM**2.LT.AMIN**2, 0 ELSE)
C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      PARAMETER (NSI=6)
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU

      REAL    STS(NSI,2)
      INTEGER IC1(2),IC2(2),JC(NFLAV,2),JC1(NFLAV,2),JC2(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      IRET=0
      AM2= (STS(4,1)+STS(4,2))**2 - (STS(3,1)+STS(3,2))**2 -
     *     (STS(2,1)+STS(2,2))**2 - (STS(1,1)+STS(1,2))**2
      AM=-AM2
      IC1(1)=NINT(STS(4+1,1))
      IC1(2)=NINT(STS(4+2,1))
      IC2(1)=NINT(STS(4+1,2))
      IC2(2)=NINT(STS(4+2,2))
      CALL IDDECO(IC1,JC1)
      CALL IDDECO(IC2,JC2)
      DO 2 NF=1,NFLAV
        JC(NF,1)=JC1(NF,1)+JC2(NF,1)
        JC(NF,2)=JC1(NF,2)+JC2(NF,2)
        IF ( NF.GT.4 .AND. (JC(NF,1).NE.0 .OR. JC(NF,2).NE.0) ) THEN
          CALL UTSTOP('UTAMST: FLAVOUR > 4                     ')
        ENDIF
 2    CONTINUE
      KEU=JC(1,1)-JC(1,2)
      KED=JC(2,1)-JC(2,2)
      KES=JC(3,1)-JC(3,2)
      KEC=JC(4,1)-JC(4,2)
      AMIN=UTAMNU(KEU,KED,KES,KEC,4)
      IF ( AM2. LT. AMIN**2 ) THEN
        IRET=3
        RETURN
      ELSE
        AM=SQRT(AM2)
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTAXIS(I,J,A1,A2,A3)

C---------------------------------------------------------------------
C  CALCULATES THE AXIS DEFINED BY THE PTLS I,J IN THE I,J CM SYSTEM
C---------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C---------------------------------------------------------------------
      A1=0.
      A2=0.
      A3=1.
      PI1=PPTL(1,I)
      PI2=PPTL(2,I)
      PI3=PPTL(3,I)
      PI4=PPTL(4,I)
      PJ1=PPTL(1,J)
      PJ2=PPTL(2,J)
      PJ3=PPTL(3,J)
      PJ4=PPTL(4,J)
      P1=PI1+PJ1
      P2=PI2+PJ2
      P3=PI3+PJ3
      P4=PI4+PJ4
      P5=SQRT(P4**2-P3**2-P2**2-P1**2)
      CALL UTLOBO(1,P1,P2,P3,P4,P5,PI1,PI2,PI3,PI4)
      CALL UTLOBO(1,P1,P2,P3,P4,P5,PJ1,PJ2,PJ3,PJ4)
      ERR=(PI1+PJ1)**2+(PI2+PJ2)**2+(PI3+PJ3)**2
      IF ( ISH.GE.90 .AND. ERR .GT. 1.E-3 ) THEN
        CALL UTMSG('UTAXIS')
        WRITE(IFCH,*)'*****  ERR=',ERR
        WRITE(IFCH,*)'PI:',PI1,PI2,PI3,PI4
        WRITE(IFCH,*)'PJ:',PJ1,PJ2,PJ3,PJ4
        CALL UTMSGF
      ENDIF
      A=SQRT( (PJ1-PI1)**2 + (PJ2-PI2)**2 + (PJ3-PI3)**2 )
      IF ( A .EQ. 0. ) RETURN
      A1=(PI1-PJ1)/A
      A2=(PI2-PJ2)/A
      A3=(PI3-PJ3)/A
      RETURN
      END
C=======================================================================

      SUBROUTINE UTCHM(ARP,ARM,II)

C-----------------------------------------------------------------------
C  CHECKS WHETHER ARP**2=0. AND ARM**2=0.
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      DOUBLE PRECISION ARM(4),ARP(4),DIFM,DIFP
      SAVE
C-----------------------------------------------------------------------
      IF(ISH.LT.90)RETURN
      DIFP=ARP(4)**2-ARP(1)**2-ARP(2)**2-ARP(3)**2
      DIFM=ARM(4)**2-ARM(1)**2-ARM(2)**2-ARM(3)**2
      IF ( ABS(DIFP) .GT. 1.D-3*ARP(4)**2   .OR.
     *     ABS(DIFM) .GT. 1.D-3*ARM(4)**2 ) THEN
        CALL UTMSG('UTCHM ')
        WRITE(IFCH,*)'*****  MASS NON ZERO  -  ',II
        WRITE(IFCH,*)'JET-MASS**2`S:    ',DIFP,DIFM
        WRITE(IFCH,*)'ENERGY**2`S:      ',ARP(4)**2,ARM(4)**2
        WRITE(IFCH,*)(SNGL(ARP(I)),I=1,4)
        WRITE(IFCH,*)(SNGL(ARM(I)),I=1,4)
        CALL UTMSGF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTCLEA(NPTL0)

C-----------------------------------------------------------------------
C  OVERWRITES ISTPTL=2 PARTICLES IN /CPTL/, REDUCES SO NPTL.
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      INTEGER NEWPTL(MXPTL)
      SAVE
C-----------------------------------------------------------------------
      ISH0=ISH
      IF ( ISHSUB/100 .EQ. 18 ) ISH=MOD(ISHSUB,100)

      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)('-',L=1,68)
        WRITE(IFCH,*)'SR UTCLEA. INITIAL.'
        WRITE(IFCH,*)('-',L=1,68)
        DO 34 N=1,NPTL
          WRITE(IFCH,116)IORPTL(N),JORPTL(N),N,IFRPTL(1,N),IFRPTL(2,N)
     *     ,IDPTL(N),SQRT(PPTL(1,N)**2+PPTL(2,N)**2),PPTL(3,N),PPTL(5,N)
     *     ,ISTPTL(N)
34      CONTINUE
116     FORMAT(1X,I6,I6,2X,I6,2X,I4,I4,I12,3(E10.2),I3,I3)
      ENDIF

      I=0
 1    I=I+1
      IF ( I .GT. NPTL ) GOTO 1000
      IF ( ISTPTL(I) .EQ. 2 ) GOTO 2
      NEWPTL(I)=I
      GOTO 1

 2    I=I-1
      J=I
 3    I=I+1
 4    J=J+1
      IF ( J .GT. NPTL ) GOTO 5
      NEWPTL(J)=0
      IF ( ISTPTL(J) .EQ. 2 ) GOTO 4
      NEWPTL(J)=I
      CALL UTREPL(I,J)
      GOTO 3

 5    NPTL=I-1
      IF ( NPTL .EQ. 0 ) THEN
        NPTL0=0
        GOTO 1000
      ENDIF

20    N0=NEWPTL(NPTL0)
      IF ( N0 .GT. 0 ) THEN
        NPTL0=N0
      ELSE
        NPTL0=NPTL0-1
        IF ( NPTL0 .GT. 0 ) GOTO 20
      ENDIF

      DO 11 K=1,NPTL
        IO=IORPTL(K)
        IF ( IO .LE. 0 ) THEN
          IORPTL(K)=IO
        ELSE
          IORPTL(K)=NEWPTL(IO)
        ENDIF

        JO=JORPTL(K)
        IF ( JO .LE. 0 ) THEN
          JORPTL(K)=JO
        ELSE
          JORPTL(K)=NEWPTL(JO)
        ENDIF

        IF1=IFRPTL(1,K)
        IF ( IF1 .LE. 0 ) THEN
          IFRPTL(1,K)=IF1
        ELSE
          IFRPTL(1,K)=NEWPTL(IF1)
        ENDIF

        IF2=IFRPTL(2,K)
        IF ( IF2 .LE. 0 ) THEN
          IFRPTL(2,K)=IF2
        ELSE
          IFRPTL(2,K)=NEWPTL(IF2)
        ENDIF
11    CONTINUE

      DO 19 K=1,NPTL
        IF ( IFRPTL(1,K).EQ.0 .AND. IFRPTL(2,K).GT.0 )
     *                                         IFRPTL(1,K)=IFRPTL(2,K)
        IF ( IFRPTL(1,K).GT.0 .AND. IFRPTL(2,K).EQ.0 )
     *                                         IFRPTL(2,K)=IFRPTL(1,K)
19    CONTINUE

1000  CONTINUE

      IF ( ISH .GE. 92 ) THEN
        WRITE(IFCH,*)('-',L=1,68)
        WRITE(IFCH,*)'SR UTCLEA. FINAL.'
        WRITE(IFCH,*)('-',L=1,68)
        DO 35 N=1,NPTL
          WRITE(IFCH,116)IORPTL(N),JORPTL(N),N,IFRPTL(1,N),IFRPTL(2,N)
     *     ,IDPTL(N),SQRT(PPTL(1,N)**2+PPTL(2,N)**2),PPTL(3,N),PPTL(5,N)
     *     ,ISTPTL(N)
35      CONTINUE
        WRITE(IFCH,*)('-',L=1,79)
      ENDIF

      ISH=ISH0
      RETURN
      END
C=======================================================================

      SUBROUTINE UTHIST(X1,X2,Y1,Y2,N,X,Y,LITY,LILO,TEXT1,TEXT2,TEXT3)

C----------------------------------------------------------------------
C  WRITES ARRAYS X,Y IN HISTO-FORMAT
C----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      REAL      X(N),Y(N)
      CHARACTER LILO*6,LITY*3,TEXT1*50,TEXT2*50,TEXT3*50
      SAVE
C----------------------------------------------------------------------
      WRITE(IFCH,100)('-',I=1,69)
      WRITE(IFCH,100)('-',I=1,69)
100   FORMAT(1X,1H',69A1,1H')
      WRITE(IFCH,101)X1,X2,0.,Y1,Y2,0.,'_____',0
     *,'___',LITY,'___','_______',LILO,0.,3,N,0.
101   FORMAT(1X,6F8.2,2X,1H',A5,1H',I11
     *,/,4X,1H',A3,1H',3X,1H',A3,A3,1H',3X,1H',A7,A6,1H',F8.2,2I6,F10.4)
      WRITE(IFCH,100)('-',I=1,69)
      WRITE(IFCH,102)0.,0.,TEXT1
      WRITE(IFCH,102)0.,0.,TEXT2
      WRITE(IFCH,102)0.,0.,TEXT3
102   FORMAT(4X,F4.2,2X,F4.2,2X,1H',A50,1H')
      WRITE(IFCH,100)('-',I=1,69)
      DO 36 K=1,N
        WRITE(IFCH,103)K,X(K),0.,Y(K),0.,0.
103     FORMAT(4X,I5,2X,5E12.4)
36    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE UTHSEA

C-----------------------------------------------------------------------
C  CREATES HISTOGRAM OF SEA STRUCTURE FUNCTION
C-----------------------------------------------------------------------
      PARAMETER (NSTRU=2049)
      COMMON /STRU/    QSEP(NSTRU),QSET(NSTRU),QVAP(NSTRU)
     *                ,QVAT(NSTRU),XCUTAR,XSTRU(NSTRU)
     *                ,IDTG
      REAL XAR(1000),YAR(1000)
      SAVE
C-----------------------------------------------------------------------
      DX=1.
      X1=-1.
      DO 3 K=1,2
        DX=DX*0.001
        X1=X1-3.
        DO 2 N=1,1000
          YAR(N)=0.
          XAR(N)=-DX*0.5+N*DX
 2      CONTINUE
        DO 1 I=1,100000
          X=UTINVT(NSTRU,XSTRU,QSEP,RANGEN()*QSEP(NSTRU))
          N=1+X/DX
          IF ( N .LE. 1000 ) YAR(N)=YAR(N)+1.
 1      CONTINUE
        CALL UTHIST(X1,X1+4,0.,5.,1000,XAR,YAR,'POC','LOGLOG'
     *        ,'XAXIS X                         $                 '
     *        ,'YAXIS COUNTS                    $                 '
     *        ,'TITLE SEA QUARK STRUCTURE FUNCTION                ')
 3    CONTINUE
      RETURN
      END
C=======================================================================

      FUNCTION UTINVT(N,X,Q,YY)

C-----------------------------------------------------------------------
C  RETURNS X WITH Y=Q(X)
C-----------------------------------------------------------------------
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CUTINV/  LUTINV
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      REAL Q(N),X(N)
      SAVE
C-----------------------------------------------------------------------
      IF ( Q(N) .EQ. 0. ) THEN
        CALL UTSTOP('UTINVT: Q(N)=0. DIMENSIONS TOO BIG      ')
      ENDIF
      Y = YY
      IF     ( Y .LT. 0. ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('UTINVT')
          WRITE(IFCH,*)'*****  Y=',Y,' < 0'
          CALL UTMSGF
        ENDIF
        Y=0.
      ELSEIF ( Y .GT. Q(N) ) THEN
        IF ( ISH .GE. 90 ) THEN
          CALL UTMSG('UTINVT')
          WRITE(IFCH,*)'*****  Y=',Y,' > ',Q(N)
          CALL UTMSGF
        ENDIF
        Y=Q(N)
      ENDIF
      LU=1
      LO=N
 1    LZ=(LO+LU)/2
      IF     ( Q(LU).LE.Y .AND. Y.LE.Q(LZ) ) THEN
        LO=LZ
      ELSEIF ( Q(LZ).LT.Y .AND. Y.LE.Q(LO) ) THEN
        LU=LZ
      ELSE
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'Q(1),Y,Q(N):',Q(1),Y,Q(N)
        WRITE(IFCH,*)'LU,LZ,LO:',LU,LZ,LO
        WRITE(IFCH,*)'Q(LU),Q(LZ),Q(LO):',Q(LU),Q(LZ),Q(LO)
        CALL UTSTOP('UTINVT: NO INTERVAL FOUND               ')
      ENDIF
      IF ( LO-LU .GE. 2 ) GOTO 1
      IF ( LO .LE. LU ) THEN
        CALL UTSTOP('UTINVT: LO.LE.LU                        ')
      ENDIF
      UTINVT=X(LU)+(Y-Q(LU))*(X(LO)-X(LU))/(Q(LO)-Q(LU))
      LUTINV=LU
      RETURN
      END
C=======================================================================

      SUBROUTINE UTKSIX(SIX,KMAX)

C-----------------------------------------------------------------------
C  RETURNS KMAX FOR SIX
C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      REAL SIX(NSI,NSIX)
      SAVE
C-----------------------------------------------------------------------
      DO 10 K=0,NSIX-1
        KMAX=K
        SIXSQR=0.
        DO 20 I=1,NSI
          SIXSQR=SIXSQR+SIX(I,K+1)**2
20      CONTINUE
        IF ( SIXSQR .LE. 1.E-5 ) RETURN
10    CONTINUE
      CALL UTSTOP('UTKSIX: DIMENSION NSIX TOO SMALL        ')
      RETURN
      END
C=======================================================================

      SUBROUTINE UTKSTR(STR,KMAX)

C-----------------------------------------------------------------------
C  RETURNS KMAX FOR STR
C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      REAL STR(NSI,NSIX+1)
      SAVE
C-----------------------------------------------------------------------
      DO 10 K=0,NSIX
        KMAX=K
        STRSQR=0.
        DO 20 I=1,NSI
          STRSQR=STRSQR+STR(I,K+1)**2
20      CONTINUE
        IF ( STRSQR .LE. 1.E-5 ) RETURN
10    CONTINUE
      CALL UTSTOP('UTKSTR: DIMENSION NSIX TOO SMALL        ')
      RETURN
      END
C=======================================================================

      SUBROUTINE UTLOBO(ISIG,P1,P2,P3,P4,P5,X1,X2,X3,X4)

C-----------------------------------------------------------------------
C  PERFORMS A LORENTZ BOOST
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      REAL BETA(4),Z(4)
      SAVE
C-----------------------------------------------------------------------
      IF ( P5 .LE. 0. ) THEN
        CALL UTMSG('UTLOBO')
        WRITE(IFCH,*)'*****  MASS <= 0.'
        WRITE(IFCH,*)'P(5): ',P1,P2,P3,P4,P5
        CALL UTMSGF
        CALL UTSTOP('UTLOBO: MASS <= 0.                      ')
      ENDIF
      Z(1)=X1
      Z(2)=X2
      Z(3)=X3
      Z(4)=X4
      BETA(1)=-P1/P5
      BETA(2)=-P2/P5
      BETA(3)=-P3/P5
      BETA(4)= P4/P5
      BP=ISIG*(Z(1)*BETA(1)+Z(2)*BETA(2)+Z(3)*BETA(3))
      AUXIL= ISIG*Z(4)+ISIG*BP/(BETA(4)+1.)
      Z(1)=Z(1)+BETA(1)*AUXIL
      Z(2)=Z(2)+BETA(2)*AUXIL
      Z(3)=Z(3)+BETA(3)*AUXIL
      Z(4)=BETA(4)*Z(4)+BP
      X1=Z(1)
      X2=Z(2)
      X3=Z(3)
      X4=Z(4)
      RETURN
      END
C=======================================================================

      SUBROUTINE UTLOB2(ISIG,P1,P2,P3,P4,P5,X1,X2,X3,X4)

C-----------------------------------------------------------------------
C  PERFORMS A LORENTZ BOOST, DOUBLE PREC.
C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION BETA(4),BP,DAUXIL,PP,P1,P2,P3,P4,P5,P5I
     *                ,XX0,X1,X10,X2,X20,X3,X30,X4,X4X,X40,Z(4)
      SAVE
C-----------------------------------------------------------------------
      IF ( ISH .GE. 90 ) THEN
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,101)SNGL(X1),SNGL(X2),SNGL(X3),SNGL(X4)
     *                  ,SNGL(X4**2-X3**2-X2**2-X1**2)
101       FORMAT(' UTLOB2:',F9.5,4F13.5)
        ENDIF
        PP=P4**2-P3**2-P2**2-P1**2
        IF ( ABS(PP-P5**2) .GT. 1.D-3*P4**2   .AND.
     *       ABS(PP-P5**2) .GT. 1.D-3        ) THEN
          CALL UTMSG('UTLOB2')
          WRITE(IFCH,*)'*****  P**2 .NE. P5**2'
          WRITE(IFCH,*)'P**2,P5**2: ',PP,P5**2
          WRITE(IFCH,*)'P: ',P1,P2,P3,P4,P5
          CALL UTMSGF
        ENDIF
        X10=X1
        X20=X2
        X30=X3
        X40=X4
      ENDIF
      XX0=X4**2-X3**2-X2**2-X1**2
      IF ( P5 .LE. 0.D0 ) THEN
        CALL UTMSG('UTLOB2')
        WRITE(IFCH,*)'*****  P5 NEGATIVE.'
        WRITE(IFCH,*)'P(5): ',P1,P2,P3,P4,P5
        CALL UTMSGF
        CALL UTSTOP('UTLOB2: P5 NEGATIVE.                    ')
      ENDIF
      Z(1)=X1
      Z(2)=X2
      Z(3)=X3
      Z(4)=X4
      P5I=1.D0/P5
      BETA(4)= P5I*P4
      BETA(1)=-P5I*P1
      BETA(2)=-P5I*P2
      BETA(3)=-P5I*P3
      BP=ISIG*( BETA(1)*Z(1)+BETA(2)*Z(2)+BETA(3)*Z(3) )
      DAUXIL = ISIG*( Z(4) + BP/(BETA(4)+1.D0) )
      Z(1)=Z(1)+BETA(1)*DAUXIL
      Z(2)=Z(2)+BETA(2)*DAUXIL
      Z(3)=Z(3)+BETA(3)*DAUXIL
      Z(4)=BETA(4)*Z(4)+BP
      X1=Z(1)
      X2=Z(2)
      X3=Z(3)
      X4=Z(4)
      IF ( ISH .GE. 93 )
     *        WRITE(IFCH,101)SNGL(X1),SNGL(X2),SNGL(X3),SNGL(X4)
     *                      ,SNGL(X4**2-X3**2-X2**2-X1**2)
      X4X=X4
      X4=SQRT(XX0+X1**2+X2**2+X3**2)
      IF ( ISH .GE. 90 ) THEN
        IF ( ISH .GE. 93 ) THEN
          WRITE(IFCH,101)SNGL(X1),SNGL(X2),SNGL(X3),SNGL(X4)
     *                ,SNGL(X4**2-X3**2-X2**2-X1**2)
          WRITE(IFCH,*)' '
        ENDIF
        IF ( ABS(X4-X4X) .GT. 1.D-2*ABS(X4)    .AND.
     *       ABS(X4-X4X) .GT. 1.D-2        )   THEN
          CALL UTMSG('UTLOB2')
          WRITE(IFCH,*)'*****  X**2_INI .NE. X**2_FIN.'
          WRITE(IFCH,*)'X1 X2 X3 X4 X**2 (INITIAL/FINAL/CORRECTED):'
          WRITE(IFCH,101)SNGL(X10),SNGL(X20),SNGL(X30),SNGL(X40)
     *                  ,SNGL(X40**2-X30**2-X20**2-X10**2)
          WRITE(IFCH,101)SNGL(X1),SNGL(X2),SNGL(X3),SNGL(X4X)
     *                  ,SNGL(X4X**2-X3**2-X2**2-X1**2)
          WRITE(IFCH,101)SNGL(X1),SNGL(X2),SNGL(X3),SNGL(X4)
     *                  ,SNGL(X4**2-X3**2-X2**2-X1**2)
          CALL UTMSGF
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTLOC(AR,N,A,L)

C-----------------------------------------------------------------------
      REAL AR(N)
      SAVE
C-----------------------------------------------------------------------
      DO 1 I=1,N
        IF ( A .LT. AR(I) ) THEN
          L=I-1
          RETURN
        ENDIF
 1    CONTINUE
      L=N
      RETURN
      END
C=======================================================================

      SUBROUTINE UTLOW(CONE)

C-----------------------------------------------------------------------
C  CONVERTS LOWER CASE CHARACTERS TO UPPER CASE CHARACTERS
C-----------------------------------------------------------------------
      CHARACTER*1 CONE
      SAVE
C-----------------------------------------------------------------------
      IF     ( CONE .EQ. 'a' ) THEN
        CONE='A'
      ELSEIF ( CONE .EQ. 'b' ) THEN
        CONE='B'
      ELSEIF ( CONE .EQ. 'c' ) THEN
        CONE='C'
      ELSEIF ( CONE .EQ. 'd' ) THEN
        CONE='D'
      ELSEIF ( CONE .EQ. 'e' ) THEN
        CONE='E'
      ELSEIF ( CONE .EQ. 'f' ) THEN
        CONE='F'
      ELSEIF ( CONE .EQ. 'g' ) THEN
        CONE='G'
      ELSEIF ( CONE .EQ. 'h' ) THEN
        CONE='H'
      ELSEIF ( CONE .EQ. 'i' ) THEN
        CONE='I'
      ELSEIF ( CONE .EQ. 'j' ) THEN
        CONE='J'
      ELSEIF ( CONE .EQ. 'k' ) THEN
        CONE='K'
      ELSEIF ( CONE .EQ. 'l' ) THEN
        CONE='L'
      ELSEIF ( CONE .EQ. 'm' ) THEN
        CONE='M'
      ELSEIF ( CONE .EQ. 'n' ) THEN
        CONE='N'
      ELSEIF ( CONE .EQ. 'o' ) THEN
        CONE='O'
      ELSEIF ( CONE .EQ. 'p' ) THEN
        CONE='P'
      ELSEIF ( CONE .EQ. 'q' ) THEN
        CONE='Q'
      ELSEIF ( CONE .EQ. 'r' ) THEN
        CONE='R'
      ELSEIF ( CONE .EQ. 's' ) THEN
        CONE='S'
      ELSEIF ( CONE .EQ. 't' ) THEN
        CONE='T'
      ELSEIF ( CONE .EQ. 'u' ) THEN
        CONE='U'
      ELSEIF ( CONE .EQ. 'v' ) THEN
        CONE='V'
      ELSEIF ( CONE .EQ. 'w' ) THEN
        CONE='W'
      ELSEIF ( CONE .EQ. 'x' ) THEN
        CONE='X'
      ELSEIF ( CONE .EQ. 'y' ) THEN
        CONE='Y'
      ELSEIF ( CONE .EQ. 'z' ) THEN
        CONE='Z'
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTLOW6(CSIX)

C-----------------------------------------------------------------------
      CHARACTER CSIX*6
      SAVE
C-----------------------------------------------------------------------
      DO 1 I=1,6
        CALL UTLOW(CSIX(I:I))
 1    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE UTMSG(TXT)

C-----------------------------------------------------------------------
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      DOUBLE PRECISION SEEDC,SEEDI
      COMMON /CSEED/   SEEDC,SEEDI
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      CHARACTER*6      TXT
      SAVE
C-----------------------------------------------------------------------
      IMSG=IMSG+1
      IF ( ISH.NE.90 .OR. ISHSUB.NE.0 ) WRITE(IFCH,*)' '
      WRITE(IFCH,*)('-',J=1,77)
      WRITE(IFCH,100)TXT,NREVT+1,IPAGE,SEEDC
100   FORMAT(1X,'***** MSG FROM ',A6,'.   EPS:',I7,I5,2X,D23.17)
      RETURN
      END
C=======================================================================

      SUBROUTINE UTMSGF

C-----------------------------------------------------------------------
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      IF ( ISH.EQ.90 .AND. ISHSUB.EQ.0 ) RETURN
      WRITE(IFCH,*)('-',J=1,77)
      WRITE(IFCH,*)' '
      RETURN
      END
C=======================================================================

      SUBROUTINE UTOVEL

C----------------------------------------------------------------------
C  FILLS ARRAY OVEL(1+I,1+J) CONTAINING THE LOGARITHM (LN) OF
C    I_OVER_J:
C  OVEL(1+I,1+J)=LOG(I!/J!/(J-I)!)      I>=0 J>=0
C----------------------------------------------------------------------
      PARAMETER (IOVMAX=100)
      PARAMETER (JOVMAX=100)
      COMMON /COVEL/   OVEL(1+IOVMAX,1+JOVMAX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C----------------------------------------------------------------------
      DO 1 I=0,IOVMAX
        OVEL(1+I,1)=0.
        OVEL(1+I,1+I)=0.
 1    CONTINUE
      DO 2 J=1,JOVMAX-1
        DO 3 I=J+1,IOVMAX
          OVEL(1+I,1+J)=OVEL(1+I,J)+LOG((I-J+1.)/J)
 3      CONTINUE
 2    CONTINUE

      IF ( ISH .GE. 90 ) THEN
        DO 5 J=1,49
          DO 5 I=J+1,50
            RELERR=ABS((EXP(OVEL(1+I,1+J))
     *           -EXP(OVEL(I,J))-EXP(OVEL(I,1+J)))/EXP(OVEL(1+I,1+J)))
            IF ( RELERR .GT. 1.E-4 ) THEN
              CALL UTMSG('UTOVEL')
              WRITE(IFCH,*)'*****  OVEL(,) VIOLATES RECURRENCE RELATION'
              WRITE(IFCH,*)EXP(OVEL(1+I,1+J))
     *                    ,EXP(OVEL(I,J))+EXP(OVEL(I,1+J))
              CALL UTMSGF
            ENDIF
 5      CONTINUE

        IF ( ISH .GE. 95 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'   EXP( OVEL(1+I,1+J) )'
          WRITE(IFCH,*)' '
          DO 6 I=0,15
            WRITE(IFCH,*)(EXP( OVEL(1+I,1+J) ),J=0,MIN(4,I))
 6        CONTINUE
        ENDIF
      ENDIF

      RETURN
      END
C=======================================================================

      SUBROUTINE UTPAGE

C-----------------------------------------------------------------------
C  INCREASES IPAGE BY 1, CHANGES ISH
C-----------------------------------------------------------------------
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CISHI/   ISHI
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      SAVE
C-----------------------------------------------------------------------
      IPAGE=IPAGE+1
      IF ( IPAGI .LE. 0 ) RETURN
      IF ( IPAGE .EQ. 1 ) ISHI=ISH
      ISH=0
      IF ( IPAGE.GE.IPAGI/10000.AND.IPAGE.LE.MOD(IPAGI,10000) ) ISH=ISHI
      RETURN
      END
C=======================================================================

      SUBROUTINE UTPART

C----------------------------------------------------------------------
C  FILLS ARRAY PARTX(K,N) CONTAINING THE NUMBER
C    OF PARTITIONS OF N INTO AT MOST K INTEGERS (K>=1 N>=1)  .
C  FILLS ARRAY PART(K,N) CONTAINING  THE NUMBER
C    OF PARTITIONS OF N INTO K INTEGERS (K>=1 N>=1)  .
C  FILLS ARRAY IPART(N,J) CONTAINING THE PARTITIONS OF N.
C----------------------------------------------------------------------
      PARAMETER (NQUAX=12)
      PARAMETER (JPAMAX=NQUAX*NQUAX*NQUAX)
      PARAMETER (KPAMAX=NQUAX)
      COMMON /CPART/   PART(KPAMAX,NQUAX)
      COMMON /CPARTA/  PARTA(NQUAX),IPART(NQUAX,JPAMAX)
      COMMON /CPARTX/  PARTX(KPAMAX,NQUAX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      INTEGER ISWAP(JPAMAX)
      SAVE
C----------------------------------------------------------------------
      WRITE(IFMT,*)'EXECUTE SR UTPART ...'
      DO 10 N=1,NQUAX
        PARTX(1,N)=1.
        PART(1,N)=1.
        DO 10 J=1,JPAMAX
          IPART(N,J)=0
10    CONTINUE

      DO 1 K=2,KPAMAX
        DO 2 N=1,NQUAX
          U=0.
          DO 3 I=0,N/K
            IF ( N-I*K .EQ. 0 ) THEN
              U=U+1.
            ELSE
              U=U+PARTX(K-1,N-I*K)
            ENDIF
 3        CONTINUE
          PARTX(K,N)=U
 2      CONTINUE
 1    CONTINUE
      DO 7 N=1,NQUAX
        PARTA(N)=PART(1,N)
        DO 7 K=2,KPAMAX
          PART(K,N)=PARTX(K,N)-PARTX(K-1,N)
          PARTA(N)=PARTA(N)+PART(K,N)
 7    CONTINUE

      IF ( ISH .GE. 93 ) THEN
121     FORMAT(1X,79A1)
        WRITE(IFCH,*)' '
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   PARTITIONS OF N INTO AT MOST K INTEGERS'
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,101)((K),K=1,10),KPAMAX
101     FORMAT(9X,'K =',10I5,3X,I5)
        WRITE(IFCH,121)('-',IC=1,79)
        DO 8 N=1,NQUAX
          WRITE(IFCH,100)N,(NINT(PARTX(K,N)),K=1,10),NINT(PARTA(N))
100       FORMAT(3X,'N = ',I2,3X,10I5,3X,I5)
 8      CONTINUE
      ENDIF

C  N=1
C  ---
      IPART(1,1)=-1
      IPART(1,2)=-1
      IPART(1,3)=1

C  N>1
C  ---
      DO 11 N=2,NQUAX
        WRITE(IFMT,*)'SR UTPART: N=',N
        IF ( ISH .GE. 95 ) THEN
          WRITE(IFCH,*)' '
          WRITE(IFCH,*)'N=',N
        ENDIF
        IPART(N,1)=-1
        IPART(N,2)=-1
        IPART(N,3)=N
        II=3

C  N=N1+N2, MULTIPLY PARTITIONS OF N1 AND N2
C  -----------------------------------------
        DO 12 N1=1,N/2
          N2=N-N1
          IF ( ISH .GE. 95 ) THEN
            WRITE(IFCH,*)' '
            WRITE(IFCH,*)'N1,N2=',N1,N2
          ENDIF
          K1MAX=-IPART(N1,1)
          K2MAX=-IPART(N2,1)
          IF ( K1MAX .LT. 0  .OR.  K2MAX .LT. 0 ) THEN
            CALL UTSTOP('UTPART: KIMAX NEGATIVE                  ')
          ENDIF
          J1=2
          DO 13 K1=1,K1MAX
            L1=-IPART(N1,J1)
            IF ( L1 .LT. 0 ) THEN
              CALL UTSTOP('UTPART: L1 NEGATIVE (1)                 ')
            ENDIF
            J2=2
            DO 15 K2=1,K2MAX
              L2=-IPART(N2,J2)
              IF ( L2 .LT. 0 ) THEN
                CALL UTSTOP('UTPART: L2 NEGATIVE (1)                 ')
              ENDIF

              IPART(N,1)=IPART(N,1)-1
              II=II+1
              IF ( II .GT. JPAMAX ) GOTO 1000
              IPART(N,II)=-L1-L2
              II0=II+1
              DO 17 I1=1,L1
                II=II+1
                IF ( II .GT. JPAMAX ) GOTO 1000
                IPART(N,II)=IPART(N1,J1+I1)
17            CONTINUE
              DO 18 I2=1,L2
                II=II+1
                IF ( II .GT. JPAMAX ) GOTO 1000
                IPART(N,II)=IPART(N2,J2+I2)
18            CONTINUE

              IF ( ISH .GE. 95 ) THEN
                WRITE(IFCH,*)' '
                WRITE(IFCH,*)'K1,K2=',K1,K2
                WRITE(IFCH,103)N,-IPART(N,1),NINT(PARTA(N))
103             FORMAT(/,3X,'N = ',I2,'   P(N) = ',I4
     *                                    ,'   P0(N) = ',I4/)
                JY=2
                DO 29 KY=1,-IPART(N,1)
                  LY=-IPART(N,JY)
                  WRITE(IFCH,102)KY,LY,(IPART(N,I)
     *                                        ,I=JY+1,JY+MIN(20,LY))
102               FORMAT(3X,I2,3X,I2,3X,20I3)
                  JY=JY+LY+1
29              CONTINUE
              ENDIF

              IPRI=0
              DO 19 I=II0,II-1
                DO 19 J=I+1,II
                  IF ( IPART(N,I) .LT. IPART(N,J) ) THEN
                    IPRI=1
                    ISTI=IPART(N,I)
                    IPART(N,I)=IPART(N,J)
                    IPART(N,J)=ISTI
                  ENDIF
19            CONTINUE
              IF ( IPRI.EQ.1 .AND. ISH.GE.95 ) THEN
                WRITE(IFCH,*)' '
                WRITE(IFCH,*)'SEQUENCE CHANGED'
              ENDIF

              JX=2
              DO 28 KX=1,-IPART(N,1)-1
                LX=-IPART(N,JX)
                IF ( LX .LT. 0 ) THEN
                 CALL UTSTOP('UTPART: LX NEGATIVE                     ')
                ENDIF
                IF ( LX .EQ. L1+L2 ) THEN
                  DO 22 L=1,LX
                    IF ( IPART(N,JX+L) .NE. IPART(N,II0-1+L) ) GOTO 23
22                CONTINUE
                  IF ( ISH .GE. 95 ) THEN
                    WRITE(IFCH,*)' '
                    WRITE(IFCH,*)'EXISTS ALREADY'
                  ENDIF
                  IPRI=1
                  DO 33 I=II0-1,II
                    IPART(N,I)=0
33                CONTINUE
                  II=II0-2
                  IPART(N,1)=IPART(N,1)+1
                  GOTO 30
23                CONTINUE
                ENDIF
                JX=JX+LX+1
28            CONTINUE
30            CONTINUE

              IF ( IPRI.EQ.1 .AND. ISH.GE.95 ) THEN
                WRITE(IFCH,103)N,-IPART(N,1),NINT(PARTA(N))
                JY=2
                DO 32 KY=1,-IPART(N,1)
                  LY=-IPART(N,JY)
                  WRITE(IFCH,102)KY,LY,(IPART(N,I)
     *                                ,I=JY+1,JY+MIN(20,LY))
                  JY=JY+LY+1
32              CONTINUE
              ENDIF

              J2=J2+L2+1
15          CONTINUE
            J1=J1+L1+1
13        CONTINUE
12      CONTINUE

C  ORDERING
C  --------
        KM=-IPART(N,1)
        IF ( KM .LT. 0 ) THEN
          CALL UTSTOP('UTPART: KM NEGATIVE                     ')
        ENDIF
        J1=2
        DO 20 K1=1,KM-1
          L1=-IPART(N,J1)
          IF ( L1 .LT. 0 ) THEN
            CALL UTSTOP('UTPART: L1 NEGATIVE (2)                 ')
          ENDIF
          J2=2
          DO 21 K2=1,KM
            L2=-IPART(N,J2)
            IF ( L2 .LT. 0 ) THEN
              CALL UTSTOP('UTPART: L2 NEGATIVE (2)                 ')
            ENDIF
            IF ( K2 .LE. K1 ) GOTO 21
            IF ( L1 .GT. L2 ) THEN
              DO 24 L=1,L2+1
                ISWAP(L)=IPART(N,J2-1+L)
24            CONTINUE
              DO 25 L=1,J2-J1
                I=J2-L
                IPART(N,I+L2+1)=IPART(N,I)
25            CONTINUE
              DO 26 L=1,L2+1
                IPART(N,J1-1+L)=ISWAP(L)
26            CONTINUE
              L1=-IPART(N,J1)
              IF ( L1 .LT. 0 ) THEN
                CALL UTSTOP('UTPART: L1 NEGATIVE (2)                 ')
              ENDIF
              IF ( ISH .GE. 95 ) THEN
                WRITE(IFCH,*)' '
                WRITE(IFCH,*)'ORDER CHANGED.      K1,K2=',K1,K2
                WRITE(IFCH,103)N,-IPART(N,1),NINT(PARTA(N))
                JY=2
                DO 31 KY=1,-IPART(N,1)
                  LY=-IPART(N,JY)
                  WRITE(IFCH,102)KY,LY,(IPART(N,I)
     *                                ,I=JY+1,JY+MIN(20,LY))
                  JY=JY+LY+1
31              CONTINUE
              ENDIF
            ENDIF
            J2=J2+L2+1
21        CONTINUE
          J1=J1+L1+1
20      CONTINUE

        IF ( ISH.GE.93 .AND. N.LE.8 ) THEN
          WRITE(IFCH,113)('-',IC=1,79),N,-IPART(N,1),('-',IC=1,79)
113       FORMAT(/,1X,79A1,/,7X,'N = ',I2,'   --->   ',I4,
     *                           ' PARTITIONS',/,1X,79A1)
          JY=2
          DO 27 KY=1,-IPART(N,1)
            LY=-IPART(N,JY)
            WRITE(IFCH,112)KY,(IPART(N,I),I=JY+1,JY+MIN(20,LY))
112         FORMAT(2X,I2,'. PARTITION:',3X,20I3)
            JY=JY+LY+1
27        CONTINUE
        ENDIF
        IF ( -IPART(N,1) .NE. NINT(PARTA(N)) ) THEN
          CALL UTSTOP('UTPART: # OF PARTITIONS WRONG           ')
        ENDIF
11    CONTINUE

      RETURN

1000  WRITE(IFCH,*)' '
      WRITE(IFCH,*)('*',J=1,79)
      WRITE(IFCH,*)'***** N=',N
      WRITE(IFCH,*)'***** JPAMAX=',JPAMAX
      WRITE(IFCH,*)('*',J=1,79)
      CALL UTSTOP('UTPART: DIMENSION JPAMAX TOO SMALL.     ')
      RETURN
      END
C=======================================================================

      FUNCTION UTPCM(A,B,C)

C-----------------------------------------------------------------------
C  CALCULATES CM MOMENTUM FOR A-->B+C
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      VAL=(A**2-B**2-C**2)**2-(2.*B*C)**2
      IF ( VAL.LT.0. .AND. VAL.GT.-1.E-4 ) THEN
        UTPCM=0.
        RETURN
      ENDIF
      UTPCM=SQRT(VAL)/(2.*A)
      RETURN
      END
C=======================================================================

      SUBROUTINE UTQUAF(FU,N,X,Q,X0,X1,X2,X3)

C-----------------------------------------------------------------------
C  RETURNS Q(I) = INTEGRAL (X(1)->X(I)) FU(X) DX
C  ACCELERATED VERSION BY      D. HECK, KFK    SEPT 20, 1993
C-----------------------------------------------------------------------
      PARAMETER (M=10)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      REAL Q(N),X(N)
      SAVE
C-----------------------------------------------------------------------
      QUOT = 1./FLOAT(M-1)
      IF ( ISH .GE. 90 ) THEN
        IF ( X1.LT.X0 .OR. X2.LT.X1 .OR. X3.LT.X2 ) THEN
          CALL UTMSG('UTQUAF')
          WRITE(IFCH,*)'   XI=',X0,X1,X2,X3
          CALL UTMSGF
        ENDIF
      ENDIF
      I1 = N/3
      I2 = 2*N/3
      FAC1 = (X1-X0)/FLOAT(I1-1)
      DO 11 I=1,I1-1
        X(I)=X0+(I-1.)*FAC1
 11   CONTINUE
      FAC2 = (X2-X1)/FLOAT(I2-I1)
      DO 12 I=I1,I2-1
        X(I)=X1+FLOAT(I-I1)*FAC2
 12   CONTINUE
      FAC3 = (X3-X2)/FLOAT(N-I2)
      DO 13 I=I2,N
        X(I)=X2+FLOAT(I-I2)*FAC3
 13   CONTINUE
      Q(1) = 0.
      Z = X(1)
      AUXIL = FU(Z)
      DO 2 I=2,N
        FACT = (X(I) - Z) * QUOT
        UTQUAD = AUXIL*0.5
        DO 3 K=2,M-1
          Z=Z+FACT
          UTQUAD=FU(Z)+UTQUAD
 3      CONTINUE
        Z = X(I)
        AUXIL = FU(Z)
        Q(I)=(AUXIL*0.5+UTQUAD)*FACT+Q(I-1)
 2    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE UTQZ(N,X,Q1,Q2,X0,X1,X2,X3)

C-----------------------------------------------------------------------
C  RETURNS Q1(I) = INTEGRAL (X(1)->X(I)) FU1(X) DX  (SEA QUARK STRUCT)
C  RETURNS Q2(I) = INTEGRAL (X(1)->X(I)) FU2(X) DX  (VAL QUARK STRUCT)
C  ACCELERATED VERSION BY      D. HECK, KFK    OCT  20, 1993
C-----------------------------------------------------------------------
      PARAMETER (M=10)
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      COMMON /CIPIO/   IPIO
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      REAL Q1(N),Q2(N),X(N)
      SAVE
C-----------------------------------------------------------------------
      IF ( ISH .GE. 90 ) THEN
        IF ( X1.LT.X0 .OR. X2.LT.X1 .OR. X3.LT.X2 ) THEN
          CALL UTMSG('UTQZ  ')
          WRITE(IFCH,*)'   XI=',X0,X1,X2,X3
          CALL UTMSGF
        ENDIF
      ENDIF
      I1 = N/3
      I2 = 2*N/3
      FAC1 = (X1-X0)/FLOAT(I1-1)
      DO 11 I=1,I1-1
        X(I)=(I-1.)*FAC1+X0
 11   CONTINUE
      FAC2 = (X2-X1)/FLOAT(I2-I1)
      DO 12 I=I1,I2-1
        X(I)=FLOAT(I-I1)*FAC2 +X1
 12   CONTINUE
      FAC3 = (X3-X2)/FLOAT(N-I2)
      DO 13 I=I2,N
        X(I)=MIN( FLOAT(I-I2)*FAC3 +X2, 0.999999 )
 13   CONTINUE
      XCUT2 = XCUT**2
      QUOT= 1./FLOAT(M-1)
      Q1(1)=0.
      Q2(1)=0.
      Z = X(1)
      DENOMI= 1./SQRT(Z**2 + XCUT2)
      IF ( IPIO .EQ. 0 ) THEN
        AUXIL1 = (1.-Z)**8.05 *DENOMI
        IF ( Z .NE. 0. ) THEN
          AUXIL2 = (1.-Z)**3.46 * Z**.419
     *                  * (2.74793064*Z+0.62452969)* DENOMI
        ELSE
          AUXIL2 = 0.
        ENDIF
        DO 3 I=2,N
          FACT = (X(I) - Z) * QUOT
          UTQUA1 = 0.5*AUXIL1
          UTQUA2 = 0.5*AUXIL2
          DO 2 K=2,M-1
            Z=Z+FACT
            DENOMI = 1./SQRT(Z**2 + XCUT2)
            UTQUA1 = DENOMI * (1.-Z)**8.05 + UTQUA1
            IF ( Z .NE. 0. ) THEN
              UTQUA2 = (1.-Z)**3.46 * Z**.419
     *                  * (2.74793064*Z+0.62452969) * DENOMI + UTQUA2
            ENDIF
 2        CONTINUE
          Z=X(I)
          DENOMI = 1./SQRT(Z**2 + XCUT2)
          AUXIL1 = DENOMI * (1.-Z)**8.05
          Q1(I) = (AUXIL1*0.5+UTQUA1) * FACT*1.265 + Q1(I-1)
          IF ( Z .NE. 0. ) THEN
            AUXIL2=(1.-Z)**3.46 * Z**.419 * (2.74793064*Z+0.62452969)
     *                  * DENOMI
            Q2(I) = (AUXIL2*0.5+UTQUA2) * FACT + Q2(I-1)
          ELSE
            AUXIL2 = 0.
            Q2(I) = FACT* UTQUA2 + Q2(I-1)
          ENDIF
 3      CONTINUE

      ELSE
        CUTLOG= LOG(XCUT)
        A0 = -5. + 6.6666667*XCUT2 - 0.53333333*XCUT2**2
        BA0= A0*XCUT
        A1 = 5. - 1.875*XCUT2
        QB = 1. - A1*XCUT2
        A2 = -3.3333333 +0.26666667*XCUT2
        A3 = 1.25
        A4 = -0.2
        AUXIL1 = (1.-Z)**5.0 * DENOMI
        IF ( Z .NE. 0. ) THEN
          AUXIL2 = (1.-Z)**0.7 * Z**0.4 * DENOMI
        ELSE
          AUXIL2 = 0.
        ENDIF
        DO 5 I=2,N
          FACT = (X(I) - Z) * QUOT
          UTQUA1 = 0.5*AUXIL1
          UTQUA2 = 0.5*AUXIL2
          DO 4 K=2,M-1
            Z=Z+FACT
            DENOMI = 1./SQRT(Z**2 + XCUT2)
            UTQUA1 = DENOMI * (1.-Z)**5.0 + UTQUA1
            IF ( Z .NE. 0. ) THEN
              UTQUA2 = DENOMI * (1.-Z)**0.7 * Z**0.4 + UTQUA2
            ENDIF
 4        CONTINUE
          Z = X(I)
          ROOT = SQRT(Z**2 + XCUT2)
          DENOMI= 1./ROOT
CC        AUXIL1=DENOMI * (1.-Z)**5.0
CC        Q1(I) = (AUXIL1*0.5+UTQUA1) * FACT*0.9 + Q1(I-1)
          Q1(I) = ( QB * ( LOG(Z+ROOT) - CUTLOG ) - BA0
     *           + ROOT * (A0+Z*(A1+Z*(A2+Z*(A3+Z*A4)))) ) *0.9
          IF ( Z .NE. 0. ) THEN
            AUXIL2 = DENOMI * (1.-Z)**0.7 * Z**0.4
            Q2(I) = (AUXIL2*0.5+UTQUA2) * FACT*0.1730725 + Q2(I-1)
          ELSE
            AUXIL2 = 0.
            Q2(I) = FACT*0.1730725 * UTQUA2 + Q2(I-1)
          ENDIF
 5      CONTINUE
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTREMB(PROJ,TARG,II)

C-----------------------------------------------------------------------
C  REMEMBERS NEVT,NSTR,NPTL,PROJ,TARG
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      PARAMETER (MXSTR=3000)
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /CREMB/   PROJRE(2,NSI,NHA),TARGRE(2,NSI,NHA)
     *                ,NEVTRE(2),NPTLRE(2),NSTRRE(2)
      COMMON /CSTR/    PSTR(5,MXSTR),ROTSTR(3,MXSTR),XORSTR(4,MXSTR)
     *                ,ICSTR(4,MXSTR),IORSTR(MXSTR),IRLSTR(MXSTR),NSTR
      REAL PROJ(NSI,NHA),TARG(NSI,NHA)
      SAVE
C-----------------------------------------------------------------------
      NEVTRE(II)=NEVT
      NSTRRE(II)=NSTR
      NPTLRE(II)=NPTL
      DO 56 M=1,NHA
        SM=0.
        DO 57 N=1,NSI
          PROJRE(II,N,M)=PROJ(N,M)
          SM=SM+PROJ(N,M)**2
57      CONTINUE
        IF ( M.GE.3 .AND. SM.LT.1.E-5 ) GOTO 58
56    CONTINUE
58    CONTINUE
      DO 53 M=1,NHA
        SM=0.
        DO 54 N=1,NSI
          TARGRE(II,N,M)=TARG(N,M)
          SM=SM+TARG(N,M)**2
54      CONTINUE
        IF ( M.GE.3 .AND. SM.LT.1.E-5 ) GOTO 55
53    CONTINUE
55    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE UTREPL(I,J)

C-----------------------------------------------------------------------
C  I IS REPLACED BY J IN /CPTL/
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      COMMON /C2PTL/   AMIPTL(MXPTL),RADPTL(MXPTL),IAAPTL(MXPTL)
      SAVE
C-----------------------------------------------------------------------
      AMIPTL(I)  =AMIPTL(J)
      IAAPTL(I)  =IAAPTL(J)
      IBPTL(1,I) =IBPTL(1,J)
      IBPTL(2,I) =IBPTL(2,J)
      IBPTL(3,I) =IBPTL(3,J)
      IBPTL(4,I) =IBPTL(4,J)
      ICLPTL(I)  =ICLPTL(J)
      IDPTL(I)   =IDPTL(J)
      IFRPTL(1,I)=IFRPTL(1,J)
      IFRPTL(2,I)=IFRPTL(2,J)
      IORPTL(I)  =IORPTL(J)
      ISTPTL(I)  =ISTPTL(J)
      JORPTL(I)  =JORPTL(J)
      NQJPTL(I)  =NQJPTL(J)
      DO 1 K=1,5
        PPTL(K,I)=PPTL(K,J)
 1    CONTINUE
      RADPTL(I)  =RADPTL(J)
      TIVPTL(1,I)=TIVPTL(1,J)
      TIVPTL(2,I)=TIVPTL(2,J)
      XORPTL(1,I)=XORPTL(1,J)
      XORPTL(2,I)=XORPTL(2,J)
      XORPTL(3,I)=XORPTL(3,J)
      XORPTL(4,I)=XORPTL(4,J)

      RETURN
      END
C=======================================================================

      SUBROUTINE UTRESM(ICP1,ICP2,ICM1,ICM2,AMP,IDPR,IADJ,IRETEN)

C-----------------------------------------------------------------------
      PARAMETER (NFLAV=6)
      INTEGER ICM(2),ICP(2),JCM(NFLAV,2),JCP(NFLAV,2)
      SAVE
C-----------------------------------------------------------------------
      ICM(1)=ICM1
      ICM(2)=ICM2
      ICP(1)=ICP1
      ICP(2)=ICP2
      CALL IDDECO(ICM,JCM)
      CALL IDDECO(ICP,JCP)
      DO 37 NF=1,NFLAV
        JCP(NF,1)=JCP(NF,1)+JCM(NF,1)
        JCP(NF,2)=JCP(NF,2)+JCM(NF,2)
37    CONTINUE
      CALL IDENCO(JCP,ICP,IRETEN)
      IDP=IDTRA(ICP,0,0,3)
      CALL IDRES(IDP,AMP,IDPR,IADJ)
      RETURN
      END
C=======================================================================

      SUBROUTINE UTREST(PROJ,TARG,II)

C-----------------------------------------------------------------------
      PARAMETER (NSI=6)
      PARAMETER (NSIX=40)
      PARAMETER (NHA=NSIX+2)
      COMMON /CEVT/    BIMEVT,COLEVT,EGYEVT,PHIEVT,PMXEVT
     *                ,KOLEVT,NEVT,NPJEVT,NTGEVT
      COMMON /CREMB/   PROJRE(2,NSI,NHA),TARGRE(2,NSI,NHA)
     *                ,NEVTRE(2),NPTLRE(2),NSTRRE(2)
      REAL PROJ(NSI,NHA),TARG(NSI,NHA)
      SAVE
C-----------------------------------------------------------------------
      NEVT=NEVTRE(II)
      NSTR=NSTRRE(II)
      NPTL=NPTLRE(II)
      DO 66 M=1,NHA
        SM=0.
        DO 67 N=1,NSI
          PROJ(N,M)=PROJRE(II,N,M)
          SM=SM+PROJRE(II,N,M)**2
67      CONTINUE
        IF ( M.GE.3 .AND. SM.LT.1.E-5 ) GOTO 68
66    CONTINUE
68    CONTINUE
      DO 63 M=1,NHA
        SM=0.
        DO 64 N=1,NSI
          TARG(N,M)=TARGRE(II,N,M)
          SM=SM+TARGRE(II,N,M)**2
64      CONTINUE
        IF ( M.GE.3 .AND. SM.LT.1.E-5 ) GOTO 65
63    CONTINUE
65    CONTINUE
      RETURN
      END
C=======================================================================

      SUBROUTINE UTROTA(ISIG,AX,AY,AZ,X,Y,Z)

C-----------------------------------------------------------------------
C  PERFORMS A ROTATION
C-----------------------------------------------------------------------
      SAVE
C-----------------------------------------------------------------------
      IF ( AZ .GE. 0. ) THEN
        RX=AX
        RY=AY
        RZ=AZ
      ELSE
        RX=-AX
        RY=-AY
        RZ=-AZ
      ENDIF
      ALP=SIGN( ABS(UTACOS(RZ/SQRT(RZ**2+RY**2))), RY )
      BET=SIGN( ABS(UTACOS(SQRT(RZ**2+RY**2)/SQRT(RZ**2+RY**2+RX**2)))
     *                                                        , RX )
      COSA=COS(ALP)
      SINA=SIN(ALP)
      COSB=COS(BET)
      SINB=SIN(BET)
      IF     ( ISIG .GT. 0 ) THEN
        XS=X*COSB-Y*SINA*SINB-Z*COSA*SINB
        YS=       Y*COSA     -Z*SINA
        ZS=X*SINB+Y*SINA*COSB+Z*COSA*COSB
      ELSEIF ( ISIG .LT. 0 ) THEN
        XS= X*COSB            +Z*SINB
        YS=-X*SINB*SINA+Y*COSA+Z*COSB*SINA
        ZS=-X*SINB*COSA-Y*SINA+Z*COSB*COSA
      ENDIF
      X=XS
      Y=YS
      Z=ZS
      RETURN
      END
C=======================================================================

      SUBROUTINE UTROT2(ISIG,AX,AY,AZ,X,Y,Z)

C-----------------------------------------------------------------------
C  PERFORMS A ROTATION, DOUBLE PREC.
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO1/   AMPRIF,AMSIAC,BMAXIM,BMINIM,CORE,CUTMSQ,CUTMSS
     *                ,DELMSS,DELREM,FCTRMX,GAUMX,OVERLP,PAREA,PDIQUA
     *                ,PHARD,PSPINL,PSPINH,PISPN,PTF,PTH,PTMX,PTQ,PUD
     *                ,PVALEN,QSEPC,QSETC,QMUST,QVAPC,QVATC,RADIAC
     *                ,RADIAS,RSTRAS,SIGJ,SIGPPI,TAUMAX,TAUMIN
     *                ,TAUMX,TAUNLL,TENSN,THEMAS,WPROJ,WTARG,WTMINI
     *                ,WTSTEP,XCUT
     *                ,IAQU,IFRADE,IOJINT,IOPBRK,IOPENT,IOPENU
     *                ,IOPTF,IOPTQ,IRESCL,IWCENT,KENTRO,KO1KO2
     *                ,LABSYS,MAXRES,NCLEAN,NCOLMX,NDECAW,NEQMN,NEQMX
     *                ,NSTTAU,NTRYMX,NUMTAU

      DOUBLE PRECISION ALP,AUXIL1,AUXIL2,AX,AY,AZ,BET,COSA,COSB
     *                ,RX,RY,RZ,SINA,SINB,X,XS,Y,YS,Z,ZS
      SAVE
C-----------------------------------------------------------------------
      IF ( AX**2.EQ.0.D0 .AND. AY**2.EQ.0.D0 .AND. AZ**2.EQ.0.D0 ) THEN
        WRITE(IFCH,*)' '
        WRITE(IFCH,*)'AX**2,AY**2,AZ**2:',AX**2,AY**2,AZ**2
        WRITE(IFCH,*)'AX,AY,AZ:',AX,AY,AZ
        CALL UTSTOP('UTROT2: ZERO VECTOR.                    ')
      ENDIF
      IF ( AZ .GE. 0.D0 ) THEN
        RX=AX
        RY=AY
        RZ=AZ
      ELSE
        RX=-AX
        RY=-AY
        RZ=-AZ
      ENDIF
      AUXIL1 = RZ**2+RY**2
      IF ( AUXIL1 .NE. 0.D0 ) THEN
        AUXIL2 = SQRT(AUXIL1)
        ALP=SIGN( ABS(ACOS(RZ/AUXIL2)), RY )
        BET=SIGN( ABS(ACOS( AUXIL2/SQRT(AUXIL1+RX**2) )), RX )
        COSA=COS(ALP)
        SINA=SIN(ALP)
        COSB=COS(BET)
        SINB=SIN(BET)
      ELSE
        COSA=0.D0
        COSB=0.D0
        SINA=1.D0
        SINB=1.D0
      ENDIF
      IF     ( ISIG .GT. 0 ) THEN
        XS=X*COSB-Y*SINA*SINB-Z*COSA*SINB
        YS=       Y*COSA     -Z*SINA
        ZS=X*SINB+Y*SINA*COSB+Z*COSA*COSB
      ELSEIF ( ISIG .LT. 0 ) THEN
        XS= X*COSB            +Z*SINB
        YS=-X*SINB*SINA+Y*COSA+Z*COSB*SINA
        ZS=-X*SINB*COSA-Y*SINA+Z*COSB*COSA
      ENDIF
      X=XS
      Y=YS
      Z=ZS
      RETURN
      END
C=======================================================================

      SUBROUTINE UTSTOP(TEXT)

C-----------------------------------------------------------------------
C  RETURNS ERROR MESSAGE AND STOPS EXECUTION
C-----------------------------------------------------------------------
      COMMON /ACCUM/   AMSAC,ILAMAS,IMSG,INOIAC,IPAGE,JERR,NAEVT,NREVT
     *                ,NRPTL,NRSTR,NTEVT
      DOUBLE PRECISION SEEDC,SEEDI
      COMMON /CSEED/   SEEDC,SEEDI
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /NEVNT/   NEVNT
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      CHARACTER STAR*1,TEXT*40
      DATA STAR/'*'/
      SAVE
C-----------------------------------------------------------------------
      DO 1 I=1,2
        IF ( I .EQ. 1 ) THEN
          IFI=IFCH
        ELSE
          IFI=IFMT
        ENDIF
        WRITE(IFI,*)' '
        WRITE(IFI,*)(STAR,J=1,77)
        WRITE(IFI,*)'***** STOP IN ',TEXT
CDH     WRITE(IFI,*)'***** CURRENT EVENT NUMBER: ',NREVT+1
        WRITE(IFI,*)'***** CURRENT EVENT NUMBER: ',NEVNT
        WRITE(IFI,*)'***** CURRENT PAGE NUMBER: ',IPAGE
        WRITE(IFI,*)'***** INITIAL SEED FOR CURRENT EVENT:',SEEDI
        WRITE(IFI,*)'***** RANDOM CALLS FOR CURRENT EVENT:',SEEDC
        WRITE(IFI,*)(STAR,J=1,77)
        WRITE(IFI,*)' '
 1    CONTINUE
      STOP
      END
C=======================================================================

      SUBROUTINE UTTAIN(I,X,Y,Z,T,N,IOPT)

C-----------------------------------------------------------------------
C  RETURNS INTERSECTION OF PTL-I-TRAJECTORY WITH TAUS-LINE.
C  N=0 IF OK, N=1 IF PTL LIVES LATER, N=2 IF EARLIER,
C    N=9 IF TIV1>TIV2, N=3,4,5 ELSE.
C  IOPT=0: FORMATION TIME CONSIDERED, IOPT=1 IF NOT
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION DD,DERR,TI1,TI2,TT,VV,VVP,VVT,XO3,XO4,ZZ,ZZA
      DATA DERR/1.D-2/
      SAVE
C-----------------------------------------------------------------------
      XO4=XORPTL(4,I)
      IF     ( IOPT .EQ. 0 ) THEN
        TI1=TIVPTL(1,I)
      ELSEIF ( IOPT .EQ. 1 ) THEN
        TI1=XO4
      ENDIF
      TI2=TIVPTL(2,I)

      IF ( TI1 .GT. TI2 ) GOTO 1009
      PPT4I = 1./PPTL(4,I)
      VV=PPTL(3,I)*PPT4I
      XO3=XORPTL(3,I)

      ZZ=XO3+(TI2-XO4)*VV
      IF ( TTAUS .LE. 0.D0 ) THEN
        TZ=TTAUS
      ELSE
        IF     ( ZZ .LE. ZZT ) THEN
          TZ=TTT+(ZZ-ZZT)*ZZT/TTT
        ELSEIF ( ZZ .GE. ZZP ) THEN
          TZ=TTP+(ZZ-ZZP)*ZZP/TTP
        ELSE
          IF ( TTAUS .GE. AINFIN ) THEN
            TZ=TTAUS
            IF ( ISH .GE. 90 ) THEN
              CALL UTMSG('UTTAIN')
              WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
              WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
              CALL UTMSGF
            ENDIF
          ELSE
C*DH        TZ=SQRT(TTAUS**2+ZZ**2)
            IF ( TI2 .LT. 0.D0 ) GOTO 1002
            IF ( TTAUS**2+ZZ**2 .GE. TI2**2 ) GOTO 1002
            GOTO 1006
          ENDIF
        ENDIF
      ENDIF
      IF ( TZ .GE. TI2 ) GOTO 1002

 1006 ZZ=XO3+(TI1-XO4)*VV
      IF ( TTAUS .GT. 0.D0 ) THEN
        IF     ( ZZ .LE. ZZT ) THEN
          TZ=TTT+(ZZ-ZZT)*ZZT/TTT
        ELSEIF ( ZZ .GE. ZZP ) THEN
          TZ=TTP+(ZZ-ZZP)*ZZP/TTP
        ELSE
          IF ( TTAUS .GE. AINFIN ) THEN
            TZ=TTAUS
            IF ( ISH .GE. 90 ) THEN
              CALL UTMSG('UTTAIN')
              WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
              WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
              CALL UTMSGF
            ENDIF
          ELSE
C*DH        TZ=SQRT(TTAUS**2+ZZ**2)
            IF ( TI1 .LT. 0.D0 ) GOTO 1007
            IF ( TTAUS**2+ZZ**2 .LE. TI1**2 ) GOTO 1001
            GOTO 1007
          ENDIF
        ENDIF
      ENDIF
      IF ( TZ .LE. TI1 ) GOTO 1001

 1007 IF ( TTAUS .LE. 0.D0 ) THEN
        TT=TTAUS
        ZZ=XO3+(TT-XO4)*VV
        IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1031
      ELSE
        ZZA=XO3-XO4*VV
        VVT=ZZT/TTT
        TT=(TTT+(ZZA-ZZT)*VVT)/(1.D0-VV*VVT)
        ZZ=XO3+(TT-XO4)*VV
        IF ( ZZ .LE. ZZT ) THEN
          IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1032
          GOTO 1000
        ENDIF
        VVP=ZZP/TTP
        TT=(TTP+(ZZA-ZZP)*VVP)/(1.D0-VV*VVP)
        ZZ=XO3+(TT-XO4)*VV
        IF ( ZZ .GE. ZZP ) THEN
          IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1033
          GOTO 1000
        ENDIF
        DD=1.D0-VV**2
        IF ( DD .EQ. 0.D0 ) THEN
          TT=-VV*(TTAUS**2+ZZA**2)*0.5D0/ZZA
        ELSE
          TT=(ZZA*VV+SQRT(ZZA**2+TTAUS**2*DD))/DD
        ENDIF
        ZZ=XO3+(TT-XO4)*VV
        IF ( TT.LT.TI1 .OR. TT.GE.TI2 ) GOTO 1034
        IF ( TT .LT. 0.D0 ) GOTO 1035
        IF ( ZZ.LE.ZZT .OR. ZZ.GE.ZZP ) GOTO 1004
        IF ( ABS(TTAUS**2-(TT+ZZ)*(TT-ZZ))
     *                         .GT. DERR*TTAUS**2 ) GOTO 1005
      ENDIF

 1000 N=0
 1011 T=TT
      Z=ZZ
      X=XORPTL(1,I)+(T-XO4)*PPTL(1,I)*PPT4I
      Y=XORPTL(2,I)+(T-XO4)*PPTL(2,I)*PPT4I
      RETURN
 1001 N=1
      RETURN
 1002 N=2
      RETURN
 1031 N=31
      GOTO 1003

 1032 N=32
      GOTO 1003
 1033 N=33
      GOTO 1003
 1034 N=34
 1003 IF ( ABS(TT-TI1) .LE. DERR*ABS(TT) ) GOTO 1000
      IF ( ABS(TT-TI2) .LE. DERR*ABS(TT) ) GOTO 1000
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('UTTAIN')
        WRITE(IFCH,*)'*****  TI1 < TT < TI2   NOT FULFILLED - ',N
        WRITE(IFCH,*)SNGL(TI1),SNGL(TT),SNGL(TI2)
        CALL UTMSGF
      ENDIF
      GOTO 1011
 1035 CONTINUE
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('UTTAIN')
        WRITE(IFCH,*)'*****  TT < 0     ( ',TT,' )'
        WRITE(IFCH,*)'VV,DD:',VV,DD
        WRITE(IFCH,*)'ZZA,TTAUS:',ZZA,TTAUS
        CALL UTMSGF
      ENDIF
      GOTO 1011
 1004 N=4
      IF ( ABS(ZZ-ZZT) .LE. DERR*ABS(ZZ) ) GOTO 1000
      IF ( ABS(ZZ-ZZP) .LE. DERR*ABS(ZZ) ) GOTO 1000
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('UTTAIN')
        WRITE(IFCH,*)'*****  ZZT < ZZ < ZZP   NOT FULFILLED'
        WRITE(IFCH,*)SNGL(ZZT),SNGL(ZZ),SNGL(ZZP)
        CALL UTMSGF
      ENDIF
      GOTO 1011
 1005 N=5
      IF ( ABS(TTAUS**2-(TT+ZZ)*(TT-ZZ)) .LE. DERR ) GOTO 1000
      IF ( ISH .GE. 90 ) THEN
        CALL UTMSG('UTTAIN')
        WRITE(IFCH,*)'*****  TTAUS**2 .NE. (TT+ZZ)*(TT-ZZ)'
        WRITE(IFCH,*)SNGL(TTAUS**2),SNGL((TT+ZZ)*(TT-ZZ))
        CALL UTMSGF
      ENDIF
      GOTO 1011
 1009 N=9
      RETURN
      END
C=======================================================================

      SUBROUTINE UTTAIX(I,TAU,ZOR,TOR,Z,T)

C-----------------------------------------------------------------------
C  RETURNS INTERSECTION Z,T OF PTL-I-TRAJECTORY WITH HYPERBOLA H.
C     H: (T-TOR)**2-(Z-ZOR)**2=TAU**2 .
C     ZOR, TOR DOUBLE PRECISION.
C-----------------------------------------------------------------------
      PARAMETER (MXPTL=70000)
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      COMMON /CPTL/    PPTL(5,MXPTL),TIVPTL(2,MXPTL),XORPTL(4,MXPTL)
     *                ,IBPTL(4,MXPTL),ICLPTL(MXPTL),IDPTL(MXPTL)
     *                ,IFRPTL(2,MXPTL),IORPTL(MXPTL),ISTPTL(MXPTL)
     *                ,JORPTL(MXPTL),NPTL,NQJPTL(MXPTL)

      DOUBLE PRECISION CC,DD,DERR,TOR,TORS,TT,TTAU,VV,ZOR,ZORS,ZZ
      SAVE
C-----------------------------------------------------------------------
      DERR=1.D-3
      TTAU=TAU
      ZORS=XORPTL(3,I)-ZOR
      TORS=XORPTL(4,I)-TOR
      VV=PPTL(3,I)/PPTL(4,I)
      IF ( ABS(VV) .GT. 1.D0 ) THEN
        IF ( ISH.GE.90 .AND. ABS(VV).GT.1.001D0 ) THEN
          CALL UTMSG('UTTAIX')
          WRITE(IFCH,*)'*****  !V! > 1'
          WRITE(IFCH,*)'V: ',VV
          WRITE(IFCH,*)'P,E: ',PPTL(3,I),PPTL(4,I)
          CALL UTMSGF
        ENDIF
        VV=SIGN(1.D0,VV)
      ENDIF
      CC=ZORS-TORS*VV
      DD=1.D0-VV**2
      IF     ( DD.EQ.0.D0 .AND. CC.EQ.0.D0 ) THEN
        IF ( TAU .EQ. 0. ) THEN
          TT=0.
        ELSE
          TT=AINFIN
        ENDIF
        ZZ=TT
        GOTO 1000
      ELSEIF ( DD .EQ. 0.D0 ) THEN
        TT=-(TTAU**2+CC**2)*0.5D0/(CC*VV)
      ELSEIF ( DD .LT. 1.D-8 ) THEN
        TT=-(TTAU**2+CC**2)*0.5D0/(CC*VV)
        IF(ISH.GE.90)THEN
          CALL UTMSG('UTTAIX')
          WRITE(IFCH,*)'*****  DD = ',DD,'    TREATED AS ZERO'
          CALL UTMSGF
        ENDIF
      ELSE
        TT=(CC*VV+SQRT(CC**2+TTAU**2*DD))/DD
      ENDIF
      ZZ=CC+TT*VV
      IF ( ISH .GE. 90 ) THEN
        IF ( ABS(TTAU**2-(TT+ZZ)*(TT-ZZ)) .GT. DERR*TTAU**2   .AND.
     *       ABS(TTAU**2-(TT+ZZ)*(TT-ZZ)) .GT. DERR           .AND.
     *       TORS**2-ZORS**2 .LT. 1.D6                 ) THEN
          CALL UTMSG('UTTAIX')
          WRITE(IFCH,*)'*****  TTAU**2 .NE. (TT+ZZ)*(TT-ZZ)'
          WRITE(IFCH,*)SNGL(TTAU**2),SNGL((TT+ZZ)*(TT-ZZ))
          WRITE(IFCH,*)'TAU,T,Z:'
          WRITE(IFCH,*)TAU,TT,ZZ
          WRITE(IFCH,*)'#,ID(PTL):',I,IDPTL(I)
          WRITE(IFCH,*)'ZOR,TOR(STR):',ZOR,TOR
          WRITE(IFCH,*)'ZORS,TORS,P,E(PTL):'
          WRITE(IFCH,*)SNGL(ZORS),SNGL(TORS),PPTL(3,I),PPTL(4,I)
          CALL UTMSGF
        ENDIF
      ENDIF
1000  Z=ZZ+ZOR
      T=TT+TOR
      RETURN
      END
C=======================================================================

      SUBROUTINE UTTAUS(Z,SZ)

C-----------------------------------------------------------------------
C  RETURNS INV Z-COORD SZ CORRESPONDING TO TTAUS AND Z.
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION ZZ
      SAVE
C-----------------------------------------------------------------------
      IF ( TTAUS .LE. 0.D0 ) THEN
        SZ=Z
        RETURN
      ENDIF
      ZZ=Z
      IF     ( ZZ .LE. ZZT ) THEN
        SZ=TTAUS*DETAT+(ZZ-ZZT)/TTAR
      ELSEIF ( ZZ .GE. ZZP ) THEN
        SZ=TTAUS*DETAP+(ZZ-ZZP)/TPRO
      ELSE
        IF ( SNGL(TTAUS) .GE. AINFIN ) THEN
          SZ=0.
          IF ( ISH .GE. 90 ) THEN
            CALL UTMSG('UTTAUS')
            WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS, SZ=0'
            WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
            CALL UTMSGF
          ENDIF
        ELSE
          TZ=SQRT(TTAUS**2+ZZ**2)
          SZ=TTAUS*0.5D0*LOG((TZ+ZZ)/(TZ-ZZ))
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE UTTAUT(Z,TZ)

C-----------------------------------------------------------------------
C  RETURNS TZ = TIME    CORRESPONDING TO TTAUS AND Z
C-----------------------------------------------------------------------
      COMMON /CNSTA/   AINFIN,PI,PIOM,PROM
      DOUBLE PRECISION DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /CTTAUS/  DETAP,DETAT,TPRO,TTAR,TTAUS,TTP,TTT
     *                ,ZPRO,ZTAR,ZZP,ZZT
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT
      DOUBLE PRECISION ZZ
      SAVE
C-----------------------------------------------------------------------
      ZZ=Z
      IF     ( ZZ .LE. ZZT ) THEN
        TZ=TTT+(ZZ-ZZT)*ZZT/TTT
      ELSEIF ( ZZ .GE. ZZP ) THEN
        TZ=TTP+(ZZ-ZZP)*ZZP/TTP
      ELSE
        IF ( TTAUS .GE. AINFIN ) THEN
          TZ=TTAUS
          IF ( ISH .GE. 90 ) THEN
            CALL UTMSG('UTTAUT')
            WRITE(IFCH,*)'*****  LARGE TTAUS; SET TZ=TTAUS'
            WRITE(IFCH,*)'TTAUS=',TTAUS,'ZZ=',ZZ
            CALL UTMSGF
          ENDIF
        ELSE
          TZ=SQRT(TTAUS**2+ZZ**2)
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
C
C     SUBROUTINE UTTIMA(TEXT)
C
C-----------------------------------------------------------------------
C  RETURNS TIME.
C-----------------------------------------------------------------------
C     COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
C     COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
C    *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
C    *                ,YHAHA,YMXIMI,YPJTL
C    *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
C    *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
C    *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
C    *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
C    *                ,MODSHO,NDECAX,NDECAY,NEVENT
C     CHARACTER*15 TEXT
C-----------------------------------------------------------------------
C     TIMAA=0.
C     TIMA=0.
C-C   CALL TIMAX(TIMAA)
C-C   CALL TIMAD(TIMA)
C     IF ( TEXT .EQ. '               ' ) RETURN
C     IF ( ISH .GE. 91 ) WRITE(IFCH,*)' '
C     WRITE(IFCH,100)TEXT,TIMA/5.,TIMAA/5.
C     IF ( ISH .GE. 91 ) WRITE(IFCH,*)' '
C100  FORMAT(1X,A15,5X,F12.5,5X,F12.5)
C     RETURN
C     END
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C
C     ENTRY UTTIMT
C
C-----------------------------------------------------------------------
C-C   CALL TIMAST(1E10)
C     RETURN
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C
C     ENTRY UTTIMX(TIMAAX)
C
C-----------------------------------------------------------------------
C     TIMAAX=0.
C-C   CALL TIMAX(TIMAAX)
C     RETURN
C     END
C=======================================================================

      SUBROUTINE UTTUCL

C----------------------------------------------------------------------
C  FILLS ARRAY TUCL(1+K,1+N) CONTAINING THE LOGARITHM (LN) OF
C  THE NUMBER OF K-TUPELS N_I WITH SUM_I N_I = N:
C    TUCL(1+K,1+N)=LOG((N+1)*(N+2)...*(N+K-1)/(K-1)!)  K>=0 N>=0 .
C  DOUBLE PRECISION TUCL.
C----------------------------------------------------------------------
      PARAMETER (KTUMAX=100)
      PARAMETER (NTUMAX=100)
      DOUBLE PRECISION TUCL
      COMMON /CTUCL/   TUCL(1+KTUMAX,1+NTUMAX)
      COMMON /FILES/   IFCH,IFDT,IFHI,IFMT,IFOP
      COMMON /PARO2/   AMPROJ,AMTARG,ANGMUE,ELEPTI,ELEPTO,ENGY
     *                ,PNLL,PNLLX,PROB(99),PROSEA,RHOPHI,TAUREA
     *                ,YHAHA,YMXIMI,YPJTL
     *                ,ICBAC(99,2),ICFOR(99,2),ICHOIC,ICLHIS,IDPM
     *                ,IDPROJ,IDTARG,IENTRO,IJPHIS,IMIHIS,IPAGI,ISH
     *                ,ISHEVT,ISHSUB,ISPALL,ISPHIS,ISTMAX,ISUP,IVI
     *                ,JPSI,JPSIFI,KUTDIQ,LAPROJ,LATARG,MAPROJ,MATARG
     *                ,MODSHO,NDECAX,NDECAY,NEVENT

      DOUBLE PRECISION DD,RELERR
      SAVE
C----------------------------------------------------------------------
      TUCL(1,1)=-100.D0
      TUCL(2,1)=0.D0
      DO 7 K=2,KTUMAX
        TUCL(1+K,1)=0.D0
 7    CONTINUE
      DO 2 N=2,NTUMAX+1
        TUCL(1,N)=-100.D0
        TUCL(2,N)=0.D0
        DO 1 K=2,KTUMAX
          DD=(N+K-2.D0)/(K-1.D0)
          TUCL(1+K,N)=TUCL(K,N)+LOG(DD)
 1      CONTINUE
 2    CONTINUE

      IF ( ISH .GE. 93 ) THEN
        DO 5 K=2,50
          DO 4 N=1,50
            RELERR=ABS((EXP(TUCL(1+K,1+N))-EXP(TUCL(1+K,N))
     *                  -EXP(TUCL(K,1+N)))/EXP(TUCL(1+K,1+N)))
            IF ( RELERR .GT. 1.D-4 ) THEN
              CALL UTMSG('UTTUCL')
              WRITE(IFCH,*)'*****  TUCL(,) VIOLATES RECURRENCE RELATION'
              WRITE(IFCH,*)SNGL(EXP(TUCL(1+K,1+N)))
     *                    ,SNGL(EXP(TUCL(1+K,N))+EXP(TUCL(K,1+N)))
              CALL UTMSGF
            ENDIF
 4        CONTINUE
 5      CONTINUE
        WRITE(IFCH,*)' '
        WRITE(IFCH,121)('-',IC=1,79)
121     FORMAT(1X,79A1)
        WRITE(IFCH,*)'   DEGENERACY OF ENERGY LEVELS OF K-DIMENSIONAL'
     *              ,' OSCILLATOR'
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   =='
     *          ,'     NUMBER OF K-TUPELS OF LENGTH N (SUM_I N_I = N)'
        WRITE(IFCH,121)('-',IC=1,79)
        WRITE(IFCH,*)'   K:','   1','   2','   3','   4'
        WRITE(IFCH,121)('-',IC=1,79)
        DO 8 N=0,50
          WRITE(IFCH,*)N,(SNGL(EXP( TUCL(1+K,1+N) )),K=1,4)
 8      CONTINUE
      ENDIF

      RETURN
      END

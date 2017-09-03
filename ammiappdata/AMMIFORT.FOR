C     AMMIFORT Version 1.0.  Hugh G. Gauch, Jr.  May 2017.
C     This FORTRAN77 program performs AMMI analysis.
C
C     Compiled with 32-bit Microsoft Developer Studio FORTRAN
C     Powerstation 4.0.
C
C     Memory requirement is about 256 megabytes for M4L = 64000000.
C     This compiler does not allow substantially more memory.
C
C     Soil and Crop Sciences,
C     619 Bradfield Hall,
C     Cornell University,
C     Ithaca, New York  14853
C     Telephone  607-255-7764
C     e-mail hgg1@cornell.edu
C
C     This software is open source and may be redistributed or
C     modified under the BSD 3-Clause License published by the
C     Regents of the University of California.
C
C     Edited 19 Aug, 8 Sep, 28 Sep, 4 Oct 2016;
C     8 Apr, 17 May 2017.
C
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8        V(32000000)
      INTEGER*4     L(64000000)
      CHARACTER*4   C(64000000)
      INTEGER*1     M(256000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      EQUIVALENCE  (V(1), Y(1), L(1), C(1), M(1))
      M4L = 64000000
C
C     Memory map for data vector with equivalenced arrays V Real*8,
C     Y Real*4, L Integer*4, and C Character*4.
C     G genotypes, E environments, R replicates, and N IPC axes.
C
C     Data            Vector     Offset             Length
C
C     SCRATCH            V       0                  KDGE*
C     TRT means          V       KDGE               GE
C     GEN means          V       KBG=KDGE+GE        G
C     ENV means          V       KBE=KBG+G          E
C     GEN scratch        V       KSCG=KBE+E         G
C     GEN IPC            V       KPGN=KSCG+G        GN
C     ENV IPC            V       KPEN=KPGN+GN       EN
C
C     GEN names          C       KNG=2*(KPEN+EN)    G
C     ENV names          C       KNE=KNG+G          E
C
C     TRT status         M       KSGE=4*(KNE+E)     GE
C     This item exists only if there are missing cells.
C
C     * KDGE = MAX(NGE, (NGER+1)/2), (6*NG)+NE)
C
C     The storage order within vectors having two or three
C     implicit subscripts is as follows.
C     For the replicated data temporarily stored in scratch Y(0),
C     the storage order is to increment REP fastest, GEN next,
C     and ENV slowest.  Hence REP 1, GEN 1, ENV 1 is entry 1;
C     REP 2, GEN 1, ENV 1 is entry 2; REP 1, GEN 2, ENV 1 is
C     entry R+1; and REP 1, GEN 1, ENV 2 is entry RG+1.
C     For the TRT vector V(KDGE), the storage order is to
C     increment GEN fastest and ENV slowest.
C     For the GEN IPC and ENV IPC vectors V(KPGN) and V(KPEN),
C     the storage order is to increment GEN (or ENV) fastest
C     and the IPC axis slowest.
C
C     Announce AMMISOFT program.
      WRITE (*, '('' AMMISOFT Version 1.0:  AMMI Analysis of '',
     1 ''Yield Trials.'', /, '' This is open source software as'',
     2 '' explained in its documentation.'', /, '' Hugh Gauch,'',
     3 '' Soil and Crop Sciences, Cornell University.'', /, /)')
      WRITE (*, '('' Computing ...'')')
C     Read input and output filenames; open files.
      OPEN (3, FILE='InHouseFilenames.txt', STATUS='OLD')
      READ (3, '(A160)')  FNDA
      READ (3, '(A160)')  FNAM
      CLOSE (3)
      OPEN (1, FILE=FNDA, STATUS='OLD')
      OPEN (2, FILE=FNAM, STATUS='NEW')
      OPEN (6, FILE='InHouseSuccessFailure.txt', STATUS='NEW')
C     Perform data analysis.
      CALL READ
      CALL COMPUT
      IF (NCA .LT. NGE) CALL IMPUTE
      CALL AMMI
      CALL OUTPUT
      CLOSE (1)
      CLOSE (2)
      STOP
      END
      SUBROUTINE READ
C
C     SCRATCH memory:
C     Data            Vector     Offset             Length
C     REP data           Y       0                  GER
C     This item exists only if the data are replicated.
C
C     Read data file and check memory requirements.
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8     V(32000000)
      INTEGER*4  L(64000000)
      CHARACTER*1   GNA, ENA, RNA, SEQA, SEQB, SEQC, FT
      CHARACTER*3   EXDS
      CHARACTER*4   C(64000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD, FMT
      LOGICAL*4  EOF
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      COMMON  /PRMS/  YMV, NDTOT, NDTRT, NDGEN, NDENV, NDGXE,
     1 NDIPC(7), NDRES, NDERR, NDBWE, NDPER, NDEFT, SSTOT,
     2 SSTRT, SSGEN, SSENV, SSGXE, SSIPC(7), SSRES, SSERR,
     3 SSBWE, SSPER, SSEFT, FG, FE, FR, FGE, FGER, GM
      EQUIVALENCE  (V(1), Y(1), L(1), C(1))
C
C     The first line of the data file (FNDA) contains two items.
C     (1)  Columns  1-72, title for the data matrix.
C     (2)  Columns 77-80, brief name for the data matrix.
C
C     The second line contains thirteen items.
C     (1)  Columns  1-20, name for the data's measurement and units.
C     (2)  Columns 21-25, scaling integer N, multiplies data by 10**N.
C     (3)  Columns 26-30, number of matrix rows (such as genotypes).
C     (4)  Columns 32-35, name for matrix rows (such as GEN).
C     (5)  Columns 36-40, number of matrix columns (environments).
C     (6)  Columns 42-45, name for matrix columns (ENV).
C     (7)  Columns 46-50, number of replications.
C     (8)  Columns 52-55, name for replications (REP).
C     (9)  Columns 58-60, order of data matrix entries, specified by
C          the first letters of the subscript names (from columns
C          32, 42, and 52), the first changing fastest and the last
C          changing slowest. For example, for order 'RGE' the second
C          matrix entry is for REP 2, GEN 1, and ENV 1.
C     (10) Columns 62-65, name for treatments (such as TRT).
C     (11) Columns 68-70, experimental design ('RCB' or other).
C     (12) Columns 71-75, missing value indicator.
C     (13) Column 80, if 'N', no format supplied;
C          if 'F', format supplied for reading matrix entries.
C
C     The optional format line, given 'F' in column 80 of line 2,
C     contains the FORTRAN variable format (format A72).
C
C     The data matrix entries then follow, in free format if no format
C     is supplied, or in fixed format if 'F'.
C
C     Row (GEN) code names and then column (ENV) names follow
C     (format 15(A4, 1X)), each list starting a new line.
C     These are optional, but include either both or else neither.
C
C     For all of the above data, enter full lines completing all
C     items called for, even if this requires blanks or zeros
C     at the end of a line.
C
C     Read title card, matrix dimensions, and format if supplied.
      READ (1, '(A72, 4X, A4)')  NMD, NMB
      NMD = TRIM(NMD)
      NMB = ADJUSTR(NMB)
      READ (1, '(A20, I5, 3(I5,1X, A4), 2X, 3A1, 1X, A4, 2X, A3,
     1 F5.0, 4X, A1)')  NMU, KSCALE, NG, NMG, NE, NME, NR,
     2 NMR, SEQA, SEQB, SEQC, NMT, EXDS, YMV, FT
      IF ((NG .LT. 3) .OR. (NE .LT. 3))  THEN
      WRITE (6, '('' ERROR:  Termination because '', A4, '' or '',
     1 A4, '' is fewer than 3.'')')  NMG, NME
      GO TO 70
      ENDIF
      IF ((EXDS .EQ. 'RCB') .AND. (NR .LT. 2))  THEN
      WRITE (6, '('' ERROR:  Termination because RCB design '',
     1 ''requires 2 or more replications.'')')
      GO TO 70
      ENDIF
      KFRCB = 0
      IF (EXDS .EQ. 'RCB')  KFRCB = 1
      KFM = -1
      IF (FT .EQ. 'N'  .OR.  FT .EQ. 'n')  KFM = 0
      IF (FT .EQ. 'F'  .OR.  FT .EQ. 'f')  KFM = 1
      IF (KFM .EQ. -1)  THEN
      WRITE (6, '('' ERROR:  Termination because format ''
     1 ''specification '', A1, '' fails to be N or F.'')')  FT
      GO TO 70
      ELSEIF (KFM .EQ. 1)  THEN
      READ (1, '(A72)')  FMT
      FMT = TRIM(FMT)
      ENDIF
C     Set several constants.
      NX = MIN(NG-2, NE-2, 7)
      NGE = NG * NE
      NGER = NGE * NR
      NRG = NR * NG
      FG = NG
      FE = NE
      FR = NR
      FGE = NGE
      FGER = NGER
C     Calculate offsets.
      KDGE = MAX(NGE, ((NGER+1)/2), (6*NG)+NE)
      KBG = KDGE + NGE
      KBE = KBG + NG
      KSCG = KBE + NE
      KPGN = KSCG + NG
      KPEN = KPGN + (NG * NX)
      KNG = 2 * (KPEN + (NE * NX))
      KNE = KNG + NG
      KSGE = 4 * (KNE + NE)
C     Check memory, apart from needing M(KSGE) if empty cells exist.
      M4D = KNE + NE
      FMEM = DFLOAT(M4D) / DFLOAT(M4L)
      IF (M4D .GT. M4L)  THEN
      WRITE (6, '('' ERROR:  Termination because required memory '',
     1 ''exceeds allocated memory'', /, ''by a factor of'', F8.3,
     2 '' or somewhat more if there are any empty cells.'')')  FMEM
      GO TO 70
      ENDIF
C     Check variable names for uniqueness of first letters.
      GNA = NMG(1:1)
      ENA = NME(1:1)
      RNA = NMR(1:1)
      IF (.NOT. (GNA .NE. ENA  .AND.  GNA .NE. RNA  .AND.
     1 ENA .NE. RNA))  THEN
      WRITE (6, '('' ERROR:  Termination because first letters '',
     1 ''of variable names are not unique.'')')
      GO TO 70
      ENDIF
C     Check validity of input order.
C     SEQA changes fastest and SEQC slowest in the input data matrix
C     order, whereas the standardized storage order in Y(base) is to
C     increment REP fastest, GEN next, and ENV slowest, which is
C     the target order of RGE.  Likewise, in V(KDGE) increment
C     GEN fastest and ENV slowest, which is the target order of GER.
      KERR = 1
      IF     (SEQA.EQ.RNA .AND. SEQB.EQ.GNA .AND. SEQC.EQ.ENA)  THEN
      KERR = 0
      MA = 1
      MB = NR
      MC = NRG
      NA = NR
      NB = NG
      NC = NE
      ELSEIF (SEQA.EQ.RNA .AND. SEQB.EQ.ENA .AND. SEQC.EQ.GNA)  THEN
      KERR = 0
      MA = 1
      MB = NRG
      MC = NR
      NA = NR
      NB = NE
      NC = NG
      ELSEIF (SEQA.EQ.GNA .AND. SEQB.EQ.ENA .AND. SEQC.EQ.RNA)  THEN
      KERR = 0
      MA = NR
      MB = NRG
      MC = 1
      NA = NG
      NB = NE
      NC = NR
      ELSEIF (SEQA.EQ.GNA .AND. SEQB.EQ.RNA .AND. SEQC.EQ.ENA)  THEN
      KERR = 0
      MA = NR
      MB = 1
      MC = NRG
      NA = NG
      NB = NR
      NC = NE
      ELSEIF (SEQA.EQ.ENA .AND. SEQB.EQ.GNA .AND. SEQC.EQ.RNA)  THEN
      KERR = 0
      MA = NRG
      MB = NR
      MC = 1
      NA = NE
      NB = NG
      NC = NR
      ELSEIF (SEQA.EQ.ENA .AND. SEQB.EQ.RNA .AND. SEQC.EQ.GNA)  THEN
      KERR = 0
      MA = NRG
      MB = 1
      MC = NR
      NA = NE
      NB = NR
      NC = NG
      ENDIF
      IF (KERR .EQ. 1)  THEN
      WRITE (6, '('' ERROR:  Termination because specified input '',
     1 ''order '', 3A1, '' not valid.'')')  GNA, ENA, RNA
      GO TO 70
      ENDIF
C     If replicated data, read into Y(base).
      IF (NR .GT. 1)  THEN
      KYB = - NR - NRG
      IF (KFM .EQ. 0)  THEN
      READ (1, *)    (((Y(KYB+MC*NCD+MB*NBD+MA*NAD), NAD=1,NA),
     1 NBD=1,NB), NCD=1,NC)
      ELSE
      READ (1, FMT)  (((Y(KYB+MC*NCD+MB*NBD+MA*NAD), NAD=1,NA),
     1 NBD=1,NB), NCD=1,NC)
      ENDIF
C     Else if unreplicated data, read into V(KDGE).
      ELSE
      IF (SEQC .NE. RNA)  THEN
      WRITE (6, '('' ERROR:  Termination because unreplicated data '',
     1 ''should list reps last in the input order.'')')
      GO TO 70
      ENDIF
      IF (KFM .EQ. 0)  THEN
      IF (SEQA .EQ. GNA)  THEN
      READ (1, *)  (V(KDGE+ND), ND=1,NGE)
      ELSE
      READ (1, *)  ((V(KDGE+NGD+((NED-1)*NG)), NED=1,NE), NGD=1,NG)
      ENDIF
      ELSE
      IF (SEQA .EQ. GNA)  THEN
      READ (1, FMT)  (V(KDGE+ND), ND=1,NGE)
      ELSE
      READ (1, FMT)  ((V(KDGE+NGD+((NED-1)*NG)), NED=1,NE), NGD=1,NG)
      ENDIF
      ENDIF
      ENDIF
C     Read variable names, or else set to blanks if end-of-file.
      IF (.NOT. EOF(1))  THEN
      READ (1, '(15(A4,1X))')  (C(KNG+ND), ND=1,NG)
      READ (1, '(15(A4,1X))')  (C(KNE+ND), ND=1,NE)
      ELSE
      DO 100 ND=1,NG
  100 C(KNG+ND) = '    '
      DO 101 ND=1,NE
  101 C(KNE+ND) = '    '
      ENDIF
      RETURN
   70 CLOSE (6)
      STOP
      END
      SUBROUTINE COMPUT
C
C     SCRATCH memory:
C     TRT status M(KSGE) created if empty cells, NCA < NGE.
C
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8        V(32000000)
      INTEGER*4     L(64000000)
      CHARACTER*4   C(64000000)
      INTEGER*1     M(256000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT, SI
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      COMMON  /PRMS/  YMV, NDTOT, NDTRT, NDGEN, NDENV, NDGXE,
     1 NDIPC(7), NDRES, NDERR, NDBWE, NDPER, NDEFT, SSTOT,
     2 SSTRT, SSGEN, SSENV, SSGXE, SSIPC(7), SSRES, SSERR,
     3 SSBWE, SSPER, SSEFT, FG, FE, FR, FGE, FGER, GM
      COMMON  /FTST/  SSS, SSE, SMS, PF, NDS, NDE, SI
      EQUIVALENCE  (V(1), Y(1), L(1), C(1), M(1))
      DYMV = YMV
C     Scale data if KSCALE not 0.
      IF (KSCALE .NE. 0)  THEN
      IF (NR .GT. 1)  THEN
      YSCALE = 10.0 ** FLOAT(KSCALE)
      DO 100 ND=1,NGER
  100 IF (Y(ND) .NE. YMV)  Y(ND) = Y(ND) * YSCALE
      ELSE
      FSCALE = 10.0D0 ** DFLOAT(KSCALE)
      DO 101 ND=1,NGE
      NSB = KDGE + ND
  101 IF (V(NSB) .NE. DYMV)  V(NSB) = V(NSB) * FSCALE
      ENDIF
      ENDIF
C     If replicated data, average Y(base) into V(KDGE);
C     count NCA and NOA.
      NCA = 0
      NOA = 0
      IF (NR .GT. 1)  THEN
      NSB = KDGE
      NSS = 0
      DO 102 ND=1,NGE
      NSB = NSB + 1
      KNT = 0
      YAV = 0.0D0
      DO 103 N2=1,NR
      NSS = NSS + 1
      IF (Y(NSS) .EQ. YMV)  GO TO 103
      NOA = NOA + 1
      KNT = KNT + 1
      YAV = YAV + Y(NSS)
  103 CONTINUE
      IF (KNT .GE. 1)  THEN
      NCA = NCA + 1
      YAV = YAV / FLOAT(KNT)
      ELSE
      YAV = YMV
      ENDIF
  102 V(NSB) = YAV
C     Else if unreplicated data, count NCA and equate NOA.
      ELSE
      DO 104 ND=1,NGE
  104 IF (V(KDGE+ND) .NE. DYMV)  NCA = NCA + 1
      NOA = NCA
      ENDIF
C     Flag KFRCB = 1 or 0, but change to -1 if NOA < NGER.
      IF (NOA .LT. NGER)  KFRCB = - KFRCB
C     If empty cells, check memory and then create M(KSGE).
      IF (NCA .LT. NGE)  THEN
      M4D = M4D + ((NGE + 3) / 4)
      IF (M4D .GT. M4L)  THEN
      FMEM = DFLOAT(M4D) / DFLOAT(M4L)
      WRITE (6, '('' ERROR:  Termination because required memory '',
     1 ''exceeds allocated memory'', /, ''by a factor of'',
     2 F8.3, ''.'')')  FMEM
      GO TO 70
      ENDIF
      DO 105 ND=1,NGE
      KST = 1
      IF (V(KDGE+ND) .EQ. DYMV)  KST = 0
  105 M(KSGE+ND) = KST
      ENDIF
C     Check that each GEN and ENV has at least 3 observations.
      IF (NCA .LT. NGE)  THEN
      KERR = 0
      DO 106 NGD=1,NG
      KT = 0
      KS = KDGE + NGD
      DO 107 NED=1,NE
  107 IF (V(KS+((NED-1)*NG)) .NE. DYMV)  KT = KT + 1
      IF (KT .GE. 3)  GO TO 106
      KERR = 1
      WRITE (6, '('' ERROR:  Fewer than 3 observations for '',
     1 A4, I6, 2X, A4)')  NMG, NGD, C(KNG+NGD)
  106 CONTINUE
      KS = KDGE
      DO 108 NED=1,NE
      KT = 0
      DO 109 NGD=1,NG
      KS = KS + 1
  109 IF (V(KS) .NE. DYMV)  KT = KT + 1
      IF (KT .GE. 3)  GO TO 108
      KERR = 1
      WRITE (6, '('' ERROR:  Fewer than 3 observations for '',
     1 A4, I6, 2X, A4)')  NME, NED, C(KNE+NED)
  108 CONTINUE
      IF (KERR .EQ. 1)  THEN
      WRITE (6, '('' ERROR:  Termination because of '',
     1 ''insufficient data.'')')
      GO TO 70
      ENDIF
      ENDIF
C     If replicated, SSERR; set SSEFT.
      IF (NR .GT. 1)  THEN
      SSERR = 0.0D0
      NDERR = - NCA
      KS = 0
      DO 110 ND=1,NGE
      AV = V(KDGE+ND)
      DO 111 NRD=1,NR
      KS = KS + 1
      IF (Y(KS) .EQ. YMV)  GO TO 111
      R = DBLE(Y(KS)) - AV
      SSERR = SSERR + (R * R)
      NDERR = NDERR + 1
  111 CONTINUE
  110 CONTINUE
      SSEFT = SSERR
      NDEFT = NDERR
      ENDIF
C     If RCB design (and no missing data), SSBWE and SSPER;
C     set SSEFT using pure error if blocks/env significant at 0.05.
      IF (KFRCB .EQ. 1)  THEN
      SSBWE = 0.0D0
      NGR = NG * NR
      FGR = NGR
      DO 112 NED=1,NE
      EAV = 0.0D0
      DO 113 NGRD=1,NGR
  113 EAV = EAV + DBLE(Y(((NED-1)*NGR)+NGRD))
      EAV = EAV / FGR
      DO 114 NRD=1,NR
      AV = 0.0D0
      DO 115 NGD=1,NG
      KS = ((NED - 1) * NGR) + ((NGD - 1) * NR) + NRD
  115 AV = AV + DBLE(Y(KS))
      AV = AV / FG
      X = AV - EAV
  114 SSBWE = SSBWE + (X * X)
  112 CONTINUE
      SSBWE = SSBWE * FG
      NDBWE = NE * (NR - 1)
      SSPER = SSERR - SSBWE
      NDPER = NDERR - NDBWE
      SSS = SSBWE
      NDS = NDBWE
      SSE = SSPER
      NDE = NDPER
      CALL FTEST
      IF (PF .LE. 0.05D0)  THEN
      SSEFT = SSPER
      NDEFT = NDPER
      ENDIF
      ENDIF
C     Initialize GEN IPCs with arbitrary trial vectors.
      IX = 20412
      IY = 28357
      IZ = 12275
      DO 116 NGD=1,NG
      IX = MOD(170*IX, 30323)
      IY = MOD(171*IY, 30269)
      IZ = MOD(172*IZ, 30307)
      RN = AMOD(FLOAT(IX)/30323.+FLOAT(IY)/30269.+FLOAT(IZ)/30307., 1.)
      R = RN - 0.5D0
      V(KPGN+NGD) = R
      V(KSCG+NGD) = R
      DO 117 NXD=1,NX
      NSUB = ((NXD - 1) * NG) + NGD
  117 V(KPGN+NSUB) = R
  116 CONTINUE
      WRITE (6, '(''Success'')')
      CLOSE (6)
      RETURN
   70 CLOSE (6)
      STOP
      END
      SUBROUTINE IMPUTE
C
C     SCRATCH memory:
C     Data            Vector     Offset             Length
C     IPC scratch        V       0                  GE
C
C     TRT status M(KSGE) exists whenever this subroutine is called.
C
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8        V(32000000)
      INTEGER*4     L(64000000)
      CHARACTER*4   C(64000000)
      INTEGER*1     M(256000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      COMMON  /PRMS/  YMV, NDTOT, NDTRT, NDGEN, NDENV, NDGXE,
     1 NDIPC(7), NDRES, NDERR, NDBWE, NDPER, NDEFT, SSTOT,
     2 SSTRT, SSGEN, SSENV, SSGXE, SSIPC(7), SSRES, SSERR,
     3 SSBWE, SSPER, SSEFT, FG, FE, FR, FGE, FGER, GM     
      EQUIVALENCE  (V(1), Y(1), L(1), C(1), M(1))
C     Initialize missing cells with AMMI0 model.
      GMU = 0.0D0
      KS = 0
      DO 100 NED=1,NE
      AV = 0.0D0
      NV = 0
      DO 101 NGD=1,NG
      KS = KS + 1
      IF (M(KSGE+KS) .EQ. 0)  GO TO 101
      X = V(KDGE+KS)
      GMU = GMU + X
      AV = AV + X
      NV = NV + 1
  101 CONTINUE
      R = AV / DFLOAT(NV)
  100 V(KBE+NED) = R
      GMU = GMU / DFLOAT(NCA)
      DO 102 NGD=1,NG
      AV = 0.0D0
      NV = 0
      DO 103 NED=1,NE
      KS = ((NED - 1) * NG) + NGD
      IF (M(KSGE+KS) .EQ. 0)  GO TO 103
      AV = AV + V(KDGE+KS)
      NV = NV + 1
  103 CONTINUE
  102 V(KBG+NGD) = AV / DFLOAT(NV)
      KS = 0
      DO 104 NED=1,NE
      DO 105 NGD=1,NG
      KS = KS + 1
  105 IF (M(KSGE+KS) .EQ. 0)
     1 V(KDGE+KS) = V(KBG+NGD) + V(KBE+NED) - GMU
  104 CONTINUE
C     Compute AMMI0 RMS; set VSFT.
      RMS = 0.0D0
      NS = KDGE
      DO 160 NED=1,NE
      R = 0.0D0
      DO 161 NGD=1,NG
      NS = NS + 1
      X = V(NS)
  161 R = R + (X * X)
  160 RMS = RMS + R
      RMS = DSQRT(RMS / FGE)
      VSFT = 1.0D-6 * RMS
C     Impute missing cells with EM-AMMI1 model.
      DO 110 NCY=1,100
      GM = 0.0D0
      KS = 0
      DO 120 NED=1,NE
      AV = 0.0D0
      DO 121 NGD=1,NG
      KS = KS + 1
  121 AV = AV + V(KDGE+KS)
      R = AV / FG
      V(KBE+NED) = R
  120 GM = GM + R
      GM = GM / FE
      DO 122 NGD=1,NG
      AV = 0.0D0
      DO 123 NED=1,NE
      KS = ((NED - 1) * NG) + NGD
  123 AV = AV + V(KDGE+KS)
  122 V(KBG+NGD) = AV / FE
      KS = 0
      DO 124 NED=1,NE
      X = V(KBE+NED) - GM
      DO 125 NGD=1,NG
      KS = KS + 1
  125 V(KS) = V(KDGE+KS) - V(KBG+NGD) - X
  124 CONTINUE
      DO 130 ITER=1,100
      EV = 0.0D0
      KS = 0
      DO 140 NED=1,NE
      X = 0.0D0
      DO 141 NGD=1,NG
      KS = KS + 1
  141 X = X + (V(KS) * V(KPGN+NGD))
      EV = EV + (X * X)
  140 V(KPEN+NED) = X
      SV = DSQRT(EV)
      DO 142 NED=1,NE
      KS = KPEN + NED
  142 V(KS) = V(KS) / SV
      EV = 0.0D0
      DO 143 NGD=1,NG
      X = 0.0D0
      DO 144 NED=1,NE
      KS = ((NED - 1) * NG) + NGD
  144 X = X + (V(KS) * V(KPEN+NED))
      EV = EV + (X * X)
  143 V(KPGN+NGD) = X
      SV = DSQRT(EV)
      SFT = 0.0D0
      DO 145 NGD=1,NG
      X = V(KPGN+NGD) / SV
      V(KPGN+NGD) = X
      KS = KSCG + NGD
      R = DABS(X - V(KS))
      IF (R .GT. SFT)  SFT = R
  145 V(KS) = X
      IF (SFT .LE. 1.0D-9)  GO TO 10
  130 CONTINUE
   10 CONTINUE
      VMIS = 0.0D0
      KS = 0
      DO 150 NED=1,NE
      DO 151 NGD=1,NG
      KS = KS + 1
      IF (M(KSGE+KS) .GT. 0)  GO TO 151
      R = V(KBG+NGD) + V(KBE+NED) - GM + (SV*V(KPGN+NGD)*V(KPEN+NED))
      X = DABS(R - V(KDGE+KS))
      IF (X .GT. VMIS)  VMIS = X
      V(KDGE+KS) = R
  151 CONTINUE
  150 CONTINUE
      IF (VMIS .LE. VSFT)  GO TO 11
  110 CONTINUE
      WRITE (2, '(''Warning:  Impute missing cells failed to '',
     1 ''converge by limit of 100 iterations.'', /,
	2 ''Largest shift'', F12.8, '' of AMMI0 RMS exceeded '',
	3 ''criterion of 0.000001.'', /, /)')  VMIS / RMS
   11 CONTINUE
      RETURN
      END
      SUBROUTINE AMMI
C
C     SCRATCH memory:
C     Data            Vector     Offset             Length
C     IPC scratch        V       0                  GE
C
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8        V(32000000)
      INTEGER*4     L(64000000)
      CHARACTER*4   C(64000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT, SI
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      COMMON  /PRMS/  YMV, NDTOT, NDTRT, NDGEN, NDENV, NDGXE,
     1 NDIPC(7), NDRES, NDERR, NDBWE, NDPER, NDEFT, SSTOT,
     2 SSTRT, SSGEN, SSENV, SSGXE, SSIPC(7), SSRES, SSERR,
     3 SSBWE, SSPER, SSEFT, FG, FE, FR, FGE, FGER, GM
      COMMON  /FTST/  SSS, SSE, SMS, PF, NDS, NDE, SI    
      EQUIVALENCE  (V(1), Y(1), L(1), C(1))
C     ANOVA PART OF AMMI.
C     Compute GEN, ENV, and grand means; SSGEN and SSENV.
      GM = 0.0D0
      DO 100 NGD=1,NG
      AV = 0.0D0
      KS = KDGE + NGD
      DO 101 NED=1,NE
  101 AV = AV + V(KS+((NED-1)*NG))
      GM = GM + AV
  100 V(KBG+NGD) = AV / FE
      GM = GM / FGE
      KS = KDGE
      DO 102 NED=1,NE
      AV = 0.0D0
      DO 103 NGD=1,NG
      KS = KS + 1
  103 AV = AV + V(KS)
  102 V(KBE+NED) = AV / FG
      SSGEN = 0.0D0
      DO 104 NGD=1,NG
      R = V(KBG+NGD) - GM
  104 SSGEN = SSGEN + (R * R)
      SSGEN = SSGEN * FE * FR
      NDGEN = NG - 1
      SSENV = 0.0D0
      DO 105 NED=1,NE
      R = V(KBE+NED) - GM
  105 SSENV = SSENV + (R * R)
      SSENV = SSENV * FG * FR
      NDENV = NE - 1
C     Put GxE into V(base); SSGXE, SSTRT, and SSTOT.
      SSGXE = 0.0D0
      KS = 0
      DO 106 NED=1,NE
      DO 107 NGD=1,NG
      KS = KS + 1
      R = V(KDGE+KS) - V(KBG+NGD) - V(KBE+NED) + GM
      SSGXE = SSGXE + (R * R)
  107 V(KS) = R
  106 CONTINUE
      SSGXE = SSGXE * FR
      NDGXE = NDGEN * NDENV
      SSTRT = SSGEN + SSENV + SSGXE
      NDTRT = NGE - 1
      SSTOT = SSTRT + SSERR
      NDTOT = NDTRT + NDERR
C     IPC PART OF AMMI.
C     Calculate IPC scores for axes 1 to NX with the power method.
      SSRES = SSGXE
      NDRES = NDGXE
      DO 120 NXD=1,NX
      KG = KPGN + ((NXD - 1) * NG)
      KE = KPEN + ((NXD - 1) * NE)
      DO 121 ITER=1,100
C     Compute ENV scores from GEN scores.
      EV = 0.0D0
      KS = 0
      DO 122 NED=1,NE
      X = 0.0D0
      DO 123 NGD=1,NG
      KS = KS + 1
  123 X = X + (V(KS) * V(KG+NGD))
      EV = EV + (X * X)
  122 V(KE+NED) = X
      SV = DSQRT(EV)
      DO 124 NED=1,NE
      KS = KE + NED
  124 V(KS) = V(KS) / SV
C     Compute GEN scores from ENV scores.
      EV = 0.0D0
      DO 125 NGD=1,NG
      X = 0.0D0
      DO 126 NED=1,NE
      KS = ((NED - 1) * NG) + NGD
  126 X = X + (V(KS) * V(KE+NED))
      EV = EV + (X * X)
  125 V(KG+NGD) = X
      SV = DSQRT(EV)
C     Test convergence and retain old GEN scores as unit vector.
      SFT = 0.0D0
      DO 127 NGD=1,NG
      KS = KG + NGD
      X = V(KS) / SV
      V(KS) = X
      KS = KSCG + NGD
      R = DABS(X - V(KS))
      IF (R .GT. SFT)  SFT = R
  127 V(KS) = X
      IF (SFT .LE. 1.0D-9)  GO TO 10
  121 CONTINUE
C     Scale scores and compute residuals; SSIPC(NXD).
   10 SSIPC(NXD) = EV * FR
      NDIPC(NXD) = NG + NE - (2 * NXD) - 1
      SSRES = SSRES - SSIPC(NXD)
      NDRES = NDRES - NDIPC(NXD)
      R = DSQRT(SV)
      DO 128 NGD=1,NG
      KS = KG + NGD
  128 V(KS) = V(KS) * R
      KS = 0
      DO 129 NED=1,NE
      NS = KE + NED
      X = V(NS) * R
      V(NS) = X
      DO 130 NGD=1,NG
      KS = KS + 1
  130 V(KS) = V(KS) - (X * V(KG+NGD))
  129 CONTINUE
  120 CONTINUE
C     If corr(GEN deviations, GEN IPC1) < 0, reverse polarity.
      X = 0.0D0
      DO 131 NGD=1,NG
  131 X = X + ((V(KBG+NGD) - GM) * V(KPGN+NGD))
      IF (X .LT. 0.0D0)  THEN
      DO 132 NGD=1,NG
  132 V(KPGN+NGD) = - V(KPGN+NGD)
      DO 133 NED=1,NE
  133 V(KPEN+NED) = - V(KPEN+NED)
      ENDIF
      RETURN
      END
      SUBROUTINE OUTPUT
C
C     Data            Vector     Offset             Length
C
C     SCRATCH memory for entire subroutine:
C     GEN scratch        V       0                  G
C
C     GEN IPC1 ranks     L       KRG=2*G            G
C     ENV IPC1 ranks     L       KRE=KRG+G          E
C
C     SCRATCH memory for winners for AMMI family:
C     GEN counts NX      L       KXG=KRE+E          GN
C     GEN counts F       L       KFG=KXG+GN         G
C
C     SCRATCH memory for ranking table:
C     GEN AMMI1 ranks    L       KR1G=KRE+E         G
C     GEN AMMIF ranks    L       KRFG=KR1G+G        G
C     GEN switch         L       KSWG=KRFG+G        G
C
C     SCRATCH memory for GEN list of means, IPC1, and IPC2.
C     GEN mean ranks     L       KRMG=KRE+E         G
C
C     SCRATCH memory for ENV list of means, IPC1, and IPC2.
C     ENV mean ranks     L       KRME=KRE+E         E
C
C     TRT status M(KSGE) exists if empty cells, NCA < NGE.
C
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8        V(32000000)
      INTEGER*4     L(64000000),    MEC(9)
      CHARACTER*4   C(64000000)
      INTEGER*1     M(256000000)
      CHARACTER*4   NMB, NMG, NME, NMR, NMT, SI
      CHARACTER*5   CL(9), CM(9)
      CHARACTER*20  NMU
      CHARACTER*160  FNDA, FNAM
      CHARACTER*72  NMD
      COMMON  /DATA/  Y(64000000)
      COMMON  /PRML/  NG, NE, NR, NGE, NGER, NX, KDGE, KBG, KBE, KSCG,
     1 KPGN, KPEN, KNG, KNE, KSGE, M4L, M4D, KSCALE, NOA, NCA, KFRCB
      COMMON  /PRMC/  NMB, NMG, NME, NMR, NMT, NMD, NMU, FNDA, FNAM
      COMMON  /PRMS/  YMV, NDTOT, NDTRT, NDGEN, NDENV, NDGXE,
     1 NDIPC(7), NDRES, NDERR, NDBWE, NDPER, NDEFT, SSTOT,
     2 SSTRT, SSGEN, SSENV, SSGXE, SSIPC(7), SSRES, SSERR,
     3 SSBWE, SSPER, SSEFT, FG, FE, FR, FGE, FGER, GM
      COMMON  /FTST/  SSS, SSE, SMS, PF, NDS, NDE, SI   
      EQUIVALENCE  (V(1), Y(1), L(1), C(1), M(1))
C     Calculate offsets.
      KRG = 2 * NG
      KRE = KRG + NG
      KXG = KRE + NE
      KFG = KXG + (NG * NX)
      KR1G = KRE + NE
      KRFG = KR1G + NG
      KSWG = KRFG + NG
      KRMG = KRE + NE
      KRME = KRE + NE
C     Print basic facts about the dataset.
      WRITE (2, '(''AMMISOFT for '', A4,
     1 '':  Data Input'')')  NMB
      WRITE (2, '(/, A72)')  NMD
      WRITE (2, '(/, ''Data input filename:'', /, A80)')  FNDA(1:80)
      IF (LEN(TRIM(FNDA)) .GE. 81)  WRITE (2, '(A80)')  FNDA(81:160)
      WRITE (2, '(/, ''AMMI output filename:'', /, A80)')  FNAM(1:80)
      IF (LEN(TRIM(FNAM)) .GE. 81)  WRITE (2, '(A80)')  FNAM(81:160)
      WRITE (2, '(/, ''Measurement  '', A20)')  NMU
      WRITE (2, '(''Number of '', A4, I9)')  NMG, NG
      WRITE (2, '(''Number of '', A4, I9)')  NME, NE
      WRITE (2, '(''Number of '', A4, I9)')  NMR, NR
      WRITE (2, '(/, ''Number of '', A4, '' with data'', I10,
     1 '' out of'', I10, '' possible.'')')  NMT, NCA, NGE
      IF (NCA .EQ. NGE)  THEN
      WRITE (2, '(''No empty cells.'')')
      ELSE
      WRITE (2, '(''EM-AMM1 imputed'', I10, '' missing cells.'',
     1 /, ''ANOVA treats all matrix cells as having all'', I4,
     2 '' replications.'')')  NGE - NCA, NR
      ENDIF
      IF (NR .GT. 1)  THEN
      WRITE (2, '(/, ''Number of observations'', I12, '' out of'',
     1 I10, '' possible.'')')  NOA, NGER
      IF (NOA .EQ. NGER)  THEN
      WRITE (2, '(''No missing data; these data are balanced.'')')
      ELSE
      WRITE (2, '(''ANOVA treats all matrix cells as having all'',
     1 I4, '' replications.'')')  NR
      ENDIF
      ENDIF
      IF (KFRCB .EQ. -1)  WRITE (2, '(/, ''RCB design ignored '',
     1 ''because of missing observations.'')')
      FMEM = DFLOAT(M4D) / DFLOAT(M4L)
      WRITE (2, '(/, ''Required memory is '', F9.4,
     1 ''% of allocated memory.'')')  100.0D0 * FMEM
      WRITE (2, '(/, ''Grand mean'', F15.5)')  GM
      IF (NR .EQ. 1)  THEN
      IF (NCA .LT. NGE)  WRITE (2,
     1 '(''(Includes imputed values.)'')')
      ELSE
      IF ((NCA .EQ. NGE) .AND. (NOA .LT. NGER))  WRITE (2,
     1 '(''(Unweighted mean.)'')')
      IF ((NCA .LT. NGE) .AND. (NOA .EQ. NCA*NR))  WRITE (2,
     1 '(''(Includes imputed values.)'')')
      IF ((NCA .LT. NGE) .AND. (NOA .LT. NCA*NR))  WRITE (2,
     1 '(''(Unweighted mean; includes imputed values.)'')')
      ENDIF
C     Print ANOVA table.
      WRITE (2, '(/, /, ''AMMISOFT for '', A4,
     1 '':  ANOVA Table'')')  NMB
      IF (NR .GT. 1)  THEN
      WRITE (2, '(/, 80(''-''), /, ''Source'', 13X, ''df'', 15X,
     1 ''SS'', 18X, ''MS'', 10X, ''Probability'', /, 80(''-''))')
      WRITE (2, '(''Total'', I17, F22.5, F20.5)')
     1 NDTOT, SSTOT, SSTOT / DBLE(NDTOT)
      SSE = SSEFT
      NDE = NDEFT
      SSS = SSTRT
      NDS = NDTRT
      CALL FTEST
      WRITE (2, '(2X, A4, I16, F22.5, F20.5, F12.7, A4)')
     1 NMT, NDS, SSS, SMS, PF, SI
      SSS = SSGEN
      NDS = NDGEN
      CALL FTEST
      WRITE (2, '(4X, A4, I14, F22.5, F20.5, F12.7, A4)')
     1 NMG, NDS, SSS, SMS, PF, SI
      SSS = SSENV
      NDS = NDENV
      CALL FTEST
      WRITE (2, '(4X, A4, I14, F22.5, F20.5, F12.7, A4)')
     1 NME, NDS, SSS, SMS, PF, SI
      SSS = SSGXE
      NDS = NDGXE
      CALL FTEST
      WRITE (2, '(4X, A1, ''x'', A1, I15, F22.5, F20.5, F12.7, A4)')
     1 NMG(1:1), NME(1:1), NDS, SSS, SMS, PF, SI
      MFR = 0
      INC = 1
      SSS = SSGXE
      EMS = SSEFT / DFLOAT(NDEFT)
      SSGXEN = DFLOAT(NDGXE) * EMS
      SSGXES = SSGXE - SSGXEN
      SSUM = 0.0D0
      SSCAP = 0.0D0
      NSN = 0
      DO 100 NXD=1,NX
      NDS = (NG - NXD) * (NE - NXD)
      IF (NXD .GT. 1)  SSS = SSS - SSIPC(NXD - 1)
      CALL FTEST
      IF (PF .LE. 0.01D0)  THEN
      MFR = MFR + INC
      ELSE
      INC = 0
      ENDIF
      R = SSIPC(NXD)
      SSUM = SSUM + R
      IF (SSGXES .GE. SSUM)  THEN
      NSN = NXD
      SSCAP = SSCAP + R
      ENDIF
      NDFX = NDIPC(NXD)
      WRITE (2, '(6X, ''IPC'', I1, I12, F22.5, F20.5, F12.7, A4)')
     1 NXD, NDFX, R, R / DFLOAT(NDFX), PF, SI
  100 CONTINUE
      SSE = SSEFT
      NDE = NDEFT
      SSS = SSRES
      NDS = NDRES
      CALL FTEST
      WRITE (2, '(6X, ''Residual'', I8, F22.5, F20.5, F12.7, A4)')
     1 NDRES, SSRES, SSRES / DFLOAT(NDRES), PF, SI
      WRITE (2, '(2X, ''Error'', I15, F22.5, F20.5)')
     1 NDERR, SSERR, SSERR / DFLOAT(NDERR)
      IF (KFRCB .EQ. 1)  THEN
      SSE = SSPER
      NDE = NDPER
      SSS = SSBWE
      NDS = NDBWE
      CALL FTEST
      WRITE (2, '(4X, ''Blocks/Env'', I8, F22.5, F20.5, F12.7,
     1 A4)')  NDS, SSS, SMS, PF, SI
      WRITE (2, '(4X, ''Pure Error'', I8, F22.5, F20.5)')
     1 NDPER, SSPER, SSPER / DFLOAT(NDPER)
      ENDIF
      WRITE (2, '(80(''-''))')
      IF (KFRCB .EQ. 1)  THEN
      IF (PF .LE. 0.05D0)  THEN
      WRITE (2, '(''F-tests use Pure Error because Blocks/Env '',
     1 ''are significant at the 0.05 level.'')')
      ELSE
      WRITE (2, '(''F-tests use Error because Blocks/Env '',
     1 ''are NOT significant at the 0.05 level.'')')
      ENDIF
      ENDIF
      ELSE
      WRITE (2, '(/, 63(''-''), /, ''Source'', 12X, ''df'', 15X,
     1 ''SS'', 18X, ''MS'', /, 63(''-''))')
      WRITE (2, '(A4, I17, F22.5, F20.5)')
     1 NMT, NDTRT, SSTRT, SSTRT / DFLOAT(NDTRT)
      WRITE (2, '(2X, A4, I15, F22.5, F20.5)')
     1 NMG, NDGEN, SSGEN, SSGEN / DFLOAT(NDGEN)
      WRITE (2, '(2X, A4, I15, F22.5, F20.5)')
     1 NME, NDENV, SSENV, SSENV / DFLOAT(NDENV)
      WRITE (2, '(2X, A1, ''x'', A1, I16, F22.5, F20.5)')
     1 NMG(1:1), NME(1:1), NDGXE, SSGXE, SSGXE / DFLOAT(NDGXE)
      DO 101 NXD=1,NX
      R = SSIPC(NXD)
      NDFX = NDIPC(NXD)
      WRITE (2, '(4X, ''IPC'', I1, I13, F22.5, F20.5)')
     1 NXD, NDFX, R, R / DFLOAT(NDFX)
  101 CONTINUE
      WRITE (2, '(4X, ''Residual'', I9, F22.5, F20.5)')
     1 NDRES, SSRES, SSRES / DFLOAT(NDRES)
      WRITE (2, '(63(''-''))')
      ENDIF
C     Print commentary on model diagnosis.
      IF (NR .GT. 1)  THEN
      WRITE (2, '(/, ''MODEL DIAGNOSIS:'')')
      WRITE (2, '(/, ''AMMI comprises a model family with AMMI0, '',
     1 ''AMMI1, AMMI2 and so on'', /, ''retaining 0, 1, 2, or more '',
     2 ''interaction principal components (IPCs)'', /, 
     3 ''before relegating higher components to a discarded '',
     4 ''residual;'', /, ''finally, the full model AMMIF equals '',
     5 ''the actual data (or averages'', /, ''over reps), so it '',
     6 ''has no residual.  Model diagnosis is required'', /,
     7 ''to determine the best AMMI model for a given dataset, '',
     8 ''based on'', /, ''statistical and practical ''
     9 ''considerations.'')')
      WRITE (2, '(/, ''FR-tests at the 0.01 level diagnose AMMI'',
     1 I1, ''.'')')  MFR
      WRITE (2, '(/, ''Estimated sums of squares for GxE signal '',
     1 ''and noise:'')')
      WRITE (2, '(A1, ''x'', A1, '' total'', F25.5)')
     1 NMG(1:1), NME(1:1), SSGXE
      R = 100.0D0 * (SSGXEN / SSGXE)
      WRITE (2, '(A1, ''x'', A1, '' noise'', F25.5, '' or'',
     1 F7.2, ''%'')')  NMG(1:1), NME(1:1), SSGXEN, R
      R = 100.0D0 * (SSGXES / SSGXE)
      WRITE (2, '(A1, ''x'', A1, '' signal'', F24.5, '' or'',
     1 F7.2, ''%'')')  NMG(1:1), NME(1:1), SSGXES, R
      WRITE (2, '(''Early IPCs selectively capture signal, '',
     1 ''and late ones noise.'')')
      KFCON = 0
      IF (MFR .EQ. NSN)  KFCON = 1
      IF (NSN .EQ. NX)  THEN
      WRITE (2, '(''Accordingly, this much signal suggests '',
     1 ''AMMI'', I1, ''.'')')  NX
      ELSE
      X = (SSGXES - SSCAP) / SSIPC(NSN+1)
      IF (X .LE. 0.3333333333D0)  THEN
      WRITE (2, '(''Accordingly, this much signal suggests '',
     1 ''AMMI'', I1, ''.'')')  NSN
      ELSEIF ((X .GT. 0.3333333D0) .AND. (X .LT. 0.6666667D0))  THEN
      WRITE (2, '(''Accordingly, this much signal suggests '',
     1 ''AMMI'', I1, '' or maybe AMMI'', I1, ''.'')')  NSN, NSN+1
      IF (MFR .EQ. NSN+1)  KFCON = 1
      ELSEIF (X .GE. 0.6666667D0)  THEN
      WRITE (2, '(''Accordingly, this much signal suggests '',
     1 ''AMMI'', I1, '' or maybe AMMI'', I1, ''.'')')  NSN+1, NSN
      IF (MFR .EQ. NSN+1)  KFCON = 1
      ENDIF
      ENDIF
      WRITE (2, '(/, ''Conclusion:'')')
      IF (KFCON .EQ. 1)  THEN
      WRITE (2, '(''The consensus diagnosis is AMMI'', I1,
     1 '' for this dataset.'')')  MFR
      IF (MFR .EQ. 0)  WRITE (2, '(''The '', A1, ''x'', A1,
     1 '' interaction is mostly noise, so AMMI analysis is '',
     2 ''inappropriate.'')')  NMG(1:1), NME(1:1)
      IF (MFR .EQ. 7)  WRITE (2, '(''(But note that AMMISOFT '',
     1 ''is limited to 7 IPCs.)'')')
      NREC = KFCON
      ELSE
      NLC = MIN(MFR, NSN)
      NUC = MAX(MFR, NSN)
      NREC = NLC
      WRITE (2, '(''Lacking consensus, parsimony favors the '',
     1 ''simpler AMMI'', I1, '' for this dataset,'', /,
     2 ''although AMMI'', I1, '' also merits consideration.'')')
     3  NLC, NUC
      IF (NLC .EQ. 0)  WRITE (2, '(''The '', A1, ''x'', A1,
     1 '' interaction might be mostly noise, so AMMI analysis '',
     2 ''is questionable.'')')  NMG(1:1), NME(1:1)
      ENDIF
      IF (NREC .GE. 2)  WRITE (2, '(''However, practical '',
     1 ''constraints that necessitate few mega-environments'', /,
     2 ''may justify selection of a model lower than AMMI'', I1,
     3 '', perhaps AMMI1.'')')  NREC
      WRITE (2, '(/, ''Perspective:'')')
      R = SSGXES / SSGEN
      WRITE (2, '(''Note that the SS for '', 2A1, ''-signal is'',
     1 F6.2, '' times that for '', A4, '' main effects.'')')
     2 NMG(1:1), NME(1:1), R, NMG
      R = SSIPC(1) / SSGEN
      IF (R .GE. 0.5D0)  THEN
      WRITE (2, '(''Hence, narrow adaptations are important '',
     1 ''for this dataset.'')')
      IF (SSGXES .GE. SSIPC(1))  WRITE (2, '(''Even just IPC1 '',
     2 ''alone is'', F6.2, '' times the '', A4,
     3 '' main effects.'')')  R, NMG
      ENDIF
      R = SSGXEN / SSGEN
      WRITE (2, '(''Also note that '', 2A1, ''-noise is'', F6.2,
     1 '' times the '', A4, '' main effects.'')')  NMG(1:1),
     2 NME(1:1), R, NMG
      WRITE (2, '(''Discarding noise improves accuracy, '',
     1 ''increases repeatability,'', /, ''simplifies '',
     2 ''conclusions, and accelerates progress.'')')
      ELSE
      WRITE (2, '(/, ''MODEL DIAGNOSIS for this dataset cannot '',
     1 ''consider FR-tests or GxE signal'', /, ''because of no '',
     2 ''replication.  Instead, consider which AMMI model '',
     3 ''yields a'', /, ''practical number of mega-environments, '',
     4 ''as well as which model parameters'', /, ''have an '',
     5 ''evident biological or ecological interpretation.'')')
      ENDIF
C     Print winners for AMMI model family; mega-environment counts.
      CALL INDEX(NG, KPGN, KRG)
      CALL INDEX(NE, KPEN, KRE)
      DO 110 NLD=1,9
  110 CL(NLD) = '-----'
      DO 111 NGXD=1,NG*NX
  111 L(KXG+NGXD) = 0
      DO 112 NGD=1,NG
  112 L(KFG+NGD) = 0
      VTEST = 0.0D0
      DO 113 NGD=1,NG
      R = V(KBG+NGD)
      IF (R .LE. VTEST)  GO TO 113
      VTEST = R
      NGMW = NGD
  113 CONTINUE
      DO 114 NED=1,NE
      DO 115 NGD=1,NG
  115 V(NGD) = V(KBG+NGD)
      DO 116 NXD=1,NX
      NXD1 = NXD - 1
      NSG = KPGN + (NXD1 * NG)
      NSE = KPEN + (NXD1 * NE) + NED
      VXE = V(NSE)
      NSC = KXG + (NXD1 * NG)
      VTEST = 0.0D0
      DO 117 NGD=1,NG
      R = V(NGD) + (V(NSG+NGD) * VXE)
      V(NGD) = R
      IF (R .LE. VTEST)  GO TO 117
      VTEST = R
      NGWIN = NGD
  117 CONTINUE
  116 L(NSC+NGWIN) = L(NSC+NGWIN) + 1
  114 CONTINUE
      DO 118 NED=1,NE
      VTEST = 0.0D0
      NOE = (NED - 1) * NG
      DO 119 NGD=1,NG
      R = V(KDGE+NOE+NGD)
      IF (R .LE. VTEST)  GO TO 119
      VTEST = R
      NGWIN = NGD
  119 CONTINUE
  118 L(KFG+NGWIN) = L(KFG+NGWIN) + 1
      DO 120 NXD=1,NX
      MCNT = 0
      NSC = KXG + ((NXD - 1) * NG)
      DO 121 NGD=1,NG
  121 IF(L(NSC+NGD) .GE. 1)  MCNT = MCNT + 1
  120 MEC(NXD+1) = MCNT
      MCNT = 0
      DO 122 NGD=1,NG
  122 IF (L(KFG+NGD) .GE. 1)  MCNT = MCNT + 1
      MEC(NX+2) = MCNT
      NX2 = NX + 2
      WRITE (2, '(/, /, ''AMMISOFT for '', A4, '':  Winners for '',
     1 ''AMMI Model Family'')')  NMB
      WRITE (2, '(/, 14(''-''), 9A5)')  (CL(NXD), NXD=1,NX2)
      WRITE (2, '(14X, ''AMMI Model Family'')')
      WRITE (2, '(14X, 9A5)')  (CL(NXD), NXD=1,NX2)
      DO 123 NXD=1,NX+1
  123 WRITE (CM(NXD), '(I5)')  NXD - 1
      CM(NX2) = '    F'
      WRITE (2, '(''Genotype      '', 9A5)')  (CM(NXD), NXD=1,NX2)
      WRITE (2, '(14(''-''), 9A5)')  (CL(NXD), NXD=1,NX2)
      NNW = 0
      DO 124 NGD=1,NG
      NGR = L(KRG+NGD)
      NLP = 0
      IF (NGR .EQ. NGMW)  THEN
      NLP = 1
      WRITE (CM(1), '(I5)')  NE
      ELSE
      CM(1) = '     '
      ENDIF
      DO 125 NXD=1,NX
      NXD1 = NXD - 1
      NXBG = L(KXG+(NXD1*NG)+NGR)
      IF (NXBG .GE. 1)  THEN
      NLP = 1
      WRITE (CM(NXD+1), '(I5)')  NXBG
      ELSE
      CM(NXD+1) = '     '
      ENDIF
  125 CONTINUE
      NXBG = L(KFG+NGR)
      IF (NXBG .GE. 1)  THEN
      NLP = 1
      WRITE (CM(NX2), '(I5)')  NXBG
      ELSE
      CM(NX2) = '     '
      ENDIF
      IF (NLP .EQ. 1)  THEN
      WRITE (2, '(I5, 2X, A4, 3X, 9A5)')  NGR, C(KNG+NGR),
     1 (CM(NXD), NXD=1,NX2)
      ELSE
      NNW = NNW + 1
      ENDIF
  124 CONTINUE
      WRITE (2, '(14(''-''), 9A5)')  (CL(NXD), NXD=1,NX2)
      WRITE (2, '(''Mega-environments 1'', 8I5)')
     1 (MEC(NXD), NXD=2,NX2)
      WRITE (2, '(14(''-''), 9A5)')  (CL(NXD), NXD=1,NX2)
      IF (NNW .GE. 1)  WRITE (2, '(''Another'', I5,
     1 '' genotypes never win, so they are not listed.'')')  NNW
      WRITE (2, '(/, ''     Mega-environments are distinguished by '',
     1 ''having different genotype'', /, ''winners.  Increasingly '',
     2 ''complex AMMI models generally have more genotype'', /,
     3 ''winners or m-es, as shown in the list at the bottom of the '',
     4 ''above table.'', /, ''This is an important reason why model '',
     5 ''diagnosis matters.'')')
      WRITE (2, '(''     Mega-environment analysis has been '',
     1 ''applied mostly to yield-trial'', /, ''data, for which '',
	2 ''larger values are better.  It is also applicable to '',
	3 /, ''other traits, such as disease resistance, but again '',
	4 ''larger values must'', /, ''indicate better genotypes.  '',
	5 ''For instance, if the original scale is '', /,
	6 ''0 for most resistant to 5 for most diseased, then '',
	7 ''subtract those values'', /, ''from 5 in order to obtain '',
	8 ''a transformed scale with 0 worst and 5 best.'')')
      WRITE (2, '(''     In the above table, the genotypes are '',
     1 ''listed in IPC1 order,'', /, ''so those at the top and '',
     2 ''bottom have opposite GxE interaction patterns.'', /,
     3 ''Often this contrast in the genotypes has an evident '',
     4 ''agricultural'', /, ''interpretation that has a '',
     5 ''corresponding contrast in the environments'', /,
     6 ''(which are listed in their IPC1 order in the following '',
     7 ''table).'', /, ''A genotype at the top of this table has '',
     8 ''positive GxE interactions with'', /, ''environments at '',
     9 ''the top of the next table and negative GxE interactions'',
     1 /, ''with environments at the bottom of the next table; '',
     2 ''and the opposite'', /, ''patterns apply to genotypes at '',
     3 ''the bottom of this table.'')')
      WRITE (2, '(''     AMMI models to the left of the best '',
     1 ''AMMI model (for optimizing'', /, ''predictive accuracy) '',
     2 ''are too simple, so they underfit real signal;'', /,
     3 ''whereas models to the right of the best one are too '',
     4 ''complex, so they'', /, ''overfit spurious noise.  At the '',
     5 ''opposite extremes, AMMI0 always captures'', /,
     6 ''no GE-signal and no GE-noise, whereas AMMIF always '',
     7 ''captures all GE-signal'', /, ''but with this also all '',
     8 ''GE-noise.  Accordingly, a parsimonious intermediate'', /,
     9 ''AMMI model is often most predictively accurate, such as '',
     1 ''AMMI1 or AMMI2,'', /, ''because missing little GE-signal '',
     2 ''while discarding much GE-noise'', /, ''is better than '',
     3 ''either of the opposite and worse problems of missing a'', /,
     4 ''huge amount of GE-signal or capturing a huge amount of '',
     5 ''GE-noise.'', /, ''This opportunity to increase predictive '',
     6 ''accuracy is another reason why '', /, ''model diagnosis '',
     7 ''matters.'')')
      WRITE (2, '(''     Practical constraints frequently limit the '',
     1 ''number of workable m-es '', /, ''to only 2 or 3 (or perhaps '',
     2 ''a few more).  This can require a lower AMMI'', /,
     3 ''model than that diagnosed solely by statistical ''
     4 ''considerations, often'', /, ''AMMI1, in order to achieve '',
     5 ''fewer m-es.  Often this lower AMMI model is'', /,
     6 ''almost as predictively accurate as the best model, and is '',
     7 ''also far more'', /, ''accurate than AMMIF, that is, the '',
     8 ''actual data (or averages over reps).'')')
	WRITE (2, '(''     The above table displays the genotype '',
     1 ''winners and the consequent'', /, ''numbers of m-es for a '',
     2 ''wide spectrum of options from AMMI0 to AMMIF.'', /,
     3 ''Thereby it illuminates the tradeoff between statistical '',
     4 ''and practical'', /, ''considerations when choosing the '',
     5 ''most appropriate member of the AMMI'', /, ''model family '',
     6 ''to use for a given dataset.'')')
C     Print ranking table for AMMI1 and AMMIF; switch IPC1 scores.
      IF (NG .GE. 5)  THEN
      WRITE (2, '(/, /, ''AMMISOFT for '', A4,
     1 '':  Ranking Table for AMMI1 and AMMIF'')')  NMB
      WRITE (2, '(/, 76(''-''), /, 30X, ''AMMI1 Ranks'', 18X,
     1 ''AMMIF Ranks'', /, 22X, 25(''-''), 4X, 25(''-''), /,
     2 ''Environment   Ratio       1    2    3    4    5'', 4X,
     3 ''    1    2    3    4    5'', /, 76(''-''))')
      NEW = -1
      NMV = 0
      DO 130 NED=1,NE
      NER = L(KRE+NED)
      KOE = (NER - 1) * NG
      DO 131 NGD=1,NG
  131 V(NGD) = V(KBG+NGD) + (V(KPGN+NGD) * V(KPEN+NER))
      CALL INDEX(NG, 0, KR1G)
      CALL INDEX(NG, KDGE+KOE, KRFG)
      NGW1 = L(KR1G+1)
      IF (NGW1 .NE. NEW)  THEN
      NMV = NMV + 1
      L(KSWG+NMV) = NGW1
      NEW = NGW1
      IF (NED .NE. 1)  WRITE (2, '()')
      ENDIF
      IF (NGW1 .EQ. NGMW)  THEN
      RATIO = 1.0D0
      ELSE
      A1AL = V(KBG+NGW1)+V(KBE+NER)-GM+(V(KPGN+NGW1)*V(KPEN+NER))
      A1ME = V(KBG+NGMW)+V(KBE+NER)-GM+(V(KPGN+NGMW)*V(KPEN+NER))
      RATIO = A1AL / A1ME
      ENDIF
  130 WRITE (2, '(I5, 2X, A4, F8.4, 3X, 5I5, 4X, 5I5)')
     1 NER, C(KNE+NER), RATIO, (L(KR1G+LD), LD=1,5),
     2 (L(KRFG+LD), LD=1,5)
      WRITE (2, '(76(''-''))')
      DO 132 MD=1,(MEC(2) - 1)
      NGL = L(KSWG+MD)
      NGR = L(KSWG+MD+1)
      SW = (V(KBG+NGR) - V(KBG+NGL)) / (V(KPGN+NGL) - V(KPGN+NGR))
  132 WRITE (2, '(''Switch from '', A4, I6, '' to '', I5, '' at '',
     1 A4, '' IPC1 score'', F20.7)')  NMG, NGL, NGR, NME, SW
      WRITE (2, '(/, ''     In the above table, the environments are '',
     1 ''listed in IPC1 order,'', /, ''so those at the top and '',
     2 ''bottom have opposite GxE interaction patterns.'')')
      WRITE (2, '(''     Ratio is the yield (or whatever the trait) '',
     1 ''for the winner within'', /, ''each environment (identified '',
     2 ''in the first column of AMMI1 ranks),'', /, ''divided by '',
     3 ''the yield for the overall winner (which is genotype '', I5,
     4 ''),'', /, ''with both yields estimated by the AMMI1 model.  '',
     5 ''Ratio automatically'', /, ''equals 1 for the overall '',
     6 ''winner.  This ratio assesses the importance'', /,
     7 ''of narrow adaptations, which are caused by GxE '',
     8 ''interactions.  When a'', /, ''5% or 10% yield increment '',
     9 ''has agricultural or economic significance,'', /,
     1 ''a ratio of 1.05 or 1.10 or more indicates that narrow '',
     2 ''adaptations offer'', /, ''substantial opportunities for '',
     3 ''yield increases, although at the cost of'', /,
     4 ''subdividing a growing region into two or more ''
     5 ''mega-environments.'')')  NGMW
      WRITE (2, '(''     M-es have different winners and they are '',
     1 ''separated by blank lines.'', /, ''Ordinarily small m-es '',
     2 ''are merged into adjacent larger m-es, especially'', /,
     3 ''when this imposes negligible loss in yield because an '',
     4 ''adjacent winner'', /, ''holds second (or third) rank in '',
     5 ''the merged environment.  Often slight'', /, ''editing of '',
     6 ''a ranking table, obtained by deleting or moving some blank'',
     7 /, ''lines, can simplify the m-e scheme considerably.  '',
     8 ''Editing may achieve'', /, ''a favorable tradeoff:  a large '',
     9 ''gain in practicality, accompanied by a'', /, ''negligible '',
     1 ''loss of yield in several of the environments.'')')
      WRITE (2, '(''     At the bottom of the above table, the '',
     1 ''environment IPC1 scores are'', /, ''listed for each switch '',
     2 ''from one winning genotype to the next.  The'', /,
     3 ''following tables list IPC1 scores for all genotypes and '',
     4 ''environments.'')')
      WRITE (2, '(''     When crafting a m-e scheme for a given '',
     1 ''dataset, bear in mind that'', /, ''narrow adaptations '',
     2 ''caused by predictable GxE interactions increase the'', /,
     3 ''number of usable m-es, whereas unpredictable GxE '',
     4 ''interactions decrease'', /, ''this number.  Often soil '',
     5 ''and management are predictable from year to'', /,
     6 ''year, whereas weather is unpredictable.'')')
      ELSE
      WRITE (2, '(/, /, ''The ranking table is not produced '',
     1 ''because there are fewer than 5 genotypes.'')')
      ENDIF 
C     Print ranked means and ranked IPC1 scores.
      CALL INDEX(NG, KBG, KRMG)
      WRITE (2, '(/, /, ''AMMISOFT for '', A4, '':  Ranked Means '',
     1 ''and IPC1 Scores'')')  NMB
      WRITE (2, '(/, 66(''-''), /, 1X, A4, ''  Code          Mean'',
     1 8X, ''|  '', A4, ''  Code      IPC1 Score'', /,
     2 66(''-''))')  NMG, NMG
      DO 140 NGD=1,NG
      IRM = L(KRMG+NGD)
      IRP = L(KRG+NGD)
  140 WRITE (2, '(I5, 2X, A4, F20.5, ''  | '', I5, 2X, A4, F20.7)')
     1 IRM, C(KNG+IRM), V(KBG+IRM), IRP, C(KNG+IRP), V(KPGN+IRP)
      WRITE (2, '(66(''-''))')
      WRITE (2, '(''Grand mean'', F21.5)')  GM
      AA = 0.0D0
      BB = 0.0D0
      AB = 0.0D0
      DO 141 NGD=1,NG
      A = V(KBG+NGD) - GM
      B = V(KPGN+NGD)
      AA = AA + (A * A)
      BB = BB + (B * B)
  141 AB = AB + (A * B)
      R = AA * BB
      IF (R .GT. 0.0D0)  THEN
      CORR = AB / DSQRT(R)
      WRITE (2, '(''Correlation between '', A4, '' means and IPC1 '',
     1 ''scores'', F11.7)')  NMG, CORR
      ENDIF
      CALL INDEX(NE, KBE, KRME)
      WRITE (2, '(/, 66(''-''), /, 1X, A4, ''  Code          Mean'',
     1 8X, ''|  '', A4, ''  Code      IPC1 Score'', /,
     2 66(''-''))')  NME, NME
      DO 142 NED=1,NE
      IRM = L(KRME+NED)
      IRP = L(KRE+NED)
  142 WRITE (2, '(I5, 2X, A4, F20.5, ''  | '', I5, 2X, A4, F20.7)')
     1 IRM, C(KNE+IRM), V(KBE+IRM), IRP, C(KNE+IRP), V(KPEN+IRP)
      WRITE (2, '(66(''-''))')
      WRITE (2, '(''Grand mean'', F21.5)')  GM
      AA = 0.0D0
      BB = 0.0D0
      AB = 0.0D0
      DO 143 NED=1,NE
      A = V(KBE+NED) - GM
      B = V(KPEN+NED)
      AA = AA + (A * A)
      BB = BB + (B * B)
  143 AB = AB + (A * B)
      R = AA * BB
      IF (R .GT. 0.0D0)  THEN
      CORR = AB / DSQRT(R)
      WRITE (2, '(''Correlation between '', A4, '' means and IPC1 '',
     1 ''scores'', F11.7)')  NME, CORR
      ENDIF
      WRITE (2, '(/, ''     The above ranked lists should be '',
     1 ''inspected carefully to determine'', /, ''whether they have '',
     2 ''a plausible biological or ecological interpretation.'', /,
     3 ''If genotype means and IPC1 scores have a small '',
     4 ''correlation, then these'', /, ''means and scores call for '',
     5 ''different biological explanations; otherwise,'', /,
     6 ''if a large correlation, then the same explanation.  The '',
     7 ''same principles'', /, ''apply to ecological interpretation '',
     8 ''of environment means and IPC1 scores.'')')
      WRITE (2, '(''     Genotype IPC1 scores and environment IPC1 '',
     1 ''scores require a coherent'', /, ''interpretation, such as '',
     2 ''drought tolerance for genotypes and rainfall for'', /,
     3 ''environments.  If either genotypes or environments are '',
     4 ''better known than'', /, ''the other, begin with whichever '',
     5 ''is more familiar.  Start by contrasting'', /, ''the top '',
     6 ''several genotypes (or environments) with the bottom '',
     7 ''several ones'', /, ''to suggest a biological (or '',
     8 ''ecological) interpretation, and then inspect'', /,
     9 ''the entire list to confirm a systematic trend and clear '',
     1 ''interpretation.'')')
      IF (NR .GT. 1)  THEN
      IF (MFR .EQ. 0)  WRITE (2, '(''     However, given the '',
     1 ''diagnosis of AMMI0 by FR-tests for this dataset,'', /,
     2 ''a compelling biological or environmental interpretation '',
     3 ''of these IPC1'', /, ''scores is rather unlikely.'')')
      ENDIF
C     If missing cells, print list of imputed values.
      IF (NCA .LT. NGE)  THEN
      WRITE (2, '(/, /, ''AMMISOFT for '', A4, '':  List of '',
     1 ''EM-AMMI1 Imputed Values'')')  NMB
      WRITE (2, '(/, 34(''-''), /, 1X, A4, 4X, A4, 8X, ''Imputed '',
     1 ''Value'', /, 34(''-''))')  NME, NMG
      KSUB = 0
      DO 150 NED=1,NE
      DO 151 NGD=1,NG
      KSUB = KSUB + 1
      IF (M(KSGE+KSUB) .EQ. 0)  WRITE (2, '(I5, 3X, I5, F21.5)')
     1 NED, NGD, V(KDGE+KSUB)
  151 CONTINUE
  150 CONTINUE
      WRITE (2, '(34(''-''))')
      WRITE (2, '(''Grand mean'', F24.5, /)')  GM
      VMINA = 1000000.0D0 * GM
      VMAXA = - GM
      VMINI = 1000000.0D0 * GM
      VMAXI = - GM
      DO 152 ND=1,NGE
      R = V(KDGE+ND)
      IF (M(KSGE+ND) .EQ. 1)  THEN
      IF (R .LT. VMINA)  VMINA = R
      IF (R .GT. VMAXA)  VMAXA = R
      ELSE
      IF (R .LT. VMINI)  VMINI = R
      IF (R .GT. VMAXI)  VMAXI = R
      ENDIF
  152 CONTINUE
      WRITE (2, '(''Actual data   '', F20.5, '' to'', F20.5)')
     1 VMINA, VMAXA
      WRITE (2, '(''Imputed values'', F20.5, '' to'', F20.5)')
     1 VMINI, VMAXI
      EC = 100.0D0 * (DFLOAT(NGE - NCA) / FGE)
      WRITE (2, '(/, ''     This dataset has'', F7.3,
     1 ''% empty cells with no data.  Inspect these'', /, ''imputed '',
     2 ''values to check for any that are especially implausible.'')')
     3 EC
      ENDIF
C     Print computer-readable AMMI parameters.
      WRITE (2, '(/, /, ''AMMISOFT for '', A4, '':  Computer-'',
     1 ''readable AMMI Parameters'')')  NMB
      WRITE (2, '(4I5, F20.5)')  NG, NE, NR, NX, GM
      DO 160 NGD=1,NG
  160 WRITE (2, '(I5, 2X, A4, F20.5, 7F20.7)')  NGD, C(KNG+NGD),
	1 V(KBG+NGD), (V(KPGN+((NXD-1)*NG)+NGD), NXD=1,NX)
      DO 161 NED=1,NE
  161 WRITE (2, '(I5, 2X, A4, F20.5, 7F20.7)')  NED, C(KNE+NED),
     1 V(KBE+NED), (V(KPEN+((NXD-1)*NE)+NED), NXD=1,NX)
      RETURN
      END
      SUBROUTINE FTEST
C     Given a source df and SS of NDS and SSS, and an error df and SS
C     of NDE and SSE, this subroutine computes the source MS of SMS,
C     F-test probability of PF, and significance annotation of SI.
C     This subroutine is adapted from William J. Kennedy, Jr. and
C     James E. Gentle, 1980, Statistical Computing, Marcel Dekker,
C     New York, pages 112-116.
C     Modified December 1990:  Limit NDS and NDE to 1000 in order
C     to avoid math overflow.
      IMPLICIT  REAL*8  (A-H, O-Z)
      CHARACTER*4  SI
      COMMON  /FTST/  SSS, SSE, SMS, PF, NDS, NDE, SI
      SMS = SSS / DFLOAT(NDS)
      EMS = SSE / DFLOAT(NDE)
      Z = SMS / EMS
      N1 = NDS
      IF (N1 .GT. 1000)  N1 = 1000
      AN1 = N1
      N2 = NDE
      IF (N2 .GT. 1000)  N2 = 1000
      AN2 = N2
      A = AN1 * Z / (AN1 * Z + AN2)
      A1 = 1.0D0 - A
      IF (A1 .LT. 0.1D-36)  A1 = 0.1D-36
      D1 = AN1 * 0.5D0
      D2 = AN2 * 0.5D0
      D3 = D1 + D2 - 1.0D0
      R = 0.0D0
      S1 = 0.0D0
      S2 = 0.0D0
      DL = 1.0D0
      XM = 1.0D0
      XK = 1.0D0
      CC = 0.25D0
      PI = 3.141592653589793D0
      N = N2
   10 MM = D2
      MM = 2 * MM
      IF (MM .NE. N)  GO TO 12
      N = D2 - 1
      IF (N .EQ. 0)  GO TO 11
      DO 100 I=1,N
      S1 = DL + S1 * R
      D2 = D2 - 1.0D0
      D3 = D3 - 1.0D0
      TEM = A1 / D2
      R = D3 * TEM
      S2 = (R + TEM) * S2
  100 CONTINUE
   11 CONTINUE
      S1 = DL + S1 * R
      DL = 0.0D0
      T = -1.0D0
      D3 = -1.0D0
      S2 = A * S2
      CC = CC + 0.5D0
      GO TO 14
   12 CONTINUE
      N = D2
      IF (N .EQ. 0)  GO TO 13
      DO 101 I=1,N
      S1 = DL + S1 * R
      D2 = D2 - 1.0D0
      D3 = D3 - 1.0D0
      TEM = A1 / D2
      R = D3 * TEM
      S2 = (R + TEM) * S2
  101 CONTINUE
   13 S1 = XK * S1
      S2 = XK * S2
      ART = DSQRT(A1)
      XM = XM * ART
      T = (XM - ART) / A1
      D3 = -0.5D0
      XK = 2.0D0 / PI
      CC = CC * 2.0D0
   14 IF (CC .GT. 0.875D0)  GO TO 15
      D2 = D1
      D3 = D2 + D3
      S2 = S1
      S1 = 0.0D0
      A1 = A
      IF (A1 .LT. 0.1D-36)  A1 = 0.1D-36
      N = N1
      GO TO 10
   15 IF (CC .LT. 1.125D0)  DL = 4.0D0 / PI * DATAN(T)
      PDFF = XM * (S2 - S1) - DL
      IF (PDFF .LT. 0.0D0)  PDFF = 0.0D0
      IF (PDFF .GT. 1.0D0)  PDFF = 1.0D0
      PF = 1.0D0 - PDFF
      SI = ' ***'
      IF (PF .LE. 0.001D0)  RETURN
      SI = ' ** '
      IF (PF .LE. 0.01D0)  RETURN
      SI = ' *  '
      IF (PF .LE. 0.05D0)  RETURN
      SI = '    '
      RETURN
      END
      SUBROUTINE INDEX(N, KBV, KBL)
C     Index N elements of vector V with base address KBV, placing
C     the indices in decending order in vector L with base address
C     KBL, using the Heapsort sorting algorithm with indirect
C     indexing through L in all references to V.
C     Hence the first index, L(KBL+1), points to V(KBV+L(KBL+1))
C     containing the largest value in vector V.
C     The contents of V(KBV) are unaltered.
C     This subroutine is adapted from William H. Press, Brian P.
C     Flannery, Saul A. Teukolsky, and William T. Vetterling,
C     1986, Numerical Recipes, Cambridge University Press,
C     Cambridge, England, pages 232-233.
      IMPLICIT  REAL*8  (A-H, O-X)
      REAL*8     V(32000000)
      INTEGER*4  L(64000000)
      COMMON  /DATA/  Y(64000000)
      EQUIVALENCE  (V(1), Y(1), L(1))
      DO 100 J=1,N
      L(KBL+J) = J
  100 CONTINUE
      K = (N / 2) + 1
      IR = N
   10 CONTINUE
      IF (K .GT. 1)  THEN
      K = K - 1
      LT = L(KBL+K)
      X = V(KBV+LT)
      ELSE
      LT = L(KBL+IR)
      X = V(KBV+LT)
      L(KBL+IR) = L(KBL+1)
      IR = IR - 1
      IF (IR .EQ. 1)  THEN
      L(KBL+1) = LT
      GO TO 30
      ENDIF
      ENDIF
      I = K
      J = K + K
   20 IF (J .LE. IR)  THEN
      IF (J .LT. IR)  THEN
      KS = KBL + J
      IF (V(KBV+L(KS)) .LT. V(KBV+L(KS+1)))  J = J + 1
      ENDIF
      IF (X .LT. V(KBV+L(KBL+J)))  THEN
      L(KBL+I) = L(KBL+J)
      I = J
      J = J + J
      ELSE
      J = IR + 1
      ENDIF
      GO TO 20
      ENDIF
      L(KBL+I) = LT
      GO TO 10
   30 KS = KBL + N + 1
      DO 101 J=1,N/2
      I = KBL + J
      K = KS - J
      LT = L(I)
      L(I) = L(K)
  101 L(K) = LT
      RETURN
      END

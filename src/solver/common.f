* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C         control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C           be checked for duplicates, and out-of-range indices.
C           Note  simple checks are always performed.
C           ICNTL(1) = 0, data checking performed.
C           Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C           ICNTL(2) = 0, input is by rows and columns in arbitrary
C           order and the output is sorted by columns.
C           ICNTL(2) = 1, the output is also row ordered
C           within each column.
C           ICNTL(2) = 2, the input is already ordered by
C           columns and is to be row ordered within each column.
C           Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C           ICNTL(3) = 0, matrix entries are ordered.
C           Otherwise, only the sparsity pattern is ordered
C           and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C           which error messages are sent. Error messages
C           can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C           which warning messages are sent. Warning
C           messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C           must be set to 0.
C           If ICNTL(6) = -1 or 1, symmetric and only the lower
C           triangular part of the reordered matrix is returned.
C           If ICNTL(6) = -2 or 2, Hermitian and only the lower
C           triangular part of the reordered matrix is returned.
C           If error checks are performed (ICNTL(1) = 0)
C           and ICNTL(6)> 1 or 2, the values of duplicate
C           entries are added together; if ICNTL(6) < -1 or -2, the
C           value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of columns in the matrix.
C NR:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of rows in the matrix.
C NE:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C            user to hold the row indices of the entries in the matrix.
C          If ICNTL(2).NE.2, the entries may be in any order.
C          If ICNTL(2).EQ.2, the entries in column J must be in
C            positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C            indices are reordered so that the entries of a single
C            column are contiguous with column J preceding column J+1, J
C            = 1, 2, ..., NC-1, with no space between columns.
C          If ICNTL(2).EQ.0, the order within each column is arbitrary;
C            if ICNTL(2) = 1 or 2, the order within each column is by
C            increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:     INTEGER array of length LJCN. Intent (INOUT).
C          If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C          to the column index of the entry
C          whose row index is held in IRN(K), K = 1, 2, ..., NE.
C          On exit, the contents of this array  will have been altered.
C          If ICNTL(2) = 2, the array is not accessed.
C LA:      INTEGER variable. Intent(IN). Defines length of array
C          A.
C A:       is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C          the I version, COMPLEX in the C version,
C          or COMPLEX"*"16 in the Z version) array of length LA.
C          Intent(INOUT).
C          If ICNTL(3).EQ.0, A(K) must be set by the user to
C          hold the value of the entry with row index IRN(K),
C          K = 1, 2, ..., NE. On exit, the array will have been
C          permuted in the same way as the array IRN.
C          If ICNTL(3).NE.0, the array is not accessed.
C LIP:     INTEGER variable. Intent(IN). Defines length of array
C          IP.
C IP:      INTEGER array of length LIP. Intent(INOUT). IP
C          need only be set by the user if ICNTL(2) = 2.
C          In this case, IP(J) holds the position in
C          the array IRN of the first entry in column J, J = 1, 2,
C          ..., NC, and IP(NC+1) is one greater than the number of
C          entries in the matrix.
C          In all cases, the array IP will have this meaning on exit
C          from the subroutine and is altered when ICNTL(2) = 2 only
C          when ICNTL(1) =  0 and there are out-of-range
C          indices or duplicates.
C LIW:     INTEGER variable. Intent(IN). Defines length of array
C          IW.
C IW:      INTEGER array of length LIW. Intent(OUT). Used by the
C          routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C          a negative value of INFO(1) is used to signal a fatal
C          error in the input data, a positive value of INFO(1)
C          indicates that a warning has been issued, and a
C          zero value is used to indicate a successful call.
C          In cases of error, further information is held in INFO(2).
C          For warnings, further information is
C          provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C          currently used and are set to zero.
C          Possible nonzero values of INFO(1):
C         -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C               Value of ICNTL(2) is given by INFO(2).
C         -2 -  NC.LE.0. Value of NC is given by INFO(2).
C         -3 -  Error in NR. Value of NR is given by INFO(2).
C         -4 -  NE.LE.0. Value of NE is given by INFO(2).
C         -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C         -6 -  LA too small. Min. value of LA is given by INFO(2).
C         -7 -  LIW too small. Value of LIW is given by INFO(2).
C         -8 -  LIP too small. Value of LIP is given by INFO(2).
C         -9 -  The entries of IP not monotonic increasing.
C        -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C        -11 -  ICNTL(6) is out of range.
C         +1 -  One or more duplicated entries. One copy of
C               each such entry is kept and, if ICNTL(3) = 0 and
C               ICNTL(6).GE.0, the values of these entries are
C               added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C               the value of the first occurrence of the entry is used.
C               Initially INFO(3) is set to zero. If an entry appears
C               k times, INFO(3) is incremented by k-1 and INFO(6)
C               is set to the revised number of entries in the
C               matrix.
C         +2 - One or more of the entries in IRN out-of-range. These
C               entries are removed by the routine.`INFO(4) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix.
C         +4 - One or more of the entries in JCN out-of-range. These
C               entries are removed by the routine. INFO(5) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix. Positive values of INFO(1) are summed so that
C               the user can identify all warnings.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. External Subroutines ..
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..

C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
      LP = ICNTL(4)
      MP = ICNTL(5)

C  Check the input data
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF

      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF

      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF

      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF

      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF

C Check workspace sufficient
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF

      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

C PART is used by MC59BD to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
      PART = 0
      IF (ICNTL6.NE.0) PART = 1

      IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59BD, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C Check for duplicates
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)

      ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
          IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF

C  Order by rows within each column
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)

      END IF

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C          for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C             (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C             (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        in arbitrary order.
C      - on exit, the entries in IRN are reordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C      - JCN(K) must be the column index of
C        the entry in IRN(K)
C      - on exit, JCN(K) is the column index for the entry with
C        row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in (IRN(K), JCN(K));
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C      - the array is used as workspace
C      - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C        of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C        of entries in matrix after removal of out-of-range entries.
C        If no data checking, KNE = NE.

C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C     ..
C     .. Executable Statements ..

C Initialise IW
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE

      KNE = 0
      IOUT = 0
      JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
      IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
C Unsymmetric
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
C Pattern only
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 130

      ELSE

C No checks
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Lower triangle ... swap if necessary
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Upper triangle ... swap if necessary
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

      IF (LA.EQ.1) THEN
C Pattern only
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE

        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF

  130 CONTINUE

      RETURN
      END
C
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C      - not set on entry.
C      - on exit,  IRN holds row indices with the row
C        indices for column 1 preceding those for column 2 and so on,
C        with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C        with indices for column 1 preceding those for column 2
C        and so on, with the order within columns is arbitrary.
C      - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in JCN(K);
C        on exit A, A(K) holds the value of the entry in IRN(K).
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C      - on entry, must be set on entry so that IW(J) points to the
C        position in JCN of the first entry in row J, J=1,...,NR, and
C        IW(NR+1) must be set to NE+1
C
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C     ..
C     .. Executable Statements ..

C  Count the number of entries in each column

      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE

      IF (LA.GT.1) THEN

        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE

C Pattern only

C  Count the number of entries in each column

        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE

      END IF

      RETURN
      END

C**********************************************************

      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        ordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C        On exit, the order within each column is by increasing
C        row indices.
C   LA - integer variable which defines the length of the array A.
C        Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in IRN(K);
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C      - on entry, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C     . .
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..

C Jump if pattern only.
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
C Next column
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE

C Pattern only.
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
C Next column
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************

      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      IF (LA.GT.1) THEN
C Matrix entries considered
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END
C***********************************************************************

      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59ED except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END

* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21BD
C     ..
C     .. Executable Statements ..
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NZ)
      INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AVAL
      INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
      IW(1,1) = 1
      IW(1,2) = LENROW(1)
      DO 10 I = 2,N
        IW(I,1) = IW(I-1,1) + LENROW(I-1)
        IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C     OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C     POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
      JJ = 1
      DO 20 I = 1,N
        IOLD = IP(I)
        IOLD = IABS(IOLD)
        LENGTH = IW(IOLD,2)
        LENROW(I) = LENGTH
        IF (LENGTH.EQ.0) GO TO 20
        IW(IOLD,1) = IW(IOLD,1) - JJ
        J2 = JJ + LENGTH - 1
        DO 15 J = JJ,J2
          IW1(J) = IOLD
   15   CONTINUE
        JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
      DO 30 I = 1,N
        IOLD = IQ(I)
        IOLD = IABS(IOLD)
        IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C     IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C     SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C     ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
      DO 200 I = 1,NZ
        IOLD = IW1(I)
        IF (IOLD.EQ.0) GO TO 200
        IPOS = I
        JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
        IF (IW(IOLD,1).EQ.0) GO TO 150
        AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C     IN ITS FINAL POSITION  .. VIZ. IPOS.
        DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C     THE CHAIN.
          NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
          IF (NEWPOS.EQ.I) GO TO 130
          A(IPOS) = A(NEWPOS)
          JNUM = ICN(NEWPOS)
          ICN(IPOS) = IW(JNUM,2)
          IPOS = NEWPOS
          IOLD = IW1(IPOS)
          IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

      END
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

      SUBROUTINE MC59A(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C         control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C           be checked for duplicates, and out-of-range indices.
C           Note  simple checks are always performed.
C           ICNTL(1) = 0, data checking performed.
C           Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C           ICNTL(2) = 0, input is by rows and columns in arbitrary
C           order and the output is sorted by columns.
C           ICNTL(2) = 1, the output is also row ordered
C           within each column.
C           ICNTL(2) = 2, the input is already ordered by
C           columns and is to be row ordered within each column.
C           Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C           ICNTL(3) = 0, matrix entries are ordered.
C           Otherwise, only the sparsity pattern is ordered
C           and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C           which error messages are sent. Error messages
C           can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C           which warning messages are sent. Warning
C           messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C           must be set to 0.
C           If ICNTL(6) = -1 or 1, symmetric and only the lower
C           triangular part of the reordered matrix is returned.
C           If ICNTL(6) = -2 or 2, Hermitian and only the lower
C           triangular part of the reordered matrix is returned.
C           If error checks are performed (ICNTL(1) = 0)
C           and ICNTL(6)> 1 or 2, the values of duplicate
C           entries are added together; if ICNTL(6) < -1 or -2, the
C           value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of columns in the matrix.
C NR:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of rows in the matrix.
C NE:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C            user to hold the row indices of the entries in the matrix.
C          If ICNTL(2).NE.2, the entries may be in any order.
C          If ICNTL(2).EQ.2, the entries in column J must be in
C            positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C            indices are reordered so that the entries of a single
C            column are contiguous with column J preceding column J+1, J
C            = 1, 2, ..., NC-1, with no space between columns.
C          If ICNTL(2).EQ.0, the order within each column is arbitrary;
C            if ICNTL(2) = 1 or 2, the order within each column is by
C            increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:     INTEGER array of length LJCN. Intent (INOUT).
C          If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C          to the column index of the entry
C          whose row index is held in IRN(K), K = 1, 2, ..., NE.
C          On exit, the contents of this array  will have been altered.
C          If ICNTL(2) = 2, the array is not accessed.
C LA:      INTEGER variable. Intent(IN). Defines length of array
C          A.
C A:       is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C          the I version, COMPLEX in the C version,
C          or COMPLEX"*"16 in the Z version) array of length LA.
C          Intent(INOUT).
C          If ICNTL(3).EQ.0, A(K) must be set by the user to
C          hold the value of the entry with row index IRN(K),
C          K = 1, 2, ..., NE. On exit, the array will have been
C          permuted in the same way as the array IRN.
C          If ICNTL(3).NE.0, the array is not accessed.
C LIP:     INTEGER variable. Intent(IN). Defines length of array
C          IP.
C IP:      INTEGER array of length LIP. Intent(INOUT). IP
C          need only be set by the user if ICNTL(2) = 2.
C          In this case, IP(J) holds the position in
C          the array IRN of the first entry in column J, J = 1, 2,
C          ..., NC, and IP(NC+1) is one greater than the number of
C          entries in the matrix.
C          In all cases, the array IP will have this meaning on exit
C          from the subroutine and is altered when ICNTL(2) = 2 only
C          when ICNTL(1) =  0 and there are out-of-range
C          indices or duplicates.
C LIW:     INTEGER variable. Intent(IN). Defines length of array
C          IW.
C IW:      INTEGER array of length LIW. Intent(OUT). Used by the
C          routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C          a negative value of INFO(1) is used to signal a fatal
C          error in the input data, a positive value of INFO(1)
C          indicates that a warning has been issued, and a
C          zero value is used to indicate a successful call.
C          In cases of error, further information is held in INFO(2).
C          For warnings, further information is
C          provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C          currently used and are set to zero.
C          Possible nonzero values of INFO(1):
C         -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C               Value of ICNTL(2) is given by INFO(2).
C         -2 -  NC.LE.0. Value of NC is given by INFO(2).
C         -3 -  Error in NR. Value of NR is given by INFO(2).
C         -4 -  NE.LE.0. Value of NE is given by INFO(2).
C         -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C         -6 -  LA too small. Min. value of LA is given by INFO(2).
C         -7 -  LIW too small. Value of LIW is given by INFO(2).
C         -8 -  LIP too small. Value of LIP is given by INFO(2).
C         -9 -  The entries of IP not monotonic increasing.
C        -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C        -11 -  ICNTL(6) is out of range.
C         +1 -  One or more duplicated entries. One copy of
C               each such entry is kept and, if ICNTL(3) = 0 and
C               ICNTL(6).GE.0, the values of these entries are
C               added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C               the value of the first occurrence of the entry is used.
C               Initially INFO(3) is set to zero. If an entry appears
C               k times, INFO(3) is incremented by k-1 and INFO(6)
C               is set to the revised number of entries in the
C               matrix.
C         +2 - One or more of the entries in IRN out-of-range. These
C               entries are removed by the routine.`INFO(4) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix.
C         +4 - One or more of the entries in JCN out-of-range. These
C               entries are removed by the routine. INFO(5) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix. Positive values of INFO(1) are summed so that
C               the user can identify all warnings.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. External Subroutines ..
      EXTERNAL MC59B,MC59C,MC59D,MC59E,MC59F
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..

C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
      LP = ICNTL(4)
      MP = ICNTL(5)

C  Check the input data
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF

      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF

      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF

      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF

      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF

C Check workspace sufficient
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF

      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

C PART is used by MC59B to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
      PART = 0
      IF (ICNTL6.NE.0) PART = 1

      IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59B, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
        CALL MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C Check for duplicates
        IF (LCHECK) CALL MC59E(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)

      ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59B(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
        IF (LCHECK) CALL MC59E(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
        CALL MC59C(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
        IF (LCHECK) THEN
          CALL MC59F(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
          IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF

C  Order by rows within each column
        CALL MC59D(NC,KNE,IRN,IP,LAA,A)

      END IF

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59A *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59A *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C          for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C             (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C             (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        in arbitrary order.
C      - on exit, the entries in IRN are reordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C      - JCN(K) must be the column index of
C        the entry in IRN(K)
C      - on exit, JCN(K) is the column index for the entry with
C        row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in (IRN(K), JCN(K));
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C      - the array is used as workspace
C      - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C        of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C        of entries in matrix after removal of out-of-range entries.
C        If no data checking, KNE = NE.

C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      REAL ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C     ..
C     .. Executable Statements ..

C Initialise IW
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE

      KNE = 0
      IOUT = 0
      JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
      IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
C Unsymmetric
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
C Pattern only
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 130

      ELSE

C No checks
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Lower triangle ... swap if necessary
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Upper triangle ... swap if necessary
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

      IF (LA.EQ.1) THEN
C Pattern only
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE

        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF

  130 CONTINUE

      RETURN
      END
C
C**********************************************************
      SUBROUTINE MC59C(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C      - not set on entry.
C      - on exit,  IRN holds row indices with the row
C        indices for column 1 preceding those for column 2 and so on,
C        with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C        with indices for column 1 preceding those for column 2
C        and so on, with the order within columns is arbitrary.
C      - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in JCN(K);
C        on exit A, A(K) holds the value of the entry in IRN(K).
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C      - on entry, must be set on entry so that IW(J) points to the
C        position in JCN of the first entry in row J, J=1,...,NR, and
C        IW(NR+1) must be set to NE+1
C
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      REAL ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C     ..
C     .. Executable Statements ..

C  Count the number of entries in each column

      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE

      IF (LA.GT.1) THEN

        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE

C Pattern only

C  Count the number of entries in each column

        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE

      END IF

      RETURN
      END

C**********************************************************

      SUBROUTINE MC59D(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        ordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C        On exit, the order within each column is by increasing
C        row indices.
C   LA - integer variable which defines the length of the array A.
C        Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in IRN(K);
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C      - on entry, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C     . .
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER IRN(NE),IP(NC)
C     ..
C     .. Local Scalars ..
      REAL ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..

C Jump if pattern only.
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
C Next column
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE

C Pattern only.
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
C Next column
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************

      SUBROUTINE MC59E(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      IF (LA.GT.1) THEN
C Matrix entries considered
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END
C***********************************************************************

      SUBROUTINE MC59F(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59E except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C     ..
C     .. Array Arguments ..
      REAL A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END

* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21B
C     ..
C     .. Executable Statements ..
      CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
     +           IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC22A(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      REAL A(NZ)
      INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C     ..
C     .. Local Scalars ..
      REAL AVAL
      INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
      IW(1,1) = 1
      IW(1,2) = LENROW(1)
      DO 10 I = 2,N
        IW(I,1) = IW(I-1,1) + LENROW(I-1)
        IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C     OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C     POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
      JJ = 1
      DO 20 I = 1,N
        IOLD = IP(I)
        IOLD = IABS(IOLD)
        LENGTH = IW(IOLD,2)
        LENROW(I) = LENGTH
        IF (LENGTH.EQ.0) GO TO 20
        IW(IOLD,1) = IW(IOLD,1) - JJ
        J2 = JJ + LENGTH - 1
        DO 15 J = JJ,J2
          IW1(J) = IOLD
   15   CONTINUE
        JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
      DO 30 I = 1,N
        IOLD = IQ(I)
        IOLD = IABS(IOLD)
        IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C     IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C     SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C     ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
      DO 200 I = 1,NZ
        IOLD = IW1(I)
        IF (IOLD.EQ.0) GO TO 200
        IPOS = I
        JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
        IF (IW(IOLD,1).EQ.0) GO TO 150
        AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C     IN ITS FINAL POSITION  .. VIZ. IPOS.
        DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C     THE CHAIN.
          NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
          IF (NEWPOS.EQ.I) GO TO 130
          A(IPOS) = A(NEWPOS)
          JNUM = ICN(NEWPOS)
          ICN(IPOS) = IW(JNUM,2)
          IPOS = NEWPOS
          IOLD = IW1(IPOS)
          IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

      END
! COPYRIGHT (c) 1999 Council for the Central Laboratory
!                    of the Research Councils
! Original date July 1999
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!         Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.6.0
! See ChangeLog for version history.
!


      SUBROUTINE MC64I(ICNTL)
      IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER ICNTL(10)
C
C  Local variables
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the defaut value, the incoming data is checked for
C     out-of-range indices and duplicates.  Setting ICNTL(4) to any
C     other will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) non-zero, if the data is
C     known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C     zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64A(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      REAL A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible.
C     The values on the diagonal are of arbitrary size. HSL subroutine
C     MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal
C     entries of the permuted matrix are one in absolute value and
C     all the off-diagonal entries are less than or equal to one in
C     absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a REAL array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >= 5N
C     JOB = 2 :  LIW >= 4N
C     JOB = 3 :  LIW >= 10N + NE
C     JOB = 4 :  LIW >= 5N
C     JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW is not used
C     JOB = 2 :  LDW >= N
C     JOB = 3 :  LDW >= NE
C     JOB = 4 :  LDW >= 2N + NE
C     JOB = 5 :  LDW >= 3N + NE
C
C DW is a REAL array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix.
C        (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C        column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C       "Algorithm 575. Permutations for a zero-free diagonal",
C       ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C       "The design and use of algorithms for permuting large
C       entries to the diagonal of sparse matrices",
C       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C       "On algorithms for permuting large entries to the diagonal
C       of sparse matrices",
C       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
      INTEGER I,J,K
      REAL FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
C External routines and functions
      EXTERNAL MC21A,MC64B,MC64R,MC64S,MC64W
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
      RINF = HUGE(RINF)

C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
C If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
        CALL MC21A(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
        CALL MC64B(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
        CALL MC64R(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
        CALL MC64S(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64W(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64W(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
C       GO TO 90
      ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF

C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PE14.4)/(14X,4(1PE14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC64B(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      REAL A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a REAL array of length N.

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      REAL CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
      REAL RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
      EXTERNAL MC64D,MC64E,MC64F


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64D(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64E(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64F(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64D(I,N,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
      IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64D(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DI

      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
C POS is index of current position of I in the tree
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64E(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      REAL DK,DR,DI

C Move last element to begin of Q
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64F(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64W (IWAY = 2).
C Move last element in the heap

      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DK,DR,DI

C Quick return, if possible
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64R(N,NE,IP,IRN,A)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      REAL A(NE)

C This subroutine sorts the entries in each column of the
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
C Local variables
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      REAL HA,KEY
C Local arrays
      INTEGER TODO(TDLEN)

      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)

C Sort array roughly with partial quicksort
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
        TD = TD - 2
        GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
C          TODO(TD-1) = FIRST
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2

  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
C There is still work to be done
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
        TD = TD - 2
        GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE

  100 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64S(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      REAL A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a REAL array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(k). The entries in each column must be
C   non-negative and ordered by decreasing value.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of
C   unmatched columns.
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know
C   a matching does not exist; in the matrix defined by the column
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      REAL BVAL,BMIN,BMAX,RINF
      EXTERNAL MC64Q,MC64U

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped.
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Compute a first maximum matching from scratch on whole matrix.
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64U(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
C Matrix is structurally singular
        BMAX = RINF
      ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute
C values of the columns
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF

C Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J.
C Set LENL(J) just after last entry in column J with value ge BMAX.
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE

C Main loop
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
            CALL MC64Q(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
C Set LEN such that all matrix entries with value lt BVAL are
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
              LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81     MOD = 1
        ELSE
C We do not have a maximum matching in IW.
          BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C          IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
          CALL MC64Q(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64U(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) =  - J
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64Q(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      REAL A(*),VAL

C This routine searches for at most XX different numerical values
C in the columns W(1:WLEN). XX>=2.
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2).
C
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      REAL SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT
C remains sorted by decreasing value.
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
C Check presence of HA in SPLIT
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
C The insertion
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
C Exit loop if XX values are found
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64U(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1

      IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C      IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
      NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1)*N

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

      DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
        ID1 = ID0 + JORD
C J is unmatched column
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row
          ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30     OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
          LAST = JORD
          GO TO 101
        ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
      LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64W(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      REAL A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   weighted matching: IPERM(I) = 0 or row I is matched to column
C   IPERM(I).
C NUM is an INTEGER variable. On exit, it contains the cardinality of
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a REAL array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and
C   V = D(N+1:2N) contains the dual column variable. If the matrix
C   is structurally nonsingular (NUM = N), the following holds:
C      U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C      U(I) = 0  if IPERM(I) = 0
C      V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      REAL CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
      REAL RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
      EXTERNAL MC64D,MC64E,MC64F


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
C Initialize U(I)
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
C Row I is not empty
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
C Don't choose cheap assignment from dense columns
        IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
C Assignment of column J to row I
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
      DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
C Continue only if column J is not empty
        IF (K1.GT.K2) GO TO 95
C       VJ = RINF
C Changes made to allow for NaNs
        I0 = IRN(K1)
        VJ = A(K1) - U(I0)
        K0 = K1
        DO 50 K = K1+1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
C Scan remaining part of assigned column JJ
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
        CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64D(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE

        DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64E(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1

C Scan column that matches with row Q0
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
            DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
C Row I is matched; do not update D(I) if DNEW is larger
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64F(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64D(I,N,Q,D,L,2)
              ENDIF
C Update tree
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

      IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END


! COPYRIGHT (c) 1999 Council for the Central Laboratory
!                    of the Research Councils
! Original date July 1999
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!         Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.6.0
! See ChangeLog for version history.
!
      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER ICNTL(10)
C
C  Local variables
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the defaut value, the incoming data is checked for
C     out-of-range indices and duplicates.  Setting ICNTL(4) to any
C     other will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) non-zero, if the data is
C     known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C     zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible.
C     The values on the diagonal are of arbitrary size. HSL subroutine
C     MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal
C     entries of the permuted matrix are one in absolute value and
C     all the off-diagonal entries are less than or equal to one in
C     absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a DOUBLE PRECISION array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >= 5N
C     JOB = 2 :  LIW >= 4N
C     JOB = 3 :  LIW >= 10N + NE
C     JOB = 4 :  LIW >= 5N
C     JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW is not used
C     JOB = 2 :  LDW >= N
C     JOB = 3 :  LDW >= NE
C     JOB = 4 :  LDW >= 2N + NE
C     JOB = 5 :  LDW >= 3N + NE
C
C DW is a DOUBLE PRECISION array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix.
C        (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C        column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C       "Algorithm 575. Permutations for a zero-free diagonal",
C       ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C       "The design and use of algorithms for permuting large
C       entries to the diagonal of sparse matrices",
C       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C       "On algorithms for permuting large entries to the diagonal
C       of sparse matrices",
C       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
C External routines and functions
      EXTERNAL MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
      RINF = HUGE(RINF)

C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
C If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
        CALL MC64RD(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
C       GO TO 90
      ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF

C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a DOUBLE PRECISION array of length N.

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
      IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI

      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
C POS is index of current position of I in the tree
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Move last element to begin of Q
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64WD (IWAY = 2).
C Move last element in the heap

      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Quick return, if possible
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)

C This subroutine sorts the entries in each column of the
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
C Local variables
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
C Local arrays
      INTEGER TODO(TDLEN)

      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)

C Sort array roughly with partial quicksort
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
        TD = TD - 2
        GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
C          TODO(TD-1) = FIRST
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2

  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
C There is still work to be done
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
        TD = TD - 2
        GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE

  100 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(k). The entries in each column must be
C   non-negative and ordered by decreasing value.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of
C   unmatched columns.
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know
C   a matching does not exist; in the matrix defined by the column
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
      EXTERNAL MC64QD,MC64UD

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped.
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Compute a first maximum matching from scratch on whole matrix.
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
C Matrix is structurally singular
        BMAX = RINF
      ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute
C values of the columns
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF

C Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J.
C Set LENL(J) just after last entry in column J with value ge BMAX.
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE

C Main loop
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
C Set LEN such that all matrix entries with value lt BVAL are
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
              LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81     MOD = 1
        ELSE
C We do not have a maximum matching in IW.
          BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C          IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) =  - J
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL

C This routine searches for at most XX different numerical values
C in the columns W(1:WLEN). XX>=2.
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2).
C
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT
C remains sorted by decreasing value.
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
C Check presence of HA in SPLIT
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
C The insertion
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
C Exit loop if XX values are found
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1

      IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C      IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
      NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1)*N

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

      DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
        ID1 = ID0 + JORD
C J is unmatched column
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row
          ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30     OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
          LAST = JORD
          GO TO 101
        ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
      LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   weighted matching: IPERM(I) = 0 or row I is matched to column
C   IPERM(I).
C NUM is an INTEGER variable. On exit, it contains the cardinality of
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a DOUBLE PRECISION array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and
C   V = D(N+1:2N) contains the dual column variable. If the matrix
C   is structurally nonsingular (NUM = N), the following holds:
C      U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C      U(I) = 0  if IPERM(I) = 0
C      V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
C Initialize U(I)
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
C Row I is not empty
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
C Don't choose cheap assignment from dense columns
        IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
C Assignment of column J to row I
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
      DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
C Continue only if column J is not empty
        IF (K1.GT.K2) GO TO 95
C       VJ = RINF
C Changes made to allow for NaNs
        I0 = IRN(K1)
        VJ = A(K1) - U(I0)
        K0 = K1
        DO 50 K = K1+1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
C Scan remaining part of assigned column JJ
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
        CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE

        DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1

C Scan column that matches with row Q0
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
            DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
C Row I is matched; do not update D(I) if DNEW is larger
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
C Update tree
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

      IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END


C COPYRIGHT (c) 2002 ENSEEIHT-IRIT, Toulouse, France and
C  Council for the Central Laboratory of the Research Councils.
C  Version 1.0.0 July 2004
C  Version 1.0.1 March 2008  Comments reflowed with length < 73
C  AUTHOR Daniel Ruiz (Daniel.Ruiz@enseeiht.fr)
C *** Copyright (c) 2004  Council for the Central Laboratory of the
C     Research Councils and Ecole Nationale Superieure
C     d'Electrotechnique, d'Electronique, d'Informatique,
C     d'Hydraulique et des Telecommunications de Toulouse.          ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC77 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***

      SUBROUTINE MC77I(ICNTL, CNTL)
C
C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC77A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL
      PARAMETER ( LICNTL=10, LCNTL=10 )
      INTEGER ICNTL(LICNTL)
      REAL CNTL(LCNTL)
C
C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the default value, the incoming data is checked for
C     out-of-range indices and duplicates, in which case the driver
C     routine will exit with an error.  Setting ICNTL(4) to any
C     other value will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) nonzero if the data is
C     known to be in range without duplicates.
C
C    ICNTL(5) has default value 0.
C     If left at the default value, it indicates that A contains some
C     negative entries, and it is necessary that their absolute values
C     be computed internally.  Otherwise, the values in the input
C     matrix A will be considered as non-negative.
C
C    ICNTL(6) has default value 0.
C     If nonzero, the input matrix A is symmetric and the user must
C     only supply the lower triangular part of A in the appropriate
C     format.  Entries in the upper triangular part of a symmetric
C     matrix will be considered as out-of-range, and are also checked
C     when ICNTL(4) is 0.
C
C    ICNTL(7) has a default value of 10.
C     It specifies the maximum number of scaling iterations that
C     may be performed.
C     Note that iteration "0", corresponding to the initial
C     normalization of the data, is always performed.
C     Restriction:  ICNTL(7) > 0
C                   (otherwise, the driver stops with an error)
C     ( ... In future release : Restriction:  ICNTL(7) >= 0 ... )
C
C    ICNTL(8) to ICNTL(15) are not currently used by MC77A/AD but are
C     set to zero in this routine.
C
C
C    CNTL(1) has a default value of 0.
C     It specifies the tolerance value when to stopping the iterations,
C     that is it is the desired value such that all row and column norms
C     in the scaled matrix lie between (1 +/- CNTL(1)).
C     If CNTL(1) is less than or equal to 0, tolerance is not checked,
C     and the algorithm will stop when the maximum number of iterations
C     given in ICNTL(7) is reached.
C
C    CNTL(2) has a default value of 1.
C     It is used in conjunction with parameter JOB set to -1,
C     to specify a REAL value for the power of norm under consideration.
C     Restriction:  CNTL(2) >= 1
C                   (otherwise, the driver stops with an error)
C
C    CNTL(3) to CNTL(10) are not currently used by MC77A/AD but are
C     set to zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 10
C     Currently unused control variables:
      DO 10 I = 8,LICNTL
        ICNTL(I) = 0
   10 CONTINUE

C Initialization of the CNTL array.
      CNTL(1) = ZERO
      CNTL(2) = ONE
C     Currently unused control variables:
      DO 20 I = 3,LCNTL
        CNTL(I) = ZERO
   20 CONTINUE

      RETURN
      END

C**********************************************************************
C***           DRIVER FOR THE HARWELL-BOEING SPARSE FORMAT          ***
C**********************************************************************
      SUBROUTINE MC77A(JOB,M,N,NNZ,JCST,IRN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCST(N+1),IRN(NNZ),IW(LIW)
      REAL A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 0.
C
C JCST is an INTEGER array of length N+1.
C   JCST(J), J=1..N, must be set by the user to the position in array
C   IRN of the first row index of an entry in column J.
C   JCST(N+1) must be set to NNZ+1.
C   The array JCST is not altered by the subroutine.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices of
C   the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restrictions:
C     The entries in A belonging to column J must be stored contiguously
C     in the positions JCST(J)..JCST(J+1)-1. The ordering of the row
C     indices within each column is unimportant.
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (REAL in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -8 : repeated entries are found. INFO(2) contains the index
C        of a column in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77N,MC77O,MC77P,MC77Q
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        INFO(2) = 7
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     Harwell-Boeing sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,M
          IW(I) = 0
    3   CONTINUE
        DO 5 J = 1,N
          DO 4 K = JCST(J),JCST(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.M) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
              INFO(1) = -9
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    5   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       Harwell-Boeing sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9021) (JCST(J),J=1,N+1)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in Harwell-Boeing SPARSE format)
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77N(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77N(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77O(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77O(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77P(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77P(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77Q(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77Q(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The Harwell-Boeing sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77A/AD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' JCST(1:N+1) = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77A/AD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77N(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (EW(J).LT.A(K)) THEN
            EW(J) = A(K)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(K)) THEN
            DW(I) = A(K)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IW(I).GT.0) THEN
                S = A(K) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77O(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            I = IRN(K)
            EW(J) = EW(J) + A(K)
            IF (JW(J).EQ.0) THEN
               JW(J) = K
            ELSE
               JW(J) = -1
            ENDIF
            DW(I) = DW(I) + A(K)
            IF (IW(I).EQ.0) THEN
               IW(I) = K
            ELSE
               IW(I) = -1
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,M
        IF (IW(K).NE.0) D(K) = SQRT(DW(K))
   40 CONTINUE
      DO 45 K=1,N
        IF (JW(K).NE.0) E(K) = SQRT(EW(K))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77P(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77Q(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IF (IJW(J).EQ.0) THEN
               IJW(J) = K
            ELSE
               IJW(J) = -1
            ENDIF
            I = IRN(K)
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
              IF (IJW(I).EQ.0) THEN
                 IJW(I) = K
              ELSE
                 IJW(I) = -1
              ENDIF
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***              DRIVER FOR THE GENERAL SPARSE FORMAT              ***
C**********************************************************************
      SUBROUTINE MC77B(JOB,M,N,NNZ,IRN,JCN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCN(NNZ),IRN(NNZ),IW(LIW)
      REAL A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C            If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 1.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices
C   of the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C JCN is an INTEGER array of NNZ.
C   JCN(J), J=1..NNZ, must be set by the user to hold the column indices
C   of the entries in the matrix.
C   The array JCN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (REAL in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -8 : repeated entries are found.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77R,MC77S,MC77T,MC77U
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     General sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 6 K = 1,NNZ
          I = IRN(K)
          J = JCN(K)
C Check for row indices that are out of range
          IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(1) = -7
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
          IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
            INFO(1) = -9
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
    6   CONTINUE
C Check for repeated row indices within a column
        DO 7 I = 1,M
          IW(I) = 0
    7   CONTINUE
        DO 9 J = 1,N
          DO 8 K = 1,NNZ
          IF (JCN(K).EQ.J) THEN
            I = IRN(K)
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = K
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),K,I,J
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
          ENDIF
    8     CONTINUE
    9   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       General sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9022) (JCN(J),J=1,NNZ)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive. The
C   default value used here 1 (e.g. check convergence every iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in general SPARSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77R(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77R(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77S(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77S(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77T(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77T(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77U(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77U(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The general sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Entry ',I8,
     &        ' has invalid row index ',I8, ' or column index ',I8)
 9007 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Duplicate entry ',I8, '   with row index ',I8/
     &        '                                 and column index ',I8)
 9008 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77B/BD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9022 FORMAT (' JCN(1:NNZ)  = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77B/BD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77R(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (EW(J).LT.A(K)) THEN
          EW(J) = A(K)
          JW(J) = I
        ENDIF
        IF (DW(I).LT.A(K)) THEN
          DW(I) = A(K)
          IW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (JW(J).GT.0) THEN
            S = A(K) / (D(I)*E(J))
            IF (EW(J).LT.S) THEN
              EW(J) = S
              JW(J) = I
            ENDIF
            IF (IW(I).GT.0) THEN
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IW(I).GT.0) THEN
              S = A(K) / (D(I)*E(J))
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77S(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          I = IRN(K)
          J = JCN(K)
          EW(J) = EW(J) + A(K)
          IF (JW(J).EQ.0) THEN
             JW(J) = K
          ELSE
             JW(J) = -1
          ENDIF
          DW(I) = DW(I) + A(K)
          IF (IW(I).EQ.0) THEN
             IW(I) = K
          ELSE
             IW(I) = -1
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).NE.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).NE.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (D(I)*E(J))
          EW(J) = EW(J) + S
          DW(I) = DW(I) + S
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77T(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (DEW(J).LT.A(K)) THEN
          DEW(J) = A(K)
          IJW(J) = I
        ENDIF
        IF (DEW(I).LT.A(K)) THEN
          DEW(I) = A(K)
          IJW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (IJW(J).GT.0) THEN
            S = A(K) / (DE(I)*DE(J))
            IF (DEW(J).LT.S) THEN
              DEW(J) = S
              IJW(J) = I
            ENDIF
            IF (IJW(I).GT.0) THEN
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IJW(I).GT.0) THEN
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77U(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          J = JCN(K)
          DEW(J) = DEW(J) + A(K)
          IF (IJW(J).EQ.0) THEN
             IJW(J) = K
          ELSE
             IJW(J) = -1
          ENDIF
          I = IRN(K)
          IF (I.NE.J) THEN
            DEW(I) = DEW(I) + A(K)
            IF (IJW(I).EQ.0) THEN
               IJW(I) = K
            ELSE
               IJW(I) = -1
            ENDIF
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (DE(I)*DE(J))
          DEW(J) = DEW(J) + S
          IF (I.NE.J)  DEW(I) = DEW(I) + S
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***             DRIVER FOR THE CASE OF DENSE MATRICES              ***
C**********************************************************************
      SUBROUTINE MC77C(JOB,M,N,A,LDA,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,LDA,LIW,LDW
      INTEGER IW(LIW)
      REAL A(LDA,*),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C A is a REAL (REAL in the D-version) array containing
C   the numerical values A(i,j), i=1..M, j=1..N, of the input matrix A.
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= (LDA*N) + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= M*(M+1)/2 + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, NNZ, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77J,MC77K,MC77L,MC77M
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check LDA
      IF (ICNTL(6).EQ.0) THEN
        IF (M.GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-M
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),'LDA < M; (LDA-M)',INFO(2)
          GO TO 99
        ENDIF
      ELSE
        IF (((M*(M+1))/2).GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-((M*(M+1))/2)
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),
     &        'LDA < M*(M+1)/2; (LDA-(M*(M+1)/2))',INFO(2)
          GO TO 99
        ENDIF
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      NNZ = (M*(M+1))/2
      IF (ICNTL(6).EQ.0)  THEN
        K = 2*(M+N)
        NNZ = LDA*N
      ENDIF
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          K = 1
          DO 5 I=1,M
C           General case :
            IF (ICNTL(6).EQ.0)
     &        WRITE(ICNTL(3),9021) I,(A(I,J),J=1,N)
C           Dense symmetric packed format :
            IF (ICNTL(6).NE.0) THEN
              WRITE(ICNTL(3),9022) I,(A(J,1),J=K,K+M-I)
              K = K+M-I+1
            ENDIF
    5     CONTINUE
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        IF (ICNTL(6).EQ.0) THEN
          DO 25 J = 1,N
            DO 20 I = 1,M
              DW(K+(J-1)*LDA+I) = ABS(A(I,J))**DP
   20       CONTINUE
   25     CONTINUE
        ELSE
          NNZ = (M*(M+1))/2
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J,1))**DP
   30     CONTINUE
        ENDIF
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          IF (ICNTL(6).EQ.0) THEN
            DO 40 J = 1,N
              DO 35 I = 1,M
                DW(K+(J-1)*LDA+I) = ABS(A(I,J))
   35         CONTINUE
   40       CONTINUE
          ELSE
            NNZ = (M*(M+1))/2
            DO 45 J = 1,NNZ
              DW(K+J) = ABS(A(J,1))
   45       CONTINUE
          ENDIF
        ENDIF
      ENDIF

C Begin the computations (input matrix in DENSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77J(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77J(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77K(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77K(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 50 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   50       CONTINUE
            RINFO(2) = ZERO
            DO 55 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   55       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77L(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77L(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77M(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77M(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9008 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77C/CD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,1:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9022 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,I:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77C/CD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77J(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (EW(J).LT.A(I,J)) THEN
            EW(J) = A(I,J)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(I,J)) THEN
            DW(I) = A(I,J)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 I=1,M
              IF (IW(I).GT.0) THEN
                S = A(I,J) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77K(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (A(I,J).GT.ZERO) THEN
            DW(I) = DW(I) + A(I,J)
            IW(I) = IW(I) + 1
            EW(J) = EW(J) + A(I,J)
            JW(J) = JW(J) + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).GT.0) .OR. (JW(K).GT.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).GT.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).GT.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77L(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      REAL THRESH,ERR
      REAL A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J)  IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  130       CONTINUE
          ELSE
            DO 140 I=J,N
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J)  IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77M(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      REAL THRESH,ERR
      REAL A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IJW(J) = IJW(J) + 1
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
            ENDIF
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
              K = K+1
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END
C COPYRIGHT (c) 2002 ENSEEIHT-IRIT, Toulouse, France and
C  Council for the Central Laboratory of the Research Councils.
C  Version 1.0.0 July 2004
C  Version 1.0.1 March 2008  Comments reflowed with length < 73
C  AUTHOR Daniel Ruiz (Daniel.Ruiz@enseeiht.fr)
C *** Copyright (c) 2004  Council for the Central Laboratory of the
C     Research Councils and Ecole Nationale Superieure
C     d'Electrotechnique, d'Electronique, d'Informatique,
C     d'Hydraulique et des Telecommunications de Toulouse.          ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC77 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***

C**********************************************************************
      SUBROUTINE MC77ID(ICNTL, CNTL)
C
C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC77A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL
      PARAMETER ( LICNTL=10, LCNTL=10 )
      INTEGER ICNTL(LICNTL)
      DOUBLE PRECISION CNTL(LCNTL)
C
C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the default value, the incoming data is checked for
C     out-of-range indices and duplicates, in which case the driver
C     routine will exit with an error.  Setting ICNTL(4) to any
C     other value will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) nonzero if the data is
C     known to be in range without duplicates.
C
C    ICNTL(5) has default value 0.
C     If left at the default value, it indicates that A contains some
C     negative entries, and it is necessary that their absolute values
C     be computed internally.  Otherwise, the values in the input
C     matrix A will be considered as non-negative.
C
C    ICNTL(6) has default value 0.
C     If nonzero, the input matrix A is symmetric and the user must
C     only supply the lower triangular part of A in the appropriate
C     format.  Entries in the upper triangular part of a symmetric
C     matrix will be considered as out-of-range, and are also checked
C     when ICNTL(4) is 0.
C
C    ICNTL(7) has a default value of 10.
C     It specifies the maximum number of scaling iterations that
C     may be performed.
C     Note that iteration "0", corresponding to the initial
C     normalization of the data, is always performed.
C     Restriction:  ICNTL(7) > 0
C                   (otherwise, the driver stops with an error)
C     ( ... In future release : Restriction:  ICNTL(7) >= 0 ... )
C
C    ICNTL(8) to ICNTL(15) are not currently used by MC77A/AD but are
C     set to zero in this routine.
C
C
C    CNTL(1) has a default value of 0.
C     It specifies the tolerance value when to stopping the iterations,
C     that is it is the desired value such that all row and column norms
C     in the scaled matrix lie between (1 +/- CNTL(1)).
C     If CNTL(1) is less than or equal to 0, tolerance is not checked,
C     and the algorithm will stop when the maximum number of iterations
C     given in ICNTL(7) is reached.
C
C    CNTL(2) has a default value of 1.
C     It is used in conjunction with parameter JOB set to -1,
C     to specify a REAL value for the power of norm under consideration.
C     Restriction:  CNTL(2) >= 1
C                   (otherwise, the driver stops with an error)
C
C    CNTL(3) to CNTL(10) are not currently used by MC77A/AD but are
C     set to zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 10
C     Currently unused control variables:
      DO 10 I = 8,LICNTL
        ICNTL(I) = 0
   10 CONTINUE

C Initialization of the CNTL array.
      CNTL(1) = ZERO
      CNTL(2) = ONE
C     Currently unused control variables:
      DO 20 I = 3,LCNTL
        CNTL(I) = ZERO
   20 CONTINUE

      RETURN
      END

C**********************************************************************
C***           DRIVER FOR THE HARWELL-BOEING SPARSE FORMAT          ***
C**********************************************************************
      SUBROUTINE MC77AD(JOB,M,N,NNZ,JCST,IRN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCST(N+1),IRN(NNZ),IW(LIW)
      DOUBLE PRECISION A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 0.
C
C JCST is an INTEGER array of length N+1.
C   JCST(J), J=1..N, must be set by the user to the position in array
C   IRN of the first row index of an entry in column J.
C   JCST(N+1) must be set to NNZ+1.
C   The array JCST is not altered by the subroutine.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices of
C   the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restrictions:
C     The entries in A belonging to column J must be stored contiguously
C     in the positions JCST(J)..JCST(J+1)-1. The ordering of the row
C     indices within each column is unimportant.
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (DOUBLE PRECISION in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (DOUBLE PRECISION in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (DOUBLE PRECISION in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -8 : repeated entries are found. INFO(2) contains the index
C        of a column in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (DOUBLE PRECISION in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77ND,MC77OD,MC77PD,MC77QD
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        INFO(2) = 7
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     Harwell-Boeing sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,M
          IW(I) = 0
    3   CONTINUE
        DO 5 J = 1,N
          DO 4 K = JCST(J),JCST(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.M) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
              INFO(1) = -9
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    5   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       Harwell-Boeing sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9021) (JCST(J),J=1,N+1)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in Harwell-Boeing SPARSE format)
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77ND(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77ND(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77OD(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77OD(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77PD(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77PD(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77QD(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77QD(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The Harwell-Boeing sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77A/AD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' JCST(1:N+1) = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77A/AD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77ND(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (EW(J).LT.A(K)) THEN
            EW(J) = A(K)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(K)) THEN
            DW(I) = A(K)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IW(I).GT.0) THEN
                S = A(K) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77OD(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            I = IRN(K)
            EW(J) = EW(J) + A(K)
            IF (JW(J).EQ.0) THEN
               JW(J) = K
            ELSE
               JW(J) = -1
            ENDIF
            DW(I) = DW(I) + A(K)
            IF (IW(I).EQ.0) THEN
               IW(I) = K
            ELSE
               IW(I) = -1
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,M
        IF (IW(K).NE.0) D(K) = SQRT(DW(K))
   40 CONTINUE
      DO 45 K=1,N
        IF (JW(K).NE.0) E(K) = SQRT(EW(K))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77PD(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77QD(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IF (IJW(J).EQ.0) THEN
               IJW(J) = K
            ELSE
               IJW(J) = -1
            ENDIF
            I = IRN(K)
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
              IF (IJW(I).EQ.0) THEN
                 IJW(I) = K
              ELSE
                 IJW(I) = -1
              ENDIF
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***              DRIVER FOR THE GENERAL SPARSE FORMAT              ***
C**********************************************************************
      SUBROUTINE MC77BD(JOB,M,N,NNZ,IRN,JCN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCN(NNZ),IRN(NNZ),IW(LIW)
      DOUBLE PRECISION A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C            If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 1.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices
C   of the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C JCN is an INTEGER array of NNZ.
C   JCN(J), J=1..NNZ, must be set by the user to hold the column indices
C   of the entries in the matrix.
C   The array JCN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (DOUBLE PRECISION in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (DOUBLE PRECISION in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (DOUBLE PRECISION in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -8 : repeated entries are found.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (DOUBLE PRECISION in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77RD,MC77SD,MC77TD,MC77UD
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     General sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 6 K = 1,NNZ
          I = IRN(K)
          J = JCN(K)
C Check for row indices that are out of range
          IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(1) = -7
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
          IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
            INFO(1) = -9
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
    6   CONTINUE
C Check for repeated row indices within a column
        DO 7 I = 1,M
          IW(I) = 0
    7   CONTINUE
        DO 9 J = 1,N
          DO 8 K = 1,NNZ
          IF (JCN(K).EQ.J) THEN
            I = IRN(K)
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = K
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),K,I,J
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
          ENDIF
    8     CONTINUE
    9   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       General sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9022) (JCN(J),J=1,NNZ)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive. The
C   default value used here 1 (e.g. check convergence every iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in general SPARSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77RD(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77RD(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77SD(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77SD(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77TD(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77TD(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77UD(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77UD(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The general sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Entry ',I8,
     &        ' has invalid row index ',I8, ' or column index ',I8)
 9007 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Duplicate entry ',I8, '   with row index ',I8/
     &        '                                 and column index ',I8)
 9008 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77B/BD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9022 FORMAT (' JCN(1:NNZ)  = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77B/BD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77RD(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (EW(J).LT.A(K)) THEN
          EW(J) = A(K)
          JW(J) = I
        ENDIF
        IF (DW(I).LT.A(K)) THEN
          DW(I) = A(K)
          IW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (JW(J).GT.0) THEN
            S = A(K) / (D(I)*E(J))
            IF (EW(J).LT.S) THEN
              EW(J) = S
              JW(J) = I
            ENDIF
            IF (IW(I).GT.0) THEN
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IW(I).GT.0) THEN
              S = A(K) / (D(I)*E(J))
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77SD(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          I = IRN(K)
          J = JCN(K)
          EW(J) = EW(J) + A(K)
          IF (JW(J).EQ.0) THEN
             JW(J) = K
          ELSE
             JW(J) = -1
          ENDIF
          DW(I) = DW(I) + A(K)
          IF (IW(I).EQ.0) THEN
             IW(I) = K
          ELSE
             IW(I) = -1
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).NE.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).NE.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (D(I)*E(J))
          EW(J) = EW(J) + S
          DW(I) = DW(I) + S
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77TD(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (DEW(J).LT.A(K)) THEN
          DEW(J) = A(K)
          IJW(J) = I
        ENDIF
        IF (DEW(I).LT.A(K)) THEN
          DEW(I) = A(K)
          IJW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (IJW(J).GT.0) THEN
            S = A(K) / (DE(I)*DE(J))
            IF (DEW(J).LT.S) THEN
              DEW(J) = S
              IJW(J) = I
            ENDIF
            IF (IJW(I).GT.0) THEN
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IJW(I).GT.0) THEN
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77UD(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          J = JCN(K)
          DEW(J) = DEW(J) + A(K)
          IF (IJW(J).EQ.0) THEN
             IJW(J) = K
          ELSE
             IJW(J) = -1
          ENDIF
          I = IRN(K)
          IF (I.NE.J) THEN
            DEW(I) = DEW(I) + A(K)
            IF (IJW(I).EQ.0) THEN
               IJW(I) = K
            ELSE
               IJW(I) = -1
            ENDIF
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (DE(I)*DE(J))
          DEW(J) = DEW(J) + S
          IF (I.NE.J)  DEW(I) = DEW(I) + S
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***             DRIVER FOR THE CASE OF DENSE MATRICES              ***
C**********************************************************************
      SUBROUTINE MC77CD(JOB,M,N,A,LDA,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,LDA,LIW,LDW
      INTEGER IW(LIW)
      DOUBLE PRECISION A(LDA,*),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C A is a REAL (DOUBLE PRECISION in the D-version) array containing
C   the numerical values A(i,j), i=1..M, j=1..N, of the input matrix A.
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= (LDA*N) + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= M*(M+1)/2 + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (DOUBLE PRECISION in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (DOUBLE PRECISION in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (DOUBLE PRECISION in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, NNZ, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77JD,MC77KD,MC77LD,MC77MD
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check LDA
      IF (ICNTL(6).EQ.0) THEN
        IF (M.GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-M
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),'LDA < M; (LDA-M)',INFO(2)
          GO TO 99
        ENDIF
      ELSE
        IF (((M*(M+1))/2).GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-((M*(M+1))/2)
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),
     &        'LDA < M*(M+1)/2; (LDA-(M*(M+1)/2))',INFO(2)
          GO TO 99
        ENDIF
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      NNZ = (M*(M+1))/2
      IF (ICNTL(6).EQ.0)  THEN
        K = 2*(M+N)
        NNZ = LDA*N
      ENDIF
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          K = 1
          DO 5 I=1,M
C           General case :
            IF (ICNTL(6).EQ.0)
     &        WRITE(ICNTL(3),9021) I,(A(I,J),J=1,N)
C           Dense symmetric packed format :
            IF (ICNTL(6).NE.0) THEN
              WRITE(ICNTL(3),9022) I,(A(J,1),J=K,K+M-I)
              K = K+M-I+1
            ENDIF
    5     CONTINUE
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        IF (ICNTL(6).EQ.0) THEN
          DO 25 J = 1,N
            DO 20 I = 1,M
              DW(K+(J-1)*LDA+I) = ABS(A(I,J))**DP
   20       CONTINUE
   25     CONTINUE
        ELSE
          NNZ = (M*(M+1))/2
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J,1))**DP
   30     CONTINUE
        ENDIF
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          IF (ICNTL(6).EQ.0) THEN
            DO 40 J = 1,N
              DO 35 I = 1,M
                DW(K+(J-1)*LDA+I) = ABS(A(I,J))
   35         CONTINUE
   40       CONTINUE
          ELSE
            NNZ = (M*(M+1))/2
            DO 45 J = 1,NNZ
              DW(K+J) = ABS(A(J,1))
   45       CONTINUE
          ENDIF
        ENDIF
      ENDIF

C Begin the computations (input matrix in DENSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77JD(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77JD(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77KD(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77KD(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 50 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   50       CONTINUE
            RINFO(2) = ZERO
            DO 55 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   55       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77LD(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77LD(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77MD(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77MD(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9008 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77C/CD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,1:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9022 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,I:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77C/CD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77JD(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (EW(J).LT.A(I,J)) THEN
            EW(J) = A(I,J)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(I,J)) THEN
            DW(I) = A(I,J)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 I=1,M
              IF (IW(I).GT.0) THEN
                S = A(I,J) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77KD(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (A(I,J).GT.ZERO) THEN
            DW(I) = DW(I) + A(I,J)
            IW(I) = IW(I) + 1
            EW(J) = EW(J) + A(I,J)
            JW(J) = JW(J) + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).GT.0) .OR. (JW(K).GT.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).GT.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).GT.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77LD(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J)  IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  130       CONTINUE
          ELSE
            DO 140 I=J,N
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J)  IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77MD(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IJW(J) = IJW(J) + 1
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
            ENDIF
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
              K = K+1
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

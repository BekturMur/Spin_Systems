CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
900   DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END
      SUBROUTINE RJBESL(X, ALPHA, NB, B, NCALC)
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C  Explanation of machine-dependent constants
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C
C  Latest modification: March 19, 1990
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
      INTEGER I,J,K,L,M,MAGX,N,NB,NBMX,NCALC,NEND,NSTART
CS    REAL               GAMMA,
      DOUBLE PRECISION  DGAMMA,
     1 ALPHA,ALPEM,ALP2EM,B,CAPP,CAPQ,CONV,EIGHTH,EM,EN,ENMTEN,ENSIG,
     2 ENTEN,FACT,FOUR,FUNC,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     3 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,TEMPC,TEST,
     4 THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,TWOPI2,X,XC,XIN,XK,XLARGE,
     5 XM,VCOS,VSIN,Z,ZERO
      DIMENSION B(NB), FACT(25)
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
CS    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
CS   1 1.935307179586476925286767E-3/
CS    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
CS    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
CS    DATA ONE30, THREE5 /130.0E0,35.0E0/
      DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535D0,6.28125D0,
     1 1.935307179586476925286767D-3/
      DATA ZERO, EIGHTH, HALF, ONE /0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO, THREE, FOUR, TWOFIV /2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30, THREE5 /130.0D0,35.0D0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/
      DATA ENTEN, ENSIG, RTNSIG /1.0D38,1.0D17,1.0D-4/
      DATA ENMTEN, XLARGE /1.2D-37,1.0D4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA FACT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     1 4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     2 8.71782912D10,1.307674368D12,2.0922789888D13,3.55687428096D14,
     3 6.402373705728D15,1.21645100408832D17,2.43290200817664D18,
     4 5.109094217170944D19,1.12400072777760768D21,
     5 2.585201673888497664D22,6.2044840173323943936D23/
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(X) = GAMMA(X)
      CONV(I) = DBLE(I)
      FUNC(X) = DGAMMA(X)
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) 
     1       .AND. (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE))  
     2   THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
            NCALC = NB
            DO 20 I=1,NB
              B(I) = ZERO
   20       CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
            IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
               TEMPA = ONE
               ALPEM = ONE + ALPHA
               HALFX = ZERO
               IF (X.GT.ENMTEN) HALFX = HALF*X
               IF (ALPHA.NE.ZERO)
     1            TEMPA = HALFX**ALPHA/(ALPHA*FUNC(ALPHA))
               TEMPB = ZERO
               IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
               B(1) = TEMPA + TEMPA*TEMPB/ALPEM
               IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
               IF (NB .NE. 1) THEN
                  IF (X .LE. ZERO) THEN
                        DO 30 N=2,NB
                          B(N) = ZERO
   30                   CONTINUE
                     ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                        TEMPC = HALFX
                        TOVER = (ENMTEN+ENMTEN)/X
                        IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                        DO 50 N=2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N))
     1                       NCALC = N-1
   50                   CONTINUE
                  END IF
               END IF
            ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
               XC = SQRT(PI2/X)
               XIN = (EIGHTH/X)**2
               M = 11
               IF (X.GE.THREE5) M = 8
               IF (X.GE.ONE30) M = 4
               XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
               T = AINT(X/(TWOPI1+TWOPI2)+HALF)
               Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
               VSIN = SIN(Z)
               VCOS = COS(Z)
               GNU = ALPHA + ALPHA
               DO 80 I=1,2
                 S = ((XM-ONE)-GNU)*((XM-ONE)+GNU)*XIN*HALF
                 T = (GNU-(XM-THREE))*(GNU+(XM-THREE))
                 CAPP = S*T/FACT(2*M+1)
                 T1 = (GNU-(XM+ONE))*(GNU+(XM+ONE))
                 CAPQ = S*T1/FACT(2*M+2)
                 XK = XM
                 K = M + M
                 T1 = T
                 DO 70 J=2,M
                   XK = XK - FOUR
                   S = ((XK-ONE)-GNU)*((XK-ONE)+GNU)
                   T = (GNU-(XK-THREE))*(GNU+(XK-THREE))
                   CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                   CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                   K = K - 2
                   T1 = T
   70            CONTINUE
                 CAPP = CAPP + ONE
                 CAPQ = (CAPQ+ONE)*(GNU*GNU-ONE)*(EIGHTH/X)
                 B(I) = XC*(CAPP*VCOS-CAPQ*VSIN)
                 IF (NB.EQ.1) GO TO 300
                 T = VSIN
                 VSIN = -VCOS
                 VCOS = T
                 GNU = GNU + TWO
   80         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
               IF (NB .GT. 2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 90 J=3,NB
                    B(J) = GNU*B(J-1)/X - B(J-2)
                    GNU = GNU + TWO
   90             CONTINUE
               END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
            ELSE
               NBMX = NB - MAGX
               N = MAGX + 1
               EN = CONV(N+N) + (ALPHA+ALPHA)
               PLAST = ONE
               P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
               TEST = ENSIG + ENSIG
               IF (NBMX .GE. 3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 130 K=NSTART,NEND
                     N = K
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN*PLAST/X - POLD
                     IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                        TOVER = ENTEN
                        P = P/TOVER
                        PLAST = PLAST/TOVER
                        PSAVE = P
                        PSAVEL = PLAST
                        NSTART = N + 1
  100                   N = N + 1
                           EN = EN + TWO
                           POLD = PLAST
                           PLAST = P
                           P = EN*PLAST/X - POLD
                        IF (P.LE.ONE) GO TO 100
                        TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                        TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))
                        TEST = TEST/ENSIG
                        P = PLAST*TOVER
                        N = N - 1
                        EN = EN - TWO
                        NEND = MIN(NB,N)
                        DO 110 L=NSTART,NEND
                           POLD = PSAVEL
                           PSAVEL = PSAVE
                           PSAVE = EN*PSAVEL/X - POLD
                           IF (PSAVE*PSAVEL.GT.TEST) THEN
                              NCALC = L - 1
                              GO TO 190
                           END IF
  110                   CONTINUE
                        NCALC = NEND
                        GO TO 190
                     END IF
  130             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
               END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  140          N = N + 1
                  EN = EN + TWO
                  POLD = PLAST
                  PLAST = P
                  P = EN*PLAST/X - POLD
               IF (P.LT.TEST) GO TO 140
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  190          N = N + 1
               EN = EN + TWO
               TEMPB = ZERO
               TEMPA = ONE/P
               M = 2*N - 4*(N/2)
               SUM = ZERO
               EM = CONV(N/2)
               ALPEM = (EM-ONE) + ALPHA
               ALP2EM = (EM+EM) + ALPHA
               IF (M .NE. 0) SUM = TEMPA*ALPEM*ALP2EM/EM
               NEND = N - NB
               IF (NEND .GT. 0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 200 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     TEMPC = TEMPB
                     TEMPB = TEMPA
                     TEMPA = (EN*TEMPB)/X - TEMPC
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        IF (N.EQ.1) GO TO 210
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                     END IF
  200             CONTINUE
               END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  210          B(N) = TEMPA
               IF (NEND .GE. 0) THEN
                  IF (NB .LE. 1) THEN
                        ALP2EM = ALPHA
                        IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                        SUM = SUM + B(1)*ALP2EM
                        GO TO 250
                     ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*TEMPA)/X - TEMPB
                        IF (N.EQ.1) GO TO 240
                        M = 2 - M
                        IF (M .NE. 0) THEN
                           EM = EM - ONE
                           ALP2EM = (EM+EM) + ALPHA
                           ALPEM = (EM-ONE) + ALPHA
                           IF (ALPEM.EQ.ZERO) ALPEM = ONE
                           SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                        END IF
                  END IF
               END IF
               NEND = N - 2
               IF (NEND .NE. 0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 230 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     B(N) = (EN*B(N+1))/X - B(N+2)
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                     END IF
  230             CONTINUE
               END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
               B(1) = TWO*(ALPHA+ONE)*B(2)/X - B(3)
  240          EM = EM - ONE
               ALP2EM = (EM+EM) + ALPHA
               IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
               SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  250          IF ((ALPHA+ONE).NE.ONE)
     1              SUM = SUM*FUNC(ALPHA)*(X*HALF)**(-ALPHA)
               TEMPA = ENMTEN
               IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
               DO 260 N=1,NB
                 IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                 B(N) = B(N)/SUM
  260          CONTINUE
            END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
         ELSE
            B(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  300 RETURN
C ---------- Last line of RJBESL ----------
      END


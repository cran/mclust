      double precision function det2mc( n, u, s)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c safeguarded computation for determinant of triangular matrix r

      implicit double precision (a-h,o-z)

      double precision                 u(n,*), s

      double precision                 zero, two 
      parameter                       (zero = 0.d0, two = 2.d0)

      double precision                 FLMAX
      common /MCLMCH/                  FLMAX
      save   /MCLMCH/

      det2mc = zero

      do k = 1, n

        q = u(k,k)*s

        if (q .eq. zero) then
          det2mc = -FLMAX
          return
        end if

        det2mc = det2mc + log(abs(q))
                 
      end do

      det2mc = two*det2mc

      return  
      end 

      double precision function detmc2( n, u)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c safeguarded computation for determinant of triangular matrix u

      implicit double precision (a-h,o-z)

      double precision                 u(n,*)

      double precision                 zero, two 
      parameter                       (zero = 0.d0, two = 2.d0)

      double precision                 FLMAX
      common /MCLMCH/                  FLMAX
      save   /MCLMCH/

      detmc2 = zero

      do k = 1, n

        q = u(k,k)

        if (q .eq. zero) then
          detmc2 = -FLMAX
          return
        end if

        detmc2 = detmc2 + log(abs(q))
                 
      end do

      detmc2 = two*detmc2

      return  
      end 
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTF2 computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U'*U  or A = L*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of DPOTF2
*
      END
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRF computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U'*U.
*
            DO 10 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L'.
*
            DO 20 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block column.
*
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOTRF
*
      END
      subroutine drnge( l, v, i, vmin, vmax)

c finds the max and min elements in absolute value of a vector

      implicit double precision (a-h,o-z)

      double precision v(*)

c----------------------------------------------------------------------------

      temp  = abs(v(1))
      vmin  = temp
      vmax  = temp

      if (l .eq. 1) return

      if (i .eq. 1) then
        do j = 2, l
          temp = abs(v(j))
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
        end do
      else
        k = 1 + i
        do j = 2, l
          temp = abs(v(k))
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
          k    = k + i
        end do
      end if

      return  
      end 
      subroutine eseee ( x, mu, Sigma, prob, n, p, G, w, z, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for constant-variance Gaussian mixtures

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

c     double precision   x(n,p),mu(p,G),Sigma(p,p),prob(G),w(p),z(n,G)
      double precision   x(n,*), mu(p,*), Sigma(p,*), prob(*), 
     *                   w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group.
c  Sigma   double  (input) (p,p) covariance matrix. Destroyed on output.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  w       double  (scratch) (p) workspace
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound condition number
c                   of sigma. On output, the loglikelihood.

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps   = sqrt(max(hood,zero))

      p1    = p + 1

      FLMIN = d1mach(1)
      FLMAX = d1mach(2)

      call dpotrf( 'L', p, Sigma, p, info)
      if (info .ne. 0) then
        hood = dble(sign(1,info))*FLMIN
        return
      end if

      call drnge( p, Sigma, p1, umin, umax)

      rc = umin/(one+umax)

      if (rc .le. eps) then
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = dble(p)*pi2log/two + detlog

      if (G .lt. 0) goto 500

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob(k) * exp(-(const+temp))
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      w(1) = umin
      w(2) = umax

      return

500   continue

      G     = -G

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = exp(-(const+temp))
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      w(1) = umin
      w(2) = umax

      return
      end
      subroutine esei ( x, mu, sigsq, prob, n, p, G, z, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step : computes conditional probabilities and loglikelihood 
c for spherical, constant-volume Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), mu(p,G), prob(G), z(n,G)
      double precision   x(n,*), mu(p,*), prob(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  mu      double  (input) (p,G) means for each group.
c  sigsq   double  (input) Sigma-squared.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound for sigsq.
c                   On output, the loglikelihood.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        d1mach
      external                d1mach

c------------------------------------------------------------------------------

      eps   = max(hood,zero)

      FLMAX = d1mach(2)

      if (sigsq .le. eps) then
        hood = FLMAX
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      if (G .lt. 0) goto 500

      hood = zero
      do i = 1, n
        sumz = zero
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          temp   = prob(k)*exp(-(const+sum/sigsq)/two)
          z(i,k) = temp
          sumz   = sumz + temp
        end do
        hood = hood + log(sumz)
        call dscal( G, (one/sumz), z(i,1), n)
      end do

      return

500   continue

      G    = -G

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sumz = zero
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          temp   = exp(-(const+sum/sigsq)/two)
          z(i,k) = temp
          sumz   = sumz + temp
        end do
        hood = hood + log(sumz)
        call dscal( G, (one/sumz), z(i,1), n)
      end do

      return
      end
      subroutine esneee( x, mu, Sigma, prob, n, p, G, w, z, hood, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for constant-variance Gaussian mixtures with Poisson noise

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

c     double precision   x(n,p),mu(p,G),Sigma(p,p),prob(G),w(p),z(n,G)
      double precision   x(n,*), mu(p,*), Sigma(p,*), prob(*), 
     *                   w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group.
c  Sigma   double  (input) (p,p) covariance matrix. Destroyed on output.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  w       double  (scratch) (p) workspace
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound condition number
c                   of sigma. On output, the loglikelihood.
c  Vinv    double  (input) estimated reciprocal hypervolume of the data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps   = sqrt(max(hood,zero))

      p1    = p + 1

      FLMIN = d1mach(1)
      FLMAX = d1mach(2)

      call dpotrf( 'L', p, Sigma, p, info)
      if (info .ne. 0) then
        hood = dble(sign(1,info))*FLMIN
        return
      end if

      call drnge( p, Sigma, p1, umin, umax)
 
      rc = umin/(one+umax)

      if (rc .le. eps) then
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = dble(p)*pi2log/two + detlog

      if (G .lt. 0) goto 500

      G1    = G + 1

      termn = prob(G1)*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob(k) * exp(-(const+temp))
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      w(1) = umin
      w(2) = umax

      return

500   continue

      G  = -G
      G1 =  G + 1

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = exp(-(const+temp))
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      w(1) = umin
      w(2) = umax

      return
      end
      subroutine esnei( x, mu, sigsq, prob, n, p, G, z, hood, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step : computes conditional probabilities and logliklihood 
c for spherical, constant-volume Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), mu(p,G), prob(G+1), z(n,G+1)
      double precision   x(n,*), mu(p,*), prob(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  mu      double  (input) (p,G) means for each group.
c  sigsq   double  (input) Sigma-squared.
c  prob    double  (input)(G+1) probabilities (not needed if equal proportions)
c  z       double  (output) (n,G+1) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound for sigsq.
c                   On output, the loglikelihood.
c  Vinv    double  (input) estimated reciprocal hypervolume of the data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps   = max(hood,zero)

      FLMAX = d1mach(2)

      if (sigsq .lt. eps) then
        hood = FLMAX
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      if (G .lt. 0) goto 100
   
      G1    =  G + 1

      termn = prob(G1)*Vinv

      hood = zero
      do i = 1, n
        sumz = termn
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          temp   = prob(k)*exp(-(const+sum/sigsq)/two)
          z(i,k) = temp
          sumz   = sumz + temp
        end do
        z(i,G1) = termn
        hood = hood + log(sumz)
        call dscal( G1, (one/sumz), z(i,1), n)
      end do

      return

100   continue

      G    = -G

      G1   = G + 1

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sumz = Vinv
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          temp   = exp(-(const+sum/sigsq)/two)
          z(i,k) = temp
          sumz   = sumz + temp
        end do
        z(i,G1) = Vinv
        hood = hood + log(sumz)
        call dscal( G1, (one/sumz), z(i,1), n)
      end do

      return
      end
      subroutine esnvi( x, mu, sigsq, prob, n, p, G, z, hood, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step : computes conditional probabilities and loglikelihood 
c for spherical, varying-volume Gaussian mixtures with Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), mu(p,G), sigsq(G), prob(G+1), z(n,G+1)
      double precision   x(n,*), mu(p,*), sigsq(*), prob(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  mu      double  (input) (p,G) means for each group.
c  sigsq   double  (input) (G) Sigma-squared.
c  prob    double  (input)(G+1) probabilities (not needed if equal proportions)
c  z       double  (output) (n,G+1) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound for sigsq.
c                   On output, the loglikelihood.
c  Vinv    double  (input) estimated reciprocal hypervolume of the data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        d1mach
      external                d1mach

c------------------------------------------------------------------------------

      eps   = max(hood,zero)

      FLMAX = d1mach(2)

      if (G .lt. 0) goto 500

      G1 = G + 1

      do k = 1, G
        sigsqk = sigsq(k)
        if (sigsqk .le. eps) then
          hood = FLMAX
          return
        end if
        const = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = prob(k)*exp(-(const+sum/sigsqk)/two)
        end do
      end do

      termn = prob(G1)*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      return

500   continue

      G  = -G
      G1 =  G + 1

      do k = 1, G
        sigsqk = sigsq(k)
        if (sigsqk .le. eps) then
          hood = FLMAX
          return
        end if
        const = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = exp(-(const+sum/sigsqk)/two)
        end do
      end do

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      return
      end
      subroutine esnvvv ( x, mu, Sigma, prob, n, p, G, w, z, hood, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for unconstrained Gaussian mixtures with Poisson noise

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

c     double precision   x(n,p),mu(p,G),Sigma(p,p,G),prob(G+1),w(p),z(n,G+1)
      double precision   x(n,*), mu(p,*), Sigma(p,p,*), prob(*), 
     *                   w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group.
c  Sigma   double  (input) (p,p,G) covariance matrices. Destroyed on output.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  w       double  (scratch) (p) workspace
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound condition number
c                   of cholesky factor of sigma. On output, the loglikelihood.
c  Vinv    double  (input) estimated reciprocal hypervolume of the data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps   = sqrt(max(hood,zero))

      p1    = p + 1

      FLMIN = d1mach(1)
      FLMAX = d1mach(2)

      if (G .lt. 0) goto 500

      G1 = G + 1

      do k = 1, G

        call dpotrf( 'L', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          hood = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, Sigma(1,1,k), p1, umin, umax)

        rc = umin/(one+umax)

        if (rc .le. eps) then
          hood = FLMAX
          return
        end if

        detlog = log(abs(Sigma(1,1,k)))
        do j = 2, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma(1,1,k), p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob(k) * exp(-(const+temp))
        end do

      end do

      termn = prob(G1)*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      return

500   continue

      G  = -G

      G1 =  G + 1

      const  = const/dble(G1)

      do k = 1, G

        call dpotrf( 'L', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          hood = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, Sigma(1,1,k), p1, umin, umax)

        rc = umin/(one+umax)

        if (rc .le. eps) then
          hood = FLMAX
          return
        end if

        detlog = log(abs(Sigma(1,1,k)))
        do j = 2, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma(1,1,k), p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = exp(-(const+temp))
        end do

      end do

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      return
      end
      subroutine esnxev ( x, mu, Sigma, prob, n, p, G, eps,
     *                    s, w, lwork, z, hood, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for prescribed-shape Gaussian mixtures with Poisson noise

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p),mu(p,G),Sigma(p,p,G),prob(G+1),
c    *                   s(p),w(lwork),z(n,G+1)
      double precision   x(n,*),mu(p,*),Sigma(p,p,*),prob(*),
     *                   s(*), w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group. On output, contains
c                   the normalized shape vector for each group.
c  Sigma   double  (input) (p,p,G) covariance matrices. On output, the first
c                   element in each matrix contains lambda for that group.
c  prob    double  (input) (G+1) probabilities
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape.
c                   On output, minimum values of these over the groups.
c  s       double  (scratch) (p)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) max(4*p,5*p-4) minumum workspace for LAPACK SVD.
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) logliklihood.
c  Vinv    double  (input) estimated reciprocal hypervolume of the data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps1   = max(eps(1),zero)
      eps2   = sqrt(max(eps(2),zero))

      p1     = p + 1

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
      EPSMIN = d1mach(3)

      vlamin = FLMAX
      rcmin  = FLMAX

      if (G .lt. 0) goto 500

      G1 = G + 1

      do k = 1, G

        call dpotrf( 'U', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood = dble(sign(1,info))*EPSMIN
          return
        end if

c zero out the lower triangle
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, Sigma(j,i,k),1)
          i = j
        end do

        call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, s, 1, smin, smax)

        rcmin = min(smin/(one+smax),rcmin)

        if (rcmin .le. eps2) then
          eps(1) = vlamin
          eps(2) = rcmin*rcmin
          hood   = -FLMAX
          return
        end if

        smax = zero
        sum  = zero
        do j = 1, p
          temp = s(j)
          smax = max(smax,temp)
          sum  = sum + log(temp*temp)
        end do
        vlam = exp(sum/dble(p))
        vlamin = min(vlamin,vlam)
        if (vlam .gt. eps1) then
          const = dble(p)*pi2log+sum
          slam = sqrt(vlam)
          call dscal( p, one/slam, s, 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w(p1), 1)
            call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
            call dgemv('N',p,p,one,Sigma(1,1,k),p,w(p1),1,zero,w,1)
            do j = 1, p
              w(j) = w(j) / s(j)
            end do
            temp   = ddot(p,w,1,w,1)/vlam
            z(i,k) = prob(k) * exp(-(const+temp)/two)
          end do
          call dscal( p, slam, s, 1)
        else
        end if
        Sigma(1,1,k) = vlam
        call dscal( p, one/smax, s, 1)
        call dcopy( p, s, 1, mu(1,k), 1)
      end do

      if (vlamin .le. eps1) then
        eps(1) = vlamin
        eps(2) = rcmin*rcmin
        hood   = FLMAX
        return
      end if

      termn = prob(G1)*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do

      eps(1) = vlamin
      eps(2) = rcmin*rcmin

      return

500   continue

      G  = -G
      G1 =  G + 1

      do k = 1, G

        call dpotrf( 'U', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*EPSMIN
          return
        end if

c zero out the lower triangle
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, Sigma(j,i,k),1)
          i = j
        end do

        call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*FLMIN
          return
        end if

        smax = zero
        sum  = zero
        do j = 1, p
          temp = s(j)
          smax = max(smax,temp)
          sum  = sum + log(temp*temp)
        end do
        vlam = exp(sum/dble(p))
        vlamin = min(vlamin,vlam)
        if (vlam .gt. eps1) then
          const = dble(p)*pi2log+sum
          slam = sqrt(vlam)
          call dscal( p, one/slam, s, 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w(p1), 1)
            call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
            call dgemv('N',p,p,one,Sigma(1,1,k),p,w(p1),1,zero,w,1)
            do j = 1, p
              w(j) = w(j) / s(j)
            end do
            temp   = ddot(p,w,1,w,1)/vlam
            z(i,k) = exp(-(const+temp)/two)
          end do
          call dscal( p, slam, s, 1)
        else
        end if
        Sigma(1,1,k) = vlam
        call dscal( p, one/smax, s, 1)
        call dcopy( p, s, 1, mu(1,k), 1)
      end do

      if (vlamin .le. eps1) then
        eps(1) = vlamin
        eps(2) = rcmin*rcmin
        hood   = FLMAX
        return
      end if

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
   
      eps(1) = vlamin
      eps(2) = rcmin*rcmin

      return
      end
      subroutine esvi ( x, mu, sigsq, prob, n, p, G, z, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step : computes conditional probabilities and loglikelihood 
c for spherical, varying-volume Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), mu(p,G), sigsq(G), prob(G), z(n,G)
      double precision   x(n,*), mu(p,*), sigsq(*), prob(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  mu      double  (input) (p,G) means for each group.
c  sigsq   double  (input) (G) Sigma-squared.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound for sigsq.
c                   On output, the loglikelihood.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        d1mach
      external                d1mach

c------------------------------------------------------------------------------

      eps    = max(hood,zero)

      FLMAX  = d1mach(2)

      if (G .lt. 0) goto 500

      do k = 1, G
        sigsqk = sigsq(k)
        if (sigsqk .le. eps) then
          hood = FLMAX
          return
        end if
        const = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = prob(k)*exp(-(const+sum/sigsqk)/two)
        end do
      end do

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      return

500   continue

      G = -G

      do k = 1, G
        sigsqk = sigsq(k)
        if (sigsqk .le. eps) then
          hood = FLMAX
          return
        end if
        const = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = exp(-(const+sum/sigsqk)/two)
        end do
      end do

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      return
      end
      subroutine esvvv ( x, mu, Sigma, prob, n, p, G,
     *                   w, z, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for unconstrained Gaussian mixtures

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

c     double precision   x(n,p),mu(p,G),Sigma(p,p,G),prob(G),w(p),z(n,G)
      double precision   x(n,*), mu(p,*), Sigma(p,p,*), prob(*), 
     *                   w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group. Detroyed on exit.
c  Sigma   double  (input) (p,p,G) covariance matrices. Destroyed on output.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  detlog  double  (scratch) (G) workspace
c  w       double  (scratch) (p) workspace
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (input/output) On input, lower bound condition number
c                   of sigma. On output, the loglikelihood.

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps    = sqrt(hood)

      p1     = p + 1

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)

      if (G .lt. 0) goto 500

      do k = 1, G

        call dpotrf( 'L', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          hood = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, Sigma(1,1,k), p1, umin, umax)

        rc = umin/(one+umax)

        if (rc .le. eps) then
          hood = FLMAX
          return
        end if

        detlog = log(abs(Sigma(1,1,k)))
        do j = 2, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma(1,1,k), p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob(k) * exp(-(const+temp))
        end do

      end do

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      return

500   continue

      G = -G

      do k = 1, G
        call dpotrf( 'L', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          hood = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, Sigma(1,1,k), p1, umin, umax)

        rc = umin/(one+umax)

        if (rc .le. eps) then
          hood = FLMAX
          return
        end if

        detlog = log(abs(Sigma(1,1,k)))
        do j = 2, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'L', 'N', 'N', p, Sigma(1,1,k), p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = exp(-(const+temp))
        end do

      end do

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      return
      end
      subroutine esxev ( x, mu, Sigma, prob, n, p, G, eps,
     *                   s, w, lwork, z, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c E-step for prescribed-shape Gaussian mixtures

      implicit double precision (a-h,o-z)

c     integer            n, p, G
      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p),mu(p,G),Sigma(p,p,G),prob(G),
c    *                   s(p),w(lwork),z(n,G)
      double precision   x(n,*),mu(p,*),Sigma(p,p,*),prob(*),
     *                   s(*), w(*), z(n,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  mu      double  (input) (p,G) means for each group. On output, contains
c                   the normalized shape vector for each group.
c  Sigma   double  (input) (p,p,G) covariance matrices. On output, the first
c                   element in each matrix contains lambda for that group.
c  prob    double  (input) (G) probabilities (not needed if equal proportions)
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape.
c                   On output, minimum values of these over the groups.
c  s       double  (scratch) (p) 
c  w       double  (scratch) (lwork)
c  lwork   integer (input) max(4*p,5*p-4) minumum workspace for LAPACK SVD.
c  z       double  (output) (n,G) Conditional probabilities.
c  hood    double  (output) loglikelihood.

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      p1     = p + 1

      eps1   = max(eps(1),zero)
      eps2   = sqrt(max(eps(2),zero))

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
      EPSMIN = d1mach(3)

      vlamin = FLMAX
      rcmin  = FLMAX

      if (G .lt. 0) goto 500

      do k = 1, G

        call dpotrf( 'U', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*EPSMIN
          return
        end if

c zero out the lower triangle
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, Sigma(j,i,k),1)
          i = j
        end do

        call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, dum, 1, dum, 1,
     *                w, lwork, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, s, 1, smin, smax)

        rc    = smin/(one+smax)

        rcmin = min(rc,rcmin)

        if (rcmin .le. eps2) then
          eps(1) =  vlamin
          eps(2) =  rcmin*rcmin
          hood   = -FLMAX
          return
        end if

        smax = zero
        sum  = zero
        do j = 1, p
          temp = s(j)
          smax = max(smax,temp)
          sum  = sum + log(temp*temp)
        end do
        vlam   = exp(sum/dble(p))
        vlamin = min(vlamin,vlam)
        if (vlam .gt. eps1) then
          const = dble(p)*pi2log+sum
          slam = sqrt(vlam)
          call dscal( p, one/slam, s, 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w(p1), 1)
            call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
            call dgemv('N',p,p,one,Sigma(1,1,k),p,w(p1),1,zero,w,1)
            do j = 1, p
              w(j) = w(j) / s(j)
            end do
            temp   = ddot(p,w,1,w,1)/vlam
            z(i,k) = prob(k) * exp(-(const+temp)/two)
          end do
          call dscal( p, slam, s, 1)
        end if
        Sigma(1,1,k) = vlam
        call dscal( p, one/smax, s, 1)
        call dcopy( p, s, 1, mu(1,k), 1)
      end do

      if (vlamin .le. eps1) then
        eps(1) = vlamin
        eps(2) = rcmin*rcmin
        hood   = FLMAX
        return
      end if

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      eps(1) = vlamin
      eps(2) = rcmin*rcmin

      return

500   continue

      G      = -G

      do k = 1, G

        call dpotrf( 'U', p, Sigma(1,1,k), p, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*EPSMIN
          return
        end if

c zero out the lower triangle
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, Sigma(j,i,k),1)
          i = j
        end do

        call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, dum, 1, dum, 1,
     *                w, lwork, info)
        if (info .ne. 0) then
          eps(1) = FLMAX
          eps(2) = FLMAX
          hood   = dble(sign(1,info))*FLMIN
          return
        end if

        call drnge( p, s, 1, smin, smax)

        rc    = smin/(one+smax)

        rcmin = min(rc,rcmin)

        if (rcmin .le. eps2) then
          eps(1) =  vlamin
          eps(2) =  rcmin*rcmin
          hood   = -FLMAX
          return
        end if

        smax = zero
        sum  = zero
        do j = 1, p
          temp = s(j)
          smax = max(smax,temp)
          sum  = sum + log(temp*temp)
        end do
        vlam = exp(sum/dble(p))
        vlamin = min(vlamin,vlam)
        if (vlam .gt. eps1) then
          const = dble(p)*pi2log+sum
          slam = sqrt(vlam)
          call dscal( p, one/slam, s, 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w(p1), 1)
            call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
            call dgemv('N',p,p,one,Sigma(1,1,k),p,w(p1),1,zero,w,1)
            do j = 1, p
              w(j) = w(j) / s(j)
            end do
            temp   = ddot(p,w,1,w,1)/vlam
            z(i,k) = exp(-(const+temp)/two)
          end do
          call dscal( p, slam, s, 1)
        end if
        Sigma(1,1,k) = vlam
        call dscal( p, one/smax, s, 1)
        call dcopy( p, s, 1, mu(1,k), 1)
      end do

      if (vlamin .le. eps1) then
        eps(1) = vlamin
        eps(2) = rcmin*rcmin
        hood   = FLMAX
        return
      end if

      hood = - dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do

      eps(1) = vlamin
      eps(2) = rcmin*rcmin

      return
      end
      subroutine hceee ( x, n, p, ic, ng, ns, io, jo, v, s, u, r)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithm in clusters share a common
c variance (shape, volume, and orientation are the same for all clusters).

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns, io(*), jo(*)

c     double precision   x(n,p), v(p), s(p*p), u(p,p), r(p,p)
      double precision   x(n,*), v(*), s(*), u(*), r(*)
c------------------------------------------------------------------------------
c
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the determinant and trace of the
c                   sum of the sample cross product matrices. Columns 3 and 4
c                   contain the merge indices if p .ge. 4
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  io,jo   integer (output [p .le. 3]) If p .lt. 3, both io and jo must be of
c                   length ns and contain the indices of the merged pairs on
c                   output. If p .eq. 3, jo must be of length ns and contains
c                   an index of each merged on output pair. Otherwise io and
c                   jo are not used and can be of length 1.
c  v       double  (scratch/output) (p) On output, algorithm breakpoints;
c                   tells where the algorithm switches from using trace
c                   to trace + det, and from trace + det to det as criterion.
c  s       double  (scratch/output) (p,p) On output the first column contains
c                   the initial trace and determinant of the sum of sample
c                   cross product matrices.
c  u,r      double  (scratch) (p,p)

      integer                 lw, q

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision       d1mach

      double precision       FLMAX
      common /MCLMCH/        FLMAX
      save   /MCLMCH/

c------------------------------------------------------------------------------

      i1 = 0
      i2 = 0

      lw = p*p

c     FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
      EPSMIN = d1mach(3)
c     EPSMCH = d1mach(4)

c     call intpr('ic', -1, ic, n)
      
c form scaled column sums
      call dscal( n*p, one/sqrt(dble(n)), x, 1)

      si = one/sqrt(dble(p))
      sj = si / dble(n)
      call dcopy( p, zero, 0, v, 1)
      do k = 1, n
        call daxpy( p, sj, x(k,1), n, v, 1)
      end do

      trc0 = zero
      call dcopy( lw, zero, 0, r, 1)
      do k = 1, n
        call dcopy( p, v, 1, s, 1)
        call daxpy( p, (-si), x(k,1), n, s, 1)
        trc0 = trc0 + ddot( p, s, 1, s, 1)
        call mclrup( (k+1), p, s, r, p)
      end do
      
      det0 = detmc2( p, r)

      DELOG = log(trc0+EPSMIN)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c     call intpr( 'ic', -1, ic, n)

      trcw = zero
      call dcopy( lw, zero, 0, r, 1)

      q = 1
      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
c update trace and Cholesky factor as if a merge
          q     = q + 2
          ni    = ic(i)
          ri    = dble(ni)
          rij   = dble(ni+1)
          sj    = sqrt(one/rij)
          si    = sqrt(ri)*sj
          call dcopy( p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
          trcw  = trcw + ddot(p, v, 1, v, 1)
          call mclrup( q, p, v, r, p)
          ic(j) = 0
          ic(i) = ic(i) + 1
          call dscal( p, si, x(i,1), n)
          call daxpy( p, sj, x(j,1), n, x(i,1), n)
c         call dcopy( p, FLMAX, 0, x(j,1), n)
c update column sum in jth row
        else
          ic(j) = 1
        end if
      end do

c     call intpr('ic', -1, ic, n)

      trc1 = trcw

      if (q .lt. p) then
        detw = -FLMAX
      else
        detw = detmc2( p, r)
      end if        

      det1 = detw

      ls =  1

      lg = ng

      l1 =  0
      l2 =  0

 100  continue

      if (q .ge. p) then
c       if (.false.) 
c    *    call intpr('PART 2 --------------------------', -1, ls, 0)
        if (detw .lt. DELOG) then
          goto 200
        else
          goto 300
        end if
      end if

      dopt = FLMAX

      do j = 2, lg
        nj = ic(j)
        rj = dble(nj)
        do i = 1, (j-1)
          ni  = ic(i)
          ri  = dble(ni)
          nij = ni + nj
          rij = dble(nij)
          si  = sqrt(ri/rij)
          sj  = sqrt(rj/rij)
          call dcopy( p, x(i,1), n, s, 1)
          call dscal( p, sj, s, 1)
          call daxpy( p, (-si), x(j,1), n, s, 1)
          tij = trcw + ddot(p, s, 1, s, 1)
          zij = max(tij,EPSMIN)
          if (zij .le. dopt) then
            dopt = zij
            nopt = nij
            siop = si
            sjop = sj
            iopt = i
            jopt = j
            call dcopy( p, s, 1, v, 1)
          end if
        end do
      end do

      trcw = dopt

      if (ls .eq. ns) goto 900

      call dscal( p, siop, x(iopt,1), n)
      call daxpy( p, sjop, x(jopt,1), n, x(iopt,1), n)

      if (jopt .ne. lg) then
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        ic(jopt) = ic(lg)
      end if

      ic(iopt) = nopt

      x(lg,1)  = detw
      x(lg,2)  = trcw
      if (p .ge. 4) then
        x(lg,3) = dble(iopt)
        x(lg,4) = dble(jopt)     
      else if (p .eq. 3) then
        x(lg,3) = dble(iopt)
        jo(ls)  = jopt
      else 
        io(ls)  = iopt
        jo(ls)  = jopt
      end if

c update the Cholesky factor

      q  =  q + 1
      call mclrup( q, p, v, r, p)

      ls = ls + 1
      lg = lg - 1

      goto 100

 200  continue

      q  =  q + 1

c     call intpr('ic', -1, ic, n)

      dopt = FLMAX
      zopt = FLMAX

      do j = 2, lg
        nj = ic(j)        
        rj = dble(nj)
        do i = 1, (j-1)
          ni  = ic(i)
          ri  = dble(ni)
          nij = ni + nj
          rij = dble(nij)
          si  = sqrt(ri/rij)
          sj  = sqrt(rj/rij)
          call dcopy( p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
          tij = trcw +  ddot(p, v, 1, v, 1)
          call dcopy( lw, r, 1, u, 1)
          call mclrup( q, p, v, u, p)
          dij = detmc2( p, u)
          if (dij .le. dopt) then
            dopt = dij
            tdet = tij
            ndet = nij
            sidt = si
            sjdt = sj
            idet = i
            jdet = j
          end if
          if (tij .eq. zero) then
            zij = -FLMAX
          else 
            zij = max(tij,EPSMIN)
            if (dij .eq. (-FLMAX)) then
              zij = log(zij)
            else if (dij .le. zero) then
              zij = log(exp(dij) + zij)
            else 
              zij = log(one + zij*exp(-dij)) + dij
            end if
          end if
          if (zij .le. zopt) then
            zopt = zij
            dijo = dij
            tijo = tij
            nopt = nij
            siop = si
            sjop = sj
            iopt = i
            jopt = j
            call dcopy( lw, u, 1, s, 1)
          end if
        end do
      end do

      if (dopt .lt. DELOG) then
        if (l1 .eq. 0) l1 = ls
        trcw = tijo
        detw = dijo
        call dcopy( lw, s, 1, r, 1)
      else 
        l2 = ls
        trcw = tdet
        detw = dopt
        siop = sidt
        sjop = sjdt
        nopt = ndet
        iopt = idet
        jopt = jdet
        call dcopy( p, x(iopt,1), n, v, 1)
        call dscal( p, sjop, v, 1)
        call daxpy( p, (-siop), x(jopt,1), n, v, 1)
        call mclrup( q, p, v, r, p)
      end if

      if (ls .eq. ns) goto 900

      call dscal( p, siop, x(iopt,1), n)
      call daxpy( p, sjop, x(jopt,1), n, x(iopt,1), n)

      if (jopt .ne. lg) then
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        ic(jopt) = ic(lg)
      end if

      ic(iopt) = nopt

      x(lg,1)  = detw
      x(lg,2)  = trcw
      if (p .ge. 4) then
        x(lg,3) = dble(iopt)
        x(lg,4) = dble(jopt)     
      else if (p .eq. 3) then
        x(lg,3) = dble(iopt)
        jo(ls)  = jopt
      else 
        io(ls)  = iopt
        jo(ls)  = jopt
      end if

      ls = ls + 1
      lg = lg - 1

      if (detw .ge. DELOG) then
c       if (.false.)
c    *    call intpr('PART 3 --------------------------', -1, ls, 0)
        goto 300
      end if

      goto 200

 300  continue

      q = q + 1

      detw = FLMAX

      do j = 2, lg
        nj = ic(j)        
        rj = dble(nj)
        do i = 1, (j-1)
          ni  = ic(i)
          ri  = dble(ni)
          nij = ni + nj
          rij = dble(nij)
          si  = sqrt(ri/rij)
          sj  = sqrt(rj/rij)
          call dcopy( p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
          call dcopy( lw, r, 1, u, 1)
          call mclrup( q, p, v, u, p)
          dij = detmc2( p, u)
          if (dij .le. detw) then
            detw = dij
            nopt = nij
            siop = si
            sjop = sj
            iopt = i
            jopt = j
            call dcopy( lw, u, 1, s, 1)
          end if
        end do
      end do

c update the trace

      call dcopy( p, x(iopt,1), n, v, 1)
      call dscal( p, sjop, v, 1)
      call daxpy( p, (-siop), x(jopt,1), n, v, 1)

      trcw = trcw + ddot( p, v, 1, v, 1)

      if (ls .eq. ns) goto 900

      call dcopy( lw, s, 1, r, 1)

      call dscal( p, siop, x(iopt,1), n)
      call daxpy( p, sjop, x(jopt,1), n, x(iopt,1), n)

      if (jopt .ne. lg) then
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        ic(jopt) = ic(lg)
      end if

      ic(iopt) = nopt

      x(lg,1)  = detw
      x(lg,2)  = trcw
      if (p .ge. 4) then
        x(lg,3) = dble(iopt)
        x(lg,4) = dble(jopt)     
      else if (p .eq. 3) then
        x(lg,3) = dble(iopt)
        jo(ls)  = jopt
      else 
        io(ls)  = iopt
        jo(ls)  = jopt
      end if

      ls = ls + 1
      lg = lg - 1

      goto 300

 900  continue

      x(lg,1) = detw
      x(lg,2) = trcw
      if (p .ge. 4) then
        if (iopt .lt. jopt) then
          x(lg,3) = dble(iopt)
          x(lg,4) = dble(jopt)     
        else
          x(lg,3) = dble(jopt)
          x(lg,4) = dble(iopt)     
        end if
      else if (p .eq. 3) then
        if (iopt .lt. jopt) then
          x(lg,3) = dble(iopt)
          jo(ls)  = jopt
        else
          x(lg,3) = dble(jopt)
          jo(ls)  = iopt
        end if
      else 
        if (iopt .lt. jopt) then
          io(ls)  = iopt
          jo(ls)  = jopt
        else
          io(ls)  = jopt
          jo(ls)  = iopt
        end if
      end if

c decode

      do k = 1, ng
       ic(k) = k
      end do

      m = ng + 1
      if (p .ge. 4) then
        l = m
        do k = 1, ns
          l      = l - 1
          i      = int(x(l,3))
          ici    = ic(i)
          j      = int(x(l,4))
          icj    = ic(j)
          if (ici .gt. icj) ic(i) = icj
          ic(j)  = ic(m - k)
          if (ici .lt. icj) then
            x(l,3) = dble(ici)
            x(l,4) = dble(icj)
          else
            x(l,3) = dble(icj)
            x(l,4) = dble(ici)
          end if
        end do
      else if (p .eq. 3) then
        l = m
        do k = 1, ns
          l      = l - 1
          i      = int(x(l,3))
          ici    = ic(i)
          j      = jo(k)
          icj    = ic(j)
          if (ici .gt. icj) ic(i) = icj
          ic(j)  = ic(m - k)
          if (ici .lt. icj) then
            x(l,3) = dble(ici)
            jo(k)  = icj
          else
            x(l,3) = dble(icj)
            jo(k)  = ici
          end if
        end do
      else
        do k = 1, ns
          i      = io(k)
          ici    = ic(i)
          j      = jo(k)
          icj    = ic(j)
          if (ici .gt. icj) ic(i) = icj
          ic(j)  = ic(m - k)
          if (ici .lt. icj) then
            io(k)  = ici
            jo(k)  = icj
          else
            io(k)  = icj
            jo(k)  = ici
          end if
        end do
      end if

      l = 2
      m = min(p,4)
      do k = ng, lg, -1
        if (k .le. l) goto 950
        call dswap( m, x(k,1), n, x(l,1), n)
        l  = l + 1
      end do
   
 950  continue

      x(1,1) = det1
      x(1,2) = trc1

      v(1) = dble(l1)
      v(2) = dble(l2)

      s(1) = det0
      s(2) = trc0

      return
      end
      subroutine hcefv ( x, n, p, ic, ng, ns, lwork, s,
     *                   v, w, t, u, a, r, nd, d)
c    *                   v, w, t, u, a, r, uu, vv, nd, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithm in which shape is fixed in advance,
c cluster volumes are equal but unknown, and orientation is allowed to vary 
c among clusters.

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns, nd

c     double precision   x(n,p), s(p), v(p), w(*), t(p,p), u(p,p), a(p,p)
c     double precision   r(p,p), d(*)
      double precision   x(n,*),s(*),v(*),w(*),t(p,*),u(p,*),a(p,*)
      double precision   r(p,*), d(*)
c     double precision   uu(p*p), vv(p*p)
c     double precision   uu(*), vv(*)
c------------------------------------------------------------------------------

c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  s       double  (input) (p) Shape matrix. It's square is the A matrix of 
c                   Banfield and Raftery.
c  v       double  (scratch) (p) 
c  t,u,a,r double  (scratch) (p*p)
c  w       double  (scratch) (max(4*p,5*p-4))
c  nd      integer (input/output) On input, the length of d. On output, the
c                   indicator from the SVD computation (abnormal if nonzero).
c  d       double  (scratch/output) max(p*p+n,((ng*(ng-1))/2,3*ns). On output
c                   the first ns elements are proportional to the change in
c                   loglikelihood associated with each merge.

      integer                 psq, pm1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        rthalf
      parameter              (rthalf = .7071067811865476d0)

      double precision        d1mach

c------------------------------------------------------------------------------

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd
     
c     FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
      EPSMIN = d1mach(3)
c     EPSMAX = d1mach(4)

      psq = p*p
      pm1 = p-1

c     lwork = max(4*p,5*p-4)
      lwopt = 0
       
      do i = 1, p
        s(i) = one/s(i)
      end do

c     call dblepr( '1/sqrt(shape)', -1, s, 2) 

c form scaled column sums
c     si = one / (sqrt(dble(n))*sqrt(dble(p)))
c     sj = si / dble(n)
c     call dcopy( p, zero, 0, v, 1)
c     do k = 1, n
c       call daxpy( p, sj, x(k,1), n, v, 1)
c     end do

c form sum of squares
c     ss = zero
c     do k = 1, p
c       call dcopy( n, x(1,k), 1, d, 1)
c       call dscal( n, si, d, 1)
c       call daxpy( n, (-one), v(k), 0, d, 1)
c       ss = ss + ddot( n, d, 1, d, 1)
c     end do

      if (ng .eq. 1) then
c       s(1) = ss
        nd   = 0
        return
      end if

c     call intpr( 'ic', -1, ic, n)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     call intpr( 'ic', -1, ic, n)

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, v, 1)
          call dscal( p, si, v, 1)
          call daxpy( p, (-sj), x(k,1), n, v, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, v, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c copy triangular factor into the rows of x
          j  = k
          m  = p
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij) = n+l
          if (l .gt. 2) then
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do i = 1, min(l-1,p)
              call dcopy( m, r(i,i), p, t(i,i), p)
              m = m - 1
            end do
            m = min(p,l-1)
            call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                    dummy, 1, dummy, 1, w, lwork, info)
c           call dgesvd( 'A', 'A', m, p, t, p, v, 
c    *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do i = 1, m
              v(i) = v(i)*s(i)
            end do
            trmij          = ddot( m, v, 1, v, 1)
            x(ic(ic(k)),1) = trmij
          else
            temp  = dnrm2( p, r, p)*s(1)
            trmij = temp*temp
          end if
          d(ll+k) = trmij
        else
          ic(k)   = 1
          d(ll+k) = zero
        end if
      end do

c     call intpr( 'ic', -1, ic, n)
c     call dblepr( 'term', -1, d(ll+1), ng)

c compute change in likelihood and determine minimum

      dopt = FLMAX

      ij = 0
      do j = 2, ng
        icj = ic(j)
        nj  = 1
        if (icj .eq. 1) then
          termj = zero
          do i = 1, (j-1)
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij   = 2
              rij   = two
              si    = rthalf
              sj    = rthalf
              sij   = rthalf
              call dcopy( p, x(i,1), n, v, 1)
              call daxpy( p, (-one), x(j,1), n, v, 1)
              call dscal( p, rthalf, v, 1)
              call dcopy( p, v, 1, u, p)  
              temp  = dnrm2(p,v,1)*s(1)
              trmij = temp*temp
              termi = zero
            else
              m  = p
              l  = ici
 110          continue
              call dcopy( m, x(l,ni), n, u(ni,ni), p)
              ni  = ni + 1
              m   = m - 1
              l   = ic(l)
              if (l .le. n) goto 110
              ni  = l - n
              termi = d(ll+i)
              ri  = dble(ni)
              nij = ni + 1
              rij = dble(nij)
              qij = one/rij
              qi  = ri*qij
              si  = sqrt(qi)
              sj  = sqrt(qij)
              sij = sj
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              call mclrup( nij, p, v, u, p)
              call dcopy( psq, zero, 0, t, 1)
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, t(k,k), p)
                m = m - 1
              end do
              m = min(p,nij-1)
              call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                      dummy, 1, dummy, 1, w, lwork, info)
c             call dgesvd( 'A', 'A', m, p, t, p, v, 
c    *                      uu, p, vv, p, w, lwork, info)
              if (info .lt. 0) then
                call intpr( 'SVD fails', -1, info, 1)
                nd = info
                return
              end if
              if (info .gt. 0) then
                call intpr( 'SVD does not converge', -1, info, 1)
                nd = info
                return
              end if
              lwopt = max(lwopt,int(w(1)))
              do k = 1, m
                v(k) = v(k)*s(k)
              end do
              trmij = ddot( m, v, 1, v, 1)
            end if
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij
            if (dij .le. dopt) then
              dopt = dij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, a(nj,nj), p)
          m  = m  - 1
          nj = nj + 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj = l - n
          termj = d(ll+j)
          rj = dble(nj)
          do i = 1, (j-1)
            m = p
            do k = 1, min(nj-1,p)
              call dcopy( m, a(k,k), p, u(k,k), p)
              m = m - 1
            end do
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij = nj + 1
              rij = dble(nij)
              qij = one/rij
              qi  = qij
              qj  = rj*qij
              si  = sqrt(qi)
              sj  = sqrt(qj)
              sij = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              termi = zero
            else
              m  = p
              l  = ici
              k  = nj + 1
 130          continue
              call dcopy( m, x(l,ni), n, v, 1)
              call mclrup( k, m, v, u(ni,ni), p)
              ni    = ni + 1
              m     = m - 1
              l     = ic(l)
              if (l .le. n) goto 130
              ni    = l - n
              termi = d(ll+i)
              ri    = dble(ni)
              nij   = ni + nj
              rij   = dble(nij)
              qij   = one/rij
              qi    = ri*qij
              qj    = rj*qij
              si    = sqrt(qi)
              sj    = sqrt(qj)
              sij   = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
            end if
            call mclrup( nij, p, v, u, p)
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do k = 1, min(nij-1,p)
              call dcopy( m, u(k,k), p, t(k,k), p)
              m = m - 1
            end do
            m = min(p,nij-1)
            call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                    dummy, 1, dummy, 1, w, lwork, info)
c           call dgesvd( 'A', 'A', m, p, t, p, v, 
c    *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do k = 1, m
              v(k) = v(k)*s(k)
            end do
            trmij = ddot( m, v, 1, v, 1)
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij
            if (dij .le. dopt) then
              dopt = dij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        end if
      end do

c     if (.false.) then
c       i  = 1
c       ij = 1
c       do j = 2, ng
c         call dblepr( 'dij', -1, d(ij), i)
c         ij = ij + i
c         i  = j
c       end do
c     end if
 
      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1,1) = dble(iopt)
          x(1,2) = dble(jopt)
        else
          x(1,1) = dble(jopt)
          x(1,2) = dble(iopt)
        end if
        d(1)   = dopt
        nd     = 0
        lwork  = lwopt
        return
      end if

      ls  = 1

 200  continue

c     call dblepr( 'col 1', -1, x(1,1), n)

      call dcopy( p, x(iopt,1), n, v, 1)
      call dscal( p, siop, v, 1)
      call daxpy( p, sjop, x(jopt,1), n, v, 1)

      if (jopt .ne. lg) then
        call mclcpy( jopt, lg, d)
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        m        = ic(jopt)
        ic(jopt) = ic(lg)
        ic(lg)   = m
      end if

      call dcopy( p, r(1,1), p, x(lg,1), n)

      if (niop .eq. 1) then

        if (njop .eq. 1) then
          ic(lg)  = n+2
        else
          l   = ic(lg)
          m   = pm1
          nij = 2
 210      continue
            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

        end if

      else

        l   = ic(iopt)
        m   = pm1
        nij = 2
 220    continue
          call dcopy( m, r(nij,nij), p, x(l,nij), n)
          nij = nij + 1
          m   = m   - 1
          k   = l
          l   = ic(l)
          if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

          if (nij .le. p .and. njop .ne. 1) then
            l     = ic(lg)
            ic(k) = l
 230        continue
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij   = nij + 1
              m     = m   - 1
              k     = l
              l     = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230
          end if

        ic(lg) = ic(iopt)
        ic(k)  = n + nopt

      end if

      ic(iopt) = lg

      if (nopt .gt. 2) then
        m      = ic(lg)
        x(m,1) = tmop
      endif   

c     if (ls. eq. 2) return

      call dcopy( p, v, 1, x(iopt,1), n)

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)
      lo     = lo - 1

      lg = lg - 1
      ld = ld - lg

c     call intpr( 'ic', -1, ic, n)

      iold  =  iopt

      dopt  = FLMAX

      ni    = nopt
      ri    = dble(ni)
      termi = tmop
      traci = trop

      ij = ((iold-1)*(iold-2))/2
      if (iold .gt. 1) then
        do j = 1, (iold-1)
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            termj = zero
          else
            m = p
            l = icj
            k = ni + 1
 310        continue
              call dcopy( m, x(l,nj), n, v, 1)
              call mclrup( k, m, v, u(nj,nj), p)
              nj = nj + 1
              m  = m - 1
              l  = ic(l)
              if (l .le. n) goto 310
            nj  = l - n            
            rj   = dble(nj)
            if (nj .gt. 2) then
              termj = x(ic(icj),1)
            else
              temp  = dnrm2(p,x(icj,1),n)*s(1)
              termj = temp*temp
            end if
            nij  = ni + nj
            rij  = dble(nij)
            qij  = one/rij
            qi   = ri*qij
            qj   = rj*qij
            si   = sqrt(qi)
            sj   = sqrt(qj)
            sij  = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
          end if
          call mclrup( nij, p, v, u, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, u(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(p,nij-1)
          call dgesvd( 'N', 'N', m, p, t, p, v,
     *                  dummy, 1, dummy, 1, w, lwork, info)
c         call dgesvd( 'A', 'A', m, p, t, p, v,
c    *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          trmij = ddot( m, v, 1, v, 1)
          dij   = trmij - (termi + termj)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = nj
            njop = ni
            siop = sj
            sjop = si
            iopt = j
            jopt = iold
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, a(k,k), p)
              m = m - 1
            end do
          end if
        end do
      end if
    
      if (iold .lt. lg) then
        i  = iold
        ij = ij + i
        do j = (iold+1), lg
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            termj = zero
          else
            m = p
            l = icj
            k = ni + 1
 410        continue
            call dcopy( m, x(l,nj), n, v, 1)
            call mclrup( k, m, v, u(nj,nj), p)
            nj = nj + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 410
            nj  = l - n
            rj  = dble(nj)
            if (nj .gt. 2) then
              termj = x(ic(icj),1)
            else
              temp  = dnrm2(p,x(icj,1),n)*s(1)
              termj = temp*temp
            end if
            nij = ni + nj
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            sij = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
          end if
          call mclrup( nij, p, v, u, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, u(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(p,nij-1)
          call dgesvd( 'N', 'N', m, p, t, p, v,
     *                  dummy, 1, dummy, 1, w, lwork, info)
c         call dgesvd( 'A', 'A', m, p, t, p, v,
c    *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          trmij = ddot( m, v, 1, v, 1)
          dij   = trmij - (termi + termj)
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = ni
            njop = nj
            siop = si
            sjop = sj
            iopt = iold
            jopt = j
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, a(k,k), p)
              m = m - 1
            end do
          end if
          ij = ij + i
          i  = j
        end do
      end if

c update d and find max

      jopt = 2
      iopt = 1

      dopt = d(1)

      if (lg .eq. 2) goto 900
        
      ij   = 1
      do i = 2, ld
        qi = d(i)
        if (qi .le. dopt) then
          ij   = i
          dopt = qi
        end if
      end do

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      if (iopt .ne. iold .and. jopt .ne. iold) then

        i   = iopt
        j   = jopt

        nj  = 1
        icj = ic(j)
        ni  = 1
        ici = ic(i)
        if (icj .eq. 1) then
          termj = zero
          if (ici .eq. 1) then
            nij = 2
            rij = two
            si  = rthalf
            sj  = rthalf
            call dcopy(p, x(i,1), n, v, 1)
            call daxpy( p, (-one), x(j,1), n, v, 1)
            call dscal( p, rthalf, v, 1)
            call dcopy( p, v, 1, r, p)
            termi = zero
          else
            m = p
            l = ici
 610        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 610
            ni  = l - n
            ri  = dble(ni)
            if (ni .gt. 2) then
              termi = x(ic(ici),1)
            else
              temp  = dnrm2(p,x(ici,1),n)*s(1)
              termi = temp*temp
            end if
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            call mclrup( nij, p, v, r, p)
          end if
        else
          m = p
          l = icj
 620      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 620
          nj  = l - n
          rj  = dble(nj)
          if (nj .gt. 2) then
            termj = x(ic(icj),1)
          else
            temp  = dnrm2(p,x(icj,1),n)*s(1)
            termj = temp*temp
          end if
          if (ici .eq. 1) then
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qj  = rj*qij
            si  = sqrt(qij)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            termi = zero
          else 
            m = p
            l = ici
            k = nj + 1
 630        continue
            call dcopy( m, x(l,ni), n, v, 1)
            call mclrup( k, m, v, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 630
            ni  = l - n
            ri  = dble(ni)
            if (ni .gt. 2) then
              termi = x(ic(ici),1)
            else
              temp  = dnrm2(p,x(ici,1),n)*s(1)
              termi = temp*temp
            end if
            nij = ni + nj
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          end if
          call mclrup( nij, p, v, r, p)
        end if

        trop = trcij
        tmop = dopt + (termi + termj)
        nopt = nij
        niop = ni
        njop = nj
        siop = si
        sjop = sj

      else
        m = p
        do k = 1, min(nopt-1,p)
          call dcopy(m, a(k,k), p, r(k,k), p)
          m = m - 1
        end do
      end if

      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 200

 900  continue

      d(lo) = dopt
      lo    = lo - 1
      d(lo) = dble(iopt)
      lo    = lo - 1
      d(lo) = dble(jopt)

c     s(1) = ss

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 1
      ld          = lo
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1,1) = si 
        x(1,2) = sj
      else
        x(1,1) = sj
        x(1,2) = si
      end if

      lg = ng + 1
      do k = 2, ns
        lo     = lo - 1
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 1
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then
          x(k,1) = dble(ici)
          x(k,2) = dble(icj)
        else
          x(k,1) = dble(icj)
          x(k,2) = dble(ici)
        end if
      end do

      ld = nd
      lo = 1
      do k = 1, ns
        si    = d(lo)
        d(lo) = d(ld)
        d(ld) = si
        ld    = ld - 1
        lo    = lo + 1
      end do

      nd    = 0

      lwork = lwopt

      return
      end
      subroutine hcei  ( x, n, p, ic, ng, ns, v, nd, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithm in which shape and orientation
c are fixed in advance, while cluster volumes are also equal but unknown.

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns, nd

c     double precision   x(n,p), v(p), d(l*(l-1)/2)
      double precision   x(n,*), v(*), d(*)
c------------------------------------------------------------------------------
c
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  v       double  (scratch)
c  nd      integer (input) The length of d.
c  d       double  (scratch/output) max(((ng-1)*(ng-2))/2,3*ns). On output
c                   the first ns elements are proportional to the change in
c                   loglikelihood associated with each merge.

      double precision        one
      parameter              (one = 1.d0)

c------------------------------------------------------------------------------

      FLMAX  = d1mach(2)

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd

c     call intpr( 'ic', -1, ic, n)
c     call intpr( 'no. of groups', -1, lg, 1)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. lg)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c     call intpr( 'ic', -1, ic, n)

      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
          ic(j) = 0
          ni    = ic(i)
          nij   = ni + 1
          ic(i) = nij
          ri    = dble(ni)
          rij   = dble(nij)
          sj    = sqrt(one/rij)
          si    = sqrt(ri)*sj
c update column sum in kth row
          call dscal( p, si, x(i,1), n)
          call daxpy( p, sj, x(j,1), n, x(i,1), n)
        else 
          ic(j) = 1
        end if
      end do

c     call intpr( 'ic', -1, ic, n)

      dopt = FLMAX

      ij   = 0
      do j = 2, lg
        nj = ic(j)
        rj = dble(nj)
        do i = 1, (j-1)
          ni    = ic(i)
          ri    = dble(ni)
          nij   = ni + nj
          rij   = dble(nij)
          si    = sqrt(ri/rij)
          sj    = sqrt(rj/rij)
          call dcopy( p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
          dij   = ddot(p, v, 1, v, 1)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            iopt = i
            jopt = j
          end if
        end do
      end do

c     if (.false.) then
c       i  = 1
c       ij = 1
c       do j = 2, ng
c         call dblepr( 'dij', -1, d(ij), i)
c         ij = ij + i
c         i  = j
c       end do
c     end if

      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1,1) = iopt
          x(1,2) = jopt
        else
          x(1,1) = jopt
          x(1,2) = iopt
        end if
        d(1)   = dopt
        return
      end if

      ls  = 1

 100  continue

      ni       = ic(iopt)
      nj       = ic(jopt)
      nij      = ni + nj 

      ic(iopt) =  nij
      ic(jopt) = -iopt

      if (jopt .ne. lg) then
        call mclswp( jopt, lg, d)
        m        = ic(jopt)
        ic(jopt) = ic(lg)
        ic(lg)   = m
      end if

      si   = dble(ni)
      sj   = dble(nj)
      sij  = dble(nij)

      dold = dopt

      iold = iopt
      jold = jopt

      iopt = -1
      jopt = -1

      dopt = FLMAX

      lg = lg - 1

      ld = ld - lg

      ii = (iold*(iold-1))/2

      if (iold .gt. 1) then
        ik = ii - iold + 1
        do j = 1, (iold - 1)
          nj    = ic(j)        
          rj    = dble(nj)
          ik    = ik + 1
          jk    = ld + j
          dij   = (rj+si)*d(ik)+(rj+sj)*d(jk)
          dij   = (dij-rj*dold)/(rj+sij)
          d(ik) = dij
        end do
      end if

      if (iold .lt. lg) then
        ik = ii + iold
        i  = iold
        do j = (iold + 1), lg
          nj    = ic(j)        
          rj    = dble(nj)
          jk    = ld + j
          dij   = (rj+si)*d(ik)+(rj+sj)*d(jk)
          dij   = (dij-rj*dold)/(rj+sij)
          d(ik) = dij
          ik    = ik + i
          i     = j
        end do
      end if

      d(lo) = dold
      lo    = lo - 1
      d(lo) = dble(iold)
      lo    = lo - 1
      d(lo) = dble(jold)
      lo    = lo - 1

c update d and find max

      jopt = 2
      iopt = 1

      dopt = d(1)

      if (lg .eq. 2) goto 900

      ij   = 1
      do i = 2, ld
        qi = d(i)
        if (qi .le. dopt) then
          ij   = i
          dopt = qi
        end if
      end do

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 100

 900  continue

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 1
      ld          = lo
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1,1) = si 
        x(1,2) = sj
      else
        x(1,1) = sj
        x(1,2) = si
      end if

      lg = ng + 1
      do k = 2, ns
        lo     = lo - 1
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 1
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then
          x(k,1) = dble(ici)
          x(k,2) = dble(icj)
        else
          x(k,1) = dble(icj)
          x(k,2) = dble(ici)
        end if
      end do

      ld = nd
      lo = 1
      do k = 1, ns
        si    = d(lo)
        d(lo) = d(ld)
        d(ld) = si
        ld    = ld - 1
        lo    = lo + 1
      end do

      return
      end

      subroutine mclswp( i, n, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Copies row n into row i and puts FLMAX in row n.
c Note : i > 1, i < n, n > 2 required.

      implicit double precision (a-h,o-z)

      double precision d(*)

      FLMAX = d1mach(2)

      i1 = i - 1
      ii = (i1*(i1-1))/2 + 1

      n1 = n - 1
      nn = (n1*(n1-1))/2 + 1

c     if (i .gt. 1) then
        call dswap( i1, d(nn), 1, d(ii), 1)
c       call dcopy( i1, FLMAX, 0, d(nn), 1)
        ii = ii + i1 + i1
        nn = nn + i
c     end if

      if (n1 .eq. i) return

      k = i

 100  continue

        temp  = d(ii)
        d(ii) = d(nn)
        d(nn) = temp

c       d(nn) = FLMAX

        ii = ii + k
 
        nn = nn + 1
        
        k  = k  + 1

        if (k .lt. n1) goto 100

c     d(nn) = FLMAX

      return
      end 
















      subroutine hcvfv ( x, n, p, ic, ng, ns, lwork, ALPHA, s, 
     *                   v, w, t, u, a, r, nd, d)
c    *                   v, w, t, u, a, r, uu, vv, nd, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithm in which shape is fixed in advance
c while volume and orientation are allowed to vary between clusters.

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns, nd

      double precision   ALPHA

c     double precision   x(n,p), s(p), v(p), w(*), t(p,p), u(p,p), a(p,p)
c     double precision   r(p,p), d(ng*(ng-1)/2)
      double precision   x(n,*),s(*),v(*),w(*),t(p,*),u(p,*),a(p,*)
      double precision   r(p,*), d(*)
c     double precision   uu(p*p), vv(p*p)
c     double precision   uu(*), vv(*)
c------------------------------------------------------------------------------

c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  ALPHA   double  (input) Additive quantity used to resolve degeneracies. 
c  s       double  (input) (p) Shape matrix. It's square is the A matrix of
c                   Banfield and Raftery.
c  v       double  (scratch) (p) 
c  t,u,a,r double  (scratch) (p*p)
c  w       double  (scratch) (max(4*p,5*p-4))
c  nd      integer (input/output) On input, the length of d. On output, the
c                   indicator from the SVD computation (abnormal if nonzero).
c  d       double  (scratch/output) max(p*p+n,((ng*(ng-1))/2,3*ns). On output
c                   the first ns elements are proportional to the change in
c                   loglikelihood associated with each merge.

      integer                 psq, pm1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        rthalf
      parameter              (rthalf = .7071067811865476d0)

      double precision        d1mach

c------------------------------------------------------------------------------

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd
     
c     FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
c     EPSMIN = d1mach(3)
      EPSMAX = d1mach(4)

      psq = p*p
      pm1 = p-1

c     lwork = max(4*p,5*p-4)

      lwopt = 0
       
      do i = 1, p
        s(i) = one/s(i)
      end do

c     call dblepr( '1/sqrt(shape)', -1, s, 2) 

      if (ng .eq. 1) then
        nd  = 0
        return
      end if

      ALPHA  = max(ALPHA,EPSMAX)

      ALFLOG = log(ALPHA)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     call intpr( 'ic', -1, ic, n)

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, v, 1)
          call dscal( p, si, v, 1)
          call daxpy( p, (-sj), x(k,1), n, v, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, v, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c copy triangular factor into the rows of x
          j  = k
          m  = p
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij) = n+l
          if (l .gt. 2) then
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do i = 1, min(l-1,p)
              call dcopy( m, r(i,i), p, t(i,i), p)
              m = m - 1
            end do
            m = min(p,l-1)
            call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                    dummy, 1, dummy, 1, w, lwork, info)
c           call dgesvd( 'A', 'A', m, p, t, p, v, 
c    *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do i = 1, m
              v(i) = v(i)*s(i)
            end do
            temp   = ddot( m, v, 1, v, 1)
c           trmij  = dble(l)*log(temp/dble(l) + ALPHA)
            trmij  = dble(l)*log((temp+ALPHA)/dble(l))
            m      = ic(ic(k))
            x(m,1) = trmij
          else
            temp   = dnrm2( p, r, p)*s(1)
            temp   = temp*temp
c           trmij  = dble(l)*log(temp/dble(l) + ALPHA)
            trmij  = dble(l)*log((temp+ALPHA)/dble(l))
          end if
          d(ll+k) = trmij
        else
          ic(k)   = 1
          d(ll+k) = ALFLOG
        end if
      end do

c     call intpr( 'ic', -1, ic, n)
c     call dblepr( 'term', -1, d(ll+1), ng)

c compute change in likelihood and determine minimum

      dopt = FLMAX

      ij = 0
      do j = 2, ng
        icj = ic(j)
        nj  = 1
        if (icj .eq. 1) then
          termj = ALFLOG
          do i = 1, (j-1)
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij   = 2
              rij   = two
              si    = rthalf
              sj    = rthalf
              sij   = rthalf
              call dcopy( p, x(i,1), n, v, 1)
              call daxpy( p, (-one), x(j,1), n, v, 1)
              call dscal( p, rthalf, v, 1)
              call dcopy( p, v, 1, u, p)  
              temp  = dnrm2(p,v,1)*s(1)
              temp  = temp*temp
c             trmij = rij*log(temp/rij + ALPHA)
              trmij = rij*log((temp + ALPHA)/rij)
              termi = ALFLOG
            else
              m  = p
              l  = ici
 110          continue
              call dcopy( m, x(l,ni), n, u(ni,ni), p)
              ni  = ni + 1
              m   = m - 1
              l   = ic(l)
              if (l .le. n) goto 110
              ni  = l - n
              termi = d(ll+i)
              ri  = dble(ni)
              nij = ni + 1
              rij = dble(nij)
              qij = one/rij
              qi  = ri*qij
              si  = sqrt(qi)
              sj  = sqrt(qij)
              sij = sj
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              call mclrup( nij, p, v, u, p)
              call dcopy( psq, zero, 0, t, 1)
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, t(k,k), p)
                m = m - 1
              end do
              m = min(p,nij-1)
              call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                      dummy, 1, dummy, 1, w, lwork, info)
c             call dgesvd( 'A', 'A', m, p, t, p, v, 
c    *                      uu, p, vv, p, w, lwork, info)
              if (info .lt. 0) then
                call intpr( 'SVD fails', -1, info, 1)
                nd = info
                return
              end if
              if (info .gt. 0) then
                call intpr( 'SVD does not converge', -1, info, 1)
                nd = info
                return
              end if
              lwopt = max(lwopt,int(w(1)))
              do k = 1, m
                v(k) = v(k)*s(k)
              end do
              temp  = ddot( m, v, 1, v, 1)
c             trmij = rij*log(temp/rij + ALPHA)
              trmij = rij*log((temp + ALPHA)/rij)
            end if
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij
            if (dij .le. dopt) then
              dopt = dij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, a(nj,nj), p)
          m  = m  - 1
          nj = nj + 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj = l - n
          termj = d(ll+j)
          rj = dble(nj)
          do i = 1, (j-1)
            m = p
            do k = 1, min(nj-1,p)
              call dcopy( m, a(k,k), p, u(k,k), p)
              m = m - 1
            end do
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij = nj + 1
              rij = dble(nij)
              qij = one/rij
              qi  = qij
              qj  = rj*qij
              si  = sqrt(qi)
              sj  = sqrt(qj)
              sij = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              termi = ALFLOG
            else
              m  = p
              l  = ici
              k  = nj + 1
 130          continue
              call dcopy( m, x(l,ni), n, v, 1)
              call mclrup( k, m, v, u(ni,ni), p)
              ni    = ni + 1
              m     = m - 1
              l     = ic(l)
              if (l .le. n) goto 130
              ni    = l - n
              termi = d(ll+i)
              ri    = dble(ni)
              nij   = ni + nj
              rij   = dble(nij)
              qij   = one/rij
              qi    = ri*qij
              qj    = rj*qij
              si    = sqrt(qi)
              sj    = sqrt(qj)
              sij   = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
            end if
            call mclrup( nij, p, v, u, p)
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do k = 1, min(nij-1,p)
              call dcopy( m, u(k,k), p, t(k,k), p)
              m = m - 1
            end do
            m = min(p,nij-1)
            call dgesvd( 'N', 'N', m, p, t, p, v, 
     *                    dummy, 1, dummy, 1, w, lwork, info)
c           call dgesvd('A', 'A', m, p, t, p, v, 
c    *                   uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do k = 1, m
              v(k) = v(k)*s(k)
            end do
            temp  = ddot( m, v, 1, v, 1)
c           trmij = rij*log(temp/rij + ALPHA)
            trmij = rij*log((temp + ALPHA)/rij)
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij
            if (dij .le. dopt) then
              dopt = dij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        end if
      end do

c     if (.false.) then
c       i  = 1
c       ij = 1
c       do j = 2, ng
c         call dblepr( 'dij', -1, d(ij), i)
c         ij = ij + i
c         i  = j
c       end do
c     end if
 
      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1,1) = iopt
          x(1,2) = jopt
        else
          x(1,1) = jopt
          x(1,2) = iopt
        end if
        d(1)   = dopt
        nd     = 0
        lwork  = lwopt
        return
      end if

      ls  = 1

 200  continue

      call dcopy( p, x(iopt,1), n, v, 1)
      call dscal( p, siop, v, 1)
      call daxpy( p, sjop, x(jopt,1), n, v, 1)

      if (jopt .ne. lg) then
        call mclcpy( jopt, lg, d)
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        m        = ic(jopt)
        ic(jopt) = ic(lg)
        ic(lg)   = m
      end if

      call dcopy( p, r(1,1), p, x(lg,1), n)

      if (niop .eq. 1) then

        if (njop .eq. 1) then
          ic(lg)  = n+2
        else
          l   = ic(lg)
          m   = pm1
          nij = 2
 210      continue
            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

        end if

      else

        l   = ic(iopt)
        m   = pm1
        nij = 2
 220    continue
          call dcopy( m, r(nij,nij), p, x(l,nij), n)
          nij = nij + 1
          m   = m   - 1
          k   = l
          l   = ic(l)
          if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

          if (nij .le. p .and. njop .ne. 1) then
            l     = ic(lg)
            ic(k) = l
 230        continue
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij   = nij + 1
              m     = m   - 1
              k     = l
              l     = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230
          end if

        ic(lg) = ic(iopt)
        ic(k)  = n + nopt

      end if

      ic(iopt) = lg

      if (nopt .gt. 2) then
        m      = ic(lg)
        x(m,1) = tmop
      endif   

c     if (ls. eq. 2) return

      call dcopy( p, v, 1, x(iopt,1), n)

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)
      lo     = lo - 1

      lg = lg - 1
      ld = ld - lg

c     call intpr( 'ic', -1, ic, n)

      iold  =  iopt

      dopt  =  FLMAX

      ni    = nopt
      ri    = dble(ni)
      termi = tmop
      traci = trop

      ij = ((iold-1)*(iold-2))/2
      if (iold .gt. 1) then
        do j = 1, (iold-1)
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            termj = ALFLOG
          else
            m = p
            l = icj
            k = ni + 1
 310        continue
              call dcopy( m, x(l,nj), n, v, 1)
              call mclrup( k, m, v, u(nj,nj), p)
              nj = nj + 1
              m  = m - 1
              l  = ic(l)
              if (l .le. n) goto 310
            nj  = l - n            
            rj   = dble(nj)
            if (nj .gt. 2) then
              termj = x(ic(icj),1)
            else
              temp  = dnrm2(p,x(icj,1),n)*s(1)
              temp  = temp*temp
c             termj = rj*log(temp/rj + ALPHA)
              termj = rj*log((temp + ALPHA)/rj)
            end if
            nij  = ni + nj
            rij  = dble(nij)
            qij  = one/rij
            qi   = ri*qij
            qj   = rj*qij
            si   = sqrt(qi)
            sj   = sqrt(qj)
            sij  = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
          end if
          call mclrup( nij, p, v, u, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, u(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(p,nij-1)
          call dgesvd( 'N', 'N', m, p, t, p, v,
     *                  dummy, 1, dummy, 1, w, lwork, info)
c         call dgesvd( 'A', 'A', m, p, t, p, v,
c    *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          temp  = ddot( m, v, 1, v, 1)
c         trmij = rij*log(temp/rij + ALPHA)
          trmij = rij*log((temp + ALPHA)/rij)
          dij   = trmij - (termi + termj)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = nj
            njop = ni
            siop = sj
            sjop = si
            iopt = j
            jopt = iold
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, a(k,k), p)
              m = m - 1
            end do
          end if
        end do
      end if
    
      if (iold .lt. lg) then
        i  = iold
        ij = ij + i
        do j = (iold+1), lg
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            termj = ALFLOG
          else
            m = p
            l = icj
            k = ni + 1
 410        continue
            call dcopy( m, x(l,nj), n, v, 1)
            call mclrup( k, m, v, u(nj,nj), p)
            nj = nj + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 410
            nj  = l - n
            rj  = dble(nj)
            if (nj .gt. 2) then
              termj = x(ic(icj),1)
            else
              temp  = dnrm2(p,x(icj,1),n)*s(1)
              temp  = temp*temp
c             termj = rj*log(temp/rj + ALPHA)
              termj = rj*log((temp + ALPHA)/rj)
            end if
            nij = ni + nj
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            sij = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
          end if
          call mclrup( nij, p, v, u, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, u(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(p,nij-1)
          call dgesvd( 'N', 'N', m, p, t, p, v,
     *                  dummy, 1, dummy, 1, w, lwork, info)
c         call dgesvd( 'A', 'A', m, p, t, p, v,
c    *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          temp  = ddot( m, v, 1, v, 1)
c         trmij = rij*log(temp/rij + ALPHA)
          trmij = rij*log((temp + ALPHA)/rij)
          dij   = trmij - (termi + termj)
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = ni
            njop = nj
            siop = si
            sjop = sj
            iopt = iold
            jopt = j
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, a(k,k), p)
              m = m - 1
            end do
          end if
          ij = ij + i
          i  = j
        end do
      end if

c update d and find max

      jopt = 2
      iopt = 1

      dopt = d(1)

      if (lg .eq. 2) goto 900

      ij   = 1
      do i = 2, ld
        qi = d(i)
        if (qi .le. dopt) then
          ij   = i
          dopt = qi
        end if
      end do

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      if (iopt .ne. iold .and. jopt .ne. iold) then

        i   = iopt
        j   = jopt

        nj  = 1
        icj = ic(j)
        ni  = 1
        ici = ic(i)
        if (icj .eq. 1) then
          termj = ALFLOG
          if (ici .eq. 1) then
            nij = 2
            rij = two
            si  = rthalf
            sj  = rthalf
            call dcopy(p, x(i,1), n, v, 1)
            call daxpy( p, (-one), x(j,1), n, v, 1)
            call dscal( p, rthalf, v, 1)
            call dcopy( p, v, 1, r, p)
            termi = ALFLOG
          else
            m = p
            l = ici
 610        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 610
            ni  = l - n
            ri  = dble(ni)
            if (ni .gt. 2) then
              termi = x(ic(ici),1)
            else
              temp  = dnrm2(p,x(ici,1),n)*s(1)
              temp  = temp*temp
c             termi = ri*log(temp/ri + ALPHA)
              termi = ri*log((temp + ALPHA)/ri)
            end if
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            call mclrup( nij, p, v, r, p)
          end if
        else
          m = p
          l = icj
 620      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 620
          nj  = l - n
          rj  = dble(nj)
          if (nj .gt. 2) then
            termj = x(ic(icj),1)
          else
            temp  = dnrm2(p,x(icj,1),n)*s(1)
            temp  = temp*temp
c           termj = rj*log(temp/rj + ALPHA)
            termj = rj*log((temp + ALPHA)/rj)
          end if
          if (ici .eq. 1) then
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qj  = rj*qij
            si  = sqrt(qij)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            termi = ALFLOG
          else 
            m = p
            l = ici
            k = nj + 1
 630        continue
            call dcopy( m, x(l,ni), n, v, 1)
            call mclrup( k, m, v, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 630
            ni  = l - n
            ri  = dble(ni)
            if (ni .gt. 2) then
              termi = x(ic(ici),1)
            else
              temp  = dnrm2(p,x(ici,1),n)*s(1)
              temp  = temp*temp
c             termi = ri*log(temp/ri + ALPHA)
              termi = ri*log((temp + ALPHA)/ri)
            end if
            nij = ni + nj
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          end if
          call mclrup( nij, p, v, r, p)
        end if

        trop = trcij
        tmop = dopt + (termi + termj)
        nopt = nij
        niop = ni
        njop = nj
        siop = si
        sjop = sj

      else
        m = p
        do k = 1, min(nopt-1,p)
          call dcopy(m, a(k,k), p, r(k,k), p)
          m = m - 1
        end do
      end if

      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 200

 900  continue

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 1
      ld          = lo
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1,1) = si 
        x(1,2) = sj
      else
        x(1,1) = sj
        x(1,2) = si
      end if

      lg = ng + 1
      do k = 2, ns
        lo     = lo - 1
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 1
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then
          x(k,1) = dble(ici)
          x(k,2) = dble(icj)
        else
          x(k,1) = dble(icj)
          x(k,2) = dble(ici)
        end if
      end do

      ld = nd
      lo = 1
      do k = 1, ns
        si    = d(lo)
        d(lo) = d(ld)
        d(ld) = si
        ld    = ld - 1
        lo    = lo + 1
      end do

      nd = 0

      lwork = lwopt

      return
      end
      subroutine hcvi  ( x, n, p, ic, ng, ns, ALPHA, v, nd, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithms in which shape and orientation 
c are fixed in advance, while volume may vary among clusters.

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns

c     double precision   x(n,p), v(p). d(*)
      double precision   x(n,*), v(*), d(*)
c------------------------------------------------------------------------------
c
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  ALPHA   double  (input) Additive quantity used to resolve degeneracies. 
c  v       double  (scratch) (p).
c  nd      integer (input) The length of d.
c  d       double  (scratch/output) max(n,((ng-1)*(ng-2))/2,3*ns). On output
c                   the first ns elements are proportional to the change in 
c                   loglikelihood associated with each merge.

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        sqrthf
      parameter              (sqrthf = .70710678118654757274D0)

      double precision        d1mach

c------------------------------------------------------------------------------

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd

c     FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
c     EPSMIN = d1mach(3)
      EPSMAX = d1mach(4)

      if (ng .eq. 1) return

      ALPHA  = max(ALPHA,EPSMAX)

      ALFLOG = log(ALPHA)

c     call intpr( 'ic', -1, ic, n) 

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. lg)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers
     
      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
c update sum of squares
          k = ic(i)
          if (k .eq. 1) then
            ic(i) = j
            ic(j) = 2
            call dscal( p, sqrthf, x(i,1), n)
            call dscal( p, sqrthf, x(j,1), n)
            call dcopy( p, x(j,1), n, v, 1)
            call daxpy( p, (-one), x(i,1), n, v, 1)
            call daxpy( p, one, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
c           x(j,1) = ddot( p, v, 1, v, 1) / two
            x(j,1) = ddot( p, v, 1, v, 1)
          else
            ic(j) = 0
            ni    = ic(k)
            ic(k) = ni + 1
            ri    = dble(ni)
            rij   = dble(ni+1)
            qj    = one/rij
            qi    = ri*qj
            si    = sqrt(qi)
            sj    = sqrt(qj)
            call dcopy( p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(i,1), n, v, 1)
c           x(k,1) = qi*x(k,1) + qj*ddot(p, v, 1, v, 1)
            x(k,1) =    x(k,1) +    ddot(p, v, 1, v, 1)
            call dscal( p, si, x(i,1), n)
            call daxpy( p, sj, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
          end if 
        else 
          ic(j) = 1
        end if
      end do
       
c store terms also so as not to recompute them

      do k = 1, ng
        i = ic(k)
        if (i .ne. 1) then
          ni     = ic(i)
          ri     = dble(ni)
c         x(i,2) = ri*log(x(i,1)+ALPHA)
          x(i,2) = ri*log((x(i,1)+ALPHA)/ri)
        end if
      end do

c     call intpr( 'ic', -1, ic, n)
c     call dblepr( 'trace', -1, x(1,1), n)
c     call dblepr( 'term', -1, x(1,2), n)

c compute change in likelihood and determine minimum

      dopt = FLMAX

      ij   = 0
      do j = 2, ng
        nj = ic(j)
        if (nj .eq. 1) then
          tracej = zero
          termj  = ALFLOG
          rj     = one
        else
          tracej = x(nj,1)
          termj  = x(nj,2)
          nj     = ic(nj)
          rj     = dble(nj)
        end if
        do i = 1, (j-1)
          ni = ic(i)
          if (ni .eq. 1) then
            tracei = zero
            termi  = ALFLOG
            ri     = one
          else 
            tracei = x(ni,1)
            termi  = x(ni,2)
            ni     = ic(ni)
            ri     = dble(ni)
          end if               
          nij = ni + nj
          rij = dble(nij)
          qij = one/rij
          qi  = ri*qij
          qj  = rj*qij
          si  = sqrt(qi)
          sj  = sqrt(qj)
          call dcopy(p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
c         trcij = (qi*tracei + qj*tracej) + qij*ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + ddot(p,v,1,v,1)
c         trmij = rij*log(trcij+ALPHA)
          trmij = rij*log((trcij+ALPHA)/rij)
          dij   = trmij - (termi + termj)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            rijo = rij
            nopt = nij
            niop = ni
            njop = nj
            siop = si
            sjop = sj
            iopt = i
            jopt = j
          end if
        end do
      end do

c     call dblepr( 'dij', -1, d, (l*(l-1))/2)

      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1,1) = dble(iopt)
          x(1,2) = dble(jopt)
        else
          x(1,1) = dble(jopt)
          x(1,2) = dble(iopt)
        end if
        d(1)   = dopt
        return
      end if

      if (niop .ne. 1) ic(ic(iopt)) = 0
      if (njop .ne. 1) ic(ic(jopt)) = 0

      ls = 1

 100  continue

c     if (.false.) then
c       ij = 1
c       jj = 1
c       do j = 2, n
c         nj = ic(j)
c         if (nj .ne. 0 .and. abs(nj) .le. n) then
c           call dblepr( 'dij', -1, d(ij), jj)
c           ij = ij + jj 
c           jj = jj + 1
c         end if
c       end do
c     end if

      call dscal( p, siop, x(iopt,1), n)
      call daxpy( p, sjop, x(jopt,1), n, x(iopt,1), n)

      if (jopt .ne. lg) then
        call mclcpy( jopt, lg, d)
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        m        = ic(jopt)
        ic(jopt) = ic(lg)
        ic(lg)   = m
      end if

      ic(iopt) =  lg
      ic(lg)   =  nopt
      x(lg,1)  =  trop
      x(lg,2)  =  tmop

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)
      lo     = lo - 1

      lg = lg - 1
      ld = ld - lg

      iold     = iopt

      iopt     = -1
      jopt     = -1

      dopt   = FLMAX

      ni     = nopt
      ri     = dble(ni)
      tracei = trop
      termi  = tmop

      ij = ((iold-1)*(iold-2))/2
      if (iold .gt. 1) then
        do j = 1, (iold - 1)
          nj = ic(j)
          if (nj .ne. 1) then
            tracej = x(nj,1)
            termj  = x(nj,2)
            nj     = ic(nj)
            rj   = dble(nj)
          else 
            tracej = zero
            termj  = ALFLOG
            rj     = one
          end if
          nij = ni + nj
          rij = dble(nij)
          qij = one/rij
          qi  = ri*qij
          qj  = rj*qij
          si  = sqrt(qi)
          sj  = sqrt(qj)
          call dcopy( p, x(iold,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
c         trcij = (qi*tracei + qj*tracej) + qij*ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + ddot(p,v,1,v,1)
c         trmij = rij*log(trcij+ALPHA)
          trmij = rij*log((trcij+ALPHA)/rij)
          dij   = trmij - (termi + termj)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            iopt = j
            jopt = iold
            nopt = nij
            niop = ni
            njop = nj
            sjop = si
            siop = sj
          end if
        end do
      end if

      if (iold .lt. lg) then
        i  = iold
        ij = ij + i
        do j = (iold + 1), lg
          nj = ic(j)
          if (nj .ne. 1) then
            tracej = x(nj,1)
            termj  = x(nj,2)
            nj     = ic(nj)
            rj     = dble(nj)
          else 
            tracej = zero
            termj  = ALFLOG
            rj     = one
          end if
          nij = ni + nj
          rij = dble(nij)
          qij = one /rij
          qi  = ri*qij
          qj  = rj*qij
          si  = sqrt(qi)
          sj  = sqrt(qj)
          call dcopy( p, x(iold,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
c         trcij = (qi*tracei + qj*tracej) + qij*ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + ddot(p,v,1,v,1)
c         trmij = rij*log(trcij+ALPHA)
          trmij = rij*log((trcij+ALPHA)/rij)
          dij   = trmij - (termi + termj)
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            iopt = iold
            jopt = j
            nopt = nij
            niop = ni
            njop = nj
            siop = si
            sjop = sj
          end if
          ij = ij + i
          i  = j
        end do
      end if

c update d and find max

      jopt = 2
      iopt = 1

      dopt = d(1)

      if (lg .eq. 2) goto 900

      ij   = 1
      dopt = d(1)
      do i = 2, ld
        qi = d(i)
        if (qi .le. dopt) then
          ij   = i
          dopt = qi
        end if
      end do

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      i   = ic(iopt)
      j   = ic(jopt)

      if (iopt .ne. iold .and. jopt .ne. iold) then
  
        if (i .ne. 1) then
          tracei = x(i,1)
          termi  = x(i,2)
          niop   = ic(i)
          ri     = dble(niop)
        else
          tracei = zero
          termi  = ALFLOG
          niop   = 1
          ri     = one
        end if

        if (j .ne. 1) then
          tracej = x(j,1)
          termj  = x(j,2)
          njop   = ic(j)
          rj     = dble(njop)
        else
          tracej = zero
          termj  = ALFLOG
          njop   = 1
          rj     = one
        end if

        nopt = niop + njop
        rij  = dble(nopt)
        qij  = one/rij
        qi   = ri*qij
        qj   = rj*qij
        siop = sqrt(qi)
        sjop = sqrt(qj)

        call dcopy( p, x(iopt,1), n, v, 1)
        call dscal( p, sjop, v, 1)
        call daxpy( p, (-siop), x(jopt,1), n, v, 1)
c       trop  = (qi*tracei + qj*tracej) + qij*ddot(p,v,1,v,1)
        trop  = (tracei + tracej) + ddot(p,v,1,v,1)
c       tmop  = rij*log(trop+ALPHA)
        tmop  = rij*log((trop+ALPHA)/rij)
      end if
     
      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 100

 900  continue

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 1
      ld          = lo
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1,1) = si 
        x(1,2) = sj
      else
        x(1,1) = sj
        x(1,2) = si
      end if

      lg = ng + 1
      do k = 2, ns
        lo     = lo - 1
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 1
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then
          x(k,1) = dble(ici)
          x(k,2) = dble(icj)
        else
          x(k,1) = dble(icj)
          x(k,2) = dble(ici)
        end if
      end do

      ld = nd
      lo = 1
      do k = 1, ns
        si    = d(lo)
        d(lo) = d(ld)
        d(ld) = si
        ld    = ld - 1
        lo    = lo + 1
      end do

      return
      end

















      subroutine hcvvv ( x, n, p, ic, ng, ns, ALPHA, BETA, 
     *                   v, u, s, r, nd, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Gaussian model-based clustering algorithm in which shape, volume, and
c orientation are allowed to vary between clusters.

      implicit double precision (a-h,o-z)

      integer            n, p, ic(n), ng, ns, nd

      double precision   ALPHA, BETA

c     double precision   x(n,p+1), v(p), u(p,p), s(p,p)
c     double precision   r(p,p), d(ng*(ng-1)/2)
      double precision   x(n,*), v(*), u(p,*), s(p,*)
      double precision   r(p,*), d(*)
c------------------------------------------------------------------------------

c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations with a column appended for scratch use.
c                   On output, the first two columns and ns rows contain the 
c                   merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ALPHA   double  (input) Additive quantity used to resolve degeneracies.
c  BETA    double  (input) Factor by which to multiply the trace which is used
c                   additively to help resolve degeneracies.
c  ng      integer (input) Number of groups in initial partition.
c  ns      integer (input) Desired number of stages of clustering.
c  v       double  (scratch) (p) 
c  u,s,r   double  (scratch) (p*p)
c  nd      integer (input) The length of d.
c  d       double  (scratch/output) max(p*p+n,((ng*(ng-1))/2,3*ns). On output
c                   the first ns elements are proportional to the change in
c                   loglikelihood associated with each merge.

      integer                 psq, pm1, pp1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        rthalf
      parameter              (rthalf = .7071067811865476d0)

      double precision        d1mach

      common /VVVMCL/         BETA0, ALPHA0, ABLOG
      save   /VVVMCL/            

      common /MCLMCH/         FLMAX
      save   /MCLMCH/  

c------------------------------------------------------------------------------

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd
     
c     FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)
c     EPSMIN = d1mach(3)
      EPSMAX = d1mach(4)

      psq = p*p
      pm1 = p-1
      pp1 = p+1

      if (ng .eq. 1) return

      ALPHA  = max(ALPHA,EPSMAX)

      BETA0  = BETA
      ALPHA0 = ALPHA

      ABLOG  = log(BETA*ALPHA)

c     call intpr( 'ic', -1, ic, n)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     call intpr( 'ic', -1, ic, n)

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, v, 1)
          call dscal( p, si, v, 1)
          call daxpy( p, (-sj), x(k,1), n, v, 1)
          trcij =    trcij +    ddot( p, v, 1, v, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, v, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c         d(ll+k) = trcij
c copy triangular factor into the rows of x
          j  = k
          m  = p
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij)    = n+l
          x(k, pp1) = zero
          if (l .ge. 2) then
            x(   k, pp1) = trcij
            trmij        = vvvtij( l, p, r, sj, trcij)
            x(ic(k),pp1) = trmij
          end if
        else
          ic(k)   = 1
c         d(ll+k) = zero
        end if
      end do

c     call intpr( 'ic', -1, ic, n)
c     call dblepr( '', -1, x(1,pp1), n)
c     call dblepr( 'trac', -1, d(ll+1), ng)
c     call dblepr( 'term', -1, term, n)

c compute change in likelihood and determine minimum

      dopt = FLMAX

      ij = 0
      do j = 2, ng
        icj = ic(j)
        nj  = 1
        if (icj .eq. 1) then
          tracj = zero
          termj = ABLOG
          do i = 1, (j-1)
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij   = 2
              rij   = two
              si    = rthalf
              sj    = rthalf
              sij   = rthalf
              call dcopy( p, x(i,1), n, v, 1)
              call daxpy( p, (-one), x(j,1), n, v, 1)
              call dscal( p, rthalf, v, 1)
c             trcij = half*ddot( p, v, 1, v, 1)
              trcij =      ddot( p, v, 1, v, 1)
              call dcopy( p, v, 1, u, p)  
c             trmij = rij*log(BETA*trcij+ALPHA)
              trmij = two*log(BETA*(trcij+ALPHA)/two)
              termi = ABLOG
            else
              m  = p
              l  = ici
 110          continue
              call dcopy( m, x(l,ni), n, u(ni,ni), p)
              ni  = ni + 1
              m   = m - 1
              l   = ic(l)
              if (l .le. n) goto 110
              ni  = l - n
c             traci = d(ll+i)
c             traci = trac(i)
c             termi = vvvtrm(i,ni,n,p,ic,x,traci)
c             termi = term(i)
              traci = x(   i , pp1)
              termi = x(ic(i), pp1)
              ri  = dble(ni)
              nij = ni + 1
              rij = dble(nij)
              qij = one/rij
              qi  = ri*qij
              si  = sqrt(qi)
              sj  = sqrt(qij)
              sij = sj
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              trcij =    traci +     ddot(p,v,1,v,1)
              call mclrup( nij, p, v, u, p)
              trmij = vvvtij( nij, p, u, sij, trcij)
            end if
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij
            if (dij .le. dopt) then
              dopt = dij
              trop = trcij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, s(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj    = l - n
c         tracj = d(ll+j)
c         termj = vvvtrm(j,nj,n,p,ic,x,tracj)
          tracj = x(    j , pp1)
          termj = x( ic(j), pp1)
          rj    = dble(nj)
          do i = 1, (j-1)
            m = p
            do k = 1, min(nj-1,p)
              call dcopy( m, s(k,k), p, u(k,k), p)
              m = m - 1
            end do
            ni  = 1
            ici = ic(i)
            if (ici .eq. 1) then
              nij = nj + 1
              rij = dble(nij)
              qij = one/rij
              qi  = qij
              qj  = rj*qij
              si  = sqrt(qi)
              sj  = sqrt(qj)
              sij = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              trcij =    tracj +     ddot(p,v,1,v,1)
              termi = ABLOG
            else
              m  = p
              l  = ici
              k  = nj + 1
 130          continue
              call dcopy( m, x(l,ni), n, v, 1)
              call mclrup( k, m, v, u(ni,ni), p)
              ni    = ni + 1
              m     = m - 1
              l     = ic(l)
              if (l .le. n) goto 130
              ni    = l - n
c             traci = d(ll+i)
c             termi = vvvtrm(i,ni,n,p,ic,x,traci)
              traci = x(   i , pp1)
              termi = x(ic(i), pp1)
              ri    = dble(ni)
              nij   = ni + nj
              rij   = dble(nij)
              qij   = one/rij
              qi    = ri*qij
              qj    = rj*qij
              si    = sqrt(qi)
              sj    = sqrt(qj)
              sij   = sqrt(qij)
              call dcopy(p, x(i,1), n, v, 1)
              call dscal( p, sj, v, 1)
              call daxpy( p, (-si), x(j,1), n, v, 1)
              trcij = (   traci +    tracj) +     ddot(p,v,1,v,1)
            end if
            call mclrup( nij, p, v, u, p)
            trmij = vvvtij( nij, p, u, sij, trcij)
            dij   = trmij - (termi + termj)
            ij    = ij + 1
            d(ij) = dij 
            if (dij .le. dopt) then
              dopt = dij
              trop = trcij
              tmop = trmij
              nopt = nij
              niop = ni
              njop = nj
              siop = si
              sjop = sj
              iopt = i
              jopt = j
              m = p
              do k = 1, min(nij-1,p)
                call dcopy( m, u(k,k), p, r(k,k), p)
                m = m - 1
              end do
            end if
          end do
        end if
      end do

c     if (.false.) then
c       i  = 1
c       ij = 1
c       do j = 2, ng
c         call dblepr( 'dij', -1, d(ij), i)
c         ij = ij + i
c         i  = j
c       end do
c     end if
 
      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1,1) = iopt
          x(1,2) = jopt
        else
          x(1,1) = jopt
          x(1,2) = iopt
        end if
        d(1)   = dopt
        return
      end if

      ls  = 1

 200  continue

      call dcopy( p, x(iopt,1), n, v, 1)
      call dscal( p, siop, v, 1)
      call daxpy( p, sjop, x(jopt,1), n, v, 1)

      if (jopt .ne. lg) then
        call mclcpy( jopt, lg, d)
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        m          = ic(jopt)
        icj        = ic(lg)
        if (icj .ne. 1) x( jopt, pp1) = x( lg, pp1)
        ic(jopt)   = icj
        ic(lg)     = m
c       term(jopt) = term(lg)
c       trac(jopt) = trac(lg)
      end if

c     term(lg)   = FLMAX
c     trac(lg)   = FLMAX

      call dcopy( p, r(1,1), p, x(lg,1), n)

      if (niop .eq. 1) then

        if (njop .eq. 1) then
          ic(lg)  = n+2
        else
          l   = ic(lg)
          m   = pm1
          nij = 2
 210      continue
            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

        end if

      else

        l   = ic(iopt)
        m   = pm1
        nij = 2
 220    continue
          call dcopy( m, r(nij,nij), p, x(l,nij), n)
          nij = nij + 1
          m   = m   - 1
          k   = l
          l   = ic(l)
          if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

          if (nij .le. p .and. njop .ne. 1) then
            l     = ic(lg)
            ic(k) = l
 230        continue
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij   = nij + 1
              m     = m   - 1
              k     = l
              l     = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230
          end if

        ic(lg) = ic(iopt)
        ic(k)  = n + nopt

      end if

c     term(iopt) = tmop
c     trac(iopt) = trop

      ic(iopt) = lg

      x(iopt, pp1) = zero
      if (nopt .ge. 2) then
        x(iopt,pp1)     = trop
        x(ic(iopt),pp1) = tmop
      endif   

      call dcopy( p, v, 1, x(iopt,1), n)

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)
      lo     = lo - 1

      lg = lg - 1
      ld = ld - lg

c     call intpr( 'ic', -1, ic, n)

      iold  =  iopt

      dopt  = FLMAX

      ni    = nopt
      ri    = dble(ni)
      termi = tmop
      traci = trop

      ij = ((iold-1)*(iold-2))/2
      if (iold .gt. 1) then
        do j = 1, (iold-1)
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            trcij =    traci +     ddot(p,v,1,v,1)
            tracj = zero
            termj = ABLOG
         else
            m = p
            l = icj
            k = ni + 1
 310        continue
              call dcopy( m, x(l,nj), n, v, 1)
              call mclrup( k, m, v, u(nj,nj), p)
              nj = nj + 1
              m  = m - 1
              l  = ic(l)
              if (l .le. n) goto 310
            nj  = l - n            
c           call vvvget(j,nj,n,p,ic,x,tracj,termj)
            tracj = x(   j ,pp1)
            termj = x(ic(j),pp1)
            rj    = dble(nj)
            nij   = ni + nj
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            qj    = rj*qij
            si    = sqrt(qi)
            sj    = sqrt(qj)
            sij   = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            trcij = (   traci +    tracj) +     ddot(p,v,1,v,1)
          end if
          call mclrup( nij, p, v, u, p)
          trmij = vvvtij( nij, p, u, sij, trcij)
          dij   = trmij - (termi + termj)
          ij    = ij + 1
          d(ij) = dij
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = nj
            njop = ni
            siop = sj
            sjop = si
            iopt = j
            jopt = iold
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, s(k,k), p)
              m = m - 1
            end do
          end if
        end do
      end if
    
      if (iold .lt. lg) then
        i  = iold
        ij = ij + i
        do j = (iold+1), lg
          m = p
          do k = 1, min(ni-1,p)
            call dcopy(m, r(k,k), p, u(k,k), p)
            m = m - 1
          end do
          nj  = 1
          icj = ic(j)
          if (icj .eq. 1) then
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            trcij =    traci +     ddot(p,v,1,v,1)
            termj = ABLOG
          else
            m = p
            l = icj
            k = ni + 1
 410        continue
            call dcopy( m, x(l,nj), n, v, 1)
            call mclrup( k, m, v, u(nj,nj), p)
            nj = nj + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 410
            nj  = l - n
c           call vvvget(j,nj,n,p,ic,x,tracj,termj)
            tracj = x(   j ,pp1)
            termj = x(ic(j),pp1)
            rj    = dble(nj)
            nij   = ni + nj
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            qj    = rj*qij
            si    = sqrt(qi)
            sj    = sqrt(qj)
            sij   = sqrt(qij)
            call dcopy(p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(iold,1), n, v, 1)
            trcij = (   traci +    tracj) +     ddot(p,v,1,v,1)
          end if
          call mclrup( nij, p, v, u, p)
          trmij = vvvtij( nij, p, u, sij, trcij)
          dij   = trmij - (termi + termj)
          d(ij) = dij
          ij    = ij + i
          i     = j
          if (dij .le. dopt) then
            dopt = dij
            trop = trcij
            tmop = trmij
            nopt = nij
            niop = ni
            njop = nj
            siop = si
            sjop = sj
            iopt = iold
            jopt = j
            m = p
            do k = 1, min(nij-1,p)
              call dcopy(m, u(k,k), p, s(k,k), p)
              m = m - 1
            end do
          end if
        end do
      end if

c update d and find max

      jopt = 2
      iopt = 1

      dopt = d(1)

      if (lg .eq. 2) goto 900
        
      ij   = 1
      do i = 2, ld
        qi = d(i)
        if (qi .le. dopt) then
          ij   = i
          dopt = qi
        end if
      end do

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      if (iopt .ne. iold .and. jopt .ne. iold) then

        i   = iopt
        j   = jopt

        nj  = 1
        icj = ic(j)
        ni  = 1
        ici = ic(i)
        if (icj .eq. 1) then
          termj = ABLOG
          if (ici .eq. 1) then
            nij = 2
            rij = two
            si  = rthalf
            sj  = rthalf
            call dcopy(p, x(i,1), n, v, 1)
            call daxpy( p, (-one), x(j,1), n, v, 1)
            call dscal( p, rthalf, v, 1)
            trcij =      ddot( p, v, 1, v, 1)
            call dcopy( p, v, 1, r, p)
            termi = ABLOG
          else
            m = p
            l = ici
 610        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 610
            ni  = l - n
c           call vvvget(i,ni,n,p,ic,x,traci,termi)
            traci = x(   i , pp1)
            termi = x(ic(i), pp1)
            ri  = dble(ni)
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            trcij =    traci +     ddot( p, v, 1, v, 1)
            call mclrup( nij, p, v, r, p)
          end if
        else
          m = p
          l = icj
 620      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 620
          nj  = l - n
c         call vvvget(j,nj,n,p,ic,x,tracj,termj)
          tracj = x(   j , pp1)
          termj = x(ic(j), pp1)
          rj  = dble(nj)
          if (ici .eq. 1) then
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qj  = rj*qij
            si  = sqrt(qij)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            trcij =    tracj +     ddot( p, v, 1, v, 1)
            termi = ABLOG
          else 
            m = p
            l = ici
            k = nj + 1
 630        continue
            call dcopy( m, x(l,ni), n, v, 1)
            call mclrup( k, m, v, r(ni,ni), p)
            ni = ni + 1
            m  = m - 1
            l  = ic(l)
            if (l .le. n) goto 630
            ni  = l - n
c           call vvvget(i,ni,n,p,ic,x,traci,termi)
            traci = x(   i , pp1)
            termi = x(ic(i), pp1)
            ri  = dble(ni)
            nij = ni + nj
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            trcij = (   traci +    tracj) +     ddot( p,v,1,v,1)
          end if
          call mclrup( nij, p, v, r, p)
        end if

        trop = trcij
        tmop = dopt + (termi + termj)
        nopt = nij
        niop = ni
        njop = nj
        siop = si
        sjop = sj

      else
        m = p
        do k = 1, min(nopt-1,p)
          call dcopy(m, s(k,k), p, r(k,k), p)
          m = m - 1
        end do
      end if

      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 200

 900  continue

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 1
      ld          = lo
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1,1) = si
        x(1,2) = sj
      else
        x(1,1) = sj
        x(1,2) = si
      end if

      lg = ng + 1
      do k = 2, ns
        lo     = lo - 1
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 1
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then 
          x(k,1) = dble(ici)
          x(k,2) = dble(icj)
        else
          x(k,1) = dble(icj)
          x(k,2) = dble(ici)
        end if
      end do

      ld = nd
      lo = 1
      do k = 1, ns
        si    = d(lo)
        d(lo) = d(ld)
        d(ld) = si
        ld    = ld - 1
        lo    = lo + 1
      end do

      return
      end

      double precision function vvvtij( l, p, r, s, trac)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

      implicit double precision (a-h,o-z)

      integer                    l, p
      double precision           r(p,*)

      double precision           zero, one
      parameter                 (zero = 0.d0, one = 1.d0)

      common /VVVMCL/            BETA, ALPHA, ABLOG
      save   /VVVMCL/            

      common /MCLMCH/            FLMAX
      save   /MCLMCH/            

      if (l .le. p) then
        vvvtij = log(BETA*(trac+ALPHA)/dble(l))
      else 
        if (trac .eq. zero) then
          vvvtij = log((ALPHA*BETA)/dble(l))
        else
          detlog = det2mc( p, r, s)
          if (detlog .eq. (-FLMAX)) then
            vvvtij = log(BETA*(trac+ALPHA)/dble(l))
          else if (detlog .le. zero) then
            vvvtij = log(exp(detlog)+BETA*(trac+ALPHA)/dble(l))
          else
            vvvtij = log(one+exp(-detlog)*(BETA*(trac+ALPHA)/dble(l))) 
     *              + detlog
          end if
        end if
      end if

      vvvtij = dble(l)*vvvtij

      return
      end
      subroutine likeee( cltree, ns, x, n, p, ic, ng, v, r, new, hood)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which
c clusters share a common variance (same shape, volume, and orientation)

      implicit double precision (a-h,o-z)

      integer            cltree(2,ns), ns, n, p, ic(n), ng, new(ns)

c     double precision   x(n,p), v(p), r(p,p), hood(ns+1)
      double precision   x(n,*), v(*), r(*), hood(*)
c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  v       double  (scratch) (p).
c  r       double  (scratch) (p*p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) Twice the loglikelihood (up to an additive 
c                   constant) for each stage.

      integer                 lw, q

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision       d1mach, detmc2

      double precision       FLMAX
      common /MCLMCH/        FLMAX
      save   /MCLMCH/

c------------------------------------------------------------------------------

      FLMAX = d1mach(2)

      dn2   = dble(n)/two

      ns1   = ns + 1

      lw    = p*p
      
c form scaled column sums
      call dscal( n*p, one/sqrt(dble(n)), x, 1)

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c     trcw = zero
      call dcopy( lw, zero, 0, r, 1)

      q = 0
      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
c update trace and Cholesky factor as if a merge
c         q     = q + 2
          q     = q + 1
          ni    = ic(i)
          ri    = dble(ni)
          rij   = dble(ni+1)
          sj    = sqrt(one/rij)
          si    = sqrt(ri)*sj
          call dcopy( p, x(i,1), n, v, 1)
          call dscal( p, sj, v, 1)
          call daxpy( p, (-si), x(j,1), n, v, 1)
c         trcw  = trcw + ddot(p, v, 1, v, 1)
          call mclrup( q+1, p, v, r, p)
          ic(j) = 0
          ic(i) = ic(i) + 1
          call dscal( p, si, x(i,1), n)
          call daxpy( p, sj, x(j,1), n, x(i,1), n)
c         call dcopy( p, FLMAX, 0, x(j,1), n)
c update column sum in jth row
        else
          ic(j) = 1
        end if
      end do

      if (q .le. p) then
        detlog  = -FLMAX
        hood(1) = -FLMAX
      else
        detlog  = detmc2( p, r)
        if (detlog .eq. -FLMAX) then
          hood(1) = -FLMAX
        else 
          hood(1) = -dn2 * detlog
        end if
      end if        

      k = 1
      do l = 2, ns1

        i   = cltree(1,k)
        j   = cltree(2,k)

        nj  = ic(j)
        rj  = dble(nj)

        ni  = ic(i)
        ri  = dble(ni)
        nij = ni + nj
        rij = dble(nij)
        si  = sqrt(ri/rij)
        sj  = sqrt(rj/rij)
        call dcopy( p, x(i,1), n, v, 1)
        call dscal( p, sj, v, 1)
        call daxpy( p, (-si), x(j,1), n, v, 1)
c       trcw = trcw + ddot(p, v, 1, v, 1)

        call dscal( p, si, x(i,1), n)
        call daxpy( p, sj, x(j,1), n, x(i,1), n)

c update the Cholesky factor

        q  =  q + 1
        call mclrup( q+1, p, v, r, p)

        ic(i)  = nij
        ic(j)  = 0

        new(k) = nij

        if (q .le. p) then
          detlog  = -FLMAX
          hood(l) = -FLMAX
        else
          detlog  = detmc2( p, r)
          if (detlog .eq. (-FLMAX)) then
            hood(l) = -FLMAX
          else 
            hood(l) = -dn2 * detlog
          end if
        end if        

        k = l
      end do

      return
      end
      subroutine likefv( cltree, ns, x, n, p, ic, ng, s, lwork, w,
     *                   v, t, r, uu, vv, new, hood, info)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which
c volume is fixed but shape and orientation may vary between clusters 

      implicit double precision (a-h,o-z)

      integer            cltree(2,ns),ns,n,p,ic(n),ng,new(ns),info

c     double precision   x(n,p),s(p),v(p),w(lwork),t(p,p),r(p,p),hood(ns+1)
      double precision   x(n,*),s(*),v(*),w(*),t(p,*),r(p,*),hood(*)

c     double precision   uu(p,p), vv(p,p)
      double precision   uu(p,*), vv(p,*)

c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  s       double  (input) (p) Shape matrix. It's square is the A matrix of
c                   Banfield and Raftery.
c  v       double  (scratch) (p).
c  w       double  (scratch) max(4*p,5*p-4).
c  t,r     double  (scratch) (p*p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) Twice the loglikelihood (up to an additive
c                   constant) for each stage.
c  info    integer (output) return code from the LAPACK svd routine

      integer                 psq, pm1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision       half
      parameter             (half = .5d0)

      double precision       rthalf
      parameter             (rthalf = .7071067811865476d0)

      double precision       d1mach

      common /MCLMCH/        FLMAX
      save   /MCLMCH/  

c------------------------------------------------------------------------------

      do i = 1, p
        s(i) = one/s(i)
      end do

c     lwork = max(4*p,5*p-4)

      lwopt = 0

      FLMAX = d1mach(2)

      dnp2   = half*dble(n*p)

      ns1   = ns + 1

      psq   = p*p
      pm1   = p-1

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     nsing = 0

      sum   = zero

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, v, 1)
          call dscal( p, si, v, 1)
          call daxpy( p, (-sj), x(k,1), n, v, 1)
c         trcij = qi*trcij + qj*ddot( p, v, 1, v, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, v, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c         term(k) = trmij
c         trac(k) = trcij
c         d(ll+k) = trcij
c copy triangular factor into the rows of x
          j  = k
          m  = p
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij) = n+l
          if (l .gt. 2) then
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do i = 1, min(l-1,p)
              call dcopy( m, r(i,i), p, t(i,i), p)
              m = m - 1
            end do
            m = min(p,l-1)
c           call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                    dummy, 1, dummy, 1, w, lwork, info)
            call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do i = 1, m
              v(i) = v(i)*s(i)
            end do
            temp = ddot( m, v, 1, v, 1)
          else 
            temp = dnrm2( p, r, p)*s(1)
            temp = temp*temp
          end if
          termij = temp
          if (l .gt. 2) x(ic(ic(k)), 1) = termij
        else
          ic(k) = 1
          termij = zero
        end if 
        sum = sum + termij
      end do

      if (sum .ne. zero) then
         hood(1) = dnp2*log(sum)
      else
        sum     =  zero
        hood(1) = -FLMAX
      end if
 
      kk = 1
      do ll = 2, ns1
    
        i   = cltree(1,kk)
        j   = cltree(2,kk)
 
        ni  = 1
        nj  = 1

        icj = ic(j)
        ici = ic(i)

        if (icj .eq. 1) then
          termj = zero
          if (ici .eq. 1) then
            termi = zero
            nij   = 2
            rij   = two
            si    = rthalf
            sj    = rthalf
            sij   = rthalf
            call dcopy( p, x(i,1), n, v, 1)
            call daxpy( p, (-one), x(j,1), n, v, 1)
            call dscal( p, rthalf, v, 1)
            call dcopy( p, v, 1, r(1,1), p)
            temp = dnrm2( p, v, 1)*s(1)
            temp = temp*temp
          else
            m  = p
            l  = ici
 110        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni  = ni + 1
            m   = m - 1
            l   = ic(l)
            if (l .le. n) goto 110
            ni    = l - n
            if (ni .gt. 2) then
              termi = x(ic(ici), 1)
            else
              temp  = dnrm2( p, r, p)*s(1)
              termi = temp*temp
            end if
            ri    = dble(ni)
            nij   = ni + 1
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            si    = sqrt(qi)
            sj    = sqrt(qij)
            sij   = sj
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            call mclrup( nij, p, v, r, p)
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do k = 1, min(p,nij-1)
              call dcopy( m, r(k,k), p, t(k,k), p)
              m = m - 1
            end do
            m = min(p,nij-1)
c           call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                    dummy, 1, dummy, 1, w, lwork, info)
            call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do k = 1, m
              v(k) = v(k)*s(k)
            end do
            temp = ddot( m, v, 1, v, 1)
          end if
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj = l - n
          if (nj .gt. 2) then
            termj = x(ic(icj), 1)
          else
            temp  = dnrm2( p, r, p)*s(1)
            termj = temp*temp
          end if
          rj  = dble(nj)
          ni  = 1
          ici = ic(i)
          if (ici .eq. 1) then
            termi = zero
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qi  = qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            sij = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          else
            m  = p
            l  = ici
            k  = nj + 1
 130        continue
            call dcopy( m, x(l,ni), n, v, 1)
            call mclrup( k, m, v, r(ni,ni), p)
            ni    = ni + 1
            m     = m - 1
            l     = ic(l)
            if (l .le. n) goto 130
            ni    = l - n
            if (ni .gt. 2) then
              termi = x(ic(ici), 1)
            else
              temp  = dnrm2( p, v, 1)*s(1)
              termi = temp*temp
            end if
            ri    = dble(ni)
            nij   = ni + nj
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            qj    = rj*qij
            si    = sqrt(qi)
            sj    = sqrt(qj)
            sij   = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          end if
          call mclrup( nij, p, v, r, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, r(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(nij-1,p)
c         call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                  dummy, 1, dummy, 1, w, lwork, info)
          call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          temp = ddot( m, v, 1, v, 1)
        end if

        termij = temp

        call dcopy( p, x(i,1), n, v, 1)
        call dscal( p, si, v, 1)
        call daxpy( p, sj, x(j,1), n, v, 1)

        call dcopy( p, r(1,1), p, x(j,1), n)

        nopt = nij

        if (ni .eq. 1) then

          if (nj .eq. 1) then
            ic(j)  = n+2
          else
            l   = icj
            m   = pm1
            nij = 2

 210        continue
            
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij = nij + 1
              m   = m   - 1
              k   = l
              l   = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

          end if

        else

          l   = ic(i)
          m   = pm1
          nij = 2

 220      continue

            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

            if (nij .le. p .and. nj .ne. 1) then

              l     = ic(j)
              ic(k) = l

 230          continue

                call dcopy( m, r(nij,nij), p, x(l,nij), n)
                nij = nij + 1
                m   = m   - 1
                k   = l
                l   = ic(l)
                if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230

            end if

          ic(j) = ic(i)
          ic(k) = n + nopt

        end if

        ic(i) = j

        if (nopt .gt. 2) then
          m      = ic(j)
          x(m,1) = termij
        endif

        call dcopy( p, v, 1, x(i,1), n)

        new(kk) = nopt

        sum = sum + (termij - (termi + termj))
        
        if (sum .ne. zero) then
          hood(ll) = -dnp2 * log(sum)
        else
          hood(ll) = -FLMAX
        end if

        kk = ll
      end do

      lwork = lwopt

      return
      end
      subroutine likei( cltree, ns, x, n, p, ic, ng, v, new, hood)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which
c are fixed in advance, and clusters have equal (but unknown) volume

      implicit double precision (a-h,o-z)

      integer            cltree(2,*), ns, n, p, ic(n), ng, new(*)

c     double precision   x(n,p), v(p), hood(ns+1)
      double precision   x(n,*), v(*), hood(*)
c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  v       double  (scratch) (p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) The loglikelihood (up to an additive 
c                   constant) for each stage.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        d1mach

c------------------------------------------------------------------------------

      FLMAX  = d1mach(2)

      ns1    = ns + 1
      dnp2   = dble(n*p)/two

c group heads should be first among rows of x

      i = 1
      j = 2
 1        continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2        continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3        continue

      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
c update sum of squares
          k = ic(i)
          if (k .eq. 1) then
            ic(i) = j
            ic(j) = 2
            call dcopy( p, x(j,1), n, v, 1)
            call daxpy( p, (-one), x(i,1), n, v, 1)
            call daxpy( p, one, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
            x(j,1) = ddot( p, v, 1, v, 1) / two
          else
            ic(j) = 0
            ni    = ic(k)
            ic(k) = ni + 1
            ri    = dble(ni)
            rij   = dble(ni+1)
            qj    = one/rij
            qi    = ri*qj
            si    = sqrt(qi)
            sj    = sqrt(qj)
            call dcopy( p, x(j,1), n, v, 1)
            call dscal( p, si, v, 1)
            call daxpy( p, (-sj), x(i,1), n, v, 1)
            x(k,1) = qi*x(k,1) + qj*ddot(p, v, 1, v, 1)
            call dscal( p, si, x(i,1), n)
            call daxpy( p, sj, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
          end if 
        else 
          ic(j) = 1
        end if
      end do

c each group now points to a scaled sum of its elements
c and has its trace stored if it has more than one element

      trw = zero

      do k = 1, ng
        i = ic(k)
        if (i .ne. 1) then
          ni    = ic(i)
          trw   = trw + x(i,1)
          ic(k) = ni
        end if
      end do

c trace is no longer stored

c     sigmsq = trw/dnp

      if (trw .eq. zero) then
        hood(1) = -FLMAX
      else
c       hood(1) = dnp*(one + log(sigmsq))
c       hood(1) = dnp*(one + (log(trw) - dnplog))  ???
        hood(1) = -dnp2*log(trw)
      end if

      k = 1
      do l = 2, ns1

        i = cltree(1,k)
        j = cltree(2,k)

        ni  = ic(i)
        ri  = dble(ni)
         
        nj  = ic(j)
        rj  = dble(nj)

        nij = ni + nj
        rij = dble(nij)
        qij = one/rij
        qi  = ri*qij
        qj  = rj*qij
        si  = sqrt(qi)
        sj  = sqrt(qj)
        call dcopy(p, x(i,1), n, v, 1)
        call dscal( p, sj, v, 1)
        call daxpy( p, (-si), x(j,1), n, v, 1)

        trw = trw + ddot(p,v,1,v,1)
 
        if (l .le. ns) then
c merge groups i and j
          call dscal( p, si, x(i,1), n)
          call daxpy( p, sj, x(j,1), n, x(i,1), n)
        end if

        ic(i)  = nij

        new(k) = nij

c       sigmsq = trw/dnp

        if (trw .eq. zero) then
          hood(l) = -FLMAX
        else
c         hood(l) = dnp*(one + log(sigmsq))
c         hood(l) = dnp*(one + (log(trw) - dnplog))  ???
          hood(l) = -dnp2*log(trw)
        end if

        k = l

      end do

      return
      end
      subroutine likvfv( cltree, ns, x, n, p, ic, ng, Vinv, s, lwork, w,
     *                   v, t, r, uu, vv, new, hood, info)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which
c shape, volume and orientation are allowed to vary between clusters 

      implicit double precision (a-h,o-z)

      integer            cltree(2,ns),ns,n,p,ic(n),ng,new(ns),info

c     double precision   x(n,p),s(p),v(p),w(lwork),t(p,p),r(p,p),hood(ns+1)
      double precision   x(n,*),s(*),v(*),w(*),t(p,*),r(p,*),hood(*)

c     double precision   uu(p,p), vv(p,p)
      double precision   uu(p,*), vv(p,*)

c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  Vinv    double  (input) Estimated reciprocal hypervolume of the data region.
c  s       double  (input) (p) Shape matrix. It's square is the A matrix of
c                   Banfield and Raftery.
c  v       double  (scratch) (p).
c  w       double  (scratch) max(4*p,5*p-4).
c  t,r     double  (scratch) (p*p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) Twice the loglikelihood (up to an additive
c                   constant) for each stage.
c  info    integer (output) return code from the LAPACK svd routine

      integer                 psq, pm1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision       half
      parameter             (half = .5d0)

      double precision       rthalf
      parameter             (rthalf = .7071067811865476d0)

      double precision       d1mach

      common /MCLMCH/        FLMAX
      save   /MCLMCH/  

c------------------------------------------------------------------------------

      do i = 1, p
        s(i) = one/s(i)
      end do

c     lwork = max(4*p,5*p-4)

      lwopt = 0

      FLMAX = d1mach(2)

      dp2   = half*dble(p)

      Vlog = -log(Vinv)

      ns1   = ns + 1

      psq   = p*p
      pm1   = p-1

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     nsing = 0

      sum   = zero

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, v, 1)
          call dscal( p, si, v, 1)
          call daxpy( p, (-sj), x(k,1), n, v, 1)
c         trcij = qi*trcij + qj*ddot( p, v, 1, v, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, v, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c         term(k) = trmij
c         trac(k) = trcij
c         d(ll+k) = trcij
c copy triangular factor into the rows of x
          j  = k
          m  = p
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij) = n+l
          if (l .gt. 2) then
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do i = 1, min(l-1,p)
              call dcopy( m, r(i,i), p, t(i,i), p)
              m = m - 1
            end do
            m = min(p,l-1)
c           call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                    dummy, 1, dummy, 1, w, lwork, info)
            call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do i = 1, m
              v(i) = v(i)*s(i)
            end do
            temp = ddot( m, v, 1, v, 1)
          else 
            temp = dnrm2( p, r, p)*s(1)
            temp = temp*temp
          end if
          if (temp .gt. zero) then
            termij = -dp2*(dble(l)*log(temp/dble(l)))
          else
c           termij = -FLMAX
            termij = -dble(l)*Vlog
c           nsing  = nsing + 1
          end if 
          if (l .gt. 2) x(ic(ic(k)), 1) = termij
        else
          ic(k) = 1
          termij = -Vlog
c         nsing  = nsing + 1
        end if 
        sum = sum + termij
      end do

c     if (nsing .eq. 0) then
c        hood(1) = dp2*sum
         hood(1) = sum
c     else
c       sum     =  zero
c       hood(1) = -FLMAX
c     end if
 
      kk = 1
      do ll = 2, ns1
    
        i   = cltree(1,kk)
        j   = cltree(2,kk)
 
        ni  = 1
        nj  = 1

        icj = ic(j)
        ici = ic(i)

        if (icj .eq. 1) then
c         termj = -FLMAX
          termj = -Vlog
          if (ici .eq. 1) then
c           termi = -FLMAX
            termi = -Vlog
            nij   = 2
            rij   = two
            si    = rthalf
            sj    = rthalf
            sij   = rthalf
            call dcopy( p, x(i,1), n, v, 1)
            call daxpy( p, (-one), x(j,1), n, v, 1)
            call dscal( p, rthalf, v, 1)
            call dcopy( p, v, 1, r(1,1), p)
            temp = dnrm2( p, v, 1)*s(1)
            temp = temp*temp
          else
            m  = p
            l  = ici
 110        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni  = ni + 1
            m   = m - 1
            l   = ic(l)
            if (l .le. n) goto 110
            ni    = l - n
            if (ni .gt. 2) then
              termi = x(ic(ici), 1)
            else
              temp  = dnrm2( p, r, p)*s(1)
              temp  = temp*temp
              if (temp .gt. zero) then
                termi = -dp2*(two*log(temp/two))
              else
                termi = -two*Vlog
              end if
            end if
            ri    = dble(ni)
            nij   = ni + 1
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            si    = sqrt(qi)
            sj    = sqrt(qij)
            sij   = sj
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
            call mclrup( nij, p, v, r, p)
            call dcopy( psq, zero, 0, t, 1)
            m = p
            do k = 1, min(p,nij-1)
              call dcopy( m, r(k,k), p, t(k,k), p)
              m = m - 1
            end do
            m = min(p,nij-1)
c           call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                    dummy, 1, dummy, 1, w, lwork, info)
            call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                    uu, p, vv, p, w, lwork, info)
            if (info .lt. 0) then
              call intpr( 'SVD fails', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            if (info .gt. 0) then
              call intpr( 'SVD does not converge', -1, info, 1)
              nd = info
              return
            end if
            lwopt = max(lwopt,int(w(1)))
            do k = 1, m
              v(k) = v(k)*s(k)
            end do
            temp = ddot( m, v, 1, v, 1)
          end if
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj = l - n
          if (nj .gt. 2) then
            termj = x(ic(icj), 1)
          else
            temp = dnrm2( p, r, p)*s(1)
            temp = temp*temp
            if (temp .gt. zero) then
              termj = -dp2*(two*log(temp/two))
            else
c             termj = -FLMAX
              termj = -two*Vlog
            end if
          end if
          rj  = dble(nj)
          ni  = 1
          ici = ic(i)
          if (ici .eq. 1) then
c           termi = -FLMAX
            termi = -Vlog
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qi  = qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            sij = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          else
            m  = p
            l  = ici
            k  = nj + 1
 130        continue
            call dcopy( m, x(l,ni), n, v, 1)
            call mclrup( k, m, v, r(ni,ni), p)
            ni    = ni + 1
            m     = m - 1
            l     = ic(l)
            if (l .le. n) goto 130
            ni    = l - n
            if (ni .gt. 2) then
              termi = x(ic(ici), 1)
            else
              temp = dnrm2( p, v, 1)*s(1)
              temp = temp*temp
              if (temp .gt. zero) then
                termi = -dp2*(two*log(temp/two))
              else
c               termi = -FLMAX
                termi = -two*Vlog
              end if
            end if
            ri    = dble(ni)
            nij   = ni + nj
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            qj    = rj*qij
            si    = sqrt(qi)
            sj    = sqrt(qj)
            sij   = sqrt(qij)
            call dcopy(p, x(i,1), n, v, 1)
            call dscal( p, sj, v, 1)
            call daxpy( p, (-si), x(j,1), n, v, 1)
          end if
          call mclrup( nij, p, v, r, p)
          call dcopy( psq, zero, 0, t, 1)
          m = p
          do k = 1, min(nij-1,p)
            call dcopy( m, r(k,k), p, t(k,k), p)
            m = m - 1
          end do
          m = min(nij-1,p)
c         call dgesvd( 'N', 'N', m, p, t, p, v, 
c    *                  dummy, 1, dummy, 1, w, lwork, info)
          call dgesvd( 'A', 'A', m, p, t, p, v, 
     *                  uu, p, vv, p, w, lwork, info)
          if (info .lt. 0) then
            call intpr( 'SVD fails', -1, info, 1)
            nd = info
            return
          end if
          lwopt = max(lwopt,int(w(1)))
          if (info .gt. 0) then
            call intpr( 'SVD does not converge', -1, info, 1)
            nd = info
            return
          end if
          do k = 1, m
            v(k) = v(k)*s(k)
          end do
          temp = ddot( m, v, 1, v, 1)
        end if

        if (temp .gt. zero) then
          termij = -dp2*(dble(nij)*log(temp/dble(nij)))
        else
          termij = -dble(nij)*Vlog
        end if

        call dcopy( p, x(i,1), n, v, 1)
        call dscal( p, si, v, 1)
        call daxpy( p, sj, x(j,1), n, v, 1)

        call dcopy( p, r(1,1), p, x(j,1), n)

        nopt = nij

        if (ni .eq. 1) then

          if (nj .eq. 1) then
            ic(j)  = n+2
          else
            l   = icj
            m   = pm1
            nij = 2

 210        continue
            
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij = nij + 1
              m   = m   - 1
              k   = l
              l   = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

          end if

        else

          l   = ic(i)
          m   = pm1
          nij = 2

 220      continue

            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

            if (nij .le. p .and. nj .ne. 1) then

              l     = ic(j)
              ic(k) = l

 230          continue

                call dcopy( m, r(nij,nij), p, x(l,nij), n)
                nij = nij + 1
                m   = m   - 1
                k   = l
                l   = ic(l)
                if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230

            end if

          ic(j) = ic(i)
          ic(k) = n + nopt

        end if

        ic(i) = j

        if (nopt .gt. 2) then
          m      = ic(j)
          x(m,1) = termij
        endif

        call dcopy( p, v, 1, x(i,1), n)

        new(kk) = nopt

c       if (termi .eq. (-FLMAX)) then
c         nsing = nsing - 1
c       else
c         sum = sum - termi
c       end if

c       if (termj .eq. (-FLMAX)) then
c         nsing = nsing - 1
c       else
c         sum = sum - termj
c       end if

c       if (termij .eq. (-FLMAX)) then
c         nsing = nsing + 1
c       else
c         sum = sum + termij
c       end if

        sum = sum + (termij - (termi + termj))
        
c       if (nsing .eq. 0) then
c         hood(ll) = dp2*sum
          hood(ll) = sum
c       else
c         hood(ll) = -FLMAX
c       end if

c       call intpr( 'll', -1, ll, 1)
c       call intpr( 'nsing', -1, nsing, 1)
c       call dblepr( 'hood', -1, hood(ll), 1)

        kk = ll
      end do

      lwork = lwopt

      return
      end
      subroutine likvi( cltree, ns, x, n, p, ic, ng, Vinv, w, new, hood)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which 
c shape and orientation are fixed in advance, and clusters have unequal volumes

      implicit double precision (a-h,o-z)

      integer            cltree(2,ns), ns, n, p, ic(n), ng, new(*)

c     double precision   x(n,p), w(p), hood(ns+1)
      double precision   x(n,*), w(*), hood(*)
c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  Vinv    double  (input) Estimated reciprocal hypervolume of the data region.
c  w       double  (scratch) (p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) The loglikelihood (up to an additive 
c                   constant) for each stage.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        d1mach

c------------------------------------------------------------------------------

      FLMAX = d1mach(2)

      dp2   = dble(p)/two

      Vlog  = -log(Vinv)

      ns1   = ns + 1

c group heads should be first among rows of x

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

c     nsing = 0
     
      do j = 1, n
        i = ic(j)
        if (i .ne. j) then
c update sum of squares
          k = ic(i)
          if (k .eq. 1) then
            ic(i) = j
            ic(j) = 2
            call dcopy( p, x(j,1), n, w, 1)
            call daxpy( p, (-one), x(i,1), n, w, 1)
            call daxpy( p, one, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
            x(j,1) = ddot( p, w, 1, w, 1) / two
          else
            ic(j) = 0
            ni    = ic(k)
            ic(k) = ni + 1
            ri    = dble(ni)
            rij   = dble(ni+1)
            qj    = one/rij
            qi    = ri*qj
            si    = sqrt(qi)
            sj    = sqrt(qj)
            call dcopy( p, x(j,1), n, w, 1)
            call dscal( p, si, w, 1)
            call daxpy( p, (-sj), x(i,1), n, w, 1)
            x(k,1) = qi*x(k,1) + qj*ddot(p, w, 1, w, 1)
            call dscal( p, si, x(i,1), n)
            call daxpy( p, sj, x(j,1), n, x(i,1), n)
c           call dcopy( p, FLMAX, 0, x(j,1), n)
          end if 
        else 
          ic(j) = 1
c         nsing = nsing + 1
        end if
      end do
       
c store terms also so as not to recompute them

      sum = zero 

      do k = 1, ng
        i = ic(k)
        if (i .ne. 1) then
          ni    = ic(i)
          traci = x(i,1)
          ri     = dble(ni)
          if (traci .ne. zero) then
            termi  = -dp2*(ri*log(traci/ri))
          else
c           x(i,2) = -FLMAX
            termi  = -ri*Vlog
c           nsing  =  nsing + 1
          end if
          x(i,2) = termi
        else
          ni    = 1
          traci = zero
          termi =  - Vlog
        end if
        sum = sum + termi
      end do

c     if (nsing .ne. 0) then

c       hood(1) = -FLMAX
c       sum     =  zero

c     else

c       do k = 1, ng
c         i    = ic(k)
c         sum  = sum + x(i,2)
c       end do

c       hood(1) = dp*(dn + sum) ???
        hood(1) = sum

c     end if

      k = 1
      do l = 2, ns1

        i = cltree(1,k)
        j = cltree(2,k)

        ni = ic(i)
        if (ni .eq. 1) then
          traci = zero
c         termi = -FLMAX
          termi = -Vlog
          ri    = one
        else 
          traci = x(ni,1)
          termi = x(ni,2)
          ni    = ic(ni)
          ri    = dble(ni)
        end if
         
        nj = ic(j)
        if (nj .eq. 1) then
          tracj = zero
c         termj = -FLMAX
          termj = -Vlog
          rj    = one
        else
          tracj = x(nj,1)
          termj = x(nj,2)
          nj    = ic(nj)
          rj    = dble(nj)
        end if

        nij = ni + nj
        rij = dble(nij)
        qij = one/rij
        qi  = ri*qij
        qj  = rj*qij
        si  = sqrt(qi)
        sj  = sqrt(qj)
        call dcopy(p, x(i,1), n, w, 1)
        call dscal( p, sj, w, 1)
        call daxpy( p, (-si), x(j,1), n, w, 1)

        tracij = (traci + tracj) + ddot(p,w,1,w,1)

        if (tracij .ne. zero) then
          termij = -dp2*(rij*log(tracij/rij))
        else
          termij = -rij*Vlog
        end if

c merge groups i and j
        call dscal( p, si, x(i,1), n)
        call daxpy( p, sj, x(j,1), n, x(i,1), n)

c       if (traci .eq. zero) then
c         nsing = nsing - 1
c       else
c         sum = sum - termi
c       end if

c       if (tracj .eq. zero) then
c         nsing = nsing - 1
c       else
c         sum = sum - termj
c       end if

c       if (tracij .eq. zero) then
c         nsing = nsing + 1
c       else
c         sum = sum + termij
c       end if

        sum = sum + (termij - (termi + termj))

        ic(i)  = j
        ic(j)  = nij
        x(j,1) = tracij
        x(j,2) = termij

        new(k)  = nij

c       if (nsing .ne. 0 .or. nbad .ne. 0) then
c         hood(l) = -FLMAX
c       else
c         hood(l) = dp*(dn + sum)  ???
          hood(l) = sum
c       end if

        k = l
      end do

      return
      end


















      subroutine likvvv( cltree, ns, x, n, p, ic, ng, Vinv, 
     *                   w, r, new, hood)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N-00014-88-K-0265 and N-00014-91-J-1074
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c likelihood and cluster size for Gaussian hierarchical clustering in which
c shape, volume and orientation are allowed to vary between clusters 

      implicit double precision (a-h,o-z)

      integer            cltree(2,ns), ns, n, p, ic(n), ng, new(ns)

c     double precision   x(n,p), w(p), r(p,p), hood(ns+1)
      double precision   x(n,*), w(*), r(p,*), hood(*)

c------------------------------------------------------------------------------
c
c  cltree  integer (input) (2 by ns) The classification tree from mclust.
c  ns      integer (input) number of columns in cltree
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, the first two columns
c                   and ns rows contain the merge indices.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  ic      integer (input) (n) Initial partitioning of the data; groups must
c                   be numbered consecutively.
c  ng      integer (input) Number of groups in initial partition.
c  Vinv    double  (input) Estimated reciprocal hypervolume of the data region.
c  w       double  (scratch) (p).
c  r       double  (scratch) (p*p).
c  new     integer (output) (ns) The size of the new cluster formed at each
c                   stage.
c  hood    double  (output) (ns+1) Twice the loglikelihood (up to an additive
c                   constant) for each stage.

      integer                 psq, pm1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision       rthalf
      parameter             (rthalf = .7071067811865476d0)

      double precision       d1mach

      common /MCLMCH/        FLMAX
      save   /MCLMCH/  

c------------------------------------------------------------------------------

      FLMAX = d1mach(2)

      dp    = dble(p)

      Vlog  = -log(Vinv)

      ns1   = ns + 1

      psq   = p*p
      pm1   = p-1

      i = 1
      j = 2
 1    continue
        icj = ic(j)
        if (icj .ne. j) goto 2
        if (j .eq. ng)  goto 3
        i = j
        j = j + 1
      goto 1

 2    continue

      k = i
      m = j + 1
      do j = m, n
        icj = ic(j)
        if (icj .gt. k) then
          k = k + 1
          call dswap( p, x(k,1), n, x(j,1), n)
          ic(j) = ic(k)
          ic(k) = icj
        end if
      end do

 3    continue

c set up pointers

      if (ng .eq. n) goto 4

      do j = n, ng+1, -1
        icj = ic(j)
        i   = ic(icj)
        ic(icj) = j
        if (i .ne. icj) then
          ic(j) = i
        else
          ic(j) = j
        end if
      end do

 4    continue

c     nsing = 0

      sum   = zero

c initialize by simulating merges       

      do k = 1, ng
        j = ic(k)
        if (j .ne. k) then
c non-singleton
          call dcopy( psq, zero, 0, r, 1)
          trcij = zero
          l     = 1
 10       continue
          m     = l + 1
          qj    = one/dble(m)
          qi    = dble(l)*qj
          si    = sqrt(qi)
          sj    = sqrt(qj)
          call dcopy( p, x(j,1), n, w, 1)
          call dscal( p, si, w, 1)
          call daxpy( p, (-sj), x(k,1), n, w, 1)
c         trcij = qi*trcij + qj*ddot( p, w, 1, w, 1)
          call dscal( p, si, x(k,1), n)
          call daxpy( p, sj, x(j,1), n, x(k,1), n)
          call mclrup( m, p, w, r, p)
          l = m
          i = ic(j)
          if (i .eq. j) goto 20
          j = i
          goto 10
 20       continue
c         term(k) = trmij
c         trac(k) = trcij
c         d(ll+k) = trcij
c copy triangular factor into the rows of x
          j  = k
          m  = p 
          do i = 1, min(l-1,p)
            j  = ic(j)
            call dcopy( m, r(i,i), p, x(j,i), n)
            m  = m - 1      
          end do
          ij = j
          if (l .ge. p) then
            do m = p, l
              icj   = ic(j)
              ic(j) = -k
              j     = icj
            end do
          end if
          ic(ij) = n+l
          if (l .gt. p) then
            detlog = detmc2( p, r)
            if (detlog .ne. (-FLMAX)) then
              termij = -dble(l)*(detlog - dp*log(dble(l)))/two
            else
c             termij = -FLMAX
              termij = -dble(l)*Vlog
c             nsing  = nsing + 1
            end if
            x(ic(ic(k)), 1) = termij
          else
c           termij = -FLMAX
            termij = -dble(l)*Vlog
c           nsing  = nsing + 1
          end if
        else
          ic(k)  = 1
          termij = -Vlog
c         nsing = nsing + 1
        end if
        sum = sum + termij
      end do

c     if (nsing .eq. 0) then
c       hood(1) = dnp + sum  ???
        hood(1) = sum
c     else
c       sum     = 0
c       hood(1) = -FLMAX
c     end if

      kk = 1
      do ll = 2, ns1

        i   = cltree(1,kk)
        j   = cltree(2,kk)

        ni  = 1
        nj  = 1

        icj = ic(j)
        ici = ic(i)

        if (icj .eq. 1) then
c         termj = -FLMAX
          termj = -Vlog
          if (ici .eq. 1) then
c           termi = -FLMAX
            termi = -Vlog
            nij   = 2
            rij   = two
            si    = rthalf
            sj    = rthalf
            sij   = rthalf
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), x(j,1), n, w, 1)
            call dscal( p, rthalf, w, 1)
            call mclrup( nij, p, w, r, p)
c           trcij = half*ddot( p, w, 1, w, 1)
c           termij = -FLMAX
            termij = -two*Vlog
          else
            m  = p
            l  = ici
 110        continue
            call dcopy( m, x(l,ni), n, r(ni,ni), p)
            ni  = ni + 1
            m   = m - 1
            l   = ic(l)
            if (l .le. n) goto 110
            ni  = l - n
            if (ni .le. p) then
c             termi = -FLMAX
              termi = -dble(ni)*Vlog
            else
              termi = x(ic(ici), 1)
            end if
            ri  = dble(ni)
            nij = ni + 1
            rij = dble(nij)
            qij = one/rij
            qi  = ri*qij
            si  = sqrt(qi)
            sj  = sqrt(qij)
            sij = sj
            call dcopy(p, x(i,1), n, w, 1)
            call dscal( p, sj, w, 1)
            call daxpy( p, (-si), x(j,1), n, w, 1)
c           trcij = qi*traci + qij*ddot(p,w,1,w,1)
            call mclrup( nij, p, w, r, p)
            if (nij .le. p) then
c             termij = -FLMAX
              termij = -rij*Vlog
            else
              detlog = detmc2( p, r)
              if (detlog .ne. (-FLMAX)) then
                termij = -rij*(detlog - dp*log(rij))/two
              else
c               termij = -FLMAX
                termij = -rij*Vlog
              endif
            end if
          end if
        else
          m = p
          l = icj
 120      continue
          call dcopy( m, x(l,nj), n, r(nj,nj), p)
          nj = nj + 1
          m  = m - 1
          l  = ic(l)
          if (l .le. n) goto 120
          nj = l - n
          if (nj .le. p) then
c           termj = -FLMAX
            termj = -dble(nj)*Vlog
          else
            termj = x(ic(icj), 1)
          end if
          rj  = dble(nj)
          ni  = 1
          ici = ic(i)
          if (ici .eq. 1) then
c           termi = -FLMAX
            termi = -Vlog
            nij = nj + 1
            rij = dble(nij)
            qij = one/rij
            qi  = qij
            qj  = rj*qij
            si  = sqrt(qi)
            sj  = sqrt(qj)
            sij = sqrt(qij)
            call dcopy(p, x(i,1), n, w, 1)
            call dscal( p, sj, w, 1)
            call daxpy( p, (-si), x(j,1), n, w, 1)
c           trcij = qj*tracj + qij*ddot(p,w,1,w,1)
          else
            m  = p
            l  = ici
            k  = nj + 1
 130        continue
            call dcopy( m, x(l,ni), n, w, 1)
            call mclrup( k, m, w, r(ni,ni), p)
            ni    = ni + 1
            m     = m - 1
            l     = ic(l)
            if (l .le. n) goto 130
            ni    = l - n
            if (ni .le. p) then
c             termi = -FLMAX
              termi = -dble(ni)*Vlog
            else
              termi = x(ic(ici), 1)
            end if
            ri    = dble(ni)
            nij   = ni + nj
            rij   = dble(nij)
            qij   = one/rij
            qi    = ri*qij
            qj    = rj*qij
            si    = sqrt(qi)
            sj    = sqrt(qj)
            sij   = sqrt(qij)
            call dcopy(p, x(i,1), n, w, 1)
            call dscal( p, sj, w, 1)
            call daxpy( p, (-si), x(j,1), n, w, 1)
c           trcij = (qi*traci + qj*tracj) + qij*ddot(p,w,1,w,1)
          end if
          call mclrup( nij, p, w, r, p)
          if (nij .le. p) then
c           termij = -FLMAX
            termij = -rij*Vlog
          else
            detlog = detmc2( p, r)
            if (detlog .ne. (-FLMAX)) then
              termij = -rij*(detlog - dp*log(rij))/two
            else
c             termij = -FLMAX
              termij = -rij*Vlog
            endif
          end if
        end if

        call dcopy( p, x(i,1), n, w, 1)
        call dscal( p, si, w, 1)
        call daxpy( p, sj, x(j,1), n, w, 1)

        call dcopy( p, r(1,1), p, x(j,1), n)

        nopt = nij

        if (ni .eq. 1) then

          if (nj .eq. 1) then
            ic(j)  = n+2
          else
            l   = icj
            m   = pm1
            nij = 2

 210        continue
            
              call dcopy( m, r(nij,nij), p, x(l,nij), n)
              nij = nij + 1
              m   = m   - 1
              k   = l
              l   = ic(l)
              if (l .le. n .and. nij .le. min(nopt-1,p)) goto 210

            ic(k) = n + nopt

          end if

        else

          l   = ic(i)
          m   = pm1
          nij = 2

 220      continue

            call dcopy( m, r(nij,nij), p, x(l,nij), n)
            nij = nij + 1
            m   = m   - 1
            k   = l
            l   = ic(l)
            if (l .le. n .and. nij .le. min(nopt-1,p)) goto 220

            if (nij .le. p .and. nj .ne. 1) then

              l     = ic(j)
              ic(k) = l

 230          continue

                call dcopy( m, r(nij,nij), p, x(l,nij), n)
                nij = nij + 1
                m   = m   - 1
                k   = l
                l   = ic(l)
                if (l .le. n .and. nij .le. min(nopt-1,p)) goto 230

            end if

          ic(j) = ic(i)
          ic(k) = n + nopt

        end if

        ic(i) = j

        if (nopt .gt. 2) then
          m      = ic(j)
          x(m,1) = termij
        endif

        call dcopy( p, w, 1, x(i,1), n)

        new(kk) = nopt

c       if (termi .eq. (-FLMAX)) then
c         nsing = nsing - 1
c       else
c         sum = sum - termi
c       end if

c       if (termj .eq. (-FLMAX)) then
c         nsing = nsing - 1
c       else
c         sum = sum - termj
c       end if

c       if (termij .eq. (-FLMAX)) then
c         nsing = nsing + 1
c       else
c         sum = sum + termij
c       end if

       sum = sum + (termij - (termi + termj))

c      if (nsing .eq. 0) then
c         hood(ll) = dnp + sum
          hood(ll) = sum
c      else
c        hood(ll) = -FLMAX
c      end if

        kk = ll
      end do

      return
      end
      subroutine mclcpy( i, n, d)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Copies row n into row i.
c Note : i > 1, i < n, n > 2 required.

      implicit double precision (a-h,o-z)

      double precision d(*)

      FLMAX = d1mach(2)

      i1 = i - 1
      ii = (i1*(i1-1))/2 + 1

      n1 = n - 1
      nn = (n1*(n1-1))/2 + 1

c     if (i .gt. 1) then
        call dcopy( i1, d(nn), 1, d(ii), 1)
c       call dcopy( i1, FLMAX, 0, d(nn), 1)
        ii = ii + i1 + i1
        nn = nn + i
c     end if

      k = i

 100  continue

        d(ii) = d(nn)
c       d(nn) = FLMAX

        ii = ii + k
 
        nn = nn + 1
        
        k  = k  + 1

        if (k .lt. n1) goto 100

c     d(nn) = FLMAX

      return
      end 
      subroutine mclrup( l, n, v, r, lr)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c rank one row update of n by n upper triangular factor r

      implicit double precision (a-h,o-z)

c     double precision v(n), r(lr,n)
      double precision v(*), r(lr,*)

      if (l .eq. 1) return

      k = l - 1
      if (k .le. n) then

        call dcopy( n, v, 1, r(k,1), lr)

        if (k .eq. 1) return

        if (n .gt. 1) then
          i = 1
          m = n
          do j = 2, k
            call drotg( r(i,i), r(k,i), cs, sn)
            m = m - 1
            call drot( m, r(i,j), lr, r(k,j), lr, cs, sn)
            i = j
          end do
        else
          call drotg( r(1,1), r(k,1), cs, sn)
        end if
  
      else

        if (n .gt. 1) then
          i = 1
          m = n
          do j = 2, n
            call drotg( r(i,i), v(i), cs, sn)
            m = m - 1
            call drot( m, r(i,j), lr, v(j), 1, cs, sn)
            i = j
          end do
        end if

        call drotg( r(n,n), v(n), cs, sn)

      end if

      return  
      end 
      subroutine mcltrw( x, n, p, u, ss)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Computes the trace of the sample cross product matrix.

      implicit double precision (a-h,o-z)

c     integer            n, p
      integer            n, p

c     double precision   x(n,p), u(p)
      double precision   x(n,*), u(*)
c------------------------------------------------------------------------------
c
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, x is overwritten.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  u       output  (scratch) (p) 
c  ss output

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

c form mean
      ss = one / dble(n)
      call dcopy( p, zero, 0, u, 1)
      do i = 1, n
        call daxpy( p, ss, x(i,1), n, u, 1)
      end do

c subtract mean and form sum of squares
      ss = zero
      do j = 1, p
        call daxpy( n, (-one), u(j), 0, x(1,j), 1)
        ss = ss + ddot(n, x(1,j), 1, x(1,j), 1)
      end do

      return
      end
      subroutine mclvol( x, n, p, u, v, w,
     *                   work, lwork, iwork, liwork, 
     *                   info)

c copyright 1996 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c Computes quantities whose product is an approximation to the hypervolume for
c of a data region via principal components.

      implicit double precision (a-h,o-z)

c     integer            n, p, iwork(liwork)
      integer            n, p, iwork(*)

c     double precision   x(n,p), u(p), v(p,p), w(p,p), work(lwork),
      double precision   x(n,*), u(*), v(p,*), w(p,p), work(*)
c------------------------------------------------------------------------------
c
c  x       double  (input/output) On input, the (n by p) matrix containing
c                   the observations. On output, x is overwritten.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  u       output  (scratch) (p) Positive quantities whose product is an
c                   approximate volume of the data region.
c  w       double  (scratch) (p*p)
c  v       double  (scratch) (p*p)
c  work    double  (scratch) (lwork)
c  lwork   integer
c  iwork   integer  (scratch) (liwork)
c  liwork  integer

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        d1mach

c------------------------------------------------------------------------------

c form mean
      s = one / dble(n)
      call dcopy( p, zero, 0, u, 1)
      do i = 1, n
        call daxpy( p, s, x(i,1), n, u, 1)
      end do

c subtract mean
      do j = 1, p
        call daxpy( n, (-one), u(j), 0, x(1,j), 1)
      end do


      if (.false.) then
c this gets the eigenvectors but x is overwritten

c get right singular vectors
        call dgesvd( 'N', 'A', n, p, x, n, u, 
     *                dummy, 1, w, p, work, lwork, info)

        if (info .lt. 0) return

        if (info .eq. 0) then
          lwork = int(work(1))
          do i = 1, p
            v(i,i) = w(i,i)
            if (i .gt. 1) then
              do j = 1, (i-1)
                v(i,j) = w(j,i)
                v(j,i) = w(i,j)
              end do
            end if
          end do
          goto 100
        end if

      end if

c form crossproduct
      call dsyrk( 'U', 'T', p, n, one, x, n, zero, w, p)

c get eigenvectors
 
      do j = 1, p
        do i = 1, j
          v(i,j) = w(i,j)
        end do
      end do

      call dsyevd( 'V', 'U', p, v, p, u, 
     *              work, lwork, iwork, liwork, info)

      if (info .lt. 0) return

      if (info .eq. 0) then
        lwork  = int(work(1))
        liwork = iwork(1)
        goto 100
      end if

      EPSMAX = d1mach(4)

      call dsyevx( 'V', 'A', 'U', p, w, p, dummy, dummy, i, i,
     *              sqrt(EPSMAX), j, u, v, p,
     *              work, lwork, iwork(p+1), iwork, info)
          
      if (info .ne. 0) return

      lwork  = int(work(1))
      liwork = -1 

 100  continue
  
      FLMAX = d1mach(2)

c form xv

c     vol = one
      do j = 1, p
        call dgemv( 'N', n, p, one, x, n, v(1,j), 1, zero, work, 1)
        cmax = -FLMAX
        cmin =  FLMAX
        do i = 1, n
          temp = work(i)
          if (temp .gt. cmax) cmax = temp
          if (temp .lt. cmin) cmin = temp
        end do
        u(j) = cmax - cmin
c       vol  = vol * (cmax - cmin)
      end do

      return

      end
      subroutine meeee ( x, z, n, p, G, eps, tol, maxi, mu, U, prob, w)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for constant-variance Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G),mu(p,G),U(p,p),prob(G),w(p)
      double precision   x(n,*), z(n,*), mu(p,*), U(p,*), prob(*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower bound on the reciprocal
c                   condition estimate of the covariance. On output,
c                   reciprocal condition estimate at the last iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for 
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  mu      double  (scratch) (p,G)
c  U       double  (scratch) (p,p)
c  prob    double  (scratch) (G) (not needed with equal proportions)
c  w       double  (scratch) (p)

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      sclfac = one/sqrt(dble(n))

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX

      iter   = 1

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

110   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
c       call dblepr( "U(*j)", -1, U(1,j), j)
        call dscal( j, sclfac, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .le. eps) goto 900

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G    = -G

c  probability = one/dble(G)

160   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .le. eps) goto 900

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

210   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
c       call dblepr( "U(*j)", -1, U(1,j), j)
        call dscal( j, sclfac, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, U, p1)
      end if

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G    = -G

c  probability = one/dble(G)

260   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, U, p1)
      end if

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = rcmin*rcmin
      maxi = iter

      return
      end
      subroutine meeev ( x, z, n, p, G, eps, tol, maxi, 
     *                   mu, U, prob, shape, s, w, lwork)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G), mu(p,G), U(p,p,G), prob(G),
c    *                   shape(p), s(p), w(lwork)
      double precision   x(n,*), z(n,*), mu(p,*), U(p,p,*), prob(*),
     *                   shape(*), s(*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower bound on the volume of
c                   the clusters and the reciprocal condition estimate of the 
c                   square root of the shape. On output, these quantities at
c                   the last iteration. Unlike the other methods, there is
c                   no provision to continue iterating for negative eps.
c  tol     double  (input/output) On input, tolerance on convergence of the
c                   loglikelihood. On output, maximum relative error for ll.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations.
c  mu      double  (scratch) (p,G)
c  U       double  (scratch) (p,p,G)
c  prob    double  (scratch) (G) (not needed with equal proportions)
c  shape   double  (scratch) (p) 
c  s       double  (scratch) (p) 
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4) workspace for LAPACK SVD.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      eps1   = eps(1)
      eps2   = sqrt(eps(2))

      tol    = max(tol,zero)

      dnp    = dble(n*p)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX
      vlamin = FLMAX

      iter   = 1

      if (G .lt. 0) goto 150

110   continue

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          vlam   = dble(sign(1,info))*FLMAX
          vlamin = dble(sign(1,info))*FLMAX
          goto 900
        end if
        do j = 1, p
          temp     = s(j)
          shape(j) = shape(j) + temp*temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        vlam   = zero
        vlamin = zero
        rc     = zero
        rcmin  = zero
        goto 900
      else 
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/dble(n)
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) goto 900

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = prob(k) * exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G    = -G

c probability = one/dble(G)

160   continue

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          vlam   = dble(sign(1,info))*FLMAX
          vlamin = dble(sign(1,info))*FLMAX
          goto 900
        end if
        do j = 1, p
          temp     = s(j)
          shape(j) = shape(j) + temp*temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        vlam   = zero
        vlamin = zero
        rc     = zero
        rcmin  = zero
        goto 900
      else 
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/dble(n)
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) goto 900

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

900   continue

      tol    = err
      eps(1) = vlam
      eps(2) = rc*rc
      maxi   = iter

      return
      end
      subroutine meei  ( x, z, n, p, G, eps, tol, maxi, y, mu, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for spherical, constant-volume Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G), y(n,G), mu(p), prob(G)
      double precision   x(n,*), z(n,*), y(n,*), mu(*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on sigsq.
c                   On output, value of sigsq at the last iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for 
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch) (n,G) (1,1 if eps .lt. 0)
c  mu      double  (scratch) (p)
c  prob    double  (scratch) (G) (not needed with equal proportions

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      dnp    = dble(n*p)

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      FLMAX  = d1mach(2)
      h      = FLMAX/two

      sigmin = FLMAX
      iter   = 1

      if (eps .le. zero) goto 200

      if (G .lt. 0) goto 150

110   continue

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          y(i,k) = sum
        end do
      end do

      sigsq  = sigsq / dnp

      sigmin = min(sigsq,sigmin)

      if (sigsq .le. eps)  goto 900

      const  = dble(p)*(pi2log+log(sigsq))

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(y(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      h    = hood
      iter = iter + 1

      goto 110

150   continue

      G    = -G

c probability = one/dble(G)

160   continue

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          y(i,k) = sum
        end do
      end do

      sigsq = sigsq / dnp

      sigmin = min(sigsq,sigmin)

      if (sigsq .le. eps) goto 900

      const  = dble(p)*(pi2log+log(sigsq))

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(y(i,k)/sigsq))/two
c         temp   = prob*exp(-temp)
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))
 
      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      h    = hood
      iter = iter + 1

      goto 160

200   continue

      if (G .lt. 0) goto 250

210   continue

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / dnp

      sigmin = min(sigsq,sigmin)

      if (sigsq .le. eps) sigsq = eps

      const  = dble(p)*(pi2log+log(sigsq))

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      h    = hood
      iter = iter + 1

      goto 210

250   continue

      G = -G

c  probability = one/dble(G)

260   continue

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / dnp

      sigmin = min(sigsq,sigmin)

      if (sigsq .le. eps) sigsq = eps

      const  = dble(p)*(pi2log+log(sigsq))

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))
 
      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      h    = hood
      iter = iter + 1

      goto 260

900   continue

      tol  = err
      eps  = sigsq
      maxi = iter

      return
      end
      subroutine meneee( x, z, n, p, G, eps, tol, maxi, mu, U, prob, 
     *                   w, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for constant-variance Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G+1),mu(p,G),U(p,p),prob(G),w(p)
      double precision   x(n,*), z(n,*), mu(p,*), U(p,*), prob(*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower bound on reciprocal 
c                   condition estimate of the covariance. On output,
c                   reciprocal condition estimate at the last iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  mu      double  (scratch) (p)
c  U       double  (scratch) (p,p)
c  prob    double  (scratch) (G)
c  w       double  (scratch) (p)
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX

      iter   = 1

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      G1 = G + 1

110   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
c       call dblepr( "U(*j)", -1, U(1,j), j)
        call dscal( j, temp, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .le. eps) goto 900

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

160   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .le. eps) goto 900

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue
 
      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      G1 = G + 1

210   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
c       call dblepr( "U(*j)", -1, U(1,j), j)
        call dscal( j, temp, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, U, p1)
      end if

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

260   continue

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, U(1,j), 1)
      end do

c condition number

      call drnge( p, U, p1, umin, umax)

      rc    = umin/(one+umax)
      rcmin = min(rcmin,rc)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, U, p1)
      end if

      detlog = log(abs(U(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = rc*rc
      maxi = iter

      return
      end
      subroutine meneev( x, z, n, p, G, eps, tol, maxi, 
     *                   mu, U, prob, shape, s, w, lwork, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G+1), mu(p,G), U(p,p,G), prob(G+1),
c    *                   shape(p), s(p), w(lwork)
      double precision   x(n,*), z(n,*), mu(p,*), U(p,p,*), prob(*),
     *                   shape(*), s(*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape.
c                   On output, minimum values of these over the iterations.
c  tol     double  (input/output) On input, tolerance on convergence of the
c                   loglikelihood. On output, maximum relative error for ll.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations.
c  mu      double  (scratch) (p,G)
c  U       double  (scratch) (p,p,G)
c  prob    double  (scratch) (G+1) (not needed with equal proportions)
c  shape   double  (scratch) (p) 
c  s       double  (scratch) (p) 
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4) workspace for LAPACK SVD.
c  Vinv    double  (input) approximate hypervolume of the data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      eps1   = max(eps(1),zero)
      eps2   = sqrt(max(eps(2),zero))

      tol    = max(tol,zero)

      dnp    = dble(n*p)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX
      vlamin = FLMAX

      iter   = 1

      if (G .lt. 0) goto 500

      G1 = G + 1

100   continue

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          vlamin = dble(sign(1,info))*FLMAX
          goto 900
        end if
        do j = 1, p
          temp     = s(j)
          shape(j) = shape(j) + temp*temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        vlam   = zero
        vlamin = zero
        rc     = zero
        rcmin  = zero
        goto 900
      else 
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/sumz
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) goto 900

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = prob(k) * exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 100

500   continue

      G  = -G

      G1 = G + 1

c probability = one/dble(G)

600   continue

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, s, 
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          vlamin = dble(sign(1,info))*FLMAX
          goto 900
        end if
        do j = 1, p
          temp     = s(j)
          shape(j) = shape(j) + temp*temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        vlam   = zero
        vlamin = zero
        rc     = zero
        rcmin  = zero
        goto 900
      else 
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/sumz
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) goto 900

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 600

900   continue

      tol    = err
      eps(1) = vlamin
      eps(2) = rcmin*rcmin
      maxi   = iter

      return
      end
      subroutine menei ( x, z, n, p, G, eps, tol, maxi, 
     *                   y, mu, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for spherical, constant-volume Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G+1), y(n,G), mu(p), prob(G)
      double precision   x(n,*), z(n,*), y(n,*), mu(*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on sigsq. 
c                   On output, value of sigsq at the final iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch)  (n,G) (1,1 only if eps .le. 0)
c  mu      double  (scratch) (p)
c  prob    double  (scratch) (G) (not needed with equal proportions)
c  Vinv    double  (input) estimated reciprocal hypervolume of data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)
      h      = FLMAX/two

      sigmin = FLMAX

      iter   = 1

      if (eps .lt. zero) goto 200

      if (G .lt. 0) goto 150

      G1 = G + 1

110   continue

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob(k) = sum/dble(n)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          y(i,k) = sum
        end do
      end do

      sigsq  = sigsq / (dble(p)*sumz)

      sigmin = min(sigmin,sigsq)

      if (sigsq .le. eps)  goto 900

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = (const+(y(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

160   continue

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          y(i,k) = sum
        end do
      end do

      sigsq = sigsq / (dble(p)*sumz)

      sigmin = min(sigmin,sigsq)

      if (sigsq .le. eps) goto 900

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = (const+(y(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      eps = -eps

      if (G .lt. 0) goto 250

      G1 = G + 1

210   continue

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob(k) = sum/dble(n)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / (dble(p)*sumz)

      sigmin = min(sigmin,sigsq)

      if (sigsq .le. eps) sigsq = eps

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

260   continue

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / (dble(p)*sumz)

      sigmin = min(sigmin,sigsq)

      if (sigsq .le. eps) sigsq = eps

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = sigsq
      maxi = iter

      return
      end
      subroutine menvev( x, z, n, p, G, eps, tol, maxi, 
     *                   mu, U, prob, vlam, shape, s, v, w, lwork, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      integer            maxi(2)
      double precision   eps(2), tol(2)

c     double precision   x(n,p), z(n,G+1), mu(p,G), U(p,p,G), prob(G+1),
c    *                   vlam(G),shape(p), s(p), v(p,G), w(max(lwork,G))
      double precision   x(n,*), z(n,*), mu(p,*), U(p,p,*), prob(*),
     *                   vlam(*), shape(*), s(*), v(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape.
c                   On output, minimum values of these over the iterations.
c  tol     double  (input/output) (2) On input, tolerance on convergence of
c                   the loglikelihood and on the inner iterations for lambda
c                   and shape. On output, maximum relative error for these.
c  maxi    integer (input/output) (2) On input, upper limit on outer/inner
c                   iterations. On output, number of outer/inner iterations.
c  mu      double  (scratch) (p,G)
c  U       double  (scratch) (p,p,G)
c  prob    double  (scratch) (G+1) (needed with equal proportions too)
c  vlam    double  (scratch) (G)
c  shape   double  (scratch) (p)
c  s       double  (scratch) (p) 
c  v       double  (scratch) (p,G)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4,G) workspace for LAPACK
c                   SVD and other things
c  Vinv    double  (input) approximate hypervolume of data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

      if (maxi1 .le. 0) return

      eps1   = max(eps(1),zero)
      eps2   = max(eps(2),zero)

      tol1   = max(tol(1),zero)
      tol2   = max(tol(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX
      vlamin = FLMAX

      errin  = FLMAX
      errim  = zero

      inmax  = 0

      iter   = 1

      if (G .lt. 0) goto 500

      G1 = G + 1

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        do j = 1, p
          temp     = v(j,k)
          temp     = temp*temp
          shape(j) = shape(j) + temp
          v(j,k)   = temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc     = zero
        rcmin  = zero
        ulam   = zero
        vlamin = zero
        goto 900
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))

      ulam = temp / dble(sumz)

      vlamin = ulam

      if (ulam .le. eps1) goto 900

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      call dcopy( G, ulam, 0, vlam, 1)

100   continue

c prob now contains n*prob

      if (maxi2 .le. 0) goto 130

      inner = 0

110   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        if (ulam .le. eps1) goto 900

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          rc    = zero
          rcmin = zero
          goto 900
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = one/exp(sum/dble(p))

        call dscal( p, temp, shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .le. eps2) goto 900

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol2) goto 120
        if (inner .lt. maxi2) goto 110

120   continue
     
      inmax = max(inner,inmax)
      errim = max(errin,errim)

130   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      call dscal( G, one/dble(n), prob, 1)

      do k = 1, G
        temp   = vlam(k)
        vlamin = min(temp,vlamin)
        if (vlamin .le. eps1) goto 900
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = prob(k)*exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      hood  = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol1 ) goto 900
      if (iter .ge. maxi1) goto 900

      iter = iter + 1
      h    = hood

      call dcopy( p, zero, 0, shape, 1)

      ulam = FLMAX

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        sum = zero
        do j = 1, p
          temp   = v(j,k)
          temp   = temp*temp
          v(j,k) = temp
          sum    = sum + temp/s(j)
        end do
        temp    = (sum/prob(k))/dble(p)
        vlam(k) = temp
        ulam    = min(temp,ulam)
        if (temp .gt. eps1) then
          do j = 1, p
            shape(j) = shape(j) + v(j,k)/temp
          end do
        end if
      end do

      vlamin = min(ulam, vlamin)

      if (ulam .le. eps1) goto 900

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc    = zero
        rcmin = zero
        goto 900
      end if 

c normalize the shape matrix
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = one/exp(sum/dble(p))

      call dscal( p, temp, shape, 1)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      goto 100

500   continue

      G  = -G

      G1 = G + 1

c probability = one/dble(G)

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        do j = 1, p
          temp     = v(j,k)
          temp     = temp*temp
          shape(j) = shape(j) + temp
          v(j,k)   = temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc     = zero
        rcmin  = zero
        ulam   = zero
        vlamin = zero
        goto 900
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))

      ulam = temp / dble(sumz)

      vlamin = ulam

      if (ulam .le. eps1) goto 900

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      call dcopy( G, ulam, 0, vlam, 1)

600   continue

c prob now contains n*prob

      if (maxi2 .le. 0) goto 630

      inner = 0

610   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        if (ulam .le. eps1) goto 900

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          rc    = zero
          rcmin = zero
          goto 900
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = one/exp(sum/dble(p))

        call dscal( p, temp, shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .le. eps2) goto 900

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol2) goto 620
        if (inner .lt. maxi2) goto 610

620   continue
     
      inmax = max(inner,inmax)
      errim = max(errin,errim)

630   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
        temp   = vlam(k)
        vlamin = min(temp,vlamin)
        if (vlamin .le. eps1) goto 900
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      hood  = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol1 ) goto 900
      if (iter .ge. maxi1) goto 900

      iter = iter + 1
      h    = hood

      call dcopy( p, zero, 0, shape, 1)

      ulam = FLMAX

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        sum = zero
        do j = 1, p
          temp   = v(j,k)
          temp   = temp*temp
          v(j,k) = temp
          sum    = sum + temp/s(j)
        end do
        temp    = (sum/prob(k))/dble(p)
        vlam(k) = temp
        ulam    = min(temp,ulam)
        if (temp .gt. eps1) then
          do j = 1, p
            shape(j) = shape(j) + v(j,k)/vlam(k)
          end do
        end if
      end do

      vlamin = min(ulam, vlamin)

      if (ulam .le. eps1) goto 900

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc    = zero
        rcmin = zero
        goto 900
      end if 

c normalize the shape matrix
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = one/exp(sum/dble(p))

      call dscal( p, temp, shape, 1)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      goto 600

900   continue

      tol(1)  = err
      if (inmax .eq. 0) then
        tol(2) = FLMAX
      else
        tol(2) = errim
      end if
      eps(1)  = vlamin
      eps(2)  = rcmin
      maxi(1) = iter
      maxi(2) = inmax

      return
      end
      subroutine menvi ( x, z, n, p, G, eps, tol, maxi, y, mu, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for spherical Gaussian mixtures with varying volumes plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G+1), y(n,G), mu(p)
      double precision   x(n,*), z(n,*), y(n,*), mu(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on sigsq. 
c                   On output, minimum value of sigsq at the final iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch) (n,G) (1,1 only if eps .lt. 0)
c  mu      double  (scratch) (p)
c  Vinv    double  (input) estimated reciprocal hypervolume of data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)
 
      h      = FLMAX/two

      sigmin = FLMAX

      iter   = 1

      if (eps .le. zero) goto 200

      if (G .lt. 0) goto 150

      G1 = G + 1

110   continue

      sigmk = FLMAX

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sumz), mu, 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          y(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        sigmk  = min(sigsqk,sigmk)
        if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            y(i,k) = probk*exp(-(const+y(i,k)/sigsqk)/two)           
          end do
        end if
      end do

      sigmin = min(sigmk,sigmin)

      if (sigmin .le. eps) goto 900

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G  = -G
      G1 = G + 1

      call dcopy( n, Vinv, 0, z(1,G1), 1)

c probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

160   continue

      sigmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        zetak = sum
c       call dblepr( "mean", -1, mu, p)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          y(i,k) = sum
        end do
        sigsqk = (sigsqk/zetak)/dble(p)
        sigmk  = min(sigsqk,sigmk)
        if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            y(i,k) = exp(-(const+y(i,k)/sigsqk)/two)           
          end do
        end if
      end do

      sigmin = min(sigmk,sigmin)

      if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      eps = -eps

      if (G .lt. 0) goto 250

      G1 = G + 1

210   continue

      sigmk = FLMAX

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sumz), mu, 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        sigmk  = min(sigsqk,sigmk)
        if (sigsqk .lt. eps) sigsqk = eps
c       if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c      end if
      end do

      sigmin = min(sigmk,sigmin)

c     if (sigmin .le. eps) goto 900

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G  = -G
      G1 = G + 1

      call dcopy( n, Vinv, 0, z(1,G1), 1)

c probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

260   continue

      sigmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        zetak = sum
c       call dblepr( "mean", -1, mu, p)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        sigmk  = min(sigsqk,sigmk)
        if (sigsqk .lt. eps) sigsqk = eps
c       if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c       end if
      end do

      sigmin = min(sigmk,sigmin)

c     if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = sigmk
      maxi = iter

      return
      end
      subroutine menvvv( x, z, n, p, G, eps, tol, maxi, 
     *                   y, mu, U, w, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for unconstrained Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G+1),y(n,G),mu(p),U(p,p),w(p)
      double precision   x(n,*), z(n,*), y(n,*), mu(*), U(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower bound on the reciprocal
c                   condition estimate of the covariances. On output, minimum 
c                    reciprocal condition estimate at the last iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch) (n,p) (1,1 if eps .lt. zero)
c  mu      double  (scratch) (p)
c  U       double  (scratch) (p,p)
c  w       double  (scratch) (p) 
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------
      if (maxi .le. 0) return

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX

      iter   = 1

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      G1 = G + 1

110   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .gt. eps) then
          detlog = log(abs(U(1,1)))
          do j = 2, p
            detlog = detlog + log(abs(U(j,j)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu, 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            y(i,k) = prob*exp(-(const+temp))
          end do
        end if
      end do

      rcmin = min(rcmk,rcmin)

      if (rcmin .le. eps) goto 900

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      hood   = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

160   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .gt. eps) then
          detlog = log(abs(U(1,1)))
          do j = 2, p
            detlog = detlog + log(abs(U(j,j)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu, 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
c           y(i,k) = prob*exp(-(const+temp))
            y(i,k) = exp(-(const+temp))
          end do
        end if
      end do

      rcmin = min(rcmk,rcmin)

      if (rcmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      G1 = G + 1

210   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, U, p1)
        end if
        detlog = log(abs(U(1,1)))
        do j = 2, p
          detlog = detlog + log(abs(U(j,j)))
        end do
        const = piterm+detlog
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob*exp(-(const+temp))
        end do
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn = (termn/dble(n))*Vinv

      hood   = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
c     h1 = -dble(n)*log(dble(G1))

260   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, U, p1)
        end if
        detlog = log(abs(U(1,1)))
        do j = 2, p
          detlog = detlog + log(abs(U(j,j)))
        end do
        const = piterm+detlog
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = exp(-(const+temp))
        end do
      end do

      hood = zero
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G1, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = rcmin*rcmin
      maxi = iter

      return
      end
      subroutine mevev ( x, z, n, p, G, eps, tol, maxi, 
     *                   mu, U, prob, vlam, shape, s, v, w, lwork)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      integer            maxi(2)
      double precision   eps(2), tol(2)

c     double precision   x(n,p), z(n,G), mu(p,G), U(p,p,G), prob(G),
c    *                   vlam(G),shape(p), s(p), v(p,G), w(max(lwork,G))
      double precision   x(n,*), z(n,*), mu(p,*), U(p,p,*), prob(*),
     *                   vlam(*), shape(*), s(*), v(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape.
c                   On output, minimum values of these at the last iteration.
c  tol     double  (input/output) (2) On input, tolerance on convergence of
c                   the loglikelihood and on the inner iterations for lambda
c                   and shape. On output, maximum relative error for these.
c  maxi    integer (input/output) (2) On input, upper limit on outer/inner
c                   iterations. On output, number of outer/inner iterations.
c  mu      double  (scratch) (p,G)
c  U       double  (scratch) (p,p,G)
c  prob    double  (scratch) (G) (needed with equal proportions too)
c  vlam    double  (scratch) (G)
c  shape   double  (scratch) (p)
c  s       double  (scratch) (p) 
c  v       double  (scratch) (p,G)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4,G) workspace for LAPACK
c                   SVD and other things

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

      if (maxi1 .le. 0) return

      eps1   = max(eps(1),zero)
      eps2   = max(eps(2),zero)

      tol1   = max(tol(1),zero)
      tol2   = max(tol(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX
      rc     = FLMAX

      errin  = FLMAX
      errim  = zero

      inmax  = 0

      iter   = 1

      if (G .lt. 0) goto 500

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        do j = 1, p
          temp     = v(j,k)
          temp     = temp*temp
          shape(j) = shape(j) + temp
          v(j,k)   = temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc     = zero
        rcmin  = zero
        ulam   = zero
        vlamin = zero
        goto 900
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))

      ulam = temp / dble(n)

      vlamin = ulam
  
      if (ulam .le. eps1) goto 900

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      call dcopy( G, ulam, 0, vlam, 1)

100   continue

c prob now contains n*prob

      if (maxi2 .le. 0) goto 130

      inner = 0

110   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        if (ulam .le. eps1) goto 900

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          rc    = zero
          rcmin = zero
          goto 900
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = one/exp(sum/dble(p))

        call dscal( p, temp, shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .le. eps2) goto 900

        errin = zero

        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol2) goto 120
        if (inner .lt. maxi2) goto 110

120   continue

      inmax = max(inner,inmax)
      errim = max(errin,errin)
        
130   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      call dscal( G, one/dble(n), prob, 1)

      do k = 1, G
        temp  = vlam(k)
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = prob(k)*exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do
 
      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol1 ) goto 900
      if (iter .ge. maxi1) goto 900

      iter = iter + 1
      h    = hood

      ulam = FLMAX

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        sum = zero
        do j = 1, p
          temp   = v(j,k)
          temp   = temp*temp
          v(j,k) = temp
          sum    = sum + temp/s(j)
        end do
        temp    = (sum/prob(k))/dble(p)
        vlam(k) = temp
        ulam    = min(temp,ulam)
        if (temp .gt. eps1) then
          do j = 1, p
            shape(j) = shape(j) + v(j,k)/temp
          end do
        end if
      end do

      vlamin = min(ulam, vlamin)

      if (ulam .le. eps1) goto 900

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc    = zero
        rcmin = zero
        goto 900
      end if 

c normalize the shape matrix
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = one/exp(sum/dble(p))

      call dscal( p, temp, shape, 1)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      goto 100

500   continue

      G    = -G

c probability = one/dble(G)

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        do j = 1, p
          temp     = v(j,k)
          temp     = temp*temp
          shape(j) = shape(j) + temp
          v(j,k)   = temp
        end do
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc     = zero
        rcmin  = zero
        ulam   = zero
        vlamin = zero
        goto 900
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))

      ulam = temp / dble(n)

      vlamin = ulam
  
      if (ulam .le. eps1) goto 900

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      call dcopy( G, ulam, 0, vlam, 1)

600   continue

c prob now contains n*prob

      if (maxi2 .le. 0) goto 630

      inner = 0

610   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        if (ulam .le. eps1) goto 900

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          rc    = zero
          rcmin = zero
          goto 900
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = one/exp(sum/dble(p))

        call dscal( p, temp, shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .le. eps2) goto 900

        errin = zero

        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol2) goto 620
        if (inner .lt. maxi2) goto 610

620   continue

      inmax = max(inner,inmax)
      errim = max(errin,errin)
        
630   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
        temp  = vlam(k)
        vlamk = min(temp,vlamk)
        if (temp .le. eps1) goto 900
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, U(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      vlamin = min(vlamk,vlamin)

      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol1 ) goto 900
      if (iter .ge. maxi1) goto 900

      iter = iter + 1
      h    = hood

      ulam = FLMAX

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j,k), w(j), cs, sn)
            call drot( p-j, U(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p,k), w(p), cs, sn)
        end do
        call dgesvd( 'N', 'O', p, p, U(1,1,k), p, v(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        sum = zero
        do j = 1, p
          temp   = v(j,k)
          temp   = temp*temp
          v(j,k) = temp
          sum    = sum + temp/s(j)
        end do
        temp    = (sum/prob(k))/dble(p)
        vlam(k) = temp
        ulam    = min(temp,ulam)
        if (temp .gt. eps1) then
          do j = 1, p
            shape(j) = shape(j) + v(j,k)/temp
          end do
        end if
      end do

      vlamin = min(ulam,vlamin)

      if (ulam .le. eps1) goto 900

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        rc    = zero
        rcmin = zero
        goto 900
      end if

c normalize the shape matrix
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = one/exp(sum/dble(p))

      call dscal( p, temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) goto 900

      goto 600

900   continue

      tol(1)  = err
      tol(2)  = errin
      eps(1)  = ulam
      eps(2)  = rc
      maxi(1) = iter
      maxi(2) = inmax

      return
      end
      subroutine mevi ( x, z, n, p, G, eps, tol, maxi, y, mu)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for spherical Gaussian mixtures with varying volumes

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G), y(n,G), mu(p)
      double precision   x(n,*), z(n,*), y(n,*), mu(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on sigsq. 
c                   On output, minimum value of sigsq at the last iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for 
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch) (n,G) (1,1 if eps .lt. zero)
c  mu      double  (scratch) (p)

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)
 
      h      = FLMAX/two

      sigmin = FLMAX

      iter   = 1

      if (eps .le. zero) goto 200

      if (G .lt. 0) goto 150

110   continue

      sigmk  = FLMAX

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sumz), mu, 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          y(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        sigmk  = min(sigsqk,sigmk)
        if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            y(i,k) = probk*exp(-(const+y(i,k)/sigsqk)/two)           
          end do
        end if
      end do

      sigmin = min(sigmk,sigmin)

      if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + y(i,k)
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G     = -G

c  probability = one/dble(G)

160   continue

      sigmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        zetak = sum
c       call dblepr( "mean", -1, mu, p)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          y(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        sigmk    = min(sigsqk,sigmk)
        if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            y(i,k) = exp(-(const+y(i,k)/sigsqk)/two)           
          end do
        end if
      end do

      sigmin = min(sigmk,sigmin)

      if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + y(i,k)
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      if (G .lt. 0) goto 250

      eps = -eps

210   continue

      sigmk  = FLMAX

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sumz), mu, 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        sigmk  = min(sigsqk,sigmk)
        sigsqk = max(sigsqk,eps)
c       if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c       end if
      end do

      sigmin = min(sigmk,sigmin)

c     if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G     = -G

c  probability = one/dble(G)

260   continue

      sigmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        zetak = sum
c       call dblepr( "mean", -1, mu, p)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/zetak)/dble(p)
        sigmk  = min(sigsqk,sigmk)
        sigsqk = max(sigsqk,eps)
c       if (sigsqk .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c       end if
      end do

      sigmin = min(sigmk,sigmin)

c     if (sigmin .le. eps) goto 900

      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = sigmk
      maxi = iter

      return
      end

      subroutine mevvv ( x, z, n, p, G, eps, tol, maxi, y, mu, U, w)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for unconstrained Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G),y(n,G),mu(p),U(p,p),w(p)
      double precision   x(n,*), z(n,*), y(n,*), mu(*), U(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower bound on the reciprocal
c                   condition estimate of the covariances. On output,
c                   minimum reciprocal condition estimate at the last
c                   iteration.
c  tol     double  (input/output) On input, tolerance on convergence of
c                   the loglikelihood. On output, maximum relative error for 
c                   the loglikelihood.
c  maxi    integer (input/output) On input, upper limit on iterations.
c                   On output, number of iterations
c  y       double  (scratch) (n,p) (1,1 if eps .lt. zero)
c  mu      double  (scratch) (p)
c  U       double  (scratch) (p,p)
c  w       double  (scratch) (p) 

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      tol    = max(tol,zero)

      FLMAX  = d1mach(2)

      h      = FLMAX/two

      rcmin  = FLMAX

      iter   = 1

      if (eps .le. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

110   continue

      rcmk   = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .gt. eps) then
          detlog = log(abs(U(1,1)))
          do j = 2, p
            detlog = detlog + log(abs(U(j,j)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu, 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            y(i,k) = prob*exp(-(const+temp))
          end do
        end if
      end do

      rcmin  = min(rcmk,rcmin)

      if (rcmin .le. eps) goto 900

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 110

150   continue

      G    = -G

c     prob = one/dble(G)

160   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .gt. eps) then
          detlog = log(abs(U(1,1)))
          do j = 2, p
            detlog = detlog + log(abs(U(j,j)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu, 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
c           y(i,k) = prob*exp(-(const+temp))
            y(i,k) = exp(-(const+temp))
          end do
        end if
      end do

      rcmin = min(rcmk,rcmin)

      if (rcmin .le. eps) goto 900

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = y(i,k)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 160

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

210   continue

      rcmk   = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        prob = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, U, p1)
        end if
        detlog = log(abs(U(1,1)))
        do j = 2, p
          detlog = detlog + log(abs(U(j,j)))
        end do
        const = piterm+detlog
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          z(i,k) = prob*exp(-(const+temp))
        end do
      end do

      rcmin  = min(rcmk,rcmin)

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 210

250   continue

      G    = -G

c     prob = one/dble(G)

260   continue

      rcmk = FLMAX

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu, 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu, 1)
        end do
        call dscal( p, (one/sum), mu, 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j), 1)
        end do
        call drnge( p, U, p1, umin, umax)
        rc   = umin/(one+umax)
        rcmk = min(rc,rcmk)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, U, p1)
        end if
        detlog = log(abs(U(1,1)))
        do j = 2, p
          detlog = detlog + log(abs(U(j,j)))
        end do
        const = piterm+detlog
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu, 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
c         z(i,k) = prob*exp(-(const+temp))
          z(i,k) = exp(-(const+temp))
        end do
      end do

      rcmin = min(rcmk,rcmin)

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
        call dscal( G, (one/sum), z(i,1), n)
      end do
      err  = abs(h-hood)/(one+abs(hood))

      if (err  .le. tol ) goto 900
      if (iter .ge. maxi) goto 900

      iter = iter + 1
      h    = hood

      goto 260

900   continue

      tol  = err
      eps  = rcmk*rcmk
      maxi = iter

      return
      end
      subroutine mseee ( x, z, n, p, G, eps, w, hood, mu, Sigma, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for constant-variance Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G),w(p),mu(p,G),Sigma(p,p),prob(G)
      double precision   x(n,*),z(n,*),w(*),mu(p,*),Sigma(p,*),prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on reciprocal 
c                   condition estimate for Sigma. If negative, a diagonal
c                   modification is used if necessary. On output, the 
c                   estimate.
c  w       double  (scratch) (p)
c  hood    double  (output) loglikelihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) Sigma.
c  prob    double  (output) (G) probabilities.

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      sclfac = one/sqrt(dble(n))

      FLMAX  = d1mach(2)

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, Sigma(1,j), 1)
      end do

c Sigma is now the Cholesky factor of Sigma

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .le. eps) then
        eps  = rc*rc
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc*rc

      return

150   continue

      G  = -G

c  probability = one/dble(G)

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, Sigma(1,j), 1)
      end do

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .le. eps) then
        eps = rc*rc
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = piterm + detlog

      hood  = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc*rc

      return

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, Sigma(1,j), 1)
      end do

c Sigma is now the Cholesky factor of Sigma

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, Sigma, p1)
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      if (rc .lt. eps) then
c undo the change
        call daxpy( p, (-one), alpha, 0, Sigma, p1)
      end if

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc*rc

      return

250   continue

      G  = -G

c  probability = one/dble(G)

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      do j = 1, p
        call dscal( j, sclfac, Sigma(1,j), 1)
      end do

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)
      rc    = rc*rc

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, Sigma, p1)
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = piterm + detlog

      hood  = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      if (rc .lt. eps) then
c undo the change
        call daxpy( p, (-one), alpha, 0, Sigma, p1)
      end if

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc

      return
      end
      subroutine mseev ( x, z, n, p, G, eps, shape, s, w, lwork,
     *                   hood, mu, Sigma, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G), shape(p), s(p), w(lwork),
c    *                   mu(p,G), Sigma(p,p,G), prob(G),
      double precision   x(n,*), z(n,*), shape(*), s(*), w(*),
     *                   mu(p,*), Sigma(p,p,*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape. On output,
c                   those values.
c  shape   double  (scratch) (p) 
c  s       double  (scratch) (p) 
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4) workspace for LAPACK SVD.
c  hood    doublr  (output) logliklihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) (p,p,G) sigmas.
c  prob    double  (output) (G) probabilities (unequal proportions only)

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps1   = max(eps(1),zero)
      eps2   = sqrt(max(eps(2),zero))

      dnp    = dble(n*p)

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)

      hood   = zero

      if (G .lt. 0) goto 500

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      if (abs(hood) .eq. FLMIN) return

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        hood   = FLMAX
        return
      else 
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/dble(n)
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) then
        eps(1) = vlam
        eps(2) = FLMAX
        hood   = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) then
        eps(1) = vlam
        eps(2) = rc*rc
        hood   = -FLMAX
        return
      end if

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam)
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = prob(k) * exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      eps(1) = vlam
      eps(2) = rc*rc

      return

500   continue

      G    = -G

c probability = one/dble(G)

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      if (abs(hood) .eq. FLMIN) return

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        hood   = FLMAX
        return
      else
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/dble(n)
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) then
        eps(1) = vlam
        eps(2) = FLMAX
        hood   = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) then
        eps(1) = vlam
        eps(2) = rc*rc
        hood   = -FLMAX
        return
      end if

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam)
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      eps(1) = vlam
      eps(2) = rc*rc

      return
      end
      subroutine msei  ( x, z, n, p, G, hood, mu, sigsq, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for spherical, constant-volume Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G), mu(p,G), prob(G)
      double precision   x(n,*), z(n,*), mu(p,*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  hood    double  (input/output) On input, lower limit on sigsq. If negative,
c                   sigsq is replaced by -hood should it fall below it. 
c                   On output, loglikelihood.
c  mu      double  (output) (p,G) means.
c  sigsq   double  (output) sigma-squared.
c  prob    double  (output) (G) probabilities.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps   = hood

      dnp   = dble(n*p)

      FLMAX = d1mach(2)

      if (eps .lt. zero) goto 200

      if (G .lt. 0) goto 150

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / dnp

      if (sigsq .le. eps) then
        hood = FLMAX
        return
      end if

      const  = dble(p)*(pi2log+log(sigsq))

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      return

150   continue

      G  = -G

c  probability = one/dble(G)

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / dnp

      if (sigsq .le. eps) then
        hood = FLMAX
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      return

200   continue
  
      eps = -eps

      if (G .lt. 0) goto 250

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / dnp

      if (sigsq .ge. eps) then
        const  = dble(p)*(pi2log+log(sigsq))

        hood   = zero
        do i = 1, n
          sum = zero
          do k = 1, G
            temp   = (const+(z(i,k)/sigsq))/two
            temp   = prob(k)*exp(-temp)
            sum    = sum + temp
          end do
          hood = hood + log(sum)
        end do
      else
        const  = dble(p)*(pi2log+log(eps))

        hood   = zero
        do i = 1, n
          sum = zero
          do k = 1, G
            temp   = (const+(z(i,k)/eps))/two
            temp   = prob(k)*exp(-temp)
            sum    = sum + temp
          end do
          hood = hood + log(sum)
        end do
      end if

      return

250   continue

      G  = -G

c  probability = one/dble(G)

      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / dnp

      if (sigsq .ge. eps) then
        const = dble(p)*(pi2log+log(sigsq))

        hood = -dble(n)*log(dble(G))
        do i = 1, n
          sum = zero
          do k = 1, G
            temp   = (const+(z(i,k)/sigsq))/two
            temp   = exp(-temp)
            sum    = sum + temp
          end do
          hood = hood + log(sum)
        end do
      else
        const = dble(p)*(pi2log+log(eps))

        hood = -dble(n)*log(dble(G))
        do i = 1, n
          sum = zero
          do k = 1, G
            temp   = (const+(z(i,k)/eps))/two
            temp   = exp(-temp)
            sum    = sum + temp
          end do
          hood = hood + log(sum)
        end do
      end if

      return
      end
      subroutine msneee( x, z, n, p, G, eps, w, 
     *                   hood, mu, Sigma, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for constant-variance Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G+1),w(p),mu(p,G),Sigma(p,p),prob(G+1)
      double precision   x(n,*),z(n,*),w(*),mu(p,*),Sigma(p,*),prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on reciprocal 
c                   condition estimate for Sigma. If negative, a diagonal
c                   modification is used if necessary. On output, the 
c                   estimate.
c  w       double  (scratch) (p)
c  hood    double  (output) loglikelihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) Sigma.
c  prob    double  (output) (G+1) probabilities.
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      FLMAX  = d1mach(2)

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      G1 = G + 1

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, Sigma(1,j), 1)
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .le. eps) then
        eps  = rc*rc
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do
  
      eps = rc*rc

      return

150   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, Sigma(1,j), 1)
      end do

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .le. eps) then
        eps  = rc*rc
        hood = FLMAX
        return
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc*rc

      return

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      G1 = G + 1

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, Sigma(1,j), 1)
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, Sigma, p1)
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = prob(k) * exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do
  
      eps = rc*rc

      return

250   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      do j = 1, p
        call dcopy( j, zero, 0, Sigma(1,j), 1)
      end do

      sumz = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p), w(p), cs, sn)
        end do
      end do

      temp = one/sqrt(sumz)
      do j = 1, p
        call dscal( j, temp, Sigma(1,j), 1)
      end do

c condition number

      call drnge( p, Sigma, p1, umin, umax)

      rc    = umin/(one+umax)

      if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
        alpha = (eps*(one+umax) - umin)/(one-eps)
        call daxpy( p, one, alpha, 0, Sigma, p1)
      end if

      detlog = log(abs(Sigma(1,1)))
      do j = 2, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const  = piterm + detlog

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp = ddot( p, w, 1, w, 1)/two
          temp   = exp(-(const+temp))
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,Sigma(1,i),1,Sigma(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      eps = rc*rc

      return
      end
      subroutine msneev( x, z, n, p, G, eps, shape, s, w, lwork,
     *                   hood, mu, Sigma, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixture with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G+1), shape(p), s(p), w(lwork),
c    *                   mu(p,G), Sigma(p,p,G), prob(G+1),
      double precision   x(n,*), z(n,*), shape(*), s(*), w(*),
     *                   mu(p,*), Sigma(p,p,*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and reciprocal condition estimate of shape. 
c                   On output, those estimates.
c  shape   double  (output) (p) shape matrix.
c  s       double  (scratch) (p) 
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4) workspace for LAPACK SVD.
c  hood    doublr  (output) logliklihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) (p,p,G) sigmas.
c  prob    double  (output) (G+1) probabilities (unequal proportions only)
c  Vinv    double  (output) approximate hypervolume of the data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps1   = max(eps(1),zero)
      eps2   = sqrt(max(eps(2),zero))

      dnp    = dble(n*p)

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)

      hood   = zero

      if (G .lt. 0) goto 500

      G1 = G + 1

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        sumz    = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else 
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      if (abs(hood) .eq. FLMIN) return

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        hood   = FLMAX
        return
      else
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/sumz
      end if

      if (vlam .le. eps1) then
        eps(1) = vlam
        eps(2) = FLMAX
        hood   = FLMAX
        return
      end if

      vlamin = min(vlam,vlamin)

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) then
        eps(1) = vlam
        eps(2) = rc*rc
        hood   = -FLMAX
        return
      end if

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam)
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = zero
      do i = 1, n
        sum = termn
        do k = 1, G
          temp   = prob(k) * exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      eps(1) = vlam
      eps(2) = rc

      return

500   continue

      G  = -G
      G1 = G + 1

c probability = one/dble(G1)

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else 
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        hood   = FLMAX
        return
      else
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))
        vlam = temp/sumz
      end if

      vlamin = min(vlam,vlamin)

      if (vlam .le. eps1) then
        eps(1) = vlam
        eps(2) = FLMAX
        hood   = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j)/temp)
      end do

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/(one+smax)
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) then
        eps(1) = vlam
        eps(2) = rc*rc
        hood   = -FLMAX
        return
      end if

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w, 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / shape(j)
          end do
          z(i,k) = ddot( p, s, 1, s, 1)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam)
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      const = dble(p)*(pi2log+log(vlam))

      hood  = -dble(n)*log(dble(G1))
      do i = 1, n
        sum = Vinv
        do k = 1, G
          temp   = exp(-(const+z(i,k)/vlam)/two)
          sum    = sum + temp
        end do
        hood = hood + log(sum)
      end do

      eps(1) = vlam
      eps(2) = rc*rc

      return
      end
      subroutine msnei ( x, z, n, p, G, hood, mu, sigsq, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for spherical, constant-volume Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G+1), mu(p,G), prob(G+1)
      double precision   x(n,*), z(n,*), mu(p,*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  hood    double  (input/output) On input, lower limit on sigsq. If negative,
c                   sigsq is replaced by -hood if it falls below it. On output,
c                   loglikelihood.
c  mu      double  (output) (p,G) means.
c  sigsq   double  (output) sigma-squared.
c  prob    double  (output) (G+1) probabilities.
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps    = hood

      FLMAX  = d1mach(2)

      if (eps .lt. zero) goto 200

      if (G .lt. 0) goto 150

      G1 = G + 1

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum/dble(n)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / (dble(p)*sumz)

      if (sigsq .le. eps) then
        hood = FLMAX
        return
      end if

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      return

150   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / (dble(p)*sumz)

      if (sigsq .lt. eps) sigsq = eps

      const  = dble(p)*(pi2log+log(sigsq))

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

200   continue

      eps = -eps

      if (G .lt. 0) goto 250

      G1 = G + 1

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum/dble(n)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq  = sigsq / (dble(p)*sumz)

      if (sigsq .le. eps) then
        hood = FLMAX
        return
      end if

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      const  = dble(p)*(pi2log+log(sigsq))

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = prob(k)*exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      return

250   continue

      G  = -G
      G1 = G + 1

c  probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      sigsq = zero

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        sumz = sumz + sum
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
          z(i,k) = sum
        end do
      end do

      sigsq = sigsq / (dble(p)*sumz)

      if (sigsq .lt. eps) sigsq = eps

      const  = dble(p)*(pi2log+log(sigsq))

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          temp   = (const+(z(i,k)/sigsq))/two
          temp   = exp(-temp)
          sum    = sum + temp
          z(i,k) = temp
        end do
        hood = hood + log(sum)
      end do

      return
      end
      subroutine msnvev( x, z, n, p, G, eps, tol, maxi,
     *                   vlam, shape, s, v, w, lwork,
     *                   hood, mu, Sigma, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G+1), mu(p,G), Sigma(p,p,G), prob(G+1),
c    *                   vlam(G),shape(p), s(p), v(p,G), w(max(lwork,G))
      double precision   x(n,*), z(n,*), mu(p,*), Sigma(p,p,*), prob(*),
     *                   vlam(*), shape(*), s(*), v(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape. On output, values.
c  tol     double  (input/output) On input tolerance for inner iteration 
c                   convergence. On output, error for inner iteration.
c  maxi    integer (input/output) On input, upper limit on number of inner 
c                   iterations. On output, number of inner iterations. 
c  vlam    double  (output) (G)
c  shape   double  (output) (p)
c  s       double  (scratch) (p) 
c  v       double  (scratch) (p,G)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4,G) workspace for LAPACK
c                   SVD and other things
c  mu      double  (output) (p,G)
c  Sigma   double  (output) (p,p,G)
c  prob    double  (output) (G+1) (needed with equal proportions too)
c  Vinv    double  (input) approximate hypervolume of the data region

      integer                 p1, G1

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps1   = max(eps(1),zero)
      eps2   = max(eps(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)

      rcmin  = FLMAX
      vlamin = FLMAX
      errin  = FLMAX

      hood   = zero

      inner  = 0

      if (G .lt. 0) goto 500

      G1 = G + 1

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, v(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = v(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              v(j,k)   = temp
            end do
          end if
        end if
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      if (abs(hood) .eq. FLMIN) then
        call dscal( G, one/dble(n), prob, 1)
        eps(1) = FLMAX
        eps(2) = FLMAX
        tol    = FLMAX
        maxi   = 0
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))
      ulam = temp / dble(n)

      if (ulam .le. eps1) then
        eps(1) = ulam
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rcmin = smin/smax

      if (rcmin .le. eps2) then
        eps(1) = ulam
        eps(2) = rcmin
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dcopy( G, ulam, 0, vlam, 1)

c prob now contains n*prob

      if (maxi .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)
        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          eps(1) =  zero
          eps(2) =  zero
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        call dscal( p, one/exp(sum/dble(p)), shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .lt. eps2) then
          eps(1) =  ulam
          eps(2) =  rc
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol ) goto 120
        if (inner .lt. maxi) goto 110

120   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      call dscal( G, one/dble(n), prob, 1)

      do k = 1, G
        temp   = vlam(k)
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w(p1), 1, 
     *                 zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = prob(k)*exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam(k))
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      hood  = zero
      do i = 1, n
        sum = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      eps(1) = ulam
      eps(2) = rcmin
 
      tol    = errin
  
      maxi   = inner

      return

500   continue

      G  = -G
      G1 = G + 1

c probability = one/dble(G)

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, v(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = v(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              v(j,k)   = temp
            end do
          end if
        end if
      end do

      if (abs(hood) .eq. FLMIN) then
        call dscal( G, one/dble(n), prob, 1)
        eps(1) = FLMAX
        eps(2) = FLMAX
        tol    = FLMAX
        maxi   = 0
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))
      ulam = temp / dble(n)

      if (ulam .le. eps1) then
        eps(1) = ulam
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rcmin = smin/(one+smax)

      if (rcmin .le. eps2) then
        eps(1) = ulam
        eps(2) = rcmin
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dcopy( G, ulam, 0, vlam, 1)

c prob now contains n*prob

      if (maxi .le. 0) goto 620

610   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)
        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          eps(1) =  zero
          eps(2) =  zero
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        call dscal( p, one/exp(sum/dble(p)), shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/(one+smax)
        rcmin = min(rc,rcmin)

        if (rc .lt. eps2) then
          eps(1) =  ulam
          eps(2) =  rc
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol ) goto 620
        if (inner .lt. maxi) goto 610

620   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
        temp   = vlam(k)
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w(p1), 1, 
     *                 zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam(k))
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      hood  = -dble(n)*log(dble(G1))
      do i = 1, n
        sum = Vinv
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      eps(1) = ulam
      eps(2) = rc
 
      tol    = errin
  
      maxi   = inner

      return
      end
      subroutine msnvi ( x, z, n, p, G, hood, mu, sigsq, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for spherical Gaussian mixtures with varying volumes and Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G+1), mu(p,G), sigsq(G), prob(G)
      double precision   x(n,*), z(n,*), mu(p,*), sigsq(*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  hood    double  (input/output) On input, lower limit on sigsq. If negative,
c                   sigsq is replaced by -hood if it falls below it. On output,
c                   loglikelihood.
c  mu      double  (output) (p,G) means.
c  sigsq   double  (output) (G) sigma-squared.
c  prob    double  (output) (G+1) probabilities.
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------
   
      eps    = hood

      FLMAX  = d1mach(2)

      sigmin = FLMAX

      if (eps .lt. zero) goto 200
 
      if (G .lt. 0) goto 150

      G1 = G + 1

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sumz), mu(1,k), 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        if (sigsqk .lt. sigmin) sigmin = sigsqk
        if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
        end if
        prob(k)  = probk
        sigsq(k) = sigsqk
      end do

      if (sigmin .le. eps) then
        hood = FLMAX
        return
      end if

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return

150   continue

      G  = -G
      G1 = G + 1

c probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        zetak = sum
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        if (sigsqk .lt. sigmin) sigmin = sigsqk
        if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)           
          end do
        end if
        sigsq(k) = sigsqk
      end do

      if (sigmin .le. eps) then
        hood = FLMAX
        return
      end if

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum  = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do
 
200   continue

      eps = -eps

      if (G .lt. 0) goto 250

      G1 = G + 1

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sumz), mu(1,k), 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/sumz)/dble(p)
        sigsq(k) = sigsqk
        probk    = sumz / dble(n)
        sigmin   = min(sigmin,sigsqk)
        if (sigsqk .lt. eps) sigsqk = eps
c       if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c       end if
        prob(k)  = probk
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      hood = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return

250   continue

      G  = -G
      G1 = G + 1

c probability = one/dble(G1)
      h1 = -dble(n)*log(dble(G1))

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        zetak = sum
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        sigsq(k) = sigsqk
        sigmin   = min(sigmin,sigsqk)
        if (sigsqk .lt. eps) sigsqk = eps
c       if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)           
          end do
c       end if
      end do

      hood = h1
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum  = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return
      end
      subroutine msnvvv( x, z, n, p, G, eps, w, 
     *                   hood, mu, Sigma, prob, Vinv)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for unconstrained Gaussian mixtures plus Poisson noise

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G+1),w(p),mu(p,G),Sigma(p,p,G),prob(G+1)
      double precision   x(n,*),z(n,*),w(*),mu(p,*),Sigma(p,p,*),prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  z       double  (input/output) (n,G+1) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on reciprocal
c                   condition estimate for covariances. If negative, a
c                   diagonal modification is used if necessary. On output,
c                   minimum estimate.
c  w       double  (scratch) (p)
c  hood    double  (output) loglikelihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) (p,p,G) Sigmas.
c  prob    double  (output) (G+1) probabilities.
c  Vinv    double  (double) estimated reciprocal hypervolume of data region

      integer                 p1, G1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      FLMAX  = d1mach(2)

      rcmin  = FLMAX

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      G1 = G + 1

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = prob(k)*exp(-(const+temp))
          end do
        end if
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      if (rcmin .le. eps) then
        eps  = rcmin*rcmin
        hood = FLMAX
        return
      end if

      hood   = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return

150   continue

      G  = -G
      G1 = G + 1

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = exp(-(const+temp))
          end do
        end if
      end do

      if (rcmin .le. eps) then
        eps  = rcmin*rcmin
        hood = FLMAX
        return
      end if

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin

      return

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      G1 = G + 1

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, Sigma(1,1,k), p1)
        end if
c       if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = prob(k)*exp(-(const+temp))
          end do
c       end if
      end do

      termn = zero
      do i = 1, n
        termn = termn + z(i,G1)
      end do
      termn    = termn/dble(n)
      prob(G1) = termn
      termn    = termn*Vinv

      hood   = zero
      do i = 1, n
        sum     = termn
        z(i,G1) = termn
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return

250   continue

      G  = -G
      G1 = G + 1

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, Sigma(1,1,k), p1)
        end if
c       if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = exp(-(const+temp))
          end do
c       end if
      end do

      hood = -dble(n)*log(dble(G1))
      do i = 1, n
        sum     = Vinv
        z(i,G1) = Vinv
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return
      end
      subroutine msvev ( x, z, n, p, G, eps, tol, maxi,
     *                   vlam, shape, s, v, w, lwork,
     *                   hood, mu, Sigma, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c EM for Gaussian mixtures with prescribed shape and constant volume

      implicit double precision (a-h,o-z)

      integer            n, p, G

      double precision   eps(2)

c     double precision   x(n,p), z(n,G), mu(p,G), Sigma(p,p,G), prob(G),
c    *                   vlam(G),shape(p), s(p), v(p,G), w(max(lwork,G))
      double precision   x(n,*), z(n,*), mu(p,*), Sigma(p,p,*), prob(*),
     *                   vlam(*), shape(*), s(*), v(p,*), w(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) (2) On input, lower bounds on lamba
c                   estimate and condition number of shape. On output, values.
c  tol     double  (input/output) On input tolerance for inner iteration 
c                   convergence. On output, error for inner iteration.
c  maxi    integer (input/output) On input, upper limit on number of inner 
c                   iterations. On output, number of inner iterations. 
c  vlam    double  (output) (G)
c  shape   double  (output) (p)
c  s       double  (scratch) (p) 
c  v       double  (scratch) (p,G)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4,G) workspace for LAPACK
c                   SVD and other things
c  mu      double  (output) (p,G)
c  Sigma   double  (output) (p,p,G)
c  prob    double  (output) (G) (needed with equal proportions too)

      integer                 p1

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps1   = max(eps(1),zero)
      eps2   = max(eps(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

      FLMIN  = d1mach(1)
      FLMAX  = d1mach(2)

      rcmin  = FLMAX

      vlamin = FLMAX
      errin  = FLMAX

      hood   = zero

      inner  = 0

      if (G .lt. 0) goto 500

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, v(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = v(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              v(j,k)   = temp
            end do
          end if
        end if
      end do

      if (abs(hood) .eq. FLMIN) then
        call dscal( G, one/dble(n), prob, 1)
        eps(1) = FLMAX
        eps(2) = FLMAX
        tol    = FLMAX
        maxi   = 0
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))
      ulam = temp / dble(n)

      if (ulam .le. eps1) then
        eps(1) = ulam
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rc    = smin/smax
      rcmin = min(rc,rcmin)

      if (rc .le. eps2) then
        eps(1) = ulam
        eps(2) = rc
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dcopy( G, ulam, 0, vlam, 1)

      vlamin = ulam

c prob now contains n*prob

      if (maxi .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        if (ulam .le. eps1) then
          eps(1) = ulam
          eps(2) = rc
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

        vlamin = min(ulam,vlamin)

        call drnge( p, shape, 1, smin, smax)
 
        if (smin .eq. zero) then
          eps(1) =  ulam
          eps(2) =  zero
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        call dscal( p, one/exp(sum/dble(p)), shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rcmin = smin/(one+smax)

        if (rcmin .lt. eps2) then
          eps(1) =  ulam
          eps(2) =  rc
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol ) goto 120
        if (inner .lt. maxi) goto 110

120   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      call dscal( G, one/dble(n), prob, 1)

      do k = 1, G
        temp   = vlam(k)
        vlamin = min(temp,vlamin)
      end do

      if (vlamin .le. eps1) then
        eps(1) =  vlamin
        eps(2) =  rc
        tol    =  errin
        maxi   =  inner
        hood   =  FLMAX
        return
      end if

      do k = 1, G
        temp  = vlam(k)
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w(p1), 1, 
     *                 zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = prob(k)*exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam(k))
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      hood  = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      eps(1) = ulam
      eps(2) = rcmin
 
      tol    = errin
  
      maxi   = inner

      return

500   continue

      G  = -G

c probability = one/dble(G)

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, Sigma(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        if (abs(hood) .ne. FLMIN) then
          call dgesvd( 'N', 'O', p, p, Sigma(1,1,k), p, v(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            hood = dble(sign(1,info))*FLMIN
          else
            do j = 1, p
              temp     = v(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              v(j,k)   = temp
            end do
          end if
        end if
      end do

      if (abs(hood) .eq. FLMIN) then
        call dscal( G, one/dble(n), prob, 1)
        eps(1) = FLMAX
        eps(2) = FLMAX
        tol    = FLMAX
        maxi   = 0
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps(1) = zero
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))
      ulam = temp / dble(n)

      if (ulam .le. eps1) then
        eps(1) = ulam
        eps(2) = zero
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call drnge( p, shape, 1, smin, smax)

      rcmin = smin/smax

      if (rcmin .le. eps2) then
        eps(1) = ulam
        eps(2) = rcmin
        tol    =  FLMAX
        maxi   = 0
        hood   = -FLMAX
        return
      end if

      call dcopy( G, ulam, 0, vlam, 1)

      vlamin = ulam

c prob now contains n*prob

      if (maxi .le. 0) goto 620

610   continue

        call dcopy( p, shape, 1, s, 1)
        call dcopy( G, vlam, 1, w, 1)

        call dcopy( p, zero, 0, shape, 1)

        ulam = FLMAX

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + v(j,k)/s(j)
          end do
          temp    = (sum/prob(k))/dble(p)
          vlam(k) = temp
          ulam    = min(temp,ulam)
          if (temp .gt. eps1) then
            do j = 1, p
              shape(j) = shape(j) + v(j,k)/temp
            end do
          end if
        end do

        inner  = inner + 1

        vlamin = min(ulam,vlamin)

        call drnge( p, shape, 1, smin, smax)
 
        if (smin .eq. zero) then
          eps(1) =  zero
          eps(2) =  zero
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        call dscal( p, one/exp(sum/dble(p)), shape, 1)

        call drnge( p, shape, 1, smin, smax)

        rc    = smin/smax
        rcmin = min(rc,rcmin)

        if (rc .lt. eps2) then
          eps(1) =  ulam
          eps(2) =  rc
          tol    =  errin
          maxi   =  inner
          hood   = -FLMAX
          return
        end if

        errin = zero
        do j = 1, p
          errin = max(abs(s(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(vlam(k)-w(k))/(one+vlam(k)), errin)
        end do
        
        if (errin .lt. tol ) goto 620
        if (inner .lt. maxi) goto 610

620   continue

      do j = 1, p
        s(j)     = shape(j)
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
        temp   = vlam(k)
        const = dble(p)*(pi2log + log(temp))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, Sigma(1,1,k), p, w(p1), 1, 
     *                 zero, w, 1)
          do j = 1, p
            w(j) = w(j) / shape(j)
          end do
          z(i,k) = exp(-(const+ddot(p,w,1,w,1)/temp)/2)
        end do
      end do

      do k = 1, G
        temp = sqrt(vlam(k))
        do i = 1, p
          call dscal( p, temp*shape(i), Sigma(i,1,k), p)
        end do
      end do

      do k = 1, G
        call dsyrk('U','T',p,p,one,Sigma(1,1,k),p,zero,x,n)
        do i = 2, p
          do j = 1, (i-1)
            x(i,j) = x(j,i)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      eps(1) = ulam
      eps(2) = rc
 
      tol    = errin
  
      maxi   = inner

      return
      end
      subroutine msvi ( x, z, n, p, G, hood, mu, sigsq, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for spherical Gaussian mixtures with varying volumes

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p), z(n,G), mu(p,G), sigsq(G), prob(G)
      double precision   x(n,*), z(n,*), mu(p,*), sigsq(*), prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  hood    double  (input/output) On input, lower limit on sigsq. If negative,
c                   sigsq is replaced by -hood should it fall below it.
c                   On output, loglikelihood.
c  mu      double  (output) (p,G) means.
c  sigsq   double  (output) (G) sigma-squared.
c  prob    double  (output) (G) probabilities.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      eps    = hood
  
      FLMAX  = d1mach(2)

      sigmin = FLMAX

      if (eps .lt. zero) goto 200

      if (G .lt. 0) goto 150

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sumz), mu(1,k), 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        if (sigsqk .lt. sigmin) sigmin = sigsqk
        if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
        end if
        prob(k)  = probk
        sigsq(k) = sigsqk
      end do

      if (sigmin .le. eps) then
        hood = FLMAX
        return
      end if

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return

150   continue

      G  = -G

c  probability = one/dble(G)

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        zetak = sum
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        if (sigsqk .lt. sigmin) sigmin = sigsqk
        if (sigmin .gt. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)           
          end do
        end if
        sigsq(k) = sigsqk
      end do

      if (sigmin .le. eps) then
        hood = FLMAX
        return
      end if

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return

200   continue

      eps = -eps

      if (G .lt. 0) goto 250

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sumz), mu(1,k), 1)
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk = (sigsqk/sumz)/dble(p)
        probk  = sumz / dble(n)
        sigmin = min(sigsqk,sigmin)
        if (sigsqk .ge. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/sigsqk)/two)           
          end do
        else
          const = dble(p)*(pi2log+log(eps))
          do i = 1, n
            z(i,k) = probk*exp(-(const+z(i,k)/eps)/two)           
          end do
        end if
        sigsq(k) = sigsqk
        prob(k)  = probk
      end do

      hood   = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return

250   continue

      G  = -G

c  probability = one/dble(G)

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        zetak = sum
        sigsqk = zero
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsqk = sigsqk + z(i,k)*sum
          z(i,k) = sum
        end do
        sigsqk   = (sigsqk/zetak)/dble(p)
        sigmin = min(sigsqk,sigmin)
        if (sigsqk .ge. eps) then
          const = dble(p)*(pi2log+log(sigsqk))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/sigsqk)/two)
          end do
        else
          const = dble(p)*(pi2log+log(eps))
          do i = 1, n
            z(i,k) = exp(-(const+z(i,k)/eps)/two)
          end do
        end if
        sigsq(k) = sigsqk
      end do

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      return
      end
      subroutine msvvv ( x, z, n, p, G, eps, w, hood, mu, Sigma, prob)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for unconstrained Gaussian mixtures

      implicit double precision (a-h,o-z)

      integer            n, p, G

c     double precision   x(n,p),z(n,G),w(p),mu(p,G),Sigma(p,p,G),prob(G)
      double precision   x(n,*),z(n,*),w(*),mu(p,*),Sigma(p,p,*),prob(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  z       double  (input/output) (n,G) Initial/final values for the
c                   conditional probabilities. 
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  G       integer (input) number of Gaussian clusters in the mixture;
c                   negative indicates equal mixing proportions
c  eps     double  (input/output) On input, lower limit on reciprocal 
c                   condition estimate for covariances. If negative, a
c                   diagonal modification is used if necessary. On output, 
c                   minimum estimate.
c  w       double  (scratch) (p)
c  hood    double  (output) loglikelihood.
c  mu      double  (output) (p,G) means.
c  Sigma   double  (output) (p,p,G) Sigmas.
c  prob    double  (output) (G) probabilities.

      integer                 p1

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c------------------------------------------------------------------------------

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      FLMAX  = d1mach(2)

      rcmin  = FLMAX

      if (eps .lt. zero) goto 200

      eps = sqrt(eps)

      if (G .lt. 0) goto 150

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = prob(k)*exp(-(const+temp))
          end do
        end if
      end do

      if (rcmin .le. eps) then
        eps  = rcmin*rcmin
        hood = FLMAX
        return
      end if
 
      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return

150   continue

      G  = -G

c probability = one/dble(G)

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = exp(-(const+temp))
          end do
        end if
      end do

      if (rcmin .le. eps) then
        eps  = rcmin*rcmin
        hood = FLMAX
        return
      end if

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return

200   continue

      eps = sqrt(-eps)

      if (G .lt. 0) goto 250

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        prob(k) = sum / dble(n)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, Sigma(1,1,k), p1)
        end if
c       if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = prob(k)*exp(-(const+temp))
          end do
c       end if
        if (rc .lt. eps) then
c restore Sigma
          call daxpy( p, (-one), alpha, 0, Sigma(1,1,k), p1)
        end if
      end do
 
      hood = zero
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return

250   continue

      G  = -G

c probability = one/dble(G)

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, Sigma(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dscal( p, sqrt(z(i,k)), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( Sigma(j,j,k), w(j), cs, sn)
            call drot( p-j, Sigma(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( Sigma(p,p,k), w(p), cs, sn)
        end do
        do j = 1, p
          call dscal( j, one/sqrt(sum), Sigma(1,j,k), 1)
        end do
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc    = umin/(one+umax)
        rcmin = min(rc,rcmin)
        if (rc .lt. eps) then
c boost the diagonal if rc falls below threshold
          alpha = (eps*(one+umax) - umin)/(one-eps)
          call daxpy( p, one, alpha, 0, Sigma(1,1,k), p1)
        end if
c       if (rc .gt. eps) then
          detlog = log(abs(Sigma(1,1,k)))
          do j = 2, p
            detlog = detlog + log(abs(Sigma(j,j,k)))
          end do
          const = piterm+detlog
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
            temp   = ddot( p, w, 1, w, 1)/two
            z(i,k) = exp(-(const+temp))
          end do
c       end if
        if (rc .lt. eps) then
c restore Sigma
          call daxpy( p, (-one), alpha, 0, Sigma(1,1,k), p1)
        end if
      end do

      hood = -dble(n)*log(dble(G))
      do i = 1, n
        sum = zero
        do k = 1, G
          sum    = sum + z(i,k)
        end do
        hood = hood + log(sum)
      end do

      do k = 1, G
        do j = 1, p
          do i = 1, j
            x(i,j) = ddot(i,Sigma(1,i,k),1,Sigma(1,j,k),1)
            if (i .ne. j) x(j,i) = x(i,j)
          end do
        end do
        do j = 1, p
          call dcopy( p, x(1,j), 1, Sigma(1,j,k), 1)
        end do
      end do

      eps = rcmin*rcmin

      return
      end
      subroutine onexev( x, n, p, s, w, lwork, mu, Sigma, hood, rcond)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c M-step for Gaussian mixtures with prescribed shape, 
c varying volume and orientation

      implicit double precision (a-h,o-z)

      integer            n, p

c     double precision   x(n,p),s(p),w(lwork),mu(p),Sigma(p,p)
      double precision   x(n,*), s(*), w(*), mu(*), Sigma(p,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on exit.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  s       double  (scratch) (p)
c  w       double  (scratch) (lwork)
c  lwork   integer (input) .ge. max(4*p,5*p-4) workspace for LAPACK SVD.
c  mu      double  (output) (p) mean.
c  Sigma   double  (output) (p,p) Sigma.
c  prob    double  (output) (G) probabilities.
c  hood    double  (output) loglikelihood.
c  rcond   double  (output) reciprocal condition estimate.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      external                d1mach
      double precision        d1mach

c---------------------------------------------------------------------------

      FLMIN = d1mach(1)
      FLMAX = d1mach(2)

      scl = (one/dble(n))
      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
        call dcopy( p, zero, 0, Sigma(1,j), 1)
      end do
   
      do i = 1, n
        call daxpy( p, (-one), mu, 1, x(i,1), n)
        j = 1
        do j1 = 2, p
          call drotg( Sigma(j,j), x(i,j), cs, sn)
          call drot( p-j, Sigma(j,j1), p, x(i,j1), n, cs, sn)
          j = j1
        end do
        call drotg( Sigma(p,p), x(i,p), cs, sn)
      end do 

      call dgesvd( 'N', 'O', p, p, Sigma, p, s, 
     *              dummy, 1, dummy, 1, w, lwork, info)
      if (info .ne. 0) then
        hood = dble(sign(1,info))*FLMIN
        return
      end if
      do j = 1, p
        temp = s(j)
        s(j) = temp*temp
      end do

      call drnge( p, s, 1, smin, smax)

      rcond = smin / (one + smax)

      if (smin .eq. zero) then
        hood = FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(s(j))
      end do
      vlmlog = sum/dble(p)
      vlam   = exp(vlmlog)

      do j = 1, p
        s(j) = sqrt(s(j)/vlam)
      end do

      vlam   = vlam/dble(n)

      temp = sqrt(vlam)
      do i = 1, p
        call dscal( p, temp*s(i), Sigma(i,1), p)
      end do

      call dsyrk('U','T',p,p,one,Sigma,p,zero,x,n)
      do i = 2, p
        do j = 1, (i-1)
          x(i,j) = x(j,i)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, Sigma(1,j), 1)
      end do

      hood = -dble(n*p) * ((pi2log+one) + log(vlam)) / two

      return
      end
      subroutine onexi( x, n, p, mu, sigsq, hood)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c loglikelihood for a single Gaussian spherical cluster

      implicit double precision (a-h,o-z)

c     integer            n, p
      integer            n, p

c     double precision   x(n,p), mu(p)
      double precision   x(n,*), mu(*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on return.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  mu      double  (ouput) (p) MLE for cluster mean.
c  sigsq   double  (output) sigma-squared for the cluster.
c  hood    double  (output) loglikelihood

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        d1mach
      external                d1mach

c------------------------------------------------------------------------------
      
      FLMAX = d1mach(2)

      dnp = dble(n*p)

      scl = (one/dble(n))

      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
      end do

      sigsq = zero

      do j = 1, p
        call daxpy( n, (-one), mu(j), 0, x(1,j), 1)
        sigsq = sigsq + ddot( n, x(1,j), 1, x(1,j), 1)
      end do
 
      sigsq = sigsq/dnp

      if (sigsq .eq. zero) then
        hood = FLMAX
      else 
        hood  = -dnp*(pi2log + (one + log(sigsq)))/two
      end if

      return
      end
      subroutine onexxx( x, n, p, mu, R, hood, rcond)

c copyright 1997 Department of Statistics, University of Washington
c funded by ONR contracts N00014-96-1-0192 and N-00014-96-1-0330
c Commercial use or distribution prohibited except by agreement with
c the University of Washington.

c loglikelihood for a single Gaussian unconstrained cluster

      implicit double precision (a-h,o-z)

c     integer            n, p
      integer            n, p

c     double precision   x(n,p), mu(p), R(p,p)
      double precision   x(n,*), mu(*), R(p,*)
c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) Matrix of observations. Destroyed on return.
c  n       integer (input) number of observations
c  p       integer (input) dimension of the data
c  mu      double  (ouput) (p) mean.
c  R       double  (output) sigma.
c  hood    double  (output) logliklihood.
c  rcond   double  (output) reciprocal condition estimate for chol of Sigma.

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        d1mach
      external                d1mach

c------------------------------------------------------------------------------

      FLMAX = d1mach(2)

      dnp = dble(n*p)

      scl = (one/dble(n))
      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
        call dcopy( p, zero, 0, R(1,j), 1)
      end do

      do i = 1, n
        call daxpy( p, (-one), mu, 1, x(i,1), n)
        j = 1
        do j1 = 2, p
          call drotg( R(j,j), x(i,j), cs, sn)
          call drot( p-j, R(j,j1), p, x(i,j1), n, cs, sn)
          j = j1
        end do
        call drotg( R(p,p), x(i,p), cs, sn)
      end do

      scl = sqrt(scl)
      do j = 1, p
        call dscal( j, scl, R(1,j), 1)
      end do

      call drnge( p, R, p+1, rmin, rmax)

      rcond = rmin / (one + rmax)
      rcond = rcond*rcond

      if (rmin .eq. zero) then
        hood = FLMAX
      else
        hood = zero
          
        do j = 1, p
          hood = hood + log(abs(R(j,j)))
        end do

        hood = -dble(n)*(hood + dble(p)*(pi2log + one)/two)
      end if

      do j = 1, p
        do i = 1, j
          x(i,j) = ddot(i,R(1,i),1,R(1,j),1)
          if (i .ne. j) x(j,i) = x(i,j)
        end do
      end do
      do j = 1, p
        call dcopy( p, x(1,j), 1, R(1,j), 1)
      end do

      return
      end

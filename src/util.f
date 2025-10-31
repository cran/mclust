*  =====================================================================
      subroutine dmvnorm ( x, mu, Sigma, n, p, w, hood, logdens)
*
*  Compute log-density of multivariate Gaussian
*
*  =====================================================================

      implicit NONE

      integer            n, p
      double precision   hood
      double precision   x(n,p), w(*), logdens(n)
      double precision   mu(p), Sigma(p,p)

      integer            info, i, j

      double precision   detlog, umin, umax, const, temp

      double precision   zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision   FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision   RTMAX
      parameter          (RTMAX = 1.340780792994260d154)

      double precision   RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision   ddot
      external           ddot

* ---------------------------------------------------------------------------


*     Cholesky factorization
      call dpotrf('U', p, Sigma, p, info)

      if (info .ne. 0) then
        w(1) = dble(info)
        hood = FLMAX
        return
      end if

      call absrng( p, Sigma, (p+1), umin, umax)

      if (umax .le. one .and. umax .ge. umin*RTMAX) then
        w(1) = zero
        hood = FLMAX
        return
      end if

      if (umax .ge. one .and. umin .le. umax*RTMIN) then
        w(1) = zero
        hood = FLMAX
        return
      end if

      detlog = zero
      do j = 1, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = dble(p)*pi2log/two + detlog

      do i = 1, n
        call dcopy(p, x(i,1), n, w, 1)
        call daxpy(p, (-one), mu(1), 1, w, 1)
        call dtrsv('U', 'T', 'N', p, Sigma, p, w, 1)
        temp = ddot(p, w, 1, w, 1)/two
        logdens(i) =  -(const+temp)
      end do

      w(1) = zero

      return
      end

*  =====================================================================
      subroutine logsumexp(x, n, k, v, lse)
*  
*  Efficiently computes log-sum-exp(x_i+v) for i = 1,...,n
*  x = matrix (n x k)
*  v = vector (k)
*  lse = output vector(n)
*
*  =====================================================================
      implicit none
      integer          :: n, k, i
      double precision :: x(n,k), v(k), lse(n), xv(k), m

      do i = 1, n
        xv = x(i,:) + v
        m = maxval(xv)
        lse(i) = m + log(sum(exp(xv - m)))
      end do

      return
      end

*  =====================================================================
      subroutine softmax(x, n, k, v, lse, z)
*  
*  Efficiently computes softmax function based on 
*    exp(x_i+v - log-sum-exp(x_i+v))   for i = 1,...,n
*  x = matrix (n x k)
*  v = vector (k)
*  z = output matrix (n x k) with rowsum(z_j) = 1 for j = 1,...,k
*
*  =====================================================================

      implicit none
      integer          :: n, k, i
      double precision :: x(n,k), v(k), xv(k), lse(n), z(n,k)

      call logsumexp(x, n, k, v, lse)

      do i = 1, n
        xv = x(i,:) + v
        z(i,:) = exp(xv - lse(i))
      end do

      return
      end

*  =====================================================================
      subroutine countf(x, n, bins, m, freq)
*  =====================================================================
C
C     Computes the counts (frequencies) of values in vector x by
C     assigning them to bins specified in vector bins.
C
C     Input:
C       x     - Integer vector of values
C       n     - Length of vector x
C       bins  - Integer vector of bin values
C       m     - Length of vector bins
C
C     Output:
C       freq  - Integer vector of frequencies for each bin value
C
      implicit none
      integer          :: n, m, i, j
      integer          :: x(n), bins(m), freq(m)

C     Initialize frequencies to zero
      do j = 1, m
         freq(j) = 0
      end do
      
C     Compute frequencies
      do i = 1, n
         do j = 1, m
            if (x(i) .eq. bins(j)) then
               freq(j) = freq(j) + 1
            end if
         end do
      end do

      return
      end

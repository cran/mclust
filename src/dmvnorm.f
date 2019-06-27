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

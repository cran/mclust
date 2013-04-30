C modified to avoid printing for calls from Fortran within R
      double precision function dgamma (x)
c jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
c jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog,
     2  dsin, dsqrt, xsml
C     external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt,
C    1  initds
      external d1mach, d9lgmc, dcsevl
c
c series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
      data gam cs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gam cs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gam cs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gam cs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gam cs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gam cs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gam cs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gam cs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gam cs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gam cs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gam cs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gam cs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gam cs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gam cs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gam cs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gam cs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gam cs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gam cs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gam cs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gam cs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gam cs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gam cs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gam cs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gam cs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gam cs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gam cs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gam cs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gam cs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gam cs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gam cs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gam cs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gam cs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gam cs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gam cs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gam cs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gam cs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gam cs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gam cs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gam cs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gam cs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gam cs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gam cs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /
c
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
c
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
c
      call d9gaml (xmin, xmax)
      xsml = exp (max (log(d1mach(1)), -log(d1mach(2)))+0.01d0)
      dxrel = sqrt (d1mach(4))
c
 10   y = abs(x)
      if (y.gt.10.d0) go to 50
c
c compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.0
c
      n = -n
C     if (x.eq.0.d0) call seteru (14hdgamma  x is 0, 14, 4, 2)
      if (x.eq.0.d0) dgamma = d1mach(2)
      if (x.eq.0.d0) return
C     if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) call seteru (
C    1  31hdgamma  x is a negative integer, 31, 4, 2)
      if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) 
     1  dgamma = -d1mach(2)
      if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) return
C     if (x.lt.(-0.5d0) .and. dabs((x-dint(x-0.5d0))/x).lt.dxrel) call
C    1  seteru (68hdgamma  answer lt half precision because x too near n
C    2egative integer, 68, 1, 1)
C     if (y.lt.xsml) call seteru (
C    1  54hdgamma  x is so close to 0.0 that the result overflows,
C    2  54, 5, 2)
      if (y.lt.xsml) dgamma = d1mach(2)
      if (y.lt.xsml) return
c
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
c
c gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
c
c gamma(x) for dabs(x) .gt. 10.0.  recall y = dabs(x).
c
C50   if (x.gt.xmax) call seteru (32hdgamma  x so big gamma overflows,
C    1  32, 3, 2)
 50   if (x.gt.xmax) dgamma = d1mach(2)
      if (x.gt.xmax) return
c
      dgamma = 0.d0
C     if (x.lt.xmin) call seteru (35hdgamma  x so small gamma underflows
C    1  , 35, 2, 0)
      if (x.lt.xmin) return
c
      dgamma = exp ((y-0.5d0)*log(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
C     if (dabs((x-dint(x-0.5d0))/x).lt.dxrel) call seteru (
C    1  61hdgamma  answer lt half precision, x too near negative integer
C    2  , 61, 1, 1)
c
      sinpiy = sin (pi*y)
C     if (sinpiy.eq.0.d0) call seteru (
C    1  31hdgamma  x is a negative integer, 31, 4, 2)
      if (sinpiy.eq.0.d0) dgamma = -d1mach(2)
      if (sinpiy.eq.0.d0) return
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end
C modified to omit priniting for calls from Fortran within R
      subroutine d9gaml (xmin, xmax)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c calculate the minimum and maximum legal bounds for x in gamma(x).
c xmin and xmax are not the only bounds, but they are the only non-
c trivial ones to calculate.
c
c             output arguments --
c xmin   dble prec minimum legal value of x in gamma(x).  any smaller
c        value of x might result in underflow.
c xmax   dble prec maximum legal value of x in gamma(x).  any larger
c        value of x might cause overflow.
c
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach,
     1  dlog
C     external d1mach, dlog
      external d1mach
c
      alnsml = log(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln+0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
C     call seteru (27hd9gaml  unable to find xmin, 27, 1, 2)
      xmin =  d1mach(2)
      xmax = -d1mach(2)
      return
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = log (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln-0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
C     call seteru (27hd9gaml  unable to find xmax, 27, 2, 2)
      xmin =  d1mach(2)
      xmax = -d1mach(2)
      return
c
 40   xmax = xmax - 0.01d0
      xmin = dmax1 (xmin, -xmax+1.d0)
c
      return
      end

      double precision function dcsevl (x, a, n)

      double precision a(n), x, twox, b0, b1, b2

      double precision d1mach
      external         d1mach
c
C     if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2)
      if (n.lt.1) dcsevl = -d1mach(2)
      if (n.lt.1) return
C     if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000,
C    1  31, 3, 2)
      if (n.gt.1000) dcsevl = d1mach(2)
      if (n.gt.1000) return
C     if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
C    1  25hdcsevl  x outside (-1,+1), 25, 1, 1)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) dcsevl = d1mach(2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) return

C added by CF to avoid uninitialized warnings
      b2 = 0
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end

      double precision function d9lgmc (x)

      double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach

      external d1mach, dcsevl, initds
c
      data algmcs(  1) / +.1666389480 4518632472 0572965082 2 d+0      /
      data algmcs(  2) / -.1384948176 0675638407 3298605913 5 d-4      /
      data algmcs(  3) / +.9810825646 9247294261 5717154748 7 d-8      /
      data algmcs(  4) / -.1809129475 5724941942 6330626671 9 d-10     /
      data algmcs(  5) / +.6221098041 8926052271 2601554341 6 d-13     /
      data algmcs(  6) / -.3399615005 4177219443 0333059966 6 d-15     /
      data algmcs(  7) / +.2683181998 4826987489 5753884666 6 d-17     /
      data algmcs(  8) / -.2868042435 3346432841 4462239999 9 d-19     /
      data algmcs(  9) / +.3962837061 0464348036 7930666666 6 d-21     /
      data algmcs( 10) / -.6831888753 9857668701 1199999999 9 d-23     /
      data algmcs( 11) / +.1429227355 9424981475 7333333333 3 d-24     /
      data algmcs( 12) / -.3547598158 1010705471 9999999999 9 d-26     /
      data algmcs( 13) / +.1025680058 0104709120 0000000000 0 d-27     /
      data algmcs( 14) / -.3401102254 3167487999 9999999999 9 d-29     /
      data algmcs( 15) / +.1276642195 6300629333 3333333333 3 d-30     /
c
      data nalgm, xbig, xmax / 0, 2*0.d0 /
c
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, 15, sngl(d1mach(3)) )
      xbig = 1.0d0/sqrt(d1mach(3))
      xmax = exp (dmin1(log(d1mach(2)/12.d0), -log(12.d0*d1mach(1))))
c
C10   if (x.lt.10.d0) 
 10   if (x.lt.10.d0) d9lgmc = d1mach(2)
      if (x.lt.10.d0) return
      if (x.ge.xmax) go to 20
c
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs,
     1  nalgm) / x
      return
c
 20   d9lgmc = 0.d0
C     call seteru (34hd9lgmc  x so big d9lgmc underflows, 34, 2, 0)
      return
c
      end

      double precision function dlngam (x)

      double precision x, dxrel, pi, sinpiy, sqpi2l, sq2pil,
     1  y, xmax, d9lgmc, d1mach
C    1  y, xmax, dgamma, d9lgmc, d1mach
      external d1mach, d9lgmc

c
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
c sq2pil = alog (sqrt(2*pi)),  sqpi2l = alog(sqrt(pi/2))
      data sqpi2l / +.2257913526 4472743236 3097614947 441 d+0    /

      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c
      data xmax, dxrel / 2*0.d0 /
c
C  added by CF to avoid uninitialized warnings
      dlngam = 0.d0

      if (xmax.ne.0.d0) go to 10
      xmax = d1mach(2)/dlog(d1mach(2))
      dxrel = dsqrt (d1mach(4))
c
 10   y = abs (x)
      if (y.gt.10.d0) go to 20
c
c dlog (dabs (dgamma(x)) ) for dabs(x) .le. 10.0
c
      dlngam = log (abs (dgamma(x)) )
      return
c
c dlog ( dabs (dgamma(x)) ) for dabs(x) .gt. 10.0
c
C20   if (y.gt.xmax) call seteru (
C    1  39hdlngam  dabs(x) so big dlngam overflows, 39, 2, 2)
 20   if (y.gt.xmax) dlngam = d1mach(2)
      if (y.gt.xmax) return
c
      if (x.gt.0.d0) dlngam = sq2pil + (x-0.5d0)*log(x) - x + d9lgmc(y)
      if (x.gt.0.d0) return
c
      sinpiy = abs (sin(pi*y))
C     if (sinpiy.eq.0.d0) call seteru (
C    1  31hdlngam  x is a negative integer, 31, 3, 2)
      if (sinpiy.eq.0.d0) dlngam = -d1mach(2)
      if (sinpiy.eq.0.d0) return
c
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
c
C     if (dabs((x-dint(x-0.5d0))*dlngam/x).lt.dxrel) call seteru (
C    1  68hdlngam  answer lt half precision because x too near negative
C    2integer, 68, 1, 1)
      return
c
      end

      function initds (dos, nos, eta)

      double precision dos(nos)

      integer          i1mach
      external         i1mach
c
C     if (nos.lt.1) call seteru (
C    1  35hinitds  number of coefficients lt 1, 35, 2, 2)
      if (nos.lt.1) initds = i1mach(9) 
c

C  added by CF to avoid uninitialized warnings
      i = 0

      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
C20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
C    1  1, 2)
 20   continue

      initds = i
c
      return
      end

      subroutine absrng( l, v, i, vmin, vmax)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      double precision v(*)

      integer          i, j, k, l

      double precision temp, vmin, vmax

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

      SUBROUTINE D2NORM ( N, X, INCX, VALUE )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), VALUE
*     ..
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*     THIS FUNCTION MODELLED AFTER DNRM2 BUT WRITTEN AS A SUBROUTINE
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      VALUE = NORM
      RETURN
*
*     End of D2NORM.
*
      END

      subroutine mclrup( l, n, v, r, lr)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer           l, n, lr

      double precision  cs, sn 

c     double precision  v(n), r(lr,n)
      double precision  v(*), r(lr,*)

      integer           i, j, k, m

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

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

      double precision   ss

c     double precision   x(n,p), u(p)
      double precision   x(n,*), u(*)

      double precision   ddot
      external           ddot

      integer                 i, j

      double precision        fac

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

c form mean
      fac = one / sqrt(dble(n))
      call dcopy( p, zero, 0, u, 1)
      do i = 1, n
        call daxpy( p, fac, x(i,1), n, u, 1)
      end do

c subtract mean and form sum of squares
      ss = zero
      do j = 1, p
        call daxpy( n, (-fac), u(j), 0, x(1,j), 1)
        ss = ss + ddot(n, x(1,j), 1, x(1,j), 1)
      end do

      return
      end

      subroutine mclvol( x, n, p, u, v, w,
     *                   work, lwork, iwork, liwork, 
     *                   info)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, lwork, liwork, info

c     integer            iwork(liwork)
      integer            iwork(*)

c     double precision   x(n,p), u(p), v(p,p), w(p,p), work(lwork),
      double precision   x(n,*), u(*), v(p,*), w(p,p), work(*)

      integer                 i, j

      double precision        temp, dummy, cmin, cmax

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157D+308)

c------------------------------------------------------------------------------

c form mean
      temp = one / dble(n)
      call dcopy( p, zero, 0, u, 1)
      do i = 1, n
        call daxpy( p, temp, x(i,1), n, u, 1)
      end do

c subtract mean
      do j = 1, p
        call daxpy( n, (-one), u(j), 0, x(1,j), 1)
      end do

c     if (.false.) then
c this gets the eigenvectors but x is overwritten

c get right singular vectors
c       call dgesvd( 'N', 'A', n, p, x, n, u, 
c    *                dummy, 1, w, p, work, lwork, info)

c       if (info .lt. 0) return

c       if (info .eq. 0) then
c         lwork = int(work(1))
c         do i = 1, p
c           v(i,i) = w(i,i)
c           if (i .gt. 1) then
c             do j = 1, (i-1)
c               v(i,j) = w(j,i)
c               v(j,i) = w(i,j)
c             end do
c           end if
c         end do
c         goto 100
c       end if

c     end if

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

c     EPSMAX = d1mach(4)

      call dsyevx( 'V', 'A', 'U', p, w, p, dummy, dummy, i, i,
     *              sqrt(EPSMAX), j, u, v, p,
     *              work, lwork, iwork(p+1), iwork, info)
          
      if (info .ne. 0) return

      lwork  = int(work(1))
      liwork = -1 

 100  continue
  
c      FLMAX = d1mach(2)

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

      subroutine sgnrng( l, v, i, vmin, vmax)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      double precision v(*)

      integer          i, j, k, l

      double precision temp, vmin, vmax

c----------------------------------------------------------------------------

      temp  = v(1)
      vmin  = temp
      vmax  = temp

      if (l .eq. 1) return

      if (i .eq. 1) then
        do j = 2, l
          temp = v(j)
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
        end do
      else
        k = 1 + i
        do j = 2, l
          temp = v(k)
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
          k    = k + i
        end do
      end if

      return  
      end 

      subroutine shapeo( TRANSP, s, O, l, m, w, info)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            TRANSP

      integer            l, m, info

c     double precision   s(l), O(l,l,m), w(l,l)
      double precision   s(*), O(l,l,*), w(l,*)

      integer                 j, k

      double precision        temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

      if (TRANSP) then
 
        do j = 1, l
          temp = sqrt(s(j))
          do k = 1, m
            call dscal( l, temp, O(j,1,k), l)
          end do
        end do

        do k = 1, m
          call dsyrk( 'U', 'T', l, l, one, O(1,1,k), l, zero, w, l)
          do j = 1, l
            call dcopy( j, w(1,j), 1, O(1,j,k), 1)
          end do
          do j = 2, l
            call dcopy( j-1, w(1,j), 1, O(j,1,k), l)
          end do
        end do

        info = 0
        return
      end if

      if (.not. TRANSP) then

        do j = 1, l
          temp = sqrt(s(j))
          do k = 1, m
            call dscal( l, temp, O(1,j,k), 1)
          end do
        end do

        do k = 1, m
          call dsyrk( 'U', 'N', l, l, one, O(1,1,k), l, zero, w, l)
          do j = 1, l
            call dcopy( j, w(1,j), 1, O(1,j,k), 1)
          end do
          do j = 2, l
            call dcopy( j-1, w(1,j), 1, O(j,1,k), l)
          end do
        end do

        info = 0
        return
      end if

      info = -1

      return
      end

      subroutine unchol ( UPPER, T, l, n, info)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            UPPER

      integer            l, n, info

c     double precision   T(abs(n), abs(n))
      double precision   T(  l   ,   *   )

       integer            i, j, k

       external           ddot
       double precision   ddot

c------------------------------------------------------------------------------

      if (UPPER) then

        do i = 2, n
          do j = 1, (i-1)
            T(i,j) = ddot( j, T(1,i), 1, T(1,j), 1)
          end do
        end do

        do k = 1, n
          T(k,k) = ddot( k, T(1,k), 1, T(1,k), 1)
        end do

        do k = 1, n-1
          call dcopy( n-k, T(k+1,k), 1, T(k,k+1), l)
        end do

        info = 0
        return
      end if

      if (.not. UPPER) then

        do i = 2, n
          do j = 1, (i-1)
            T(j,i) = ddot( j, T(i,1), l, T(j,1), l)
          end do
        end do
        
        do k = 1, n
          T(k,k) = ddot( k, T(k,1), l, T(k,1), l)
        end do

        do k = 2, n
          call dcopy( k-1, T(1,k), 1, T(k,1), l)
        end do

        return
        end if

      info = -1
      return
      end

      subroutine wardsw( i, n, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             i, n

      double precision    d(*)

      integer             i1, n1, ii, nn, k

      double precision    temp

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157D+308)

*-----------------------------------------------------------------------------

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

      subroutine es1e ( x, mu, sigsq, pro, n, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, G

      double precision   sigsq, hood, Vinv

c     double precision   x(n), mu(G), pro(G[+1]), z(n,G[+1])
      double precision   x(*), mu(*), pro(  *  ), z(n,  *  )

      integer                 i, k, nz

      double precision        temp, const, muk, prok, tmin, tmax, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (sigsq .le. zero) then
        hood  = FLMAX
        return
      end if

      const = pi2log + log(sigsq)

      do k = 1, G
        muk  = mu(k)
c       prok = pro(k)
        do i = 1, n
          temp   = x(i) - muk
c         z(i,k) = prok*exp(-(const+(temp*temp)/sigsq)/two)
          if (sigsq .lt. one .and. 
     *        abs(temp) .ge. sqrt(sigsq)*RTMAX) then
            hood = FLMAX
            return
          end if 
          z(i,k) = -(const+(temp*temp)/sigsq)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       temp = zero
c       do k = 1, nz
c         temp = temp + z(i,k)
c       end do
c       hood = hood + log(temp)
c       call dscal( nz, (one/temp), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if   
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine hc1e  ( x, n, ic, ng, ns, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, ic(n), ng, ns, nd

c     double precision    x(n), d(ng*(ng-1)/2)
      double precision    x(*), d(*)

      integer             lg, ld, ll, lo, ls
      integer             i, j, k, m 
      integer             ni, nj, nij, iopt, jopt, iold, jold
      integer             ij, ici, icj, ii, ik, jk

      double precision    ri, rj, rij, si, sj, sij
      double precision    temp, dij, dopt, dold

      external            wardsw

      double precision    one
      parameter          (one = 1.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      iopt   = 0
      jopt   = 0

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
          k     = k + 1
c         call dswap( p, x(k,1), n, x(j,1), n)
          temp  = x(k)
          x(k)  = x(j)
          x(j)  = temp
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
c         call dscal( p, si, x(i,1), n)
c         call daxpy( p, sj, x(j,1), n, x(i,1), n)
          x(i) = si*x(i) + sj*x(j)
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
c         call dcopy( p, x(i,1), n, v, 1)
c         call dscal( p, sj, v, 1)
c         call daxpy( p, (-si), x(j,1), n, v, 1)
c         dij   = ddot(p, v, 1, v, 1)
          temp  = sj*x(i) - si*x(j)
          dij   = temp*temp
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
          x(1)  = dble(iopt)
          ic(1) = jopt
        else
          x(1)  = dble(jopt)
          ic(1) = iopt
        end if
        d(1)    = dopt
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
        call wardsw( jopt, lg, d)
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
        si = d(i)
        if (si .le. dopt) then
          ij   = i
          dopt = si
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
        x(1)  = si 
        d(ld) = sj
      else
        x(1)  = sj
        d(ld) = si
      end if
      ld = ld - 1

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
          x(k)  = dble(ici)
          d(ld) = dble(icj)
        else
          x(k)  = dble(icj)
          d(ld) = dble(ici)
        end if
        ld = ld - 1
      end do

      ld = nd
      lo = nd - 1
      do k = 1, ns
        ic(k) = int(d(lo))
        lo    = lo - 1
        ld    = ld - 1
        d(ld) = d(lo)
        lo    = lo - 1
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

      subroutine me1e ( EQPRO, x, n, G, Vinv, z, maxi, tol, eps, 
     *                  mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, G, maxi

      double precision    Vinv, eps, tol

c     double precision    x(n), z(n,G[+1]), mu(G), sigsq, pro(G[+1])
      double precision    x(*), z(n,  *  ), mu(*), sigsq, pro(  *  )

      integer                 nz, iter, k, i

      double precision        hold, hood, err, prok, tmin, tmax, ViLog
      double precision        const, sum, sumz, smu, temp, term, zsum
      double precision        rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)
      end if
 
      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX 
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter  + 1

      sumz  = zero
      sigsq = zero

      zsum  = one

      do k = 1, G
        sum = zero
        smu = zero
        do i = 1, n
          temp = z(i,k)   
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
        sumz   = sumz + sum
        if (.not. EQPRO) pro(k) = sum / dble(n)
        zsum   = min(sum,zsum)
        if (sum .gt. rteps) then
          smu    = smu / sum
          mu(k)  = smu
          do i = 1, n
            temp   = x(i) - smu
            temp   = temp*temp
            sigsq  = sigsq + z(i,k)*temp
            z(i,k) = temp
          end do
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      if (Vinv .le. zero) then
        sigsq  = sigsq / dble(n)
      else
        sigsq  = sigsq / sumz
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
         
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      const = pi2log + log(sigsq)

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
c         z(i,k) = temp*exp(-(const+(z(i,k)/sigsq))/two)           
          z(i,k) = -(const+(z(i,k)/sigsq))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine me1ep ( EQPRO, x, n, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, G, maxi

      double precision    pshrnk, pmu, pscale, pdof

      double precision    Vinv, eps, tol

c     double precision    x(n), z(n,G[+1]), mu(G), sigsq, pro(G[+1])
      double precision    x(*), z(n,  *  ), mu(*), sigsq, pro(  *  )

      integer                 nz, iter, k, i

      double precision        hold, hood, err, prok, tmin, tmax, ViLog
      double precision        const, sum, sumz, smu, temp, term, zsum
      double precision        pmupmu, cgam, cmu, rmu, rgam, rteps
      
      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      double precision        dlngam
      external                dlngam

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)
      end if
 
      eps    = max(eps,zero)
      tol    = max(tol,zero)

      rteps  = sqrt(eps)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX 
      err    = FLMAX

      pmupmu = pmu*pmu

      iter   = 0

100   continue

      iter  = iter  + 1

      sigsq = zero

      zsum  = one

      do k = 1, G
        sumz = zero
        smu  = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          smu  = smu + temp*x(i)
        end do
        if (.not. EQPRO) pro(k)   = sumz / dble(n)
        zsum = min(zsum,sumz)
        if (sumz .gt. rteps) then
          smu = smu/sumz
          sum = zero
          do i = 1, n
            term   = x(i) - smu
            term   = term*term
            sum    = sum + z(i,k)*term
          end do
          term  = (pshrnk*sumz)/(pshrnk+sumz)
          temp  = (pmupmu + smu*smu) - two*pmu*smu
          sigsq = sigsq + sum + term*temp
          term  = sumz/(pshrnk+sumz)
          temp  = pshrnk/(pshrnk+sumz)
          mu(k) = term*smu + temp*pmu
        end if
      end do 

      if (zsum .le. rteps) then
        tol  = zsum
        eps  = -FLMAX
        maxi = iter
        return
      end if

      sigsq = (pscale + sigsq)/(pdof + dble(n+G) + two)

c     if (Vinv .le. zero) then
c       sigsq  = sigsq / dble(n)
c     else
c       sigsq  = sigsq / sumz
c     end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
         
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      const = pi2log + log(sigsq)

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
          term = x(i) - mu(k)
          z(i,k) = -(const+((term*term)/sigsq))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      cmu  = dble(G)*(pi2log-log(pshrnk))/two

      sum  = zero
      do k = 1, G
        temp = pmu - mu(k)
        temp = temp*temp
        sum  = sum - (pshrnk/sigsq)*temp
      end do

      term = log(sigsq)
      rmu  = (sum - dble(G)*term)/two

      temp =  pdof/two
      cgam =  temp*log(pscale/two) - dlngam(temp)

      rgam = -(temp+one)*term - (pscale/sigsq)/two

      pdof = (cmu+cgam) + (rmu+rgam)

      return
      end

      subroutine ms1e ( x, z, n, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, G

c     double precision    x(n), z(n,G), mu(G), sigsq, pro(G)
      double precision    x(*), z(n,*), mu(*), sigsq, pro(*)

      integer                 i, k

      double precision        sum, smu, sumz, temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      sumz  = zero
 
      sigsq = zero

      do k = 1, G
        sum = zero
        smu = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
        sumz   = sumz + sum
        pro(k) = sum / dble(n)
        if (sigsq .gt. one .or. smu .le. sum*FLMAX) then
          smu    = smu  / sum
           mu(k)  = smu 
          if (sigsq .ne. FLMAX) then 
            do i = 1, n
              temp = abs(x(i) - smu)
              sigsq  = sigsq + z(i,k)*(temp*temp)
            end do
          end if
        else 
          mu(k) = FLMAX 
          sigsq = FLMAX
        end if
      end do

c sumz .eq. n when no noise
      if (sigsq .ne. FLMAX) sigsq  = sigsq / sumz

      return
      end

      subroutine ms1ep ( x, z, n, G,
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, G

      double precision    pshrnk, pmu, pscale, pdof

c     double precision    x(n), z(n,G), mu(G), sigsq, pro(G)
      double precision    x(*), z(n,*), mu(*), sigsq, pro(*)

      integer                 k, i

      double precision        pmupmu
      double precision        sum, sumz, smu, temp, term
      
      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      pmupmu = pmu*pmu 

      sigsq = zero

      do k = 1, G
        sumz = zero
        smu  = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          smu  = smu + temp*x(i)
        end do
        pro(k)   = sumz / dble(n)
        if (sumz .gt. one .or. smu .lt. sumz*FLMAX) then
          smu = smu/sumz
          sum = zero
          term  = sumz/(pshrnk+sumz)
          temp  = pshrnk/(pshrnk+sumz)
          mu(k) = term*smu + temp*pmu
          if (sigsq .ne. FLMAX) then          
            do i = 1, n
              term = abs(x(i) - smu)
              sum = sum + z(i,k)*(term*term)
            end do
            term  = (pshrnk*sumz)/(pshrnk+sumz)
            temp  = (pmupmu + smu*smu) - two*pmu*smu
            sigsq = sigsq + sum + term*temp
          end if
        else
          mu(k) = FLMAX
          sigsq = FLMAX
        end if
      end do 

      if (sigsq .ne. FLMAX) then
        temp = pdof + dble(n) + two
        if (pshrnk .gt. zero) temp = temp + dble(G)
        sigsq = (pscale + sigsq)/temp
      end if

      return
      end

      subroutine eseee ( CHOL, x, mu, Sigma, pro, n, p, G, Vinv,
     *                   w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     character          CHOL
      logical            CHOL

c     integer            n, p, G
      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), w(p), z(n,G[+1])
      double precision   x(n,*), w(*), z(n,  *  )

c     double precision   mu(p,G), Sigma(p,p), pro(G[+1])
      double precision   mu(p,*), Sigma(p,*), pro(  *  )

      integer                 info, i, j, k, nz

      double precision        rteps, detlog, prok, tmin, tmax
      double precision        umin, umax, const, temp, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------


c     if (CHOL .eq. 'N') then
      if (.not. CHOL) then

c Cholesky factorization
        call dpotrf( 'U', p, Sigma, p, info)

        if (info .ne. 0) then
          w(1) = dble(info)
          hood = FLMAX
          return
        end if

      end if

      call absrng( p, Sigma, (p+1), umin, umax)

c     rc   = umin/(one+umax)

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

      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
c         z(i,k) = prok*exp(-(const+temp))
          z(i,k) = -(const+temp)
        end do
      end do

      w(1) = zero
      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          w(1) = zero
          hood = FLMAX
          return
        end if 
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      w(1) = zero
 
      return
      end

      double precision function detmc2( n, u)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer                          k, n
 
      double precision                 q

      double precision                 u(n,*)

      double precision                 zero, two
      parameter                       (zero = 0.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

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

      subroutine meeee ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps,
     *                   mu, U, pro, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p), pro(G)
      double precision   mu(p,*), U(p,*), pro(*)

      integer                 nz, p1, iter, i, j, k, j1

      double precision        piterm, sclfac, sumz, sum, zsum
      double precision        cs, sn, umin, umax, rc, detlog, rteps
      double precision        const, hold, hood, err, temp, term
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

c zero out the lower triangle
      i = 1
      do j = 2, p
        call dcopy( p-i, zero, 0, U(j,i), 1)
        i = j
      end do

      iter = 0

100   continue

      iter = iter + 1

      do j = 1, p
        call dcopy( j, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      zsum = one
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum / dble(n)
        zsum = min(zsum,sum)
        if (sum .gt. rteps) then
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
        end if
      end do
 
      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        hood = eps  
        maxi =  iter
        return
      end if

      if (Vinv .le. zero) then
        sclfac = one/sqrt(dble(n))
      else
        sclfac = one/sqrt(sumz)
      end if

      do j = 1, p
        call dscal( j, sclfac, U(1,j), 1)
      end do

c condition number

      call absrng( p, U, p1, umin, umax)

      rc = umin/(one+umax)

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if   

      end if

      if (rc .le. rteps) then
        tol  = err
        eps  = FLMAX
        hood = eps
        maxi = iter
        return
      end if

      detlog = zero
      do j = 1, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const = piterm + detlog

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp * exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine meeeep( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps,
     *                   mu, U, pro, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p), pro(G)
      double precision   mu(p,*), U(p,*), pro(*)

      integer                 nz, p1, iter, i, j, k, j1

      double precision        piterm, sclfac, sumz, sum, zsum
      double precision        cs, sn, umin, umax, rc, detlog, rteps
      double precision        const, hold, hood, err, temp, term
      double precision        prok, tmin, tmax, ViLog
      double precision        cmu, cgam, rmu, rgam

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        twolog
      parameter              (twolog  = 0.6931471805599453d0)

      double precision        pilog
      parameter              (pilog = 1.144729885849400d0)

      double precision        pi2log      
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot, dlngam
      double precision        ddot, dlngam

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero 

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      sclfac = one/sqrt(dble(n))

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      iter   = iter + 1

c copy pscale to U
      do j = 1, p
        call dcopy( p, pscale(1,j), 1, U(1,j), 1)
      end do

      sumz = zero
      zsum = one
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum / dble(n)
        zsum = min(zsum,sum)
        if (sum .gt. rteps) then     
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
          call dcopy( p, pmu, 1, w, 1)  
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          const = sum + pshrnk
          temp  = (sum*pshrnk)/const
          call dscal( p, sqrt(temp), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
          call dscal( p, sum/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = pdof + dble(p) + one
      if (pshrnk .gt. zero) term = term + dble(G)

      if (Vinv .le. zero) then
        sclfac = one/sqrt(term+dble(n))
      else 
        sclfac = one/sqrt(term+dble(sumz))
      end if

      do j = 1, p
        call dscal( j, sclfac, U(1,j), 1)
      end do

      if (Vinv .gt. zero) then
        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if   
      end if

c condition number

      call absrng( p, U, p1, umin, umax)

      rc = umin/(one+umax)

      if (rc .le. rteps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      detlog = zero
      do j = 1, p
        detlog = detlog + log(abs(U(j,j)))
      end do

      const = piterm + detlog

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp * exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol   = err
      eps   = hood
      maxi  = iter

      if (pshrnk .gt. zero) then
        cmu = dble(p)*(log(pshrnk) - pi2log)/two

        rmu = zero
        do k = 1, G
          call daxpy( p, (-one), mu(1,k), 1, pmu, 1)
          call dtrsv('U','T','N',p,U,p,pmu,1)
          rmu = rmu + ddot( p, pmu, 1, pmu, 1)
        end do

        sum  = zero
        term = zero
        temp = zero
        do j = 1, p
          call dcopy( p, pscale(j,1), p, pmu, 1)
c         call dtrsv('U','T','N', p, U, p, pmu, 1)
          i = p-j+1
c         call dtrsv('U','T','N', i, U(j,j),i,pmu(j),1)
          call dtrsv('U','T','N', i, U(j,j),p,pmu(j),1)
          sum  = sum + ddot(i, pmu(j), 1, pmu(j), 1)
          temp = temp + log(abs(pscale(j,j)))
          term = term + dlngam((pdof+one-dble(j))/two)
        end do

        rmu   = -(detlog+pshrnk*rmu/two)

        const = -dble(p)*(pdof*twolog+(dble(p)-one)*pilog/two)
        cgam  = (const/two-pdof*temp) - term

        rgam  = -((pdof+dble(p)+one)*detlog + sum/two)

        pdof  = (dble(G)*cmu+rmu) + (cgam+rgam)

      else
        pdof  = FLMAX
      end if

      return
      end

      subroutine mseee ( x, z, n, p, G, w, mu, U, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p), pro(G)
      double precision   mu(p,*), U(p,*), pro(*)

c------------------------------------------------------------------------------
c
c  x       double    (input) (n,p) matrix of observations.
c  z       double    (input) (n,G) conditional probabilities. 
c  n       integer   (input) number of observations.
c  p       integer   (input) dimension of the data.
c  G       integer   (input) number of Gaussian clusters in the mixture.
c  w       double    (scratch) (p) 
c  mu      double    (output) (p,G) mean for each group.
c  U       double    (output) (p,p) upper triangular Cholesky factor of the 
c         common covariance matrix for the groups: transpose(U) * U = Sigma.
c  pro     double    (output) (G) mixing proportions (ignore result if equal).

      integer                 i, j, k, j1

      double precision        sum, sumz, zsum, temp, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      do j = 1, p
        call dcopy( p, zero, 0, U(1,j), 1)
      end do

      sumz = zero
      zsum = one

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum / dble(n)
        if (sum .gt. one .or. one .gt. sum*FLMAX) then
          zsum = min(zsum,sum)
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
        else 
          zsum = zero
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (zsum .eq. zero) return

c sumz .eq. n when no noise
      do j = 1, p
        call dscal( j, one/sqrt(sumz), U(1,j), 1)
      end do

      return
      end

      subroutine mseeep( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   w, mu, U, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

c     double precision   mu(p,G), U(p,p), pro(G)
      double precision   mu(p,*), U(p,*), pro(*)

      integer                 i, j, k, j1

      double precision        sclfac, const, temp
      double precision        sum, sumz, zsum, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      if (pshrnk .le. zero) pshrnk = zero
 
      do j = 1, p
        call dcopy( p, pscale(1,j), 1, U(1,j), 1)
      end do

      sumz = zero
      zsum = one

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum / dble(n)
        if (sum .ge. one .or. one .gt. sum*FLMAX) then
          zsum = min(zsum,sum)
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
          call dcopy( p, pmu, 1, w, 1)  
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          const = sum + pshrnk
          temp  = (sum*pshrnk)/const
          call dscal( p, sqrt(temp), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( U(j,j), w(j), cs, sn)
            call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( U(p,p), w(p), cs, sn)
          call dscal( p, sum/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        else
          zsum = zero 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (zsum .eq. zero) return

      temp   = pdof+dble(n+p+1)
      if (pshrnk .gt. zero) temp = temp + dble(G)
      sclfac = one/sqrt(temp)
      
      do j = 1, p
        call dscal( j, sclfac, U(1,j), 1)
      end do

      return
      end

      subroutine eseei ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   scale, hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), shape(p), pro(G[+1])
      double precision   mu(p,*), shape(*), pro(  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        smin, smax, prok

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (scale .le. zero) then
        hood = FLMAX
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do

      const = dble(p)*(pi2log+log(scale)) 


      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            if (abs(temp) .ge. shape(j)*FLMAX) then
              hood = FLMAX
              return
            end if
            temp = temp/shape(j)
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      if (pro(1) .lt. zero) return 

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp = z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine meeei ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, scale

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), shape(p), pro(G[+1])
      double precision    mu(p,*), shape(*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, zsum
      double precision    const, hold, hood, err, smin, smax
      double precision    prok, tmin, tmax, ViLog, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      zsum = one

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        sum = zero
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum/dble(n)
        zsum = min(zsum,sum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum = sum + temp*temp
            end do
            shape(j) = shape(j) + sum
          end do
        end if 
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) then
        scale = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        scale = FLMAX
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if 
 
      if (temp .gt. SMALOG) then 
        temp = exp(temp)
      else
        temp = zero
      end if

      if (Vinv .le. zero) then
        scale = temp/dble(n)
      else 
        scale = temp/sumz
      end if  

      if (temp .le. eps) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps .or. scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      const = dble(p)*(pi2log+log(scale))

      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j)             
          end do        
c         z(i,k) = prok*exp(-(const+(sum/scale))/two)
          z(i,k) = -(const+(sum/scale))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine meeeip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

      double precision    Vinv, eps, tol, scale

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), shape(p), pro(G[+1])
      double precision    mu(p,*), shape(*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, zsum
      double precision    const, hold, hood, err, smin, smax
      double precision    prok, tmin, tmax, ViLog, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps) 

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      call dcopy( p, pscale, 0, shape, 1)

      sumz  = zero
      zsum  = one

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum/dble(n)
        zsum = min(zsum,sum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          term  = pshrnk + sum
          const = (pshrnk*sum)/term
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum = sum + (temp*temp)
            end do
            shape(j) = shape(j) + sum
            temp     = pmu(j) - mu(j,k)
            shape(j) = shape(j) + const*(temp*temp)
          end do
        end if
      end do

      if (zsum .le. rteps) then
        tol   =  zsum
        eps   = -FLMAX
        maxi  =  iter
        return
      end if

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) then
        scale = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        scale = FLMAX
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if 
 
      if (temp .gt. SMALOG) then 
        temp = exp(temp)
      else
        temp = zero
      end if

      term = pdof + one
      if (pshrnk .gt. zero) term = term + one

      if (Vinv .le. zero) then
        scale = temp/(term + dble(n))
      else 
        scale = temp/(term + sumz)
      end if

      if (temp .le. eps) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps .or. scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      const = dble(p)*(pi2log+log(scale))

      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j)             
          end do        
c         z(i,k) = prok*exp(-(const+(sum/scale))/two)
          z(i,k) = -(const+(sum/scale))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine mseei ( x, z, n, p, G, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G)
      double precision   x(n,*), z(n,*)

c     double precision   mu(p,G), scale, shape(p), pro(G)
      double precision   mu(p,*), scale, shape(*), pro(*)

      integer             i, j, k

      double precision    sum, sumz, temp, smin, smax

      double precision    zero, one
      parameter          (zero = 0.d0, one = 1.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      sumz = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        sumz   = sumz + sum
        pro(k) = sum/dble(n)
        if (sum .gt. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
         else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      call dcopy( p, zero, 0, shape, 1)

      do j = 1, p
        sum = zero
        do i = 1, n
          do k = 1, G
            if (mu(1,k) .eq. FLMAX) then
              scale =  FLMAX
              return
            end if
            temp = z(i,k)*(x(i,j) - mu(j,k))
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
        end do
        shape(j) = shape(j) + sum
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do      
      temp = sum/dble(p)
  
      if (temp .gt. BIGLOG) then
        scale = FLMAX
        call dcopy( p, FLMAX, 0, shape, 1)
        return
      end if 
  
      if (temp .gt. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if

      if (sumz .lt. one .and. temp .ge. sumz*FLMAX) then
        scale = FLMAX
        call dcopy( p, FLMAX, 0, shape, 1)
        return
      end if
  
      scale = temp/sumz

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        return
      end if

      call dscal( p, one/temp, shape, 1)
      
      return
      end

      subroutine mseeip( x, z, n, p, G,
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, G

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G)
      double precision    x(n,*), z(n,*)

c     double precision    mu(p,G), scale, shape(p), pro(G[+1])
      double precision    mu(p,*), scale, shape(*), pro(  *  )

      integer             i, j, k

      double precision    sum, sumz, temp, term
      double precision    const, smin, smax

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      call dcopy( p, pscale, 0, shape, 1)

      sumz  = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        sumz   = sumz + sum
        pro(k) = sum/dble(n)
        if (sum .gt. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          term  = pshrnk + sum
          const = (pshrnk*sum)/term
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
            end do
            shape(j) = shape(j) + sum
            temp     = pmu(j) - mu(j,k)
            shape(j) = shape(j) + const*(temp*temp)
          end do
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .ge. BIGLOG) then
        scale = FLMAX
        call dcopy( p, FLMAX, 0, shape, 1)
        return
      end if 

      if (temp .gt. SMALOG) then
        smin = exp(temp)
      else
        smin = zero
      end if 

      term  = pdof + sumz + two
      if (pshrnk .gt. zero) term = term + dble(G)
      scale = smin/term

      if (smin .lt. one .and. one .ge. smin*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        return
      end if

      call dscal( p, one/smin, shape, 1)

      return
      end

      subroutine eseev ( x, mu, scale, shape, O, pro, n, p, G, 
     *                   Vinv, v, w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p, G
      integer            n, p, G

      double precision   scale, Vinv, hood

c     double precision   x(n,p), v(p), w(p), z(n,G[+1])
      double precision   x(n,*), v(*), w(*), z(n,  *  )

c     double precision   mu(p,G), shape(p), O(p,p,G), pro(G[+1])
      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

      integer                 i, j, k, nz

      double precision        const, temp, tmin, tmax
      double precision        smin, smax, prok, eps, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (scale .le. zero) then
        hood = FLMAX
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)
  
      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do
        
      const = dble(p)*(pi2log + log(scale))

      do k = 1, G

c       prok = pro(k)

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, 
     *                 w, 1, zero, v, 1)
          do j = 1, p
            if (shape(j) .lt. one .and. 
     *          abs(v(j)) .ge. shape(j)*FLMAX) then
              hood = FLMAX
              return
            end if 
            v(j) = v(j)/shape(j)
          end do
          temp   = ddot( p, v, 1, v, 1)
c         z(i,k) = prok*exp(-(const+temp)/two)
          z(i,k) = -(const+temp)/two
        end do

      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine meeev ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps,
     *                   lwork, mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi, lwork

      double precision	 Vinv, eps, tol, scale

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

      integer                 nz, p1, iter, i, j, k, l, j1, info

      double precision        dnp, dummy, temp, term, rteps
      double precision        sumz, sum, smin, smax, cs, sn
      double precision        const, rc, hood, hold, err
      double precision        prok, tmin, tmax, ViLog, zsum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      p1     = p + 1

      dnp    = dble(n*p)

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      iter = iter + 1

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero
      zsum = one

      l = 0

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        zsum = min(zsum,sum)
        if (.not. EQPRO) pro(k) = sum / dble(n)
        if (sum .ge. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else 
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      if (l .ne. 0 .or. zsum .lt. rteps) then
        lwork = l        
c       w(1)  = FLMAX
        tol   = err
        if (l .ne. 0) then
          eps =  FLMAX
        else
          eps = -FLMAX
        end if
        maxi  = iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. zero) then
        lwork = 0
c       w(1)  = smin
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if
 
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp  = sum/dble(p)

      if (temp .gt. BIGLOG) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if 

      if (temp .gt. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if 
     
      if (Vinv .le. zero) then
        scale = temp/dble(n)
      else
        scale = temp/sumz
      end if

      if (temp .le. eps) then
        lwork = 0
c       w(1)  = temp
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call sgnrng( p, shape, 1, smin, smax)
      
      if (smin .le. eps) then
        lwork = 0
c       w(1)  = -smin
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if
      
      if (scale .le. eps) then
c       w(1)  = -scale
        lwork = 0
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        w(j) = temp*sqrt(shape(j))
      end do

      call absrng( p, w, 1, smin, smax)

      rc = smin / (one + smax)
      
      if (smin .le. rteps) then
c       w(1)  = -smin
        lwork = 0
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      const = dble(p)*(pi2log + log(scale))/two

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, w(p1), 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / w(j)
          end do
          sum    = ddot( p, s, 1, s, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      lwork = 0

c     w(1)  = rc

      tol   = err
      eps   = hood
      maxi  = iter

      return
      end

      subroutine meeevp( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps,
     *                   lwork, mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi, lwork

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   Vinv, eps, tol, scale

c     double precision   x(n,p), z(n,G[+1]), w(lwork), s(p)
      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

c     double precision   mu(p,G), shape(p), O(p,p,G), pro(G[+1])
      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

      integer                 nz, p1, iter, i, j, k, l, j1, info

      double precision        dnp, dummy, temp, term, rteps
      double precision        sumz, sum, smin, smax, cs, sn
      double precision        const, rc, hood, hold, err
      double precision        prok, tmin, tmax, ViLog, zsum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      p1    = p + 1

      dnp   = dble(n*p)

      eps   = max(eps,zero)
      rteps = sqrt(eps)

      tol   = max(tol,zero)

      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      call dcopy( p, zero, 0, shape, 1)

      zsum  = one
      sumz  = zero

      l     = 0

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum / dble(n)
        zsum = min(zsum,sum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          term  = sum+pshrnk
          const = (sum*pshrnk)/term
          call dscal( p, sqrt(const), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( O(j,j,k), w(j), cs, sn)
            call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( O(p,p,k), w(p), cs, sn)
          call dscal( p, sum/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else 
            do j = 1, p
              temp     = s(j)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      if (l .ne. 0 .or. zsum .le. rteps) then
        lwork = l        
c       w(1)  = FLMAX
        tol   = err
        if (l .ne. 0) then
          eps =  FLMAX
        else
          eps = -FLMAX
        end if
        maxi  = iter
        return
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
c       w(1)  = smin
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if
 
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if 

      if (temp .gt. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if 

      term = pdof + dble(p) + one
      if (pshrnk .gt. zero) term = term + one 
      if (Vinv .le. zero) then
        scale = temp/(term + dble(n))
      else
        scale = temp/(term + sumz)
      end if

      if (temp .le. eps) then
c       w(1)  = temp
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call sgnrng( p, shape, 1, smin, smax)
      
      if (smin .le. eps) then
c       w(1)  = -smin
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if
      
      if (scale .le. eps) then
c       w(1)  = -scale
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        w(j) = temp*sqrt(shape(j))
      end do

      call sgnrng( p, w, 1, smin, smax)

      rc = smin / (one + smax)
      
      if (smin .le. rteps) then
c       w(1)  = -smin
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      const = dble(p)*(pi2log + log(scale))/two

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, w(p1), 1, zero, s, 1)
          do j = 1, p
            s(j) = s(j) / w(j)
          end do
          sum    = ddot( p, s, 1, s, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      lwork = 0

c     w(1)  = rc

      tol   = err
      eps   = hood
      maxi  = iter

      return
      end
      subroutine mseev ( x, z, n, p, G, w, lwork, 
     *                   mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G, lwork

      double precision   scale

c     double precision   x(n,p), z(n,G), w(lwork)
      double precision   x(n,*), z(n,*), w(  *  )

c     double precision   shape(p), O(p,p,G), mu(p,G), pro(G)
      double precision   shape(*), O(p,p,*), mu(p,*), pro(*)

      integer                 i, j, k, j1, l, info

      double precision        dummy, sum, sumz, temp
      double precision        cs, sn, smin, smax

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        BIGLOG
      parameter              (BIGLOG =  709.d0)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call dcopy( p, zero, 0, shape, 1)

      l     = 0

      sumz  = zero
      scale = zero

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum / dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k), 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else if (scale .ne. FLMAX) then
            do j = 1, p
              temp     = z(j,k)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        else
          scale = FLMAX          
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (scale .eq. FLMAX  .or. l .ne. 0) then
        lwork = l       
        if (l .ne. 0) then
          scale =  FLMAX
        else
          scale = -FLMAX
        end if
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      lwork = 0

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = FLMAX
        return
      end if
 
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        scale = FLMAX
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if 

      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if 

      if (temp .ge. sumz*FLMAX) then
        scale = FLMAX 
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      scale = temp/sumz

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      call dscal( p, one/temp, shape, 1)

      return
      end

      subroutine mseevp( x, z, n, p, G,
     *                   pshrnk, pmu, pscale, pdof,
     *                   w, lwork, mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G, lwork

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   scale

c     double precision   x(n,p), z(n,G), w(lwork)
      double precision   x(n,*), z(n,*), w(  *  )

c     double precision   mu(p,G), shape(p), O(p,p,G), pro(G)
      double precision   mu(p,*), shape(*), O(p,p,*), pro(*)

      integer            p1, i, j, k, l, j1, info

      double precision   dummy, temp, term, const
      double precision   sumz, sum, smin, smax, cs, sn

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   FLMAX
      parameter         (FLMAX = 1.7976931348623157d308)

      double precision   BIGLOG
      parameter         (BIGLOG =  709.d0)

      double precision   SMALOG
      parameter         (SMALOG = -708.d0)

      external           ddot
      double precision   ddot

c------------------------------------------------------------------------------

      if (pshrnk .gt. zero) pshrnk = zero

      p1     = p + 1

      call dcopy( p, zero, 0, shape, 1)

      l      = 0

      sumz   = zero
      scale  = zero

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum / dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          term  = sum+pshrnk
          const = (sum*pshrnk)/term
          call dscal( p, sqrt(const), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( O(j,j,k), w(j), cs, sn)
            call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( O(p,p,k), w(p), cs, sn)
          call dscal( p, sum/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k), 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else if (scale .ne. FLMAX) then
            do j = 1, p
              temp     = z(j,k)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        else
          scale = FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (scale .eq. FLMAX  .or. l .ne. 0) then
        lwork = l        
        if (l .ne. 0) then
          scale =  FLMAX
        else
          scale = -FLMAX
        end if
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      lwork = 0

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = FLMAX
        return
      end if
 
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        scale = FLMAX
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if

      term = pdof + dble(p) + one
      if (pshrnk .gt. zero) term = term + one 
      scale = temp/(term + sumz)

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        return
      end if

      call dscal( p, one/temp, shape, 1)

      return
      end

      subroutine eseii ( x, mu, sigsq, pro, n, p, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   sigsq, hood, Vinv

c     double precision   x(n,p), mu(p,G), pro(G[+1]), z(n,G[+1])
      double precision   x(n,*), mu(p,*), pro(  *  ), z(n,  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, prok, tmin, tmax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (sigsq .le. zero) then
        hood = FLMAX
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum/sigsq)/two)
          if (sigsq .lt. one .and. sum .ge. sigsq*FLMAX) then
            hood = FLMAX
            return
          end if
          z(i,k) = -(const+sum/sigsq)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if 
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine hceii ( x, n, p, ic, ng, ns, v, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, ic(n), ng, ns, nd

c     double precision    x(n,p), v(p), d(ng*(ng-1)/2)
      double precision    x(n,*), v(*), d(*)

      integer             lg, ld, ll, lo, ls
      integer             i, j, k, m
      integer             ni, nj, nij, iopt, jopt, iold, jold
      integer             ij, ici, icj, ii, ik, jk

      double precision    ri, rj, rij, si, sj, sij
      double precision    dij, dopt, dold

      external            wardsw

      double precision    one
      parameter          (one = 1.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    ddot
      external            ddot

c------------------------------------------------------------------------------

      iopt   = 0

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
        call wardsw( jopt, lg, d)
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
        si = d(i)
        if (si .le. dopt) then
          ij   = i
          dopt = si
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

      subroutine meeii ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, sigsq

c     double precision    x(n,p), z(n,G[+1]), mu(p,G), pro(G[+1])
      double precision    x(n,*), z(n,  *  ), mu(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, prok, tmax, tmin, rteps
      double precision    const, hold, hood, err, dnp, ViLog, zsum

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG
      parameter          (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      dnp = dble(n*p)

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      sigsq = zero

      sumz  = zero
      zsum  = one

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum/dble(n)
        zsum = min(sum,zsum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = abs(x(i,j) - mu(j,k))
              if (temp .gt. RTMIN) sum  = sum + temp*temp
            end do
            if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)
     *        sigsq  = sigsq + z(i,k)*sum
            z(i,k) = sum
          end do
        else
          sigsq = FLMAX 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      if (Vinv .le. zero) then
        sigsq  = sigsq / dnp
      else 
        sigsq = sigsq / (dble(p)*sumz)
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      do k = 1, G
c       temp = pro(k)
        do i = 1, n
c         z(i,k) = temp*exp(-(const+(z(i,k)/sigsq))/two)
          z(i,k) = -(const+(z(i,k)/sigsq))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine meeiip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

      double precision    Vinv, eps, tol, sigsq

c     double precision    x(n,p), z(n,G[+1]), mu(p,G), pro(G[+1])
      double precision    x(n,*), z(n,  *  ), mu(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumk, sumz, temp, term, tmax, tmin
      double precision    const, hold, hood, err, dnp, ViLog, prok
      double precision    pmupmu, cmu, cgam, rmu, rgam, zsum, rteps 

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

      double precision    ddot, dlngam
      external            ddot, dlngam

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      dnp = dble(n*p)

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      rteps  = sqrt(eps)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

      pmupmu = ddot( p, pmu, 1, pmu, 1)

100   continue

      iter  = iter + 1

      sigsq = zero

      sumz = zero
      zsum = one

      do k = 1, G
        sumk = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumk = sumk + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sumk
        if (.not. EQPRO) pro(k) = sumk/dble(n)
        zsum = min(zsum,sumk)
        if (sumk .gt. rteps) then         
          call dscal( p, (one/sumk), mu(1,k), 1)
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = x(i,j) - mu(j,k)
              sum  = sum + temp*temp
            end do
            sigsq  = sigsq + z(i,k)*sum
          end do
          temp  = pmupmu + ddot( p, mu(1,k), 1, mu(1,k), 1)
          temp  = temp - two*ddot( p, mu(1,k), 1, pmu, 1)
          const = sumk+pshrnk
          sigsq = sigsq + ((pshrnk*sumk)/const)*temp
          call dscal( p, (sumk/const), mu(1,k), 1)
          call daxpy(p, (pshrnk/const), pmu, 1, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = zero
      if (Vinv .le. zero) then

        sigsq  = sigsq / (pdof + dble((n+G)*p) + two)

      else 

        sigsq = sigsq / (pdof + (sumz+dble(G))*dble(p) + two)

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      const = dble(p)*(pi2log+log(sigsq))

      do i = 1, n
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = -(const+(sum/sigsq))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMLOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      cmu   = dble(p)*(log(pshrnk)-pi2log)/two

      const = pdof/two
      cgam  = const*log(pscale/two)-dlngam(const)

      rmu   = zero
      do k = 1, G
        temp = pmupmu + ddot( p, mu(1,k), 1, mu(1,k), 1)
        temp = temp - two*ddot( p, mu(1,k), 1, pmu, 1)
        rmu  = rmu + (pshrnk*temp)/sigsq
      end do

      term = log(sigsq)
      rmu  = -(rmu + dble(p)*term)/two

      rgam = -(const+one)*term - (pscale/sigsq)/two

      pdof  = (dble(G)*cmu+cgam) + (rmu+rgam)

      return
      end

      subroutine mseii ( x, z, n, p, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), mu(p,G), sigsq, pro(G)
      double precision   x(n,*), z(n,*), mu(p,*), sigsq, pro(*)

      integer                 i, j, k

      double precision        sum, sumz, temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      sumz  = zero
      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum/dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then 
          call dscal( p, (one/sum), mu(1,k), 1)
          if (sigsq .ne. FLMAX) then
            do i = 1, n
              sum = zero
              do j = 1, p
                temp = abs(x(i,j) - mu(j,k))
                if (temp .gt. RTMIN) sum = sum + temp*temp
              end do
              if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)
     *          sigsq = sigsq + z(i,k)*sum
            end do
          end if
        else
          sigsq = FLMAX  
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if      
      end do

c sumz .eq. n when no noise
      if (sigsq .ne. FLMAX) sigsq = sigsq / (sumz*dble(p))

      return
      end

      subroutine mseiip( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, G

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G), mu(p,G), sigsq, pro(G)
      double precision    x(n,*), z(n,*), mu(p,*), sigsq, pro(*)

      integer             i, j, k

      double precision    sum, sumz, zsum, pmupmu
      double precision    const, temp, dnp

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    ddot
      external            ddot

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      dnp = dble(n*p)

      pmupmu = ddot( p, pmu, 1, pmu, 1)

      sumz = zero
      sigsq = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        pro(k) = sum/dble(n)
        if (sum .gt. one .or. one .le. sum*FLMAX) then         
          call dscal( p, (one/sum), mu(1,k), 1)
          temp  = pmupmu + ddot( p, mu(1,k), 1, mu(1,k), 1)
          temp  = temp - two*ddot( p, mu(1,k), 1, pmu, 1)
          const = sum+pshrnk
          call dscal( p, (sum/const), mu(1,k), 1)
          call daxpy(p, (pshrnk/const), pmu, 1, mu(1,k), 1)
          if (sigsq .ne. FLMAX) then
            sigsq = sigsq + ((pshrnk*sum)/const)*temp
            do i = 1, n
              sum = zero
              do j = 1, p
                temp = abs(x(i,j) - mu(j,k))
                if (temp .gt. RTMIN) sum  = sum + temp*temp
              end do
              if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)
     *            sigsq  = sigsq + z(i,k)*sum
            end do
          end if
        else
          sigsq = FLMAX  
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (sigsq .eq. FLMAX) return

      temp  = pdof + sumz*dble(p) + two
      if (pshrnk .gt. zero) temp = temp + dble(G*p)
      sigsq = sigsq / temp

      return
      end

      subroutine esevi ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   scale, hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), shape(p,G), pro(G[+1])
      double precision   mu(p,*), shape(p,*), pro(  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        smin, smax, prok

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (scale .le. zero) then
        hood = FLMAX
        return
      end if

      do k = 1, G
        call sgnrng( p, shape(1,k), 1, smin, smax)
        if (smin .eq. zero) then
          hood = FLMAX
          return
        end if
      end do

      temp = sqrt(scale)

      do k = 1, G
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
      end do

      const  = dble(p)*(pi2log+log(scale))

      do k = 1, G
c       prok   = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            if (shape(j,k) .lt. one .and. 
     *          abs(temp) .ge. shape(j,k)*FLMAX) then
              hood = FLMAX
              return
            end if
            temp = temp/shape(j,k)
            if (abs(temp) .ge. RTMAX) then
              hood = FLMAX
              return
            end if 
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine meevi ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, scale

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), shape(p,G), pro(G[+1])
      double precision    mu(p,*), shape(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, epsmin
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog, zsum, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      sumz  = zero
      zsum  = one
       
      do k = 1, G
        call dcopy( p, zero, 0, shape(1,k), 1)
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum /dble(n)
        zsum = min(sum,zsum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              if (sqrt(z(i,k))*abs(temp) .gt. RTMIN)
     *          sum = sum + z(i,k)*(temp*temp)
            end do
            shape(j,k) = shape(j,k) + sum
          end do
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol   =  zsum
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      scale  = zero
      epsmin = FLMAX
      do k = 1, G
        call sgnrng(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .gt. zero) then
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp = sum/dble(p)
          if (temp .gt. BIGLOG) then
            scale = FLMAX 
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
            tol   = err
            eps   = FLMAX
            maxi  = iter
            return
          end if
          if (temp .gt. SMALOG) then
            temp = exp(temp)
          else
            temp = zero
          end if
          scale  = scale + temp
          epsmin = min(temp,epsmin)
          if (temp .lt. eps) then
            scale = FLMAX 
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
            tol   = err
            eps   = FLMAX
            maxi  = iter
            return
          end if
          call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      term = zero
      if (Vinv .gt. zero) then
        scale = scale /sumz
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else
        scale = scale /dble(n)
      end if

      if (scale .le. eps) then
        tol   = epsmin
        eps   = FLMAX
        maxi  = iter
        return
      end if

      do k = 1, G

        call sgnrng( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = smin
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      const = dble(p)*(pi2log + log(scale))

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j,k)
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/scale))/two)
          z(i,k) = -(const+(sum/scale))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp = z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine meevip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

      double precision    Vinv, eps, tol, scale

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), shape(p,G), pro(G[+1])
      double precision    mu(p,*), shape(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, epsmin, zsum
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)
 
      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      sumz  = zero
      zsum  = one

      do k = 1, G
        call dcopy( p, pscale, 0, shape(1,k), 1)
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        if (.not. EQPRO) pro(k) = sum /dble(n)
        zsum = min(sum,zsum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          term  = pshrnk + sum
          const = (pshrnk*sum)/term
          do j = 1, p
            do i = 1, n
              temp       = x(i,j) - mu(j,k)
              if (abs(temp)*sqrt(z(i,k)) .gt. RTMIN)
     *          shape(j,k) = shape(j,k) + z(i,k)*(temp*temp)
            end do
            temp       = pmu(j) - mu(j,k)
            shape(j,k) = shape(j,k) + const*(temp*temp)
          end do
          call dscal( p, sum/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if 
      end do

      if (zsum .le. rteps) then
        tol   =  zsum
        eps   = -FLMAX
        maxi  =  iter
        return
      end if

      scale  = zero
      epsmin = FLMAX
      do k = 1, G
        call sgnrng(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .gt. zero) then
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp = sum/dble(p)
          if (temp .gt. BIGLOG) then
            scale = FLMAX 
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
            tol   = err
            eps   = FLMAX
            maxi  = iter
            return
          end if
          if (temp .gt. SMALOG) then
            temp = exp(temp)
          else
            temp = zero
          end if
          scale  = scale + temp
          epsmin = min(temp,epsmin)
          if (temp .lt. eps) then
            scale = FLMAX 
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
            tol   = err
            eps   = FLMAX
            maxi  = iter
            return
          end if
          call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      term = pdof + one
      if (Vinv .le. zero) then
        term = term + dble(n)
      else
        term = term + sumz
      end if 
      if (pshrnk .gt. zero) term = term + one

      scale  = scale/term

      if (Vinv .gt. zero) then
        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      end if

      if (scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call sgnrng( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = err
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      const = dble(p)*(pi2log + log(scale))

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j,k)
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/scale))/two)
          z(i,k) = -(const+(sum/scale))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine msevi ( x, z, n, p, G, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G)
      double precision   x(n,*), z(n,*)

c     double precision   mu(p,G), scale, shape(p,G), pro(G)
      double precision   mu(p,*), scale, shape(p,*), pro(*)

      integer                 i, j, k

      double precision        smin, smax
      double precision        sum, sumz, temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      scale = zero
      sumz  = zero

      do k = 1, G
        call dcopy( p, zero, 0, shape(1,k), 1)
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        sumz   = sumz + sum
        pro(k) = sum/dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
        else
          scale = FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if
      end do

      if (scale .eq. FLMAX) return

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      scale  = zero
     
      do k = 1, G

        call sgnrng(p, shape(1,k), 1, smin, smax)

        if (smin .le. zero) then
          scale = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100 
        end if

        sum   = zero
        do j = 1, p
          sum = sum + log(shape(j,k))
        end do      
        temp   = sum/dble(p)

        if (temp .ge. BIGLOG) then
          scale = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if
 
        if (temp .ge. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if
    
        if (scale .ne. FLMAX) scale  = scale + temp
     
        if (temp .lt. one .and. one .ge. temp*FLMAX) then
          scale = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        call dscal( p, one/temp , shape(1,k), 1)

100   continue
  
      end do
      
      if (sumz .lt. one .and. one .ge. sumz*FLMAX) then
        scale = FLMAX
        return 
      end if

      scale = scale/sumz

      return
      end

      subroutine msevip( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, G

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G)
      double precision    x(n,*), z(n,*)

c     double precision    mu(p,G), scale, shape(p,G), pro(G)
      double precision    mu(p,*), scale, shape(p,*), pro(*)

      integer             i, j, k

      double precision    sum, sumz, temp, term
            double precision    smin, smax, const

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      sumz = zero
      scale = zero

      do k = 1, G
        call dcopy( p, pscale, 0, shape(1,k), 1)
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz   = sumz + sum
        pro(k) = sum /dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          term  = pshrnk + sum
          const = (pshrnk*sum)/term
          do j = 1, p
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              if (abs(temp)*sqrt(z(i,k)) .gt. RTMIN)
     *          shape(j,k) = shape(j,k) + z(i,k)*(temp*temp)
            end do
            temp       = pmu(j) - mu(j,k)
            shape(j,k) = shape(j,k) + const*(temp*temp)
          end do
          call dscal( p, sum/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
        else
          scale = FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if 
      end do

      if (scale .eq. FLMAX) return

      scale  = zero
      do k = 1, G

        call sgnrng(p, shape(1,k), 1, smin, smax)

        if (smin .le. zero) then
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j,k))
        end do
        temp = sum/dble(p)

        if (temp .ge. BIGLOG) then
          scale = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        if (temp .ge. SMALOG) then 
          temp = exp(temp)
        else
          temp = zero
        endif

        if (scale .ne. FLMAX) scale  = scale + temp

        if (temp .le. one .and. one .ge. temp*FLMAX) then
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if
 
        call dscal( p, one/temp, shape(1,k), 1)

100   continue
 
      end do

      term = pdof + sumz + two
      if (pshrnk .gt. zero) term = term + dble(G)

      scale  = scale/term

      return
      end

      subroutine es1v  ( x, mu, sigsq, pro, n, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, G

      double precision   hood, Vinv

c     double precision   x(n), mu(G), sigsq(G), pro(G[+1]), z(n,G[+1])
      double precision   x(*), mu(*), sigsq(*), pro(  *  ), z(n,  *  )

      integer                 i, k, nz

      double precision        temp, const, tmin, tmax, sum
      double precision        muk, sigsqk, prok, sigmin

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call sgnrng( G, sigsq, 1, sigmin, temp)  

      if (sigmin .le. zero) then
        hood = FLMAX
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        muk    = mu(k)
        sigsqk = sigsq(k)
        const  = pi2log + log(sigsqk)
        do i = 1, n
          temp   = x(i) - muk
c         z(i,k) = prok*exp(-(const+(temp*temp)/sigsqk)/two)
          if (sigsqk .lt. one .and. 
     *        abs(temp) .ge. sqrt(sigsqk)*RTMAX) then
            hood = FLMAX
            return
          end if 
          z(i,k) = -(const+(temp*temp)/sigsqk)/two
        end do
      end do
 
      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       temp = zero
c       do k = 1, nz
c         temp = temp + z(i,k)
c       end do
c       hood = hood + log(temp)
c       call dscal( nz, (one/temp), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if  
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine hc1v  ( x, n, ic, ng, ns, ALPHA, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer             n, ic(n), ng, ns, nd
      integer             n, ic(*), ng, ns, nd

c     double precision    x(n), ALPHA, d(ng*(ng-1)/2)
      double precision    x(*), ALPHA, d(*)

      integer             lg, ld, ll, lo, ls, i, j, k, m
      integer             ni, nj, nij, nopt, niop, njop  
      integer             ij, ici, icj, iopt, jopt, iold

      double precision    ALFLOG
      double precision    qi, qj, qij, ri, rj, rij, si, sj
      double precision    tracei, tracej, trcij, trop
      double precision    termi, termj, trmij, tmop
      double precision    temp, dij, dopt, siop, sjop

      double precision    zero, one
      parameter          (zero = 0.d0, one = 1.d0)

      double precision    sqrthf
      parameter          (sqrthf = .70710678118654757274d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

c------------------------------------------------------------------------------

c     call dblepr( 'x', -1, x, n) 
c     call intpr( 'n', -1, n, 1) 
c     call intpr( 'ic', -1, ic, n) 
c     call intpr( 'ng', -1, ng, 1) 
c     call intpr( 'ns', -1, ns, 1) 
c     call dblepr( 'alpha', -1, alpha, 1) 
c     call intpr( 'nd', -1, nd, 1) 

      iopt = 0
      jopt = 0
      niop = 0
      njop = 0
      nopt = 0
      siop = 0
      sjop = 0
      tmop = 0.d0
      trop = 0.d0

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd

      if (ng .eq. 1) return

      ALPHA  = max(ALPHA,EPSMAX)

      ALFLOG = log(ALPHA)

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
c         call dswap( p, x(k,1), n, x(j,1), n)
          temp  = x(k)
          x(k)  = x(j)
          x(j)  = temp
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
c           call dscal( p, sqrthf, x(i,1), n)
c           call dscal( p, sqrthf, x(j,1), n)
c           call dcopy( p, x(j,1), n, v, 1)
c           call daxpy( p, (-one), x(i,1), n, v, 1)
c           call daxpy( p, one, x(j,1), n, x(i,1), n)
c           x(j,1) = ddot( p, v, 1, v, 1)
            temp   = sqrthf*(x(j)  - x(i))
            x(i)   = sqrthf*(x(j)  + x(i))
            x(j)   = temp*temp
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
c           call dcopy( p, x(j,1), n, v, 1)
c           call dscal( p, si, v, 1)
c           call daxpy( p, (-sj), x(i,1), n, v, 1)
c           x(k,1) =    x(k,1) +    ddot(p, v, 1, v, 1)
c           call dscal( p, si, x(i,1), n)
c           call daxpy( p, sj, x(j,1), n, x(i,1), n)
            temp = si*x(j) - sj*x(i)
            x(k) = x(k) + temp*temp
            x(i) = si*x(i) + sj*x(j)
          end if 
        else 
          ic(j) = 1
        end if
      end do
       
c store terms also so as not to recompute them

      do k = 1, ng
        i = ic(k)
        if (i .ne. 1) then
          ni        = ic(i)
          ri        = dble(ni)
          d(nd-k+1) = ri*log((x(i)+ALPHA)/ri)
        end if
      end do

c     call intpr( 'ic', -1, ic, n)

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
          tracej = x(nj)
          nj     = ic(nj)
          rj     = dble(nj)
          termj  = d(nd-j+1)
        end if
        do i = 1, (j-1)
          ni = ic(i)
          if (ni .eq. 1) then
            tracei = zero
            termi  = ALFLOG
            ri     = one
          else 
            tracei = x(ni)
            ni     = ic(ni)
            ri     = dble(ni)
            termi  = d(nd-i+1)
          end if               
          nij = ni + nj
          rij = dble(nij)
          qij = one/rij
          qi  = ri*qij
          qj  = rj*qij
          si  = sqrt(qi)
          sj  = sqrt(qj)
c         call dcopy(p, x(i,1), n, v, 1)
c         call dscal( p, sj, v, 1)
c         call daxpy( p, (-si), x(j,1), n, v, 1)
          temp  = sj*x(i) - si*x(j)
c         trcij = (tracei + tracej) + ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + temp*temp
          trmij = rij*log((trcij+ALPHA)/rij)
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
          end if
        end do
      end do

c     call dblepr( 'dij', -1, d, (ng*(ng-1))/2)

      if (ns .eq. 1) then
        if (iopt .lt. jopt) then
          x(1)  = dble(iopt)
          ic(1) = jopt
        else
          x(1)  = dble(jopt)
          ic(1) = iopt
        end if
        d(1)    = dopt
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

c     call dscal( p, siop, x(iopt,1), n)
c     call daxpy( p, sjop, x(jopt,1), n, x(iopt,1), n)
      x(iopt) = siop*x(iopt)+sjop*x(jopt)

      if (jopt .ne. lg) then
        call wardsw( jopt, lg, d)
c       call dcopy( p, x(lg,1), n, x(jopt,1), n)
        x(jopt)  = x(lg)
        m        = ic(jopt)
        ic(jopt) = ic(lg)
        ic(lg)   = m
      end if

      ic(iopt) =  lg
c     ic(lg)   =  nopt
c     x(lg,1)  =  trop
      x(lg)    =  trop
c     x(lg,2)  =  tmop

      d(lo)  = dopt
      lo     = lo - 1
      ic(lg) =  lo
      d(lo)  = tmop
      lo     = lo - 1
      d(lo)  = dble(nopt)
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
c           tracej = x(nj,1)
            tracej = x(nj)
            k      = ic(nj)
            termj  = d(k)
            nj     = int(d(k-1))
            rj     = dble(nj)
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
c         call dcopy( p, x(iold,1), n, v, 1)
c         call dscal( p, sj, v, 1)
c         call daxpy( p, (-si), x(j,1), n, v, 1)
          temp  = sj*x(iold)-si*x(j)
c         trcij = (tracei + tracej) + ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + temp*temp
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
c           tracej = x(nj,1)
            tracej = x(nj)
            k      = ic(nj)
            termj  = d(k)
            nj     = int(d(k-1))
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
c         call dcopy( p, x(iold,1), n, v, 1)
c         call dscal( p, sj, v, 1)
c         call daxpy( p, (-si), x(j,1), n, v, 1)
          temp  = sj*x(iold) - si*x(j)
c         trcij = (tracei + tracej) + ddot(p,v,1,v,1)
          trcij = (tracei + tracej) + temp*temp
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
          tracei = x(i)
          ici    = ic(i)
          termi  = d(ici)
          niop   = int(d(ici-1))
          ri     = dble(niop)
        else
          tracei = zero
          termi  = ALFLOG
          niop   = 1
          ri     = one
        end if

        if (j .ne. 1) then
c         tracej = x(j,1)
          tracej = x(j)
          icj    = ic(j)
          termj  = d(icj)
          njop   = int(d(icj-1))
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
c       call dcopy( p, x(iopt,1), n, v, 1)
c       call dscal( p, sjop, v, 1)
c       call daxpy( p, (-siop), x(jopt,1), n, v, 1)
        temp  = sjop*x(iopt)-siop*x(jopt)
c       trop  = (tracei + tracej) + ddot(p,v,1,v,1)
        trop  = (tracei + tracej) + temp*temp
        tmop  = rij*log((trop+ALPHA)/rij)
      end if
     
      ls = ls + 1

      if (ls .eq. ns) goto 900

      goto 100

 900  continue

      d(lo)  = dopt
      lo     = lo - 1
      d(lo)  = tmop
      lo     = lo - 1
      d(lo)  = dble(nopt)
      lo     = lo - 1
      d(lo)  = dble(iopt)
      lo     = lo - 1
      d(lo)  = dble(jopt)

      do i = 1, ng
        ic(i) = i
      end do

      lo          = nd - 3
      ld          = nd - 1
      si          = d(lo)
      lo          = lo - 1
      sj          = d(lo)
      lo          = lo - 1
      ic(int(sj)) = ng

      if (si .lt. sj) then
        x(1)  = si 
        d(ld) = sj
      else
        x(1)  = sj
        d(ld) = si
      end if
      ld = ld - 1

      lg = ng + 1
      do k = 2, ns
        d(ld)  = d(lo)
        ld     = ld - 1
        lo     = lo - 3
        i      = int(d(lo))
        ici    = ic(i)
        lo     = lo - 1
        j      = int(d(lo))
        lo     = lo - 1
        icj    = ic(j)
        if (ici .gt. icj) ic(i) = icj
        ic(j)  = ic(lg-k)
        if (ici .lt. icj) then
          x(k)  = dble(ici)
          d(ld) = dble(icj)
        else
          x(k)  = dble(icj)
          d(ld) = dble(ici)
        end if
        ld = ld - 1
      end do

      ld = nd
      lo = nd - 1
      do k = 1, ns
        ic(k) = int(d(lo))
        lo    = lo - 1
        ld    = ld - 1
        d(ld) = d(lo)
        lo    = lo - 1
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

      subroutine me1v ( EQPRO, x, n, G, Vinv, z, maxi, tol, eps,
     *                  mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical              EQPRO

      integer              n, G, maxi

      double precision     Vinv, eps, tol

c     double precision     x(n), z(n,G[+1]), mu(G), sigsq(G), pro(G[+1])
      double precision     x(*), z(n,  *  ), mu(*), sigsq(*), pro(  *  )

      integer                 nz, iter, k, i

      double precision        hold, hood, err, sum, smu, zsum
      double precision        const, temp, term, sigmin, sigsqk
      double precision        prok, tmin, tmax, ViLog, rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      rteps  = sqrt(eps)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      iter   = iter + 1

      zsum   = one

      do k = 1, G
        sum = zero
        smu = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
        if (.not. EQPRO) pro(k)   = sum / dble(n)
        zsum = min(sum,zsum)
        if (sum .gt. rteps) then
          smu    = smu / sum
          mu(k)  = smu
          sigsqk = zero
          do i = 1, n
            temp   = x(i) - smu
            temp   = temp*temp
            sigsqk = sigsqk + z(i,k)*temp
            z(i,k) = temp
          end do
          sigsq(k) = sigsqk / sum
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
        
      end if

      sigmin = FLMAX
      do k = 1, G
        sigmin = min(sigmin,sigsq(k))
      end do

      if (sigmin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G
        sigsqk = sigsq(k)
        const  = pi2log + log(sigsqk)
        do i = 1, n
c         z(i,k) = temp*exp(-(const+(z(i,k)/sigsqk))/two)           
          z(i,k) = -(const+(z(i,k)/sigsqk))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine me1vp ( EQPRO, x, n, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c       http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical              EQPRO

      integer              n, G, maxi

      double precision     pshrnk, pmu, pscale, pdof

      double precision     Vinv, eps, tol

c     double precision     x(n), z(n,G[+1]), mu(G), sigsq(G), pro(G[+1])
      double precision     x(*), z(n,  *  ), mu(*), sigsq(*), pro(  *  )

      integer                 nz, iter, k, i

      double precision        hold, hood, err, pmupmu
      double precision        sumz, sum, smu, zsum, rteps
      double precision        const, temp, term, sigmin, sigsqk
      double precision        prok, tmin, tmax, ViLog
      double precision        cmu, cgam, rmu, rgam

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        three
      parameter              (three = 3.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      double precision        dlngam
      external                dlngam

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      rteps  = sqrt(eps)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX
 
      pmupmu = pmu*pmu

      iter   = 0

100   continue

      iter   = iter + 1

      zsum   = one

      do k = 1, G
        sumz = zero
        smu  = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          smu  = smu + temp*x(i)
        end do
        if (.not. EQPRO) pro(k)   = sumz / dble(n)
        zsum = min(zsum,sumz) 
        if (sumz .gt. rteps) then
          smu    = smu/sumz
          sum    = zero
          do i = 1, n
            term   = abs(x(i) - smu)
            if (term .ge. eps .or. sqrt(z(i,k))*term .gt. RTMIN) 
     *                 sum = sum + z(i,k)*(term*term)
          end do
          term     = (pshrnk*sumz)/(pshrnk+sumz)
          temp     = (pmupmu + smu*smu) - two*pmu*smu
          sigsq(k) = (pscale + sum + term*temp)/(pdof+sumz+three)
          term     = sumz/(pshrnk+sumz)
          temp     = pshrnk/(pshrnk+sumz)
          mu(k)    = term*smu + temp*pmu
        end if
      end do

      if (zsum .le. rteps) then 
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
        
      end if

      sigmin = FLMAX
      do k = 1, G
        sigmin = min(sigmin,sigsq(k))
      end do

      if (sigmin .le. eps) then 
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G
        sigsqk = sigsq(k)
        const  = pi2log + log(sigsqk)
        do i = 1, n
          term   = abs(x(i) - mu(k))
          if (term .gt. RTMIN) then
            z(i,k) = -(const+((term*term)/sigsqk))/two
          else 
            z(i,k) = -const/two
          end if   
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol .and. iter .lt. maxi) goto 100

      tol   = err
      eps   = hood
      maxi  = iter

      cmu   = dble(G)*(pi2log-log(pshrnk))/two

      const = pdof/two
      cgam  = dble(G)*(const*log(pscale/two) - dlngam(const))

      rmu   = zero
      rgam  = zero
      do k = 1, G
        temp = pmu - mu(k)
        temp = temp*temp
        term = log(sigsq(k))
        rmu  = rmu + (term + (pshrnk/sigsq(k))*temp)
        rgam = rgam + ((pdof+3.d0)*term + pscale/sigsq(k))
      end do
      rmu   = -rmu /two    
      rgam  = -rgam/two
  
      pdof  = (cmu+cgam) + (rmu+rgam)

      return
      end

      subroutine ms1v ( x, z, n, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, G

c     double precision   x(n), z(n,G), mu(G), sigsq(G), pro(G)
      double precision   x(*), z(n,*), mu(*), sigsq(*), pro(*)

      integer                 i, k
     
      double precision        sum, smu, temp, sigsqk

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      do k = 1, G
        sum = zero
        smu = zero      
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
        pro(k) = sum / dble(n)
        if (sum .gt. one .or. smu .le. sum*FLMAX) then
          smu    = smu / sum
          mu(k)  = smu 
          sigsqk = zero
          do i = 1, n
            temp = abs(x(i) - smu)
            sigsqk = sigsqk + z(i,k)*(temp*temp)
          end do
          sigsq(k) = sigsqk / sum
        else 
          mu(k)    = FLMAX  
          sigsq(k) = FLMAX
        end if 
      end do

      return
      end

      subroutine ms1vp ( x, z, n, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c       http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer              n, G

      double precision     pshrnk, pmu, pscale, pdof

c     double precision     x(n), z(n,G), mu(G), sigsq(G), pro(G)
      double precision     x(*), z(n,*), mu(*), sigsq(*), pro(*)

      integer              k, i

      double precision     pmupmu
      double precision     sumz, sum, smu
      double precision     temp, term

      double precision     zero, one, two
      parameter           (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision     FLMAX
      parameter           (FLMAX = 1.7976931348623157d308)

      double precision     RTMIN
      parameter           (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero 

      pmupmu = pmu*pmu

      do k = 1, G
        sumz = zero
        smu  = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          smu  = smu + temp*x(i)
        end do
        pro(k) = sumz / dble(n)
        if (sumz .gt. one .or. smu .le. sumz*FLMAX) then
          smu   = smu/sumz
          term  = sumz/(pshrnk+sumz)
          temp  = pshrnk/(pshrnk+sumz)
          mu(k) = term*smu + temp*pmu
          sum   = zero
          do i = 1, n
            term   = abs(x(i) - smu)
            sum = sum + z(i,k)*(term*term)
          end do
          term = (pshrnk*sumz)/(pshrnk+sumz)
          temp = (pmupmu + smu*smu) - two*pmu*smu
          if (pshrnk .gt. zero) then
            sigsq(k) = (pscale + sum + term*temp)/(pdof+sumz+3.d0)
          else
            sigsq(k) = (pscale + sum + term*temp)/(pdof+sumz+two)
          end if
        else 
          mu(k)    = FLMAX
          sigsq(k) = FLMAX
        end if
      end do

      return
      end

      subroutine esvei ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), scale(G), shape(p), pro(G[+1])
      double precision   mu(p,*), scale(*), shape(*), pro(  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        smin, smax, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            if (shape(j) .lt. one .and.
     *          abs(temp) .ge. shape(j)*FLMAX) then
              hood = FLMAX
              return
            end if
            temp = temp/shape(j)
            if (abs(temp) .ge. RTMAX) then 
              hood = FLMAX
              return
            end if
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum/scalek)/two)
          if (scalek .lt. one .and.
     *        sum .ge. scalek*FLMAX) then
            hood = FLMAX
            return
          end if
          z(i,k) = -(const+sum/scalek)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and.
     *      one .le. sum*FLMAX) then
          hood = FLMAX
          return
        end if  
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine mevei ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi(2)

      double precision    Vinv, eps, tol(2)

c     double precision    x(n,p), z(n,G[+1]), scl(G), shp(p), w(p,G)
      double precision    x(n,*), z(n,  *  ), scl(*), shp(*), w(p,*)

c     double precision    mu(p,G), scale(G), shape(p), pro(G[+1])
      double precision    mu(p,*), scale(*), shape(*), pro(  *  )

      integer             nz, i, j, k
      integer             iter, maxi1, maxi2, inner, inmax

      double precision    tol1, tol2, sum, temp, term, tmin, tmax
      double precision    prok, scalek, smin, smax, const, zsum
      double precision    hold, hood, err, errin, dnp, ViLog, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG 
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      maxi1 = maxi(1)
      maxi2 = max(maxi(2),0)

      if (maxi1 .le. 0) return

      dnp   = dble(n*p)
    
      inmax = 0

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps   = max(eps,zero)
      tol1  = max(tol(1),zero)
      tol2  = max(tol(2),zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX
      errin = FLMAX

c start with shape and scale equal to 1

      call dcopy(p, one, 0, shape, 1) 

      call dcopy(G, one, 0, scale, 1) 

      iter  = 0

100   continue

      inner = 0

      zsum = one

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sum
        zsum   = min(zsum,sum)
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              sum  = sum + z(i,k)*(temp*temp)
            end do
            w(j,k) = sum
          end do
        end if
      end do

      call dscal( G, dble(p), pro, 1)

      if (zsum .le. rteps) then
        eps     = -FLMAX
        tol(1)  =  zsum
        tol(2)  =  errin
        maxi(1) =  iter
        maxi(2) =  max(inner,inmax)
        return
      end if

      if (maxi2 .le. 0) goto 120

110   continue

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      inner = inner + 1

c scale estimate

      call dcopy( G, scale, 1, scl, 1)

      do k = 1, G
        sum = zero
        do j = 1, p
          sum = sum + w(j,k)/shape(j)
        end do
        scale(k) = sum/pro(k)
      end do

      call sgnrng(G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

c shape estimate

      call dcopy( p, shape, 1, shp, 1)

      do j = 1, p
        sum = zero
        do k = 1, G
          sum = sum + w(j,k)/scale(k)
        end do
        shape(j) = sum
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      if (temp .gt. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if

      if (temp .le. eps) then
        eps     = temp
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      call dscal( p, one/temp, shape, 1)

      errin = zero

      do k = 1, G
        errin = max(errin, abs(scl(k)-scale(k))/(one + scale(k)))
      end do

      do j = 1, p
        errin = max(errin, abs(shp(j)-shape(j))/(one + shape(j)))
      end do

      if (errin .gt. tol2 .and. inner .le. maxi2) goto 110

120   continue

      iter = iter + 1

      inmax = max(inner, inmax)

      if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)

      term = zero
      if (Vinv .gt. zero) then
     
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j)
          end do
c         z(i,k) = prok*exp(-(const+sum/scalek)/two)
          z(i,k) = -(const+sum/scalek)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol1 .and. iter .lt. maxi1) goto 100

      tol(1)  = err
      tol(2)  = errin
      eps     = hood
      maxi(1) = iter
      maxi(2) = inmax

      return
      end

      subroutine meveip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi(2)

      double precision    Vinv, eps, tol(2)

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G[+1]), scl(G), shp(p), w(p,G)
      double precision    x(n,*), z(n,  *  ), scl(*), shp(*), w(p,*)

c     double precision    mu(p,G), scale(G), shape(p), pro(G[+1])
      double precision    mu(p,*), scale(*), shape(*), pro(  *  )

      integer             nz, i, j, k
      integer             iter, maxi1, maxi2, inner, inmax

      double precision    tol1, tol2, sum, temp, term, tmin, tmax
      double precision    prok, scalek, smin, smax, const, sumz
      double precision    hold, hood, err, errin, dnp, ViLog, zsum
      double precision    rteps
    
      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero
 
      maxi1 = maxi(1)
      maxi2 = max(maxi(2),0)

      if (maxi1 .le. 0) return

      dnp   = dble(n*p)
    
      inmax = 0

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps   = max(eps,zero)
      tol1  = max(tol(1),zero)
      tol2  = max(tol(2),zero)
 
      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX
      errin = FLMAX

c start with shape and scale equal to 1

      call dcopy(p, one, 0, shape, 1) 

      call dcopy(G, one, 0, scale, 1) 

      iter  = 0

100   continue

      inner = 0

      zsum  = one

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sumz   = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz
        zsum   = min(zsum,sumz)
        if (sumz .gt. rteps) then
          term   = pshrnk + sumz
          const  = (pshrnk*sumz)/term
          call dscal( p, (one/sumz), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              sum  = sum + z(i,k)*(temp*temp)
            end do
            temp   = pmu(j) - mu(j,k)
            w(j,k) = pscale + sum + const*(temp*temp)
          end do
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1) 
        end if
      end do

      call dscal( G, dble(p), pro, 1)

      if (zsum .le. rteps) then
        eps     = -FLMAX
        tol(1)  =  zsum
        tol(2)  =  errin
        maxi(1) =  iter
        maxi(2) =  max(inner,inmax)
        return
      end if

      if (maxi2 .le. 0) goto 120

110   continue

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      inner = inner + 1

c scale estimate

      call dcopy( G, scale, 1, scl, 1)

      temp = pdof + two
      if (pshrnk .gt. zero) temp = temp + one
      do k = 1, G
        sum = zero
        do j = 1, p
          sum = sum + w(j,k)/shape(j)
        end do
        scale(k) = sum/(pro(k)+temp)
      end do

      call sgnrng(G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

c shape estimate

      call dcopy( p, shape, 1, shp, 1)

      do j = 1, p
        sum = zero
        do k = 1, G
          sum = sum + w(j,k)/scale(k)
        end do
        shape(j) = sum
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      if (temp .gt. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if

      if (temp .le. eps) then
        eps     = temp
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = max(inner,inmax)
        return
      end if

      call dscal( p, one/temp, shape, 1)

      errin = zero

      do k = 1, G
        errin = max(errin, abs(scl(k)-scale(k))/(one + scale(k)))
      end do

      do j = 1, p
        errin = max(errin, abs(shp(j)-shape(j))/(one + shape(j)))
      end do

      if (errin .gt. tol2 .and. inner .le. maxi2) goto 110

120   continue

      iter = iter + 1

      inmax = max(inner, inmax)

      if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)

      term = zero
      if (Vinv .gt. zero) then
     
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j)
          end do
c         z(i,k) = prok*exp(-(const+sum/scalek)/two)
          z(i,k) = -(const+sum/scalek)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol1 .and. iter .lt. maxi1) goto 100

      tol(1)  = err
      tol(2)  = errin
      eps     = hood
      maxi(1) = iter
      maxi(2) = inmax

      return
      end

      subroutine msvei ( x, z, n, p, G, maxi, tol, 
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G, maxi

      double precision   tol

c     double precision   x(n,p), z(n,G), scl(G), shp(p), w(p,G)
      double precision   x(n,*), z(n,*), scl(*), shp(*), w(p,*)

c     double precision   mu(p,G), scale(G), shape(p), pro(G)
      double precision   mu(p,*), scale(*), shape(*), pro(*)

      integer                 i, j, k, inner

      double precision        sum, temp, smin, smax, err

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      tol = max(tol,zero)
      err = FLMAX

c start with the equal volume and shape estimate

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sum
        if (sum .gt. one .or. one .lt. sum*FLMAX) then
          err = min(err,sum)
          call dscal( p, (one/sum), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              temp = temp*temp
              temp = z(i,k)*temp
              sum  = sum + temp
            end do
            w(j,k) = sum
          end do
        else
          err = -FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (err .lt. zero) then
        call dscal( G, one/dble(n), pro, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        tol  = FLMAX
        maxi = 0
        return
      end if

      call dcopy( p, one, 0, shape, 1)

      call dcopy( G, one, 0, scale, 1)

      call dscal( G, dble(p), pro, 1)

      inner = 0
      err  = FLMAX

100   continue

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) goto 200

      inner = inner + 1

c scale estimate

      call dcopy( G, scale, 1, scl, 1)

      do k = 1, G
        sum = zero       
        do j = 1, p
          if (shape(j) .gt. one .or. 
     *        w(j,k) .lt. shape(j)*FLMAX) then
            sum = sum + w(j,k)/shape(j)
          else
            scale(k) = FLMAX
            goto 110
          end if
        end do
        scale(k) = sum/pro(k)
110     continue
      end do

      call sgnrng(G, scale, 1, smin, smax)

      if (smin .le. zero .or. smax .eq. FLMAX) goto 200

c shape estimate

      call dcopy( p, shape, 1, shp, 1)

      do j = 1, p
        sum = zero
        do k = 1, G
          if (scale(k) .gt. one .or. w(j,k) .lt. scale(k)*FLMAX) then
            sum = sum + w(j,k)/scale(k)
          else
            shape(j) = FLMAX
            goto 120
          end if
        end do
        shape(j) = sum
120     continue
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero .or. smax .eq. FLMAX) goto 200

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do      
      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        call dcopy( G, FLMAX, 0, scale, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200 
      end if 
  
      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if 

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200
      end if

      call dscal( p, one/temp, shape, 1)

      err = zero
      
      do k = 1, G
        err = max(err, abs(scl(k) - scale(k))/(one + scale(k)))        
      end do
      
      do j = 1, p
        err = max(err, abs(shp(j) - shape(j))/(one + shape(j)))       
      end do
 
      if (err .gt. tol .and. inner .le. maxi) goto 100

200   continue

      call dscal( G, one/dble(n*p), pro, 1)
    
      tol  = err
      maxi = inner

      return
      end

      subroutine msveip( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   maxi, tol,
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, G, maxi

      double precision    tol

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G), scl(G), shp(p), w(p,G)
      double precision    x(n,*), z(n,*), scl(*), shp(*), w(p,*)

c     double precision    mu(p,G), scale(G), shape(p), pro(G)
      double precision    mu(p,*), scale(*), shape(*), pro(*)

      integer             i, j, k, inner

      double precision    sum, temp, term, err
      double precision    smin, smax, const, sumz
    
      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------
    
      tol   = max(tol,zero)

      err   = FLMAX

c start with shape and scale equal to 1

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz
        if (sumz .gt. one .or. one .lt. sumz*FLMAX) then
          err   = min(err,sumz)
          term  = pshrnk + sumz
          const = (pshrnk*sumz)/term
          call dscal( p, (one/sumz), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = x(i,j) - mu(j,k)
              sum  = sum + z(i,k)*(temp*temp)
            end do
            temp   = pmu(j) - mu(j,k)
            w(j,k) = pscale + sum + const*(temp*temp)
          end do
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1) 
        else
          err = -FLMAX
          call dcopy( p, FLMAX, 0,  mu(1,k), 1)
        end if
      end do

      if (err .lt. zero) then
        call dscal( G, one/dble(n), pro, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        tol  = FLMAX
        maxi = 0
        return
      end if

      call dcopy(p, one, 0, shape, 1) 

      call dcopy(G, one, 0, scale, 1) 

      call dscal( G, dble(p), pro, 1)

      if (maxi .le. 0) return

      inner = 0
      err   = FLMAX

100   continue

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) goto 200

      inner = inner + 1

c scale estimate

      call dcopy( G, scale, 1, scl, 1)

      do k = 1, G
        sum = zero
        do j = 1, p
          if (shape(j) .ge. one .or. 
     *        w(j,k) .le. shape(j)*FLMAX) then
             sum = sum + w(j,k)/shape(j)
          else 
            scale(k) = FLMAX
            goto 110
          end if
        end do
        temp  = pdof + pro(k) + two
        if (pshrnk .gt. zero) temp = temp + one
        scale(k) = sum/temp
110     continue
      end do

      call sgnrng(G, scale, 1, smin, smax)

      if (smin .le. zero .or. smax .ge. FLMAX) then
        call dcopy( G, FLMAX, 0, scale, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200
      end if

c shape estimate

      call dcopy( p, shape, 1, shp, 1)

      do j = 1, p
        sum = zero
        do k = 1, G
          if (scale(k) .gt. w(j,k) .or. 
     *        w(j,k) .lt. scale(k)*FLMAX) then
            sum = sum + w(j,k)/scale(k)
          else
            shape(j) = FLMAX
            goto 120
          end if
        end do
        shape(j) = sum
120     continue
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero .or. smax .ge. FLMAX) then
        call dcopy( G, FLMAX, 0, scale, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = sum/dble(p)

      if (temp .ge. BIGLOG) then
        call dcopy( G, FLMAX, 0, scale, 1)
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200
      end if

      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else
        temp = zero
      end if

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200         
      end if

      call dscal( p, one/temp, shape, 1)

      err = zero

      do k = 1, G
        err = max(err, abs(scl(k)-scale(k))/(one + scale(k)))
      end do

      do j = 1, p
        err = max(err, abs(shp(j)-shape(j))/(one + shape(j)))
      end do

      if (err .gt. tol .and. inner .le. maxi) goto 100

200   continue
 
      call dscal( G, one/dble(n*p), pro, 1)
    
      tol  = err
      maxi = inner

      return
      end

      subroutine esvev ( x, mu, scale, shape, O, pro, n, p, G, 
     *                   Vinv, v, w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p, G
      integer            n, p, G

      double precision   Vinv, hood

c     double precision   x(n,p), z(n,G[+1]), mu(p,G), pro(G[+1])
      double precision   x(n,*), z(n,  *  ), mu(p,*), pro(  *  )

c     double precision   v(p), w(p)
      double precision   v(*), w(*)

c     double precision   scale(G), shape(p), O(p,p,G)
      double precision   scale(*), shape(*), O(p,p,*)

      integer                 i, j, k, nz

      double precision        const, temp, tmin, tmax
      double precision        smin, smax, scalek, prok, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G

        scalek = scale(k)
        
        const = dble(p)*(pi2log+log(scalek))

c       prok  = pro(k)

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, 
     *                 w, 1, zero, v, 1)
          do j = 1, p
            if (shape(j) .lt. one .and. 
     *          abs(v(j)) .ge. shape(j)*FLMAX) then
              hood = FLMAX
              return
            end if
            v(j) = v(j)/shape(j)
          end do
          temp   = ddot( p, v, 1, v, 1)
          if (scalek .lt. one .and. temp .ge. scalek*FLMAX) then
            hood = FLMAX
            return
          end if
          temp   = temp/scalek
c         z(i,k) = prok*exp(-(const+temp)/two)
          z(i,k) = -(const+temp)/two
        end do

      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine mevev ( EQPRO, x, n, p, G, Vinv, z, 
     *                   maxi, tol, eps, lwork,
     *                   mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi(2), lwork

      double precision   Vinv, eps, tol(2)

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

      double precision   mu(p,*), pro(  *  )

      double precision   scale(*), shape(*), O(p,p,*)

      integer                 maxi1, maxi2, p1, inmax, iter
      integer                 nz, i, j, k, l, j1, info, inner

      double precision        tol1, tol2, dnp, term, rteps, ViLog
      double precision        errin, smin, smax, sumz, tmin, tmax
      double precision        cs, sn, dummy, hold, hood, err, zsum
      double precision        const, temp, sum, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

      if (maxi1 .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol1   = max(tol(1),zero)
      tol2   = max(tol(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

c     FLMAX  = d1mach(2)

      hold   = FLMAX/two
      hood   = FLMAX

      err    = FLMAX
      errin  = FLMAX

      inmax  = 0

      iter   = 0

100   continue

      sumz = zero
      zsum = one

      l    = 0
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        zsum = min(zsum,sum) 
        pro(k) = sum
        if (sum .ge. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else 
            do j = 1, p
              temp     = z(j,k)
              z(j,k)   = temp*temp
            end do
          end if
        end  if
      end do

      iter = iter + 1

      if (l .ne. 0 .or. zsum .lt. rteps) then

        if (Vinv .ge. zero) then
          term = zero
          do i = 1, n
            term = term + z(i,nz)
          end do
          temp    = term / dble(n)
          pro(nz) = temp

          call dcopy( n, ViLog, 0, z(1,nz), 1)

          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else
          if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
        end if
        lwork   = l
c       w(1)    = FLMAX
        tol(1)  = err
        tol(2)  = errin
        if (l .ne. 0) then
          eps =  FLMAX
        else
          eps = -FLMAX
        end if
        maxi(1) = -1
        maxi(2) = -1
        return
      end if

      if (iter .eq. 1) then
        call dcopy( p, zero, 0, shape, 1)
        do j = 1, p
          sum = zero
          do k = 1, G
            sum = sum + z(j,k)
          end do
          shape(j) = sum
        end do

        call sgnrng( p, shape, 1, smin, smax)

        if (smin .le. zero) then
          if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
          if (Vinv .ge. zero) then
            term = zero
            do i = 1, n
              term = term + z(i,nz)
            end do
            temp    = term / dble(n)
            pro(nz) = temp

            call dcopy( n, ViLog, 0, z(1,nz), 1)

            if (EQPRO) then
              temp = (one - pro(nz))/dble(G)
              call dcopy( G, temp, 0, pro, 1)
            end if
          else if (EQPRO) then
            call dcopy( G, one/dble(G), 0, pro, 1)
          end if
          lwork   = 0
c         w(1)    = smin
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp  = sum/dble(p)
 
        if (temp .gt. BIGLOG) then
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if 

        if (temp .gt. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if 

        if (Vinv .le. zero) then
          call dcopy (G, temp/dble(n), 0, scale, 1)
        else
          call dcopy (G, temp/sumz, 0, scale, 1)
        end if

        if (temp .le. eps) then
          if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
          if (Vinv .gt. zero) then
            term = zero
            do i = 1, n
              term = term + z(i,nz)
            end do
            temp    = term / dble(n)
            pro(nz) = temp

            call dcopy( n, ViLog, 0, z(1,nz), 1)

            if (EQPRO) then
              temp = (one - pro(nz))/dble(G)
              call dcopy( G, temp, 0, pro, 1)
            end if
          else 
           if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
          end if
          lwork   = 0
c         w(1)    = temp
c         w(2)    = zero
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if
  
        call dscal( p, one/temp, shape, 1)

      end if

c inner iteration to estimate scale and shape
c pro now contains n*pro

      inner = 0
      errin = zero

      if (maxi2 .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + z(j,k)/w(j)
          end do
          temp     = sum/(pro(k)*dble(p))
          scale(k) = temp
          if (temp .le. eps) then
            lwork   = 0
c           w(1)    = temp
            tol(1)  = err
            tol(2)  = errin
            eps     = FLMAX
            maxi(1) = iter
            maxi(2) = max(inner,inmax)
            return
          end if
          do j = 1, p
            shape(j) = shape(j) + z(j,k)/temp
          end do
        end do

        inner = inner + 1

        call sgnrng( p, shape, 1, smin, smax)

        if (smin .le. zero) then
           if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
           if (Vinv .gt. zero) then
             term = zero
             do i = 1, n
               term = term + z(i,nz)
             end do
             temp    = term / dble(n)
             pro(nz) = temp

             call dcopy( n, ViLog, 0, z(1,nz), 1)

             if (EQPRO) then
               temp = (one - pro(nz))/dble(G)
               call dcopy( G, temp, 0, pro, 1)
             end if
          else 
            if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
          end if
          lwork   = 0
c         w(1)    = smin
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = iter
          maxi(2) = max(inner,inmax)
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp  = sum/dble(p)
 
        if (temp .gt. BIGLOG) then
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if 

        if (temp .gt. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if 

        if (temp .le. eps) then
           if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
           if (Vinv .gt. zero) then
             term = zero
             do i = 1, n
               term = term + z(i,nz)
             end do
             temp    = term / dble(n)
             pro(nz) = temp

             call dcopy( n, ViLog, 0, z(1,nz), 1)
             if (EQPRO) then
               temp = (one - pro(nz))/dble(G)
               call dcopy( G, temp, 0, pro, 1)
             end if
          else 
            if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
          end if
          lwork   = 0
c         w(1)    = temp
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
          maxi(1) = iter
          maxi(2) = max(inner,inmax)
        end if

        call dscal( p, one/temp, shape, 1)

        errin = zero
        do j = 1, p
          errin = max(abs(w(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(scale(k)-w(p+k))/(one+scale(k)), errin)
        end do

        if (errin .ge. tol2 .and. inner .lt. maxi2) goto 110

120   continue

      inmax = max(inner,inmax)

      smin = smin/temp
      smax = smax/temp

      if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)

      term = zero
      if (Vinv .gt. zero) then
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else 
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      call sgnrng( p, shape, 1, smin, smax)
      
      if (smin .le. eps) then
        lwork   = 0
c       w(1)    = -smin
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        lwork   = 0
c       w(1)    = -smin
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      do j = 1, p
        s(j) = sqrt(shape(j))
      end do

      call sgnrng( p, s, 1, smin, smax)
      
      if (smin .le. rteps) then
        lwork   = 0
c       w(1)    = -smin
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log + log(scalek))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / s(j)
          end do
          sum    = ddot(p,w,1,w,1)/scalek
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol1 .and. iter .lt. maxi1) goto 100

c     smin  = sqrt(smin)
c     smax  = sqrt(smax)

c     rcmin = FLMAX
c     do k = 1, G
c       temp = sqrt(scale(k))
c       rcmin = min(rcmin,(temp*smin)/(one+temp*smax))
c      end do

      lwork   = 0
     
c     w(1)    = rcmin

      tol(1)  = err
      tol(2)  = errin

      eps     = hood

      maxi(1) = iter
      maxi(2) = inmax

      return
      end

      subroutine mevevp( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, lwork,
     *                   mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi(2), lwork

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   Vinv, eps, tol(2)

c     double precision   x(n,p), z(n,G[+1]), w(lwork), s(p)
      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

c     double precision   mu(p,G), pro(G[+1])
      double precision   mu(p,*), pro(  *  )

c     double precision   scale(G), shape(p), O(p,p,G)
      double precision   scale(*), shape(*), O(p,p,*)

      integer                 maxi1, maxi2, p1, inmax, iter
      integer                 nz, i, j, k, l, j1, info, inner

      double precision        tol1, tol2, dnp, term, rteps, ViLog
      double precision        errin, smin, smax, sumz, tmin, tmax
      double precision        cs, sn, dummy, hold, hood, err, zsum 
      double precision        const, temp, sum, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

      if (maxi1 .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol1   = max(tol(1),zero)
      tol2   = max(tol(2),zero)

      p1     = p + 1

      dnp    = dble(n*p)

c     FLMAX  = d1mach(2)

      hold   = FLMAX/two
      hood   = FLMAX

      err    = FLMAX
      errin  = FLMAX

      inmax  = 0
      inner  = 0

      iter   = 0

100   continue

      zsum = one
      l    = 0
      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, O(1,j,k), 1)
        end do
        sumz = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz / dble(n)
        zsum   = min(zsum,sumz)
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          term  = sumz+pshrnk
          const = (sumz*pshrnk)/term
          call dscal( p, sqrt(const), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( O(j,j,k), w(j), cs, sn)
            call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( O(p,p,k), w(p), cs, sn)
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k),
     *                 dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else 
            do j = 1, p
              temp   = z(j,k)
              z(j,k) = temp*temp
            end do
          end if
        end if
      end do

      iter = iter + 1

      if (l .ne. 0 .or. zsum .le. rteps) then
        lwork   = l
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        if (l .ne. 0) then          
          eps =  FLMAX
        else 
          eps = -FLMAX
        end if
        goto 200
      end if

      if (iter .eq. 1) then

        call dcopy( p, zero, 0, shape, 1)
        do j = 1, p
          sum = zero
          do k = 1, G
            sum = sum + z(j,k)
          end do
          shape(j) = sum
        end do

        call sgnrng( p, shape, 1, smin, smax)

        if (smin .le. zero) then
          eps     = FLMAX
          goto 200
          return
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        temp  = sum/dble(p)

        if (temp .gt. BIGLOG) then
          eps  = FLMAX
          goto 200
          return
        end if 

        if (temp .ge. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if 

        do k = 1, G
          scale(k) = temp / (pro(k)*dble(n))
        end do

        if (temp .le. eps) then
          eps = FLMAX
          goto 200
          return
        end if
  
        call dscal( p, one/temp, shape, 1)

      end if

      inner = 0
      errin = zero

      if (maxi2 .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)
 
        call sgnrng( p+G, w, 1, smin, smax)

        if (smin .le. zero) then
          call dcopy( p, FLMAX, 0, shape, 1)
          call dcopy( G, FLMAX, 0, scale, 1)
          goto 200
        end if

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            if (w(j) .le. z(j,k) .and. z(j,k) .lt. w(j)*rteps) then
              call dcopy( p, FLMAX, 0, shape, 1) 
              call dcopy( G, FLMAX, 0, scale, 1)
              goto 200
            end if
            sum = sum + z(j,k)/w(j)
          end do
          temp     = sum/(pro(k)*dble(n*p))
          scale(k) = temp
          do j = 1, p
            if (temp .le. z(j,k) .and. z(j,k) .lt. temp*rteps) then
              call dcopy( p, FLMAX, 0, shape, 1) 
              call dcopy( G, FLMAX, 0, scale, 1)
              goto 200
            end if
            shape(j) = shape(j) + z(j,k)/temp
          end do
        end do

        inner = inner + 1

        call sgnrng( p, shape, 1, smin, smax)

        if (smin .le. zero) then
          eps  = FLMAX
          goto 200
          return
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        temp  = sum/dble(p)

        if (temp .gt. BIGLOG) then
          eps  = FLMAX
          goto 200
          return
        end if 

        if (temp .gt. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if 

        if (temp .le. eps) then
          eps = FLMAX
          goto 200
        end if

        call dscal( p, one/temp, shape, 1)

        errin = zero
        do j = 1, p
          errin = max(abs(w(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(scale(k)-w(p+k))/(one+scale(k)), errin)
        end do

        if (errin .ge. tol2 .and. inner .lt. maxi2) goto 110

120   continue

      inmax = max(inner,inmax)

      smin = smin/temp
      smax = smax/temp

      term = zero
      if (Vinv .gt. zero) then
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else 
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      call sgnrng( p, shape, 1, smin, smax)
      
      if (smin .le. eps) then
        eps = FLMAX
        goto 200 
        return
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        goto 200
        return
      end if

      do j = 1, p
        s(j) = sqrt(shape(j))
      end do

      call sgnrng( p, s, 1, smin, smax)
      
      if (smin .le. rteps) then
        eps = FLMAX
        goto 200
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log + log(scalek))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, w(p1), 1, zero, w, 1)
          do j = 1, p
            w(j) = w(j) / s(j)
          end do
          sum    = ddot(p,w,1,w,1)/scalek
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol1 .and. iter .lt. maxi1) goto 100

c     smin  = sqrt(smin)
c     smax  = sqrt(smax)

c     rcmin = FLMAX
c     do k = 1, G
c       temp = sqrt(scale(k))
c       rcmin = min(rcmin,(temp*smin)/(one+temp*smax))
c      end do
     
c     w(1)    = rcmin

      lwork   = 0

      eps     = hood

200   continue

      tol(1)  = err
      tol(2)  = errin

      maxi(1) = iter
      maxi(2) = inmax

      return
      end

      subroutine msvev ( x, z, n, p, G, w, lwork, maxi, tol, 
     *                   mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G, maxi, lwork

      double precision   tol

c     double precision   x(n,p), z(n,G), w(max(4*p,5*p-4,p+G))
      double precision   x(n,*), z(n,*), w(*)

c     double precision   scale(G), shape(p), O(p,p,G), mu(p,G), pro(G)
      double precision   scale(*), shape(*), O(p,p,*), mu(p,*), pro(*)

      integer                 p1, i, j, k, j1, inner, info

      double precision        dnp, err, dummy
      double precision        temp, sum, smin, smax, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX =  1.7976931348623157d308)

      double precision        BIGLOG
      parameter              (BIGLOG =  709.d0)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      tol   = max(tol,zero)

      p1    = p + 1

      err   = FLMAX
      
      inner = 0

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, zero, 0, O(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sum
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          if (lwork .gt. 0) then
            do i = 1, n
              call dcopy( p, x(i,1), n, w, 1)
              call daxpy( p, (-one), mu(1,k), 1, w, 1)
              call dscal( p, sqrt(z(i,k)), w, 1)
              j = 1
              do j1 = 2, p
                call drotg( O(j,j,k), w(j), cs, sn)
                call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
                j = j1
              end do
              call drotg( O(p,p,k), w(p), cs, sn)
            end do
            call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k),
     *                  dummy, 1, dummy, 1, w, lwork, info)
            if (info .ne. 0) then
              inner = info
            else
              do j = 1, p
                temp     = z(j,k)
                temp     = temp*temp
                shape(j) = shape(j) + temp
                z(j,k)   = temp
              end do
            end if
          end if
        else
          err = zero
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

c inner iteration estimates scale and shape
c pro now contains n*pro
   
      if (inner .ne. 0 .or. err .eq. zero) then
        lwork = inner
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      lwork = 0

      call sgnrng( p, shape, 1, smin, smax)
  
      if (smin .eq. zero) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do

      temp = sum/dble(p)

      if (temp .gt. BIGLOG) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else 
        temp = zero
      end if

      do k = 1, G
        scale(k) = temp / (pro(k)*dble(n))
      end do
     
      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        goto 200
      end if  

      call dscal( p, one/temp, shape, 1)

c iteration to estimate scale and shape
c pro now contains n*pro

      if (maxi .le. 0) goto 200
      
100   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)

        call absrng( p, w, 1, smin, smax)

        if (smin .le. one .and. one .ge. smin*FLMAX) then
          call dcopy( p, FLMAX, 0, shape, 1)
          call dcopy( G, FLMAX, 0, scale, 1)
          goto 200
        end if

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + z(j,k)/w(j)
          end do
          temp     = (sum/pro(k))/dble(p)
          scale(k) = temp
          if (temp .lt. one .and. one .ge. temp*FLMAX) then
            call dcopy( p, FLMAX, 0, shape, 1) 
            goto 200
          end if
          do j = 1, p
            shape(j) = shape(j) + z(j,k)/temp
          end do
        end do

        inner  = inner + 1

        call sgnrng( p, shape, 1, smin, smax)
 
        if (smin .le. zero) then
          call dcopy( p, FLMAX, 0, shape, 1) 
          goto 200
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        temp = sum/dble(p)
 
        if (temp .ge. BIGLOG) then
          call dcopy( p, FLMAX, 0, shape, 1)
          goto 200
        end if

        if (temp .ge. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if

        if (temp .lt. one .and. one .ge. temp*FLMAX) then
          call dcopy( p, FLMAX, 0, shape, 1)
          goto 200
        end if

        call dscal( p, one/temp, shape, 1)

        err = zero
        do j = 1, p
          err = max(abs(w(j)-shape(j))/(one+shape(j)), err)
        end do

        do k = 1, G
          err = max(abs(scale(k)-w(p+k))/(one+scale(k)), err)
        end do
        
        if (err .ge. tol .and. inner .lt. maxi) goto 100

200   continue

      call dscal( G, one/dble(n), pro, 1)

      tol  = err
      maxi = inner

      return
      end

      subroutine msvevp( x, z, n, p, G,
     *                   pshrnk, pmu, pscale, pdof,
     *                   w, lwork, maxi, tol, 
     *                   mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G, maxi, lwork

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   tol

c     double precision   x(n,p), z(n,G), w(lwork)
      double precision   x(n,*), z(n,*), w(  *  )

c     double precision   mu(p,G), pro(G)
      double precision   mu(p,*), pro(*)

c     double precision   scale(G), shape(p), O(p,p,G)
      double precision   scale(*), shape(*), O(p,p,*)

      integer            p1, i, j, k, l, j1, inner, info

      double precision   sum, term, temp, err, smin, smax
      double precision   sumz, cs, sn, dummy, const

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   FLMAX
      parameter         (FLMAX = 1.7976931348623157d308)

      double precision   SMALOG, BIGLOG
      parameter         (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (pshrnk .le. zero) pshrnk = zero

      tol   = max(tol,zero)

      p1    = p + 1

      err   = FLMAX

      inner = 0
      l     = 0 

      call dcopy( p, zero, 0, shape, 1)

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, O(1,j,k), 1)
        end do
        sumz = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz / dble(n)
        if (sumz .ge. one .or. one .lt. sumz*FLMAX) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( O(j,j,k), w(j), cs, sn)
              call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( O(p,p,k), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          term  = sumz+pshrnk
          const = (sumz*pshrnk)/term
          call dscal( p, sqrt(const), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( O(j,j,k), w(j), cs, sn)
            call drot( p-j, O(j,j1,k), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( O(p,p,k), w(p), cs, sn)
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, z(1,k),
     *                 dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
          else 
            do j = 1, p
              temp     = z(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              z(j,k)   = temp
            end do
          end if
        else
          err = zero 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (l .ne. 0 .or. err .eq. zero) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      call sgnrng( p, shape, 1, smin, smax)

      if (smin .le. zero) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do

      temp  = sum/dble(p)

      if (temp .gt. BIGLOG) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if

      if (temp .ge. SMALOG) then
        temp = exp(temp)
      else 
        temp = zero
      end if

      do k = 1, G
        scale(k) = temp / (pro(k)*dble(n))
      end do

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        call dcopy( G, FLMAX, 0, scale, 1)
        goto 200
      end if  
  
      call dscal( p, one/temp, shape, 1)

      if (maxi .le. 0) goto 200

100   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)

        call sgnrng( p+G, w, 1, smin, smax)

        if (smin .le. zero) then
          call dcopy( p, FLMAX, 0, shape, 1)
          call dcopy( G, FLMAX, 0, scale, 1)
          goto 200
        end if

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            if (w(j) .le. z(j,k) .and. z(j,k) .ge. w(j)*FLMAX) then
              call dcopy( p, FLMAX, 0, shape, 1) 
              call dcopy( G, FLMAX, 0, scale, 1)
              goto 200
            end if
            sum = sum + z(j,k)/w(j)
          end do
          temp     = sum/(pro(k)*dble(n*p))
          scale(k) = temp
          do j = 1, p
            if (temp .le. z(j,k) .and. z(j,k) .ge. temp*FLMAX) then
              call dcopy( p, FLMAX, 0, shape, 1) 
              call dcopy( G, FLMAX, 0, scale, 1)
              goto 200
            end if
            shape(j) = shape(j) + z(j,k)/temp
          end do
        end do

        inner = inner + 1

        call sgnrng( p, shape, 1, smin, smax)

        if (smin .le. zero) then
          call dcopy( p, FLMAX, 0, shape, 1) 
          goto 200
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        temp = sum/dble(p)

        if (temp .ge. BIGLOG) then
          call dcopy( p, FLMAX, 0, shape, 1)
          goto 200
        end if

        if (temp .ge. SMALOG) then
          temp = exp(temp)
        else
          temp = zero
        end if

        if (temp .lt. one .and. one .ge. temp*FLMAX) then
          call dcopy( p, FLMAX, 0, shape, 1)
          goto 200
        end if

        call dscal( p, one/temp, shape, 1)

        err = zero
        do j = 1, p
          err = max(abs(w(j)-shape(j))/(one+shape(j)), err)
        end do

        do k = 1, G
          err = max(abs(scale(k)-w(p+k))/(one+scale(k)), err)
        end do

        if (err .ge. tol .and. inner .lt. maxi) goto 100

200   continue

      lwork = l
      tol   = err
      maxi  = inner

      return
      end

      subroutine esvii ( x, mu, sigsq, pro, n, p, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), sigsq(G), pro(G[+1])
      double precision   mu(p,*), sigsq(*), pro(  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        prok, sigsqk, sigmin

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call sgnrng( G, sigsq, 1, sigmin, temp)    

      if (sigmin .le. zero) then
        hood  = FLMAX
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        sigsqk = sigsq(k)
        const  = dble(p)*(pi2log+log(sigsq(k)))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            if (abs(temp) .ge. RTMAX) then
              hood = FLMAX
              return  
            end if
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum/sigsqk)/two)
          if (sigsqk .lt. one .and. sum .ge. sigsqk*FLMAX) then
            hood = FLMAX
            return
          end if  
          z(i,k) = -(const+sum/sigsqk)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end 

      subroutine hcvii ( x, n, p, ic, ng, ns, ALPHA, v, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, ic(n), ng, ns, nd

c     double precision    x(n,p), v(p). d(*), ALPHA
      double precision    x(n,*), v(*), d(*), ALPHA

      integer             lg, ld, ll, lo, ls, i, j, k, m
      integer             ni, nj, nij, nopt, niop, njop
      integer             ij, ici, icj, iopt, jopt, iold

      double precision    ALFLOG
      double precision    qi, qj, qij, ri, rj, rij, si, sj
      double precision    tracei, tracej, trcij, trop
      double precision    termi, termj, trmij, tmop
      double precision    dij, dopt, siop, sjop

      double precision    zero, one
      parameter          (zero = 0.d0, one = 1.d0)

      double precision    sqrthf
      parameter          (sqrthf = .70710678118654757274d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

      double precision    ddot
      external            ddot

c------------------------------------------------------------------------------

      iopt = 0
      niop = 0
      njop = 0
      nopt = 0
      tmop = 0.d0
      trop = 0.d0

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd

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
        call wardsw( jopt, lg, d)
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
            rj     = dble(nj)
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

      subroutine mevii ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), sigsq(G), pro(G[+1])
      double precision    mu(p,*), sigsq(*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sumz, sum, temp, const, term, zsum
      double precision    sigmin, sigsqk, hold, hood, err
      double precision    prok, tmin, tmax, ViLog, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG
      parameter          (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      zsum  = one

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        if (.not. EQPRO) pro(k) = sumz / dble(n)
        zsum = min(sumz,zsum)
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          sigsqk = zero
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = abs(x(i,j) - mu(j,k))
              if (temp .gt. RTMIN) sum  = sum + temp*temp
            end do
            if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)   
     *          sigsqk = sigsqk + z(i,k)*sum
            z(i,k) = sum
          end do
          sigsq(k) = (sigsqk/sumz)/dble(p)
        else
          sigsq(k) = FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call sgnrng( G, sigsq, 1, sigmin, temp)

      if (sigmin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G
c       temp   = pro(k)
        sigsqk = sigsq(k)
        const  = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
c         z(i,k) = temp*exp(-(const+z(i,k)/sigsqk)/two)           
          z(i,k) = -(const+z(i,k)/sigsqk)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        const = zero - tmax
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) + const
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)-const)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine meviip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

      double precision    Vinv, eps, tol

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), sigsq(G), pro(G[+1])
      double precision    mu(p,*), sigsq(*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sumz, sum, temp, const, term, zsum
      double precision    sigmin, sigsqk, hold, hood, err
      double precision    prok, tmin, tmax, ViLog, rteps
      double precision    pmupmu, cmu, cgam, rmu, rgam

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    SMALOG
      parameter          (SMALOG = -708.d0)

      double precision    ddot, dlngam
      external            ddot, dlngam 

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      rteps  = sqrt(eps)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0
 
      pmupmu = ddot(p,pmu,1,pmu,1)

100   continue

      iter   = iter + 1

      zsum   = one

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        if (.not. EQPRO) pro(k) = sumz / dble(n)
        zsum = min(sumz,zsum)
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          sigsqk = pscale
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = abs(x(i,j) - mu(j,k))
              if (temp .gt. RTMIN) sum  = sum + temp*temp
            end do
            if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)
     *          sigsqk = sigsqk + z(i,k)*sum
          end do
          temp   = pmupmu + ddot(p, mu(1,k), 1, mu(1,k), 1)
          temp   = temp - two*ddot(p,mu(1,k),1,pmu,1) 
          const  = sumz+pshrnk
          sigsqk = sigsqk + ((sumz*pshrnk)/const) * temp
c         sigsq(k) = sigsqk/(pdof+(sumz+one)*dble(p)+two)
          temp   = pdof+sumz*dble(p)+two
          if (pshrnk .gt. zero) temp = temp + dble(p)
          sigsq(k) = sigsqk/temp
          call dscal( p, sumz/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        else
          sigsq(k) = FLMAX
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call sgnrng( G, sigsq, 1, sigmin, temp)

      if (sigmin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G
        sigsqk = sigsq(k)
        const  = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          z(i,k) = -(const+sum/sigsqk)/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        const = zero - tmax
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) + const
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)-const)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      if (pshrnk .gt. zero) then 
        cmu   = dble(p)*(log(pshrnk)-pi2log)/two

        const = pdof/two
        cgam  = const*log(pscale/two)-dlngam(const)

        rmu   = zero
        rgam  = zero
        do k = 1, G
          term = log(sigsq(k))
          temp = pmupmu + ddot( p, mu(1,k), 1, mu(1,k), 1)
          temp = temp - two*ddot( p, mu(1,k), 1, pmu, 1)
          rmu  = rmu + (pshrnk*temp)/sigsq(k)
          rgam = rgam + ((pdof+3.d0)*term - (pscale/sigsq(k)))
        end do

        rmu  = -rmu /two
        rgam = -rgam/two

        pdof  = (dble(G)*cmu+rmu) + (dble(G)*cgam+rgam)
      else
        pdof  = FLMAX
      end if

      return
      end

      subroutine msvii ( x, z, n, p, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), mu(p,G), sigsq(G), pro(G)
      double precision   x(n,*), z(n,*), mu(p,*), sigsq(*), pro(*)
      
      integer                 i, j, k
     
      double precision        sum, sumz, temp, sigsqk

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision       FLMAX
      parameter             (FLMAX = 1.7976931348623157d308)

      double precision       RTMIN
      parameter             (RTMIN = 1.49166814624d-154)

c------------------------------------------------------------------------------

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz / dble(n)
        if (sumz .ge. one .or. one .le. sumz*FLMAX) then  
          call dscal( p, (one/sumz), mu(1,k), 1)
          sigsqk = zero
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = abs(x(i,j) - mu(j,k))
              if (temp .gt. RTMIN) sum  = sum + temp*temp
            end do
            if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN) 
     *          sigsqk = sigsqk + z(i,k)*sum
          end do
          temp = sumz*dble(p)
          if (temp .ge. one .or. sigsqk .le. temp*FLMAX) then  
            sigsq(k) = sigsqk/temp
          else
            sigsq(k) = FLMAX
          end if
        else
          sigsq(k) = FLMAX 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if   
      end do

      return
      end

      subroutine msviip( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer             n, p, G

c     double precision    pshrnk, pmu(p), pscale, pdof
      double precision    pshrnk, pmu(*), pscale, pdof

c     double precision    x(n,p), z(n,G)
      double precision    x(n,*), z(n,*)

c     double precision    mu(p,G), sigsq(G), pro(G)
      double precision    mu(p,*), sigsq(*), pro(*)

      integer             i, j, k

      double precision    sumz, sum, temp
      double precision    sigsqk, const, pmupmu

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    RTMIN
      parameter          (RTMIN = 1.49166814624d-154)

      double precision    ddot
      external            ddot

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      pmupmu = ddot(p,pmu,1,pmu,1)

      do k = 1, G
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz / dble(n)
        if (sumz .ge. one .or. one .lt. sumz*FLMAX) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          sigsqk = pscale
          do i = 1, n
            sum = zero
            do j = 1, p
              temp = abs(x(i,j) - mu(j,k))
              if (temp .gt. RTMIN) sum  = sum + temp*temp
            end do
            if (sqrt(z(i,k))*sqrt(sum) .gt. RTMIN)
     *          sigsqk = sigsqk + z(i,k)*sum
          end do
          temp     = pmupmu + ddot(p, mu(1,k), 1, mu(1,k), 1)
          temp     = temp - two*ddot(p,mu(1,k),1,pmu,1) 
          const    = sumz+pshrnk
          sigsqk   = sigsqk + ((sumz*pshrnk)/const) * temp
          temp     = pdof+sumz*dble(p)+two
          if (pshrnk .gt. zero) temp = temp + dble(p)
          sigsq(k) = sigsqk/temp
          call dscal( p, sumz/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        else 
          sigsq(k) = FLMAX 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

      return
      end

      subroutine esvvi ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), scale(G), shape(p,G), pro(G[+1])
      double precision   mu(p,*), scale(*), shape(p,*), pro(  *  )

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        smin, smax, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. zero) then
        hood = FLMAX
        return
      end if  

      do k = 1, G
        call sgnrng( p, shape(1,k), 1, smin, smax)
        if (smin .le. zero) then
          hood = FLMAX
          return
        end if    
        temp = sqrt(scale(k))
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
      end do

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            if (shape(j,k) .lt. one .and.
     *         abs(temp) .ge. shape(j,k)*FLMAX) then 
              hood = FLMAX
              return
            end if 
            temp = temp/shape(j,k)
            if (abs(temp) .gt. RTMIN) sum = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          hood = FLMAX
          return
        end if
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine mevvi ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), scale(*), shape(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sum, temp, term, scalek, epsmin
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog, zsum, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      tol   = max(tol,zero)
      eps   = max(eps,zero)
      rteps = sqrt(eps)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      zsum  = one

      do k = 1, G
        call dcopy( p, zero, 0, shape(1,k), 1)
        call dcopy( p, zero, 0, mu(1,k), 1)
        sum = zero
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        pro(k) = sum
c pro(k) now contains n_k
        zsum = min(zsum,sum) 
        if (sum .gt. rteps) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum  = sum + temp*temp
            end do
            shape(j,k) = shape(j,k) + sum
          end do
        end if
      end do

      if (zsum .le. rteps) then
        call dscal( G, one/dble(n), pro, 1)
        tol  =  zsum
        eps  = -FLMAX
        maxi =  iter
        return
      end if

      epsmin = FLMAX
      do k = 1, G
        call sgnrng(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .le. zero) then
          scale(k) = zero
        else
          temp = zero
          do j = 1, p
            temp = temp + log(shape(j,k))
          end do
          temp = temp/dble(p)
          if (temp .gt. BIGLOG) then
            call dscal( G, one/dble(n), pro, 1)
            tol  =  zsum 
            eps  =  FLMAX
            maxi =  iter
            return
          end if
          if (temp .gt. SMALOG) then
            temp = exp(temp)
          else
            temp = zero
          end if
          scale(k) = temp/pro(k)
          epsmin   = min(temp,epsmin)
          if (temp .le. eps) then
            call dscal( G, one/dble(n), pro, 1)
            tol  =  zsum 
            eps  =  FLMAX
            maxi =  iter
            return
          end if
          call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      if (.not. EQPRO) then
        call dscal( G, one/dble(n), pro, 1)
      else if (Vinv .le. zero) then
        call dscal( G, one/dble(G), pro, 1)
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (epsmin .le. eps) then
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call sgnrng( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = err
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      do k = 1, G
        scalek = scale(k)
        const = dble(p)*(pi2log + log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j,k)
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/scalek))/two)
          z(i,k) = -(const+(sum/scalek))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine mevvip( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

c     double precision   pshrnk, pmu(p), pscale, pdof
      double precision   pshrnk, pmu(*), pscale, pdof

      double precision    Vinv, eps, tol

c     double precision    x(n,p), z(n,G[+1])
      double precision    x(n,*), z(n,  *  )

c     double precision    mu(p,G), scale(G), shape(p,G), pro(G[+1])
      double precision    mu(p,*), scale(*), shape(p,*), pro(  *  )

      integer             nz, iter, i, j, k

      double precision    sumz, sum, temp, term, scalek, epsmin
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog, zsum, rteps

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

      rteps = sqrt(eps) 

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      zsum  = one 

      do k = 1, G
        call dcopy( p, pscale, 0, shape(1,k), 1)
        call dcopy( p, zero, 0, mu(1,k), 1)
        sumz = zero
        do i = 1, n
          temp   = z(i,k)
          sumz   = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        pro(k) = sumz
        zsum   = min(zsum,sumz)
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          term   = pshrnk+sumz
          const  = (pshrnk*sumz)/term
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum  = sum + temp*temp
            end do
            shape(j,k) = shape(j,k) + sum
            temp       = pmu(j) - mu(j,k)
            shape(j,k) = shape(j,k) + const*(temp*temp)
          end do
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
        end if
      end do

      if (zsum .le. rteps) then
        call dscal( G, one/dble(n), pro, 1)
        tol  =  zsum 
        eps  = -FLMAX
        maxi =  iter
        return
      end if

c pro(k) now contains n_k

      epsmin = FLMAX
      term   = pdof+two
      if (pshrnk .gt. zero) term = term + one 
      do k = 1, G
        call sgnrng(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .eq. zero) then
          scale(k) = zero
        else
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp = sum/dble(p)
          if (temp .gt. BIGLOG) then
            call dscal( G, one/dble(n), pro, 1)
            tol  =  zsum 
            eps  =  FLMAX
            maxi =  iter
            return
          end if
          if (temp .gt. SMALOG) then
            temp = exp(temp)
          else
            temp = zero
          end if
c pro(k) contains n_k
          scale(k) = temp/(pro(k)+term)
          epsmin   = min(temp,epsmin)
          if (temp .le. eps) then
            call dscal( G, one/dble(n), pro, 1)
            tol  =  zsum 
            eps  =  FLMAX
            maxi =  iter
            return
          end if
          call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      if (.not. EQPRO) then
        call dscal( G, one/dble(n), pro, 1)
      else if (Vinv .le. zero) then
        call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      if (epsmin .le. eps) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call sgnrng( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call sgnrng( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = err
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      do k = 1, G
        scalek = scale(k)
        const = dble(p)*(pi2log + log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j,k)
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/scalek))/two)
          z(i,k) = -(const+(sum/scalek))/two
        end do
      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine msvvi ( x, z, n, p, G, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G)
      double precision   x(n,*), z(n,*)

c     double precision   mu(p,G), scale(G), shape(p,G), pro(G)
      double precision   mu(p,*), scale(*), shape(p,*), pro(*)

      integer                 i, j, k

      double precision        sum, temp, smin, smax

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      do k = 1, G
        call dcopy( p, zero, 0, shape(1,k), 1)
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        pro(k) = sum
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
      end do

c pro(k) now contains n_k

      do k = 1, G
        if (mu(1,k) .ne. FLMAX) then
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum  = sum + temp*temp
            end do
            shape(j,k) = shape(j,k) + sum
          end do
        else
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if
      end do

      do k = 1, G 

        call sgnrng(p, shape(1,k), 1, smin, smax)

        if (smin .le. zero) then
          scale(k) = zero
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        if (smax .eq. FLMAX) then 
          scale(k) = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j,k))
        end do      
        temp = sum/dble(p)

        if (temp .gt. BIGLOG) then
          scale(k) = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        if (temp .lt. SMALOG) then
          temp     = zero
          scale(k) = zero
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        temp = exp(temp)
        if (pro(k) .lt. one .and. temp .ge. pro(k)*FLMAX) then
          scale(k) = FLMAX
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if

        scale(k) = temp/pro(k)
        if (temp .lt. one .and. one .ge. temp*FLMAX) then
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
          goto 100
        end if
            
        call dscal( p, one/temp, shape(1,k), 1)

 100    continue

      end do

      call dscal( G, one/dble(n), pro, 1)

      return
      end

      subroutine msvvip( x, z, n, p, G, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   pshrnk, pmu(p), pscale, pdof
      double precision   pshrnk, pmu(*), pscale, pdof

c     double precision   x(n,p), z(n,G)
      double precision   x(n,*), z(n,*)

c     double precision   mu(p,G), scale(G), shape(p,G), pro(G)
      double precision   mu(p,*), scale(*), shape(p,*), pro(*)

      integer             i, j, k

      double precision    sumz, sum, temp, term
      double precision    smin, smax, const

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMALOG, BIGLOG
      parameter          (SMALOG = -708.d0, BIGLOG = 709.d0)

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      do k = 1, G
        call dcopy( p, pscale, 0, shape(1,k), 1)
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sumz   = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
          z(i,k) = sqrt(temp)
        end do
        pro(k) = sumz
        if (sumz .ge. one .or. one .le. sumz*FLMAX) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          term   = pshrnk+sumz
          const  = (pshrnk*sumz)/term
          do j = 1, p
            sum = zero
            do i = 1, n
              temp = z(i,k)*(x(i,j) - mu(j,k))
              sum  = sum + temp*temp
            end do
            shape(j,k) = shape(j,k) + sum
            temp       = pmu(j) - mu(j,k)
            shape(j,k) = shape(j,k) + const*(temp*temp)
          end do
          call dscal( p, sumz/term, mu(1,k), 1)
          call daxpy( p, pshrnk/term, pmu, 1, mu(1,k), 1)
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        end if
      end do

c pro(k) now contains n_k

      do k = 1, G
        call sgnrng(p, shape(1,k), 1, smin, smax)
        if (smin .le. zero) then
          scale(k) = zero
          call dcopy( p, FLMAX, 0, shape(1,k), 1)
        else if (smax .eq. FLMAX) then
          scale(k) = FLMAX
        else
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp = sum/dble(p)
          if (temp .gt. BIGLOG) then
            scale(k) = FLMAX
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
          else if (temp .lt. SMALOG) then
            temp     = zero
            scale(k) = zero
            call dcopy( p, FLMAX, 0, shape(1,k), 1)
          else
            temp     = exp(temp)
c pro(k) contains n_k
            term = pro(k) + pdof + two
            if (pshrnk .gt. zero) term = term + one
            scale(k) = temp/term
            if (temp .ge. one .or. one .le. temp*FLMAX) then
              call dscal( p, one/temp, shape(1,k), 1)
            else
              call dcopy( p, FLMAX, 0, shape(1,k), 1)
            end if
          end if
        end if
      end do
 
      call dscal( G, one/dble(n), pro, 1)

      return
      end

      subroutine esvvv ( CHOL, x, mu, Sigma, pro, n, p, G, Vinv, 
     *                   w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     character          CHOL
      logical            CHOL

c     integer            n, p, G
      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), w(p), z(n,G[+1])
      double precision   x(n,*), w(*), z(n,  *  )

c     double precision   mu(p,G), Sigma(p,p,G), pro(G[+1])
      double precision   mu(p,*), Sigma(p,p,*), pro(  *  )

      integer                 nz, p1, info, i, j, k

      double precision        const, detlog, temp, prok, tmin, tmax
      double precision        umin, umax, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      p1    = p + 1

c     if (CHOL .eq. 'N') then
      if (.not. CHOL) then

        do k = 1, G

          call dpotrf( 'U', p, Sigma(1,1,k), p, info)

          w(1) = dble(info)

          if (info .ne. 0) then
            hood = FLMAX
            return
          end if
       
        end do

      end if

      do k = 1, G
        
        call absrng( p, Sigma(1,1,k), p1, umin, umax)

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

      end do

      do k = 1, G

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

c       prok  = pro(k)

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
          temp   = ddot( p, w, 1, w, 1)/two
c         z(i,k) = prok*exp(-(const+temp))
          z(i,k) = -(const+temp)
        end do

      end do

      w(1) = zero
      if (pro(1) .lt. zero) return

      nz = G
      if (Vinv .gt. zero) then
        nz = nz + 1
c       call dcopy( n, pro(nz)*Vinv, 0, z(1,nz), 1)
        call dcopy( n, log(Vinv), 0, z(1,nz), 1)
      end if

c     hood = zero
c     do i = 1, n
c       sum = zero
c       do k = 1, nz
c         sum = sum + z(i,k)
c       end do
c       hood = hood + log(sum)
c       call dscal( nz, (one/sum), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood = hood + (log(sum)+tmax)
        if (sum .lt. one .and. one .ge. sum*FLMAX) then
          w(1) = zero
          hood = FLMAX
          return
        end if 
        call dscal( nz, (one/sum), z(i,1), n)
      end do

      w(1) = zero

      return
      end

       subroutine hcvvv ( x, n, p, ic, ng, ns, ALPHA, BETA, 
     *                    v, u, s, r, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, ic(n), ng, ns, nd

      double precision   ALPHA, BETA

c     double precision   x(n,p+1), v(p), u(p,p), s(p,p)
c     double precision   r(p,p), d(ng*(ng-1)/2)
      double precision   x(n,*), v(*), u(p,*), s(p,*)
      double precision   r(p,*), d(*)

      integer            psq, pm1, pp1
      integer            i, j, k, l, m, ij, iold
      integer            lg, ld, ll, lo, ls
      integer            ici, icj, ni, nj, nij
      integer            nopt, niop, njop, iopt, jopt

      double precision   trcij, trmij, trop, tmop
      double precision   traci, tracj, termi, termj
      double precision   qi, qj, qij, si, sj, sij, ri, rj, rij
      double precision   dij, dopt, siop, sjop

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   rthalf
      parameter         (rthalf = .7071067811865476d0)

      external           ddot, vvvtij
      double precision   ddot, vvvtij

      double precision   BETA0, ALPHA0, ABLOG
      common /VVVMCL/    BETA0, ALPHA0, ABLOG
      save   /VVVMCL/            

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

c------------------------------------------------------------------------------

      iopt = 0
      niop = 0
      nopt = 0
      tmop = 0.d0
      trop = 0.d0

      lg     =  ng
      ld     = (ng*(ng-1))/2
      ll     =  nd-ng
      lo     =  nd

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
        call wardsw( jopt, lg, d)
        call dcopy( p, x(lg,1), n, x(jopt,1), n)
        m          = ic(jopt)
        icj        = ic(lg)
        if (icj .ne. 1) x( jopt, pp1) = x( lg, pp1)
        ic(jopt)   = icj
        ic(lg)     = m
      end if

      if (niop .eq. 1) then
        ic(iopt) = lg
      else 
        l = ic(iopt)
        do k = 1, min(niop-1,p)
          m = l
          l = ic(l)
        end do
        if (l .lt. n) call intpr("l .lt. n", 8, l, 1)
        ic(m) = lg
      end if

      l = ic(iopt)
      do k = 1, min(nopt-1,p)
        call dcopy( p, r(1,1), p, x(l,1), n)
        m = l
        l = ic(l)
      end do
      ic(m) = nopt + n        

c     call intpr('ic', 2, ic, n)

c     term(iopt) = tmop
c     trac(iopt) = trop

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

      iold  =  iopt

      dopt  = FLMAX

      ni    = nopt
      ri    = dble(ni)
      termi = tmop
      traci = trop

      ij = ((iold-1)*(iold-2))/2
      if (iold .gt. 1) then
        do j = 1, (iold-1)
          call dcopy(psq, zero, 0, u, 1)
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
              k  = k  + 1
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
          call dcopy(psq, zero, 0, u, 1)
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
            k  = k + 1
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

c     call dblepr("d", 1, d, nd)
c     call dblepr("d", 1, d, ld)

      if (ij .gt. 1) then
        do i = 2, ij
          iopt = iopt + 1
          if (iopt .ge. jopt) then
            jopt = jopt + 1
            iopt = 1
          end if
        end do
      end if

      do k = 1, p
        call dcopy( p, zero, 0, r(1,k), 1)
      end do

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

        l = ic(iopt)
        if (l .ne. 1) then
710       continue
          if (l .le. n) then
            l = ic(l)
            goto 710
          end if
          niop = l-n
        else
          niop = 1
        end if

        l = ic(jopt)
        if (l .ne. 1) then
720       continue
          if (l .le. n) then
            l = ic(l)
            goto 720
          end if
          njop = l-n
        else
          njop = 1
        end if

        nopt = niop + njop        
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

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer                    l, p

      double precision           r(p,*), s, trac

      double precision           detlog

      double precision           zero, one
      parameter                 (zero = 0.d0, one = 1.d0)

      external                   det2mc
      double precision           det2mc

      double precision           BETA, ALPHA, ABLOG
      common /VVVMCL/            BETA, ALPHA, ABLOG
      save   /VVVMCL/            

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      if (l .le. p) then
        vvvtij = log(BETA*(trac+ALPHA)/dble(l))
      else 
        if (trac .eq. zero) then
          vvvtij = log((ALPHA*BETA)/dble(l))
        else
          detlog = det2mc( p, r, s)
          if (detlog .eq. -FLMAX) then
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

      double precision function det2mc( n, u, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer                          k, n

      double precision                 q, s
      double precision                 u(n,*)

      double precision                 zero, two
      parameter                       (zero = 0.d0, two = 2.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      det2mc = zero

      do k = 1, n

        q = u(k,k)*s

        if (abs(q) .le. zero) then
          det2mc = -FLMAX
          return
        end if

        det2mc = det2mc + log(abs(q))

      end do

      det2mc = two*det2mc

      return
      end

      subroutine mevvv ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, U, pro, w, S)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p,G), pro(G), S(p,p)
      double precision   mu(p,*), U(p,p,*), pro(*), S(p,*)

      integer                 nz, p1, iter, i, j, k, l, j1

      double precision        piterm, hold, rcmin, rteps
      double precision        temp, term, cs, sn, umin, umax
      double precision        sumz, sum, detlog, const, hood, err
      double precision        prok, tmin, tmax, ViLog, zsum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

c zero out the lower triangle
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, S(j,i), 1)
          i = j
        end do
        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do
      end do

      iter = 0

100   continue

      iter = iter + 1

      zsum = one
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do
        do j = 1, p
          call dcopy( j, zero, 0, S(1,j), 1)
        end do
        call dcopy( p, zero, 0, mu(1,k), 1)
        sumz = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        if (.not. EQPRO) pro(k) = sumz / dble(n)
        zsum = min(sumz,zsum) 
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( S(j,j), w(j), cs, sn)
              call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( S(p,p), w(p), cs, sn)
          end do
          do j = 1, p
            call dscal( j, one/sqrt(sumz), S(1,j), 1)
          end do
        else
          call dcopy( p, FLMAX, 0,  z(1,k), 1) 
        end if
        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do
      end do

      if (zsum .le. rteps) then
        tol  = zsum
        eps  = -FLMAX
        maxi = iter
        return
      end if

      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      rcmin = FLMAX
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,K)
          end do
        end do
        call absrng( p, S, p1, umin, umax)
        rcmin = min(umin/(one+umax),rcmin)
      end do

      if (rcmin .le. rteps) then
        tol  = rcmin
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,K)
          end do
        end do

c       temp = pro(k)

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(S(j,j)))
        end do

        const = piterm+detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, S, p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do

      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol .and. iter .lt. maxi) goto 100

c     w(1) = rcmin

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine mevvvp( EQPRO, x, n, p, G, Vinv, 
     *                   pshrnk, pmu, pscale, pdof,
     *                   z, maxi, tol, eps, mu, U, pro, w, S)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   Vinv, eps, tol

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p,G), pro(G), S(p,p)
      double precision   mu(p,*), U(p,p,*), pro(*), S(p,*)

      integer                 nz, p1, iter, i, j, k, l, j1

      double precision        piterm, hold, rcmin, rteps
      double precision        temp, term, cs, sn, umin, umax
      double precision        sum, sumz, detlog, const, hood, err
      double precision        prok, tmin, tmax, ViLog
      double precision        cmu, cgam, rmu, rgam, zsum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        twolog
      parameter              (twolog  = 0.6931471805599453d0)

      double precision        pilog
      parameter              (pilog = 1.144729885849400d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot, dlngam
      double precision        ddot, dlngam

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter = 0

100   continue

      iter = iter + 1

      zsum = one
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, S(1,j), 1)
        end do
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        if (.not. EQPRO) pro(k) = sumz / dble(n)
        zsum = min(sumz,zsum)
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( S(j,j), w(j), cs, sn)
              call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( S(p,p), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          const = sumz+pshrnk
          temp  = (sumz*pshrnk)/const
          call dscal( p, sqrt(temp), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( S(j,j), w(j), cs, sn)
            call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( S(p,p), w(p), cs, sn)
          do j = 1, p
            temp = pdof+sumz+dble(p)+two
            call dscal( j, one/sqrt(temp), S(1,j), 1)
          end do
          call dscal( p, sumz/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        else
          call dcopy( p, FLMAX, 0,  z(1,k), 1) 
        end if
        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do
      end do

      if (zsum .le. rteps) then
        tol  = zsum
        eps  = -FLMAX
        maxi = iter
        return
      end if

      term = zero
      if (Vinv .gt. zero) then

        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      rcmin = FLMAX
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,K)
          end do
        end do
        call absrng( p, S, p1, umin, umax)
        rcmin = min(umin/(one+umax),rcmin)
      end do

      if (rcmin .le. rteps) then
        tol  = rcmin
        eps  = FLMAX
        maxi = iter
        return
      end if

      rmu  = zero
      rgam = zero 
      do k = 1, G

c       temp = pro(k)
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(S(j,j)))
        end do

        rmu  = rmu - detlog
        rgam = rgam - (pdof+dble(p)+one)*detlog

        const = piterm+detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, S, p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do

      end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do
        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol .and. iter .lt. maxi) goto 100

c     w(1) = rcmin
 
      tol  = err
      eps  = hood
      maxi = iter

      cmu  = dble(p)*(log(pshrnk) - pi2log)/two

      sum = zero
      do k = 1, G
        do j = 1, p
          do l = 1, p
             S(l,j) = U(l,j,k)
          end do
        end do
        call daxpy( p, (-one), mu(1,k), 1, pmu, 1)
        call dtrsv('U','T','N',p, S, p, pmu, 1)
        sum = sum + ddot( p, pmu, 1, pmu, 1)
      end do

      rmu = rmu - pshrnk*sum/two

      sum    = zero
      term   = zero
      temp   = zero
      do j = 1, p
        call dcopy( p, pscale(j,1), p, pmu, 1)
c       call dtrsv('U','T','N', p, U, p, pmu, 1)
        i = p-j+1
c       call dtrsv('U','T','N', i, U(j,j,k),i,pmu(j),1)
        call dtrsv('U','T','N', i, S(j,j), p, pmu(j), 1)
        sum    = sum + ddot(i, pmu(j), 1, pmu(j), 1)
        temp   = temp + log(abs(pscale(j,j)))
        term   = term + dlngam((pdof+one-dble(j))/two)
      end do

      rgam  = rgam - sum/two

      const = -dble(p)*(pdof*twolog+(dble(p)-one)*pilog/two)
      cgam  = (const-pdof*temp)/two-term

      pdof  = (dble(G)*cmu+rmu) + (dble(G)*cgam+rgam)

      return
      end

      subroutine msvvv ( x, z, n, p, G, w, mu, U, pro, S)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p,G), pro(G), S(p,p)
      double precision   mu(p,*), U(p,p,*), pro(*), S(p,*)

      integer                 i, j, k, l, j1

      double precision        sum, temp, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
c         call dcopy( j, zero, 0, U(1,j,k), 1)
          call dcopy( j, zero, 0, S(1,j), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sum / dble(n)
        if (sum .ge. one .or. one .lt. sum*FLMAX) then
          call dscal( p, (one/sum), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( S(j,j), w(j), cs, sn)
              call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( S(p,p), w(p), cs, sn)
          end do
          temp = sqrt(sum)
          if (temp .ge. one .or. one .lt. temp*FLMAX) then
            do j = 1, p
              call dscal( j, one/temp, S(1,j), 1)
            end do
          else
            do j = 1, p
              call dcopy( j, FLMAX, 0, S(1,j), 1)
            end do
          end if
        else
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do
      end do

      return
      end

      subroutine msvvvp( x, z, n, p, G,
     *                   pshrnk, pmu, pscale, pdof,
     *                   w, mu, U, pro, S)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n, p, G

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p,G), pro(G), S(p,p)
      double precision   mu(p,*), U(p,p,*), pro(*), S(p,*)

c------------------------------------------------------------------------------
c
c  x       double  (input) (n,p) matrix of observations.
c  z       double  (input/output) (n,G[+1]) conditional probabilities. 
c  n       integer (input) number of observations.
c  p       integer (input) dimension of the data.
c  G       integer (input) number of Gaussian clusters in the mixture.
c  mu      double  (output) (p,G) mean for each group.
c  U       double  (output) (p,p,G)
c  pro     double  (output) (G) mixing proportions (used even if equal).
c  w       double  (scratch) (max(p,G))

      integer                 i, j, k, l, j1

      double precision        sumz, temp, cs, sn, const

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do
        do j = 1, p
          call dcopy( p, pscale(1,j), 1, S(1,j))
        end do
        sumz = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        pro(k) = sumz / dble(n)
        if (sumz .ge. one .or. one .lt. sumz*FLMAX) then
          call dscal( p, (one/sumz), mu(1,k), 1)
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( S(j,j), w(j), cs, sn)
              call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( S(p,p), w(p), cs, sn)
          end do
          call dcopy( p, pmu, 1, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          const = sumz+pshrnk
          temp  = (sumz*pshrnk)/const
          call dscal( p, sqrt(temp), w, 1)
          j = 1
          do j1 = 2, p
            call drotg( S(j,j), w(j), cs, sn)
            call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
            j = j1
          end do
          call drotg( S(p,p), w(p), cs, sn)
          temp = pdof+sumz+dble(p)+one
          if (pshrnk .gt. zero) temp = temp + one
          do j = 1, p
            call dscal( j, one/sqrt(temp), S(1,j), 1)
          end do
          call dscal( p, sumz/const, mu(1,k), 1)
          call daxpy( p, pshrnk/const, pmu, 1, mu(1,k), 1)
        else 
          call dcopy( p, FLMAX, 0, mu(1,k), 1)
        end if
        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do
      end do

      return
      end

      subroutine mvn1d ( x, n, mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n
      integer            n

      double precision   mu,  sigsq, hood

c     double precision   x(n)
      double precision   x(*)

c------------------------------------------------------------------------------
c
c  x       double  (input) (n) matrix of observations (destroyed).
c  n       integer (input) number of observations.
c  mu      double  (output) mean.
c  sigsq   double  (output) variance.
c  hood    double  (output) loglikelihood.

      double precision        dn

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot
      external                ddot

c------------------------------------------------------------------------------
      
      dn    = dble(n)

      mu    = ddot( n, one/dn, 0, x, 1)

      sigsq = zero

      call daxpy( n, (-one), mu, 0, x, 1)
      sigsq = ddot( n, x, 1, x, 1)/dn
 
      if (sigsq .eq. zero) then
        hood = FLMAX
      else 
        hood = -dn*(pi2log + (one + log(sigsq)))/two
      end if

      return
      end

      subroutine mvn1p ( x, n, pshrnk, pmu, pscale, pdof,
     *                   mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      integer            n

      double precision   pshrnk, pmu, pscale, pdof
      double precision   mu, sigsq, hood

c     double precision   x(n)
      double precision   x(*)

      integer                 i

      double precision   dn, scl, const, term, temp, xbar
      double precision   cmu, cgam, rmu, rgam

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   pi2log
      parameter         (pi2log = 1.837877066409345d0)

      double precision   FLMAX
      parameter         (FLMAX = 1.7976931348623157d308)

      double precision        ddot, dlngam
      external                ddot, dlngam

c------------------------------------------------------------------------------

      if (pshrnk .lt. zero) pshrnk = zero

      dn    = dble(n)

      xbar  = ddot( n, one/dn, 0, x, 1) 
      const = pshrnk + dn
      mu    = (dn/const)*xbar + (pshrnk/const)*pmu
      sigsq = zero
      do i = 1, n
        temp  = xbar - x(i)
        sigsq = sigsq + temp*temp
      end do     
      temp  = xbar - pmu
      sigsq = sigsq + pscale + dn*(pshrnk/const)*(temp*temp)
      temp = pdof + dn + two 
      if (pshrnk .gt. zero) temp = temp + one
      sigsq = sigsq / temp
      
      if (sigsq .eq. zero) then
        hood = FLMAX
      else 
        call daxpy( n, (-one), mu, 0, x, 1)
        temp = ddot( n, x, 1, x, 1)
        if (sigsq .lt. one .and. temp .ge. sigsq*FLMAX) then
          hood = FLMAX
          return
        end if
        temp = temp/sigsq
        hood = -(dn*(pi2log + log(sigsq)) + temp)/two
      end if

      if (pshrnk .gt. zero) then
        cmu    = (pi2log-log(pshrnk))/two
        term   =  pdof/two
        cgam   =  term*log(pscale/two) - dlngam(term)
        temp   =  pmu - mu
        const  =  log(sigsq)
        rmu    = -(const - (pshrnk/sigsq)*(temp*temp))/two
        rgam   = -(term+one)*const - (pscale/sigsq)/two

        pdof   = (cmu+rmu) + (cgam+rgam)
      else
        pdof   = FLMAX
      end if

      return
      end

      subroutine mnxiip( x, n, p, pshrnk, pmu, pscale, pdof, 
     *                   mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

c     double precision   pshrnk, pmu(p), pscale, pdof
      double precision   pshrnk, pmu(*), pscale, pdof

      double precision   sigsq, hood

c     double precision   x(n,p), mu(p)
      double precision   x(n,*), mu(*)

      integer                 i, j

      double precision        dnp, scl, temp, term, sum
      double precision        dmudmu, pmupmu, cmu, cgam, rmu, rgam

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot, dlngam
      external                ddot, dlngam

c------------------------------------------------------------------------------
      
      dnp = dble(n*p)

      scl = one/dble(n)
      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
      end do

      sum = zero
      do i = 1, n
        do j = 1, p
          temp = x(i,j) - mu(j)
          sum  = sum + temp*temp
        end do
      end do

      pmupmu = ddot(p,pmu,1,pmu,1)
      dmudmu = ddot(p,mu,1,mu,1)
      temp   = dmudmu + pmupmu
      temp   = temp - two*ddot(p,mu,1,pmu,1)

      term   = pshrnk + dble(n)
      scl    = (pshrnk*dble(n))/term
      sigsq  = pscale + scl*temp + sum

      temp   = pdof + dble(n*p) + two
      if (pshrnk .gt. zero) temp = temp + dble(p)
      sigsq  = sigsq/temp

      call dscal( p, dble(n)/term, mu, 1)
      call daxpy( p, pshrnk/term, pmu, 1, mu, 1)

      if (sigsq .eq. zero) then
        hood  = FLMAX
      else 
        sum  = zero
        do i = 1, n
          do j = 1, p
            temp = x(i,j) - mu(j)
            sum  = sum + temp*temp
          end do
        end do
        hood  = -(sum/sigsq + dnp*(pi2log + log(sigsq)))/two
      end if

      if (pshrnk .gt. zero) then 
        dmudmu = ddot(p,mu,1,mu,1)

        cmu    = dble(p)*(log(pshrnk)-pi2log)/two
        temp   = (dmudmu+pmupmu) - two*ddot(p,pmu,1,mu,1)
        term   = log(sigsq)
        rmu    = -(dble(p)*term + (pshrnk*temp)/sigsq)/two

        temp   = pdof/two     
        cgam   = temp*log(pscale/two) - dlngam(temp)
        rgam   = -(temp+one)*term - pscale/(two*sigsq)

        pdof   = (cmu+rmu) + (cgam+rgam)
      else
        pdof   = FLMAX
      end if

      return
      end

      subroutine mvnxii( x, n, p, mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

      double precision   sigsq, hood

c     double precision   x(n,p), mu(p)
      double precision   x(n,*), mu(*)

      integer                 j

      double precision        dnp, scl

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot
      external                ddot

c------------------------------------------------------------------------------
      
      dnp = dble(n*p)

      scl = one/dble(n)
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
        hood  = FLMAX
      else 
        hood  = -dnp*(pi2log + (one + log(sigsq)))/two
      end if

      return
      end

      subroutine mnxxip( x, n, p, pshrnk, pmu, pscale, pdof, 
     *                   mu, scale, shape, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

c     double precision   pshrnk, pmu(p), pscale, pdof
      double precision   pshrnk, pmu(*), pscale, pdof

      double precision   scale, hood

c     double precision   x(n,p), mu(p), shape(p)
      double precision   x(n,*), mu(*), shape(*)

      integer            i, j

      double precision   sum, temp, smin, smax
      double precision   term, const, scl

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision   pi2log
      parameter         (pi2log = 1.837877066409345d0)

      double precision   FLMAX
      parameter         (FLMAX = 1.7976931348623157d308)

      double precision   SMALOG, BIGLOG
      parameter         (SMALOG = -708.d0, BIGLOG = 709.d0)

      double precision   ddot
      external           ddot

c------------------------------------------------------------------------------

      temp = one/dble(n)
      do j = 1, p
        mu(j)    = ddot( n, temp, 0, x(1,j), 1)
        shape(j) = zero
      end do

      do j = 1, p
        sum = zero
        do i = 1, n
          temp = x(i,j) - mu(j)
          sum  = sum + temp*temp
        end do
        shape(j) = shape(j) + sum
      end do

      term = pshrnk + dble(n)
      scl  = (pshrnk*dble(n))/term
      do j = 1, p
        temp     = pmu(j) - mu(j)
        shape(j) = shape(j) + scl*(temp*temp) + pscale
      end do

      call dscal( p, dble(n)/term, mu, 1)
      call daxpy( p, pshrnk/term, pmu, 1, mu, 1)

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        scale = zero
        hood  = FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp  = sum/dble(p)

      if (temp .ge. BIGLOG) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        scale = FLMAX
        hood  = FLMAX
       return
      end if

      if (temp .le. SMALOG) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        scale = zero
        hood  = FLMAX
       return
      end if

      temp  = exp(temp)

      term = pdof + dble(n) + two
      if (pshrnk .gt. zero) term = term + one

      scale = temp/term

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1) 
        hood = FLMAX
        return
      end if 
 
      call dscal( p, one/temp, shape, 1)

      const = dble(p)*(pi2log+log(scale))

      hood  = zero
      do i = 1, n
        sum = zero
        do j = 1, p
          temp = x(i,j) - mu(j)
          sum  = sum + (temp*temp)/shape(j)             
        end do    
        hood = hood - (const+(sum/scale))/two
      end do

c log posterior computation not yet available
      pdof = FLMAX 

      return
      end

      subroutine mvnxxi( x, n, p, mu, scale, shape, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

      double precision   scale, hood

c     double precision   x(n,p), mu(p), shape(p)
      double precision   x(n,*), mu(*), shape(*)

      integer                 i, j

      double precision        dn, scl, sum, temp, smin, smax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG, BIGLOG
      parameter              (SMALOG = -708.d0, BIGLOG = 709.d0)

      double precision        ddot
      external                ddot

c------------------------------------------------------------------------------

      dn  = dble(n)

      scl = one/dn
      do j = 1, p
        mu(j)    = ddot( n, scl, 0, x(1,j), 1)
        shape(j) = zero
      end do

      do j = 1, p
        sum = zero
        do i = 1, n
          temp = x(i,j) - mu(j)
          sum  = sum + temp*temp
        end do
        shape(j) = shape(j) + sum
      end do

      call sgnrng(p, shape, 1, smin, smax)

      if (smin .le. zero) then
        call dcopy( p, FLMAX, 0, shape, 1)
        scale = zero
        hood  = FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do

      temp  = sum/dble(p)

      if (temp .gt. BIGLOG) then
        call dcopy( p, FLMAX, 0, shape, 1)
        scale = FLMAX
        hood  = FLMAX
        return
      end if

      if (temp .lt. SMALOG) then
        call dcopy( p, FLMAX, 0, shape, 1)
        scale = zero
        hood  = FLMAX
        return
      end if

      temp  = exp(temp)

      scale = temp/dn

      if (temp .lt. one .and. one .ge. temp*FLMAX) then
        call dcopy( p, FLMAX, 0, shape, 1)
        hood = FLMAX
        return
      end if 

      call dscal( p, one/temp, shape, 1)

      hood  = -dble(n*p)*(one + pi2log + log(scale))/two

      return
      end

      subroutine mnxxxp( x, n, p, w,
     *                   pshrnk, pmu, pscale, pdof,
     *                   mu, U, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

c     double precision   pshrnk, pmu(p), pscale(p,p), pdof
      double precision   pshrnk, pmu(*), pscale(p,*), pdof

      double precision   hood

c     double precision   x(n,p), w(p), mu(p), U(p,p)
      double precision   x(n,*), w(*), mu(*), U(p,*)

      integer                 i, j, j1

      double precision        dnp, scl, detlog, sum, term, temp
      double precision        umin, umax, cs, sn, const
      double precision        cmu, cgam, rmu, rgam

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        twolog
      parameter              (twolog  = 0.6931471805599453d0)

      double precision        pilog
      parameter              (pilog = 1.144729885849400d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot, dlngam
      external                ddot, dlngam

c------------------------------------------------------------------------------

      dnp = dble(n*p)

      scl = one/dble(n)
      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
        call dcopy( p, pscale(1,j), 1, U(1,j), 1)
      end do
c mu contains ybar; U contains Cholesky factor of inverse Wishart scale

      do i = 1, n
        call dcopy( p, x(i,1), n, w, 1)
        call daxpy( p, (-one), mu, 1, w, 1)
        j = 1
        do j1 = 2, p 
          call drotg( U(j,j), w(j), cs, sn)
          call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
          j = j1
        end do
        call drotg( U(p,p), w(p), cs, sn)
      end do

      call dcopy( p, pmu, 1, w, 1)
      call daxpy( p, (-one), mu, 1, w, 1)
      term = (pshrnk*dble(n))/(pshrnk+dble(n))
      call dscal( p, sqrt(term), w, 1)
      j = 1
      do j1 = 2, p
        call drotg( U(j,j), w(j), cs, sn)
        call drot( p-j, U(j,j1), p, w(j1), 1, cs, sn)
        j = j1
      end do
      call drotg( U(p,p), w(p), cs, sn)

      scl = pdof + dble(n+p+1)
      if (pshrnk .gt. zero) scl = scl + one
      scl = one/sqrt(scl)

      do j = 1, p
        call dscal( j, scl, U(1,j), 1)
      end do

      term = pshrnk + dble(n)

      call dscal( p, dble(n)/term, mu, 1)
      call daxpy( p, pshrnk/term, pmu, 1, mu, 1)

      call absrng( p, U, p+1, umin, umax)

c     rcond = umin / (one + umax)

      if (umin .eq. zero) then
        hood  = FLMAX
        return
      end if 

      detlog  = zero
      do j = 1, p
        detlog = detlog + log(abs(U(j,j)))
      end do
      const = dble(n)*(detlog + dble(p)*pi2log/two)
      sum = zero
      do i = 1, n
        call dcopy( p, x(i,1), n, w, 1)
        call daxpy( p, (-one), mu, 1, w, 1)
        call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
        sum = sum + ddot(p, w, 1, w, 1)
      end do
      hood = -(const+sum/two)

      cmu = dble(p)*(log(pshrnk) - pi2log)/two

      call dcopy( p, pmu, 1, w, 1)
      call daxpy( p, (-one), mu, 1, w, 1)
      call dtrsv( 'U', 'T', 'N', p, U, p, w, 1)
      temp = ddot( p, w, 1, w, 1)

      sum    = zero
      term   = zero
      do j = 1, p
        term = term + dlngam((pdof+dble(1-j))/two)
        call dcopy( p, pscale(j,1), p, pmu, 1)
c       call dtrsv('U','T','N', p, U, p, pmu, 1)
        i = p-j+1
c       call dtrsv('U','T','N', i, U(j,j),i,pmu(j),1)
        call dtrsv('U','T','N', i, U(j,j),p,pmu(j),1)
        sum  = sum + ddot(i, pmu(j), 1, pmu(j), 1)
      end do

      if (pshrnk .gt. zero) then 
        rmu   = -(detlog+pshrnk*temp/two)

        const = -dble(p)*(pdof*twolog+(dble(p)-one)*pilog/two)
        cgam  = (const/two-pdof*detlog) - term

        rgam = -((pdof+dble(p)+one)*detlog + sum/two)

        pdof  = (cmu+cgam) + (rmu+rgam)
      else
        pdof = FLMAX
      end if

      return
      end

      subroutine mvnxxx( x, n, p, mu, U, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

c     integer            n, p
      integer            n, p

      double precision   hood

c     double precision   x(n,p), mu(p), U(p,p)
      double precision   x(n,*), mu(*), U(p,*)

      integer                 i, j, j1

      double precision        dn, dnp, scl
      double precision        umin, umax, cs, sn

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot
      external                ddot

c------------------------------------------------------------------------------

      dn  = dble(n)
      dnp = dble(n*p)

      scl = one/dn
      do j = 1, p
        mu(j) = ddot( n, scl, 0, x(1,j), 1)
        call dcopy( p, zero, 0, U(1,j), 1)
      end do

      do i = 1, n
        call daxpy( p, (-one), mu, 1, x(i,1), n)
        j = 1
        do j1 = 2, p
          call drotg( U(j,j), x(i,j), cs, sn)
          call drot( p-j, U(j,j1), p, x(i,j1), n, cs, sn)
          j = j1
        end do
        call drotg( U(p,p), x(i,p), cs, sn)
      end do

      scl = sqrt(scl)
      do j = 1, p
        call dscal( j, scl, U(1,j), 1)
      end do

      call absrng( p, U, p+1, umin, umax)

c     rcond = umin / (one + umax)

      if (umin .eq. zero) then
        hood  = FLMAX
      else
        hood  = zero
        do j = 1, p
          hood = hood + log(abs(U(j,j)))
        end do
        hood = -dn*(hood + dble(p)*(pi2log + one)/two)
      end if
c
c     do j = 1, p
c       do i = 1, j
c         x(i,j) = ddot(i,U(1,i),1,U(1,j),1)
c         if (i .ne. j) x(j,i) = x(i,j)
c        end do
c     end do
c     do j = 1, p
c       call dcopy( p, x(1,j), 1, U(1,j), 1)
c     end do

      return
      end

c Luca: add to check if compile ok

      subroutine hceee ( x, n, p, ic, ng, ns, io, jo, v, s, u, r)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

c Gaussian model-based clustering algorithm in clusters share a common
c variance (shape, volume, and orientation are the same for all clusters).

      implicit NONE

      integer            n, p, ic(n), ng, ns, io(*), jo(*)

c     double precision   x(n,p), v(p), s(p,p), u(p,p), r(p,p)
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

      integer                 q, i, j, k, l, m, i1, i2, l1, l2
      integer                 ni, nj, nij, lw, ls, lg, ici, icj
      integer                 nopt, iopt, jopt, idet, jdet, ndet

      double precision        DELOG
      double precision        ri, rj, rij, dij, tij, zij
      double precision        trc0, trc1, trcw, det0, det1, detw
      double precision        si, sj, siop, sjop, sidt, sjdt
      double precision        dopt, zopt, dijo, tijo, tdet

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      external                ddot, detmc2
      double precision        ddot, detmc2

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    EPSMIN
      parameter          (EPSMIN = 1.1102230246251565d-16)

c------------------------------------------------------------------------------

      i1 = 0
      i2 = 0

      lw = p*p

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

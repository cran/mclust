c     mclust.f for R: changed passing of characters from R to fortran
c     to passing of integers: emeee, emeev, emvev, emvvv, shapeo, unchol,
c     denvvv, deneee
      subroutine wardsw( i, n, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      subroutine mclrup( l, n, v, r, lr)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      double precision function detmc2( n, u)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      double precision function det2mc( n, u, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      subroutine hc1e  ( x, n, ic, ng, ns, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer             n, ic(n), ng, ns, nd

      double precision    x(*), d(*)

c------------------------------------------------------------------------------

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

      subroutine hc1v  ( x, n, ic, ng, ns, ALPHA, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c     integer             n, ic(n), ng, ns, nd
      integer             n, ic(*), ng, ns, nd

      double precision    x(*), ALPHA, d(*)

c------------------------------------------------------------------------------

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

      subroutine hceee ( x, n, p, ic, ng, ns, io, jo, v, s, u, r)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, ic(n), ng, ns, io(*), jo(*)

      double precision   x(n,*), v(*), s(*), u(*), r(*)

c------------------------------------------------------------------------------

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

      subroutine hceii ( x, n, p, ic, ng, ns, v, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer             n, p, ic(n), ng, ns, nd

      double precision    x(n,*), v(*), d(*)

*------------------------------------------------------------------------------

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

      subroutine hcvii ( x, n, p, ic, ng, ns, ALPHA, v, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer             n, p, ic(n), ng, ns, nd

      double precision    x(n,*), v(*), d(*), ALPHA

c------------------------------------------------------------------------------

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

      subroutine hcvvv ( x, n, p, ic, ng, ns, ALPHA, BETA, 
     *                   v, u, s, r, nd, d)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, ic(n), ng, ns, nd

      double precision   ALPHA, BETA

      double precision   x(n,*), v(*), u(p,*), s(p,*)
      double precision   r(p,*), d(*)

c------------------------------------------------------------------------------

      integer                 psq, pm1, pp1
      integer                 i, j, k, l, m, ij, iold
      integer                 lg, ld, ll, lo, ls
      integer                 ici, icj, ni, nj, nij
      integer                 nopt, niop, njop, iopt, jopt

      double precision        trcij, trmij, trop, tmop
      double precision        traci, tracj, termi, termj
      double precision        qi, qj, qij, si, sj, sij, ri, rj, rij
      double precision        dij, dopt, siop, sjop

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        rthalf
      parameter              (rthalf = .7071067811865476d0)

      external                ddot, vvvtij
      double precision        ddot, vvvtij

      double precision        BETA0, ALPHA0, ABLOG
      common /VVVMCL/         BETA0, ALPHA0, ABLOG
      save   /VVVMCL/            

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

c------------------------------------------------------------------------------

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

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      subroutine em1e ( EQPRO, x, n, G, Vinv, mu, sigsq, pro,
     *                  maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical            EQPRO

      integer             n, G, maxi

      double precision    tol, eps, Vinv

      double precision    x(*), z(n,  *  ), mu(*), sigsq, pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, iter, k, i

      double precision        hold, hood, err, tmin, tmax, ViLog
      double precision        const, sum, sumz, temp, muk, prok

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)

      eps   = max(eps,zero)
      tol   = max(tol,zero)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      iter  = iter + 1

      const = pi2log+log(sigsq)

      do k = 1, G
        muk  = mu(k)
c       prok = pro(k)
        do i = 1, n
          temp   = x(i) - muk
c         z(i,k) = prok*exp(-(const+(temp*temp)/sigsq)/two)
          z(i,k) = -(const+(temp*temp)/sigsq)/two
        end do
      end do
      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      sumz  = zero

      sigsq = zero

      do k = 1, G
        muk = zero
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          muk  = muk + temp*x(i)
        end do
        sumz   = sumz + sum
        muk    = muk / sum
        mu(k)  = muk
        if (.not. EQPRO) pro(k) = sum / dble(n)
        do i = 1, n
          temp   = x(i) - muk
          temp   = temp*temp
          sigsq  = sigsq + z(i,k)*temp
        end do
      end do

      if (Vinv .le. zero) then
        sigsq  = sigsq / dble(n)
      else
        sigsq  = sigsq / sumz
      end if

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine em1v ( EQPRO, x, n, G, Vinv, mu, sigsq, pro,
     *                  maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical              EQPRO

      integer              n, G, maxi

      double precision     Vinv, tol, eps

      double precision     x(*), z(n,  *  )

      double precision     mu(*), sigsq(*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, iter, k, i

      double precision        hold, hood, err, const, sum
      double precision        temp, sigmin, muk, sigsqk
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      err    = FLMAX

      do k = 1, G
        if (sigsq(k) .le. eps) then
          tol = FLMAX
          eps = FLMAX
          maxi = 0
        end if
      end do

      iter   = 0

100   continue

      iter   = iter + 1

      do k = 1, G
        muk    = mu(k)
c       prok   = pro(k)
        sigsqk = sigsq(k)
        const  = pi2log + log(sigsqk)
        do i = 1, n
          temp   = x(i) - muk
c         z(i,k) = prok*exp(-(const+((temp*temp)/sigsqk))/two)           
          z(i,k) = -(const+((temp*temp)/sigsqk))/two
        end do
      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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
        hood  = hood + (log(sum) + tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      sigmin = FLMAX

      do k = 1, G
        sum = zero
        muk = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          muk  = muk + temp*x(i)
        end do
        muk    = muk / sum
        mu(k)  = muk
        sigsqk = zero
        do i = 1, n
          temp   = x(i) - muk
          temp   = temp*temp
          sigsqk = sigsqk + z(i,k)*temp
          z(i,k) = temp
        end do
        sigsqk   = sigsqk / sum
        sigmin   = min(sigmin,sigsqk)
        sigsq(k) = sigsqk
        if (.not. EQPRO) pro(k)   = sum / dble(n)
      end do

      if (sigmin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emeee ( CHOL, EQPRO, x, n, p, G, Vinv, 
     *                   mu, U, pro, maxi, tol, eps, w, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

c      character          CHOL
      integer            CHOL

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, tol, eps

      double precision   x(n,*), z(n,  *  ), w(*)

      double precision   mu(p,*), U(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, info, p1, iter, i, j, k, j1

      double precision        piterm, sclfac, sum, sumz
      double precision        temp, cs, sn, umin, umax, rc, rteps
      double precision        const, hold, hood, err, detlog
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

c     FLMAX = d1mach(2)

c     if (CHOL .eq. 'N') then
      if (CHOL .eq. 0) then


c Cholesky factorization
        call dpotrf( 'U', p, U, p, info)

        if (info .ne. 0) then
c         w(1) = FLMAX
          w(1) = dble(info)
          tol  = FLMAX
          eps  = FLMAX
          return
        end if

      end if

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

      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

c condition number

      call drnge( p, U, p1, umin, umax)

      rc = umin/(one+umax)

      if (rc .le. rteps) then
c       w(1) = rc
        w(1) = zero
        eps  = FLMAX
        tol  = err
        maxi = iter
        return
      end if

      iter = iter + 1

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

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy (G, temp, 0, pro, 1)
        end if

      end if

      sumz = zero

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
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

      if (Vinv .le. zero) then
        do j = 1, p
          call dscal( j, sclfac, U(1,j), 1)
        end do
      else 
        do j = 1, p
          call dscal( j, one/sqrt(sumz), U(1,j), 1)
        end do
      end if

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      call drnge( p, U, p1, umin, umax)

      rc = umin/(one+umax)

c     w(1) = rc
      w(1) = zero

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emeei ( EQPRO, x, n, p, G, Vinv, 
     *                   mu, scale, shape, pro, maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    tol, eps, Vinv, scale

      double precision    x(n,*), z(n,*), mu(p,*), shape(*), pro(*)

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, temp, sumz, const
      double precision    hold, hood, err, smin, smax
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz    = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      if (scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      const = dble(p)*(pi2log + log(scale))

      iter = iter + 1

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + (temp*temp)/shape(j)
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/scale))/two)
          z(i,k) = -(const+(sum/scale))/two
        end do
      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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
	hood = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call dcopy( p, zero, 0, shape, 1)

      sumz  = zero
      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum/dble(n)
        do j = 1, p
          sum = zero
          do i = 1, n
            temp = x(i,j) - mu(j,k)
            sum  = sum + z(i,k)*temp*temp
          end do
          shape(j) = shape(j) + sum
        end do
      end do

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp  = exp(sum/dble(p))

      if (Vinv .le. zero) then
        scale  = temp / dble(n)
      else
        scale  = temp / sumz
      end if

      if (temp .le. eps) then
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end
      subroutine emeev ( SIGMA, EQPRO, x, n, p, G, Vinv, mu,
     *                   scale, shape, O, pro, maxi, tol, eps,
     *                   lwork, w, z, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            SIGMA
c      character          SIGMA

      logical            EQPRO

      integer            n, p, G, maxi, lwork

      double precision	 Vinv, eps, tol, scale

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, p1, iter, i, j, k, l, j1, info

      double precision        dnp, dummy, temp, rteps
      double precision        sumz, sum, smin, smax, cs, sn
      double precision        const, rc, hood, hold, err
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

c     FLMAX  = d1mach(2)

      p1 = p + 1

      if (SIGMA .ne. 0) then

c Cholesky factorization and singular value decomposition

        call dcopy( p, zero, 0, shape, 1)

        temp = zero

        l = 0
        do k = 1, G
          if (SIGMA .ne. 1) then
            call dpotrf( 'U', p, O(1,1,k), p, info)
            if (info .ne. 0) then
              temp = dble(info)
              goto 10
            end if
          end if
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, s, 
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
            goto 10
          end if
          do j = 1, p
            temp     = s(j)
            shape(j) = shape(j) + temp*temp
          end do
        end do

   10   continue

        if (temp .ne. zero .or. l .ne. 0) then
          lwork = l
c         w(1)  = FLMAX
c         w(2)  = temp
          tol   = FLMAX
          eps   = FLMAX
          err   = FLMAX
          return
        end if

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          lwork = 0
c         w(1)  = smin
c         w(2)  = zero
          tol   = err
          eps   = FLMAX
          maxi  = iter
          return
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp  = exp(sum/dble(p))

        if (Vinv .le. zero) then
          scale = temp/dble(n)
        else
          scale = temp/sumz
        end if

        if (temp .le. eps) then
          lwork = 0
c         w(1)  = temp
c         w(2)  = zero
          tol   = err
          eps   = FLMAX
          maxi  = iter
          return
        end if

        call dscal( p, one/temp, shape, 1)

      end if

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)

      dnp    = dble(n*p)

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      temp = sqrt(scale)
      do j = 1, p
        w(j) = temp*sqrt(shape(j))
      end do

      call drnge( p, w, 1, smin, smax)
      
      rc = smin / (one+smax)

      if (rc .le. rteps) then
        lwork = 0
c       w(1)  = rc
c       w(2)  = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      iter = iter + 1

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

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      call dcopy( p, zero, 0, shape, 1)

      sumz = zero

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
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
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
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          l = info
        else 
          do j = 1, p
            temp     = s(j)
            shape(j) = shape(j) + temp*temp
          end do
        end if
      end do

 110  continue

      if (l .ne. 0) then
        lwork = l        
c       w(1)  = FLMAX
c       w(2)  = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        lwork = 0
c       w(1)  = smin
c       w(2)  = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp  = exp(sum/dble(p))

      if (Vinv .le. zero) then
        scale = temp/dble(n)
      else
        scale = temp/sumz
      end if

      if (temp .le. eps) then
        lwork = 0
c       w(1)  = temp
c       w(2)  = zero
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      lwork = 0

c     w(1)  = rc
c     w(2)  = zero

      tol   = err
      eps   = hood
      maxi  = iter

      return
      end

      subroutine emeii ( EQPRO, x, n, p, G, Vinv, mu, sigsq, pro,
     *                   maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    tol, eps, Vinv, sigsq

      double precision    x(n,*), z(n,*), mu(p,*), pro(*)

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    dnp, sum, temp, sumz, tmin, tmax
      double precision    const, hold, hood, err, prok, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      dnp    = dble(n*p)

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      if (sigsq .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      iter = iter + 1

      const  = dble(p)*(pi2log+log(sigsq))

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
c         z(i,k) = pro(k)*exp(-(const+(sum/sigsq))/two)
          z(i,k) = -(const+(sum/sigsq))/two
        end do
      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
        end do
      end do

      if (Vinv .le. zero) then
        sigsq  = sigsq / dnp
      else
        sigsq  = sigsq / (dble(p)*sumz)
      end if

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emevi ( EQPRO, x, n, p, G, Vinv, 
     *                   mu, scale, shape, pro, maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    scale, tol, eps, Vinv

      double precision    x(n,*), z(n,*)

      double precision    mu(p,*), shape(p,*), pro(*)

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, const, epsmin
      double precision    hold, hood, err, smin, smax
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz    = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      if (scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call drnge( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = err
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      iter = iter + 1

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

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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
        hood = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      sumz = zero
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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      epsmin = FLMAX
      scale  = zero
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .ne. zero) then
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp   = exp(sum/dble(p))
          scale  = scale + temp
          epsmin = min(temp,epsmin)
          if (temp .gt. eps)
     *      call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      if (Vinv .gt. zero) then
        scale = scale / sumz
        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else 
        scale = scale / dble(n)
      end if

      if (epsmin .le. eps) then
        scale = zero
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emvei ( EQPRO, x, n, p, G, Vinv, mu, scale, shape, pro,
     *                   maxi, tol, eps, z, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi(2)

      double precision    Vinv, eps, tol(2)

      double precision    x(n,*), z(n,  *  ), scl(*), shp(*), w(p,*)

      double precision    mu(p,*), scale(*), shape(*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, i, j, k
      integer             iter, maxi1, maxi2, inner, inmax

      double precision    tol1, tol2, sum, temp, tmax, tmin
      double precision    prok, scalek, smin, smax, const
      double precision    hold, hood, err, errin, dnp, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

      if (maxi1 .le. 0) return

      if (Vinv .gt. zero) then
        nz    = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      eps    = max(eps,zero)

      tol1   = max(tol(1),zero)
      tol2   = max(tol(2),zero)

      dnp    = dble(n*p)

c     FLMAX  = d1mach(2)

      hold   = FLMAX/two
      hood   = FLMAX

      err    = FLMAX
      errin  = FLMAX

      inmax  = 0

      iter   = 0

100   continue

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        eps     = hood
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        maxi(1) = iter
        maxi(2) = inmax
        tol(1)  = err
        tol(2)  = errin
        eps     = hood
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

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
        do j = 1, p
          sum = zero
          do i = 1, n
            temp = x(i,j) - mu(j,k)
            sum  = sum + z(i,k)*(temp*temp)
          end do
          w(j,k) = sum
        end do
      end do

      call dscal( G, dble(p), pro, 1)

c inner iteration to estimate scale and shape
c prob now contains n*prob

      iter = iter + 1

      inner = 0

      if (maxi2 .le. 0) goto 120

110   continue

      call drnge(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
        if (Vinv .gt. zero) then
          sum = zero
          do i = 1, n
            sum = sum + z(i,nz)
          end do
          pro(nz) = sum / dble(n)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else 
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
        eps     = smin
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

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
        if (Vinv .gt. zero) then
          sum = zero
          do i = 1, n
            sum = sum + z(i,nz)
          end do
          pro(nz) = sum / dble(n)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else 
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
        eps     = smin
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

      call drnge(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
        if (Vinv .gt. zero) then
          sum = zero
          do i = 1, n
            sum = sum + z(i,nz)
          end do
          pro(nz) = sum / dble(n)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else 
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
        eps     = smin
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
      temp  = exp(sum/dble(p))

      if (temp .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
        if (Vinv .gt. zero) then
          sum = zero
          do i = 1, n
            sum = sum + z(i,nz)
          end do
          pro(nz) = sum / dble(n)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else 
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
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

      inmax = max(inner, inmax)

      if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else 
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      if (err  .gt. tol1 .and. iter .lt. maxi1) goto 100

      maxi(1) = iter
      maxi(2) = inmax

      tol(1)  = err
      tol(2)  = errin

      eps     = hood

      return
      end

      subroutine emvev ( SIGMA, EQPRO, x, n, p, G, Vinv, mu, 
     *                   scale, shape, O, pro, maxi, tol, eps, 
     *                   lwork, w, z, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

c      character          SIGMA
      integer            SIGMA

      logical            EQPRO

      integer            n, p, G, maxi(2), lwork

      double precision   Vinv, eps, tol(2)

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(p, *)

      double precision   mu(p,*), pro(  *  )

      double precision   scale(*), shape(*), O(p,p,*)

c------------------------------------------------------------------------------

      integer                 maxi1, maxi2, nz, p1, inmax, iter
      integer                 i, j, k, l, j1, info, inner

      double precision        tol1, tol2, dnp, rteps, ViLog
      double precision        rcmin, errin, smin, smax, tmin, tmax
      double precision        cs, sn, dummy, hold, hood, err
      double precision        const, temp, sum, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------
     
      maxi1  = maxi(1)
      maxi2  = maxi(2)

c     FLMAX  = d1mach(2)

      if (maxi1 .le. 0) return

      if (SIGMA .ne. 0) then

c Cholesky factorization and singular value decomposition

        temp = zero

        l = 0
        do k = 1, G
          
          if (SIGMA .ne. 1) then
            call dpotrf( 'U', p, O(1,1,k), p, info)
            if (info .ne. 0) then
              l = info
              goto 10
            end if
          end if
          call dgesvd( 'N', 'O', p, p, O(1,1,k), p, s,
     *                  dummy, 1, dummy, 1, w, lwork, info)
          if (info .ne. 0) then
            l = info
            goto 10
          end if

          do j = 1, p
            temp     = s(j,1)
            shape(j) = shape(j) + temp*temp
          end do

          call drnge( p, s, 1, smin, smax)
          if (smin .eq. zero) then
            scale(k) = zero
          else 
            sum = zero
            do j = 1, p
              sum = sum + log(s(j,1))
            end do
            scale(k) = exp(sum)
          end if

        end do

   10   continue

        if (l .ne. 0) then
          lwork   = l
          maxi(1) = -1
          maxi(2) = -1
          tol(1)  = FLMAX
          tol(2)  = FLMAX
          eps     = FLMAX
          return
        end if

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          lwork   = 0
c         w(1)    = smin
          tol(1)  = FLMAX
          tol(2)  = FLMAX
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if

        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp  = exp(sum/dble(p))

        if (temp .le. eps) then
          lwork   = 0
c         w(1)    = temp
c         w(2)    = zero
          tol(1)  = FLMAX
          tol(2)  = FLMAX
          eps     = FLMAX
          maxi(1) = -1
          maxi(2) = -1
          return
        end if

        call dscal( p, one/temp, shape, 1)

      end if

      if (Vinv .gt. zero) then
        nz    = G + 1
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

      hold   = FLMAX/two
      hood   = FLMAX

      err    = FLMAX
      errin  = FLMAX

      inmax  = 0

      iter   = 0

100   continue

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        lwork   = 0
c       w(1)    = smin
c       w(2)    = zero
        maxi(1) = iter
        maxi(2) = inmax
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        return
      end if

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        lwork   = 0
        w(1)    = -smin
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      do j = 1, p
        s(j,1) = sqrt(shape(j))
      end do

      call drnge( p, s, 1, smin, smax)

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
        const = dble(p)*(pi2log + log(scalek))
        do i = 1, n
          call dcopy( p, x(i,1), n, w(p1), 1)
          call daxpy( p, (-one), mu(1,k), 1, w(p1), 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, w(p1), 1, 
     *                 zero, w, 1)
          do j = 1, p
            w(j) = w(j) / s(j,1)
          end do
          sum    = ddot(p,w,1,w,1)/scalek
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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
        hood = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n)
      end do
      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      iter = iter + 1

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
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
        call dgesvd( 'N', 'O', p, p, O(1,1,k), p, s(1,k),
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          l = info
        else
          do j = 1, p
            temp   = s(j,k)
            s(j,k) = temp*temp
          end do
        end if
      end do

      if (l .ne. 0) then
        if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
        if (Vinv .gt. zero) then
          sum = zero
          do i = 1, n
            sum = sum + z(i,nz)
          end do
          pro(nz) = sum / dble(n)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else 
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
        lwork   = l
c       w(1)    = FLMAX
c       w(2)    = zero
        maxi(1) = iter
        maxi(2) = inner
        tol(1)  = err
        tol(2)  = errin
        eps     = FLMAX
        return
      end if

c inner iteration to estimate scale and shape
c prob now contains n*prob

      inner = 0

      if (maxi2 .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + s(j,k)/w(j)
          end do
          temp     = (sum/pro(k))/dble(p)
          scale(k) = temp
          if (temp .le. eps) then
            if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
            if (Vinv .gt. zero) then
              sum = zero
              do i = 1, n
                sum = sum + z(i,nz)
              end do
              pro(nz) = sum / dble(n)
              if (EQPRO) then
                temp = (one - pro(nz))/dble(G)
                call dcopy( G, temp, 0, pro, 1)
              end if
            else 
              if (EQPRO) call dscal( G, one/dble(G), pro, 1)
            end if
            lwork   = 0
c           w(1)    = temp
c           w(2)    = zero
            maxi(1) = iter
            maxi(2) = max(inner, inmax)
            tol(1)  = err
            tol(2)  = errin
            eps     = FLMAX
            return
          end if
          do j = 1, p
            shape(j) = shape(j) + s(j,k)/temp
          end do
        end do

        inner  = inner + 1

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
          if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
          if (Vinv .gt. zero) then
            sum = zero
            do i = 1, n
              sum = sum + z(i,nz)
            end do
            pro(nz) = sum / dble(n)
            if (EQPRO) then
              temp = (one - pro(nz))/dble(G)
              call dcopy( G, temp, 0, pro, 1)
            end if
          else 
            if (EQPRO) call dscal( G, one/dble(G), pro, 1)
          end if
          lwork   = 0
c         w(1)    = smin
c         w(2)    = zero
          maxi(1) = iter
          maxi(2) = max(inner,inmax)
          tol(1)  = err
          tol(2)  = errin
          eps     = hood
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do
        temp = exp(sum/dble(p))

        if (temp .le. eps) then
          if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)
          if (Vinv .gt. zero) then
            sum = zero
            do i = 1, n
              sum = sum + z(i,nz)
            end do
            pro(nz) = sum / dble(n)
            if (EQPRO) then
              temp = (one - pro(nz))/dble(G)
              call dcopy( G, temp, 0, pro, 1)
            end if
          else 
            if (EQPRO) call dscal( G, one/dble(G), pro, 1)
          end if
          lwork   = 0
c         w(1)    = temp
c         w(2)    = zero
          maxi(1) = iter
          maxi(2) = max(inner,inmax)
          tol(1)  = err
          tol(2)  = errin
          eps     = FLMAX
        end if

        call dscal( p, one/temp, shape, 1)
        smin = smin/temp
        smax = smax/temp

        errin = zero
        do j = 1, p
          errin = max(abs(w(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(scale(k)-w(p+k))/(one+scale(k)), errin)
        end do

        if (errin .ge. tol2 .and. inner .lt. maxi2) goto 110

120   continue

      inmax = max(inner, inmax)

      if (.not. EQPRO) call dscal( G, one/dble(n), pro, 1)

      if (Vinv .gt. zero) then
        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)
        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if
      else 
        if (EQPRO) call dscal( G, one/dble(G), pro, 1)
      end if

      if (err  .gt. tol1 .and. iter .lt. maxi1) goto 100

      lwork = 0

      smin  = sqrt(smin)
      smax  = sqrt(smax)

      rcmin = FLMAX
      do k = 1, G
        temp = sqrt(scale(k))
        rcmin = min(rcmin,(temp*smin)/(one+temp*smax))
      end do

      lwork   = 0
     
c     w(1)    = rcmin
c     w(2)    = zero

      maxi(1) = iter
      maxi(2) = inmax

      tol(1)  = err
      tol(2)  = errin

      eps     = hood

      return
      end

      subroutine emvii ( EQPRO, x, n, p, G, Vinv, mu, sigsq, pro, 
     *                   maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, tol, eps

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), sigsq(*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sumz, sum, temp, const, err, tmin, tmax
      double precision    sigmin, sigsqk, prok, hold, hood, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX
      sigmin = FLMAX

      iter   = 0

100   continue

      if (sigmin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      iter   = iter + 1

      do k = 1, G
c       prok   = pro(k)
        sigsqk = sigsq(k)
        const  = dble(p)*(pi2log+log(sigsqk))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+(sum/sigsqk))/two)           
          z(i,k) = -(const+(sum/sigsqk))/two
        end do
      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      sigmin = FLMAX

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
        end do
        sigsqk   = (sigsqk/sumz)/dble(p)
        sigsq(k) = sigsqk
        sigmin   = min(sigsqk,sigmin)
        if (.not. EQPRO) pro(k) = sumz / dble(n)
      end do

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emvvi ( EQPRO, x, n, p, G, Vinv, 
     *                   mu, scale, shape, pro, maxi, tol, eps, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    tol, eps, Vinv

      double precision    x(n,*), z(n,*)

      double precision    mu(p,*), scale(*), shape(p,*), pro(*)

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, temp, const, epsmin
      double precision    hold, hood, err, smin, smax, scalek
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps    = max(eps,zero)
      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call drnge( p, shape(1,k), 1, smin, smax)

        if (smin .le. eps) then
          tol  = err
          eps  = FLMAX
          maxi = iter
          return
        end if

      end do

      iter = iter + 1

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

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      epsmin = FLMAX
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .eq. zero) then
          scale(k) = zero
        else
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp     = exp(sum/dble(p))
          scale(k) = temp/pro(k)
          epsmin   = min(temp,epsmin)
          if (temp .gt. eps)
     *      call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      if (.not. EQPRO) then
        call dscal( G, one/dble(n), pro, 1)
      else if (Vinv .le. zero) then
        call dscal( G, one/dble(G), pro, 1)
      end if

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

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

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine emvvv ( CHOL, EQPRO, x, n, p, G, Vinv, mu, U, pro,
     *                   maxi, tol, eps, w, z) 

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

c      character          CHOL
      integer            CHOL

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

      double precision   x(n,*), z(n,*), w(*)

      double precision   mu(p,*), U(p,p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 nz, info, p1, iter, i, j, k, j1

      double precision        piterm, hold, rcmin
      double precision        temp, cs, sn, umin, umax, rteps
      double precision        sum, detlog, const, hood, err
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

c     FLMAX = d1mach(2)

      if (CHOL .eq. 0) then

c Cholesky factorization

        temp  = zero
        do k = 1, G
          call dpotrf( 'U', p, U(1,1,k), p, info)
          if (info .ne. 0) temp = dble(info)
        end do

        if (temp .ne. zero) then
c         w(1)  = FLMAX
          w(1)  = temp
          tol   = FLMAX
          eps   = FLMAX
          maxi  = -1
          return
        end if

      end if

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      if (EQPRO) call dcopy( nz, one/dble(nz), 0, pro, 1)

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

      hold   = FLMAX/two
      hood   = FLMAX

      err    = FLMAX

      iter   = 0

100   continue

      rcmin  = FLMAX

      do k = 1, G
        call drnge( p, U(1,1,k), p1, umin, umax)
        rcmin = min(umin/(one+umax),rcmin)
      end do

      if (rcmin .le. rteps) then
c       w(1) = rcmin
        w(1) = zero
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      iter  = iter + 1

      do k = 1, G

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(U(j,j,k)))
        end do

        const = piterm+detlog

c       temp = pro(k)
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U(1,1,k), p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do

      end do

      if (Vinv .gt. zero) call dcopy( n, ViLog, 0, z(1,nz), 1)

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

      if (Vinv .gt. zero) then

        sum = zero
        do i = 1, n
          sum = sum + z(i,nz)
        end do
        pro(nz) = sum / dble(n)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)

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

        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j,k), 1)
        end do

      end do

      if (err  .gt. tol .and. iter .lt. maxi) goto 100

c     w(1) = rcmin
      w(1) = zero

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine me1e ( EQPRO, x, n, G, Vinv, z, maxi, tol, eps, 
     *                  mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, G, maxi

      double precision    Vinv, eps, tol

      double precision    x(*), z(n,  *  ), mu(*), sigsq, pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, iter, k, i

      double precision        hold, hood, err, prok, tmin, tmax, ViLog
      double precision        const, sum, sumz, smu, temp, term

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX 
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter  + 1

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
        smu    = smu / sum
        mu(k)  = smu
        if (.not. EQPRO) pro(k) = sum / dble(n)
        do i = 1, n
          temp   = x(i) - smu
          temp   = temp*temp
          sigsq  = sigsq + z(i,k)*temp
          z(i,k) = temp
        end do
      end do

      if (Vinv .le. zero) then
        sigsq  = sigsq / dble(n)
      else
        sigsq  = sigsq / sumz
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
        eps  = hood
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

      return
      end

      subroutine me1v ( EQPRO, x, n, G, Vinv, z, maxi, tol, eps,
     *                  mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical              EQPRO

      integer              n, G, maxi

      double precision     Vinv, eps, tol

      double precision     x(*), z(n,  *  ), mu(*), sigsq(*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, iter, k, i

      double precision        hold, hood, err, sum, smu
      double precision        const, temp, term, sigmin, sigsqk
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      iter   = iter + 1

      do k = 1, G
        sum = zero
        smu = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
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
        if (.not. EQPRO) pro(k)   = sum / dble(n)
      end do

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

      return
      end

      subroutine meeee ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps,
     *                   mu, U, pro, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p), pro(G)
      double precision   mu(p,*), U(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 nz, p1, iter, i, j, k, j1

      double precision        piterm, sclfac, sumz, sum
      double precision        cs, sn, umin, umax, rc, detlog, rteps
      double precision        const, hold, hood, err, temp, term
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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

      sclfac = one/sqrt(dble(n))

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

      iter   = 0

100   continue

      iter = iter + 1

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
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

      term = zero
      if (Vinv .le. zero) then
        do j = 1, p
          call dscal( j, sclfac, U(1,j), 1)
        end do
      else
        do j = 1, p
          call dscal( j, one/sqrt(sumz), U(1,j), 1)
        end do
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

      call drnge( p, U, p1, umin, umax)

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

      return
      end

      subroutine meeei ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, scale

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), shape(*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term
      double precision    const, hold, hood, err, smin, smax
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

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

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      call dcopy( p, zero, 0, shape, 1)

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum/dble(n)
      end do

      do j = 1, p
        sum = zero
        do i = 1, n
          do k = 1, G
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + (temp*temp)
          end do
        end do
        shape(j) = shape(j) + sum
      end do

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      smin  = exp(sum/dble(p))

      term = zero
      if (Vinv .le. zero) then
        scale = smin/dble(n)
      else 
        scale = smin/sumz

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

      if (smin .le. eps) then
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      call dscal( p, one/smin, shape, 1)

      call drnge(p, shape, 1, smin, smax)

      if (min(scale,smin) .le. eps) then
        tol   = err
        eps   = FLMAX
        maxi  = iter
        return
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

      return
      end

      subroutine meeev ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps,
     *                   lwork, mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi, lwork

      double precision	 Vinv, eps, tol, scale

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, p1, iter, i, j, k, l, j1, info

      double precision        dnp, dummy, temp, term, rteps
      double precision        sumz, sum, smin, smax, cs, sn
      double precision        const, rc, hood, hold, err
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
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
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          l = info
        else 
          do j = 1, p
            temp     = s(j)
            shape(j) = shape(j) + temp*temp
          end do
        end if
      end do

      if (l .ne. 0) then
        lwork = l        
c       w(1)  = FLMAX
        tol   = err
        eps   = FLMAX
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

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
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
      temp  = exp(sum/dble(p))

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

      call drnge( p, shape, 1, smin, smax)
      
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

      call drnge( p, w, 1, smin, smax)

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

      lwork = 0

c     w(1)  = rc

      tol   = err
      eps   = hood
      maxi  = iter

      return
      end

      subroutine meeii ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, sigsq

      double precision    x(n,*), z(n,  *  ), mu(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, prok, tmax, tmin
      double precision    const, hold, hood, err, dnp, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

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

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      sigsq = zero

      sumz = zero

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum/dble(n)
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

      term = zero
      if (Vinv .le. zero) then
        sigsq  = sigsq / dnp
      else 
        sigsq = sigsq / (dble(p)*sumz)

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

      return
      end

      subroutine meevi ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol, scale

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), shape(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, sumz, temp, term, epsmin
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

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

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

      sumz = zero
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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum /dble(n)
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      scale  = zero
      epsmin = FLMAX
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .ne. zero) then
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp   = exp(sum/dble(p))
          scale  = scale + temp
          epsmin = min(temp,epsmin)
          if (temp .gt. eps)
     *      call dscal( p, one/temp, shape(1,k), 1)
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

      if (epsmin .le. eps) then
        tol   = err
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      if (scale .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call drnge( p, shape(1,k), 1, smin, smax)

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

      return
      end

      subroutine mevei ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi(2)

      double precision    Vinv, eps, tol(2)

      double precision    x(n,*), z(n,  *  ), scl(*), shp(*), w(p,*)

      double precision    mu(p,*), scale(*), shape(*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, i, j, k
      integer             iter, maxi1, maxi2, inner, inmax

      double precision    tol1, tol2, sum, temp, term, tmin, tmax
      double precision    prok, scalek, smin, smax, const
      double precision    hold, hood, err, errin, dnp, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

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

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp   = z(i,k)
          sum    = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
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
      end do

      call dscal( G, dble(p), pro, 1)

      if (maxi2 .le. 0) goto 120

110   continue

      call drnge(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
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
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
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

      call drnge(G, scale, 1, smin, smax)

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
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
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
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

      call drnge(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
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
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
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
      temp  = exp(sum/dble(p))

      if (temp .le. eps) then
        if (.not. EQPRO) call dscal( G, one/dnp, pro, 1)
        if (Vinv .gt. zero) then
          term = zero
          do i = 1, n
            term = term + z(i,nz)
          end do
          temp    = term / dble(n)
          pro(nz) = temp
          term    = temp * Vinv
          call dcopy( n, term, 0, z(1,nz), 1)
          if (EQPRO) then
            temp = (one - pro(nz))/dble(G)
            call dcopy( G, temp, 0, pro, 1)
          end if
        else
          if (EQPRO) call dscal( G, one/dble(G), pro, 1)
        end if
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

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps     = FLMAX
        tol(1)  = err
        tol(2)  = errin
        maxi(1) = iter
        maxi(2) = inmax
        return
      end if

      call drnge( p, shape, 1, smin, smax)

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

      if (err .gt. tol1 .and. iter .lt. maxi1) goto 100

      tol(1)  = err
      tol(2)  = errin
      eps     = hood
      maxi(1) = iter
      maxi(2) = inmax

      return
      end

      subroutine mevev ( EQPRO, x, n, p, G, Vinv, z, 
     *                   maxi, tol, eps, lwork,
     *                   mu, scale, shape, O, pro, w, s)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi(2), lwork

      double precision   Vinv, eps, tol(2)

      double precision   x(n,*), z(n,  *  ), w(  *  ), s(*)

      double precision   mu(p,*), pro(  *  )

      double precision   scale(*), shape(*), O(p,p,*)

c------------------------------------------------------------------------------

      integer                 maxi1, maxi2, p1, inmax, iter
      integer                 nz, i, j, k, l, j1, info, inner

      double precision        tol1, tol2, dnp, term, rteps, ViLog
      double precision        errin, smin, smax, sumz, tmin, tmax
      double precision        cs, sn, dummy, hold, hood, err
      double precision        const, temp, sum, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
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
     *                dummy, 1, dummy, 1, w, lwork, info)
        if (info .ne. 0) then
          l = info
        else 
          do j = 1, p
            temp     = z(j,k)
            z(j,k)   = temp*temp
          end do
        end if
      end do

      iter = iter + 1

      if (l .ne. 0) then
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
        eps     = FLMAX
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

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
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
        temp  = exp(sum/dble(p))

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
          temp     = (sum/pro(k))/dble(p)
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

        call drnge( p, shape, 1, smin, smax)

        if (smin .eq. zero) then
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
        temp = exp(sum/dble(p))

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

      call drnge( p, shape, 1, smin, smax)
      
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

      call drnge( G, scale, 1, smin, smax)

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

      call drnge( p, s, 1, smin, smax)
      
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

      subroutine mevii ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), sigsq(*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sumz, sum, temp, const, term
      double precision    sigmin, sigsqk, hold, hood, err
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

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

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

      iter   = 0

100   continue

      iter   = iter + 1

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
        sigsq(k) = (sigsqk/sumz)/dble(p)
        if (.not. EQPRO) pro(k) = sumz / dble(n)
      end do

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
            if (temp .ge. SMLOG) then
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

      subroutine mevvi ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical             EQPRO

      integer             n, p, G, maxi

      double precision    Vinv, eps, tol

      double precision    x(n,*), z(n,  *  )

      double precision    mu(p,*), scale(*), shape(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer             nz, iter, i, j, k

      double precision    sum, temp, term, scalek, epsmin
      double precision    hold, hood, err, smin, smax, const
      double precision    prok, tmin, tmax, ViLog

      double precision    zero, one, two
      parameter          (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision    pi2log
      parameter          (pi2log = 1.837877066409345d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157d308)                  

      double precision    SMLOG
      parameter          (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
      end if

      eps   = max(eps,zero)
      tol   = max(tol,zero)

c     FLMAX = d1mach(2)
      hold  = FLMAX/two
      hood  = FLMAX
      err   = FLMAX

      iter  = 0

100   continue

      iter  = iter + 1

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      epsmin = FLMAX
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .eq. zero) then
          scale(k) = zero
        else
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do
          temp     = exp(sum/dble(p))
          scale(k) = temp/pro(k)
          epsmin   = min(temp,epsmin)
          if (temp .gt. eps)
     *      call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do

      if (.not. EQPRO) then
        call dscal( G, one/dble(n), pro, 1)
      else if (Vinv .le. zero) then
        call dscal( G, one/dble(G), pro, 1)
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
        eps   = -FLMAX
        maxi  = iter
        return
      end if

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

        call drnge( p, shape(1,k), 1, smin, smax)

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

      return
      end

      subroutine mevvv ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, U, pro, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

      double precision   x(n,*), z(n,*), w(*)

      double precision   mu(p,*), U(p,p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 nz, p1, iter, i, j, k, j1

      double precision        piterm, hold, rcmin, rteps
      double precision        temp, term, cs, sn, umin, umax
      double precision        sum, detlog, const, hood, err
      double precision        prok, tmin, tmax, ViLog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

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
        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, U(j,i,k), 1)
          i = j
        end do
      end do

      iter   = 0

100   continue

      iter = iter + 1

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j,k), 1)
        end do
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        if (.not. EQPRO) pro(k) = sum / dble(n)
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
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j,k), 1)
        end do
      end do

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
        call drnge( p, U(1,1,k), p1, umin, umax)
        rcmin = min(umin/(one+umax),rcmin)
      end do

      if (rcmin .le. rteps) then
c       w(1) = rcmin
        tol  = err
        eps  = FLMAX
        maxi = iter
        return
      end if

      do k = 1, G

c       temp = pro(k)

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(U(j,j,k)))
        end do

        const = piterm+detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, U(1,1,k), p, w, 1)
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

      if (err .gt. tol .and. iter .lt. maxi) goto 100

c     w(1) = rcmin

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

      subroutine es1e ( x, mu, sigsq, pro, n, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, G

      double precision   sigsq, hood, Vinv

      double precision   x(*), mu(*), pro(  *  ), z(n,  *  )

c------------------------------------------------------------------------------

      integer                 i, k, nz

      double precision        temp, const, muk, prok, tmin, tmax, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (sigsq .le. max(hood,zero)) then
c       FLMAX = d1mach(2)
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
          z(i,k) = -(const+(temp*temp)/sigsq)/two
        end do
      end do

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

      return
      end

      subroutine es1v  ( x, mu, sigsq, pro, n, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, G

      double precision   hood, Vinv

c     double precision   x(n), mu(G), sigsq(G), pro(G[+1]), z(n,G[+1])
      double precision   x(*), mu(*), sigsq(*), pro(  *  ), z(n,  *  )

c------------------------------------------------------------------------------

      integer                 i, k, nz

      double precision        temp, const, tmin, tmax, sum
      double precision        muk, sigsqk, prok, sigmin

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      sigmin = FLMAX

      do k = 1, G
        sigmin = min(sigmin,sigsq(k))
      end do

      if (sigmin .le. max(hood,zero)) then
        hood  = FLMAX
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
          z(i,k) = -(const+(temp*temp)/sigsqk)/two
        end do
      end do

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

      return
      end

      subroutine eseee ( CHOL, x, mu, Sigma, pro, n, p, G, Vinv,
     *                   w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      character          CHOL

c     integer            n, p, G
      integer            n, p, G

      double precision   hood, Vinv

c     double precision   x(n,p), w(p), z(n,G[+1])
      double precision   x(n,*), w(*), z(n,  *  )

c     double precision   mu(p,G), Sigma(p,p), pro(G[+1])
      double precision   mu(p,*), Sigma(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 info, i, j, k, nz

      double precision        eps, rteps, detlog, prok, tmin, tmax
      double precision        umin, umax, const, temp, rc, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)
      
      eps   = max(hood,zero)
      rteps = sqrt(eps)

      if (CHOL .eq. 'N') then

c Cholesky factorization
        call dpotrf( 'U', p, Sigma, p, info)

        w(1) = dble(info)

        if (info .ne. 0) then
          hood  = FLMAX
          return
        end if

      end if

      call drnge( p, Sigma, (p+1), umin, umax)

      rc   = umin/(one+umax)

      if (rc .le. rteps) then
        w(1)  = zero
c       FLMAX = d1mach(2)
        hood  = FLMAX
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

      w(1) = zero

      return
      end

      subroutine eseei ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   scale, hood, Vinv

c     double precision   x(n,p), z(n,G[+1])
      double precision   x(n,*), z(n,  *  )

c     double precision   mu(p,G), shape(p), pro(G[+1])
      double precision   mu(p,*), shape(*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        eps, rteps, smin, smax, prok

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      eps   = max(hood, zero)
      rteps = sqrt(eps)

      if (scale .le. eps) then
c       FLMAX = d1mach(2)
        hood  = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. rteps) then
c       FLMAX = d1mach(2)
        hood  = FLMAX
        return
      end if

      const = dble(p)*(pi2log+log(scale)) 

      do k = 1, G
c       prok = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

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

      return
      end

      subroutine eseev ( x, mu, scale, shape, O, pro, n, p, G, 
     *                   Vinv, v, w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

c     integer            n, p, G
      integer            n, p, G

      double precision   scale, Vinv, hood

      double precision   x(n,*), v(*), w(*), z(n,  *  )

      double precision   mu(p,*), shape(*), O(p,p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        const, temp, rteps, tmin, tmax
      double precision        smin, smax, prok, eps, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      eps    = max(hood, zero)
      rteps  = sqrt(eps)
  
      if (scale .le. eps) then
        hood = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do

      call drnge( p, shape, 1, smin, smax)
  
      if (smin .le. rteps) then
        hood = FLMAX
        return
      end if
        
      const = dble(p)*(pi2log + log(scale))

      do k = 1, G

c       prok = pro(k)

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, 
     *                 w, 1, zero, v, 1)
          do j = 1, p
            v(j) = v(j)/shape(j)
          end do
          temp   = ddot( p, v, 1, v, 1)
c         z(i,k) = prok*exp(-(const+temp)/two)
          z(i,k) = -(const+temp)/two
        end do

      end do

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

      return
      end

      subroutine eseii ( x, mu, sigsq, pro, n, p, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   sigsq, hood, Vinv

c     double precision   x(n,p), mu(p,G), pro(G[+1]), z(n,G[+1])
      double precision   x(n,*), mu(p,*), pro(  *  ), z(n,  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, prok, tmin, tmax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

      if (sigsq .le. max(hood,zero)) then
c       FLMAX = d1mach(2)
        hood  = FLMAX
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
          z(i,k) = -(const+sum/sigsq)/two
        end do
      end do

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

      return
      end

      subroutine esevi ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   scale, hood, Vinv

      double precision   x(n,*), z(n,  *  )

      double precision   mu(p,*), shape(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, eps, tmin, tmax
      double precision        smin, smax, prok, rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      eps   = max(hood, zero)
      rteps = sqrt(eps)
    
      if (scale .le. eps) then
        hood  = FLMAX
        return
      end if

      temp = sqrt(scale)
      do k = 1, G
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
        call drnge( p, shape(1,k), 1, smin, smax)
        if (smin .le. rteps) then
          hood  = FLMAX
          return
        end if
      end do

      const  = dble(p)*(pi2log+log(scale))

      do k = 1, G
c       prok   = pro(k)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j,k)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

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

      return
      end

      subroutine estep2( n, ncolz, pro, z, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, ncolz

      double precision   pro(  *  ), z(n,  *  ), hood

c------------------------------------------------------------------------------

      integer                 i, k

      double precision        temp, tmin, tmax, sum, prok

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157D+308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     hood = zero
c     do i = 1, n
c       temp = zero
c       do k = 1, ncolz
c         temp = temp + z(i,k)
c       end do
c       hood = hood + log(temp)
c       call dscal( ncolz, (one/temp), z(i,1), n)
c     end do

      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX
        do k = 1, ncolz
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
        do k = 1, ncolz
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
        call dscal( ncolz, (one/sum), z(i,1), n)
      end do

      return
      end

      subroutine esvei ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

      double precision   x(n,*), z(n,  *  )

      double precision   mu(p,*), scale(*), shape(*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, eps, tmin, tmax
      double precision        smin, smax, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      eps   = max(hood, zero)

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        hood  = FLMAX
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        hood  = FLMAX
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
            temp = (x(i,j) - mu(j,k))/shape(j)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum/scalek)/two)
          z(i,k) = -(const+sum/scalek)/two
        end do
      end do

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

      return
      end

      subroutine esvev ( x, mu, scale, shape, O, pro, n, p, G, 
     *                   Vinv, v, w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

c     integer            n, p, G
      integer            n, p, G

      double precision   Vinv, hood

      double precision   x(n,*), z(n,  *  ), mu(p,*), pro(  *  )

      double precision   v(*), w(*)

      double precision   scale(*), shape(*), O(p,p,*)

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        const, temp, eps, tmin, tmax
      double precision        smin, smax, scalek, prok, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      eps    = max(hood, zero)

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        hood = FLMAX
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
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
            v(j) = v(j)/shape(j)
          end do
          temp   = ddot( p, v, 1, v, 1)/scalek
c         z(i,k) = prok*exp(-(const+temp)/two)
          z(i,k) = -(const+temp)/two
        end do

      end do

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

      smin = FLMAX
      smax = zero
      do k = 1, G
        scalek = sqrt(scale(k))
        do j = 1, p
          temp = scalek*shape(j)
          smin = min(temp, smin)
          smax = min(temp, smax)
        end do
      end do

      w(1) = smin / (one+smax)

      return
      end

      subroutine esvii ( x, mu, sigsq, pro, n, p, G, Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

      double precision   x(n,*), z(n,  *  )

      double precision   mu(p,*), sigsq(*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmin, tmax
      double precision        prok, sigsqk, sigmin

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      sigmin = FLMAX

      do k = 1, G
        sigmin = min(sigmin,sigsq(k))
      end do

      if (sigmin .le. max(hood,zero)) then
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
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum/sigsqk)/two)
          z(i,k) = -(const+sum/sigsqk)/two
        end do
      end do

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

      return
      end 

      subroutine esvvi ( x, mu, scale, shape, pro, n, p, G, 
     *                   Vinv, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   hood, Vinv

      double precision   x(n,*), z(n,  *  )

      double precision   mu(p,*), scale(*), shape(p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 i, j, k, nz

      double precision        sum, temp, const, eps, tmin, tmax
      double precision        smin, smax, prok, scalek, rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      eps   = max(hood, zero)
      rteps = sqrt(eps)

      do k = 1, G
        temp = sqrt(scale(k))
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
      end do

      do k = 1, G
       
        call drnge( p, shape(1,k), 1, smin, smax)

        if (smin .le. rteps) then
          hood  = FLMAX
          return
        end if

      end do

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        hood  = FLMAX
        return
      end if

      do k = 1, G
c       prok   = pro(k)
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j,k)
            sum  = sum + temp*temp
          end do
c         z(i,k) = prok*exp(-(const+sum)/two)
          z(i,k) = -(const+sum)/two
        end do
      end do

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

      return
      end

      subroutine esvvv ( CHOL, x, mu, Sigma, pro, n, p, G, Vinv, 
     *                   w, hood, z)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      character          CHOL

      integer            n, p, G

      double precision   hood, Vinv

      double precision   x(n,*), w(*), z(n,  *  )

      double precision   mu(p,*), Sigma(p,p,*), pro(  *  )

c------------------------------------------------------------------------------

      integer                 nz, p1, info, i, j, k

      double precision        const, detlog, temp, prok, tmin, tmax
      double precision        umin, umax, rcmin, eps, rteps, sum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMLOG
      parameter              (SMLOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      p1    = p + 1

c     FLMAX = d1mach(2)

      eps   = max(hood,zero)
      rteps = sqrt(eps)

      if (CHOL .eq. 'N') then

        do k = 1, G

          call dpotrf( 'U', p, Sigma(1,1,k), p, info)

          w(1) = dble(info)

          if (info .ne. 0) then
            hood = FLMAX
            return
          end if
       
        end do

      end if

      rcmin = FLMAX
      do k = 1, G
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rcmin = min(rcmin,umin/(one+umax))
      end do

      if (rcmin .le. rteps) then
        w(1) = zero
        hood = FLMAX
        return
      end if

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

      w(1) = zero

      return
      end

      subroutine ms1e ( x, z, n, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer             n, G

      double precision    x(*), z(n,*), mu(*), sigsq, pro(*)

c------------------------------------------------------------------------------

      integer                 i, k

      double precision        sum, smu, sumz, temp

      double precision        zero
      parameter              (zero = 0.d0)

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
	smu    = smu  / sum
	mu(k)  = smu 
        pro(k) = sum / dble(n)
        do i = 1, n
          temp   = x(i) - smu
          temp   = temp*temp
          sigsq  = sigsq + z(i,k)*temp
        end do
      end do

c sumz .eq. n when no noise
      sigsq  = sigsq / sumz

      return
      end

      subroutine ms1v ( x, z, n, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, G

      double precision   x(*), z(n,*), mu(*), sigsq(*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, k
     
      double precision        sum, smu, temp, sigsqk

      double precision        zero
      parameter              (zero = 0.d0)

c------------------------------------------------------------------------------

      do k = 1, G
        sum = zero
        smu = zero      
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          smu  = smu + temp*x(i)
        end do
        smu    = smu / sum
        mu(k)  = smu 
        sigsqk = zero
        do i = 1, n
          temp   = x(i) - smu
          temp   = temp*temp
          sigsqk = sigsqk + z(i,k)*temp
        end do
        sigsq(k) = sigsqk / sum
        pro(k)   = sum / dble(n)
      end do

      return
      end

      subroutine mseee ( x, z, n, p, G, w, mu, U, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   x(n,*), z(n,*), w(*)

      double precision   mu(p,*), U(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k, j1

      double precision        sum, sumz, temp, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

      do j = 1, p
        call dcopy( p, zero, 0, U(1,j), 1)
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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum / dble(n)
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

c sumz .eq. n when no noise
      do j = 1, p
        call dscal( j, one/sqrt(sumz), U(1,j), 1)
      end do

      return
      end

      subroutine mseei ( x, z, n, p, G, eps, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   eps

      double precision   x(n,*), z(n,*)

      double precision   mu(p,*), scale, shape(*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, sumz, temp, smin, smax

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum/dble(n)
      end do

      call dcopy( p, zero, 0, shape, 1)

      do j = 1, p
        sum = zero
        do i = 1, n
          do k = 1, G
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
        end do
        shape(j) = shape(j) + sum
      end do

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        eps   = smin
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do      
      temp  = exp(sum/dble(p))

      scale = temp/sumz

      if (temp .le. eps) then
        eps = temp
        return
      end if

      call dscal( p, one/temp, shape, 1)

      eps = temp
      
      return
      end

      subroutine mseev ( x, z, n, p, G, w, lwork, eps, 
     *                   mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G, lwork

      double precision   eps, scale

      double precision   x(n,*), z(n,*), w(  *  )

      double precision   shape(*), O(p,p,*), mu(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k, j1, l, info
      double precision        dummy, sum, sumz, temp, const
      double precision        cs, sn, smin, smax

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

      call dcopy( p, zero, 0, shape, 1)

      l = 0

      sumz = zero

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum / dble(n)
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
            l = info
          else 
            do j = 1, p
              temp     = z(j,k)
              shape(j) = shape(j) + temp*temp
            end do
          end if
        end if
      end do

      lwork = l

      if (lwork .ne. 0) return

      call drnge( p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        eps   = smin
        scale = zero
        return
      end if
        
      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      const = exp(sum/dble(p))

c sumz .eq. n when no noise
      scale = const/sumz

      if (const .le. eps) then
        eps = const
        return
      end if

      call dscal( p, one/const, shape, 1)

      eps = const

      return
      end

      subroutine mseii ( x, z, n, p, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   x(n,*), z(n,*), mu(p,*), sigsq, pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, sumz, temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum/dble(n)
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          sigsq  = sigsq + z(i,k)*sum
        end do
      end do

c sumz .eq. n when no noise
      sigsq  = sigsq / (sumz*dble(p))

      return
      end

      subroutine msevi ( x, z, n, p, G, eps, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   eps, scale

      double precision   x(n,*), z(n,*)

      double precision   mu(p,*), shape(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        smin, smax
      double precision        sum, sumz, temp, epsmin

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      sumz = zero
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
        sumz = sumz + sum
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum/dble(n)
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      epsmin = FLMAX
      scale  = zero
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .ne. zero) then
          sum   = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do      
          temp   = exp(sum/dble(p))
          scale  = scale + temp
          epsmin = min(temp,epsmin)
          if (temp .gt. eps) 
     *      call dscal( p, one/temp, shape(1,k), 1)
        end if
      end do
      
      if (epsmin .eq. zero) then
        scale = zero
      else
        scale = scale/sumz
      end if

      eps = epsmin
      
      return
      end

      subroutine msvei ( x, z, n, p, G, maxi, tol, eps,
     *                   mu, scale, shape, pro, scl, shp, w)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G, maxi

      double precision   eps, tol

      double precision   x(n,*), z(n,*), scl(*), shp(*), w(p,*)

      double precision   mu(p,*), scale(*), shape(*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k, iter

      double precision        sum, temp, smin, smax, err

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps   = max(eps,zero)
      tol   = max(tol,zero)

c     FLMAX = d1mach(2)
      err   = FLMAX

c start with the equal volume and shape estimate

      do k = 1, G
        sum = zero
        call dcopy( p, zero, 0, mu(1,k), 1)
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
        do j = 1, p
          sum = zero
          do i = 1, n
            temp     = x(i,j) - mu(j,k)
            temp     = temp*temp
            temp     = z(i,k)*temp
            sum      = sum + temp
          end do
          w(j,k)   = sum
        end do
      end do

      iter = 0

      do j = 1, p
        sum = zero       
        do k = 1, G
          sum = sum + w(j,k)
        end do
        shape(j) = sum
      end do

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        call dscal( G, one/dble(n), pro, 1)
        eps   = smin
        tol  = err
        maxi = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do      
      temp  = exp(sum/dble(p))

      call dcopy( G, temp/dble(n), 0, scale, 1)

      if (temp .le. eps) then
        call dscal( G, one/dble(n), pro, 1)
        eps  = temp
        tol  = err
        maxi = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      call dscal( G, dble(p), pro, 1)

100   continue

      call drnge(p, shape, 1, smin, smax)

      if (smin .le. eps) then
        call dscal( G, one/dble(n*p), pro, 1)
        eps  = smin
        tol  = err
        maxi = iter
        return
      end if

      iter = iter + 1

c scale estimate

      call dcopy( G, scale, 1, scl, 1)

      do k = 1, G
        sum = zero       
        do j = 1, p
          sum = sum + w(j,k)/shape(j)
        end do
        scale(k) = sum/pro(k)
      end do

      call drnge(G, scale, 1, smin, smax)

      if (smin .le. eps) then
        call dscal( G, one/dble(n*p), pro, 1)
        eps  = smin
        tol  = err
        maxi = iter
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

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        call dscal( G, one/dble(n*p), pro, 1)
        eps   = smin
        tol  = err
        maxi = iter
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do      
      temp  = exp(sum/dble(p))

      if (temp .le. eps) then
        call dscal( G, one/dble(n*p), pro, 1)
        eps  = temp
        tol  = err
        maxi = iter
        return
      end if

      call dscal( p, one/temp, shape, 1)

      err = zero
      
      do k = 1, G
        err = max(err, abs(scl(k) - scale(k))/(one + scale(k)))        
      end do
      
      do j = 1, p
        err = max(err, abs(shp(j) - shape(j))/(one + shape(j)))       
      end do
 
      if (err .gt. tol .and. iter .le. maxi) goto 100

      call dscal( G, one/dble(n*p), pro, 1)

      call drnge(G, scale, 1, smin, smax)

      eps  = smin

      call drnge(p, shape, 1, smin, smax)
    
      eps  = min(eps,smin)
      tol  = err
      maxi = iter

      return
      end

      subroutine msvev ( x, z, n, p, G, w, lwork, maxi, tol, eps,
     *                   mu, scale, shape, O, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G, maxi, lwork

      double precision   eps, tol

      double precision   x(n,*), z(n,*), w(*)

      double precision   scale(*), shape(*), O(p,p,*), mu(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 p1, i, j, k, j1, inner, info

      double precision        dnp, epsmin, errin, dummy
      double precision        temp, sum, smin, smax, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX =  1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps    = max(eps,zero)
      tol    = max(tol,zero)

      p1     = p + 1

      dnp    = dble(n*p)

c     FLMAX  = d1mach(2)

      epsmin = FLMAX

      errin  = FLMAX

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
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
            lwork = 0
          else
            do j = 1, p
              temp     = z(j,k)
              temp     = temp*temp
              shape(j) = shape(j) + temp
              z(j,k)   = temp
            end do
          end if
        end if
      end do
  
      lwork = info
 
      if (info .ne. 0) then
        call dscal( G, one/dble(n), pro, 1)
        maxi = -1
        eps  = -one
        tol  = -one
      end if

      call drnge( p, shape, 1, smin, smax)
  
      if (smin .eq. zero) then
        call dscal( G, one/dble(n), pro, 1)
        maxi = -1
        tol  = -one
        eps  = smin
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp = exp(sum/dble(p))
 
      call dcopy( G, temp/dble(n), 0, scale, 1)
     
      if (temp .le. eps) then
        call dscal( G, one/dble(n), pro, 1)
        maxi = -1
        tol  = -one
        eps  = temp
        return
      end if 

      epsmin = min(epsmin,temp)

      call dscal( p, one/temp, shape, 1)

c iteration to estimate scale and shape
c pro now contains n*pro

      inner  = 0

      if (maxi .le. 0) goto 120

110   continue

        call dcopy( p, shape, 1, w    , 1)
        call dcopy( G, scale, 1, w(p1), 1)

        call dcopy( p, zero, 0, shape, 1)

        do k = 1, G
          sum = zero
          do j = 1, p
            sum = sum + z(j,k)/w(j)
          end do
          temp     = (sum/pro(k))/dble(p)
          scale(k) = temp
          if (temp .le. eps) then
            maxi = inner
            tol  = errin
            eps  = temp
            return
          end if
          epsmin = min(epsmin,temp)
          do j = 1, p
            shape(j) = shape(j) + z(j,k)/temp
          end do
        end do

        inner  = inner + 1

        call drnge( p, shape, 1, smin, smax)
 
        if (smin .eq. zero) then
          call dscal( G, one/dble(n), pro, 1)
          maxi = inner
          tol  = errin
          eps  = smin
          return
        end if

c normalize the shape matrix
        sum = zero
        do j = 1, p
          sum = sum + log(shape(j))
        end do

        temp = exp(sum/dble(p))

        if (temp .le. eps) then
          call dscal( G, one/dble(n), pro, 1)
          maxi = inner
          tol  = errin
          eps  = temp
          return
        end if

        epsmin = min(epsmin,temp)

        call dscal( p, one/temp, shape, 1)

        errin = zero
        do j = 1, p
          errin = max(abs(w(j)-shape(j))/(one+shape(j)), errin)
        end do

        do k = 1, G
          errin = max(abs(scale(k)-w(p+k))/(one+scale(k)), errin)
        end do
        
        if (errin .ge. tol .and. inner .lt. maxi) goto 110

120   continue

      call dscal( G, one/dble(n), pro, 1)

      eps  = epsmin
      tol  = errin
      maxi = inner

      return
      end

      subroutine msvii ( x, z, n, p, G, mu, sigsq, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   x(n,*), z(n,*), mu(p,*), sigsq(*), pro(*)

c------------------------------------------------------------------------------
      
      integer                 i, j, k
     
      double precision        sum, sumz, temp, sigsqk

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

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
        end do
        sigsq(k) = (sigsqk/sumz)/dble(p)
        pro(k)   = sumz / dble(n)
      end do

      return
      end

      subroutine msvvi ( x, z, n, p, G, eps, mu, scale, shape, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   eps

      double precision   x(n,*), z(n,*)

      double precision   mu(p,*), scale(*), shape(p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, temp, smin, smax, epsmin

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

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
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum
      end do

c pro(k) now contains n_k

      do j = 1, p
        do k = 1, G
          sum = zero
          do i = 1, n
            temp = z(i,k)*(x(i,j) - mu(j,k))
            sum  = sum + temp*temp
          end do
          shape(j,k) = shape(j,k) + sum
        end do
      end do

      epsmin = FLMAX
      do k = 1, G
        call drnge(p, shape(1,k), 1, smin, smax)
        epsmin = min(smin,epsmin)
        if (smin .eq. zero) then
          scale(k) = zero
        else 
          sum = zero
          do j = 1, p
            sum = sum + log(shape(j,k))
          end do      
          temp     = exp(sum/dble(p))
          scale(k) = temp/pro(k)
          epsmin   = min(temp,epsmin)
          if (temp .gt. eps) 
     *      call dscal( p, one/temp, shape(1,k), 1)
        end if
        pro(k) = pro(k)/dble(n)
      end do
      
      eps = epsmin
      
      return
      end

      subroutine msvvv ( x, z, n, p, G, w, mu, U, pro)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   x(n,*), z(n,*), w(*)

      double precision   mu(p,*), U(p,p,*), pro(*)

c------------------------------------------------------------------------------

      integer                 i, j, k, j1

      double precision        sum, temp, cs, sn

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

      do k = 1, G
        call dcopy( p, zero, 0, mu(1,k), 1)
        do j = 1, p
          call dcopy( j, zero, 0, U(1,j,k), 1)
        end do
        sum = zero
        do i = 1, n
          temp = z(i,k)
          sum  = sum + temp
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1)
        end do
        call dscal( p, (one/sum), mu(1,k), 1)
        pro(k) = sum / dble(n)
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
        do j = 1, p
          call dscal( j, one/sqrt(sum), U(1,j,k), 1)
        end do
      end do

      return
      end

      subroutine mvn1d ( x, n, mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c     integer            n
      integer            n

      double precision   mu,  sigsq, hood

      double precision   x(*)

c------------------------------------------------------------------------------

      double precision        dn, scl

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot
      external                ddot

c------------------------------------------------------------------------------
      
      dn = dble(n)

      scl = one/dble(n)
      mu  = ddot( n, scl, 0, x, 1)

      sigsq = zero

      call daxpy( n, (-one), mu, 0, x, 1)
      sigsq = ddot( n, x, 1, x, 1)
 
      sigsq = sigsq/dn

      if (sigsq .eq. zero) then
        hood  = FLMAX
      else 
        hood  = -dn*(pi2log + (one + log(sigsq)))/two
      end if

      return
      end

      subroutine mvnxii( x, n, p, mu, sigsq, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p

      double precision   sigsq, hood

      double precision   x(n,*), mu(*)

c------------------------------------------------------------------------------

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

      subroutine mvnxxi( x, n, p, mu, scale, shape, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c     integer            n, p
      integer            n, p

      double precision   scale, hood

      double precision   x(n,*), mu(*), shape(*)

c------------------------------------------------------------------------------

      integer                 i, j

      double precision        sum, temp, smin, smax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        ddot
      external                ddot

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
          sum = sum + temp*temp
        end do
        shape(j) = shape(j) + sum
      end do

      call drnge(p, shape, 1, smin, smax)

      if (smin .eq. zero) then
        scale = zero
        hood  = FLMAX
        return
      end if

      sum = zero
      do j = 1, p
        sum = sum + log(shape(j))
      end do
      temp  = exp(sum/dble(p))

      scale = temp/dble(n)

      call dscal( p, one/temp, shape, 1)

      hood  = -dble(n*p)*(one + pi2log + log(scale))/two

      return
      end

      subroutine mvnxxx( x, n, p, mu, U, hood)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p

      double precision   hood

      double precision   x(n,*), mu(*), U(p,*)

c------------------------------------------------------------------------------

      integer                 i, j, j1

      double precision        dnp, scl
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

      dnp = dble(n*p)

      scl = one/dble(n)
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

      call drnge( p, U, p+1, umin, umax)

c     rcond = umin / (one + umax)

      if (umin .eq. zero) then
        hood  = FLMAX
      else
        hood  = zero
        do j = 1, p
          hood = hood + log(abs(U(j,j)))
        end do
        hood = -dble(n)*(hood + dble(p)*(pi2log + one)/two)
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

      subroutine drnge( l, v, i, vmin, vmax)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

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

      subroutine mcltrw( x, n, p, u, ss)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit double precision (a-h,o-z)

      integer            n, p

      double precision   x(n,*), u(*)

c------------------------------------------------------------------------------

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
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit double precision (a-h,o-z)

c     integer            n, p, iwork(liwork)
      integer            n, p, iwork(*)

      double precision   x(n,*), u(*), v(p,*), w(p,p), work(*)

c------------------------------------------------------------------------------

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

      double precision    EPSMAX
      parameter          (EPSMAX = 2.2204460492503131d-16)

      double precision    FLMAX
      parameter          (FLMAX = 1.7976931348623157D+308)

c     double precision        d1mach

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

      subroutine shapeo( TRANSP, s, O, l, m, w, info)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c      character          TRANSP
      integer            TRANSP

      integer            l, m, info

      double precision   s(*), O(l,l,*), w(l,*)

c------------------------------------------------------------------------------

      integer                 j, k

      double precision        temp

      double precision        zero, one
      parameter              (zero = 0.d0, one = 1.d0)

c------------------------------------------------------------------------------

      if (TRANSP .eq. 1) then
 
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

      if (TRANSP .eq. 0) then

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

      subroutine unchol ( UPLO, T, l, n, info)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c      character          UPLO
      integer          UPLO

      integer            l, n, info

c     double precision   T(abs(n), abs(n))
      double precision   T(  l   ,   *   )

c------------------------------------------------------------------------------

       integer            i, j, k

       external           ddot
       double precision   ddot

c------------------------------------------------------------------------------

      if (UPLO .eq. 1) then

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

      if (UPLO .eq. 0) then

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

      subroutine den1e ( x, mu, sigsq, n, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, G

      double precision   sigsq, eps

      double precision   x(*), mu(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, k

      double precision        sumz, temp, const

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps   = max(eps,zero)

      if (sigsq .le. eps) then
        eps = FLMAX
        return
      end if

      eps = sigsq

      const = pi2log+log(sigsq)

      do i = 1, n
        sumz = zero
        do k = 1, G
          temp      = x(i) - mu(k)
          temp      = temp*temp
          dens(i,k) = -(const+temp/sigsq)/two
        end do
      end do

      return
      end

      subroutine den1v ( x, mu, sigsq, n, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, G

      double precision   eps

      double precision   x(*), mu(*), sigsq(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, k

      double precision        temp, const, muk, sigsqk

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps  = max(eps,zero)

      temp = sigsq(1)

      if (G .gt. 1) then
        do k = 2, G
          temp = min(temp,sigsq(k))
        end do
      end if

      if (temp .le. eps) then
        eps = FLMAX
        return
      end if

      eps = temp

      do k = 1, G
        sigsqk = sigsq(k)
        const  = pi2log+log(sigsqk)
        muk    = mu(k)
        do i = 1, n
          temp      = x(i) - muk
          dens(i,k) = -(const+(temp*temp)/sigsqk)/two
        end do
      end do

      return
      end

      subroutine deneee( CHOL, x, mu, Sigma, n, p, G, w, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            CHOL
c      character          CHOL

c     integer            n, p, G
      integer            n, p, G

      double precision   eps

      double precision   x(n,*), mu(p,*), Sigma(p,*)

      double precision   w(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 info, i, j, k

      double precision        detlog, const, temp
      double precision        umin, umax, rc

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      if (CHOL .eq. 0) then

c Cholesky factorization
        call dpotrf( 'U', p, Sigma, p, info)

        w(1) = dble(info)

        if (info .ne. 0) then
          w(1)  = FLMAX
          return
        end if

      end if

      call drnge( p, Sigma, (p+1), umin, umax)

      eps = max(eps,zero)
      eps = sqrt(eps)

      rc  = umin/(one+umax)

      if (rc .le. eps) then
        eps = rc
        return
      end if

      eps = rc

      detlog = zero
      do j = 1, p
        detlog = detlog + log(abs(Sigma(j,j)))
      end do

      const = dble(p)*pi2log/two + detlog

      do k = 1, G
        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma, p, w, 1)
          temp      = ddot( p, w, 1, w, 1)/two
          dens(i,k) = -(const+temp)
        end do
      end do

      w(1) = zero

      return
      end

      subroutine deneei ( x, mu, scale, shape, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   scale, eps

      double precision   x(n,*), mu(p,*), shape(*)

      double precision   dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, temp, const, smin, smax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps = max(eps,zero)

      if (scale .le. eps) then
        eps = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. sqrt(eps)) then
c       FLMAX = d1mach(2)
        eps   = FLMAX
        return
      end if

      eps = smin*smin

      const = dble(p)*(pi2log+log(scale)) 

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum)/two
        end do
      end do

      return
      end

      subroutine deneev( x, mu, scale, shape, O, n, p, G, 
     *                   v, w, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   scale, eps

      double precision   x(n,*), mu(p,*), shape(*), O(p,p,*)

      double precision   v(*), w(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        const, temp, smin, smax

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      eps    = max(eps, zero)
  
      if (scale .le. eps) then
        eps = FLMAX
        return
      end if

      temp = sqrt(scale)
      do j = 1, p
        shape(j) = temp*sqrt(shape(j))
      end do

      call drnge( p, shape, 1, smin, smax)
  
      if (smin .le. sqrt(eps)) then
        eps = FLMAX
        return
      end if

      eps = smin*smin
        
      const = dble(p)*(pi2log + log(scale))/two

      do k = 1, G

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, 
     *                 w, 1, zero, v, 1)
          do j = 1, p
            v(j) = v(j)/shape(j)
          end do
          temp      = ddot( p, v, 1, v, 1)/two
          dens(i,k) = -(const+temp)
        end do

      end do

      return
      end

      subroutine deneii( x, mu, sigsq, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   sigsq, eps

      double precision   x(n,*), mu(p,*), dens(n, *)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, sumz, temp, const

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps   = max(eps,zero)

      if (sigsq .le. eps) then
        eps = FLMAX
        return
      end if

      eps = sigsq

      const = dble(p)*(pi2log+log(sigsq))

      do i = 1, n
        sumz = zero
        do k = 1, G
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum/sigsq)/two
        end do
      end do

      return
      end

      subroutine denevi( x, mu, scale, shape, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   scale, eps

      double precision   x(n,*), mu(p,*), shape(p,*)

      double precision   dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, temp, const
      double precision        smin, smax, rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      eps   = max(eps, zero)
      rteps = sqrt(eps)
    
      if (scale .le. eps) then
        eps = FLMAX
        return
      end if

      temp = sqrt(scale)
      do k = 1, G
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
        call drnge( p, shape(1,k), 1, smin, smax)
        if (smin .le. rteps) then
          eps = FLMAX
          return
        end if
      end do

      const  = dble(p)*(pi2log+log(scale))

      do k = 1, G
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j,k)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum)/two
        end do
      end do

      return
      end

      subroutine denvei( x, mu, scale, shape, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   eps

      double precision   x(n,*), mu(p,*), scale(*), shape(*)

      double precision   dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, temp, const, smin, smax, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)
 
      eps   = max(eps,zero)

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        return
      end if

      eps = smin

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        return
      end if

      do j = 1, p
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum/scalek)/two
        end do
      end do

      return
      end

      subroutine denvev( x, mu, scale, shape, O, n, p, G, 
     *                   v, w, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

c     integer            n, p, G
      integer            n, p, G

      double precision   eps

      double precision   x(n,*), mu(p,*)

      double precision   scale(*), shape(*), O(p,p,*)

      double precision   v(*), w(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        const, temp, smin, smax, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

c     FLMAX  = d1mach(2)

      eps    = max(eps, zero)

      call drnge( p, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        return
      end if

      call drnge( p, shape, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        return
      end if

      eps = smin

      do j = 1, p
        shape(j) = sqrt(shape(j))
      end do

      do k = 1, G

        scalek = scale(k)
        
        const = dble(p)*(pi2log+log(scalek))

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dgemv( 'N', p, p, one, O(1,1,k), p, 
     *                 w, 1, zero, v, 1)
          do j = 1, p
            v(j) = v(j)/shape(j)
          end do
          temp      = ddot( p, v, 1, v, 1)/scalek
          dens(i,k) = -(const+temp)/two
        end do

      end do

      return
      end

      subroutine denvii( x, mu, sigsq, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   sigsq(*), eps

      double precision   x(n,*), mu(p,*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, sigsqk, temp, const

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

      eps  = max(eps,zero)

      temp = sigsq(1)

      if (G .gt. 1) then
        do k = 2, G
          temp = min(temp,sigsq(k))
        end do
      end if

      if (temp .le. eps) then
        eps = FLMAX
        return
      end if

      eps = temp

      do k = 1, G
        sigsqk = sigsq(k)
        const  = dble(p)*(pi2log+log(sigsq(k)))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = x(i,j) - mu(j,k)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum/sigsqk)/two
        end do
      end do

      return
      end
      subroutine denvvi( x, mu, scale, shape, n, p, G, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      integer            n, p, G

      double precision   eps

      double precision   x(n,*), mu(p,*), scale(*), shape(p,*)

      double precision   dens(n,*)

c------------------------------------------------------------------------------

      integer                 i, j, k

      double precision        sum, temp, const
      double precision        smin, smax, scalek, rteps

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

c------------------------------------------------------------------------------

c     FLMAX = d1mach(2)

      eps   = max(eps, zero)
      rteps = sqrt(eps)

      do k = 1, G
        temp = sqrt(scale(k))
        do j = 1, p
          shape(j,k) = temp*sqrt(shape(j,k))
        end do
      end do

      do k = 1, G
       
        call drnge( p, shape(1,k), 1, smin, smax)

        if (smin .le. rteps) then
          eps = FLMAX
          return
        end if

      end do

      call drnge( G, scale, 1, smin, smax)

      if (smin .le. eps) then
        eps = FLMAX
        return
      end if

      eps = smin

      do k = 1, G
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum = zero
          do j = 1, p
            temp = (x(i,j) - mu(j,k))/shape(j,k)
            sum  = sum + temp*temp
          end do
          dens(i,k) = -(const+sum)/two
        end do
      end do

      return
      end

      subroutine denvvv( CHOL, x, mu, Sigma, n, p, G, w, eps, dens)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.

      implicit NONE

      character          CHOL

c     integer            n, p, G
      integer            n, p, G

      double precision   eps

      double precision   x(n,*), mu(p,*), Sigma(p,p,*)

      double precision   w(*), dens(n,*)

c------------------------------------------------------------------------------

      integer                 p1, info, i, j, k

      double precision        const, detlog, temp
      double precision        umin, umax, rc

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      p1    = p + 1

c     FLMAX = d1mach(2)

      if (CHOL .eq. 'N') then

        do k = 1, G

          call dpotrf( 'U', p, Sigma(1,1,k), p, info)

          w(1) = dble(info)

          if (info .ne. 0) then
            w(1) = FLMAX
            return
          end if
       
        end do

      end if

      eps = max(eps,zero)
      eps = sqrt(eps)

      rc   = FLMAX
      do k = 1, G
        call drnge( p, Sigma(1,1,k), p1, umin, umax)
        rc  = min(rc,umin/(one+umax))
      end do

      if (rc .le. eps) then
        eps = rc
        return
      end if
 
      eps = rc

      do k = 1, G

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(Sigma(j,j,k)))
        end do

        const = dble(p)*pi2log/two + detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, Sigma(1,1,k), p, w, 1)
          temp      = ddot( p, w, 1, w, 1)/two
          dens(i,k) = -(const+temp)
        end do

      end do

      w(1) = zero

      return
      end
c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt
c Distribution of MCLUST is prohibited except by agreement with the 
c University of Washington.
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

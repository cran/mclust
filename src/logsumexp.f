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

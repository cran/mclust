*  =====================================================================
      subroutine transpose(X, p)
*
*  Compute transpose of a matrix 
*
*  =====================================================================
      implicit NONE
      integer :: p, i, j
      double precision :: X(p,p), temp
      
      do j = 2, p
        do i = 1, j-1
            temp = X(i,j)
            X(i,j) = X(j,i)
            X(j,i) = temp
        end do
      end do
      
      return
      end

*  =====================================================================
      subroutine crossprodf(X, Y, n, p, q, XTY)
*
*  Given matrices X and Y of dimension (n x p) and (n x q) computes 
*  the matrix of cross-product, i.e. X' Y
*
*  =====================================================================
      implicit NONE
      integer             n, p, q
      double precision    X(n,p), Y(n,q), XTY(p,q)
*     Compute X'Y using DGEMM blas subroutine
      call DGEMM('T', 'N', p, q, n, 1.d0, X, n, Y, n, 0.d0, XTY, p)
      end

* ======================================================================
      subroutine covwf ( X, Z, n, p, G, M, S, W )
*
*   Given data matrix X(n x p) and weight matrix Z(n x G) computes
*   weighted means M(p x G), weighted covariance matrices S(p x p x G)
*   and weighted scattering matrices W(p x p x G)
*
* ======================================================================

      implicit none

      integer :: n, p, G
      double precision :: X(n,p), Z(n,G)

      double precision :: M(p,G), S(p,p,G), W(p,p,G)

      integer :: j, k
      double precision :: sumZ(G), temp(n,p)

*     compute X'Z using BLAS
      call dgemm('T', 'N', p, G, n, 1.d0, X, n, Z, n, 0.d0, M, p)

*     compute row sums of Z
      sumZ = sum(Z, DIM = 1)

      do k = 1,G

*         compute means
          call dscal(p, (1.d0/sumZ(k)), M(:,k), 1)     
        
              do j = 1,p                  
*                 compute sqrt(Z) * (X - M)
                  temp(:,j) = sqrt(Z(:,k)) * (X(:,j) - M(j,k))
              end do
        
*         compute scattering matrix    
          call dgemm('T', 'N', p, p, n, 1.d0, temp, n, temp, n,
     *               0.d0, W(:,:,k), p)

*         compute covariance matrix
          S(:,:,k) = W(:,:,k)/sumZ(k)

      end do
  
      return 
      end

************************************************************************
**** EVV model
************************************************************************        
* ======================================================================
      subroutine msevv (x,z, n,p,G, mu,O,U,scale,shape,pro, lwork,info, 
     *                      eps)
*       Maximization step for model EEV
* ======================================================================

      implicit none
      
      integer :: n, p, G
      
      double precision :: x(n,p), z(n,G)
      
      double precision :: mu(p,G), O(p,p,*), U(p,p,*), pro(G)
      double precision :: scale(G), shape(p,G)
      
      double precision :: sumz(G)
      integer :: i, j, k, info, lwork, dummy, l
          
      double precision :: temp(p), wrk(lwork), eps
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
*        
*      double precision :: BIGLOG
*      parameter (BIGLOG =  709.d0)
*        
*      double precision :: SMALOG
*      parameter (SMALOG = -708.d0)
        
*-----------------------------------------------------------------------

*       colsums of z
      sumz = sum(z, dim = 1)
      
*       a priori probabilities
      pro = sumz / dble(n)
*      pro = sumz / sum(sumz)
*       if there is noise sum(sumz) does not sum to n. See help(mstep)
      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
*            U(:,:,k) = U(:,:,k) + 
*     *                  spread(temp, dim = 2, ncopies = p)*
*     *                  spread(temp, dim = 1, ncopies = p)
*       outer product, Press et al. (1992), p. 970
            call dger(p, p, 1.d0, temp, 1, temp, 1, U(:,:,k), p)
*       more efficient            
        end do
*       U contains the weighted scattering matrix          

        O(:,:,k) = U(:,:,k)
*        call dgesvd('O', 'N', p, p, O(:,:,k), p, shape(:,k),
*     *              dummy, 1, dummy, 1, wrk, lwork, info)
        call dgesvd('N', 'O', p, p, O(:,:,k), p, shape(:,k),
     *              dummy, 1, dummy, 1, wrk, lwork, info)
*       O now contains eigenvectors of the scattering matrix
*       ##### NOTE: O is transposed
*       shape contains the eigenvalues    
     
*       check if dgesvd converged (info == 0)
        if (info .ne. 0) then            
            l = info
        else
           scale(k) = exp( sum( log(shape(:,k)) ) )**(1.d0/p)
           call dscal(p*p, 1.d0/scale(k), U(:,:,k), 1)
           call dscal(p, 1.d0/scale(k), shape(:,k), 1)
*       now U is the matrix Ck (Celeux, Govaert 1995, p.787)
*       and shape is the proper scaled shape (matrix A)           
        end if      
        
      end do
      
*       check very small eigenvalues (singular covariance) 
      if (minval(shape) .le. sqrt(eps) .or. 
     *          minval(scale) .le. sqrt(eps)) then
        shape = FLMAX
        scale = FLMAX
        return
      end if
      
      scale(1) = sum(scale) / sum(sumz)
        
      return
      end
      

* ======================================================================
      subroutine esevv (x,z, n,p,G,Gnoise, mu,O,scale,shape,pro, Vinv, 
     *                 loglik, eps)
*       Expectation step for model EVV
* ======================================================================
      
      implicit none
      
      integer :: n, p, G, Gnoise
      double precision :: x(n,p), z(n,Gnoise)
      double precision :: mu(p,G), O(p,p,G), scale, shape(p,G)
      double precision :: Vinv, pro(Gnoise)
      
      double precision :: temp1(p), temp2(p), temp3
      
      integer :: i, k, j
      double precision :: const, logdet, loglik, eps
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)      
      
      external :: ddot
      double precision :: ddot

      
*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------
      
*       check very small eigenvalues (singular covariance) 
      if (minval(shape) .le. sqrt(eps) .or. scale .le. sqrt(eps)) then
        loglik = FLMAX
        return
      end if
      
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O(:,:,k), p, temp1, 1, 0.d0, temp2, 1)
            temp2 = temp2/sqrt(scale*shape(:,k))
            temp3 = ddot(p, temp2, 1, temp2, 1)
*       temp3 contains the mahalanobis distance
*            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
            z(i,k) = const - logdet/2.d0 - temp3/2.d0
*       help(cdens) --> The densities are not scaled by mixing proportions            
        end do
*       z contains the log-density log(N(x|theta_k)) 
            
      end do
      
      if ( pro(1) .lt. 0.d0 ) return
*       cdens function

*       noise component
      if (Vinv .gt. 0.d0) then
        call dcopy( n, log(Vinv), 0, z(:,Gnoise), 1)
      end if 
*       now column Gnoise of z contains log(Vinv)           
      
      do i = 1,n
      
        z(i,:) = z(i,:) + log( pro(:) )
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1) )

*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
        
      end do
      
      return
      end
      
      
* ======================================================================
      subroutine meevv (x,z, n,p,G,Gnoise, mu,O,U,scale,shape,pro,Vinv,
     *                  loglik, eqpro,itmax,tol,eps,
     *                  niterout,errout,lwork,info)
*       Maximization-expectation algorithm for model EVV
* ======================================================================
      
      implicit none
      
      logical :: eqpro
      integer :: n, p, G, Gnoise
      double precision :: x(n,p), z(n,Gnoise)
      double precision :: mu(p,G), O(p,p,G),scale(G),shape(p,G)
      double precision :: Vinv, pro(Gnoise)
      double precision :: U(p,p,G), sumz(Gnoise)
      
      double precision :: temp1(p), temp2(p), temp3, scsh(p)
*      double precision :: temp(*)
      
      integer :: i, j, k, info, lwork, dummy, l, itmax, niterout
      double precision :: tol, eps, errout, rteps
      double precision :: const, logdet, loglik, lkprev, wrk(lwork)
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot

*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------

      l = 0
      rteps = sqrt(eps)
      niterout = 0
      errout = FLMAX
      lkprev = FLMAX/2
      loglik = FLMAX
      
      const = (-dble(p)/2.d0)*log2pi
      
*       WHILE loop using goto statement
100   continue
      
      niterout = niterout + 1

      sumz = sum(z, dim = 1)
      if ( eqpro ) then
        if ( Vinv .gt. 0 ) then
            pro(Gnoise) = sumz(Gnoise) / dble(n)
            pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
            sumz = pro * dble(n)
        else 
            pro = 1 / dble(G)
            sumz = pro * dble(n)
        end if
      else
        pro = sumz / dble(n)
      end if
      
*       re-initialise U      
      call dcopy(p*p*G, 0.d0, 0, U, 1)
      
*       M step..........................................................      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1, 1, temp1, 1, U(:,:,k), p)         
        end do
*       U contains the weighted scattering matrix          

        O(:,:,k) = U(:,:,k)
        call dgesvd('N', 'O', p, p, O(:,:,k), p, shape(:,k),
     *              dummy, 1, dummy, 1, wrk, lwork, info)
*       O now contains eigenvectors of the scattering matrix
*       ##### NOTE: O is transposed
*       shape contains the eigenvalues    
     
*       check if dgesvd converged (info == 0)
        if (info .ne. 0) then            
            l = info
            return
        else
           scale(k) = exp( sum( log(shape(:,k)) ) )**(1.d0/dble(p))
           call dscal(p*p, 1.d0/scale(k), U(:,:,k), 1)
           call dscal(p, 1.d0/scale(k), shape(:,k), 1)
*       now U is the matrix Ck (Celeux, Govaert 1995, p.787)
*       and shape is the proper scaled shape (matrix A)           
        end if      
    
      end do
      
      if ( Vinv .gt. 0.d0 ) then
        scale(1) = sum(scale) / sum(sumz(1:G))
      else
        scale(1) = sum(scale)/dble(n)
      end if
*       if noise lambda = num/sum_{k=1}^{G} n_k; pag. 787 Celeux, Govaert      
*       ................................................................
            
*       check very small eigenvalues (singular covariance) 
      if (minval(shape) .le. rteps .or. minval(scale) .le. rteps) then
        loglik = FLMAX
        return
      end if
      
*       E step..........................................................
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale(1)) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed
        scsh = sqrt(scale(1)*shape(:,k))
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O(:,:,k), p, temp1, 1, 0.d0, temp2, 1)
            temp2 = temp2/scsh
            temp3 = ddot(p, temp2, 1, temp2, 1)
*       temp3 contains the mahalanobis distance
            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
        end do
*       z contains the log-density log(N(x|theta_k)) + log(p_k) 
            
      end do
      
*       noise component
      if (Vinv .gt. 0.d0) then
*        call dcopy( n, log(Vinv) + log(pro(Gnoise)), 0, z(:,Gnoise), 1)
        z(:,Gnoise) = log(Vinv) + log( pro(Gnoise) )
      end if 
*       now column Gnoise of z contains log(Vinv) + log(p_0)
*       with p_0 the proportion of noise       
      
      loglik = 0.d0
      do i = 1,n
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1) )
        
*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
      end do
*       ................................................................      
      
      errout = abs(loglik - lkprev)/(1.d0 + abs(loglik))
*      errout = abs(loglik - lkprev)
      lkprev = loglik
*      temp(niterout) = loglik

* Chris F (June 2015): pro should not be computed in the E-step
*     sumz = sum(z, dim = 1)
*     if ( eqpro ) then
*       if ( Vinv .gt. 0 ) then
*           pro(Gnoise) = sumz(Gnoise) / dble(n)
*           pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
*           sumz = pro * dble(n)
*       else 
*           pro = 1 / dble(G)
*           sumz = pro * dble(n)
*       end if
*     else
*       pro = sumz / dble(n)
*     end if
      
*       check if empty components
*      if ( minval(pro) .lt. rteps ) then
      if ( any(sumz .lt. rteps, 1) ) then
        loglik = -FLMAX
        return
      end if
      
*       WHILE condition      
      if ( errout .gt. tol .and. niterout .lt. itmax ) goto 100  
      
      return
      end
  
 
************************************************************************
**** VEE model
************************************************************************
* ======================================================================
      subroutine msvee (x,z, n,p,G, mu,U,C,scale,pro, lwork,info,
     *                  itmax,tol, niterin,errin)
*       Maximization step for model VEE
* ======================================================================

      implicit none
      
      integer :: n, p, G
      
      double precision :: x(n,p), z(n,G)
      
      double precision :: mu(p,G), U(p,p,G), C(p,p), pro(G)
*       ### NOTE: shape and orientation parameters are computed in R 
      double precision :: scale(G)
      
      double precision :: sumz(G)
      integer :: i, j, k, info, lwork, l
*     integer :: dummy          
      double precision :: temp1(p), temp2(p,p), temp3 
      double precision :: wrk(lwork), tol, errin, trgt, trgtprev
      integer :: itmax, niterin
    
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot
        
*-----------------------------------------------------------------------

*       colsums of z
      sumz = sum(z, dim = 1)
      
*       a priori probabilities
      pro = sumz / dble(n)
*      pro = sumz / sum(sumz)
*       if there is noise sum(sumz) does not sum to n. See help(mstep)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1, 1, temp1, 1, U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       check if U is positive definite (see help of dpotrf)
*       (through Choleski is more efficient)
        temp2 = U(:,:,k)
        call dpotrf('U', p, temp2, p, info)
        if ( info .ne. 0 ) then
            if ( info .lt. 0) then
                l = info
                return
            else if ( info .gt. 0 ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
        
      end do
      
*       covariance matrix components estimation      
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2
      
*       WHILE loop using goto statement
100   continue

      niterin = niterin + 1
      
*       initialise C      
      call dcopy(p*p, 0.d0, 0, C, 1)
*       ### NOTE: scale is initialised in R
      
      do k = 1,G
        C = C + U(:,:,k)/scale(k)
      end do
*       C contains the numerator of matrix C in pag.785, Celeux, Govaert

      temp2 = C
      call dsyev('N', 'U', p, temp2, p, temp1, wrk, lwork, info)      
      temp1 = temp1(p:1:-1)
*       temp1 contains the (decreasing) ordered eigenvalues of C

*       check if dsyev converged or illegal value
      if ( info .ne. 0 ) then
        l = info
        return
      end if      
      
      temp3 = exp( sum(log(temp1)) )**(1/dble(p))
*       temp3 is the denominator of C
      
      C = C/temp3
*       C is now the actual matrix C of pag.785      
      
*       compute the inverse of C via Choleski
      temp2 = C
      call dpotrf('U', p, temp2, p, info)
        if ( info .ne. 0 ) then
            if ( info .lt. 0) then
                l = info
                return
            else if ( info .gt. 0 ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
      call dpotri('U', p, temp2, p, info)
        if ( info .ne. 0 ) return
      do j = 2,p
        do k = 1,(j-1)
            temp2(j,k) = temp2(k,j)
        end do
      end do
*       temp2 is now the inverse of C      
      
      scale = 0.d0
      do k = 1,G
        do j = 1,p
            scale(k) = scale(k) + ddot(p, U(j,:,k), 1, temp2(:,j), 1) 
        end do
        scale(k) = scale(k) / (dble(p)*sumz(k))
      end do
*       scale contains now the lambdas (pag.784 of Celeux, Govaert)      
      
*       evaluate target function      
*      trgt = dble(n)*dble(p) + dble(p)*SUM(log(scale)*sumz)
      trgt = sum(sumz)*dble(p) + dble(p)*SUM(log(scale)*sumz)
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt
      
*       WHILE condition      
      if ( errin .gt. tol .and. niterin .lt. itmax ) goto 100  
        
      return
      end


* ======================================================================
      subroutine esvee (x,z, n,p,G,Gnoise, mu,O,scale,shape,pro, Vinv,
     *                  loglik, eps)
*       Expectation step for model VEE
* ======================================================================
      
      implicit none
      
      integer :: n,p,G,Gnoise
      double precision :: x(n,p), z(n,Gnoise), pro(Gnoise), Vinv
      double precision :: mu(p,G), O(p,p), scale(G), shape(p)
      
      double precision :: temp1(p), temp2(p), temp3
      
      integer :: i, k, j
      double precision :: const, logdet, loglik, eps
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)            
      
      external :: ddot
      double precision :: ddot

      
*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------
      
*       check very small eigenvalues (cannot compute E step) 
      if ( minval(shape) .le. sqrt(eps) .or. 
     *                  minval(scale) .le. sqrt(eps) ) then
        loglik = FLMAX
        return
      end if
      
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j)) + log(scale(k)) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O, p, temp1, 1, 0.d0, temp2, 1)
            temp2 = temp2/sqrt(scale(k)*shape)
            temp3 = ddot(p, temp2, 1, temp2, 1)
*       temp3 contains the mahalanobis distance
*            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
            z(i,k) = const - logdet/2.d0 - temp3/2.d0
*       help(cdens) --> The densities are not scaled by mixing proportions            
        end do
*       z contains the log-density log(N(x|theta_k)) 
            
      end do
      
      if ( pro(1) .lt. 0.d0 ) return
*       cdens function     

*       noise component
      if (Vinv .gt. 0.d0) then
        call dcopy( n, log(Vinv), 0, z(:,Gnoise), 1)
      end if 
*       now column Gnoise of z contains log(Vinv)           
      
      do i = 1,n
      
        z(i,:) = z(i,:) + log( pro(:) )
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
        
        z(i,:) = exp( z(i,:) - temp1(1) )
        
*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
        
      end do
      
      return
      end
      
      
* ======================================================================
      subroutine mevee ( x,z, n,p,G,Gnoise, mu,C,U,scale,shape,pro,Vinv,
     *                  loglik, eqpro,itmaxin,tolin,itmaxout,tolout,eps, 
     *                  niterin,errin,niterout,errout,lwork,info )
*       Maximization-expectation algorithm for model VEE
* ======================================================================
      
      implicit none
      
      logical :: eqpro
      integer :: n,p,G,Gnoise
      double precision :: x(n,p), z(n,Gnoise), pro(Gnoise), Vinv
      double precision :: mu(p,G), C(p,p), scale(G), shape(p)
      double precision :: U(p,p,G), sumz(Gnoise)
      
      double precision :: temp1(p), temp2(p,p), temp3, temp4(p)
      
      integer :: i, j, k, info, lwork, l
*			integer :: dummy
      integer :: itmaxin, itmaxout, niterin, niterout
      double precision :: tolin, tolout, errin, errout, eps, rteps
      double precision :: const, logdet, loglik, lkprev, wrk(lwork)
      double precision :: trgt, trgtprev
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot

*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------

      l = 0
      rteps = sqrt(eps)
      niterout = 0
      errout = FLMAX
      lkprev = FLMAX/2
      loglik = FLMAX
      
      const = (-dble(p)/2.d0)*log2pi
      
*       WHILE loop for EM algorithm
100   continue
    
      niterout = niterout + 1
      
      sumz = sum(z, dim = 1)
      if ( eqpro ) then
        if ( Vinv .gt. 0 ) then
            pro(Gnoise) = sumz(Gnoise) / dble(n)
            pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
            sumz = pro * dble(n)
        else 
            pro = 1 / dble(G)
            sumz = pro * dble(n)
        end if
      else
        pro = sumz / dble(n)
      end if
      
*       re-initialise U      
      call dcopy(p*p*G, 0.d0, 0, U, 1)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1, 1, temp1, 1, U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       check if U is positive definite (see help of dpotrf)
*       (through Choleski is more efficient)
        temp2 = U(:,:,k)
        call dpotrf('U', p, temp2, p, info)
        if ( info .ne. 0 ) then
            if ( info .lt. 0) then
                l = info
                return
            else if ( info .gt. 0 ) then
                info = 0
                loglik = FLMAX
                return
            end if
        end if
        
      end do
      
*       M step..........................................................       
*       covariance matrix components estimation      
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2
*       initialise scale        
      call dcopy(G, 1.d0, 0, scale, 1)
           
*       WHILE loop for M step
110   continue

      niterin = niterin + 1
      
*       initialise C      
      call dcopy(p*p, 0.d0, 0, C, 1)
      
      do k = 1,G
        C = C + U(:,:,k)/scale(k)
      end do
*       C contains the numerator of matrix C in pag.785, Celeux, Govaert

      temp2 = C
      call dsyev('N', 'U', p, temp2, p, temp1, wrk, lwork, info)      
      temp1 = temp1(p:1:-1)
*       temp1 contains the (decreasing) ordered eigenvalues of C
      
*       check if dsyev converged or illegal value
      if ( info .ne. 0 ) then
        l = info
        return
      end if      

      
      temp3 = exp( sum(log(temp1)) )**(1/dble(p))
*       temp3 is the denominator of C
      
      C = C/temp3
*       C is now the actual matrix C of pag.785      
      
*       compute the inverse of C via Choleski
      temp2 = C      
      call dpotrf('U', p, temp2, p, info)
        if ( info .ne. 0 ) then
            if ( info .lt. 0) then
                l = info
                return
            else if ( info .gt. 0 ) then
                info = 0
                loglik = FLMAX
                return
            end if
        end if
      call dpotri('U', p, temp2, p, info)
        if ( info .ne. 0 ) return
      do j = 2,p
        do k = 1,(j-1)
            temp2(j,k) = temp2(k,j)
        end do
      end do
*       temp2 is now the inverse of C      
      
      scale = 0.d0
      do k = 1,G
        do j = 1,p
            scale(k) = scale(k) + ddot(p, U(j,:,k), 1, temp2(:,j), 1) 
        end do
        scale(k) = scale(k) / (dble(p)*sumz(k))
      end do
*       scale contains now the lambdas (pag.784 of Celeux, Govaert)      
      
*       evaluate target function      
*      trgt = dble(n)*dble(p) + dble(p)*SUM(log(scale)*sumz)
       trgt = sum(sumz(1:G))*dble(p) + dble(p)*SUM(log(scale)*sumz(1:G))
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt
      
*       WHILE condition for M step     
      if ( errin .gt. tolin .and. niterin .lt. itmaxin ) goto 110 
*       ................................................................

*       eigenvalues of C      
      shape = temp1 / temp3
      
*       check very small eigenvalues (singular covariance)     
      if (minval(shape) .le. rteps .or. minval(scale) .le. rteps) then
        loglik = FLMAX
        return
      end if
      
*       E step..........................................................
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + log(shape(j)) + log(scale(k))
        end do
    
*       compute mahalanobis distance for each observation
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp4, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 temp2, p, temp1, 1, 0.d0, temp4, 1)
            temp4 = temp4/scale(k)
            temp3 = ddot(p, temp4, 1, temp1, 1)
*       temp3 contains the mahalanobis distance
            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
*            z(i,k) = const - logdet/2.d0 - temp3/2.d0          
        end do
*       z contains the log-density log(N(x|theta_k)) + log(p_k) 
      end do
      
*      if ( pro(1) .lt. 0.d0 ) return
*       cdens function   

*       noise component
      if (Vinv .gt. 0.d0) then
        z(:,Gnoise) = log(Vinv) + log( pro(Gnoise) )
      end if 
*       now column Gnoise of z contains log(Vinv) + log(p_0)
*       with p_0 the proportion of noise       
      
      loglik = 0.d0
      do i = 1,n
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
        
        z(i,:) = exp( z(i,:) - temp1(1) )
*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
      end do
*       ................................................................

      errout = abs(loglik - lkprev)/(1.d0 + abs(loglik))
*      errout = abs(loglik - lkprev)
      lkprev = loglik
*      temp(niterout) = loglik

* Chris F (June 2015): pro should not be computed in the E-step
*     sumz = sum(z, dim = 1)
*     if ( eqpro ) then
*       if ( Vinv .gt. 0 ) then
*           pro(Gnoise) = sumz(Gnoise) / dble(n)
*           pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
*           sumz = pro * dble(n)
*       else 
*           pro = 1 / dble(G)
*           sumz = pro * dble(n)
*       end if
*     else
*       pro = sumz / dble(n)
*     end if
      
*       check if empty components
      if ( minval(sumz) .lt. rteps ) then
        loglik = -FLMAX
        return
      end if
      
*       WHILE condition EM     
      if ( errout .gt. tolout .and. niterout .lt. itmaxout ) goto 100

      return
      end
      

************************************************************************
**** EVE model
************************************************************************
* ======================================================================
      subroutine mseve (x,z, n,p,G, mu,U,O,scale,shape,pro, lwork,info,
     *                  itmax,tol, niterin,errin, eps)
*       Maximization step for model EVE
* ======================================================================

      implicit none
      
      integer :: n, p, G
      
      double precision :: x(n,p), z(n,G)
      
      double precision :: mu(p,G), U(p,p,G), pro(G), O(p,p)
      double precision :: scale, shape(p,G)
      
      double precision :: sumz(G), omega(G)
      integer :: i, j, k, info, lwork
          
      double precision :: temp1(p,p), temp2(p,p), temp3(p,p), temp4(p) 
      double precision :: wrk(lwork), tol, errin, trgt, trgtprev, eps
      integer :: itmax, niterin
    
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot
        
*-----------------------------------------------------------------------

*       colsums of z
      sumz = sum(z, dim = 1)
      
*       a priori probabilities
      pro = sumz / dble(n)
*      pro = sumz / sum(sumz)
*       if there is noise sum(sumz) does not sum to n. See help(mstep)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1(:,1), 1, temp1(:,1), 1, 
     *                 U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       compute the eigenvalues of U to be stored in omega
        temp2 = U(:,:,k)
        call dsyev('N', 'U', p, temp2, p, temp1(:,1), wrk, lwork, info)
*       now temp1 contains all the eigenvalues of U

*       check if dsyev converged and positive definite
        if ( info .ne. 0 ) then
            return
        else
            if ( minval(temp1(:,1)) .lt. sqrt(eps) ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
        
        omega(k) = temp1(p,1)      
        
      end do
*       omega contains the largest eigenvalue of each scattering matrix
   
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2
      
*       covariance matrix components estimation
*       we consider algorithm MM 1 and MM 2 of Browne, McNicholas 2013
*       with a modification in computing the orientation matrix in the MM 2 step

*      shape (matrix A) and orientation (matrix D) initialised in R        
*      shape = matrix(1, p,G)
*      O = diag(p)      
      
*       WHILE loop using goto statement
100   continue

*       ### NOTE: O is transposed

      niterin = niterin + 1
      temp2 = 0.d0
      temp3 = 0.d0
*       temp3 will contain matrix F      
      
*       Algorithm MM 1 ......................................
      do k = 1,G
        
        do j = 1,p
*            temp1(j,:) = O(:,j) / shape(j,k)
            temp1(j,:) = O(j,:) / shape(j,k)
        end do
*       temp1 contains inv(A)t(D)
        
        call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, U(:,:,k),p, 
     *              0.d0, temp2,p )
*       temp2 contains inv(A) %*% t(D) %*% W    
        
        temp1 = temp2 - omega(k)*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do
      
*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return
      
*       NOTE: t(O) = t( R %*% t(P) ) = P %*% t(R)      
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       Algorithm MM 2 ......................................
*      call dgemm( 'T','T', p,p,p, 1.d0, temp2,p, temp1,p, 
*     *              0.d0, O,p )
      call transpose(O, p)
*       O contains matrix D of Browne, McNicholas
  
*       Algorithm MM 2
      temp1 = 0.d0
      temp3 = 0.d0
      do k = 1,G
        
      call dgemm( 'N','N', p,p,p, 1.d0, U(:,:,k),p, O,p, 
     *            0.d0, temp1,p )
*       temp1 contains W %*% D
        
        do j = 1,p
            temp2(:,j) = temp1(:,j) / shape(j,k)
        end do
*       temp2 contains W %*% D %*% inv(A)
        
        temp1 = temp2 - maxval( 1/shape(:,k) )*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do

*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return

*       NOTE: t(O) = R %*% t(P)      
*      call dgemm( 'T','T', p,p,p, 1.d0, temp2,p, temp1,p, 
*     *              0.d0, O,p )
      O = 0.d0
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
*       NOTE: we compute the TRANSPOSED of the matrix in the output in the paper
      call transpose(O, p)
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       compute shape (matrix A) and target function
      trgt = 0.d0
      do k = 1,G
      
        temp1 = 0.d0
        call dgemm( 'N','N', p,p,p, 1.d0, O,p, U(:,:,k),p, 
     *              0.d0, temp1,p )
*       temp1 contains t(D) %*% W
        do j = 1,p
            shape(j,k) = ddot(p, temp1(j,:), 1, O(j,:), 1)
        end do
        shape(:,k) = shape(:,k)/
     *                  exp( sum( log(shape(:,k)) ) )**(1.d0/dble(p))
*       now shape contains matrix A of Celeux, Govaert pag. 785

*       check positive values
        if ( minval(shape(:,k)) .lt. sqrt(eps) ) then
            info = 0
            shape = FLMAX
            return
        end if

        temp4(1) = 0.d0
        do j = 1,p
*            temp2(:,j) = O(:,j) * 1.d0/shape(j,k)
            temp2(:,j) = O(j,:) * 1.d0/shape(j,k)
            temp4(1) = temp4(1) + ddot(p, temp1(j,:), 1, temp2(:,j), 1)
        end do
        trgt = trgt + temp4(1)
        
      end do
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt      

*       WHILE condition      
      if ( errin .gt. tol .and. niterin .lt. itmax ) goto 100 

      scale = trgt / ( sum(sumz)*dble(p) )
        
      return
      end


* ======================================================================
      subroutine eseve (x,z, n,p,G,Gnoise, mu,O,scale,shape,pro, Vinv, 
     *                 loglik, eps)
*       Expectation step for model EVE
* ======================================================================
      
      implicit none
      
      integer :: n, p, G, Gnoise
      double precision :: x(n,p), z(n,Gnoise)
      double precision :: mu(p,G), O(p,p), scale, shape(p,G)
      double precision :: Vinv, pro(Gnoise)
      
      double precision :: temp1(p), temp2(p), temp3, temp4(n)
      
      integer :: i, k, j
      double precision :: const, logdet, loglik, eps
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)      
      
      external :: ddot
      double precision :: ddot

      
*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------
      
*       check very small eigenvalues (singular covariance) 
      if (minval(shape) .le. sqrt(eps) .or. scale .le. sqrt(eps)) then
        loglik = FLMAX
        return
      end if
      
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O, p, temp1, 1, 0.d0, temp2, 1)
            temp2 = temp2/sqrt(scale*shape(:,k))
            temp3 = ddot(p, temp2, 1, temp2, 1)
            temp4(i) = temp3
*       temp3 contains the mahalanobis distance
*            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
            z(i,k) = const - logdet/2.d0 - temp3/2.d0
*       help(cdens) --> The densities are not scaled by mixing proportions            
        end do
*       z contains the log-density log(N(x|theta_k)) 
            
      end do
      
      if ( pro(1) .lt. 0.d0 ) return
*       cdens function

*       noise component
      if (Vinv .gt. 0.d0) then
        call dcopy( n, log(Vinv), 0, z(:,Gnoise), 1)
      end if 
*       now column Gnoise of z contains log(Vinv)           
      
      do i = 1,n
      
        z(i,:) = z(i,:) + log( pro )
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1) )

*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
        
      end do
      
      return
      end
      
      
* ======================================================================
      subroutine meeve ( x,z, n,p,G,Gnoise, mu,O,U,scale,shape,pro,Vinv,
     *                  loglik, eqpro,itmaxin,tolin,itmaxout,tolout,eps, 
     *                  niterin,errin,niterout,errout,lwork,info )
*       Maximization-expectation algorithm for model EVE
* ======================================================================
      
      implicit none
      
      logical :: eqpro
      integer :: n,p,G,Gnoise
      double precision :: x(n,p), z(n,Gnoise), pro(Gnoise), Vinv
      double precision :: mu(p,G), O(p,p), scale, shape(p,G)
      double precision :: U(p,p,G), sumz(Gnoise), omega(G)
      
      double precision :: temp1(p,p), temp2(p,p), temp3(p,p), temp4(p)
      
      integer :: i, j, k, info, lwork
      integer :: itmaxin, itmaxout, niterin, niterout
      double precision :: tolin, tolout, errin, errout, eps, rteps
      double precision :: const, logdet, loglik, lkprev, wrk(lwork)
      double precision :: trgt, trgtprev
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot

*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------

      rteps = sqrt(eps)
      niterout = 0
      errout = FLMAX
      lkprev = FLMAX/2
      loglik = FLMAX
      
      const = (-dble(p)/2.d0)*log2pi
      
*       WHILE loop for EM algorithm  
100   continue
    
      niterout = niterout + 1
      
      sumz = sum(z, dim = 1)
      if ( eqpro ) then
        if ( Vinv .gt. 0 ) then
            pro(Gnoise) = sumz(Gnoise) / dble(n)
            pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
            sumz = pro * dble(n)
        else 
            pro = 1 / dble(G)
            sumz = pro * dble(n)
        end if
      else
        pro = sumz / dble(n)
      end if
      
*       re-initialise U      
      call dcopy(p*p*G, 0.d0, 0, U, 1)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1(:,1), 1, temp1(:,1), 1, 
     *                 U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       compute the eigenvalues of U to be stored in omega
        temp2 = U(:,:,k)
        call dsyev('N', 'U', p, temp2, p, temp1(:,1), wrk, lwork, info)
*       now temp1 contains all the eigenvalues of U

*       check if dsyev converged and positive definite
        if ( info .ne. 0 ) then
            return
        else
            if ( minval(temp1(:,1)) .lt. rteps ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
        
        omega(k) = temp1(p,1)      
        
      end do
*       omega contains the largest eigenvalue of each scattering matrix 

*       M step..........................................................
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2
      
*       covariance matrix components estimation
*       we consider algorithm MM 1 and MM 2 of Browne, McNicholas 2013
*       with a modification in computing the orientation matrix in the MM 2 step

*       shape (matrix A) and orientation (matrix D) initialised in R        
*       shape = matrix(1, p,G)
*       O = diag(p)             
*       ### NOTE: we don't re-initialize shape and orientation at each 
*                 outer iteration of the EM algorithm
      
*       WHILE loop for M step
110   continue

*       ### NOTE: O is transposed

      niterin = niterin + 1

      temp2 = 0.d0
      temp3 = 0.d0
*       temp3 will contain matrix F      
      
*       Algorithm MM 1 ......................................
      do k = 1,G
        
        do j = 1,p
*            temp1(j,:) = O(:,j) / shape(j,k)
            temp1(j,:) = O(j,:) / shape(j,k)
        end do
*       temp1 contains inv(A)t(D)
        
        call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, U(:,:,k),p, 
     *              0.d0, temp2,p )
*       temp2 contains inv(A) %*% t(D) %*% W    
        
        temp1 = temp2 - omega(k)*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do
      
*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return
      
*       NOTE: t(O) = t( R %*% t(P) ) = P %*% t(R)      
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       Algorithm MM 2 ......................................
*      call dgemm( 'T','T', p,p,p, 1.d0, temp2,p, temp1,p, 
*     *              0.d0, O,p )
      call transpose(O, p)
*       O contains matrix D of Browne, McNicholas
  
*       Algorithm MM 2
      temp1 = 0.d0
      temp3 = 0.d0
      do k = 1,G
        
      call dgemm( 'N','N', p,p,p, 1.d0, U(:,:,k),p, O,p, 
     *            0.d0, temp1,p )
*       temp1 contains W %*% D
        
        do j = 1,p
            temp2(:,j) = temp1(:,j) / shape(j,k)
        end do
*       temp2 contains W %*% D %*% inv(A)
        
        temp1 = temp2 - maxval( 1/shape(:,k) )*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do

*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return

*       NOTE: t(O) = R %*% t(P)      
*      call dgemm( 'T','T', p,p,p, 1.d0, temp2,p, temp1,p, 
*     *              0.d0, O,p )
      O = 0.d0
*       NOTE: we compute the TRANSPOSED of the matrix in the output in the paper
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
      call transpose(O, p)
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       compute shape (matrix A) and target function
      trgt = 0.d0
      do k = 1,G
      
        temp1 = 0.d0
        call dgemm( 'N','N', p,p,p, 1.d0, O,p, U(:,:,k),p, 
     *              0.d0, temp1,p )
*       temp1 contains t(D) %*% W
        do j = 1,p
            shape(j,k) = ddot(p, temp1(j,:), 1, O(j,:), 1)
        end do
        shape(:,k) = shape(:,k)/
     *                  exp( sum( log(shape(:,k)) ) )**(1.d0/dble(p))
*       now shape contains matrix A of Celeux, Govaert pag. 785

*       check positive values
        if ( minval(shape(:,k)) .lt. rteps ) then
            info = 0
            loglik = FLMAX
            return
        end if

        temp4(1) = 0.d0
        do j = 1,p
*            temp2(:,j) = O(:,j) * 1.d0/shape(j,k)
            temp2(:,j) = O(j,:) * 1.d0/shape(j,k)
            temp4(1) = temp4(1) + ddot(p, temp1(j,:), 1, temp2(:,j), 1)
        end do
        trgt = trgt + temp4(1)
        
      end do
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt      

*       WHILE condition M step     
      if ( errin .gt. tolin .and. niterin .lt. itmaxin ) goto 110 

      scale = trgt / ( sum(sumz(1:G))*dble(p) )
*       ................................................................       

*       E step..........................................................
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2(:,1), 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O, p, temp1(:,1), 1, 0.d0, temp2(:,1), 1)
            temp2(:,1) = temp2(:,1)/sqrt(scale*shape(:,k))
            temp3(1,1) = ddot(p, temp2(:,1), 1, temp2(:,1), 1)
*       temp3 contains the mahalanobis distance
            z(i,k) = const - logdet/2.d0 - temp3(1,1)/2.d0 + log(pro(k))
*            z(i,k) = const - logdet/2.d0 - temp3(1,1)/2.d0         
        end do
*       z contains the log-density log(N(x|theta_k)) + log(p_k) 
            
      end do

*       noise component
      if (Vinv .gt. 0.d0) then
        z(:,Gnoise) = log(Vinv) + log( pro(Gnoise) )
      end if 
*       now column Gnoise of z contains log(Vinv) + log(p_0)          
      
      loglik = 0.d0
      do i = 1,n
*       Numerical Recipes pag.844
        temp3(1,1) = maxval(z(i,:))
        temp1(1,1) = temp3(1,1) + log( sum(exp(z(i,:) - temp3(1,1))) )
        loglik = loglik + temp1(1,1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1,1) )

*       re-normalize probabilities
        temp3(1,1) = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3(1,1), z(i,:), 1 )
      end do
*       ................................................................        

      errout = abs(loglik - lkprev)/(1.d0 + abs(loglik))
      lkprev = loglik

* Chris F (June 2015): pro should not be computed in the E-step
*     sumz = sum(z, dim = 1)
*     if ( eqpro ) then
*       if ( Vinv .gt. 0 ) then
*           pro(Gnoise) = sumz(Gnoise) / dble(n)
*           pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
*           sumz = pro * dble(n)
*       else 
*           pro = 1 / dble(G)
*           sumz = pro * dble(n)
*       end if
*     else
*       pro = sumz / dble(n)
*     end if
      
*       check if empty components
      if ( minval(sumz) .lt. rteps ) then
        loglik = -FLMAX
        return
      end if
      
*       WHILE condition EM     
      if ( errout .gt. tolout .and. niterout .lt. itmaxout ) goto 100

      return
      end
      

************************************************************************
**** VVE model
************************************************************************
* ======================================================================
      subroutine msvve (x,z, n,p,G, mu,U,O,scale,shape,pro, lwork,info,
     *                  itmax,tol, niterin,errin, eps)
*       Maximization step for model VVE
* ======================================================================

      implicit none
      
      integer :: n, p, G
      
      double precision :: x(n,p), z(n,G)
      
      double precision :: mu(p,G), U(p,p,G), pro(G), O(p,p)
      double precision :: scale(G), shape(p,G)
      
      double precision :: sumz(G), omega(G)
      integer :: i, j, k, info, lwork
          
      double precision :: temp1(p,p), temp2(p,p), temp3(p,p), temp4(p) 
      double precision :: wrk(lwork), tol, errin, trgt, trgtprev, eps
      integer :: itmax, niterin
    
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot
        
*-----------------------------------------------------------------------

*       colsums of z
      sumz = sum(z, dim = 1)
      
*       a priori probabilities
      pro = sumz / dble(n)
*      pro = sumz / sum(sumz)
*       if there is noise sum(sumz) does not sum to n. See help(mstep)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1(:,1), 1, temp1(:,1), 1, 
     *                 U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       compute the eigenvalues of U to be stored in omega
        temp2 = U(:,:,k)
        call dsyev('N', 'U', p, temp2, p, temp1(:,1), wrk, lwork, info)
*       now temp1 contains all the eigenvalues of U

*       check if dsyev converged and positive definite
        if ( info .ne. 0 ) then
            return
        else
            if ( minval(temp1(:,1)) .lt. sqrt(eps) ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
        
        omega(k) = temp1(p,1)      
        
      end do
*       omega contains the largest eigenvalue of each scattering matrix
   
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2

*       covariance matrix components estimation
*       we consider algorithm MM 1 and MM 2 of Browne, McNicholas 2013
*       with a modification in computing the orientation matrix in the MM 2 step

*      shape (matrix A) and orientation (matrix D) initialised in R        
*      shape = matrix(1, p,G)
*      O = diag(p)      
      
*       WHILE loop using goto statement
100   continue

*       ### NOTE: O is transposed

      niterin = niterin + 1

      temp2 = 0.d0
      temp3 = 0.d0
*       temp3 will contain matrix F      
      
*       Algorithm MM 1 ......................................      
      do k = 1,G
        
        do j = 1,p
*            temp1(j,:) = O(:,j) / shape(j,k)
            temp1(j,:) = O(j,:) / shape(j,k)
        end do
*       temp1 contains inv(A)t(D)
        
        call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, U(:,:,k),p, 
     *              0.d0, temp2,p )
*       temp2 contains inv(A) %*% t(D) %*% W    
        
        temp1 = temp2 - omega(k)*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do
      
*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return
      
*       NOTE: t(P %*% t(R)) = R %*% t(P)      
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       Algorithm MM 2 ......................................
      call transpose(O, p)
*       O contains matrix D of Browne, McNicholas
  
*       Algorithm MM 2
      temp1 = 0.d0
      temp3 = 0.d0
      do k = 1,G
        
      call dgemm( 'N','N', p,p,p, 1.d0, U(:,:,k),p, O,p, 
     *            0.d0, temp1,p )
*       temp1 contains W %*% D
        
        do j = 1,p
            temp2(:,j) = temp1(:,j) / shape(j,k)
        end do
*       temp2 contains W %*% D %*% inv(A)
        
        temp1 = temp2 - maxval( 1/shape(:,k) )*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do

*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return

*       NOTE: t(O) = R %*% t(P)      
      O = 0.d0
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
      call transpose(O, p)
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       compute shape (matrix A) and target function
      trgt = 0.d0
      do k = 1,G
      
        temp1 = 0.d0
        call dgemm( 'N','N', p,p,p, 1.d0, O,p, U(:,:,k),p, 
     *              0.d0, temp1,p )
*       temp1 contains t(D) %*% W
        do j = 1,p
            shape(j,k) = ddot(p, temp1(j,:), 1, O(j,:), 1)
        end do
*        shape(:,k) = shape(:,k)/
*     *                  exp( sum( log(shape(:,k)) ) )**(1.d0/dble(p))
        shape(:,k) = shape(:,k)/sumz(k)
*       now shape contains matrix A (scale*A) of Celeux, Govaert pag. 785

*       compute scale parameter and shape matrix A
        scale(k) = exp( sum( log(shape(:,k)) ) )**(1/dble(p))
        shape(:,k) = shape(:,k)/scale(k)

*       check positive values
        if ( minval(shape(:,k)) .lt. sqrt(eps) ) then
            info = 0
            shape = FLMAX
            return
        end if

        temp4(1) = 0.d0
        do j = 1,p
*            temp2(:,j) = O(:,j) * 1.d0/shape(j,k)
            temp2(:,j) = O(j,:) * 1.d0/shape(j,k)
            temp4(1) = temp4(1) + ddot(p, temp1(j,:), 1, temp2(:,j), 1)
        end do
        trgt = trgt + temp4(1)
        
      end do
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt      

*       WHILE condition      
      if ( errin .gt. tol .and. niterin .lt. itmax ) goto 100 
        
      return
      end


* ======================================================================
      subroutine esvve (x,z, n,p,G,Gnoise, mu,O,scale,shape,pro, Vinv, 
     *                 loglik, eps)
*       Expectation step for model VVE
* ======================================================================
      
      implicit none
      
      integer :: n, p, G, Gnoise
      double precision :: x(n,p), z(n,Gnoise)
      double precision :: mu(p,G), O(p,p), scale(G), shape(p,G)
      double precision :: Vinv, pro(Gnoise)
      
      double precision :: temp1(p), temp2(p), temp3, temp4(n)
      
      integer :: i, k, j
      double precision :: const, logdet, loglik, eps
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)      
      
      external :: ddot
      double precision :: ddot

      
*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------
      
*       check very small eigenvalues (singular covariance) 
      if ( minval(shape) .le. sqrt(eps) .or. 
     *               minval(scale) .le. sqrt(eps) ) then
        loglik = FLMAX
        return
      end if
      
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale(k)) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1 = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2, 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O, p, temp1, 1, 0.d0, temp2, 1)
            temp2 = temp2/sqrt(scale(k)*shape(:,k))
            temp3 = ddot(p, temp2, 1, temp2, 1)
            temp4(i) = temp3
*       temp3 contains the mahalanobis distance
*            z(i,k) = const - logdet/2.d0 - temp3/2.d0 + log(pro(k))
            z(i,k) = const - logdet/2.d0 - temp3/2.d0
*       help(cdens) --> The densities are not scaled by mixing proportions            
        end do
*       z contains the log-density log(N(x|theta_k)) 
            
      end do
      
      if ( pro(1) .lt. 0.d0 ) return
*       cdens function

*       noise component
      if (Vinv .gt. 0.d0) then
        call dcopy( n, log(Vinv), 0, z(:,Gnoise), 1)
      end if 
*       now column Gnoise of z contains log(Vinv)           
      
      do i = 1,n
      
        z(i,:) = z(i,:) + log( pro )
*       Numerical Recipes pag.844
        temp3 = maxval(z(i,:))
        temp1(1) = temp3 + log( sum(exp(z(i,:) - temp3)) )
        loglik = loglik + temp1(1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1) )

*       re-normalize probabilities
        temp3 = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3, z(i,:), 1 ) 
        
      end do
      
      return
      end
      
      
* ======================================================================
      subroutine mevve ( x,z, n,p,G,Gnoise, mu,O,U,scale,shape,pro,Vinv,
     *                  loglik, eqpro,itmaxin,tolin,itmaxout,tolout,eps, 
     *                  niterin,errin,niterout,errout,lwork,info)
*       Maximization-expectation algorithm for model VVE
* ======================================================================
      
      implicit none
      
      logical :: eqpro
      integer :: n,p,G,Gnoise
      double precision :: x(n,p), z(n,Gnoise), pro(Gnoise), Vinv
      double precision :: mu(p,G), O(p,p), scale(G), shape(p,G)
      double precision :: U(p,p,G), sumz(Gnoise), omega(G)
      
      double precision :: temp1(p,p), temp2(p,p), temp3(p,p), temp4(p)
      
      integer :: i, j, k, info, lwork
      integer :: itmaxin, itmaxout, niterin, niterout
      double precision :: tolin, tolout, errin, errout, eps, rteps
      double precision :: const, logdet, loglik, lkprev, wrk(lwork)
      double precision :: trgt, trgtprev
      
      double precision :: log2pi
      parameter (log2pi = 1.837877066409345d0)
      
      double precision :: FLMAX
      parameter (FLMAX = 1.7976931348623157d308)
      
      external :: ddot
      double precision :: ddot

*      double precision :: smalog
*      parameter (smalog = -708.d0)
      
*-----------------------------------------------------------------------

      rteps = sqrt(eps)
      niterout = 0
      errout = FLMAX
      lkprev = FLMAX/2
      loglik = FLMAX
      
      const = (-dble(p)/2.d0)*log2pi
      
*       WHILE loop for EM algorithm  
100   continue
    
      niterout = niterout + 1
      
      sumz = sum(z, dim = 1)
      if ( eqpro ) then
        if ( Vinv .gt. 0 ) then
            pro(Gnoise) = sumz(Gnoise) / dble(n)
            pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
            sumz = pro * dble(n)
        else 
            pro = 1 / dble(G)
            sumz = pro * dble(n)
        end if
      else
        pro = sumz / dble(n)
      end if
      
*       re-initialise U      
      call dcopy(p*p*G, 0.d0, 0, U, 1)
      
*       compute weighted scattering matrix and means      
      do k = 1,G
      
        do j = 1,p
            mu(j,k) = sum(x(:,j)*z(:,k))/sumz(k)
        end do
        
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) ) * sqrt(z(i,k))
            call dger(p, p, 1.d0, temp1(:,1), 1, temp1(:,1), 1, 
     *                 U(:,:,k), p)
        end do
*       U contains the weighted scattering matrix  

*       compute the eigenvalues of U to be stored in omega
        temp2 = U(:,:,k)
        call dsyev('N', 'U', p, temp2, p, temp1(:,1), wrk, lwork, info)
*       now temp1 contains all the eigenvalues of U

*       check if dsyev converged and positive definite
        if ( info .ne. 0 ) then
            return
        else
            if ( minval(temp1(:,1)) .lt. rteps ) then
                info = 0
                scale = FLMAX
                return
            end if
        end if
        
        omega(k) = temp1(p,1)      
        
      end do
*       omega contains the largest eigenvalue of each scattering matrix 

*       M step..........................................................
      niterin = 0
      errin = FLMAX
      trgt = FLMAX
      trgtprev = FLMAX/2
      
*       covariance matrix components estimation
*       we consider algorithm MM 1 and MM 2 of Browne, McNicholas 2013
*       with a modification in computing the orientation matrix in the MM 2 step

*       shape (matrix A) and orientation (matrix D) initialised in R        
*       shape = matrix(1, p,G)
*       O = diag(p)             
*       ### NOTE: we don't re-initialize shape and orientation at each 
*                 outer iteration of the EM algorithm
      
*       WHILE loop for M step
110   continue

*       ### NOTE: O is transposed

      niterin = niterin + 1

      temp2 = 0.d0
      temp3 = 0.d0
*       temp3 will contain matrix F      
      
*       Algorithm MM 1 ......................................      
      do k = 1,G
        
        do j = 1,p
*            temp1(j,:) = O(:,j) / shape(j,k)
            temp1(j,:) = O(j,:) / shape(j,k)
        end do
*       temp1 contains inv(A)t(D)
        
        call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, U(:,:,k),p, 
     *              0.d0, temp2,p )
*       temp2 contains inv(A) %*% t(D) %*% W    
        
        temp1 = temp2 - omega(k)*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do
      
*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return
      
*       NOTE: t(P %*% t(R)) = R %*% t(P)      
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
*       O contains TRANSPOSED orientation (matrix D of Browne, McNicholas)
*       .....................................................

*       Algorithm MM 2 ......................................
      call transpose(O, p)
*       O contains matrix D of Browne, McNicholas
  
*       Algorithm MM 2
      temp1 = 0.d0
      temp3 = 0.d0
      do k = 1,G
        
      call dgemm( 'N','N', p,p,p, 1.d0, U(:,:,k),p, O,p, 
     *            0.d0, temp1,p )
*       temp1 contains W %*% D
        
        do j = 1,p
            temp2(:,j) = temp1(:,j) / shape(j,k)
        end do
*       temp2 contains W %*% D %*% inv(A)
        
        temp1 = temp2 - maxval( 1/shape(:,k) )*temp1
        temp3 = temp3 + temp1
*       temp3 contains the matrix F          
        
      end do

*       compute matrices P and R where svd(F) = P %*% B %*% t(R)
      call dgesvd('A','A', p,p, temp3,p, temp4, temp1,p, temp2,p,
     *              wrk, lwork, info)
*       now temp1 contains matrix P, temp2 contains matrix t(R)
*       temp4 contains the singular values

*       check if dgesvd converged
      if ( info .ne. 0 ) return

*       NOTE: t(O) = R %*% t(P)      
      O = 0.d0
      call dgemm( 'N','N', p,p,p, 1.d0, temp1,p, temp2,p, 
     *              0.d0, O,p )
      call transpose(O, p)
*       O contains TRANSPOSED matrix D of Browne, McNicholas
*       .....................................................

*       compute shape (matrix A) and target function
      trgt = 0.d0
      do k = 1,G
      
        temp1 = 0.d0
        call dgemm( 'N','N', p,p,p, 1.d0, O,p, U(:,:,k),p, 
     *              0.d0, temp1,p )
*       temp1 contains t(D) %*% W
        do j = 1,p
            shape(j,k) = ddot(p, temp1(j,:), 1, O(j,:), 1)
        end do
*        shape(:,k) = shape(:,k)/
*     *                  exp( sum( log(shape(:,k)) ) )**(1.d0/dble(p))
        shape(:,k) = shape(:,k)/sumz(k)
*       now shape contains matrix A (scale*A) of Celeux, Govaert pag. 785

*       compute scale parameter and shape matrix A
        scale(k) = exp( sum( log(shape(:,k)) ) )**(1/dble(p))
        shape(:,k) = shape(:,k)/scale(k)

*       check positive values
        if (minval(shape(:,k)) .lt. rteps .or.
     *                      scale(k) .lt. rteps) then
            info = 0
            loglik = FLMAX
            return
        end if

        temp4(1) = 0.d0
        do j = 1,p
*            temp2(:,j) = O(:,j) * 1.d0/shape(j,k)
            temp2(:,j) = O(j,:) * 1.d0/shape(j,k)
            temp4(1) = temp4(1) + ddot(p, temp1(j,:), 1, temp2(:,j), 1)
        end do
        trgt = trgt + temp4(1)
        
      end do
      
*       error
      errin = abs(trgt - trgtprev)/(1.d0 + abs(trgt))
      
      trgtprev = trgt      

*       WHILE condition M step     
      if ( errin .gt. tolin .and. niterin .lt. itmaxin ) goto 110
      
*      do k = 1,G
*        scale(k) = exp( sum( log(shape(:,k)) ) )**(1/dble(p))
*        shape(:,k) = shape(:,k)/scale(k) 
*      end do
*       ................................................................       

*       E step..........................................................
      const = (-dble(p)/2.d0)*log2pi
      
      do  k = 1,G
        
        logdet = 0.d0
        do j = 1,p        
            logdet = logdet + ( log(shape(j,k)) + log(scale(k)) )
        end do

*       compute mahalanobis distance for each observation
*       ##### NOTE: O is transposed   
        do i = 1,n
            temp1(:,1) = ( x(i,:) - mu(:,k) )
            call dcopy(p, 0.d0, 0, temp2(:,1), 1)
            call dgemv('N', p, p, 1.d0, 
     *                 O, p, temp1(:,1), 1, 0.d0, temp2(:,1), 1)
            temp2(:,1) = temp2(:,1)/sqrt(scale(k)*shape(:,k))
            temp3(1,1) = ddot(p, temp2(:,1), 1, temp2(:,1), 1)
*       temp3 contains the mahalanobis distance
            z(i,k) = const - logdet/2.d0 - temp3(1,1)/2.d0 + log(pro(k))
*            z(i,k) = const - logdet/2.d0 - temp3(1,1)/2.d0         
        end do
*       z contains the log-density log(N(x|theta_k)) + log(p_k) 
            
      end do

*       noise component
      if (Vinv .gt. 0.d0) then
        z(:,Gnoise) = log(Vinv) + log( pro(Gnoise) )
      end if 
*       now column Gnoise of z contains log(Vinv) + log(p_0)          
      
      loglik = 0.d0
      do i = 1,n
*       Numerical Recipes pag.844
        temp3(1,1) = maxval(z(i,:))
        temp1(1,1) = temp3(1,1) + log( sum(exp(z(i,:) - temp3(1,1))) )
        loglik = loglik + temp1(1,1)
*       ##### NOTE: do we need to check if (z - zmax) is too small?
        
        z(i,:) = exp( z(i,:) - temp1(1,1) )

*       re-normalize probabilities
        temp3(1,1) = sum( z(i,:) )
        call dscal( Gnoise, 1.d0/temp3(1,1), z(i,:), 1 )
      end do
*       ................................................................        

      errout = abs(loglik - lkprev)/(1.d0 + abs(loglik))
      lkprev = loglik

* Chris F (June 2015): pro should not be computed in the E-step
*     sumz = sum(z, dim = 1)
*     if ( eqpro ) then
*       if ( Vinv .gt. 0 ) then
*           pro(Gnoise) = sumz(Gnoise) / dble(n)
*           pro(1:G) = ( 1 - pro(Gnoise) ) / dble(G)
*           sumz = pro * dble(n)
*       else 
*           pro = 1 / dble(G)
*           sumz = pro * dble(n)
*       end if
*     else
*       pro = sumz / dble(n)
*     end if
      
*       check if empty components
      if ( minval(sumz) .lt. rteps ) then
        loglik = -FLMAX
        return
      end if
      
*       WHILE condition EM     
      if ( errout .gt. tolout .and. niterout .lt. itmaxout ) goto 100

      return
      end

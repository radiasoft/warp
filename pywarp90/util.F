#include "top.h"
c************************************************************************** 
c@(#) File UTIL.M, version $Revision: 1.1 $ $Date: 2001/06/30 00:53:59 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c Contains various utility routines that are part of the TOP package.
c This file should rarely be changed, if ever, except to add to it.
c sdot was replaced by sdotlocal to aboid naming problems.
c
c Currently contains:
c     From LINPACK:
c        ssifa(a,lda,n,kpvt,info)                 subroutine
c        ssidi(a,lda,n,kpvt,det,inert,work,job)   subroutine
c        ssisl(a,lda,n,kpvt,b)                    subroutine
c        sswap (n,sx,incx,sy,incy)                subroutine
c        saxpy(n,sa,sx,incx,sy,incy)              subroutine
c        scopy(n,sx,incx,sy,incy)                 subroutine
c        isamax(n,sx,incx)                        integer function
c        sdotlocal(n,sx,incx,sy,incy)             real function
c     From "Numerical Recipes"
c        svbksb(u,w,v,m,n,mp,np,b,x,tmp)          subroutine 
c        svdcmp(a,m,n,mp,np,w,v,tmp)              subroutine
c        svdfit(x,y,ndata,ndatap,a,basis,ma,map,  subroutine 
c               u,w,v,tmp,tol,chisq)
c
c
c***************************************************************************

c** from netlib, Fri Jul 27 13:53:38 EDT 1990 ***
      subroutine ssifa(a,lda,n,kpvt,info)
      integer(ISZ):: lda,n,kpvt(1),info
      real(kind=8):: a(lda,1)
c
c     ssifa factors a real symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow ssifa by ssisl.
c     to compute  inverse(a)*c , follow ssifa by ssisl.
c     to compute  determinant(a) , follow ssifa by ssidi.
c     to compute  inertia(a) , follow ssifa by ssidi.
c     to compute  inverse(a) , follow ssifa by ssidi.
c
c     on entry
c
c        a       real(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that ssisl or ssidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas saxpy,sswap,isamax
c     fortran abs,max,sqrt
c
c     internal variables
c
      real(kind=8):: ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real(kind=8):: absakk,alpha,colmax,rowmax
      integer(ISZ):: imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,isamax
      logical(ISZ):: swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = abs(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = isamax(k-1,a(1,k),1)
         colmax = abs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = max(rowmax,abs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = isamax(imax-1,a(1,imax),1)
               rowmax = max(rowmax,abs(a(jmax,imax)))
   50       continue
            if (abs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (max(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call sswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call saxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call sswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0e0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call saxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call saxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      integer(ISZ) function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      real(kind=8):: sx(1),smax
      integer(ISZ):: i,incx,ix,n
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
c** from netlib, Fri Jul 27 13:53:22 EDT 1990 ***
      subroutine ssidi(a,lda,n,kpvt,det,inert,work,job)
      integer(ISZ):: lda,n,job
      real(kind=8):: a(lda,1),work(1)
      real(kind=8):: det(2)
      integer(ISZ):: kpvt(1),inert(3)
c
c     ssidi computes the determinant, inertia and inverse
c     of a real symmetric matrix using the factors from ssifa.
c
c     on entry
c
c        a       real(lda,n)
c                the output from ssifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from ssifa.
c
c        work    real(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    real(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  ssico  has set rcond .eq. 0.0
c        or  ssifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas saxpy,scopy,sdotlocal,sswap
c     fortran abs,iabs,mod
c
c     internal variables.
c
      real(kind=8):: akkp1,sdotlocal,temp
      real(kind=8):: ten,d,t,ak,akp1
      integer(ISZ):: j,jb,k,km1,ks,kstep
      logical(ISZ):: noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
   20    continue
         t = 0.0e0
         do 130 k = 1, n
            d = a(k,k)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = abs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0e0) go to 30
                  t = abs(a(k,k+1))
                  d = (d/t)*a(k+1,k+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0e0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0e0) inert(1) = inert(1) + 1
               if (d .lt. 0.0e0) inert(2) = inert(2) + 1
               if (d .eq. 0.0e0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0e0) go to 110
   70             if (abs(det(1)) .ge. 1.0e0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0e0
                  go to 70
   80             continue
   90             if (abs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0e0
                  go to 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
  150    if (k .gt. n) go to 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               a(k,k) = 1.0e0/a(k,k)
               if (km1 .lt. 1) go to 170
                  call scopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = sdotlocal(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + sdotlocal(km1,work,1,a(1,k),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = abs(a(k,k+1))
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0e0)
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call scopy(km1,a(1,k+1),1,work,1)
                  do 190 j = 1, km1
                     a(j,k+1) = sdotlocal(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1) + sdotlocal(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + sdotlocal(km1,a(1,k),1,a(1,k+1),1)
                  call scopy(km1,a(1,k),1,work,1)
                  do 200 j = 1, km1
                     a(j,k) = sdotlocal(j,a(1,j),1,work,1)
                     call saxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k) + sdotlocal(km1,work,1,a(1,k),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               call sswap(ks,a(1,ks),1,a(1,k),1)
               do 230 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) go to 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
      subroutine sswap (n,sx,incx,sy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real(kind=8):: sx(1),sy(1),stemp
      integer(ISZ):: i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real(kind=8):: sx(1),sy(1),sa
      integer(ISZ):: i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end
      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real(kind=8):: sx(1),sy(1)
      integer(ISZ):: i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
      real(kind=8) function sdotlocal(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real(kind=8):: sx(1),sy(1),stemp
      integer(ISZ):: i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdotlocal = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdotlocal = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdotlocal = stemp
      return
      end

c***************************************************************************
c** from netlib, Mon Oct 21 15:51:27 EDT 1991 ***
      subroutine ssisl(a,lda,n,kpvt,b)
      integer(ISZ):: lda,n,kpvt(1)
      real(kind=8):: a(lda,1),b(1)
c
c     ssisl solves the real symmetric system
c     a * x = b
c     using the factors computed by ssifa.
c
c     on entry
c
c        a       real(lda,n)
c                the output from ssifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from ssifa.
c
c        b       real(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  ssico  has set rcond .eq. 0.0
c        or  ssifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call ssifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call ssisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas saxpy,sdotlocal
c     fortran iabs
c
c     internal variables.
c
      real(kind=8):: ak,akm1,bk,bkm1,sdotlocal,denom,temp
      integer(ISZ):: k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call saxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call saxpy(k-2,b(k),a(1,k),1,b(1),1)
               call saxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + sdotlocal(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + sdotlocal(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + sdotlocal(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
c***********************************************************************
      subroutine svdfit(x,y,ndata,ndatap,a,basis,ma,map,u,w,v,tmp,
     *                  tol,chisq)
      integer(ISZ):: ndata,ndatap,map,ma
      real(kind=8):: tol,chisq
      real(kind=8):: x(ndatap),y(ndatap),a(map),basis(ndatap,map),
     *     u(ndatap,map),w(map),v(map,map), tmp(ndatap)
c     
c     Solution of linear least-squares matrix problem by the use of 
c     singular value decompositions.   
c     Taken from Press et. al "Numerical Recipies, The Art of Scientific 
c     Computing," Cambridge University Press, copyright 1986, 
c     3rd printing 1988, pg 518.  Program has been extensively modified.  
c     Changes include: explicit input of the basis matrix and the singular 
c     value decomposition of the basis matrix, physical and logical 
c     dimensions passed for all arrays, elimination of variances for each 
c     data point, explicit input arguments for workspace vector tmp passed 
c     to svbksb and tolerance input for singular values. 
c
c     Given a set of ndata points x(i), y(i), use chi-squared minimization 
c     to determine the ma coefficients of a of the fitting function 
c     y = \sum_i a_i basis_i(x).  Here, the sum over i is over the ma 
c     fit coefficients and basis_i(x) denotes the ith basis function.  
c     These basis functions are input as a ndata x ma matrix basis where 
c     each row corresponds to the ma basis functions basis_i(x) evaluated 
c     at one of the ndata x data points.  The least-squares problem is 
c     solved using the singular value decomposition of the matrix basis, 
c     defined as basis = u.w.transpose(v),as provided by the routine svdcmp.  
c     u, w, and v are input.  ma and ndata, and map and ndatap are the 
c     logical and physical dimensions of these matrices, respectively, 
c     as indicated below.  It is necessary that ndata >= ma.  The program 
c     returns the values for the ma fit parameters a, and the chi-squared 
c     error measure chisq.  
c
      real(kind=8):: wmax,thresh,sum
      integer(ISZ):: i,j

      wmax=0.
      do 13 j=1,ma
        if(w(j).gt.wmax)wmax=w(j)
13    continue
      thresh=tol*wmax
      do 14 j=1,ma
        if(w(j).lt.thresh)w(j)=0.
14    continue
      call svbksb(u,w,v,ndata,ma,ndatap,map,y,a,tmp)
      chisq=0.
      do 16 i=1,ndata
        sum=0.
        do 15 j=1,ma
          sum=sum+a(j)*basis(i,j)
15      continue
        chisq=chisq+(y(i)-sum)**2
16    continue
      return
      end
c***********************************************************************
      subroutine svdcmp(a,m,n,mp,np,w,v,tmp)
      integer(ISZ):: m,n,mp,np
      real(kind=8):: a(mp,np),w(np),v(np,np),tmp(mp)
c     
c     Singular value decomposition of a matrix 
c     Taken from Press et. al "Numerical Recipies, The Art of Scientific 
c     Computing," Cambridge University Press, copyright 1986, 
c     3rd printing 1988, pg 60.  Program error trap is modified, and 
c     the argument list has been expanded to include work space vector TMP.   
c
c     Given a matrix a, with logical dimensions m by n and physical 
c     dimensions mp by np, this routine computes its singular value 
c     decomposition, a = u.w.transpose(v).  The matrix u replaces a on 
c     output.  The diagonal matrix of singular values is output as a 
c     vector w.  The matrix v (not the transpose of v) is output as v.  
c     m must be greater or equal to n; if it is smaller, then a should be 
c     filled up to a square with zero rows.  
c
      real(kind=8):: g,scale,anorm,s,f,h,c,x,y,z
      integer(ISZ):: i,j,l,its,nm,k

      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        tmp(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if (i.le.m) then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if (scale.ne.0.0) then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            if (i.ne.n) then
              do 15 j=l,n
                s=0.0
                do 13 k=i,m
                  s=s+a(k,i)*a(k,j)
13              continue
                f=s/h
                do 14 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
            endif
            do 16 k= i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if ((i.le.m).and.(i.ne.n)) then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if (scale.ne.0.0) then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              tmp(k)=a(i,k)/h
19          continue
            if (i.ne.m) then
              do 23 j=l,m
                s=0.0
                do 21 k=l,n
                  s=s+a(j,k)*a(i,k)
21              continue
                do 22 k=l,n
                  a(j,k)=a(j,k)+s*tmp(k)
22              continue
23            continue
            endif
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(tmp(i))))
25    continue
      do 32 i=n,1,-1
        if (i.lt.n) then
          if (g.ne.0.0) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=tmp(i)
        l=i
32    continue
      do 39 i=n,1,-1
        l=i+1
        g=w(i)
        if (i.lt.n) then
          do 33 j=l,n
            a(i,j)=0.0
33        continue
        endif
        if (g.ne.0.0) then
          g=1.0/g
          if (i.ne.n) then
            do 36 j=l,n
              s=0.0
              do 34 k=l,m
                s=s+a(k,i)*a(k,j)
34            continue
              f=(s/a(i,i))*g
              do 35 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
          endif
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if ((abs(tmp(l))+anorm).eq.anorm)  go to 2
            if ((abs(w(nm))+anorm).eq.anorm)  go to 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*tmp(i)
            if ((abs(f)+anorm).ne.anorm) then
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.0/h
              c= (g*h)
              s=-(f*h)
              do 42 j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
42            continue
            endif
43        continue
2         z=w(k)
          if (l.eq.k) then
            if (z.lt.0.0) then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            go to 3
          endif
          if (its.eq.30) then 
            print*,"Routine svdcmp - No convergence in 30 iterations"
            call kaboom (0) 
          endif 
          x=w(l)
          nm=k-1
          y=w(nm)
          g=tmp(nm)
          h=tmp(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=sqrt(f*f+1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=tmp(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            tmp(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 nm=1,n
              x=v(nm,j)
              z=v(nm,i)
              v(nm,j)= (x*c)+(z*s)
              v(nm,i)=-(x*s)+(z*c)
45          continue
            z=sqrt(f*f+h*h)
            w(j)=z
            if (z.ne.0.0) then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 nm=1,m
              y=a(nm,j)
              z=a(nm,i)
              a(nm,j)= (y*c)+(z*s)
              a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
          tmp(l)=0.0
          tmp(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      end
c***********************************************************************
      subroutine svbksb(u,w,v,m,n,mp,np,b,x,tmp)
      integer(ISZ):: m,n,mp,np
      real(kind=8):: u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(mp)
c     
c     Singular value "backsubstitution" routine for the solution of matrix 
c     problems with singular value matrix decompositions. 
c     Taken from Press et. al "Numerical Recipies, The Art of Scientific 
c     Computing," Cambridge University Press, copyright 1986, 
c     3rd printing 1988, pg 57.  Program argument list has been expanded 
c     to include work space vector tmp.   
c
c     Solves the matrix problem a.x = b for a vector x, where a is 
c     specified by the arrays u, w, v as returned by the singular value 
c     matrix decomposition routine svdcmp.  m and n are the logical 
c     dimensions of a, and will be equal for square matrices.  mp and np
c     are the physical dimensions of a.  b is the input right-hand side.  
c     x is the output solution vector.  No input quantities are destroyed, 
c     so the routine may be called sequentially with different b's.  m 
c     must be greater or equal to n; see svdcmp.
c
      integer(ISZ):: i,j,jj
      real(kind=8):: s
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      end
c***********************************************************************

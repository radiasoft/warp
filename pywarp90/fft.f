#include "top.h"
c************************************************************************** 
c@(#) File FFT.M, version $Revision: 1.1 $ $Date: 2001/06/30 00:53:58 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c Contains FFT routines that are part of the TOP package.
c This file should rarely be changed, if ever, accept to add to it.
c
c Currently contains:
c     FFT routines:
c        vrpft2                                   subroutine
c        vrpfti2                                  subroutine
c        vcpft                                    subroutine
c
c***************************************************************************
c========================================================================
c    The following three routines, vrpft2, vrpfti2, vcpft are the FFT
c    routines used in the field solvers, f3d and frz.
c========================================================================

      subroutine vrpft2 (a, b, n, incp, lenv, lfd)
      integer(ISZ):: n,incp,lenv,lfd
      real(kind=8):: a(lfd,1), b(lfd,1)

c Vectorized version of Langdon's RPFT2
c Call after VCPFT
c Written by Alex Friedman, LLNL, April 1989.

      integer(ISZ):: inc,ninc,iv,lp,lm
      real(kind=8):: rp,rm,ip,im

      inc = incp
      ninc = n * inc
      do 100 iv = 1, lenv
         a(1,iv) = 2. * a(1,iv)
         b(1,iv) = 2. * b(1,iv)
  100 continue
      lp = inc
      lm = ninc - lp
      if (lp .ge. lm) go to 2
  1   do 200 iv = 1, lenv
         rp = a(lp+1,iv)
         rm = a(lm+1,iv)
         ip = b(lp+1,iv)
         im = b(lm+1,iv)
         a(lp+1,iv) = rm + rp
         b(lm+1,iv) = rm - rp
         b(lp+1,iv) = ip + im
         a(lm+1,iv) = ip - im
  200 continue
      lp = lp + inc
      lm = ninc - lp
      if (lp .lt. lm) go to 1
  2   if (lp .gt. ninc) return
      do 300 iv = 1, lenv
         a(lp+1,iv) = 2. * a(lp+1,iv)
         b(lp+1,iv) = 2. * b(lp+1,iv)
  300 continue

      return
      end
c======================================================================
      subroutine vrpfti2 (a, b, n, incp, lenv, lfd)
      integer(ISZ):: n, incp, lenv, lfd
      real(kind=8):: a(lfd,1), b(lfd,1)

c Vectorized version of Langdon's RPFTI2
c Call before VCPFT(-1)
c Inverts transform vcpft(+1),vrpft2 except arrays have been
c   multiplied by 2*n.
c Written by Alex Friedman, LLNL, April 1989.

      integer(ISZ):: inc,ninc,lp,lm,iv
      real(kind=8):: ca,sb,cb,sa

      inc = incp
      ninc = n * inc
      lp = inc
      lm = ninc - lp
      if (lp .ge. lm) return
  3   do 100 iv = 1, lenv
         ca = a(lp+1,iv)
         sb = b(lm+1,iv)
         cb = b(lp+1,iv)
         sa = a(lm+1,iv)
         a(lp+1,iv) = ca - sb
         a(lm+1,iv) = ca + sb
         b(lp+1,iv) = cb + sa
         b(lm+1,iv) = cb - sa
  100 continue
      lp = lp + inc
      lm = ninc - lp
      if (lp .lt. lm) go to 3

      return
      end
c======================================================================
      subroutine vcpft(r,i,n,incp,signp,lenv,lfd)
      integer(ISZ):: n,signp,lenv,lfd
      real(kind=8):: r(lfd,1), i(lfd,1)
c----------------------------------------------------------------------
c A vectorized version of CPFT (Langdon's).  Written by Dale Nielsen.
c----------------------------------------------------------------------
c  r      real part of data vector.
c  i      imag part of data vector.
c  n      number of elements (=1,2,4,8...131072).
c  inc    spacing in memory of data (usually 1, but see below).
c  sign   its sign will be sign of argument in transform exponential.
c  lenv   the length of the vector to be transformed
c  lfd    length of first dimension in data arrays
c
c    on entry arrays r and i contain the sequence to be transformed.
c  on exit they contain the transform. input and output sequences are
c  both in natural order (i.e. not bit-reversed scrambled).
c
c    a call to cpft with sign=+1, followed by another call with the
c  first 4 parameters the same and sign=-1, will leave r and i with
c  their original contents times n. the same is true if first sign=-1,
c  and next sign=+1.
c----------------------------------------------------------------------

      integer(ISZ):: span, rc,is,inc,ninc,it,k1,k0,iv,n1,n2,ij,ji,incp
      real(kind=8):: sines(17), r0,r1,i0,i1,t,sgn,c,s,c0,s0
      save sines
      data sines(1)/0./

      if ( sines(1).ne.1. ) then
c--------------------------------------------- calculate table of sines
        sines(1) = 1.
        t = atan(1.)
        do 1 is = 2,17
         sines(is) = sin(t)
         t = t/2.
1       continue
      endif

      if ( n.eq.1 ) return

c  set up various indices.

      inc  = incp
      sgn  = signp
      ninc = n * inc
      span = ninc
      it = n/2
      do 3 is=1,17
        if ( it.eq.1 ) go to 12
        it = it/2
3     continue

c***********************************************************************
c    there are 2 inner loops which run over the n/(2*span) replications
c  of transforms of length (2*span).
c  one loop is for arbitrary
c  rotation factor angle. the other takes care of the special case in
c  which the angle is zero so that no complex multiplication is needed.
c  this is more efficient than testing and branching inside the inner
c  loop, as is often done. the other special case in which no complex
c  multiply is needed is angle=pi (i.e. factor=i); this is not handled
c  specially. these measures are most helpful for small n.
c
c    the organization of the recursion is that of sande (ref. (3),
c  pp. 566-568). that is, the data is in normal order to start and
c  scrambled afterward, and the exponential rotation ('twiddle') factor
c  angles are used in ascending order during each recursion level.
c  all the sines and cosines needed are generated from a short table
c  using a stable multiple-angle recursion (ref. (1), p651 and ref. (2),
c  pp. 179-180). this method is economical in storage and time, and
c  yields accuracy comparable to good library sin-cos routines.
c  angles between 0 and pi are needed. the recursion is used for
c  angles up to pi/2; larger angles are obtained by reflection in the
c  imaginary axis (angle:=pi-angle). these pairs of angles are used
c  one right after the other.
c
c    for simplicity, commentary below applies to inc=1 case.
c***********************************************************************

10    t = s+(s0*c-c0*s)
      c = c-(c0*c+s0*s)
      s = t

c----------------------------------------------------- replication loop
11    k1 = k0+span

      do 5 iv = 1,lenv
        r0 = r(k0+1,iv)
        r1 = r(k1+1,iv)
        i0 = i(k0+1,iv)
        i1 = i(k1+1,iv)
        r(k0+1,iv) = r0+r1
        i(k0+1,iv) = i0+i1
        r0 = r0-r1
        i0 = i0-i1
        r(k1+1,iv) = c*r0-s*i0
        i(k1+1,iv) = s*r0+c*i0
5     continue

      k0 = k1+span
      if ( k0.lt.ninc ) go to 11
      k1 = k0-ninc
      c = -c
      k0 = span-k1
      if ( k1.lt.k0 ) go to 11
      k0 = k0+inc
      if ( k0.lt.k1 ) go to 10

c  recursion to next level.
12    continue
      span = span/2
      k0 = 0

c  angle = 0 loop.
13    k1 = k0+span

      do 6 iv=1,lenv
        r0 = r(k0+1,iv)
        r1 = r(k1+1,iv)
        i0 = i(k0+1,iv)
        i1 = i(k1+1,iv)
        r(k0+1,iv) = r0+r1
        i(k0+1,iv) = i0+i1
        r(k1+1,iv) = r0-r1
        i(k1+1,iv) = i0-i1
6     continue

      k0 = k1+span
      if ( k0.lt.ninc ) go to 13

      if ( span.eq.inc ) then
c------------------------------------------------------------- finished
        go to 20
      else
c---------------------------------------------- prepare non-zero angles
        c0 = 2.*sines(is)**2
        is = is-1
        s = sign( sines(is),sgn )
        s0 = s
        c = 1.-c0
        k0 = inc
        go to 11
      endif

c----------------------------------------------------------------------
c    arrays r and i now contain transform, but stored in 'reverse-
c    binary' order. the re-ordering is done by pair exchanges.
c
c    once again, commentary applies to inc=1 case.
c  indices are:
c    ij:=0,1,2...n/2-1 ( a simple counter).
c    ji:=reversal of ij.
c    rc:=reversal of 0,2,4...n/2 (incremented n/4 times).
c  rc is incremented thusly: starting with the next-to-leftmost bit,
c  change each bit up to and including first 0. (the actual coding is
c  done so as to work for any inc>0 with equal efficiency.)
c    for all exchanges ij fits one of these cases:
c      (1) 1st and last bits are 0 (ij,ji even and <n/2), and ij<=ji.
c      (2) one's complement of case (1) (both odd and >n/2).
c      (3) 1st bit 0, last bit 1 (ij odd and <n/2, ji>n/2).
c    the code from label even down to odd is entered with ij even and
c  <=ji. first time thru the complements are done -case (2). second
c  time thru gets case (1). thus a pair of elements both in the first
c  half of the sequence, and another pair in the 2nd half, are
c  exchanged. the condition ij<ji prevents a pair from being exchanged
c  twice.
c    the code from label odd down to increv does case (3).
c***********************************************************************

20    n1 = ninc-inc
      n2 = ninc/2
      rc = 0
      ij = 0
      ji = 0
      if ( n2.eq.inc ) return
      go to 22

c----------------------------------------------------------------- even
21    ij = n1-ij
      ji = n1-ji

      do 7 iv  =  1,lenv
        t = r(ij+1,iv)
        r(ij+1,iv) = r(ji+1,iv)
        r(ji+1,iv) = t
        t = i(ij+1,iv)
        i(ij+1,iv) = i(ji+1,iv)
        i(ji+1,iv) = t
7     continue

      if ( ij.gt.n2 ) go to 21

c------------------------------------------------------------------ odd
22    ij = ij+inc
      ji = ji+n2

      do 8 iv = 1,lenv
        t = r(ij+1,iv)
        r(ij+1,iv) = r(ji+1,iv)
        r(ji+1,iv) = t
        t = i(ij+1,iv)
        i(ij+1,iv) = i(ji+1,iv)
        i(ji+1,iv) = t
8     continue
      it = n2

c--------------------------------------- increment the reversed counter
23    it = it/2
      rc = rc-it
      if ( rc.ge.0 ) go to 23
      rc = rc+2*it
      ji = rc
      ij = ij+inc

      if ( ij.lt.ji ) go to 21
      if ( ij.lt.n2 ) go to 22

      return
      end

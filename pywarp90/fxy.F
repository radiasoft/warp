#include "top.h"
c=======================================================================
c@(#) File FXY.M, version $Revision: 3.0 $, $Date: 2001/04/11 17:32:10 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  written by David P. Grote
c  It is the source file for the package FXY of the PIC code WARPxy,
c  but it may be useful by itself.  It contains:
c
c    Capacity matrix solver which works in xy space
c
c=======================================================================
      subroutine fxyinit
      use FXYversion

c  Called at first reference to package (not nec. a "run" etc.).


      call fxyvers (STDOUT)

      return
      end
c=============================================================================
      subroutine fxyvers (iout)
      use FXYversion
      integer     :: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"Fieldsolver FXY",versfxy)
      return
      end
c=============================================================================
      subroutine fxygen()

      return
      end
c=======================================================================
      subroutine fxyfin()
      use FXYversion

c  Deallocates dynamic storage used by test driver


      call gfree ("FXYvars")

      return
      end
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c                                                                          =
c Capacity matrix solver                                                   =
c                                                                          =
c===========================================================================
c=======================================================================
      subroutine cmsetxy(phi,kxsq,kysq,attx,atty,filt,xlen,ylen,nx,ny,dx,dy,
     &                   xmmin,ymmin,scrtch1,scrtch2,l2symtry,l4symtry)
      use CapMatxy
      integer     :: nx,ny
      real(kind=8):: phi(0:nx,0:ny)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1)
      real(kind=8):: filt(5,2),scrtch1(0:nx,0:ny),scrtch2(*)
      real(kind=8):: xlen,ylen,dx,dy,xmmin,ymmin
      logical     :: l2symtry,l4symtry

c Sets up capacity matrix for xy fieldsolver. Note that the capacity matrix
c is calculated using field solves which do not do any filtering.
c (Experience has shown that when filtering is used, the inverse capacity
c matrix is very ill-conditioned and so the capacity matrix becomes garbled.)


      real(kind=8):: wx,wy,tt1,tt2,cx,cy
      integer     :: ix,iy,i,j,info

c     --- If the number of conductor points is zero, print a warning
c     --- and return without doing anything
      if (ncxy == 0) then
        call remark("NOTICE: the capacity matrix field solver is turned on")
        call remark("        but no conductors have been defined.")
        return
      endif

c     --- Make a copy of phi since it will be scrambled below.
      call copyarry(phi,scrtch1,(nx+1)*(ny+1))

c find capacity matrix

c     --- find inverse capacity matrix
      do i=1,ncxy

c       --- start with zero phi
c       call zeroarry(phi,(nx+1)*(ny+1))
        phi = 0.

c       --- set phi with 1 Coulomb at conductor points (charge is scaled by
c       --- 2 or 4 at a symmetric boundary)
        ix = int((xcond(i)-xmmin)/dx)
        wx =     (xcond(i)-xmmin)/dx  - ix
        iy = int((ycond(i)-ymmin)/dy)
        wy =     (ycond(i)-ymmin)/dy  - iy
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy .eq. 0) cy = 2.
        if (l4symtry .and. ix .eq. 0) cx = 2.
        if (l4symtry .and. iy .eq. 0) cy = 2.
        phi(ix  ,iy  ) = (1.-wx)*(1.-wy)*cx*cy
        phi(ix+1,iy  ) =     wx *(1.-wy)*cy
        phi(ix  ,iy+1) = (1.-wx)*    wy *cx
        phi(ix+1,iy+1) =     wx *    wy

c       --- Solve for fields, piece-by-piece so that the filtering can be
c       --- skipped.
        call vpois2d(2,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)
        call vpois2d(4,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)
        call vpois2d(5,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

c       --- fill matrix with phi
        do j=1,ncxy
          ix = int((xcond(j)-xmmin)/dx)
          wx =     (xcond(j)-xmmin)/dx  - ix
          iy = int((ycond(j)-ymmin)/dy)
          wy =     (ycond(j)-ymmin)/dy  - iy
          cmatxy(j,i) = phi(ix  ,iy  )*(1.-wx)*(1.-wy) +
     &                  phi(ix+1,iy  )*    wx *(1.-wy) +
     &                  phi(ix  ,iy+1)*(1.-wx)*    wy  +
     &                  phi(ix+1,iy+1)*    wx *    wy
        enddo

      enddo

c     --- invert to get capacity matrix (checking for singular matrix)

c     --- These routines don't seem to give as good a result as those below.
c     call ssifa(cmatxy,ncxy,ncxy,kpvtxy,info)
c     if (info .ne. 0) then
c       print*,"ERROR: Capacity matrix is singular - it is likely that"
c       print*,"two or more conductor points are within the same grid cell."
c       call kaboom(0)
c     endif
c     call ssidi(cmatxy,ncxy,ncxy,kpvtxy,tt1,tt2,scrtch2,001)

c     --- These routines are from LAPACK.
#if WORDSIZE == 64
      call ssytrf("U",ncxy,cmatxy,ncxy,kpvtxy,scrtch2,(1+nx)*(1+ny),info)
#else
      call dsytrf("U",ncxy,cmatxy,ncxy,kpvtxy,scrtch2,(1+nx)*(1+ny),info)
#endif
      if (info .ne. 0) then
        print*,"ERROR: Capacity matrix is singular - it is likely that"
        print*,"two or more conductor points are within the same grid cell."
        print*,info
        call kaboom(0)
      endif
c     --- pcond is passed as workspace.
#if WORDSIZE == 64
      call ssytri("U",ncxy,cmatxy,ncxy,kpvtxy,pcond,info)
#else
      call dsytri("U",ncxy,cmatxy,ncxy,kpvtxy,pcond,info)
#endif

c     --- Fill in lower half (which is the same as the upper half)
      do j=2,ncxy
        do i=1,j-1
          cmatxy(j,i) = cmatxy(i,j)
        enddo
      enddo

c     --- Return the original saved phi.
      call copyarry(scrtch1,phi,(nx+1)*(ny+1))

      return
      end
c============================================================================== 
      subroutine capmatxyf(iwhich,phi,kxsq,kysq,attx,atty,filt,xlen,ylen,nx,ny,
     &                    dx,dy,xmmin,ymmin,scrtch1,scrtch2,l2symtry,l4symtry)
      use CapMatxy
      integer     :: iwhich,nx,ny
      real(kind=8):: phi(0:nx,0:ny)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1)
      real(kind=8):: filt(5,2),scrtch1(0:nx,0:ny),scrtch2(*)
      real(kind=8):: xlen,ylen,dx,dy,xmmin,ymmin
      logical     :: l2symtry,l4symtry

c Applies the capacity matrix for xy code electrostatic field solver.
c Multiplies phi by cap to get induced rho.
c Assumes that beam's rho has been set, and copied into phi.
c When filtering is being done, care is taken so that only the beam's charge
c density is filtered and not the induced charge on the conductor.
c Note that fieldsolve is set up now so that the rho array is not used, since
c phi (which holds rho initially) is saved after the first FFT. This allows
c the capacity matrix to be used in bends.


      real(kind=8):: wx,wy,cx,cy
      integer     :: ix,iy,i,j
      logical     :: lattenuate

c     --- Initialize the arrays for poisson solve and the capacity matrix.
      if (iwhich == 1 .or. iwhich == 0) then

        call vpois2d(1,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

        call cmsetxy(phi,kxsq,kysq,attx,atty,filt,xlen,ylen,nx,ny,dx,dy,
     &               xmmin,ymmin,scrtch1,scrtch2,l2symtry,l4symtry)

      elseif (iwhich > 1) then

c       --- Make call to do the specialized action requested.
        call vpois2d(iwhich,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

      endif

c     --- If a full field solve was not desired, return.
      if (iwhich > 0) return

c     --- Do the first field solve (note that rho must have been copied into
c     --- phi for this to work)
c     --- This is done piece-by-piece so that when filtering is done, the
c     --- filtered source can be saved.

c     --- First transform phi (actually the source).
      call vpois2d(2,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

c     --- Check if filtering is turned on. If so, filter.
      lattenuate = .false.
      do j=1,2
        do i=1,5
          if (filt(i,j) .ne. 0) lattenuate = .true.
        enddo
      enddo
      if (lattenuate) then
        call vpois2d(3,phi,phi,kxsq,kysq,attx,atty,filt,
     &               xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)
      endif

c     --- Now, save the transformed and possibly filtered source.
      call copyarry(phi,scrtch1,(nx+1)*(ny+1))

c     --- Finish the first field solve: divide by k-squared and untransform.
      call vpois2d(4,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)
      call vpois2d(5,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

c     --- Extract phi error from FFT solution and zero qcond
      do i=1,ncxy
        ix = int((xcond(i)-xmmin)/dx)
        wx =     (xcond(i)-xmmin)/dx  - ix
        iy = int((ycond(i)-ymmin)/dy)
        wy =     (ycond(i)-ymmin)/dy  - iy
        pcond(i) = vcond(i) - (phi(ix  ,iy  )*(1.-wx)*(1.-wy) +
     &                         phi(ix+1,iy  )*    wx *(1.-wy) +
     &                         phi(ix  ,iy+1)*(1.-wx)*    wy  +
     &                         phi(ix+1,iy+1)*    wx *    wy)
        qcond(i) = 0.
      enddo

c     --- Multiply pcond by capacity matrix to get induced charge
      do j=1,ncxy
        do i=1,ncxy
          qcond(i) = qcond(i) + pcond(j)*cmatxy(i,j)
        enddo
      enddo

c     --- Zero out phi (the source is added back in later).
c     call zeroarry(phi,(nx+1)*(ny+1))
      phi = 0.

c     --- Deposit induced charge onto grid
      do i=1,ncxy
        ix = int((xcond(i)-xmmin)/dx)
        wx =     (xcond(i)-xmmin)/dx  - ix
        iy = int((ycond(i)-ymmin)/dy)
        wy =     (ycond(i)-ymmin)/dy  - iy
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy .eq. 0) cy = 2.
        if (l4symtry .and. ix .eq. 0) cx = 2.
        if (l4symtry .and. iy .eq. 0) cy = 2.
        phi(ix  ,iy  ) = phi(ix  ,iy  ) + qcond(i)*(1.-wx)*(1.-wy)*cx*cy
        phi(ix+1,iy  ) = phi(ix+1,iy  ) + qcond(i)*    wx *(1.-wy)*cy
        phi(ix  ,iy+1) = phi(ix  ,iy+1) + qcond(i)*(1.-wx)*    wy *cx
        phi(ix+1,iy+1) = phi(ix+1,iy+1) + qcond(i)*    wx *    wy
      enddo

c     --- Redo the field solve, again piece-by-piece so that the already
c     --- transformed and filtered source may be added and so that the
c     --- filtering can be skipped. Note that the initial transform is
c     --- only applied to the induced charge.
      call vpois2d(2,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

c     --- Sum the induced and transformed beam charge.
      do i=0,nx
        do j=0,ny
          phi(i,j) = scrtch1(i,j) + phi(i,j)
        enddo
      enddo

c     --- Complete the field solve.
      call vpois2d(4,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)
      call vpois2d(5,phi,phi,kxsq,kysq,attx,atty,filt,
     &             xlen,ylen,nx,ny,scrtch2,0,l2symtry,l4symtry)

      return
      end
c============================================================================== 

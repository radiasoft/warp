#include "top.h"
c=======================================================================
c@(#) File F3D.M, version $Revision: 3.136 $, $Date: 2009/01/05 18:27:25 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  written by Alex Friedman and David P. Grote
c  It is the source file for the package F3D of the PIC code WARP3d,
c  but it may be useful by itself.  It contains:
c
c     1) Dale Nielsen's vectorized version of Langdon's cpft, VCPFT
c     2) Vectorized two at a time real(kind=8):: transforms: VRPFT2, VRPFTI2
c     3) Vectorized sine transforms, in X and Y: VSFTX, VSFTY
c     4) A vectorized 3d sine-sine-periodic Poisson solver: VPOIS3D
c     5) A Python test driver (invoke via PACKAGE F3D;GENERATE)
c        The arrays and scalars used are local to package F3D; some have the
c        same names as those in the associated particle code package W6,
c        but they are <<not>> the same variables.  After calling this driver
c        and examining the results, the user should type FINISH to deallocate.
c
c  Call chain for the vector 3d s-s-p fieldsolve is (roughly):
c  VPOIS3D - VSFTX,VSFTY,VCPFT,VRPFT2,VRPFTI2
c
c     6) Capacity matrix solver which works in kz space
c     7) Capacity matrix solver which works in real(kind=8):: space
c     8) PSOR field solver including internal conductors and sub-grid scale
c        placement of boudaries
c
c=======================================================================
      subroutine f3dinit
      use F3Dversion

c  Called at first reference to package (not nec. a "run" etc.).


      call f3dvers (STDOUT)

      return
      end
c=============================================================================
      subroutine f3dvers (iout)
      use F3Dversion
      integer(ISZ):: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"Fieldsolver F3D",versf3d)
      return
      end
c=============================================================================
      subroutine f3dgen ()
      use Constant
      use F3Dvars

c  Invoked by the GENERATE command, it tests the fieldsolver.

      real(kind=8):: dx,dy,dz,cx,cy,cz,aa,sumsq,sqm,errmax,time,delsq,err
      integer(ISZ):: nx2,ny2,nz2,ix,iy,iz,kx,ky,kz,itest
      real(kind=8):: wranf,wtimeoff

      nx = 64
      ny = 64
      nz = 64
      call gallot("F3Dvars",0)

c  Check for improper input; we're using power-of-2 fft's

      if (  ( nx .ne. 2**int( .5 + log(1.*nx) / log(2.) ) )
     & .or. ( ny .ne. 2**int( .5 + log(1.*ny) / log(2.) ) )
     & .or. ( nz .ne. 2**int( .5 + log(1.*nz) / log(2.) ) ) ) then
         print*," ---> Bad nx, ny, or nz ... must be powers of 2"
         print*,"      nx =",nx,"  ny =",ny,"  nz =",nz
         return
      endif

c  Echo system size to the standard output (generally the user's terminal)

      print*," ***  testing 3d sine-sine-periodic transform  ***"
      print*,"      nx =",nx,"  ny =",ny,"  nz =",nz

c  Set constants: pi, grid spacing

      dx = lx / nx
      dy = ly / ny
      dz = lz / nz
      nx2 = nx / 2
      ny2 = ny / 2
      nz2 = nz / 2

c  Check 3d forward transform on known modes

      cx =    pi/nx
      cy =    pi/ny
      cz = 2.*pi/nz
      do ix = 1, nx-1
        do iy = 1, ny-1
          do iz = 0, nz-1
            b(ix,iy,iz)=0.
            b(ix,iy,iz)=b(ix,iy,iz)+   sin(cx*ix*1.)*sin(cy*iy)*cos(cz*iz*3.)
            b(ix,iy,iz)=b(ix,iy,iz)+2.*sin(cx*ix*2.)*sin(cy*iy*3.)
            b(ix,iy,iz)=b(ix,iy,iz)+3.*sin(cx*ix*3.)*sin(cy*iy)*sin(cz*iz*1.)
            b(ix,iy,iz)=b(ix,iy,iz)+4.*sin(cx*ix*1.)*sin(cy*iy)*(-1.)**iz
          enddo
        enddo
      enddo
c     --- Call the solver asking for just a forward transform
      call vpois3d (2,b,b,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &              lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,0,.false.,.false.,
     &              2,2,0)
c     --- Subtract off correct mode amplitudes.
      b(1,1,nz2-3) = b(1,1,nz2-3) - lx*ly*lz * .125
      b(2,3,nz2  ) = b(2,3,nz2  ) - lx*ly*lz * .5
      b(3,1,nz2+1) = b(3,1,nz2+1) - lx*ly*lz * .375
      b(1,1,0    ) = b(1,1,0    ) - lx*ly*lz * 1.
      call remark("Incorrect modes, if any: ")
      do ix = 1, nx-1
        do iy = 1, ny-1
          do iz = 0, nz-1
            kx = ix
            ky = iy
            kz = iz-nz2
            aa = b(ix,iy,iz)
            if( abs(aa) > 1.e-6*lx*ly*lz ) print*,kx, ky, kz, aa
          enddo
        enddo
      enddo

c  check 3d transform and inverse on 1-point-source and random data.

      do itest = 0,1
         sumsq=0.
         do ix = 0, nx-1
           do iy = 0, ny-1
             do iz = 0, nz
               b(ix,iy,iz) = 0.
               if (itest==1 .and. ix > 0 .and. iy > 0) then
                  if (iz < nz) b(ix,iy,iz) = wranf() - .5
               elseif (ix==nx2 .and. iy==ny2 .and. iz==nz2) then
                  b(ix,iy,iz) = 1.
               endif
               bsav(ix,iy,iz) = b(ix,iy,iz)
               sumsq = sumsq + dx*dy*dz*b(ix,iy,iz)**2
             enddo
           enddo
         enddo
c        --- Call the solver asking for just a forward transform
         call vpois3d (2,b,b,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,0,.false.,.false.,
     &                2,2,0)
c        --- sum of squares from fourier mode amplitudes.
         sqm = 0.
         do ix = 1, nx-1
         do iy = 1, ny-1
            sqm = sqm + ( b(ix,iy,nz2)**2 + b(ix,iy,0)**2 )
            do iz = 1, nz2-1
               sqm = sqm + 2. * ( b(ix,iy,nz2+iz)**2 + b(ix,iy,nz2-iz)**2 )
            enddo
         enddo
         enddo
         sqm = sqm * 4. / (lx * ly * lz)
         err = sumsq - sqm
         print*,"Sums of squares: sumsq =",sumsq," sqm =",sqm," err =",err
c        --- Call the solver asking for just an inverse transform
         call vpois3d (5,b,b,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,0,.false.,.false.,
     &                2,2,0)
         errmax = 0.
         do iz = 0,nz-1
           do ix = 1,nx-1
             do iy = 1,ny-1
              if (abs(b(ix,iy,iz)-bsav(ix,iy,iz)) > 1.e-5) then
                print*,ix,iy,iz,bsav(ix,iy,iz),b(ix,iy,iz)
              endif
              errmax = max( errmax,abs( b(ix,iy,iz) - bsav(ix,iy,iz) ) )
           enddo
           enddo
         enddo
         print*,"Maximum error in inverse =",errmax
      enddo

c  Check the fieldsolver on random data

c     --- Dummy call to get k-squared set up and show it works more than once.
      call vpois3d (0,b,b,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &              lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,0,.false.,.false.,
     &              2,2,0)
      do ix = 0, nx-1
      do iy = 0, ny-1
      do iz = 0, nz
         b(ix,iy,iz) = 0.
         if (ix > 0.and.iy > 0.and.iz < nz) b(ix,iy,iz)=wranf()-.5
c*****   if (ix==nx2.and.iy==ny2.and.iz==nz2) b(ix,iy,iz) = 1.
         bsav(ix,iy,iz) = b(ix,iy,iz)
      enddo
      enddo
      enddo
      call wtimeon
c     --- call the fieldsolver, asserting the k-squared arrays are already set
      call vpois3d (-1,b,b,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &              lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,0,.false.,.false.,
     &              2,2,0)
      time = wtimeoff()
c     --- check using 7 point poisson difference equation.
      errmax = 0.
      do ix = 1, nx-1
      do iy = 1, ny-1
         delsq=(b(ix+1,iy  ,0)-2.*b(ix,iy,0)+b(ix-1,iy  ,0   ))/dx**2
     &        +(b(ix  ,iy+1,0)-2.*b(ix,iy,0)+b(ix  ,iy-1,0   ))/dy**2
     &        +(b(ix  ,iy  ,1)-2.*b(ix,iy,0)+b(ix  ,iy  ,nz-1))/dz**2
         bsav(ix,iy,0) = bsav(ix,iy,0) + bsav(ix,iy,nz)
     &                   + delsq*eps0
         errmax = max(errmax, abs(bsav(ix,iy,0)))
         do iz = 1, nz-1
            delsq=( b(ix+1,iy,iz)-2.*b(ix,iy,iz)+b(ix-1,iy,iz) )/dx**2
     &           +( b(ix,iy+1,iz)-2.*b(ix,iy,iz)+b(ix,iy-1,iz) )/dy**2
     &           +( b(ix,iy,iz+1)-2.*b(ix,iy,iz)+b(ix,iy,iz-1) )/dz**2
            bsav(ix,iy,iz) = bsav(ix,iy,iz) + delsq*eps0
            errmax = max(errmax, abs(bsav(ix,iy,iz)))
         enddo
      enddo
      enddo
      print*,"poisson time=",time," milliseconds"
      print*,"check by 7 point formula; residual max=",errmax

      return
      end
c=======================================================================
      subroutine f3dfin ()
      use F3Dvars

c  Deallocates dynamic storage used by test driver


      call gfree ("F3Dvars")

      return
      end
c=======================================================================
      subroutine vp3x (iwhich)
      use F3Dvars
      integer(ISZ):: iwhich

c  Python interface to VPOIS3D, using variables from the database for
c  package F3D.  The user's application package should contain a
c  similar subroutine for convenience, using its own database variables.
c  To avoid confusion, it should have a different name, say "vp3"


      call vpois3d (iwhich,b,b,kxsq,kysq,kzsq,attx,atty,attz,
     &              filt,lx,ly,lz,nx,ny,nz,nz,work,work2d,zwork,
     &              0,.false.,.false.,2,0,0)

      return
      end
c=======================================================================
c=============================================================================
c=============================================================================
      subroutine vpois3d(iwhich,a,ak,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                   lx,ly,lzfull,nx,ny,nzlocal,nz,w,xywork,zwork,
     &                   ibc,l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Constant
      use Timers
      use GlobalVars
#ifdef MPIPARALLEL
      use Parallel
      use Transpose_work_space
#endif
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz,ibc
      real(kind=8)::  a(-1:nx+1, -1:ny+1, 0:nzlocal)
      real(kind=8):: ak(-1:nx+1, -1:ny+1, -nzlocal/2:nzlocal/2)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz), filt(5,3)
      real(kind=8):: lx, ly, lzfull
      real(kind=8):: w(*), xywork(0:nx,0:ny), zwork(2,0:nx,0:nz)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy
      
c-----------------------------------------------------------------------------
c               Vectorized 3d Poisson solver:
c               Sine transform in x and y, periodic in z
c               Alex Friedman, LLNL, May 1989.
c
c          Returned value: 0 if ok, 1 if error
c
c          IWHICH = ... -1: full fieldsolve; assumes kxsq etc. and attx etc.
c                           have been set already.  This is equivalent to a
c                           sequence of calls with iwhich = 2,3,4,5.
c                        0: full fieldsolve; kxsq etc. will be set herein.
c                           This is equivalent to a sequence of calls with
c                           iwhich = 1,2,3,4,5. (THIS IS THE SIMPLEST USAGE.)
c                        1: set kxsq etc. and attx etc., return.
c                        2: forward transform (no k^2 factor), return.
c                        3: apply k-space filter using attx, etc., return.
c                        4: apply inverse-k^2 factor (rho to phi), return.
c                        5: inverse transform (no k^2 factor), return.
c                        8: apply reciprocal of filter in k-space, return.
c                        9: apply k^2 factor (turns phi into rho), return.
c                        10: forward transform in z (no k^2 factor), return.
c                        11: inverse transform in z (no k^2 factor), return.
c                        12: forward transform in x&y (no k^2 factor), return.
c                        13: inverse transform in x&y (no k^2 factor), return.
c                        14: tridiagonal solve along z, return
c              A, AK ... Array to transform (pass it to this routine twice)
c                        Use of both a and ak simplifies subscripting of the
c                        Fourier coefficients.
c                        See array declarations above these comments.
c   KXSQ, KYSQ, KZSQ ... Discrete analogs to eps0*kx^2 etc.
c                        See array declarations below these comments.
c   ATTX, ATTY, ATTZ ... Attenuation factors for each mode number in x, y, z.
c                        See array declarations below these comments.
c          FILT(i,j) ... Spatial filtering coefficients; array should be
c                        dimensioned (5,3) in calling routine.
c                        Here j is the axis (1=x, 2=y, 3=z),
c                        and there are four coefficients for each axis:
c                        filt(1,j): coefficient a1, of sin^2 term (boosting)
c                        filt(2,j): coefficient a2, of tan^4 term (smoothing)
c                        filt(3,j): cutoff k dx, k dy, or k dz, e.g. pi/2
c                        filt(4,j): exponent N (sharpness of cutoff) e.g. 8.
c                        filt(5,j): W factor in binomial filter, e.g. 0.5
c                        (A value of .5 gives the familiar 1/4:1/2:1/4 filter)
c                        See Birdsall and Langdon, "Plasma Physics via
c                        Computer Simulation," Appendices B and C.
c               EPS0 ... Epsilon-0; use 1/4Pi for CGS, 1 for rationalized.
c         LX, LY, LZ ... System lengths
c         NX, NY, NZ ... Number of mesh points (mesh runs from 0 to NX, etc.)
c                        At present these must be powers of two.
c                  W ... Workspace, length at least  NX + NY - 4
c             XYWORK ... Transverse work array, dimensined (0:nx,0:ny)
c              ZWORK ... Longitudinal work array, dimensioned (2,0:nx,0:nz)
c                IBC ... Unused
c           L2SYMTRY ... Switch for transverse 2-fold symmetry
c           L4SYMTRY ... Switch for transverse 4-fold symmetry
c
c  Note that for x (y is similar), k dx = pi j/nx, j = 1, ..., nx-1.
c        The longest wavelength is 2 nx dx.
c        The j=nx (odd-even) term is absent because the boundary pts are zero.
c  For the (periodic) z direction, k dz = 2 pi j/nz, j = 0, ..., nz/2.
c        The longest wavelength (excepting the k=0 mode) is nz dz = lz.
c
c  This routine was inspired by Langdon's scalar sine-periodic solver "sppois"
c
c  The tranposed logical flag is only used in the parallel version, and only
c  for the options used with iwhich <= 0. If a transpose is needed, this
c  routine will always take the inverse as well so that the data returned
c  is never transposed. Also, it is always assumed that the input data
c  is not transposed.
c-----------------------------------------------------------------------------

      integer(ISZ):: icmx,icpx,icmy,icpy,icox,icoy,nz2,ikxmin,ikymin
      integer(ISZ):: ikx,iky,ikz,ix,iy,iz,izmin,izmax
      real(kind=8):: klast, kdxb2, kdyb2, kdzb2, kdx, kdy, kdz
      real(kind=8):: dx,dy,dz,pib2,norm
      integer(ISZ):: i,j
      logical(ISZ):: lattenuate,doperiodic,dodirichlet,doneumann
      logical(ISZ):: periodicsolve,dirichletsolve
#ifdef MPIPARALLEL
      integer(ISZ):: nx_tran,ny_tran,ikxm,ikym
      logical(ISZ):: transposed
#endif


#ifdef MPIPARALLEL
c  The input data is always assumed to be not transposed
      transposed = .false.
#endif

#ifdef ESSL

c  The ESSL routines themselves have no limit on the number of points.
c  Because of the way vpftz works, ny must be even (this is easily fixable).
c  Because of the way the cosqx and cosqy routines work, nx and ny must
c  be even (i.e. when symmetry is turned on).
      if (mod(ny,2) == 1) then
        call remark("NY must be even")
        return
      endif
      if (l4symtry .and. (mod(nx,2) == 1 .or. mod(ny,2) == 1)) then
        call remark("NX and NY must be even")
        return
      endif

#else

c  Error checking; so far, only possible error is non-power-of-2 nx,ny,nz
c  Note that nz may be a non-power-of-2 with the options 1, 12, 13, and 14.
c  This allows a non-power-of-2 for the tridiag in z field solver.

      if (  ( nx .ne. 2**int( .5 + log(1.*nx) / log(2.) ) )
     & .or. ( ny .ne. 2**int( .5 + log(1.*ny) / log(2.) ) )
     & .or.(( nz .ne. 2**int( .5 + log(1.*nz) / log(2.) ) ) .and.
     &      ( iwhich .ne. 1 .and. iwhich .ne. 12 .and. iwhich .ne. 13 .and.
     &        iwhich .ne. 14 ) ) ) then
         call kaboom("NX and NY must be powers of 2, NZ also unless using tridiag solver")
         return
      endif
#endif

c  Set pointers into workspace needed by x and y sine transforms

      icmx  = 1
      icpx  = icmx  + nx/2-1
      icmy  = icpx  + nx/2-1
      icpy  = icmy  + ny/2-1
      icox  = icpy  + ny/2-1
      icoy  = icox  + nx
c     inext = icoy  + ny     is first word beyond workspace

c  Set useful constants

      nz2 = nz / 2
      dx = lx / nx
      dy = ly / ny
      dz = lzfull / nz
      pib2 = .5 * pi

c     --- Setup convenient logicals
      doperiodic = (bound0 == periodic) .or. (boundnz == periodic)
      dodirichlet = (bound0 == dirichlet) .or. (boundnz == dirichlet)
      doneumann = (bound0 == neumann) .or. (boundnz == neumann)

      if (doneumann) then
        call kaboom("FFT solver does not support Neumann boundary conditions yet")
        return
      endif

c     --- These are used to select which type of longitudinal solver to do.
c     --- Periodic uses FFT's, Dirichlet used tridiagonal
      periodicsolve = doperiodic .and. (iwhich <= 0)
      dirichletsolve = dodirichlet .and. (iwhich <= 0)

      izmin = 0
      if (dodirichlet) then
c       --- This only applies to the tridiagonal solver, which is the only
c       --- solver here which can apply Dirichlet boundary conditions.
c       --- Sine the tridiag solve is done in k-space, the transforms must
c       --- be done over the full array.
c       --- This is only used in the sections 12 and 13 when doing the
c       --- transverse FFT's.
        izmax = nzlocal
      else
        izmax = nzlocal - 1
      endif

c  Minimum index of kxsq and kysq that is used
      ikxmin = 1
      if (l4symtry) ikxmin = 0
      ikymin = 1
      if (l2symtry .or. l4symtry) ikymin = 0

c ----------------------------------------------------------------------------
c  For poisson equation, kxsq(ikx) is a discrete
c  analog to eps0 kx^2, etc.  If the user has requested it,
c  we set up ksq arrays for a seven point scheme now; the coding for kxsq
c  etc. is arranged so that it will vectorize on "most" compilers.
c  Also, compute attenuation factors as functions of mode numbers in x, y, z.
c ----------------------------------------------------------------------------

      if (iwhich == 0 .or. iwhich == 1) then

c        --- compute x direction rho-to-phi coefficients ---
         do ikx = ikxmin, nx-1
            kxsq(ikx) = ikx
         enddo
         if (l4symtry) then
           do ikx = 0, nx-1
              kxsq(ikx) = 2./dx**2*(1.-cos(pib2*(2.*kxsq(ikx)+1.)/nx)) * eps0
           enddo
         else
           do ikx = 1, nx-1
              kxsq(ikx) = (2./dx*sin(pib2*kxsq(ikx)/nx))**2 * eps0
           enddo
         endif
c        --- compute y direction coefficients ---
         do iky = ikymin, ny-1
            kysq(iky) = iky
         enddo
         if (l2symtry .or. l4symtry) then
           do iky = 0, ny-1
              kysq(iky) = 2./dy**2*(1.-cos(pib2*(2.*kysq(iky)+1.)/ny)) * eps0
           enddo
         else
           do iky = 1, ny-1
              kysq(iky) = (2./dy*sin(pib2*kysq(iky)/ny))**2 * eps0
           enddo
         endif
c        --- compute z direction coefficients (note differences from x,y) ---
         do ikz = 0, nz2
            kzsq(nz2+ikz) = (2./dz*sin(pi*ikz/nz))**2 * eps0
            kzsq(nz2-ikz) = kzsq(nz2+ikz)
         enddo

c        --- spatial filtering in x; binomial filter, or unity
         do ikx = ikxmin, nx-1
            kdx = pi * ikx / nx
            attx(ikx) = (1.+2.*filt(5,1)*cos(kdx)) / (1.+2.*filt(5,1))
         enddo
         if (filt(1,1) .ne. 0. .or. filt(2,1) .ne. 0.) then
c           --- compute first form of spatial filtering
            do ikx = ikxmin, nx-1
               kdxb2 = pi * ikx / (2.*nx)
               attx(ikx) = attx(ikx) * (exp( filt(1,1)*(sin(kdxb2))**2
     &          - filt(2,1)*(tan(kdxb2))**4 ))**2
            enddo
         endif
         if (filt(3,1) .ne. 0.) then
c           --- compute second form of spatial filtering
            klast = filt(3,1) * nx / pi
            do ikx = ikxmin, nx-1
               attx(ikx) = attx(ikx) * exp(-(ikx/klast)**filt(4,1))
            enddo
         endif

c        --- spatial filtering in y; binomial filter, or unity
         do iky = ikymin, ny-1
            kdy = pi * iky / ny
            atty(iky) = (1.+2.*filt(5,2)*cos(kdy)) / (1.+2.*filt(5,2))
         enddo
         if (filt(1,2) .ne. 0. .or. filt(2,2) .ne. 0.) then
c           --- compute first form of spatial filtering
            do iky = ikymin, ny-1
               kdyb2 = pi * iky / (2.*ny)
               atty(iky) = atty(iky) * (exp( filt(1,2)*(sin(kdyb2))**2
     &          - filt(2,2)*(tan(kdyb2))**4 ))**2
            enddo
         endif
         if (filt(3,2) .ne. 0.) then
c           --- compute second form of spatial filtering
            klast = filt(3,2) * ny / pi
            do iky = ikymin, ny-1
               atty(iky) = atty(iky) * exp(-(iky/klast)**filt(4,2))
            enddo
         endif

c        --- spatial filtering in z; binomial filter, or unity
         do ikz = 0, nz2
           kdz = 2. * pi * ikz / nz
           attz(nz2+ikz) = (1.+2.*filt(5,3)*cos(kdz)) / (1.+2.*filt(5,3))
         enddo
         if (filt(1,3) .ne. 0. .or. filt(2,3) .ne. 0.) then
c          --- compute first form of spatial filtering
           do ikz = 0, nz2
              kdzb2 = pi * ikz / nz
              attz(nz2+ikz) = attz(nz2+ikz) * (exp( filt(1,3)*(sin(kdzb2))**2
     &         - filt(2,3)*(tan(kdzb2))**4 ))**2
           enddo
         endif
         if (filt(3,3) .ne. 0.) then
c          --- compute second form of spatial filtering
           klast = filt(3,3) * nz / (2. * pi)
           do ikz = 0, nz2
              attz(nz2+ikz) = attz(nz2+ikz) * exp(-(ikz/klast)**filt(4,3))
           enddo
         endif
         do ikz = 1, nz2
           attz(nz2-ikz) = attz(nz2+ikz)
         enddo

      endif

c ----------------------------------------------------------------------------
c  Set up the tables needed by the sine transforms (call them w/ ISETUP=1)
c ----------------------------------------------------------------------------

      if (iwhich <= 0  .or. iwhich == 2  .or. iwhich == 5 .or.
     &    iwhich == 10 .or. iwhich == 11 .or.
     &    iwhich == 12 .or. iwhich == 13) then

        if (l4symtry) then
          call cosqx (a, xywork, w(icox), nx, ny, 1,-1)
          call cosqy (a, xywork, w(icoy), nx, ny, 1,-1)
        elseif (l2symtry) then
          call vsftx (a, xywork, w(icpx), w(icmx), nx, ny, 1)
          call cosqy (a, xywork, w(icoy), nx, ny, 1,-1)
        else
          call vsftx (a, xywork, w(icpx), w(icmx), nx, ny, 1)
          call vsfty (a, xywork, w(icpy), w(icmy), nx, ny, 1)
        endif

      endif

!$OMP PARALLEL PRIVATE(xywork,zwork)

c ----------------------------------------------------------------------------
c  Do the forward transform in x and y
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 2 .or. iwhich == 12) then

c  take z-sequence of sine-sine transforms in the x-y plane.
c  We should multitask this loop, being careful about scratch
c  which is (nx+1)*(ny+1) per task.  For now, we use the plane
c  at nz for the single task's scratch, and so need no other.

        if  (periodicsolve .or. iwhich == 2) then
c         --- Include normalization factor for the z transform
          norm = .5*dz
        else
          norm = 1.
        endif

      if (l4symtry) then
c       factor of .5dx for cosqx and .5*dy for cosqy
        norm = .25*dx*dy*norm
!$OMP DO
        do iz = izmin,izmax
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call cosqx(a(-1,-1,iz),xywork,w(icox),nx,ny,0,-1)
           call cosqy(a(-1,-1,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
!$OMP END DO
      else if (l2symtry) then
c       factor of .5*dy for cosqy
        norm = .5*dx*dy*norm
!$OMP DO
        do iz = izmin,izmax
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(-1,-1,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call cosqy(a(-1,-1,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
!$OMP END DO
      else
        norm = dx*dy*norm
!$OMP DO
        do iz = izmin,izmax
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(-1,-1,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call vsfty(a(-1,-1,iz),xywork,w(icpy),w(icmy),nx,ny,0)
        enddo
!$OMP END DO
      endif

      endif

c ----------------------------------------------------------------------------
c  Do the forward transform in z
c ----------------------------------------------------------------------------

      if (periodicsolve .or. iwhich == 2 .or. iwhich == 10) then

c  Here we perform the periodic transforms in z, two y's at a time.
c  We vectorize over x.  No scratch space is needed.

        if (periodicsolve .or. iwhich == 2) then
c         --- Normalization factor include above in the transverse transform
          norm = 0.
        else
          norm = .5*dz
        endif

#ifdef MPIPARALLEL
c  For parallel code, transpose is used to redistribute phi so each axial
c  line is on a single processor.
        call warptranspose(nx,ny,nzlocal,a,nx_tran,ny_tran)
        transposed = .true.

        if (ny_tran == 0) then
          call kaboom("vpois3d: having nslaves > ny is not supported now")
c         --- If there is only one y line, pass nx as zero so that
c         --- transform done on 2 y lines actually transforms on 2 x's
c         --- next to each other.
          call vpftz(0,nx_tran,nz,phi_trnsps,norm,0,1,zwork)
        else
c         --- pass normal arguments
          call vpftz(nx_tran,ny_tran,nz,phi_trnsps,norm,1,0,zwork)
        endif

        if (iwhich == 2 .or. iwhich == 10 .or. ny <= nslaves) then
c         --- Transpose back if not doing the full field solve
c         --- Note that if ny <= nslaves, transform back as well since the
c         --- ordering of the x and y lines are modified and not compatible
c         --- with the attenuation and rhotophi calls. This is a lazy fix,
c         --- a better way would likely be to write appropriate copy loops
c         --- in vpftz to handle that case, since data is being copied
c         --- anyway, it wouldn't be much more work to reorganize it.
          call warptransposei(nx,ny,nzlocal,a,nx_tran,ny_tran)
          transposed = .false.
        endif
#else
c  Make call for serial code.
        call vpftz(nx,ny,nzlocal,a,norm,1,1,zwork)
#endif

      endif

c ----------------------------------------------------------------------------
c  Apply k-space filter; attenuate mode-by-mode
c ----------------------------------------------------------------------------

      if (periodicsolve .or. iwhich == 3) then

c       --- Check if attenuation is needed
c       --- This check is done everytime just in case filt is changed by
c       --- the user.
        lattenuate = .false.
        do j=1,3
          do i=1,5
            if (filt(i,j) .ne. 0) lattenuate = .true.
          enddo
        enddo

        if (lattenuate) then

#ifdef MPIPARALLEL
c  For parallel code, data may be arranged differently then serial version.

          if (transposed) then
c           --- Find location of transverse data relative to original dimensions
            ix = mod(my_index*nx_tran,nx)
            iy = int(my_index*ny/real(nzprocs))
c           --- Check if work includes boundary points.  If so, use
c           --- precalculate min's.
            if (ix == 0) then
              ikxm = ikxmin
            else
              ikxm = 0
            endif
            if (iy == 0) then
              ikym = ikymin
            else
              ikym = 0
            endif
c           --- attenuate transposed array
            call attenuate(nx_tran,ny_tran,nz,phi_trnsps,attx(ix),atty(iy),attz,
     &                     ikxm,ikym,1,0,1)
          else
            iz = my_index*nzlocal
            call attenuate(nx,ny,nzlocal,a,attx,atty,attz(iz),ikxmin,ikymin,1,1,1)
          endif
#else
c  Normal call for serial version.
          call attenuate(nx,ny,nzlocal,a,attx,atty,attz,ikxmin,ikymin,1,1,1)
#endif

        endif

      endif

c ----------------------------------------------------------------------------
c  Convert rhok to phik.
c ----------------------------------------------------------------------------

      if (periodicsolve .or. iwhich == 4) then

#ifdef MPIPARALLEL
c  For parallel code, data may be arranged differently then serial version.
        if (transposed) then
c         --- Find location of transverse data relative to original dimensions.
          ix = mod(my_index*nx_tran,nx)
          iy = int(my_index*ny/real(nzprocs))
c         --- Check if work includes boundary points.  If so, use precalculated
c         --- min's.
          if (ix == 0) then
            ikxm = ikxmin
          else
            ikxm = 0
          endif
          if (iy == 0) then
            ikym = ikymin
          else
            ikym = 0
          endif
c         --- make conversion with transposed array
          call rhotophi(nx_tran,ny_tran,nz,phi_trnsps,kxsq(ix),kysq(iy),kzsq,
     &                  ikxm,ikym,1,0,1)
        else
          iz = my_index*nzlocal
          call rhotophi(nx,ny,nzlocal,a,kxsq,kysq,kzsq(iz),ikxmin,ikymin,1,1,1)
        endif
#else
c  Call for Serial version
        call rhotophi(nx,ny,nzlocal,a,kxsq,kysq,kzsq,ikxmin,ikymin,1,1,1)
#endif

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform
c ----------------------------------------------------------------------------

      if (periodicsolve .or. iwhich == 5 .or. iwhich == 11) then

c  Inverse transform and shift in z.

        if (periodicsolve .or. iwhich == 5) then
c         --- Normalization factor is include below in the transverse transform
          norm = 0.
        else
          norm = 1./lzfull
        endif

#ifdef MPIPARALLEL
c  For parallel code, inverse transpose is used to redistribute phi into its
c  original format.
        if (.not. transposed) then
          call warptranspose(nx,ny,nzlocal,a,nx_tran,ny_tran)
          transposed = .true.
        endif
        if (ny_tran == 0) then
c         --- If there is only one y line, pass nx as zero so that
c         --- transform done on 2 y lines actually transforms on 2 x's
c         --- next to each other.
          call vpftzi(0,nx_tran,nz,phi_trnsps,norm,0,1,zwork)
        else
c         --- pass normal arguments
          call vpftzi(nx_tran,ny_tran,nz,phi_trnsps,norm,1,0,zwork)
        endif

c       --- Rearrange phi again back to original arrangement.
        call warptransposei(nx,ny,nzlocal,a,nx_tran,ny_tran)
        transposed = .false.

#else
c  Call for serial version
        call vpftzi(nx,ny,nzlocal,a,norm,1,1,zwork)
#endif

      endif

c ----------------------------------------------------------------------------
c Tridiagonal solve in z, including kxsq and kysq on left hand side
c ----------------------------------------------------------------------------

      if (dirichletsolve .or. iwhich == 14) then

        norm = dz**2/eps0

#ifdef MPIPARALLEL
c  For parallel code, transpose is used to redistribute phi so each axial
c  line is on a single processor.
c       --- Rearrange phi array so each axial line is on a single processor.
        call warptranspose(nx,ny,nzlocal,a,nx_tran,ny_tran)

c       --- Find location of transverse data relative to original dimensions.
        ix = mod(my_index*nx_tran,nx)
        iy = int(my_index*ny/real(nzprocs))
c       --- Check if work includes boundary points.  If so, use precalculated
c       --- min's.
        if (ix == 0) then
          ikxm = ikxmin
        else
          ikxm = 0
        endif
        if (iy == 0) then
          ikym = ikymin
        else
          ikym = 0
        endif
        call tridiag(nx_tran,ny_tran,nz,phi_trnsps,norm,
     &               kxsq(ix),kysq(iy),kzsq,ikxm,ikym,0,zwork)

c       --- Rearrange phi again back to original arrangement.
        call warptransposei(nx,ny,nzlocal,a,nx_tran,ny_tran)
#else
c  Serial version
        call tridiag(nx,ny,nzlocal,a,norm,kxsq,kysq,kzsq,ikxmin,ikymin,1,zwork)
#endif

      endif
c ----------------------------------------------------------------------------
c  Do the inverse transform in x and y only
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 5 .or. iwhich == 13) then

        if (periodicsolve .or. iwhich == 5) then
c         --- Include the normalization term for the z transform
          norm = 1./lzfull
        else
          norm = 1.
        endif

        if (l4symtry) then
c         factor of .25/lx for cosqx and .25/ly for cosqy
          norm = .0625/(lx*ly)*norm
!$OMP DO
          do iz = izmin,izmax
c            --- Invoke the inverse transform routines for the current z-slice
             call cosqx(a(-1,-1,iz),xywork,w(icox),nx,ny,0,1)
             call cosqy(a(-1,-1,iz),xywork,w(icoy),nx,ny,0,1)
c            --- Normalize the array
             do iy = 0, ny-1
               do ix = 0, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
          enddo
!$OMP END DO
        else if (l2symtry) then
c         factor of 2./lx for vsftx and .25/ly for cosqy
          norm = .5/(lx*ly)*norm
!$OMP DO
          do iz = izmin,izmax
c            --- Invoke the inverse transform routines for the current z-slice
             call vsftx ( a(-1,-1,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
             call cosqy(a(-1,-1,iz),xywork,w(icoy),nx,ny,0,1)
c            --- Normalize the array
             do iy = 0, ny-1
               do ix = 0, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
          enddo
!$OMP END DO
        else
          norm = 4. / (lx * ly)*norm
!$OMP DO
          do iz = izmin,izmax
c            --- Invoke the transform routines
             call vsfty ( a(-1,-1,iz), xywork, w(icpy), w(icmy), nx, ny, 0)
             call vsftx ( a(-1,-1,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
c            --- Normalize the array
             do iy = 1, ny-1
               do ix = 1, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
           enddo
!$OMP END DO
        endif

c  Enforce periodic bc's in z (now done in w3d)

      endif

c ----------------------------------------------------------------------------
c  The usual fieldsolve has been completed.  Utility operations follow.
c ----------------------------------------------------------------------------

c ----------------------------------------------------------------------------
c  Un-filter in k space
c ----------------------------------------------------------------------------

      if (iwhich == 8) then

#ifdef MPIPARALLEL
c  For parallel code, at this point, the data will never be transposed.
        iz = my_index*nzlocal
        call unattenuate(nx,ny,nzlocal,a,attx,atty,attz(iz),ikxmin,ikymin,1,1,1)
#else
c  Call for Serial version
        call unattenuate(nx,ny,nzlocal,a,attx,atty,attz,ikxmin,ikymin,1,1,1)
#endif

      endif

c ----------------------------------------------------------------------------
c  Multiply by k^2 factor, in contrast to the usual division.
c  This is useful if one wants to recreate rho from phi.
c ----------------------------------------------------------------------------

      if (iwhich == 9) then

#ifdef MPIPARALLEL
c  For parallel code, at this point, the data will never be transposed.
        iz = my_index*nzlocal
        call phitorho(nx,ny,nzlocal,a,kxsq,kysq,kzsq(iz),ikxmin,ikymin,1,1,1)
#else
c  Call for Serial version
#endif
        call phitorho(nx,ny,nzlocal,a,kxsq,kysq,kzsq,ikxmin,ikymin,1,1,1)

      endif


!$OMP END PARALLEL

c ----------------------------------------------------------------------------
c  End of VPOIS3D
c ----------------------------------------------------------------------------

      return
      end
c=============================================================================

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c=======================================================================
c The following two routines use the vcpft routine for the FFT's. That
c routine runs faster when the data is in transposed order and the following
c versions were kept.  Original unvectoized routine is in appendix to 
c  Birdsall/Langdon.
c The difference in speed is around 15% using the Intel compiler on Linux
c and about a factor of 2 on the IBM SP.
c=============================================================================
      subroutine vpftx(nx,ny,nz,a,norm,esx,esy,xywork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,xywork(2,0:nx,0:ny)

c Routine does the fft over x for entire array.  Outer loop is in z
c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c 
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iz = 0, nz-2, 2
c       --- shift k-space origin to middle and apply normalization factor
c       --- if non-zero.
        if (norm == 0) then
          do ix = 1, nx-esx, 2
            do iy = 0, ny-esy
              xywork(1,iy,ix-1) = +a(ix-1,iy,iz  )
              xywork(2,iy,ix-1) = +a(ix-1,iy,iz+1)
              xywork(1,iy,ix  ) = -a(ix  ,iy,iz  )
              xywork(2,iy,ix  ) = -a(ix  ,iy,iz+1)
            enddo
          enddo
        else
          do ix = 1, nx-esx, 2
            do iy = 0, ny-esy
              xywork(1,iy,ix-1) = +a(ix-1,iy,iz  )*norm
              xywork(2,iy,ix-1) = +a(ix-1,iy,iz+1)*norm
              xywork(1,iy,ix  ) = -a(ix,  iy,iz  )*norm
              xywork(2,iy,ix  ) = -a(ix,  iy,iz+1)*norm
            enddo
          enddo
        endif

c       --- do the transforms

        call vcpft (xywork(1,0,0), xywork(2,0,0), nx, 2*(ny+1),-1,ny+(1-esy),2)
        call vrpft2(xywork(1,0,0), xywork(2,0,0), nx, 2*(ny+1),   ny+(1-esy),2)

        do ix = 0, nx-1
          do iy = 0, ny-esy
            a(ix,iy  ,iz) = xywork(1,iy,ix)
            a(ix,iy,iz+1) = xywork(2,iy,ix)
          enddo
        enddo

      enddo
!$OMP END DO

      return
      end

c=============================================================================
      subroutine vpftxi(nx,ny,nz,a,norm,esx,esy,xywork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,xywork(2,0:nx,0:ny)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iz = 0, nz-2, 2

        do ix = 0, nx-esx
          do iy = 0, ny-esy
            xywork(1,iy,ix) = a(ix,iy,iz)
            xywork(2,iy,ix) = a(ix,iy,iz+1)
          enddo
        enddo

c       --- do the transforms

        call vrpfti2(xywork(1,0,0),xywork(2,0,0),nx, 2*(ny+1),   ny+(1-esy),2)
        call vcpft  (xywork(1,0,0),xywork(2,0,0),nx, 2*(ny+1),+1,ny+(1-esy),2)

c       --- shift k-space origin back and apply normalization factor
c       --- if non-zero.
        if (norm == 0.) then
          do ix = 1, nx-esx, 2
            do iy = 0, ny-esy
              a(ix-1,iy,iz  ) = +xywork(1,iy,ix-1)
              a(ix-1,iy,iz+1) = +xywork(2,iy,ix-1)
              a(ix,  iy,iz  ) = -xywork(1,iy,ix  )
              a(ix,  iy,iz+1) = -xywork(2,iy,ix  )
            enddo
          enddo
        else
          do ix = 1, nx-esx, 2
            do iy = 0, ny-esy
              a(ix-1,iy,iz  ) = +xywork(1,iy,ix-1)*norm
              a(ix-1,iy,iz+1) = +xywork(2,iy,ix-1)*norm
              a(ix,  iy,iz  ) = -xywork(1,iy,ix  )*norm
              a(ix,  iy,iz+1) = -xywork(2,iy,ix  )*norm
            enddo
          enddo
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
c Periodic forward and reverse transforms in y
c=============================================================================
      subroutine vpfty(nx,ny,nz,a,norm,esx,esy,xywork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,xywork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iz = 0, nz-2, 2
c       --- shift k-space origin to middle and apply normalization factor
c       --- if non-zero.
        if (norm == 0) then
          do iy = 1, ny-esy, 2
            do ix = 0, nx-esx
              xywork(1,ix,iy-1) = +a(ix,iy-1,iz  )
              xywork(2,ix,iy-1) = +a(ix,iy-1,iz+1)
              xywork(1,ix,iy  ) = -a(ix,iy  ,iz  )
              xywork(2,ix,iy  ) = -a(ix,iy  ,iz+1)
            enddo
          enddo
        else
          do iy = 1, ny-esy, 2
            do ix = 0, nx-esx
              xywork(1,ix,iy-1) = +a(ix,iy-1,iz  )*norm
              xywork(2,ix,iy-1) = +a(ix,iy-1,iz+1)*norm
              xywork(1,ix,iy  ) = -a(ix,iy  ,iz  )*norm
              xywork(2,ix,iy  ) = -a(ix,iy  ,iz+1)*norm
            enddo
          enddo
        endif

c       --- do the transforms

        call vcpft (xywork(1,0,0),xywork(2,0,0),ny, 2*(nx+1),-1,nx+(1-esx),2)
        call vrpft2(xywork(1,0,0),xywork(2,0,0),ny, 2*(nx+1),   nx+(1-esx),2)

        do iy = 0, ny-esy
          do ix = 0, nx-esx
            a(ix,iy,iz  ) = xywork(1,ix,iy)
            a(ix,iy,iz+1) = xywork(2,ix,iy)
          enddo
        enddo

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine vpftyi(nx,ny,nz,a,norm,esx,esy,xywork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,xywork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iz = 0, nz-2, 2

        do iy = 0, ny-esy
          do ix = 0, nx-esx
            xywork(1,ix,iy) = a(ix,iy,iz  )
            xywork(2,ix,iy) = a(ix,iy,iz+1)
          enddo
        enddo

c       --- do the transforms

        call vrpfti2(xywork(1,0,0),xywork(2,0,0), ny, 2*(nx+1),   nx+(1-esx),2)
        call vcpft  (xywork(1,0,0),xywork(2,0,0), ny, 2*(nx+1),+1,nx+(1-esx),2)

c       --- shift k-space origin back and apply normalization factor
c       --- if non-zero.
        if (norm == 0.) then
          do iy = 1, ny-esy, 2
            do ix = 0, nx-esx
              a(ix,iy-1,iz  ) = +xywork(1,ix,iy-1)
              a(ix,iy-1,iz+1) = +xywork(2,ix,iy-1)
              a(ix,iy  ,iz  ) = -xywork(1,ix,iy  )
              a(ix,iy  ,iz+1) = -xywork(2,ix,iy  )
            enddo
          enddo
        else
          do iy = 1, ny-esy, 2
            do ix = 0, nx-esx
              a(ix,iy-1,iz  ) = +xywork(1,ix,iy-1)*norm
              a(ix,iy-1,iz+1) = +xywork(2,ix,iy-1)*norm
              a(ix,iy  ,iz  ) = -xywork(1,ix,iy  )*norm
              a(ix,iy  ,iz+1) = -xywork(2,ix,iy  )*norm
            enddo
          enddo
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=============================================================================
c=============================================================================
c=============================================================================
#ifdef ESSL
c=======================================================================
c The following two subroutines use the dcft library routine to do the
c FFT's. That routine is faster when the data passed into it is in normal
c order, that is the data to be transformed is the first index, vectorized
c over the second index. However, when vcpft is used, it is faster when the
c data is in transposed order, that is the data to be transformed is the
c second index, vectorized over the first index. Both version were kept.
c The difference in speed is around 15% using the Intel compiler on Linux
c and about a factor of 2 on the IBM SP.
c Using the dcft routines, the code runs about 30-40% faster than with vcpft.
c=============================================================================
      subroutine vpftz(nx,ny,nz,a,norm,esx,esy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,zwork(2,0:nz,0:nx)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

      integer(ISZ):: naux1,naux2
      integer(ISZ):: maxnaux1,maxnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
      data maxnaux1/0/,maxnaux2/0/
      save maxnaux1,maxnaux2,aux1,aux2
 
c Setup and make call to the ESSL routine for sine-transforms

c     --- Calculate how much space is needed. Only allocate if it is
c     --- different from last time (deallocating if needed).
      if (nz <= 2048) then
        naux1 = 20000
      else
        naux1 = 20000+2.28*nz
      endif
      if (nz <= 2048) then
        naux2 = 20000
      else
        naux2 = 20000+2.28*nz
      endif
      if (nz >= 252) naux2 = naux2 + (2*nz+256)*(min(64, nx+(1-esx)))

      if (naux1 > maxnaux1) then
        maxnaux1 = naux1
        if (allocated(aux1)) deallocate(aux1)
        allocate(aux1(maxnaux1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vpftz: allocation error ",allocerror,
     &           ": could not allocate aux1 to shape ",maxnaux1
          call kaboom("vpftz: allocation error")
          return
        endif
      endif
      if (naux2 > maxnaux2) then
        maxnaux2 = naux2
        if (allocated(aux2)) deallocate(aux2)
        allocate(aux2(maxnaux2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vpftz: allocation error ",allocerror,
     &           ": could not allocate aux2 to shape ",maxnaux2
          call kaboom("vpftz: allocation error")
          return
        endif
      endif

      call dcft(1_4,zwork,1_4,int(nz+1,4),zwork,1_4,int(nz+1,4),
     &          int(nz,4),int(nx+(1-esx),4),int(+1,4),1.,
     &          aux1,int(naux1,4),aux2,int(naux2,4))
 
!$OMP DO FIRSTPRIVATE(norm)
      do iy = 0, ny-2+(1-esy), 2
c       --- shift k-space origin to middle and apply normalization factor
c       --- if non-zero.
        if (norm == 0) then
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              zwork(1,iz-1,ix) = +a(ix,iy  ,iz-1)
              zwork(2,iz-1,ix) = +a(ix,iy+1,iz-1)
              zwork(1,iz  ,ix) = -a(ix,iy  ,iz  )
              zwork(2,iz  ,ix) = -a(ix,iy+1,iz  )
            enddo
          enddo
        else
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              zwork(1,iz-1,ix) = +a(ix,iy  ,iz-1)*norm
              zwork(2,iz-1,ix) = +a(ix,iy+1,iz-1)*norm
              zwork(1,iz  ,ix) = -a(ix,iy  ,iz  )*norm
              zwork(2,iz  ,ix) = -a(ix,iy+1,iz  )*norm
            enddo
          enddo
        endif

c       --- do the transforms

        call dcft(0_4,zwork,1_4,int(nz+1,4),zwork,1_4,
     &            int(nz+1,4),int(nz,4),int(nx+(1-esx),4),int(+1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))
        call vrpft2(zwork(1,0,0), zwork(2,0,0), nz, 2,   nx+(1-esx),2*(nz+1))

        do iz = 0, nz-1
          do ix = 0, nx-esx
            a(ix,iy  ,iz) = zwork(1,iz,ix)
            a(ix,iy+1,iz) = zwork(2,iz,ix)
          enddo
        enddo

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine vpftzi(nx,ny,nz,a,norm,esx,esy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,zwork(2,0:nz,0:nx)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

      integer(ISZ):: naux1,naux2
      integer(ISZ):: maxnaux1,maxnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
      data maxnaux1/0/,maxnaux2/0/
      save maxnaux1,maxnaux2,aux1,aux2
 
c Setup and make call to the ESSL routine for sine-transforms

c     --- Calculate how much space is needed. Only allocate if it is
c     --- different from last time (deallocating if needed).
      if (nz <= 2048) then
        naux1 = 20000
      else
        naux1 = 20000+2.28*nz
      endif
      if (nz <= 2048) then
        naux2 = 20000
      else
        naux2 = 20000+2.28*nz
      endif
      if (nz >= 252) naux2 = naux2 + (2*nz+256)*(min(64, nx+(1-esx)))

      if (naux1 > maxnaux1) then
        maxnaux1 = naux1
        if (allocated(aux1)) deallocate(aux1)
        allocate(aux1(maxnaux1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vpftzi: allocation error ",allocerror,
     &           ": could not allocate aux1 to shape ",maxnaux1
          call kaboom("vpftzi: allocation error")
          return
        endif
      endif
      if (naux2 > maxnaux2) then
        maxnaux2 = naux2
        if (allocated(aux2)) deallocate(aux2)
        allocate(aux2(maxnaux2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vpftzi: allocation error ",allocerror,
     &           ": could not allocate aux2 to shape ",maxnaux2
          call kaboom("vpftzi: allocation error")
          return
        endif
      endif

      call dcft(1_4,zwork,1_4,int(nz+1,4),zwork,1_4,int(nz+1,4),
     &          int(nz,4),int(nx+(1-esx),4),int(-1,4),1.,
     &          aux1,int(naux1,4),aux2,int(naux2,4))

!$OMP DO FIRSTPRIVATE(norm)
      do iy = 0, ny-2+(1-esy), 2

        do iz = 0, nz-1
          do ix = 0, nx-esx
            zwork(1,iz,ix) = a(ix,iy  ,iz)
            zwork(2,iz,ix) = a(ix,iy+1,iz)
          enddo
        enddo

c       --- do the transforms

        call vrpfti2(zwork(1,0,0), zwork(2,0,0), nz, 2,   nx+(1-esx),2*(nz+1))
        call dcft(0_4,zwork,1_4,int(nz+1,4),zwork,1_4,
     &            int(nz+1,4),int(nz,4),int(nx+(1-esx),4),int(-1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))

c       --- shift k-space origin back and apply normalization factor
c       --- if non-zero.
        if (norm == 0.) then
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              a(ix,iy  ,iz-1) = +zwork(1,iz-1,ix)
              a(ix,iy+1,iz-1) = +zwork(2,iz-1,ix)
              a(ix,iy  ,iz  ) = -zwork(1,iz  ,ix)
              a(ix,iy+1,iz  ) = -zwork(2,iz  ,ix)
            enddo
          enddo
        else
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              a(ix,iy  ,iz-1) = +zwork(1,iz-1,ix)*norm
              a(ix,iy+1,iz-1) = +zwork(2,iz-1,ix)*norm
              a(ix,iy  ,iz  ) = -zwork(1,iz  ,ix)*norm
              a(ix,iy+1,iz  ) = -zwork(2,iz  ,ix)*norm
            enddo
          enddo
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
#else
c=======================================================================
c The following two routines use the vcpft routine for the FFT's. That
c routine runs faster when the data is in transposed order and the following
c versions were kept.
c The difference in speed is around 15% using the Intel compiler on Linux
c and about a factor of 2 on the IBM SP.
c=============================================================================
      subroutine vpftz(nx,ny,nz,a,norm,esx,esy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,zwork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iy = 0, ny-2+(1-esy), 2
c       --- shift k-space origin to middle and apply normalization factor
c       --- if non-zero.
        if (norm == 0) then
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              zwork(1,ix,iz-1) = +a(ix,iy  ,iz-1)
              zwork(2,ix,iz-1) = +a(ix,iy+1,iz-1)
              zwork(1,ix,iz  ) = -a(ix,iy  ,iz  )
              zwork(2,ix,iz  ) = -a(ix,iy+1,iz  )
            enddo
          enddo
        else
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              zwork(1,ix,iz-1) = +a(ix,iy  ,iz-1)*norm
              zwork(2,ix,iz-1) = +a(ix,iy+1,iz-1)*norm
              zwork(1,ix,iz  ) = -a(ix,iy  ,iz  )*norm
              zwork(2,ix,iz  ) = -a(ix,iy+1,iz  )*norm
            enddo
          enddo
        endif

c       --- do the transforms

        call vcpft (zwork(1,0,0), zwork(2,0,0), nz, 2*(nx+1),-1,nx+(1-esx),2)
        call vrpft2(zwork(1,0,0), zwork(2,0,0), nz, 2*(nx+1),   nx+(1-esx),2)

        do iz = 0, nz-1
          do ix = 0, nx-esx
            a(ix,iy  ,iz) = zwork(1,ix,iz)
            a(ix,iy+1,iz) = zwork(2,ix,iz)
          enddo
        enddo

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine vpftzi(nx,ny,nz,a,norm,esx,esy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm,zwork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.

      integer(ISZ):: ix,iy,iz

!$OMP DO FIRSTPRIVATE(norm)
      do iy = 0, ny-2+(1-esy), 2

        do iz = 0, nz-1
          do ix = 0, nx-esx
            zwork(1,ix,iz) = a(ix,iy  ,iz)
            zwork(2,ix,iz) = a(ix,iy+1,iz)
          enddo
        enddo

c       --- do the transforms

        call vrpfti2(zwork(1,0,0), zwork(2,0,0), nz, 2*(nx+1),   nx+(1-esx),2)
        call vcpft  (zwork(1,0,0), zwork(2,0,0), nz, 2*(nx+1),+1,nx+(1-esx),2)

c       --- shift k-space origin back and apply normalization factor
c       --- if non-zero.
        if (norm == 0.) then
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              a(ix,iy  ,iz-1) = +zwork(1,ix,iz-1)
              a(ix,iy+1,iz-1) = +zwork(2,ix,iz-1)
              a(ix,iy  ,iz  ) = -zwork(1,ix,iz  )
              a(ix,iy+1,iz  ) = -zwork(2,ix,iz  )
            enddo
          enddo
        else
          do iz = 1, nz-1, 2
            do ix = 0, nx-esx
              a(ix,iy  ,iz-1) = +zwork(1,ix,iz-1)*norm
              a(ix,iy+1,iz-1) = +zwork(2,ix,iz-1)*norm
              a(ix,iy  ,iz  ) = -zwork(1,ix,iz  )*norm
              a(ix,iy+1,iz  ) = -zwork(2,ix,iz  )*norm
            enddo
          enddo
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
#endif
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine attenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,esx,esy,esz)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy,esz
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

!$OMP DO
      do ikz = 0, nz-esz
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*attx(ikx)*atty(iky)*attz(ikz)
          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine unattenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,
     &                       esx,esy,esz)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy,esz
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz-1)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

!$OMP DO
      do ikz = 0, nz-esz
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)/(attx(ikx)*atty(iky)*attz(ikz))
          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine rhotophi(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,esx,esy,esz)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy,esz
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz)
      real(kind=8):: kxsq(0:nx-esx), kysq(0:ny-esy), kzsq(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz
      real(kind=8):: oneok

!$OMP DO
      do ikz = 0, nz-esz
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            oneok = 1./(kxsq(ikx) + kysq(iky) + kzsq(ikz))
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*oneok
          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine phitorho(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,esx,esy,esz)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy,esz
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz)
      real(kind=8):: kxsq(0:nx-esx), kysq(0:ny-esy), kzsq(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

!$OMP DO
      do ikz = 0, nz-esz
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*(kxsq(ikx) + kysq(iky) + kzsq(ikz))
          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine tridiag(nx,ny,nz,a,norm,kxsq,kysq,kzsq,ikxmin,ikymin,esy,
     &                   zwork)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esy
      real(kind=8):: a(-1:nx+1,-1:ny+1,0:nz),norm
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: zwork(0:nx,0:nz)

c     --- Use a(nx,,0) to hold diagonal matrix elements.
c     --- Use a(nx,,1:nz-1) and a(nx,,nz) as temp arrays in tridiag.
c     --- The lines iz=-1 and nz+1 should be used as scratch arrays to allow
c     --- for Neumann or Periodic boundary conditions in z, where the
c     --- 2-d scratch array would be phi(nx,,0:nz), but phi is declared
c     --- as phi(,,0:nz) so the iz=-1 and nz+1 elements are not accessible.
c     --- Off-diagonal terms are all the constant (-1.).
c     --- I know that this looks messy and nearly undecipherable, but I
c     --- wanted to do this without creating extra temporary arrays,
c     --- without using the scrtch array which is holding data needed
c     --- for the transverse FFTs, and without using the planes at ix=iy=0
c     --- since those planes are used when there are transverse symmetries.
c     --- (signed DPG)

      integer(ISZ):: ix,iy,iz

!$OMP DO
      do iy=ikymin,ny-esy
        do iz=1,nz-1
          do ix=ikxmin,nx-1
c           --- the RHS is the transverse FFTs of rho times dz**2/eps0
c           --- This multiply could actually be done within the tridiag
c           --- solving loops below, making them somewhat more complicated.
c           --- That saves a few percent in the run time.
            a(ix,iy,iz) = a(ix,iy,iz)*norm
          enddo
        enddo
        do ix=ikxmin,nx-1
c         --- diagonal matrix elements
          zwork(ix,0) = 2. + (kxsq(ix)+kysq(iy))*norm
c         --- set the end points using Dirichlet boundary conditions.
          a(ix,iy,1) = a(ix,iy,1) + a(ix,iy,0)
          a(ix,iy,nz-1) = a(ix,iy,nz-1) + a(ix,iy,nz)
        enddo

c       --- do tridiag matrix solve
        do ix = ikxmin,nx-1
          zwork(ix,nz) = zwork(ix,0)
          a(ix,iy,1) = a(ix,iy,1)/zwork(ix,nz)
        enddo
        do iz = 2,nz-1
          do ix=ikxmin,nx-1
            zwork(ix,iz) = -1./zwork(ix,nz)
            zwork(ix,nz) = zwork(ix,0) - (-1.)*zwork(ix,iz)
            a(ix,iy,iz) = (a(ix,iy,iz) - (-1.)*a(ix,iy,iz-1))/zwork(ix,nz)
          enddo
        enddo
        do iz = nz-2,1,-1
          do ix=ikxmin,nx-1
            a(ix,iy,iz) = a(ix,iy,iz) - zwork(ix,iz+1)*a(ix,iy,iz+1)
          enddo
        enddo

      enddo
!$OMP END DO

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine vpois2d(iwhich,a,ak,kxsq,kysq,attx,atty,filt,lx,ly,nx,ny,
     &                   w,xywork,ibc,l2symtry,l4symtry)
      use Constant
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: iwhich,nx,ny,ibc
      real(kind=8)::  a(-1:nx+1, -1:ny+1)
      real(kind=8):: ak(-1:nx+1, -1:ny+1)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), filt(5,2)
      real(kind=8):: lx, ly, w(*), xywork(0:nx,0:ny)
      logical(ISZ):: l2symtry,l4symtry
      
c-----------------------------------------------------------------------------
c               Vectorized 2d Poisson solver:
c               Sine transform in x and y
c               Alex Friedman, LLNL, June 1989.
c
c          Returned value: 0 if ok, 1 if error
c
c          IWHICH = ... -1: full fieldsolve; assumes kxsq etc. and attx etc.
c                           have been set already.  This is equivalent to a
c                           sequence of calls with iwhich = 2,3,4,5.
c                        0: full fieldsolve; kxsq etc. will be set herein.
c                           This is equivalent to a sequence of calls with
c                           iwhich = 1,2,3,4,5. (THIS IS THE SIMPLEST USAGE.)
c                        1: set kxsq etc. and attx etc., return.
c                        2: forward transform (no k^2 factor), return.
c                        3: apply k-space filter using attx, etc., return.
c                        4: apply inverse-k^2 factor (rho to phi), return.
c                        5: inverse transform (no k^2 factor), return.
c                        8: apply reciprocal of filter in k-space, return.
c                        9: apply k^2 factor (turns phi into rho), return.
c              A, AK ... Array to transform (pass it to this routine twice)
c                        Use of both a and ak for consistency with the
c                        3d package.
c                        See array declarations above these comments.
c         KXSQ, KYSQ ... Discrete analogs to eps0*kx^2 etc.
c                        See array declarations below these comments.
c         ATTX, ATTY ... Attenuation factors for each mode number in x, y.
c                        See array declarations below these comments.
c          FILT(i,j) ... Spatial filtering coefficients; array should be
c                        dimensioned (5,2) in calling routine.
c                        Here j is the axis (1=x, 2=y),
c                        and there are four coefficients for each axis:
c                        filt(1,j): coefficient a1, of sin^2 term (boosting)
c                        filt(2,j): coefficient a2, of tan^4 term (smoothing)
c                        filt(3,j): cutoff k dx or k dy, e.g. pi/2
c                        filt(4,j): exponent N (sharpness of cutoff) e.g. 8.
c                        filt(5,j): W factor in binomial filter, e.g. 0.5
c                        (A value of .5 gives the familiar 1/4:1/2:1/4 filter)
c                        See Birdsall and Langdon, "Plasma Physics via
c                        Computer Simulation," Appendices B and C.
c               EPS0 ... Epsilon-0; use 1/4Pi for CGS, 1 for rationalized.
c             LX, LY ... System lengths
c             NX, NY ... Number of mesh points (mesh runs from 0 to NX, etc.)
c                        At present these must be powers of two.
c                  W ... Workspace, length at least  (NX+2)(NY+2)-7
c                IBC ... Boundary condition switch, char*8 (for future use)
c
c  Note that for x (y is similar), k dx = pi j/nx, j = 1, ..., nx-1.
c        The longest wavelength is 2 nx dx.
c        The j=nx (odd-even) term is absent because the boundary pts are zero.
c
c  This routine was inspired by Langdon's scalar sine-periodic solver "sppois"
c-----------------------------------------------------------------------------
c The parallel version works by having subgroups of processors that
c cooperate on the calculation. Each processor on the group works on one
c chunk of the array, and then a gatherall is done to redistribute the
c result among all of the processors in the group. The coding is organized
c so that during a complete field solve, the number of gathers is
c minimized.  Note - if only one operation is done, i.e.  iwhich = 2, 3,
c or 4, then the data returned has not been gathered, so none of the
c processors has the complete result. The one exception is the inverse
c FFT, iwhich = 5, which always calls the gather since it is the last
c operation anyway.
c The number of processors in each group is set by the variable
c xynppgroup (it defaults to 1 which means the no parallelization is done).
c The variable needs to be optimized through experimentation.
c-----------------------------------------------------------------------------

      integer(ISZ):: iwxy,icmx,icpx,icmy,icpy,ikx,iky,ix,iy
      integer(ISZ):: icox,icoy,ikxmin,ikymin
      integer(ISZ):: ixstart,ixend,iystart,iyend
      real(kind=8):: klast, kdxb2, kdyb2, kdx, kdy
      real(kind=8):: dx,dy,pib2,norm
      real(kind=8):: atemp(0:ny,0:nx)
      real(kind=8):: w1(-1:3,-1:ny+1),w2(-1:3,-1:ny+1)
      w1 = 0.
      w2 = 0.

c  Error checking; so far, only possible error is non-power-of-2 nx,ny

      if (  ( nx .ne. 2**int( .5 + log(1.*nx) / log(2.) ) )
     & .or. ( ny .ne. 2**int( .5 + log(1.*ny) / log(2.) ) ) ) then
         call kaboom("NX and NY must be powers of 2")
         return
      endif

c  Set pointers into workspace needed by x and y sine transforms

      iwxy  = 1
c     icmx  = iwxy  + (nx+1)*(ny+1)
      icmx  = iwxy
      icpx  = icmx  + nx/2-1
      icmy  = icpx  + nx/2-1
      icpy  = icmy  + ny/2-1
      icox  = icpy  + ny/2-1
      icoy  = icox  + nx
c     inext = icoy  + ny is first word beyond workspace

c  Set useful constants

      dx = lx / nx
      dy = ly / ny
      pib2 = .5 * pi

c  Minimum index of kxsq and kysq that is used
      ikxmin = 1
      if (l4symtry) ikxmin = 0
      ikymin = 1
      if (l2symtry .or. l4symtry) ikymin = 0

c  Start and end of loops
      ixstart = 0
      ixend = nx-1
      iystart = 0
      iyend = ny-1
#ifdef MPIPARALLEL
      if (xynppgroup > nzprocs) xynppgroup = nzprocs
      if (nx > xynppgroup .and. xynppgroup > 1) then
        ixstart = mod(my_index,xynppgroup)*nx/xynppgroup
        ixend = ixstart + nx/xynppgroup - 1
      endif
      if (ny > xynppgroup .and. xynppgroup > 1) then
        iystart = mod(my_index,xynppgroup)*ny/xynppgroup
        iyend = iystart + ny/xynppgroup - 1
      endif
#endif

c ----------------------------------------------------------------------------
c  For poisson equation, kxsq(ikx) is a discrete
c  analog to eps0 kx^2, etc.  If the user has requested it,
c  we set up ksq arrays for a seven point scheme now; the coding for kxsq
c  etc. is arranged so that it will vectorize on "most" compilers.
c  Also, compute attenuation factors as functions of mode numbers in x, y.
c ----------------------------------------------------------------------------

      if (iwhich == 0 .or. iwhich == 1) then

c        --- compute x direction rho-to-phi coefficients ---
         do ikx = ikxmin, nx-1
           kxsq(ikx) = ikx
         enddo
         if (l4symtry) then
           do ikx = 0, nx-1
             kxsq(ikx) = 2./dx**2*(1.-cos(pib2*(2.*kxsq(ikx)+1.)/nx)) * eps0
           enddo
         else
           do ikx = 1, nx-1
             kxsq(ikx) = (2./dx*sin(pib2*kxsq(ikx)/nx))**2 * eps0
           enddo
         endif
c        --- compute y direction coefficients ---
         do iky = ikymin, ny-1
           kysq(iky) = iky
         enddo
         if (l2symtry .or. l4symtry) then
           do iky = 0, ny-1
             kysq(iky) = 2./dy**2*(1.-cos(pib2*(2.*kysq(iky)+1.)/ny)) * eps0
           enddo
         else
           do iky = 1, ny-1
              kysq(iky) = (2./dy*sin(pib2*kysq(iky)/ny))**2 * eps0
           enddo
         endif

c        --- spatial filtering in x; binomial filter, or unity
         do ikx = ikxmin, nx-1
            kdx = pi * ikx / nx
            attx(ikx) = (1.+2.*filt(5,1)*cos(kdx)) / (1.+2.*filt(5,1))
         enddo
         if (filt(1,1) .ne. 0. .or. filt(2,1) .ne. 0.) then
c           --- compute first form of spatial filtering
            do ikx = ikxmin, nx-1
               kdxb2 = pi * ikx / (2.*nx)
               attx(ikx) = attx(ikx) * (exp( filt(1,1)*(sin(kdxb2))**2
     &          - filt(2,1)*(tan(kdxb2))**4 ))**2
            enddo
         endif
         if (filt(3,1) .ne. 0.) then
c           --- compute second form of spatial filtering
            klast = filt(3,1) * nx / pi
            do ikx = ikxmin, nx-1
               attx(ikx) = attx(ikx) * exp(-(ikx/klast)**filt(4,1))
            enddo
         endif

c        --- spatial filtering in y; binomial filter, or unity
         do iky = ikymin, ny-1
            kdy = pi * iky / ny
            atty(iky) = (1.+2.*filt(5,2)*cos(kdy)) / (1.+2.*filt(5,2))
         enddo
         if (filt(1,2) .ne. 0. .or. filt(2,2) .ne. 0.) then
c           --- compute first form of spatial filtering
            do iky = ikymin, ny-1
               kdyb2 = pi * iky / (2.*ny)
               atty(iky) = atty(iky) * (exp( filt(1,2)*(sin(kdyb2))**2
     &          - filt(2,2)*(tan(kdyb2))**4 ))**2
            enddo
         endif
         if (filt(3,2) .ne. 0.) then
c           --- compute second form of spatial filtering
            klast = filt(3,2) * ny / pi
            do iky = ikymin, ny-1
               atty(iky) = atty(iky) * exp(-(iky/klast)**filt(4,2))
            enddo
         endif

      endif

c ----------------------------------------------------------------------------
c  Set up the tables needed by the sine transforms (call them w/ ISETUP=1)
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 2 .or. iwhich == 5) then

        if (l4symtry) then
          call cosqx (a, xywork, w(icox), nx, ny, 1,-1)
          call cosqy (a, xywork, w(icoy), nx, ny, 1,-1)
        elseif (l2symtry) then
          call vsftx (a, xywork, w(icpx), w(icmx), nx, ny, 1)
          call cosqy (a, xywork, w(icoy), nx, ny, 1,-1)
        else
          call vsftx (a, xywork, w(icpx), w(icmx), nx, ny, 1)
          call vsfty (a, xywork, w(icpy), w(icmy), nx, ny, 1)
        endif

      endif

!$OMP PARALLEL PRIVATE(xywork)

c ----------------------------------------------------------------------------
c  Do the forward transform
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 2) then

        if (l4symtry) then
c         --- factor of .5dx for cosqx and .5*dy for cosqy
          norm = 0.25*dx*dy
        else if (l2symtry) then
c         --- factor of .5*dy for cosqy
          norm = 0.5*dx*dy
        else
          norm = dx*dy
        endif

c       --- Invoke the transform routines for the current z-slice
        if (l4symtry) then
!$OMP DO
          do iy = iystart,iyend,2
            a(:,iy:iy+1) = a(:,iy:iy+1)*norm
            call cosqx(a(:,iy-1:iy+2),xywork(0,iy),w(icox),nx,2,0,-1)
          enddo
!$OMP END DO NOWAIT
        else
!$OMP DO
          do iy = iystart,iyend,2
            a(:,iy:iy+1) = a(:,iy:iy+1)*norm
            call vsftx(a(:,iy-1:iy+2),xywork(0,iy),w(icpx),w(icmx),nx,2,0)
          enddo
!$OMP END DO NOWAIT
        endif

#ifdef MPIPARALLEL
        if (ny > xynppgroup .and. xynppgroup > 1)
     &    call parallelgatherall(a,(1+nx)*ny/xynppgroup,xynppgroup,1)
#endif
!$OMP FLUSH (a)

        if (l4symtry .or. l2symtry) then
!$OMP DO
          do ix = ixstart,ixend,2
c           --- The arrays are copied into temporary arrays to get around
c           --- a bug in the intel compiler, version 8.
            w1(0:2,0:ny) = a(ix:ix+2,0:ny)
            w2(0:2,0:ny) = xywork(ix:ix+2,0:ny)
            call cosqy(w1,w2,w(icoy),2,ny,0,-1)
            a(ix:ix+2,0:ny) = w1(0:2,0:ny)
            xywork(ix:ix+2,0:ny) = w2(0:2,0:ny)
          enddo
!$OMP END DO NOWAIT
        else
!$OMP DO
          do ix = ixstart,ixend,2
            call vsfty(a(ix-1:ix+3,:),xywork(ix:ix+2,:),w(icpy),w(icmy),2,ny,0)
          enddo
!$OMP END DO NOWAIT
        endif
!$OMP FLUSH (a)

      endif

c ----------------------------------------------------------------------------
c  Apply k-space filter; attenuate mode-by-mode
c ----------------------------------------------------------------------------

      if ((iwhich <= 0 .or. iwhich == 3) .and. maxval(abs(filt)) > 0.) then

!$OMP DO
        do iky = ikymin, ny-1
          do ikx = ixstart,ixend
           ak(ikx,iky) = ak(ikx,iky) * attx(ikx) * atty(iky)
          enddo
        enddo
!$OMP END DO NOWAIT

      endif

c ----------------------------------------------------------------------------
c  Convert rhok to phik.  We should multitask this loop.
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 4) then

!$OMP DO
        do iky = ikymin, ny-1
          do ikx = ixstart,ixend
           ak(ikx,iky) = ak(ikx,iky)/(kxsq(ikx) + kysq(iky))
          enddo
        enddo
!$OMP END DO NOWAIT
!$OMP FLUSH (a)

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 5) then

        if (l4symtry) then
c         --- factor of .25/lx for cosqx and .25/ly for cosqy
          norm = .0625/(lx*ly)
        else if (l2symtry) then
c         --- factor of 2./lx for vsftx and .25/ly for cosqy
          norm = .5/(lx*ly)
        else
          norm = 4./(lx*ly)
        endif

        if (l4symtry .or. l2symtry) then
c         --- Invoke the inverse transform routines for the current z-slice
!$OMP DO
          do ix = ixstart,ixend,2
c           --- The arrays are copied into temporary arrays to get around
c           --- a bug in the intel compiler, version 8.
            w1(0:2,0:ny) = a(ix:ix+2,0:ny)
            w2(0:2,0:ny) = xywork(ix:ix+2,0:ny)
            call cosqy(w1,w2,w(icoy),2,ny,0,1)
            a(ix:ix+2,0:ny) = w1(0:2,0:ny)
            xywork(ix:ix+2,0:ny) = w2(0:2,0:ny)
          enddo
!$OMP END DO NOWAIT
        else
!$OMP DO
          do ix = ixstart,ixend,2
            call vsfty(a(ix-1:ix+3,:),xywork(ix:ix+2,:),w(icpy),w(icmy),2,ny,0)
          enddo
!$OMP END DO NOWAIT
        endif

#ifdef MPIPARALLEL
        if (nx > xynppgroup .and. xynppgroup > 1) then
          atemp = transpose(a)
          call parallelgatherall(atemp,(1+ny)*nx/xynppgroup,xynppgroup,1)
          a = transpose(atemp)
        endif
#endif

!$OMP FLUSH (a)

        if (l4symtry) then
!$OMP DO
          do iy = iystart,iyend,2
            call cosqx(a(:,iy-1:iy+2),xywork(0,iy),w(icox),nx,2,0,1)
            a(0:nx,iy:iy+1) = a(0:nx,iy:iy+1)*norm
          enddo
!$OMP END DO NOWAIT
        else
!$OMP DO
          do iy = iystart,iyend,2
            call vsftx(a(:,iy-1:iy+2),xywork(0,iy),w(icpx),w(icmx),nx,2,0)
            a(0:nx,iy:iy+1) = a(0:nx,iy:iy+1)*norm
          enddo
!$OMP END DO NOWAIT
        endif

!$OMP FLUSH (a)
!$OMP SINGLE
        a(nx,:) = 0.
        a(:,ny) = 0.
        if (.not. l4symtry) then
          a(0,:) = 0.
        endif
        if (.not. l2symtry .and. .not. l4symtry) then
          a(:,0) = 0.
        endif
!$OMP END SINGLE

#ifdef MPIPARALLEL
        if (ny > xynppgroup .and. xynppgroup > 1)
     &    call parallelgatherall(a,(1+nx)*ny/xynppgroup,xynppgroup,1)
#endif

      endif

c ----------------------------------------------------------------------------
c  The usual fieldsolve has been completed.  Utility operations follow.
c ----------------------------------------------------------------------------

c ----------------------------------------------------------------------------
c  Un-filter in k space
c ----------------------------------------------------------------------------

      if (iwhich == 8) then

!$OMP DO
      do iky = ikymin, ny-1
        do ikx = ikxmin, nx-1
          ak(ikx,iky) = ak(ikx,iky) / (attx(ikx) * atty(iky))
        enddo
      enddo
!$OMP END DO

      endif

c ----------------------------------------------------------------------------
c  Multiply by k^2 factor, in contrast to the usual division.
c  This is useful if one wants to recreate rho from phi.
c ----------------------------------------------------------------------------

      if (iwhich == 9) then

!$OMP DO
      do iky = ikymin, ny-1
        do ikx = ikxmin, nx-1
          ak(ikx,iky) = ak(ikx,iky) * ( kxsq(ikx) + kysq(iky) )
        enddo
      enddo
!$OMP END DO

      endif

c ----------------------------------------------------------------------------
c  End of VPOIS2D
c ----------------------------------------------------------------------------

!$OMP END PARALLEL

      return
      end
c======================================================================
      subroutine vsftx (a, w, cp, cm, nx, ny, isetup)
      use Constant
      integer(ISZ):: nx,ny,isetup
      real(kind=8):: a(-1:nx+1,-1:ny+1), w(0:nx,0:ny)
      real(kind=8):: cp(nx/2-1), cm(nx/2-1)

c Many sine transforms in x; nx must be a power of two
c Uses vectorized routine VCPFT
c If ISETUP = 1, does only a setup and returns
c This routine is its own inverse, except for a factor of 2/nx
c In comparison with Langdon's routine MSTX,
c    r is stored in w( ,iy+1);    i is stored in w( ,iy).

c -----------------------------------------------------------------
#ifndef ESSL
c Written by Alex Friedman, LLNL, April 1989.


      real(kind=8):: aa,bb
      integer(ISZ):: nx2,ix,iy,ix2
      save nx2

c Initialize
      if (isetup == 1) then
         nx2 = nx / 2
         do ix = 1, nx2-1
            cp(ix) = 1. / ( 8.*sin(pi/nx * ix) ) + .25
            cm(ix) = 1. / ( 8.*sin(pi/nx * ix) ) - .25
         enddo
         return
      endif

c Set up for transform
      do iy = 0, ny-1
         w(0,iy) = 2. * a(1,iy)
      enddo
      do iy = 0, ny-2, 2
         do ix = 1, nx2-1
            ix2 = 2 * ix
            aa = a(ix2,iy)
            bb = a(ix2+1,iy+1) - a(ix2-1,iy+1)
            w(ix   ,iy+1) = bb + aa
            w(nx-ix,iy+1) = bb - aa
            aa = a(ix2+1,iy) - a(ix2-1,iy)
            bb = a(ix2,iy+1)
            w(ix   ,iy) = aa - bb
            w(nx-ix,iy) = aa + bb
         enddo
      enddo
      do iy = 0, ny-1
         w(nx2,iy) = -2. * a(nx-1,iy)
      enddo

c Do vector complex periodic transform
      call vcpft ( w(0,1), w(0,0), nx, 1, 1, ny/2, 2*(nx+1) )

c Separate sine coefficients
      do iy = 0, ny-1
         do ix = 1, nx2-1
            aa = w(ix   ,iy)
            bb = w(nx-ix,iy)
            a(ix   ,iy) = cp(ix)*aa + cm(ix)*bb
            a(nx-ix,iy) = cm(ix)*aa + cp(ix)*bb
         enddo
      enddo
      do iy = 0, ny-1
         a(nx2,iy) = .25 * w(nx2,iy)
      enddo

#else
c -----------------------------------------------------------------

      integer(ISZ):: naux1,naux2
      integer(ISZ):: maxnaux1,maxnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
      data naux1/0/,naux2/0/
      data maxnaux1/0/,maxnaux2/0/
      save aux1,aux2,maxnaux1,maxnaux2

c Setup and make call to the ESSL routine for sine-transforms

c     --- Calculate how much space is needed. Only allocate if it is
c     --- different from last time (deallocating if needed).
      if (2*nx <= 16384) then
        naux1 = 50000
      else
        naux1 = 20000+.60*nx*2
      endif
      if (2*nx <= 16384) then
        naux2 = 20000
      else
        naux2 = 20000+.64*nx*2
      endif

      if (naux1 > maxnaux1) then
        maxnaux1 = naux1
        if (allocated(aux1)) deallocate(aux1)
        allocate(aux1(maxnaux1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vsftx: allocation error ",allocerror,
     &           ": could not allocate aux1 to shape ",maxnaux1
          call kaboom("vsftx: allocation error")
          return
        endif
      endif
      if (naux2 > maxnaux2) then
        maxnaux2 = naux2
        if (allocated(aux2)) deallocate(aux2)
        allocate(aux2(maxnaux2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vsftx: allocation error ",allocerror,
     &           ": could not allocate aux2 to shape ",maxnaux2
          call kaboom("vsftx: allocation error")
          return
        endif
      endif

c     --- XXX Note that this will be inefficient since a copy of 'a' will
c     --- be made.
      call dsinf(1_4,a(0:nx,0:ny),1_4,int(1+nx,4),a(0:nx,0:ny),
     &           1_4,int(1+nx,4),int(2*nx,4),int(ny,4),
     &           1.,aux1,int(naux1,4),aux2,int(naux2,4))

      if (isetup == 1) return

      call dsinf(0_4,a(0:nx,0:ny),1_4,int(1+nx,4),a(0:nx,0:ny),
     &           1_4,int(1+nx,4),int(2*nx,4),int(ny,4),
     &           1.,aux1,int(naux1,4),aux2,int(naux2,4))

#endif

      return
      end
c======================================================================
      subroutine vsfty (a, w, cp, cm, nx, ny, isetup)
      use Constant
      integer(ISZ):: nx,ny,isetup
      real(kind=8):: a(-1:nx+1,-1:ny+1), w(0:nx,0:ny)
      real(kind=8):: cp(ny/2-1), cm(ny/2-1)

c Many sine transforms in y; ny must be a power of two
c Uses vectorized routine VCPFT
c If ISETUP = 1, does only a setup and returns
c This routine is its own inverse, except for a factor of 2/nx
c In comparison with Langdon's routine MSTX,
c    r is stored in w(ix+1, );    i is stored in w(ix, ).

c -----------------------------------------------------------------
#ifndef ESSL

c Written by Alex Friedman, LLNL, April 1989.


      real(kind=8):: aa,bb
      integer(ISZ):: ny2,ix,iy,iy2
      save ny2

c Initialize
      if (isetup == 1) then
         ny2 = ny / 2
         do iy = 1, ny2-1
            cp(iy) = 1. / ( 8.*sin(pi/ny * iy) ) + .25
            cm(iy) = 1. / ( 8.*sin(pi/ny * iy) ) - .25
         enddo
         return
      endif

c Set up for transform
      do ix = 0, nx-1
         w(ix,0) = 2. * a(ix,1)
      enddo
      do iy = 1, ny2-1
         iy2 = 2 * iy
         do ix = 0, nx-2, 2
            aa = a(ix,iy2)
            bb = a(ix+1,iy2+1) - a(ix+1,iy2-1)
            w(ix+1,iy   ) = bb + aa
            w(ix+1,ny-iy) = bb - aa
            aa = a(ix,iy2+1) - a(ix,iy2-1)
            bb = a(ix+1,iy2)
            w(ix,iy   ) = aa - bb
            w(ix,ny-iy) = aa + bb
         enddo
      enddo
      do ix = 0, nx-1
         w(ix,ny2) = -2. * a(ix,ny-1)
      enddo

c Do vector complex periodic transform
      call vcpft ( w(1,0), w(0,0), ny, nx+1, 1, nx/2, 2 )

c Separate sine coefficients
      do iy = 1, ny2-1
         do ix = 0, nx-1
            aa = w(ix,iy   )
            bb = w(ix,ny-iy)
            a(ix,iy   ) = cp(iy)*aa + cm(iy)*bb
            a(ix,ny-iy) = cm(iy)*aa + cp(iy)*bb
         enddo
      enddo
      do ix = 0, nx-1
         a(ix,ny2) = .25 * w(ix,ny2)
      enddo

#else
c -----------------------------------------------------------------

      integer(ISZ):: naux1,naux2
      integer(ISZ):: maxnaux1,maxnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
      data naux1/0/,naux2/0/
      data maxnaux1/0/,maxnaux2/0/
      save aux1,aux2,maxnaux1,maxnaux2

c Setup and make call to the ESSL routine for sine-transforms

c     --- Calculate how much space is needed. Only allocate if it is
c     --- different from last time (deallocating if needed).
      if (2*ny <= 16384) then
        naux1 = 50000
      else
        naux1 = 20000+.60*nx*2
      endif
      if (2*ny <= 16384) then
        naux2 = 20000
      else
        naux2 = 20000+.64*nx*2
      endif
      if (2*ny >= 252) naux2 = naux2 + (ny + 257)*min(128,1+nx)

      if (naux1 > maxnaux1) then
        maxnaux1 = naux1
        if (allocated(aux1)) deallocate(aux1)
        allocate(aux1(maxnaux1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vsfty: allocation error ",allocerror,
     &           ": could not allocate aux1 to shape ",maxnaux1
          call kaboom("vsfty: allocation error")
          return
        endif
      endif
      if (naux2 > maxnaux2) then
        maxnaux2 = naux2
        if (allocated(aux2)) deallocate(aux2)
        allocate(aux2(maxnaux2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vsfty: allocation error ",allocerror,
     &           ": could not allocate aux2 to shape ",maxnaux2
          call kaboom("vsfty: allocation error")
          return
        endif
      endif

c     --- XXX Note that this will be inefficient since a copy of 'a' will
c     --- be made.
      call dsinf(1_4,a(0:nx,0:ny),int(1+nx,4),1_4,a(0:nx,0:ny),
     &           int(1+nx,4),1_4,int(2*ny,4),int(nx,4),
     &           1.,aux1,int(naux1,4),aux2,int(naux2,4))

      if (isetup == 1) return
 
      call dsinf(0_4,a(0:nx,0:ny),int(1+nx,4),1_4,a(0:nx,0:ny),
     &           int(1+nx,4),1_4,int(2*ny,4),int(nx,4),
     &           1.,aux1,int(naux1,4),aux2,int(naux2,4))

#endif

      return
      end
c=======================================================================
#ifdef ESSL
c=======================================================================
c The following two subroutines use the dcft library routine to do the
c FFT's. That routine is faster when the data passed into it is in normal
c order, that is the data to be transformed is the first index, vectorized
c over the second index. However, when vcpft is used, it is faster when the
c data is in transposed order, that is the data to be transformed is the
c second index, vectorized over the first index. Both version were kept.
c The difference in speed is around 15% using the Intel compiler on Linux
c and about a factor of 2 on the IBM SP.
c Using the dcft routines, the code runs about 30-40% faster than with vcpft.
c=======================================================================
      subroutine cosqx(a,w,c,nx,ny,isetup,isign)
      integer(ISZ):: nx,ny,isetup,isign
      real(kind=8):: a(-1:nx+1,-1:ny+1),w(2,0:nx,ny/2),c(nx-1)

c  This routine does Fourier transforms with only odd cosines.

c  This routine is based off the routines cosqi,cosqf1,cosqb1 fftpack of
c  the NETLIB collection.  The files were obtained in July of 1994.
c
c  Changes were made to the routines to allow vectorization over the
c  direction perpendicular to the transforms and to use vcpft to transform
c  two sets of data at once.  Also, the output of vcpft is organized
c  differently than the FFT routine used by NETLIB.  The loops were all
c  rearranged for optimization.  The FFT's are actually performed in the
c  scratch arrays.  This avoids an additional loop over the arrays to
c  copy the data back to the original array.

      real(kind=8):: pih,dt,fk
      integer(ISZ):: k,kc,ix,iy,ii,jj

      integer(ISZ):: naux1,naux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
 
c     --- Calculate how much space is needed.
      naux1 = 20000
      if (nx > 2048) naux1 = naux1 + 2.28*nx
      naux2 = 20000
      if (nx > 2048) naux2 = naux2 + 2.28*nx

c     --- Allocate the work space needed by dcft
c     --- Note that since this happens inside a OpenMP parallel construct, each thread
c     --- must allocate the work space each time.
      allocate(aux1(naux1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"cosqx: allocation error ",allocerror,
     &         ": could not allocate aux1 to shape ",naux1
        call kaboom("cosqx: allocation error")
        return
      endif
      allocate(aux2(naux2),stat=allocerror)
      if (allocerror /= 0) then
        print*,"cosqx: allocation error ",allocerror,
     &         ": could not allocate aux2 to shape ",naux2
        call kaboom("cosqx: allocation error")
        return
      endif

c     --- Initialize the auxiliary arrays
      call dcft(1_4,w,1_4,int(nx+1,4),w,int(1,4),int(nx+1,4),
     &          int(nx,4),int(ny/2,4),int(-isign,4),1.,
     &          aux1,int(naux1,4),aux2,int(naux2,4))

c     do setup
      if (isetup == 1) then
        pih = 2.*atan(1.)
        dt = pih/nx
        fk = 0.
        do k=1,nx/2-1
           fk = fk+1.
           c(k)    = cos(fk*dt) + cos((nx-fk)*dt)
           c(nx-k) = cos(fk*dt) - cos((nx-fk)*dt)
        enddo
        c(nx/2) = cos(nx*0.5*dt)
        return
      endif
   
c     do forward transform
      if (isign == -1) then

        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          do k=1,nx/2-1
            kc = nx-k
            w(ii,k ,jj)  = c(k)*a(k,iy) - c(kc)*a(kc,iy)
            w(ii,kc,jj) = c(kc)*a(k,iy) + c(k)*a(kc,iy)
          enddo
        enddo
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          w(ii,0   ,jj) = a(0,iy)
          w(ii,nx/2,jj) = c(nx/2)*(a(nx/2,iy) + a(nx/2,iy))
        enddo
        call dcft(0_4,w,1_4,int(nx+1,4),w,1_4,int(nx+1,4),
     &            int(nx,4),int(ny/2,4),int(+1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))
        call vrpft2(w(1,0,1), w(2,0,1), nx, 2,   ny/2,2*(nx+1))
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          a(0,iy) = w(ii,0,jj)
          a(nx-1,iy) = w(ii,nx/2,jj)
        enddo
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          do ix=1,nx/2-1
            a(2*ix-1,iy) = w(ii,ix,jj) - w(ii,nx-ix,jj)
            a(2*ix,iy)   = w(ii,ix,jj) + w(ii,nx-ix,jj)
          enddo
        enddo

c     do backward transform
      elseif (isign == 1) then

        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          w(ii,0   ,jj) = a(0,iy) + a(0,iy)
          w(ii,nx/2,jj) = a(nx-1,iy) + a(nx-1,iy)
        enddo
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          do ix=1,nx/2-1
            w(ii,ix   ,jj) = a(2*ix-1,iy) + a(2*ix,iy)
            w(ii,nx-ix,jj) = a(2*ix,iy) - a(2*ix-1,iy)
          enddo
        enddo
        call vrpfti2(w(1,0,1), w(2,0,1), nx, 2,   ny/2,2*(nx+1))
        call dcft(0_4,w,1_4,int(nx+1,4),w,1_4,int(nx+1,4),
     &            int(nx,4),int(ny/2,4),int(-1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          a(0,iy) = w(ii,0,jj)+w(ii,0,jj)
          a(nx/2,iy) = c(nx/2)*(w(ii,nx/2,jj) + w(ii,nx/2,jj))
        enddo
        do iy=0,ny-1
          ii = mod(iy,2) + 1
          jj = iy/2 + 1
          do k=1,nx/2-1
            kc = nx-k
            a(k,iy) = c(kc)*w(ii,kc,jj) + c(k)*w(ii,k,jj)
            a(kc,iy) = c(k)*w(ii,kc,jj) - c(kc)*w(ii,k,jj)
          enddo
        enddo

      endif

      deallocate(aux1)
      deallocate(aux2)

      return
      end
c=======================================================================
      subroutine cosqy(a,w,c,nx,ny,isetup,isign)
      integer(ISZ):: nx,ny,isetup,isign
      real(kind=8):: a(-1:nx+1,-1:ny+1),w(2,0:ny,nx/2),c(ny-1)

c  This routine does Fourier transforms with only odd cosines.

c  This routine is based off the routines cosqi,cosqf1,cosqb1 fftpack of
c  the NETLIB collection.  The files were obtained in July of 1994.
c
c  Changes were made to the routines to allow vectorization over the
c  direction perpendicular to the transforms and to use vcpft to transform
c  two sets of data at once.  Also, the output of vcpft is organized
c  differently than the FFT routine used by NETLIB.  The loops were all
c  rearranged for optimization.  The FFT's are actually performed in the
c  scratch arrays.  This avoids an additional loop over the arrays to
c  copy the data back to the original array.

      real(kind=8):: pih,dt,fk
      integer(ISZ):: k,kc,ix,iy,ii,jj

      integer(ISZ):: naux1,naux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      integer(ISZ):: allocerror
 
c Setup and make call to the ESSL routine for Fourier-transforms

c     --- Calculate how much space is needed.
      naux1 = 20000
      if (ny > 2048) naux1 = naux1 + 2.28*ny
      naux2 = 20000
      if (ny > 2048) naux2 = naux2 + 2.28*ny

c     --- Allocate the work space needed by dcft
      allocate(aux1(naux1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"cosqy: allocation error ",allocerror,
     &         ": could not allocate aux1 to shape ",naux1
        call kaboom("cosqy: allocation error")
        return
      endif
      allocate(aux2(naux2),stat=allocerror)
      if (allocerror /= 0) then
        print*,"cosqy: allocation error ",allocerror,
     &         ": could not allocate aux2 to shape ",naux2
        call kaboom("cosqy: allocation error")
        return
      endif

c     --- Initialize the auxiliary arrays
      call dcft(1_4,w,1_4,int(ny+1,4),w,1_4,int(ny+1,4),
     &          int(ny,4),int(nx/2,4),int(-isign,4),1.,
     &          aux1,int(naux1,4),aux2,int(naux2,4))

c     do setup
      if (isetup == 1) then
        pih = 2.*atan(1.)
        dt = pih/ny
        fk = 0.
        do k=1,ny/2-1
           fk = fk+1.
           c(k)    = cos(fk*dt) + cos((ny-fk)*dt)
           c(ny-k) = cos(fk*dt) - cos((ny-fk)*dt)
        enddo
        c(ny/2) = cos(ny*0.5*dt)
        return
      endif
   
c     do forward transform
      if (isign == -1) then

        do k=1,ny/2-1
          kc = ny-k
          do ix=0,nx-1
            ii = mod(ix,2) + 1
            jj = ix/2 + 1
            w(ii,k ,jj) = c(k)*a(ix,k) - c(kc)*a(ix,kc)
            w(ii,kc,jj) = c(kc)*a(ix,k) + c(k)*a(ix,kc)
          enddo
        enddo
        do ix=0,nx-1
          ii = mod(ix,2) + 1
          jj = ix/2 + 1
          w(ii,0,jj) = a(ix,0)
          w(ii,ny/2,jj) = c(ny/2)*(a(ix,ny/2) + a(ix,ny/2))
        enddo
        call dcft(0_4,w,1_4,int(ny+1,4),w,1_4,int(ny+1,4),
     &            int(ny,4),int(nx/2,4),int(+1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))
        call vrpft2(w(1,0,1), w(2,0,1), ny, 2,   nx/2,2*(ny+1))
        do ix=0,nx-1
          ii = mod(ix,2) + 1
          jj = ix/2 + 1
          a(ix,0) = w(ii,0,jj)
          a(ix,ny-1) = w(ii,ny/2,jj)
        enddo
        do iy=1,ny/2-1
          do ix=0,nx-1
            ii = mod(ix,2) + 1
            jj = ix/2 + 1
            a(ix,2*iy-1) = w(ii,iy,jj) - w(ii,ny-iy,jj)
            a(ix,2*iy)   = w(ii,iy,jj) + w(ii,ny-iy,jj)
          enddo
        enddo

c     do backward transform
      elseif (isign == 1) then

        do ix=0,nx-1
          ii = mod(ix,2) + 1
          jj = ix/2 + 1
          w(ii,0,jj) = a(ix,0) + a(ix,0)
          w(ii,ny/2,jj) = a(ix,ny-1) + a(ix,ny-1)
        enddo
        do iy=1,ny/2-1
          do ix=0,nx-1
            ii = mod(ix,2) + 1
            jj = ix/2 + 1
            w(ii,iy,jj) = a(ix,2*iy-1) + a(ix,2*iy)
            w(ii,ny-iy,jj) = a(ix,2*iy) - a(ix,2*iy-1)
          enddo
        enddo
        call vrpfti2(w(1,0,1), w(2,0,1), ny, 2,   nx/2,2*(ny+1))
        call dcft(0_4,w,1_4,int(ny+1,4),w,1_4,int(ny+1,4),
     &            int(ny,4),int(nx/2,4),int(-1,4),1.,
     &            aux1,int(naux1,4),aux2,int(naux2,4))
        do ix=0,nx-1
          ii = mod(ix,2) + 1
          jj = ix/2 + 1
          a(ix,0) = w(ii,0,jj)+w(ii,0,jj)
          a(ix,ny/2) = c(ny/2)*(w(ii,ny/2,jj)+w(ii,ny/2,jj))
        enddo
        do k=1,ny/2-1
          kc = ny-k
          do ix=0,nx-1
            ii = mod(ix,2) + 1
            jj = ix/2 + 1
            a(ix,k) = c(kc)*w(ii,kc,jj) + c(k)*w(ii,k,jj)
            a(ix,kc) = c(k)*w(ii,kc,jj) - c(kc)*w(ii,k,jj)
          enddo
        enddo

      endif

      deallocate(aux1)
      deallocate(aux2)

      return
      end
c===========================================================================
#else
c=======================================================================
c The following two routines use the vcpft routine for the FFT's. That
c routine runs faster when the data is in transposed order and the following
c versions were kept.
c The difference in speed is around 15% using the Intel compiler on Linux
c and about a factor of 2 on the IBM SP.
c=======================================================================
      subroutine cosqx(a,w,c,nx,ny,isetup,isign)
      integer(ISZ):: nx,ny,isetup,isign
      real(kind=8):: a(-1:nx+1,-1:ny+1),w(0:ny-1,0:nx),c(nx-1)

c  This routine does Fourier transforms with only odd cosines.

c  This routine is based off the routines cosqi,cosqf1,cosqb1 fftpack of
c  the NETLIB collection.  The files were obtained in July of 1994.
c
c  Changes were made to the routines to allow vectorization over the
c  direction perpendicular to the transforms and to use vcpft to transform
c  two sets of data at once.  Also, the output of vcpft is organized
c  differently than the FFT routine used by NETLIB.  The loops were all
c  rearranged for optimization.  The FFT's are actually performed in the
c  scratch arrays.  This avoids an additional loop over the arrays to
c  copy the data back to the original array.

      real(kind=8):: pih,dt,fk
      integer(ISZ):: k,kc,ix,iy

c     do setup
      if (isetup == 1) then
        pih = 2.*atan(1.)
        dt = pih/nx
        fk = 0.
        do k=1,nx/2-1
           fk = fk+1.
           c(k)    = cos(fk*dt) + cos((nx-fk)*dt)
           c(nx-k) = cos(fk*dt) - cos((nx-fk)*dt)
        enddo
        c(nx/2) = cos(nx*0.5*dt)
        return
      endif
   
c     do forward transform
      if (isign == -1) then

        do iy=0,ny-1
          do k=1,nx/2-1
            kc = nx-k
            w(iy,k)  = c(k)*a(k,iy) - c(kc)*a(kc,iy)
            w(iy,kc) = c(kc)*a(k,iy) + c(k)*a(kc,iy)
          enddo
        enddo
        do iy=0,ny-1
          w(iy,0) = a(0,iy)
          w(iy,nx/2) = c(nx/2)*(a(nx/2,iy) + a(nx/2,iy))
        enddo
        call vcpft (w(0,0), w(1,0), nx, ny,-1,ny/2,2)
        call vrpft2(w(0,0), w(1,0), nx, ny,   ny/2,2)
        do iy=0,ny-1
          a(0,iy) = w(iy,0)
          a(nx-1,iy) = w(iy,nx/2)
        enddo
        do iy=0,ny-1
          do ix=1,nx/2-1
            a(2*ix-1,iy) = w(iy,ix) - w(iy,nx-ix)
            a(2*ix,iy)   = w(iy,ix) + w(iy,nx-ix)
          enddo
        enddo

c     do backward transform
      elseif (isign == 1) then

        do iy=0,ny-1
          w(iy,0) = a(0,iy) + a(0,iy)
          w(iy,nx/2) = a(nx-1,iy) + a(nx-1,iy)
        enddo
        do iy=0,ny-1
          do ix=1,nx/2-1
            w(iy,ix) = a(2*ix-1,iy) + a(2*ix,iy)
            w(iy,nx-ix) = a(2*ix,iy) - a(2*ix-1,iy)
          enddo
        enddo
        call vrpfti2(w(0,0), w(1,0), nx, ny,   ny/2,2)
        call vcpft  (w(0,0), w(1,0), nx, ny,+1,ny/2,2)
        do iy=0,ny-1
          a(0,iy) = w(iy,0)+w(iy,0)
          a(nx/2,iy) = c(nx/2)*(w(iy,nx/2) + w(iy,nx/2))
        enddo
        do iy=0,ny-1
          do k=1,nx/2-1
            kc = nx-k
            a(k,iy) = c(kc)*w(iy,kc) + c(k)*w(iy,k)
            a(kc,iy) = c(k)*w(iy,kc) - c(kc)*w(iy,k)
          enddo
        enddo

      endif

      return
      end
c=======================================================================
      subroutine cosqy(a,w,c,nx,ny,isetup,isign)
      integer(ISZ):: nx,ny,isetup,isign
      real(kind=8):: a(-1:nx+1,-1:ny+1),w(0:nx-1,0:ny),c(ny-1)

c  This routine does Fourier transforms with only odd cosines.

c  This routine is based off the routines cosqi,cosqf1,cosqb1 fftpack of
c  the NETLIB collection.  The files were obtained in July of 1994.
c
c  Changes were made to the routines to allow vectorization over the
c  direction perpendicular to the transforms and to use vcpft to transform
c  two sets of data at once.  Also, the output of vcpft is organized
c  differently than the FFT routine used by NETLIB.  The loops were all
c  rearranged for optimization.  The FFT's are actually performed in the
c  scratch arrays.  This avoids an additional loop over the arrays to
c  copy the data back to the original array.

      real(kind=8):: pih,dt,fk
      integer(ISZ):: k,kc,ix,iy

c     do setup
      if (isetup == 1) then
        pih = 2.*atan(1.)
        dt = pih/ny
        fk = 0.
        do k=1,ny/2-1
           fk = fk+1.
           c(k)    = cos(fk*dt) + cos((ny-fk)*dt)
           c(ny-k) = cos(fk*dt) - cos((ny-fk)*dt)
        enddo
        c(ny/2) = cos(ny*0.5*dt)
        return
      endif
   
c     do forward transform
      if (isign == -1) then

        do k=1,ny/2-1
          kc = ny-k
          do ix=0,nx-1
            w(ix,k)  = c(k)*a(ix,k) - c(kc)*a(ix,kc)
            w(ix,kc) = c(kc)*a(ix,k) + c(k)*a(ix,kc)
          enddo
        enddo
        do ix=0,nx-1
          w(ix,0) = a(ix,0)
          w(ix,ny/2) = c(ny/2)*(a(ix,ny/2) + a(ix,ny/2))
        enddo
        call vcpft (w(0,0), w(1,0), ny, nx,-1,nx/2,2)
        call vrpft2(w(0,0), w(1,0), ny, nx,   nx/2,2)
        do ix=0,nx-1
          a(ix,0) = w(ix,0)
          a(ix,ny-1) = w(ix,ny/2)
        enddo
        do iy=1,ny/2-1
          do ix=0,nx-1
            a(ix,2*iy-1) = w(ix,iy) - w(ix,ny-iy)
            a(ix,2*iy)   = w(ix,iy) + w(ix,ny-iy)
          enddo
        enddo

c     do backward transform
      elseif (isign == 1) then

        do ix=0,nx-1
          w(ix,0) = a(ix,0) + a(ix,0)
          w(ix,ny/2) = a(ix,ny-1) + a(ix,ny-1)
        enddo
        do iy=1,ny/2-1
          do ix=0,nx-1
            w(ix,iy) = a(ix,2*iy-1) + a(ix,2*iy)
            w(ix,ny-iy) = a(ix,2*iy) - a(ix,2*iy-1)
          enddo
        enddo
        call vrpfti2(w(0,0), w(1,0), ny, nx,   nx/2,2)
        call vcpft  (w(0,0), w(1,0), ny, nx,+1,nx/2,2)
        do ix=0,nx-1
          a(ix,0) = w(ix,0)+w(ix,0)
          a(ix,ny/2) = c(ny/2)*(w(ix,ny/2)+w(ix,ny/2))
        enddo
        do k=1,ny/2-1
          kc = ny-k
          do ix=0,nx-1
            a(ix,k) = c(kc)*w(ix,kc) + c(k)*w(ix,k)
            a(ix,kc) = c(k)*w(ix,kc) - c(kc)*w(ix,k)
          enddo
        enddo

      endif

      return
      end
c===========================================================================
#endif
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c                                                                          =
c General capacity matrix solver                                           =
c                                                                          =
c===========================================================================
c=======================================================================
      subroutine cmset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                   xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                   dx,dy,dz,xmmin,ymmin,zmmin,
     &                   scrtch,xywork,zwork,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
      use CapMat3d
      integer(ISZ):: nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Sets up capacity matrix for 3d fieldsolver. Note that the capacity matrix
c is calculated using field solves which do not do any filtering.
c (Experience has shown that when filtering is used, the inverse capacity
c matrix is very ill-conditioned and so the capacity matrix becomes garbled.)

      real(kind=8),allocatable:: zerofilt(:,:)
      integer(ISZ):: allocerror
      real(kind=8):: wx,wy,wz,cx,cy
c     real(kind=8):: tt1,tt2
      integer(ISZ):: ix,iy,iz,i,j,info

c find capacity matrix

      if (nc3d == 0) return

c     --- find inverse capacity matrix
      do i=1,nc3d

c       --- start with zero phi
c       call zeroarry(phi,(nx+1)*(ny+1)*(1+nzlocal))
        phi = 0.

c       --- set phi with 1 Coulomb at conductor points (charge is scaled by
c       --- 2 or 4 at a symmetric boundary)
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        iz = int((zcond3d(i)-zmmin)/dz)
        wz =     (zcond3d(i)-zmmin)/dz  - iz
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy == 0) cy = 2.
        if (l4symtry .and. ix == 0) cx = 2.
        if (l4symtry .and. iy == 0) cy = 2.
        phi(ix  ,iy  ,iz  ) = (1.-wx)*(1.-wy)*(1.-wz)*cx*cy
        phi(ix+1,iy  ,iz  ) =     wx *(1.-wy)*(1.-wz)*cy
        phi(ix  ,iy+1,iz  ) = (1.-wx)*    wy *(1.-wz)*cx
        phi(ix+1,iy+1,iz  ) =     wx *    wy *(1.-wz)
        phi(ix  ,iy  ,iz+1) = (1.-wx)*(1.-wy)*    wz *cx*cy
        phi(ix+1,iy  ,iz+1) =     wx *(1.-wy)*    wz *cy
        phi(ix  ,iy+1,iz+1) = (1.-wx)*    wy *    wz *cx
        phi(ix+1,iy+1,iz+1) =     wx *    wy *    wz 

c       --- Solve for fields, with the filtering turned off
        allocate(zerofilt(size(filt,1),size(filt,2)),stat=allocerror)
        if (allocerror /= 0) then
          print*,"cmset3d: allocation error ",allocerror,
     &           ": could not allocate zerofilt to shape ",
     &           size(filt,1),size(filt,2)
          call kaboom("cmset3d: allocation error")
          return
        endif
        zerofilt = 0.
        call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &                  zerofilt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                  scrtch,xywork,zwork,0,l2symtry,l4symtry,
     &                  bound0,boundnz,boundxy)
        deallocate(zerofilt)

c       --- fill matrix with phi
        do j=1,nc3d
          ix = int((xcond3d(j)-xmmin)/dx)
          wx =     (xcond3d(j)-xmmin)/dx  - ix
          iy = int((ycond3d(j)-ymmin)/dy)
          wy =     (ycond3d(j)-ymmin)/dy  - iy
          iz = int((zcond3d(j)-zmmin)/dz)
          wz =     (zcond3d(j)-zmmin)/dz  - iz
          cmat3d(j,i,0) = phi(ix  ,iy  ,iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &                    phi(ix+1,iy  ,iz  )*    wx *(1.-wy)*(1.-wz) +
     &                    phi(ix  ,iy+1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &                    phi(ix+1,iy+1,iz  )*    wx *    wy *(1.-wz) +
     &                    phi(ix  ,iy  ,iz+1)*(1.-wx)*(1.-wy)*    wz  +
     &                    phi(ix+1,iy  ,iz+1)*    wx *(1.-wy)*    wz  +
     &                    phi(ix  ,iy+1,iz+1)*(1.-wx)*    wy *    wz  +
     &                    phi(ix+1,iy+1,iz+1)*    wx *    wy *    wz 
        enddo

      enddo

c     --- invert to get capacity matrix (checking for singular matrix)

c     --- These routines don't seem to give as good a result as those below.
c     call dsifa(cmat3d,nc3d,nc3d,kpvt3d,info)
c     if (info .ne. 0) then
c       call kaboom("cmset3d: ERROR: Capacity matrix is singular - it is likely that
c    & two or more conductor points are within the same grid cell.")
c     endif
c     call dsidi(cmat3d,nc3d,nc3d,kpvt3d,tt1,tt2,scrtch,001)

c     --- These routines are from LAPACK.
#if WORDSIZE == 64
      call ssytrf("U",nc3d,cmat3d,nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),info)
#else
#ifdef CYGWIN
      call dsytrf_("U",nc3d,cmat3d,nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),info)
#else
      call dsytrf("U",nc3d,cmat3d,nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),info)
#endif
#endif
      if (info .ne. 0) then
c       print*,info
        call kaboom("cmset3d: ERROR: Capacity matrix is singular - it is likely that
     & two or more conductor points are within the same grid cell.")
      endif
c     --- pcond3d is passed as workspace.
#if WORDSIZE == 64
      call ssytri("U",nc3d,cmat3d,nc3d,kpvt3d,pcond3d,info)
#else
#ifdef CYGWIN
      call dsytri_("U",nc3d,cmat3d,nc3d,kpvt3d,pcond3d,info)
#else
      call dsytri("U",nc3d,cmat3d,nc3d,kpvt3d,pcond3d,info)
#endif
#endif

c     --- Fill in lower half (which is the same as the upper half)
      do j=2,nc3d
        do i=1,j-1
          cmat3d(j,i,0) = cmat3d(i,j,0)
        enddo
      enddo

      return
      end
c============================================================================== 
      subroutine capmat3df(iwhich,phi,rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                    filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,dx,dy,dz,
     &                    xmmin,ymmin,zmmin,scrtch,xywork,zwork,
     &                    l2symtry,l4symtry,bound0,boundnz,boundxy)
      use CapMat3d
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Applies the capacity matrix for 3d code electrostatic field solver.
c Multiplies phi by cap to get induced rho.
c Assumes that beam's rho has been set, and copied into phi.
c When filtering is being done, field solve is inexact since filtering is
c applied to the induced charge.


      real(kind=8):: wx,wy,wz,cx,cy
      real(kind=8):: wx0,wy0,wz0
      integer(ISZ):: ix,iy,iz,i,j

c     --- Initialize the arrays for poisson solve and the capacity matrix.
      if (iwhich == 1 .or. iwhich == 0) then

        call vpois3d(1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &               scrtch,xywork,zwork,0,l2symtry,l4symtry,
     &               bound0,boundnz,boundxy)

        call cmset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &               xlen,ylen,zlen,nx,ny,nzlocal,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &               scrtch,xywork,zwork,l2symtry,l4symtry,
     &               bound0,boundnz,boundxy)
        if (iwhich == 1) then
          phi(0:nx,0:ny,0:nzlocal) = rho
        endif

      elseif (iwhich > 1) then

c       --- Make call to do the specialized action requested.
        call vpois3d(iwhich,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &               scrtch,xywork,zwork,0,l2symtry,l4symtry,
     &               bound0,boundnz,boundxy)

      endif

c     --- If a full field solve was not desired, return.
      if (iwhich > 0) return

c     --- Do the first field solve (note that rho must have been copied into
c     --- phi for this to work)
      call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &             scrtch,xywork,zwork,0,l2symtry,l4symtry,
     &             bound0,boundnz,boundxy)

c     --- Extract phi error from FFT solution and zero qcond
      do i=1,nc3d
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        iz = int((zcond3d(i)-zmmin)/dz)
        wz =     (zcond3d(i)-zmmin)/dz  - iz
        pcond3d(i,0)=vcond3d(i)-(phi(ix  ,iy  ,iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &                           phi(ix+1,iy  ,iz  )*    wx *(1.-wy)*(1.-wz) +
     &                           phi(ix  ,iy+1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &                           phi(ix+1,iy+1,iz  )*    wx *    wy *(1.-wz) +
     &                           phi(ix  ,iy  ,iz+1)*(1.-wx)*(1.-wy)*    wz  +
     &                           phi(ix+1,iy  ,iz+1)*    wx *(1.-wy)*    wz  +
     &                           phi(ix  ,iy+1,iz+1)*(1.-wx)*    wy *    wz  +
     &                           phi(ix+1,iy+1,iz+1)*    wx *    wy *    wz )
        qcond3d(i,0) = 0.
      enddo

c     --- Multiply pcond by capacity matrix to get induced charge
      do j=1,nc3d
        do i=1,nc3d
          qcond3d(i,0) = qcond3d(i,0) + pcond3d(j,0)*cmat3d(i,j,0)
        enddo
      enddo

c     --- Recopy rho into phi.
      phi(0:nx,0:ny,0:nzlocal) = rho

c     --- Deposit induced charge onto grid
      do i=1,nc3d
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        wx0 = 1. - wx
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        wy0 = 1. - wy
        iz = int((zcond3d(i)-zmmin)/dz)
        wz =     (zcond3d(i)-zmmin)/dz  - iz
        wz0 = 1. - wz
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy == 0) cy = 2.
        if (l4symtry .and. ix == 0) cx = 2.
        if (l4symtry .and. iy == 0) cy = 2.
        phi(ix  ,iy  ,iz  )=phi(ix  ,iy  ,iz  )+qcond3d(i,0)*wx0*wy0*wz0*cx*cy
        phi(ix+1,iy  ,iz  )=phi(ix+1,iy  ,iz  )+qcond3d(i,0)*wx *wy0*wz0*cy
        phi(ix  ,iy+1,iz  )=phi(ix  ,iy+1,iz  )+qcond3d(i,0)*wx0*wy *wz0*cx
        phi(ix+1,iy+1,iz  )=phi(ix+1,iy+1,iz  )+qcond3d(i,0)*wx *wy *wz0
        phi(ix  ,iy  ,iz+1)=phi(ix  ,iy  ,iz+1)+qcond3d(i,0)*wx0*wy0*wz *cx*cy
        phi(ix+1,iy  ,iz+1)=phi(ix+1,iy  ,iz+1)+qcond3d(i,0)*wx *wy0*wz *cy
        phi(ix  ,iy+1,iz+1)=phi(ix  ,iy+1,iz+1)+qcond3d(i,0)*wx0*wy *wz *cx
        phi(ix+1,iy+1,iz+1)=phi(ix+1,iy+1,iz+1)+qcond3d(i,0)*wx *wy *wz
      enddo

c     --- Redo the field solve, now including both the beam space-charge and
c     --- the induced charge in the conductor.
      call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

      return
      end
c============================================================================== 
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c                                                                          =
c General capacity matrix solver in kz space                               =
c                                                                          =
c===========================================================================
c=======================================================================
      subroutine cmkzset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                     xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                     dx,dy,dz,xmmin,ymmin,zmmin,
     &                     scrtch,xywork,zwork,l2symtry,l4symtry,
     &                     bound0,boundnz,boundxy)
      use CapMat3d
      integer(ISZ):: nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Sets up capacity matrix for 3d fieldsolver. Note that the capacity matrix
c is calculated using field solves which do not do any filtering.
c (Experience has shown that when filtering is used, the inverse capacity
c matrix is very ill-conditioned and so the capacity matrix becomes garbled.)


      real(kind=8):: wx,wy,cx,cy
c     real(kind=8):: tt1,tt2
      integer(ISZ):: ix,iy,iz,i,j,info

c find capacity matrix

      if (nc3d == 0) return

c     --- Make sure capacity matrix is fully allocated.
      nc3dz = nzlocal
#ifdef MPIPARALLEL
      nc3dz2 = nzlocal - 1
#else
      nc3dz2 = nzlocal/2
#endif
      call gchange("CapMat3d",0)

c     --- find inverse capacity matrix
      do i=1,nc3d

c       --- start with zero phi
c       call zeroarry(phi,(nx+1)*(ny+1)*(1+nzlocal))
        phi = 0.

c       --- set phi with 1 Coulomb at conductor points (charge is scaled by
c       --- 2 or 4 at a symmetric boundary)
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy == 0) cy = 2.
        if (l4symtry .and. ix == 0) cx = 2.
        if (l4symtry .and. iy == 0) cy = 2.
        do iz=0,nc3dz2
          phi(ix  ,iy  ,iz) = (1.-wx)*(1.-wy)*cx*cy
          phi(ix+1,iy  ,iz) =     wx *(1.-wy)*cy
          phi(ix  ,iy+1,iz) = (1.-wx)*    wy *cx
          phi(ix+1,iy+1,iz) =     wx *    wy
        enddo

c       --- Solve for fields, transforming only in x and y and skipping
c       --- filtering.
        call vpois3d(12,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)
        call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)
        call vpois3d(13,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)

c       --- fill matrix with phi
        do j=1,nc3d
          ix = int((xcond3d(j)-xmmin)/dx)
          wx =     (xcond3d(j)-xmmin)/dx  - ix
          iy = int((ycond3d(j)-ymmin)/dy)
          wy =     (ycond3d(j)-ymmin)/dy  - iy
          do iz=0,nc3dz2
            cmat3d(j,i,iz) = phi(ix  ,iy  ,iz)*(1.-wx)*(1.-wy) +
     &                       phi(ix+1,iy  ,iz)*    wx *(1.-wy) +
     &                       phi(ix  ,iy+1,iz)*(1.-wx)*    wy  +
     &                       phi(ix+1,iy+1,iz)*    wx *    wy
          enddo
        enddo

      enddo

c     --- invert to get capacity matrices (checking for singular matrix)
      do iz=0,nc3dz2

c       --- These routines don't seem to give as good a result as those below.
c       call dsifa(cmat3d(1,1,iz),nc3d,nc3d,kpvt3d,info)
c       if (info .ne. 0) then
c         call kaboom("cmkzset3d: ERROR: Capacity matrix is singular - it is likely that
c    & two or more conductor points are within the same grid cell.")
c       endif
c       call dsidi(cmat3d(1,1,iz),nc3d,nc3d,kpvt3d,tt1,tt2,scrtch,001)

c       --- These routines are from LAPACK.
#if WORDSIZE == 64
        call ssytrf("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),
     &              info)
#else
#ifdef CYGWIN
        call dsytrf_("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),
     &              info)
#else
        call dsytrf("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),
     &              info)
#endif
#endif
        if (info .ne. 0) then
c         print*,info
          call kaboom("cmkzset3d: ERROR: Capacity matrix is singular - it is likely that
     & two or more conductor points are within the same grid cell.")
        endif
c       --- pcond3d is passed as workspace.
#if WORDSIZE == 64
        call ssytri("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,pcond3d,info)
#else
#ifdef CYGWIN
        call dsytri_("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,pcond3d,info)
#else
        call dsytri("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,pcond3d,info)
#endif
#endif

c     --- Fill in lower half (which is the same as the upper half)
        do j=2,nc3d
          do i=1,j-1
            cmat3d(j,i,iz) = cmat3d(i,j,iz)
          enddo
        enddo

c     --- End of loop over nc3dz
      enddo

      return
      end
c============================================================================== 
      subroutine capmatkz3d(iwhich,phi,rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                      filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,dx,dy,dz,
     &                      xmmin,ymmin,zmmin,scrtch,xywork,zwork,
     &                      l2symtry,l4symtry,bound0,boundnz,boundxy)
      use CapMat3d
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Applies the capacity matrix for 3d code electrostatic field solver.
c Multiplies phi by cap to get induced rho.
c Assumes that beam's rho has been set, and copied into phi.
c When filtering is being done, field solve is inexact since filtering is
c applied to the induced charge.
c For the serial version, advantage is taken by realizing the the capacity
c matrix for the kz plane iz is the same as that for kz plane nzlocal-iz. The
c amount of storage then for the matrix then is halved.
c For the parallel version though, the matrices for all planes is calculated
c for simplicity. This avoids having to pass the matrices around, and does not
c slow the initialization down any. (Granted that the matrix inversion time
c could be cut is half if careful load balancing were done.)


      real(kind=8):: wx,wy,cx,cy
      real(kind=8):: wx0,wy0
      integer(ISZ):: ix,iy,iz,i,j,iiz

c     --- Initialize the arrays for poisson solve and the capacity matrix.
      if (iwhich == 1 .or. iwhich == 0) then

        call vpois3d(1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &              filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &              l2symtry,l4symtry,bound0,boundnz,boundxy)

        call cmkzset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                 xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                 dx,dy,dz,xmmin,ymmin,zmmin,
     &                 scrtch,xywork,zwork,l2symtry,l4symtry,
     &                 bound0,boundnz,boundxy)
        if (iwhich == 1) then
          phi(0:nx,0:ny,0:nzlocal) = rho
        endif

      elseif (iwhich > 1) then

c       --- Make call to do the specialized action requested.
        call vpois3d(iwhich,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &               scrtch,xywork,zwork,0,l2symtry,l4symtry,
     &               bound0,boundnz,boundxy)

      endif

c     --- If a full field solve was not desired, return.
      if (iwhich > 0) return

c     --- Do the first field solve (note that rho must have been copied into
c     --- phi for this to work). The inverse transform in z is skipped.
      call vpois3d( 2,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 3,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d(13,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

c     --- Extract phi error from FFT solution and zero qcond3d
      do i=1,nc3d
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        do iz=0,nc3dz-1
          pcond3d(i,iz) = -(phi(ix  ,iy  ,iz)*(1.-wx)*(1.-wy) +
     &                      phi(ix+1,iy  ,iz)*    wx *(1.-wy) +
     &                      phi(ix  ,iy+1,iz)*(1.-wx)*    wy  +
     &                      phi(ix+1,iy+1,iz)*    wx *    wy)
          qcond3d(i,iz) = 0.
        enddo
      enddo

c     --- Multiply pcond3d by capacity matrix to get induced charge.
c     --- Note the clever setting of iiz. For the serial code, this does
c     --- the correct thing, since the capacity matrix for plane nzlocal-iz is the
c     --- same as the matrix for plane iz.  For the parallel code, it also
c     --- does the correct thing, i.e. nothing, since iiz is always < nc3dz2.
      do iz=0,nc3dz-1
        iiz = iz
        if (iiz > nc3dz2) iiz = nc3dz - iz
        do j=1,nc3d
          do i=1,nc3d
            qcond3d(i,iz) = qcond3d(i,iz) + pcond3d(j,iz)*cmat3d(i,j,iiz)
          enddo
        enddo
      enddo

c     --- Recopy rho into phi.
      phi(0:nx,0:ny,0:nzlocal) = rho

c     --- Transform in z
      call vpois3d(10,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

c     --- Deposit induced charge onto grid
      do i=1,nc3d
        ix = int((xcond3d(i)-xmmin)/dx)
        wx =     (xcond3d(i)-xmmin)/dx  - ix
        wx0 = 1. - wx
        iy = int((ycond3d(i)-ymmin)/dy)
        wy =     (ycond3d(i)-ymmin)/dy  - iy
        wy0 = 1. - wy
        cx = 1.
        cy = 1.
        if (l2symtry .and. iy == 0) cy = 2.
        if (l4symtry .and. ix == 0) cx = 2.
        if (l4symtry .and. iy == 0) cy = 2.
        do iz=0,nc3dz-1
          phi(ix  ,iy  ,iz) = phi(ix  ,iy  ,iz) + qcond3d(i,iz)*wx0*wy0*cx*cy
          phi(ix+1,iy  ,iz) = phi(ix+1,iy  ,iz) + qcond3d(i,iz)*wx *wy0*cy
          phi(ix  ,iy+1,iz) = phi(ix  ,iy+1,iz) + qcond3d(i,iz)*wx0*wy *cx
          phi(ix+1,iy+1,iz) = phi(ix+1,iy+1,iz) + qcond3d(i,iz)*wx *wy
        enddo
      enddo

c     --- Finish the field solve, now including both the beam space-charge and
c     --- the induced charge in the conductor.
      call vpois3d(12,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 5,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

      return
      end
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
c                                                                          =
c Capacity matrix solver with infinite axial pipe: the capacity matrix is  =
c applied in kz space                                                      =
c                                                                          =
c===========================================================================
      subroutine pipe3df(iwhich,pipeshpe,rho,phi,kxsq,kysq,kzsq,attx,atty,attz, 
     &                  filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,
     &                  l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Pipe3d
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      character(*):: pipeshpe
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy
 
c  External interface to capacity matrix field solver.
c  Passes external and local variables to field solver routines.
c  If iwhich is '0' or '1', also allocates local variables


c  Allocate local variables
      pipenz = nzlocal
#ifdef MPIPARALLEL
      pipenz2 = nzlocal - 1
#else
      pipenz2 = nzlocal/2
#endif
      if (iwhich == 0 .or. iwhich == 1) then
        if (pipen == 0) pipen = 2*nx
        pipe8th = pipen/8 + 1
        call gallot("Pipe3d",0)
      endif

c  Call field solver
      call vpipe3d(iwhich,pipeshpe,rho,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,
     &             piperadi,pipen,pipenz,pipenz2,pipe8th,pipex,pipey,
     &             pipephi,pipeq,cap3d,kpvt,l2symtry,l4symtry,
     &             bound0,boundnz,boundxy)

      return
      end
c=======================================================================
      subroutine pipest3d(pipeshpe,cap,phi,kxsq,kysq,kzsq,attx,atty,attz, 
     &                    filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                    scrtch,xywork,zwork,
     &                    piperadi,pipen,pipenz2,pipe8th,pipex,pipey,cap3d,kpvt,
     &                    l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Constant
      integer(ISZ):: nx,ny,nzlocal,nz
      character(*):: pipeshpe
      real(kind=8):: cap(2*(nx+ny),*)
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen, piperadi
      integer(ISZ):: pipen,pipenz2,pipe8th
      real(kind=8):: pipex(pipen+1), pipey(pipen+1), cap3d(pipen,pipe8th,0:pipenz2)
      integer(ISZ):: kpvt(pipen)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c sets up stuff for capacity matrix
c CAP is a work array, it must be at least 4*(nx+ny)**2 in size
c PHI can be used since both CAP and PHI are not used at the same time


      real(kind=8):: yymax,xx,yy,tt1,tt2
      integer(ISZ):: ix,iy,iz,i,j,info

      if (piperadi == 0) piperadi = nx/2 - 1.

c Find points on first eigth of pipe
      if (pipeshpe == "hyp") then
c --- points on the quadrupole, i.e., on a hyperbola
        if (piperadi == 0) piperadi = nx/3
        yymax = sqrt((nx*0.5-1.)**2 - piperadi**2)
        do i=1,pipe8th
          pipey(i) = yymax*((i-1.)/(pipe8th-1.))**2 + ny*0.5
          pipex(i) = sqrt(piperadi**2 + (pipey(i) - ny*0.5)**2) + nx*0.5
        enddo
      else
c --- circle is the default
        if (piperadi == 0) piperadi = nx*0.5- 1.
        do i=1,pipe8th
          pipex(i) = piperadi*cos((i-1.)*pi*0.25/(pipe8th-1.)) + nx*0.5
          pipey(i) = piperadi*sin((i-1.)*pi*0.25/(pipe8th-1.)) + ny*0.5
        enddo
      endif

c     --- now use that info to find all of the points.
c     --- this insures eight fold symmetry.
      do i=2,pipe8th
        pipex(i+1*pipe8th-1) =      pipey(pipe8th-i+1)
        pipey(i+1*pipe8th-1) =      pipex(pipe8th-i+1)
        pipex(i+2*pipe8th-2) = ny - pipey(i)
        pipey(i+2*pipe8th-2) =      pipex(i)
        pipex(i+3*pipe8th-3) = nx - pipex(pipe8th-i+1)
        pipey(i+3*pipe8th-3) =      pipey(pipe8th-i+1)
        pipex(i+4*pipe8th-4) = nx - pipex(i)
        pipey(i+4*pipe8th-4) = ny - pipey(i)
        pipex(i+5*pipe8th-5) = ny - pipey(pipe8th-i+1)
        pipey(i+5*pipe8th-5) = nx - pipex(pipe8th-i+1)
        pipex(i+6*pipe8th-6) =      pipey(i)
        pipey(i+6*pipe8th-6) = nx - pipex(i)
        pipex(i+7*pipe8th-7) =      pipex(pipe8th-i+1)
        pipey(i+7*pipe8th-7) = ny - pipey(pipe8th-i+1)
      enddo

c find capacity matrix for each kz
c     --- find inverse capacity matrix for all kz
c     --- calculate stuff for first eighth (using phi for convenience)
c     --- storing it in cap3d temporarily
      do i=1,pipe8th
c       call zeroarry(phi,(nx+1)*(ny+1)*(nzlocal+1))
        phi = 0.
c       --- set phi with z transform of sinusoidal line charge for all kz
        ix = pipex(i)
        xx = pipex(i) - ix
        iy = pipey(i)
        yy = pipey(i) - iy
        do iz=0,pipenz2
          phi(ix  ,iy  ,iz) = (1.-xx)*(1.-yy)
          phi(ix+1,iy  ,iz) =     xx *(1.-yy)
          phi(ix  ,iy+1,iz) = (1.-xx)*    yy 
          phi(ix+1,iy+1,iz) =     xx *    yy 
        enddo
c       --- do transform in x and y
        call vpois3d(12,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,
     &               xywork,zwork,0,l2symtry,l4symtry,bound0,boundnz,boundxy)
c       --- divide by k squared
        call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,
     &               xywork,zwork,0,l2symtry,l4symtry,bound0,boundnz,boundxy)
c       --- inverse transform in x and y to find phi(kz)
        call vpois3d(13,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,
     &               xywork,zwork,0,l2symtry,l4symtry,bound0,boundnz,boundxy)
c       --- fill matrix with phi(kz)
        do j=1,pipen
          ix = pipex(j)
          xx = pipex(j) - ix
          iy = pipey(j)
          yy = pipey(j) - iy
          do iz=0,pipenz2
            cap3d(j,i,iz) = phi(ix  ,iy  ,iz) * (1.-xx)*(1.-yy) +
     &                      phi(ix+1,iy  ,iz) *     xx *(1.-yy) +
     &                      phi(ix  ,iy+1,iz) * (1.-xx)*    yy  +
     &                      phi(ix+1,iy+1,iz) *     xx *    yy 
          enddo
        enddo
      enddo

c     --- Now find the capacity matrix one kz at a time
      do iz=0,pipenz2
c       --- fill cap with inverse capacity matrix
        do i=0,pipe8th-2
          do j=1,pipen
            cap(j,i+1)         = cap3d(j,i+1,iz)
          enddo
          do j=1,2*pipe8th-1
            cap(j,i+1*pipe8th  ) = cap3d(2*pipe8th-j,pipe8th-i,iz)
          enddo
          do j=2*pipe8th,pipen
            cap(j,i+1*pipe8th  ) = cap3d(pipen-j+2*pipe8th,pipe8th-i,iz)
          enddo
          do j=1,2*pipe8th-2
            cap(j,i+2*pipe8th-1) = cap3d(pipen-2*pipe8th+2+j,i+1,iz)
          enddo
          do j=2*pipe8th-1,pipen
            cap(j,i+2*pipe8th-1) = cap3d(j-2*pipe8th+2,i+1,iz)
          enddo
          do j=1,4*pipe8th-3
            cap(j,i+3*pipe8th-2) = cap3d(4*pipe8th-2-j,pipe8th-i,iz)
          enddo
          do j=4*pipe8th-2,pipen
            cap(j,i+3*pipe8th-2) = cap3d(pipen+4*pipe8th-2-j,pipe8th-i,iz)
          enddo
          do j=1,4*pipe8th-4
            cap(j,i+4*pipe8th-3) = cap3d(pipen-4*pipe8th+4+j,i+1,iz)
          enddo
          do j=4*pipe8th-3,pipen
            cap(j,i+4*pipe8th-3) = cap3d(j-4*pipe8th+4,i+1,iz)
          enddo
          do j=1,6*pipe8th-5
            cap(j,i+5*pipe8th-4) = cap3d(6*pipe8th-4-j,pipe8th-i,iz)
          enddo
          do j=6*pipe8th-4,pipen
            cap(j,i+5*pipe8th-4) = cap3d(pipen+6*pipe8th-4-j,pipe8th-i,iz)
          enddo
          do j=1,6*pipe8th-6
            cap(j,i+6*pipe8th-5) = cap3d(pipen-6*pipe8th+6+j,i+1,iz)
          enddo
          do j=6*pipe8th-5,pipen
            cap(j,i+6*pipe8th-5) = cap3d(j-6*pipe8th+6,i+1,iz)
          enddo
          cap(1,i+7*pipe8th-6) = cap3d(1,pipe8th-i,iz)
          do j=2,8*pipe8th-8
            cap(j,i+7*pipe8th-6) = cap3d(8*pipe8th-6-j,pipe8th-i,iz)
          enddo
        enddo

c       --- invert to get capacity matrix
        call dsifa(cap,2*(nx+ny),pipen,kpvt,info)
        call dsidi(cap,2*(nx+ny),pipen,kpvt,tt1,tt2,scrtch,001)

c       --- store matrix in 3-D array (cap only contains upper half)
        do i=1,pipe8th
          do j=i,pipen
            cap3d(j,i,iz) = cap(i,j)
          enddo
        enddo
        do i=1,pipe8th-1
          do j=i+1,pipe8th
            cap3d(i,j,iz) = cap(i,j)
          enddo
        enddo

      enddo

      return
      end
c============================================================================== 
      subroutine vpipe3d(iwhich,pipeshpe,rho,phi,kxsq,kysq,kzsq,attx,atty,attz, 
     &                   filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,
     &                   piperadi,pipen,pipenz,pipenz2,pipe8th,pipex,pipey,
     &                   pipephi,pipeq,cap3d,kpvt,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      character(*):: pipeshpe
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(2,0:nx,0:nz)
      real(kind=8):: xlen,ylen,zlen,piperadi
      integer(ISZ):: pipen,pipenz,pipenz2,pipe8th
      real(kind=8):: pipex(pipen+1),pipey(pipen+1),cap3d(pipen,pipe8th,0:pipenz2)
      integer(ISZ):: kpvt(pipen)
      real(kind=8):: pipephi(pipen,0:pipenz-1),pipeq(pipen,0:pipenz-1)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c applies capacity matrix, for round pipe electrostatic solver.
c multiplies phi by cap to get rho.
c assumes that rho has been set, and copied into phi.

      real(kind=8):: xx,yy
      integer(ISZ):: ix,iy,iz,iiz,i,j

c initialize capacity matrix and arrays for poisson solve
      if (iwhich == 1 .or. iwhich == 0) then
        call vpois3d(1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)
        call pipest3d(pipeshpe,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                scrtch,xywork,zwork,
     &                piperadi,pipen,pipenz2,pipe8th,pipex,pipey,cap3d,kpvt,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)
        if (iwhich == 0) then
c         --- Copy rho back into phi since phi was scrambled during the
c         --- calculation of the matrix.
          phi(0:nx,0:ny,0:nzlocal) = rho
        endif
      endif

      if (iwhich == 1) return

c do first partial field solve (without z transform at end)
      call vpois3d( 2,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 3,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d(13,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

c apply cap matrix for each kz to get induced rho(kz)
c     --- extract phi error from FFT solution and zero pipeq
      do i=1,pipen
        ix = pipex(i)
        xx = pipex(i) - ix
        iy = pipey(i)
        yy = pipey(i) - iy
        do iz=0,pipenz-1
          pipephi(i,iz) = - (phi(ix  ,iy  ,iz)*(1.-xx)*(1.-yy) +
     &                       phi(ix+1,iy  ,iz)*    xx *(1.-yy) +
     &                       phi(ix  ,iy+1,iz)*(1.-xx)*    yy  +
     &                       phi(ix+1,iy+1,iz)*    xx *    yy)
          pipeq(i,iz) = 0.
        enddo
      enddo

c     --- mpy pipephit by cap matrix (using reduced storage)
      do iz=0,pipenz-1
        iiz = iz
        if (iiz > pipenz2) iiz = pipenz - iz
        do j=0,pipe8th-2
          do i=1,pipen
            pipeq(i,iz) = pipeq(i,iz) + pipephi(j+1,iz)*cap3d(i,j+1,iiz)
          enddo
          do i=1,2*pipe8th-1
            pipeq(i,iz) = pipeq(i,iz) +
     &               pipephi(j+1*pipe8th,iz)*cap3d(2*pipe8th-i,pipe8th-j,iiz)
          enddo
          do i=2*pipe8th,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &           pipephi(j+1*pipe8th,iz)*cap3d(pipen-i+2*pipe8th,pipe8th-j,iiz)
          enddo
          do i=1,2*pipe8th-2
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+2*pipe8th-1,iz)*cap3d(pipen-2*pipe8th+2+i,j+1,iiz)
          enddo
          do i=2*pipe8th-1,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &                   pipephi(j+2*pipe8th-1,iz)*cap3d(i-2*pipe8th+2,j+1,iiz) 
          enddo
          do i=1,4*pipe8th-3
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+3*pipe8th-2,iz)*cap3d(4*pipe8th-2-i,pipe8th-j,iiz)
          enddo
          do i=4*pipe8th-2,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &       pipephi(j+3*pipe8th-2,iz)*cap3d(pipen+4*pipe8th-2-i,pipe8th-j,iiz)
          enddo
          do i=1,4*pipe8th-4
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+4*pipe8th-3,iz)*cap3d(pipen-4*pipe8th+4+i,j+1,iiz)
          enddo
          do i=4*pipe8th-3,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &                   pipephi(j+4*pipe8th-3,iz)*cap3d(i-4*pipe8th+4,j+1,iiz) 
          enddo
          do i=1,6*pipe8th-5
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+5*pipe8th-4,iz)*cap3d(6*pipe8th-4-i,pipe8th-j,iiz)
          enddo
          do i=6*pipe8th-4,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &       pipephi(j+5*pipe8th-4,iz)*cap3d(pipen+6*pipe8th-4-i,pipe8th-j,iiz)
          enddo
          do i=1,6*pipe8th-6
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+6*pipe8th-5,iz)*cap3d(pipen-6*pipe8th+6+i,j+1,iiz)
          enddo
          do i=6*pipe8th-5,pipen
            pipeq(i,iz) = pipeq(i,iz) +
     &                   pipephi(j+6*pipe8th-5,iz)*cap3d(i-6*pipe8th+6,j+1,iiz) 
          enddo
          pipeq(1,iz) = pipeq(1,iz) +
     &                 pipephi(j+7*pipe8th-6,iz)*cap3d(1,pipe8th-j,iiz)
          do i=2,8*pipe8th-8
            pipeq(i,iz) = pipeq(i,iz) +
     &             pipephi(j+7*pipe8th-6,iz)*cap3d(8*pipe8th-6-i,pipe8th-j,iiz)
          enddo
        enddo
      enddo

c     --- set up for field solve
      phi(0:nx,0:ny,0:nzlocal) = rho

c     --- transform in z
      call vpois3d(10,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

c     --- deposit induced rho(kz) onto grid
      do i=1,pipen
        ix = pipex(i)
        xx = pipex(i) - ix
        iy = pipey(i)
        yy = pipey(i) - iy
        do iz=0,pipenz-1
          phi(ix  ,iy  ,iz) = phi(ix  ,iy  ,iz) + pipeq(i,iz)*(1.-xx)*(1.-yy)
          phi(ix+1,iy  ,iz) = phi(ix+1,iy  ,iz) + pipeq(i,iz)*    xx *(1.-yy)
          phi(ix  ,iy+1,iz) = phi(ix  ,iy+1,iz) + pipeq(i,iz)*(1.-xx)*    yy 
          phi(ix+1,iy+1,iz) = phi(ix+1,iy+1,iz) + pipeq(i,iz)*    xx *    yy 
        enddo
      enddo

c     --- complete the field solve
      call vpois3d(12,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 4,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)
      call vpois3d( 5,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nzlocal,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)

      return
      end
c============================================================================== 
c============================================================================== 
c============================================================================== 
c============================================================================== 
c
c Capacity matrix solver for 3-d structures
c
c============================================================================== 
      subroutine setcap3d(ncndpts,quadx,quady,quadz,quadcap,kpvt,nx,ny,nzlocal,
     &                    nz,
     &                    phi,work,xywork,zwork,kxsq,kysq,kzsq,attx,atty,attz,
     &                    xlen,ylen,zlen,filt,l2symtry,l4symtry,
     &                    bound0,boundnz,boundxy)
      integer(ISZ):: ncndpts,nx,ny,nzlocal,nz
      real(kind=8):: quadx(ncndpts,4), quady(ncndpts,4), quadz(ncndpts)
      real(kind=8):: quadcap(4*ncndpts,4*ncndpts)
      integer(ISZ):: kpvt(4*ncndpts)
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: work(nx+ny-4),xywork(*),zwork(2,0:nx,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3)
      real(kind=8):: xlen, ylen, zlen
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c  Calculates the capacity matrix for a single lens, which consists of a block
c  of four conductors.
c  input:
c  ncndpts   number of points on one conductor, total number of points is
c            four times ncndpts
c  quadx, quady, quadz   list of conductor points, second dimension is number
c            of conductor, ordered like    1
c                                         3 4
c                                          2
c            quadz is the same for each of the four conductors
c  output:
c  quadcap   on return is the inverse matrix in reduced form (upper triangular)
c  kpvt      on return is the pivot vector created during reduction
c
c  rest of variables are input and have the standard meanings.

      real(kind=8),allocatable:: zerofilt(:,:)
      integer(ISZ):: allocerror
      real(kind=8):: xx,yy,zz,dx,dy,dz
      integer(ISZ):: ix,iy,iz,i,j,ic,info

      dx = xlen/nx
      dy = ylen/ny
      dz = zlen/nz

c  calculate inverse cap matrix

c  first calculate capacities for one conductor
      do i=1,ncndpts
c       call zeroarry(phi,(nx+1)*(ny+1)*(nzlocal+1))
        phi = 0.
c  deposit unit charge on grid
        ix = quadx(i,1)/dx + nx/2
        xx = quadx(i,1)/dx + nx/2 - ix
        iy = quady(i,1)/dy + ny/2
        yy = quady(i,1)/dy + ny/2 - iy
        iz = quadz(i)/dz + nzlocal/2 - 4
        zz = quadz(i)/dz + nzlocal/2 - 4 - iz
        phi(ix  ,iy  ,iz  ) = (1. - xx)*(1. - yy)*(1. - zz)
        phi(ix+1,iy  ,iz  ) = (     xx)*(1. - yy)*(1. - zz)
        phi(ix  ,iy+1,iz  ) = (1. - xx)*(     yy)*(1. - zz)
        phi(ix+1,iy+1,iz  ) = (     xx)*(     yy)*(1. - zz)
        phi(ix  ,iy  ,iz+1) = (1. - xx)*(1. - yy)*(     zz)
        phi(ix+1,iy  ,iz+1) = (     xx)*(1. - yy)*(     zz)
        phi(ix  ,iy+1,iz+1) = (1. - xx)*(     yy)*(     zz)
        phi(ix+1,iy+1,iz+1) = (     xx)*(     yy)*(     zz)
c  do 3-d field solve, turning filtering off
        allocate(zerofilt(size(filt,1),size(filt,2)),stat=allocerror)
        if (allocerror /= 0) then
          print*,"setcap3d: allocation error ",allocerror,
     &           ": could not allocate zerofilt to shape ",
     &           size(filt,1),size(filt,2)
          call kaboom("setcap3d: allocation error")
          return
        endif
        zerofilt = 0.
        call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &               zerofilt,xlen,ylen,zlen,nx,ny,nzlocal,nz,work,
     &               xywork,zwork,0,l2symtry,l4symtry,bound0,boundnz,boundxy)
        deallocate(zerofilt)
c  fill i'th column of inverse cap matrix
        do ic = 1,4
          do j=1,ncndpts
            ix = quadx(j,ic)/dx + nx/2
            xx = quadx(j,ic)/dx + nx/2 - ix
            iy = quady(j,ic)/dy + ny/2
            yy = quady(j,ic)/dy + ny/2 - iy
            iz = quadz(j)/dz + nzlocal/2 - 4
            zz = quadz(j)/dz + nzlocal/2 - 4 - iz
            quadcap(i,j+(ic-1)*ncndpts) =
     &                     phi(ix  ,iy  ,iz  )*(1. - xx)*(1. - yy)*(1. - zz) +
     &                     phi(ix+1,iy  ,iz  )*(     xx)*(1. - yy)*(1. - zz) +
     &                     phi(ix  ,iy+1,iz  )*(1. - xx)*(     yy)*(1. - zz) +
     &                     phi(ix+1,iy+1,iz  )*(     xx)*(     yy)*(1. - zz) +
     &                     phi(ix  ,iy  ,iz+1)*(1. - xx)*(1. - yy)*(     zz) +
     &                     phi(ix+1,iy  ,iz+1)*(     xx)*(1. - yy)*(     zz) +
     &                     phi(ix  ,iy+1,iz+1)*(1. - xx)*(     yy)*(     zz) +
     &                     phi(ix+1,iy+1,iz+1)*(     xx)*(     yy)*(     zz)
          enddo
        enddo
      enddo
          
c  then fill in rest of matrix using symmetries of the conductor arrangement
        do i=1,ncndpts
          do j=1,ncndpts
            quadcap(i+  ncndpts,j          ) = quadcap(i,j+  ncndpts)
            quadcap(i+  ncndpts,j+  ncndpts) = quadcap(i,j          )
            quadcap(i+  ncndpts,j+2*ncndpts) = quadcap(i,j+3*ncndpts)
            quadcap(i+  ncndpts,j+3*ncndpts) = quadcap(i,j+2*ncndpts)
            quadcap(i+2*ncndpts,j          ) = quadcap(i,j+2*ncndpts)
            quadcap(i+2*ncndpts,j+  ncndpts) = quadcap(i,j+3*ncndpts)
            quadcap(i+2*ncndpts,j+2*ncndpts) = quadcap(i,j          )
            quadcap(i+2*ncndpts,j+3*ncndpts) = quadcap(i,j+  ncndpts)
            quadcap(i+3*ncndpts,j          ) = quadcap(i,j+3*ncndpts)
            quadcap(i+3*ncndpts,j+  ncndpts) = quadcap(i,j+2*ncndpts)
            quadcap(i+3*ncndpts,j+2*ncndpts) = quadcap(i,j+  ncndpts)
            quadcap(i+3*ncndpts,j+3*ncndpts) = quadcap(i,j          )
          enddo
        enddo

c  now reduce to allow inversion
      call dsifa(quadcap,4*ncndpts,4*ncndpts,kpvt,info)

      return
      end
c==========================================================================
      subroutine capfs(iwhich,nx,ny,nzlocal,nz,phi,rho,
     &                 xlen,ylen,zlen,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                 work,workz,xywork,zwork,xmmax,pipeshpe,
     &                 l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Capvars
      use GlobalVars
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal), rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: work(nx+ny-4), workz(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: xlen, ylen, zlen, filt(5,3), xmmax
      character(*):: pipeshpe
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c  Capacity matrix field solve routine
c  input:
c  pipeshpe   shape of conductors, default is circlular, "hyp" is hyperbolic
c  rest of variables have the standard meanings

      real(kind=8):: xx,yy,zz,dx,dy,dz,dxi,dyi,dzi
      integer(ISZ):: iend,iquad,ic,iq,ix,iy,iz

      dx = xlen/nx
      dy = ylen/ny
      dz = zlen/nz
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c  initialization
      if (iwhich == 0 .or. iwhich == 1) then
c    intialize vpois3d vars
        call vpois3d (1, phi, phi, kxsq, kysq, kzsq,
     &                attx, atty, attz, filt, xlen, ylen, zlen,
     &                nx, ny, nzlocal,nz, work, xywork, zwork, 0,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)
c    set up conductor point arrays and gchange array sizes
        call cndctr3d(nx,ny,xmmax,dx,dy,quadradi,quadcent,
     &                ncndpts, quadx, quady, quadz, quadv, nzpts, nendquad,
     &                nzquad, qdslclen, pipeshpe, loadquad)
        call gchange("Capvars",0)
c    set up capacity matrix
        call setcap3d(ncndpts,quadx,quady,quadz,quadcap,kkkk,nx,ny,nzlocal,nz,
     &                phi,work,xywork,zwork,kxsq,kysq,kzsq,
     &                attx,atty,attz,xlen,ylen,
     &                zlen,filt,l2symtry,l4symtry,bound0,boundnz,boundxy)
        if (iwhich == 0) then
c         --- Copy rho back into phi since phi was scrambled during the
c         --- calculation of the matrix.
          phi(0:nx,0:ny,0:nzlocal) = rho
        endif
      endif

c  return if only being initialized
      if (iwhich == 1) return

c  do 3-d field solve
      call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,work,xywork,zwork,
     &             0,l2symtry,l4symtry,bound0,boundnz,boundxy)

c  zero induced charge array
c     call zeroarry(quadrho,numends*4*ncndpts)
      quadrho = 0.

c loop over the quadrupole focusing lenses
      do iend=1,numends

c    find z's of conducting points relative to grid
        do iquad=1,ncndpts
          workz(iquad) = (quadz(iquad) + quadend(iend))*dzi
          if (workz(iquad) < 0) workz(iquad) = 0.
          if (workz(iquad) >= nzlocal) workz(iquad) = nzlocal - 1.
        enddo

c    apply 3-d cap matrix to the four conductors
        do ic = 1,4

c    get (volt - vshift - phi) on conductor
          do iq=1,ncndpts
            ix = quadx(iq,ic)*dxi + nx/2
            xx = quadx(iq,ic)*dxi + nx/2 - ix
            iy = quady(iq,ic)*dyi + ny/2
            yy = quady(iq,ic)*dyi + ny/2 - iy
            iz = workz(iq)
            zz = workz(iq) - iz
            quadrho(iq,ic,iend) = quadvlt(iend)*quadv(iq,ic) + vshift(iend) -
     &                    phi(ix  ,iy  ,iz  )*(1. - xx)*(1. - yy)*(1. - zz) -
     &                    phi(ix+1,iy  ,iz  )*(     xx)*(1. - yy)*(1. - zz) -
     &                    phi(ix  ,iy+1,iz  )*(1. - xx)*(     yy)*(1. - zz) -
     &                    phi(ix+1,iy+1,iz  )*(     xx)*(     yy)*(1. - zz) -
     &                    phi(ix  ,iy  ,iz+1)*(1. - xx)*(1. - yy)*(     zz) -
     &                    phi(ix+1,iy  ,iz+1)*(     xx)*(1. - yy)*(     zz) -
     &                    phi(ix  ,iy+1,iz+1)*(1. - xx)*(     yy)*(     zz) -
     &                    phi(ix+1,iy+1,iz+1)*(     xx)*(     yy)*(     zz)
          enddo
        enddo

c     solve for induced charge by solving phi = matrix*rho
        call dsisl(quadcap,4*ncndpts,4*ncndpts,kkkk,quadrho(1,1,iend))

      enddo

c  reload rho back into phi
      phi(0:nx,0:ny,0:nzlocal) = rho

c  if periodic in z, zero out iz=nzlocal plane so that plane will only have
c  image charges (plane 0 will have both image charges and beam charge)
c     if (bound0==periodic) then
        do ix=0,nx
          do iy=0,ny
            phi(ix,iy,nzlocal) = 0.
          enddo
        enddo
c     endif

c  load induced charges onto phi
      do iend=1,numends

c    find z's relative to grid and zero quadrho that are off the grid
        do iquad=1,ncndpts
          workz(iquad) = (quadz(iquad) + quadend(iend))*dzi
          if (workz(iquad) < 0) then
            workz(iquad) = 0.
            quadrho(iquad,1,iend) = 0.
            quadrho(iquad,2,iend) = 0.
            quadrho(iquad,3,iend) = 0.
            quadrho(iquad,4,iend) = 0.
          endif
          if (workz(iquad) >= nzlocal) then
            workz(iquad) = nzlocal - 1.
            quadrho(iquad,1,iend) = 0.
            quadrho(iquad,2,iend) = 0.
            quadrho(iquad,3,iend) = 0.
            quadrho(iquad,4,iend) = 0.
          endif
        enddo

c    load charge onto each of the four conductors in this lens
        do ic = 1,4
          do iq=1,ncndpts
            ix = quadx(iq,ic)*dxi + nx/2
            xx = quadx(iq,ic)*dxi + nx/2 - ix
            iy = quady(iq,ic)*dyi + ny/2
            yy = quady(iq,ic)*dyi + ny/2 - iy
            iz = workz(iq)
            zz = workz(iq) - iz
            phi(ix  ,iy  ,iz  ) = phi(ix  ,iy  ,iz  ) +
     &                            (1.-xx)*(1.-yy)*(1.-zz)*quadrho(iq,ic,iend)
            phi(ix+1,iy  ,iz  ) = phi(ix+1,iy  ,iz  ) +
     &                            (   xx)*(1.-yy)*(1.-zz)*quadrho(iq,ic,iend)
            phi(ix  ,iy+1,iz  ) = phi(ix  ,iy+1,iz  ) +
     &                            (1.-xx)*(   yy)*(1.-zz)*quadrho(iq,ic,iend)
            phi(ix+1,iy+1,iz  ) = phi(ix+1,iy+1,iz  ) +
     &                            (   xx)*(   yy)*(1.-zz)*quadrho(iq,ic,iend)
            phi(ix  ,iy  ,iz+1) = phi(ix  ,iy  ,iz+1) +
     &                            (1.-xx)*(1.-yy)*(   zz)*quadrho(iq,ic,iend)
            phi(ix+1,iy  ,iz+1) = phi(ix+1,iy  ,iz+1) +
     &                            (   xx)*(1.-yy)*(   zz)*quadrho(iq,ic,iend)
            phi(ix  ,iy+1,iz+1) = phi(ix  ,iy+1,iz+1) +
     &                            (1.-xx)*(   yy)*(   zz)*quadrho(iq,ic,iend)
            phi(ix+1,iy+1,iz+1) = phi(ix+1,iy+1,iz+1) +
     &                            (   xx)*(   yy)*(   zz)*quadrho(iq,ic,iend)
          enddo
        enddo
      enddo

c Make phi periodic by adding planes iz=0 and nzlocal and placing result back into
c both planes.  (Note that plane 0 has beam charge and additional image charge 
c while plane nzlocal has only image charge.)
c     if (bound0==periodic) then
        do ix=0,nx
          do iy=0,ny
            phi(ix,iy,0) = phi(ix,iy,0) + phi(ix,iy,nzlocal)
            phi(ix,iy,nzlocal) = phi(ix,iy,0)
          enddo
        enddo
c     endif

c  finished with conductors, so...
c  now do full 3-d field solve, with induced charges
      call vpois3d(-1,phi,phi,kxsq,kysq,kzsq,
     &              attx,atty,attz,filt,xlen,ylen,zlen,
     &              nx,ny,nzlocal,nz,work,xywork,zwork,
     &              0,l2symtry,l4symtry,bound0,boundnz,boundxy)

      return
      end
c==========================================================================
      subroutine cndctr3d(nx,ny,xmmax,dx,dy,quadradi,quadcent, ncndpts,
     &                    quadx, quady, quadz, quadv, nzpts, nendquad,nzquad,
     &                    qdslclen, pipeshpe, loadquad)
      integer(ISZ):: nx,ny,ncndpts,nzpts,nendquad,nzquad
      real(kind=8):: quadx(ncndpts,4), quady(ncndpts,4), quadz(ncndpts)
      real(kind=8):: quadv(ncndpts,4)
      real(kind=8):: xmmax,dx,dy,quadradi,quadcent,qdslclen
      character(*):: pipeshpe
      logical(ISZ):: loadquad

c  Calculate set of points representing the conductor
c  Only calculate one quarter, rest are found by symmetry.
c  (Same with user inputed points if loadquad is false)

      real(kind=8):: theta0,theta,rad
      real(kind=8):: piperadi,yymax,qcmge
      integer(ISZ):: i,j,iend,iquad,ixmin,ixmax,iymin,iymax,iystep,iq,ncndptst,ix,iy

      if (loadquad) then

c       --- find points on surface of first slice of one conductor
        if (pipeshpe == "hyp") then
c       --- hyperbolic conductors
          piperadi = (quadcent - quadradi)
c         --- find max y where quad interceps first row inside grid
          yymax = (xmmax*0.5-dx)**2 - piperadi**2
          if (yymax > 0.) then
            yymax = sqrt(yymax)
          else
            nzpts = 0
            call remark("CNDCTR3D WARNING: quads not within grid, number of points is zero")
          endif
          do i=1,nzpts/2
            quady(i,1) = yymax*(i*1./(nzpts*0.5))**2
            quadx(i,1) = sqrt(piperadi**2 + (quady(i,1) - xmmax*0.5)**2)
            quadz(i) = 0.
          enddo
          do i=nzpts/2+1,nzpts
            quady(i,1) = - yymax*((i-nzpts*0.5)*1./(nzpts*0.5))**2
            quadx(i,1) = sqrt(piperadi**2 + (quady(i,1) - xmmax*0.5)**2)
            quadz(i) = 0.
          enddo
        else
c       --- round conductors are the default
c          --- find angle where quad interceps first row inside grid
          qcmge =  (quadcent - xmmax + dx)/quadradi
          if (abs(qcmge) <= 1.) then
            theta0 = acos(qcmge)
          else if (qcmge > 1.) then
c           --- quads outside grid
            nzpts = 0
            call remark("CNDCTR3D WARNING: quads not within grid, number of points is zero")
          else if (qcmge < -1) then
c           --- quads completely inside grid
c           --- set to Pi
            theta0 = 4. * atan(1.)
          endif
          do i=1,nzpts
            theta = theta0 - 2.*theta0*(i-1.)/(nzpts-1.)
            quadx(i,1) = ( - quadradi*cos(theta) + quadcent)
            quady(i,1) = quadradi*sin(theta)
            quadz(i) = 0.
          enddo
        endif

c       --- fill in rest of z surface slices in this conductor
          do i=2,nzquad
            do j=1,nzpts
            quadx(j+(i-1)*nzpts,1) = quadx(j,1)
            quady(j+(i-1)*nzpts,1) = quady(j,1)
            quadz(j+(i-1)*nzpts) = (i-1.)*qdslclen
          enddo
        enddo

c       --- find points on ends of quad
c       --- takes points on mesh within bounding circle
c       --- not quite right for hyperbolic quads
c       --- guaranteed to have points centered about ny/2
        iend = 0
        iquad = nzquad*nzpts
c       --- find rectangular bounds of semi-circle
        ixmin = (quadcent - quadradi)/dx + 1 + nx/2
        ixmax = nx - 1
        iymin = ny/2
        iymax = min(quady(1,1)/dy + ny/2 - 1,ny-1.)
c       --- distance between points in y, increases with x to give nice distrbtn
        iystep = 1
        do ix = ixmin,ixmax
          iystep = iystep + 1
          do iy = iymin,iymax,iystep
c           --- radius of current points from quad center
            rad = ((ix-nx/2)*dx - quadcent)**2 + ((iy-ny/2)*dy)**2
c           --- if within quad, use it, use twice if off axis to give symmetric
c           --- distribution
            if (rad < (quadradi-1.1*dx)**2) then
              if (iend < nendquad) then
                iend = iend + 1
                quadx(iquad + iend,1) = (ix - nx/2)*dx
                quady(iquad + iend,1) = (iy - ny/2)*dy
                quadz(iquad + iend) = 0.
              endif
              if (iend < nendquad .and. iy .ne. ny/2) then
                iend = iend + 1
                quadx(iquad + iend,1) = (ix - nx/2)*dx
                quady(iquad + iend,1) = (ny - iy - ny/2)*dy
                quadz(iquad + iend) = 0.
              endif
            endif
          enddo
        enddo
  
c       --- reset nendquad and ncndpts
        ncndptst = ncndpts - 2*nendquad + 2*iend
        nendquad = iend

c       --- fill in other end of conductor
        iquad = nzquad*nzpts + 1
        do i=iquad,iquad+nendquad-1
          quadx(i+nendquad,1) = quadx(i,1)
          quady(i+nendquad,1) = quady(i,1)
          quadz(i+nendquad) = (nzquad - 1.)*qdslclen
        enddo

c       --- set sign of voltage
        do iq=1,ncndptst
          quadv(iq,1) = 1.
        enddo

c     --- end of point loader
      else
        ncndptst = ncndpts
      endif

c  fill in other three conductors
      do iq=1,ncndptst
        quadx(iq,2) = - quadx(iq,1)
        quady(iq,2) = - quady(iq,1)
        quadv(iq,2) =   quadv(iq,1)

        quadx(iq,3) =   quady(iq,1)
        quady(iq,3) =   quadx(iq,1)
        quadv(iq,3) = - quadv(iq,1)

        quadx(iq,4) = - quady(iq,1)
        quady(iq,4) = - quadx(iq,1)
        quadv(iq,4) = - quadv(iq,1)
      enddo

c  final value of ncndpts
      ncndpts = ncndptst

      return
      end
c=============================================================================
      subroutine findqdnd(nquad,quadzs,quadze,quadde,quadvx,quadvy,
     &                    zmmin,zgrid,nzlocal,dz,
     &                    numends,nendmax,quadend,quadvlt,
     &                    vshift,quadcent,quadradi)
      integer(ISZ):: nquad,nzlocal,numends,nendmax
      real(kind=8):: quadzs(0:nquad), quadze(0:nquad), quadde(0:nquad)
      real(kind=8):: quadvx(0:nquad), quadvy(0:nquad)
      real(kind=8):: quadvlt(nendmax), quadend(nendmax), vshift(nendmax)
      real(kind=8):: zmmin,zgrid,dz,quadcent,quadradi

c  Find starts of quads in grid, and put quad voltage into quadvlt
c  Also sets array vshift
c  dependent upon quad variables

      real(kind=8):: offset
      integer(ISZ):: iz,iq,iqprev

      numends = 0

c     --- If there are no quads, then just return
      if (nquad < 0) return

c  get rest of quads within mesh
      iq = 0
      iqprev = -1
      do iz=1,nzlocal-1
        call getquadid(zmmin+iz*dz+zgrid,offset,iq,1)
        if (iq .ne. iqprev .and. numends < nendmax) then
          numends = numends + 1
          quadend(numends) = (quadzs(iq) + offset - zgrid - zmmin)
          quadvlt(numends) = (quadvx(iq) - quadvy(iq))*0.5
          if (quadvx(iq) .ne. 0. .or. quadvy(iq) .ne. 0.)
     &      vshift(numends) = (quadvx(iq) + quadvy(iq))*0.5
        endif
        iqprev = iq
      enddo

      return
      end
c=============================================================================
      subroutine vcap3d(iwhich,rho,phi,kxsq,kysq,kzsq, 
     &                  attx,atty,attz,filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                  scrtch,xywork,zwork,xmmax,zmmin,zgrid,pipeshpe,
     &                  l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Capvars
      use Lattice
      integer(ISZ):: iwhich,nx,ny,nzlocal,nz
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*), xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen
      real(kind=8):: xmmax,zmmin,zgrid
      character(*):: pipeshpe
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy
 
c  External interface to capacity matrix field solver.
c  Passes external and local variables to field solver routines.
c  If iwhich is '0' or '1', also allocates local variables


      real(kind=8):: dz,offset
      integer(ISZ):: iq

      dz = zlen/nz

c  find quad starts
      call findqdnd(nquad,quadzs,quadze,quadde,quadvx,quadvy,
     &              zmmin,zgrid,nzlocal,dz,
     &              numends,nendmax,quadend,quadvlt,vshift,quadcent,quadradi)

c  Allocate local variables
      if (iwhich == 0 .or. iwhich == 1) then
        if (nquad >= 0 .and. loadquad) then
          if (nzpts == 0) nzpts = 16
          if (nendquad == 0) nendquad = 20
          iq = 0
          call getquadid(zmmin+zgrid,offset,iq,1)
          nzquad =  ((quadze(iq) - quadzs(iq))/dz + 1.5)
          qdslclen = (quadze(iq) - quadzs(iq))/(nzquad-1.)
          ncndpts = nzquad*nzpts + 2*nendquad
          nendmax = numends + 2
        endif
        call gchange("Capvars",0)
      endif

c  Call field solver, which may also do initialization
      call capfs(iwhich,nx,ny,nzlocal,nz,phi,rho,xlen,ylen,zlen,kxsq,kysq,kzsq,
     &           attx,atty,attz,filt,scrtch,scrtch,xywork,zwork,
     &           xmmax,pipeshpe,l2symtry,l4symtry,bound0,boundnz,boundxy)

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================

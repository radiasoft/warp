#include "top.h"
c=======================================================================
c@(#) File F3D.M, version $Revision: 3.38 $, $Date: 2002/04/30 13:56:27 $
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
c     5) A BASIS test driver (invoke via PACKAGE F3D;GENERATE)
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
c  (The only BASIS/MPPL dependent coding is use of STDOUT, return (value),
c  and the next few lines).

c  Bring in the common block; allocate dynamic arrays ("group allot")


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
     &              lx,ly,lz,nx,ny,nz,work,work2d,zwork,0,.false.,.false.)
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
     &                 lx,ly,lz,nx,ny,nz,work,work2d,zwork,0,.false.,.false.)
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
     &                 lx,ly,lz,nx,ny,nz,work,work2d,zwork,0,.false.,.false.)
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
     &              lx,ly,lz,nx,ny,nz,work,work2d,zwork,0,.false.,.false.)
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
     &              lx,ly,lz,nx,ny,nz,work,work2d,zwork,0,.false.,.false.)
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

c  BASIS interface to VPOIS3D, using variables from the database for
c  package F3D.  The user's application package should contain a
c  similar subroutine for convenience, using its own database variables.
c  To avoid confusion, it should have a different name, say "vp3"


      call vpois3d (iwhich,b,b,kxsq,kysq,kzsq,attx,atty,attz,
     &              filt,lx,ly,lz,nx,ny,nz,work,work2d,zwork,ibc,
     &              .false.,.false.)

      return
      end
c=======================================================================
c=============================================================================
c=============================================================================
      subroutine vpois3d(iwhich,a,ak,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                   lx,ly,lz,nx,ny,nz,w,xywork,zwork,ibc,
     &                   l2symtry,l4symtry)
      use Constant
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: iwhich,nx,ny,nz,ibc
      real(kind=8)::  a(0:nx, 0:ny, 0:nz)
      real(kind=8):: ak(0:nx, 0:ny, -nz/2:nz/2)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz), filt(5,3)
      real(kind=8):: lx, ly, lz, w(*), xywork(0:nx,0:ny), zwork(2,0:nx,0:nz)
      logical(ISZ):: l2symtry,l4symtry
      
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
c                IBC ... Boundary condition switch, char*8 (for future use)
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
      integer(ISZ):: ikx,iky,ikz,ix,iy,iz
      real(kind=8):: klast, kdxb2, kdyb2, kdzb2, kdx, kdy, kdz
      real(kind=8):: dx,dy,dz,pib2,norm
      integer(ISZ):: fullnz
      real(kind=8):: fulllz
      integer(ISZ):: i,j
      logical(ISZ):: lattenuate
#ifdef MPIPARALLEL
      integer(ISZ):: nx_tran,ny_tran,ikxm,ikym
      logical(ISZ):: transposed
#endif

      fullnz = nz
      fulllz = lz

#ifdef MPIPARALLEL
c  Value of nz and lz for full array across all processors.
      fullnz = nz*nslaves
      fulllz = lz*nslaves

c  The input data is always assumed to be not transposed
      transposed = .false.
#endif

#ifdef ESSL

c  The ESSL routines themselves have no limit on the number of points.
c  Because of the way vsftz works, ny must be even (this is easily fixable).
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
     & .or.(( fullnz .ne. 2**int( .5 + log(1.*fullnz) / log(2.) ) ) .and.
     &      ( iwhich .ne. 1 .and. iwhich .ne. 12 .and. iwhich .ne. 13 .and.
     &        iwhich .ne. 14 ) ) ) then
         call remark("NX and NY must be powers of 2, NZ also unless using tridiag solver")
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

      nz2 = fullnz / 2
      dx = lx / nx
      dy = ly / ny
      dz = lz / nz
      pib2 = .5 * pi

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
            kzsq(nz2+ikz) = (2./dz*sin(pi*ikz/fullnz))**2 * eps0
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
           kdz = 2. * pi * ikz / fullnz
           attz(nz2+ikz) = (1.+2.*filt(5,3)*cos(kdz)) / (1.+2.*filt(5,3))
         enddo
         if (filt(1,3) .ne. 0. .or. filt(2,3) .ne. 0.) then
c          --- compute first form of spatial filtering
           do ikz = 0, nz2
              kdzb2 = pi * ikz / fullnz
              attz(nz2+ikz) = attz(nz2+ikz) * (exp( filt(1,3)*(sin(kdzb2))**2
     &         - filt(2,3)*(tan(kdzb2))**4 ))**2
           enddo
         endif
         if (filt(3,3) .ne. 0.) then
c          --- compute second form of spatial filtering
           klast = filt(3,3) * fullnz / (2. * pi)
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
          call cosqx (a(0,0,0), xywork, w(icox), nx, ny, 1,-1)
          call cosqy (a(0,0,0), xywork, w(icoy), nx, ny, 1,-1)
        elseif (l2symtry) then
          call vsftx (a(0,0,0), xywork, w(icpx), w(icmx), nx, ny, 1)
          call cosqy (a(0,0,0), xywork, w(icoy), nx, ny, 1,-1)
        else
          call vsftx (a(0,0,0), xywork, w(icpx), w(icmx), nx, ny, 1)
          call vsfty (a(0,0,0), xywork, w(icpy), w(icmy), nx, ny, 1)
        endif

      endif

c ----------------------------------------------------------------------------
c  Do the forward transform
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 2) then

c  take z-sequence of sine-sine transforms in the x-y plane.
c  We should multitask this loop, being careful about scratch
c  which is (nx+1)*(ny+1) per task.  For now, we use the plane
c  at nz for the single task's scratch, and so need no other.

      if (l4symtry) then
c       factor of .5dx for cosqx and .5*dy for cosqy
        norm = .125*dx*dy*dz
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call cosqx(a(0,0,iz),xywork,w(icox),nx,ny,0,-1)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
      else if (l2symtry) then
c       factor of .5*dy for cosqy
        norm = .25*dx*dy*dz
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(0,0,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
      else
        norm = .5 * dx * dy * dz
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(0,0,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call vsfty(a(0,0,iz),xywork,w(icpy),w(icmy),nx,ny,0)
        enddo
      endif

c  Here we perform the periodic transforms in z, two y's at a time.
c  We vectorize over x.  No scratch space is needed.

#ifdef MPIPARALLEL
c  For parallel code, transpose is used to redistribute phi so each axial
c  line is on a single processor.
        call transpose(nx,ny,nz,a,nx_tran,ny_tran)
        transposed = .true.
        if (ny_tran == 1) then
c         --- If ny_tran = 1, i.e. there is only one y line, pass nx as
c         --- zero so that transform done on 2 y lines actually
c         --- transforms on 2 x's next to each other.
          call vsftz(0,nx_tran,fullnz,a,0.,0,1,nx_tran+1,zwork)
        else
c         --- Otherwise pass normal arguments to vsftz
          call vsftz(nx_tran,ny_tran-1,fullnz,a,0.,1,0,(nx_tran+1)*ny_tran,
     &               zwork)
        endif
c       --- Transpose back if only doing forward transform.
        if (iwhich == 2) then
          call transposei(nx,ny,nz,a,nx_tran,ny_tran)
          transposed = .false.
        endif
#else
c  For the serial code, vsftz is called with normal arguments.
        call vsftz(nx,ny,nz,a,0.,1,1,(nx+1)*(ny+1),zwork)
#endif

      endif

c ----------------------------------------------------------------------------
c  Do the forward transform in x and y only
c ----------------------------------------------------------------------------

      if (iwhich == 12) then

c  take z-sequence of sine-sine transforms in the x-y plane.
c  We should multitask this loop, being careful about scratch
c  which is (nx+1)*(ny+1) per task.  For now, we use the plane
c  at nz for the single task's scratch, and so need no other.

      if (l4symtry) then
c       factor of .5*dx for cosqx and .5*dy for cosqy
        norm = .25*dx*dy
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call cosqx(a(0,0,iz),xywork,w(icox),nx,ny,0,-1)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
      else if (l2symtry) then
c       factor of .5*dy for cosqy
        norm = .5*dx*dy
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(0,0,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,-1)
        enddo
      else
        norm = dx * dy
        do iz = 0, nz-1
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
c          --- Invoke the transform routines for the current z-slice
           call vsftx(a(0,0,iz),xywork,w(icpx),w(icmx),nx,ny,0)
           call vsfty(a(0,0,iz),xywork,w(icpy),w(icmy),nx,ny,0)
        enddo
      endif

      endif

c ----------------------------------------------------------------------------
c  Do the forward transform in z only
c ----------------------------------------------------------------------------

      if (iwhich == 10) then

        norm = .5 * dz

#ifdef MPIPARALLEL
c  For parallel code, transpose is used to redistribute phi so each axial
c  line is on a single processor.
        call transpose(nx,ny,nz,a,nx_tran,ny_tran)

        if (ny_tran == 1) then
c         --- If there is only one y line, pass nx as zero so that
c         --- transform done on 2 y lines actually transforms on 2 x's
c         --- next to each other.
          call vsftz(0,nx_tran,fullnz,a,norm,0,1,nx_tran+1,zwork)
        else
c         --- pass normal arguments
          call vsftz(nx_tran,ny_tran-1,fullnz,a,norm,1,0,(nx_tran+1)*ny_tran,
     &               zwork)
        endif

        call transposei(nx,ny,nz,a,nx_tran,ny_tran)
#else
c  Make call for serial code.
        call vsftz(nx,ny,nz,a,norm,1,1,(nx+1)*(ny+1),zwork)
#endif

      endif

c ----------------------------------------------------------------------------
c  Apply k-space filter; attenuate mode-by-mode
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 3) then

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
            iy = int(my_index*ny/real(nslaves))
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
            call attenuate(nx_tran,ny_tran-1,fullnz,a,attx(ix),atty(iy),attz,
     &                     ikxm,ikym,1,0)
          else
            iz = my_index*nz
            call attenuate(nx,ny,nz,a,attx,atty,attz(iz),ikxmin,ikymin,1,1)
          endif
#else
c  Normal call for serial version.
          call attenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,1,1)
#endif

        endif

      endif

c ----------------------------------------------------------------------------
c  Convert rhok to phik.
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 4) then

#ifdef MPIPARALLEL
c  For parallel code, data may be arranged differently then serial version.
        if (transposed) then
c         --- Find location of transverse data relative to original dimensions.
          ix = mod(my_index*nx_tran,nx)
          iy = int(my_index*ny/real(nslaves))
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
          call rhotophi(nx_tran,ny_tran-1,fullnz,a,kxsq(ix),kysq(iy),kzsq,
     &                  ikxm,ikym,1,0)
        else
          iz = my_index*nz
          call rhotophi(nx,ny,nz,a,kxsq,kysq,kzsq(iz),ikxmin,ikymin,1,1)
        endif
#else
c  Call for Serial version
        call rhotophi(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,1,1)
#endif

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 5) then

c  Inverse transform and shift in z.

#ifdef MPIPARALLEL
c  For parallel code, inverse transpose is used to redistribute phi into its
c  original format.
        if (.not. transposed) then
          call transpose(nx,ny,nz,a,nx_tran,ny_tran)
          transposed = .true.
        endif
        if (ny_tran == 1) then
c         --- If there is only one y line, pass nx as zero so that
c         --- transform done on 2 y lines actually transforms on 2 x's
c         --- next to each other.
          call vsftzi(0,nx_tran,fullnz,a,0.,0,1,nx_tran+1,zwork)
        else
c         --- pass normal arguments
          call vsftzi(nx_tran,ny_tran-1,fullnz,a,0.,1,0,(nx_tran+1)*ny_tran,
     &                zwork)
        endif

c       --- Rearrange phi again back to original arrangement.
        call transposei(nx,ny,nz,a,nx_tran,ny_tran)
        transposed = .false.

#else
c  Call for serial version
        call vsftzi(nx,ny,nz,a,0.,1,1,(nx+1)*(ny+1),zwork)
#endif

c  inverse transform and normalization in xy.  We should multitask this loop.

      if (l4symtry) then
c       factor of .25/lx for cosqx and .25/ly for cosqy
        norm = .0625/(lx*ly*fulllz)
        do iz = 0, nz-1
c          --- Invoke the inverse transform routines for the current z-slice
           call cosqx(a(0,0,iz),xywork,w(icox),nx,ny,0,1)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,1)
c          --- set boundary elements to zero
           do iy = 0, ny
             a(nx,iy,iz) = 0.
           enddo
           do ix = 0, nx-1
             a(ix,ny,iz) = 0.
           enddo
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
        enddo
      else if (l2symtry) then
c       factor of 2./lx for vsftx and .25/ly for cosqy
        norm = .5/(lx*ly*fulllz)
        do iz = 0, nz-1
c          --- Invoke the inverse transform routines for the current z-slice
           call vsftx ( a(0,0,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
           call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,1)
c          --- set boundary elements to zero
           do iy = 0, ny
             a(0 ,iy,iz) = 0.
             a(nx,iy,iz) = 0.
           enddo
           do ix = 1, nx-1
             a(ix,ny,iz) = 0.
           enddo
c          --- Normalize the array
           do iy = 0, ny-1
             do ix = 0, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
        enddo
      else
        norm = 4. / (lx * ly * fulllz)
        do iz = 0, nz-1
c          --- Invoke the transform routines
           call vsfty ( a(0,0,iz), xywork, w(icpy), w(icmy), nx, ny, 0)
           call vsftx ( a(0,0,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
c          --- Normalize the array, and set boundary elements to zero
           do ix = 0, nx-1
             a(ix,0 ,iz) = 0.
             a(ix,ny,iz) = 0.
           enddo
           do iy = 1, ny-1
             a(0 ,iy,iz) = 0.
             a(nx,iy,iz) = 0.
             do ix = 1, nx-1
               a(ix,iy,iz) = a(ix,iy,iz) * norm
             enddo
           enddo
         enddo
      endif

c  Enforce periodic bc's in z (now done in w3d)

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform in x and y only
c ----------------------------------------------------------------------------

      if (iwhich == 13) then

c  inverse transform and normalization in xy.  We should multitask this loop.

        if (l4symtry) then
c         factor of .25/lx for cosqx and .25/ly for cosqy
          norm = .0625/(lx*ly)
          do iz = 0, nz-1
c            --- Invoke the inverse transform routines for the current z-slice
             call cosqx(a(0,0,iz),xywork,w(icox),nx,ny,0,1)
             call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,1)
c            --- set boundary elements to zero
             do iy = 0, ny
               a(nx,iy,iz) = 0.
             enddo
             do ix = 0, nx-1
               a(ix,ny,iz) = 0.
             enddo
c            --- Normalize the array
             do iy = 0, ny-1
               do ix = 0, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
          enddo
        else if (l2symtry) then
c         factor of .25/ly for cosqy
          norm = .5/(lx*ly)
          do iz = 0, nz-1
c            --- Invoke the inverse transform routines for the current z-slice
             call vsftx ( a(0,0,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
             call cosqy(a(0,0,iz),xywork,w(icoy),nx,ny,0,1)
c            --- set boundary elements to zero
             do iy = 0, ny
               a(0 ,iy,iz) = 0.
               a(nx,iy,iz) = 0.
             enddo
             do ix = 1, nx-1
               a(ix,ny,iz) = 0.
             enddo
c            --- Normalize the array
             do iy = 0, ny-1
               do ix = 0, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
          enddo
        else
          norm = 4. / (lx * ly)
          do iz = 0, nz-1
c            --- Invoke the transform routines
             call vsfty ( a(0,0,iz), xywork, w(icpy), w(icmy), nx, ny, 0)
             call vsftx ( a(0,0,iz), xywork, w(icpx), w(icmx), nx, ny, 0)
c            --- Normalize the array, and set boundary elements to zero
             do ix = 0, nx-1
               a(ix,0 ,iz) = 0.
               a(ix,ny,iz) = 0.
             enddo
             do iy = 1, ny-1
               a(0 ,iy,iz) = 0.
               a(nx,iy,iz) = 0.
               do ix = 1, nx-1
                 a(ix,iy,iz) = a(ix,iy,iz) * norm
               enddo
             enddo
           enddo

        endif

c  Enforce periodic bc's in z (resets nz plane to clear out garbage left over
c  from that plane being used as scratch space) (now done in w3d)

      endif
c ----------------------------------------------------------------------------
c  Do the inverse transform only in z
c ----------------------------------------------------------------------------

      if (iwhich == 11) then

c  Inverse transform and shift in z.

        norm = 1./fulllz

#ifdef MPIPARALLEL
c  For parallel code, inverse transpose is used to redistribute phi into its
c  original format.
        call transpose(nx,ny,nz,a,nx_tran,ny_tran)

        if (ny_tran == 1) then
c         --- If there is only one y line, pass nx as zero so that
c         --- transform done on 2 y lines actually transforms on 2 x's
c         --- next to each other.
          call vsftzi(0,nx_tran,fullnz,a,norm,0,1,nx_tran+1,zwork)
        else
c         --- pass normal arguments
          call vsftzi(nx_tran,ny_tran-1,fullnz,a,norm,1,0,(nx_tran+1)*ny_tran,
     &                zwork)
        endif

c       --- Rearrange phi again back to original arrangement.
        call transposei(nx,ny,nz,a,nx_tran,ny_tran)

#else
c  Serial version
        call vsftzi(nx,ny,nz,a,norm,1,1,(nx+1)*(ny+1),zwork)
#endif

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
        iz = my_index*nz
        call unattenuate(nx,ny,nz,a,attx,atty,attz(iz),ikxmin,ikymin,1,1)
#else
c  Call for Serial version
        call unattenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,1,1)
#endif

      endif

c ----------------------------------------------------------------------------
c  Multiply by k^2 factor, in contrast to the usual division.
c  This is useful if one wants to recreate rho from phi.
c ----------------------------------------------------------------------------

      if (iwhich == 9) then

#ifdef MPIPARALLEL
c  For parallel code, at this point, the data will never be transposed.
        iz = my_index*nz
        call phitorho(nx,ny,nz,a,kxsq,kysq,kzsq(iz),ikxmin,ikymin,1,1)
#else
c  Call for Serial version
#endif
        call phitorho(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,1,1)

      endif

c ----------------------------------------------------------------------------
c Tridiagonal solve in z, including kxsq and kysq on left hand side
c ----------------------------------------------------------------------------

      if (iwhich == 14) then

        norm = dz**2/eps0

#ifdef MPIPARALLEL
c  For parallel code, transpose is used to redistribute phi so each axial
c  line is on a single processor.
c       --- Rearrange phi array so each axial line is on a single processor.
        call transpose(nx,ny,nz,a,nx_tran,ny_tran)

c       --- Find location of transverse data relative to original dimensions.
        ix = mod(my_index*nx_tran,nx)
        iy = int(my_index*ny/real(nslaves))
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
        call tridiag(nx_tran,ny_tran-1,fullnz,nz2,a,norm,
     &               kxsq(ix),kysq(iy),kzsq,ikxm,ikym,0,zwork)

c       --- Rearrange phi again back to original arrangement.
        call transposei(nx,ny,nz,a,nx_tran,ny_tran)
#else
c  Serial version
        call tridiag(nx,ny,nz,nz2,a,norm,kxsq,kysq,kzsq,ikxmin,ikymin,1,zwork)
#endif

      endif

c ----------------------------------------------------------------------------
c  End of VPOIS3D
c ----------------------------------------------------------------------------

      return
      end
c=============================================================================
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
      subroutine vsftz(nx,ny,nz,a,norm,esx,esy,nxy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy,nxy
      real(kind=8):: a(0:nx,0:ny,0:nz),norm,zwork(2,0:nz,0:nx)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.
c nxy is the size of the transverse planes.

      integer(ISZ):: ix,iy,iz

      integer(ISZ):: naux1,naux2
      integer(ISZ):: prevnaux1,prevnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      save naux1,naux2,aux1,aux2
 
c Setup and make call to the ESSL routine for sine-transforms

c     if (isetup == 1) then
c       --- Calculate how much space is needed. Only allocate if it is
c       --- different from last time (deallocating if needed).
        prevnaux1 = naux1
        prevnaux2 = naux2
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

        if (prevnaux1 > 0 .and. prevnaux1 .ne. naux1) deallocate(aux1)
        if (prevnaux2 > 0 .and. prevnaux2 .ne. naux2) deallocate(aux2)
        if (prevnaux1 .ne. naux1) allocate(aux1(naux1))
        if (prevnaux2 .ne. naux2) allocate(aux2(naux2))
        call dcft(1,zwork,1,nz+1,zwork,1,nz+1,nz,nx+(1-esx),+1,1.,
     &            aux1,naux1,aux2,naux2)
c       return
c     endif
 
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

        call dcft(0,zwork,1,nz+1,zwork,1,nz+1,nz,nx+(1-esx),+1,1.,
     &            aux1,naux1,aux2,naux2)
        call vrpft2(zwork(1,0,0), zwork(2,0,0), nz, 2,   nx+(1-esx),2*(nz+1))

        do iz = 0, nz-1
          do ix = 0, nx-esx
            a(ix,iy  ,iz) = zwork(1,iz,ix)
            a(ix,iy+1,iz) = zwork(2,iz,ix)
          enddo
        enddo

      enddo

      return
      end
c=============================================================================
      subroutine vsftzi(nx,ny,nz,a,norm,esx,esy,nxy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy,nxy
      real(kind=8):: a(0:nx,0:ny,0:nz),norm,zwork(2,0:nz,0:nx)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.
c nxy is the size of the transverse planes.

      integer(ISZ):: ix,iy,iz

      integer(ISZ):: naux1,naux2
      integer(ISZ):: prevnaux1,prevnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      save naux1,naux2,aux1,aux2
 
c Setup and make call to the ESSL routine for sine-transforms

c     if (isetup == 1) then
c       --- Calculate how much space is needed. Only allocate if it is
c       --- different from last time (deallocating if needed).
        prevnaux1 = naux1
        prevnaux2 = naux2
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

        if (prevnaux1 > 0 .and. prevnaux1 .ne. naux1) deallocate(aux1)
        if (prevnaux2 > 0 .and. prevnaux2 .ne. naux2) deallocate(aux2)
        if (prevnaux1 .ne. naux1) allocate(aux1(naux1))
        if (prevnaux2 .ne. naux2) allocate(aux2(naux2))
        call dcft(1,zwork,1,nz+1,zwork,1,nz+1,nz,nx+(1-esx),-1,1.,
     &            aux1,naux1,aux2,naux2)
c       return
c     endif

      do iy = 0, ny-2+(1-esy), 2

        do iz = 0, nz-1
          do ix = 0, nx-esx
            zwork(1,iz,ix) = a(ix,iy  ,iz)
            zwork(2,iz,ix) = a(ix,iy+1,iz)
          enddo
        enddo

c       --- do the transforms

        call vrpfti2(zwork(1,0,0), zwork(2,0,0), nz, 2,   nx+(1-esx),2*(nz+1))
        call dcft(0,zwork,1,nz+1,zwork,1,nz+1,nz,nx+(1-esx),-1,1.,
     &            aux1,naux1,aux2,naux2)

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
      subroutine vsftz(nx,ny,nz,a,norm,esx,esy,nxy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy,nxy
      real(kind=8):: a(0:nx,0:ny,0:nz),norm,zwork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.
c nxy is the size of the transverse planes.

      integer(ISZ):: ix,iy,iz

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

      return
      end
c=============================================================================
      subroutine vsftzi(nx,ny,nz,a,norm,esx,esy,nxy,zwork)
      integer(ISZ):: nx,ny,nz,esx,esy,nxy
      real(kind=8):: a(0:nx,0:ny,0:nz),norm,zwork(2,0:nx,0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.
c
c esx and esy are the amount of extra space at the ends of the dimensions.
c nxy is the size of the transverse planes.

      integer(ISZ):: ix,iy,iz

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

      return
      end
c=============================================================================
#endif
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine attenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,esx,esy)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy
      real(kind=8):: a(0:nx,0:ny,0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

      do ikz = 0, nz-1
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*attx(ikx)*atty(iky)*attz(ikz)
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine unattenuate(nx,ny,nz,a,attx,atty,attz,ikxmin,ikymin,esx,esy)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy
      real(kind=8):: a(0:nx,0:ny,0:nz-1)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

      do ikz = 0, nz-1
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)/(attx(ikx)*atty(iky)*attz(ikz))
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine rhotophi(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,esx,esy)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy
      real(kind=8):: a(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz
      real(kind=8):: oneok

      do ikz = 0, nz-1
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            oneok = 1./(kxsq(ikx) + kysq(iky) + kzsq(ikz))
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*oneok
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine phitorho(nx,ny,nz,a,kxsq,kysq,kzsq,ikxmin,ikymin,esx,esy)
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,esx,esy
      real(kind=8):: a(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)

c This is put into a subroutine for the parallel field solver since the
c grid dimensions passed are not the true grid dimensions.

      integer(ISZ):: ikx,iky,ikz

      do ikz = 0, nz-1
        do iky = ikymin, ny-esy
          do ikx = ikxmin, nx-esx
            a(ikx,iky,ikz) = a(ikx,iky,ikz)*(kxsq(ikx) + kysq(iky) + kzsq(ikz))
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine tridiag(nx,ny,nz,nz2,a,norm,kxsq,kysq,kzsq,ikxmin,ikymin,esy,
     &                   zwork)
      integer(ISZ):: nx,ny,nz,nz2,ikxmin,ikymin,esy
      real(kind=8):: a(0:nx,0:ny,0:nz),norm
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
      real(kind=8)::  a(0:nx, 0:ny)
      real(kind=8):: ak(0:nx, 0:ny)
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

c  Error checking; so far, only possible error is non-power-of-2 nx,ny

      if (  ( nx .ne. 2**int( .5 + log(1.*nx) / log(2.) ) )
     & .or. ( ny .ne. 2**int( .5 + log(1.*ny) / log(2.) ) ) ) then
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
          call cosqx (a(0,0), xywork, w(icox), nx, ny, 1,-1)
          call cosqy (a(0,0), xywork, w(icoy), nx, ny, 1,-1)
        elseif (l2symtry) then
          call vsftx (a(0,0), xywork, w(icpx), w(icmx), nx, ny, 1)
          call cosqy (a(0,0), xywork, w(icoy), nx, ny, 1,-1)
        else
          call vsftx (a(0,0), xywork, w(icpx), w(icmx), nx, ny, 1)
          call vsfty (a(0,0), xywork, w(icpy), w(icmy), nx, ny, 1)
        endif

      endif

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

c       --- Normalize the array
        do iy = iystart,iyend
          a(:,iy) = a(:,iy)*norm
        enddo

c       --- Invoke the transform routines for the current z-slice
        if (l4symtry) then
          do iy = iystart,iyend,2
            call cosqx(a(0,iy),xywork(0,iy),w(icox),nx,2,0,-1)
          enddo
        else
          do iy = iystart,iyend,2
            call vsftx(a(0,iy),xywork(0,iy),w(icpx),w(icmx),nx,2,0)
          enddo
        endif

#ifdef MPIPARALLEL
        if (ny > xynppgroup .and. xynppgroup > 1)
     &    call parallelgatherall(a,(1+nx)*ny/xynppgroup,xynppgroup,1)
#endif

        if (l4symtry .or. l2symtry) then
          do ix = ixstart,ixend,2
            call cosqy(a(ix:ix+2,0:ny),xywork(ix:ix+2,0:ny),w(icoy),2,ny,0,-1)
          enddo
        else
          do ix = ixstart,ixend,2
            call vsfty(a(ix:ix+2,0:ny),xywork(ix:ix+2,0:ny),w(icpy),w(icmy),2,ny,0)
          enddo
        endif

      endif

c ----------------------------------------------------------------------------
c  Apply k-space filter; attenuate mode-by-mode
c ----------------------------------------------------------------------------

      if ((iwhich <= 0 .or. iwhich == 3) .and. maxval(abs(filt)) > 0.) then

        do iky = ikymin, ny-1
          do ikx = ixstart,ixend
           ak(ikx,iky) = ak(ikx,iky) * attx(ikx) * atty(iky)
          enddo
        enddo

      endif

c ----------------------------------------------------------------------------
c  Convert rhok to phik.  We should multitask this loop.
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 4) then

        do iky = ikymin, ny-1
          do ikx = ixstart,ixend
           ak(ikx,iky) = ak(ikx,iky)/(kxsq(ikx) + kysq(iky))
          enddo
        enddo

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform
c ----------------------------------------------------------------------------

      if (iwhich <= 0 .or. iwhich == 5) then

        if (l4symtry .or. l2symtry) then
c         --- Invoke the inverse transform routines for the current z-slice
          do ix = ixstart,ixend,2
            call cosqy(a(ix:ix+2,0:ny),xywork(ix:ix+2,0:ny),w(icoy),2,ny,0,1)
          enddo
        else
          do ix = ixstart,ixend,2
            call vsfty(a(ix:ix+2,0:ny),xywork(ix:ix+2,0:ny),w(icpy),w(icmy),2,ny,0)
          enddo
        endif

#ifdef MPIPARALLEL
        if (nx > xynppgroup .and. xynppgroup > 1) then
          atemp = transpose(a)
          call parallelgatherall(atemp,(1+ny)*nx/xynppgroup,xynppgroup,1)
          a = transpose(atemp)
        endif
#endif

        if (l4symtry) then
          do iy = iystart,iyend,2
            call cosqx(a(0,iy),xywork(0,iy),w(icox),nx,2,0,1)
          enddo
        else
          do iy = iystart,iyend,2
            call vsftx(a(0,iy),xywork(0,iy),w(icpx),w(icmx),nx,2,0)
          enddo
        endif

        if (l4symtry) then
c         --- factor of .25/lx for cosqx and .25/ly for cosqy
          norm = .0625/(lx*ly)
        else if (l2symtry) then
c         --- factor of 2./lx for vsftx and .25/ly for cosqy
          norm = .5/(lx*ly)
        else
          norm = 4./(lx*ly)
        endif

        a(nx,:) = 0.
        a(:,ny) = 0.
        if (.not. l4symtry) then
          a(0,:) = 0.
        endif
        if (.not. l2symtry .and. .not. l4symtry) then
          a(:,0) = 0.
        endif

        do iy = iystart,iyend
          a(:,iy) = a(:,iy)*norm
        enddo

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

      do iky = ikymin, ny-1
        do ikx = ikxmin, nx-1
          ak(ikx,iky) = ak(ikx,iky) / (attx(ikx) * atty(iky))
        enddo
      enddo

      endif

c ----------------------------------------------------------------------------
c  Multiply by k^2 factor, in contrast to the usual division.
c  This is useful if one wants to recreate rho from phi.
c ----------------------------------------------------------------------------

      if (iwhich == 9) then

      do iky = ikymin, ny-1
        do ikx = ikxmin, nx-1
          ak(ikx,iky) = ak(ikx,iky) * ( kxsq(ikx) + kysq(iky) )
        enddo
      enddo

      endif

c ----------------------------------------------------------------------------
c  End of VPOIS2D
c ----------------------------------------------------------------------------

      return
      end
c======================================================================
      subroutine vsftx (a, w, cp, cm, nx, ny, isetup)
      use Constant
      integer(ISZ):: nx,ny,isetup
      real(kind=8):: a(0:nx,0:ny), w(0:nx,0:ny)
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
      integer(ISZ):: prevnaux1,prevnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      save naux1,naux2,aux1,aux2

c Setup and make call to the ESSL routine for sine-transforms

      if (isetup == 1) then
c       --- Calculate how much space is needed. Only allocate if it is
c       --- different from last time (deallocating if needed).
        prevnaux1 = naux1
        prevnaux2 = naux2
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
        if (prevnaux1 > 0 .and. prevnaux1 .ne. naux1) deallocate(aux1)
        if (prevnaux2 > 0 .and. prevnaux2 .ne. naux2) deallocate(aux2)
        if (prevnaux1 .ne. naux1) allocate(aux1(naux1))
        if (prevnaux2 .ne. naux2) allocate(aux2(naux2))
        call dsinf(1,a,1,1+nx,a,1,1+nx,2*nx,ny,1.,aux1,naux1,aux2,naux2)
        return
      endif
 
      call dsinf(0,a,1,1+nx,a,1,1+nx,2*nx,ny,1.,aux1,naux1,aux2,naux2)

#endif

      return
      end
c======================================================================
      subroutine vsfty (a, w, cp, cm, nx, ny, isetup)
      use Constant
      integer(ISZ):: nx,ny,isetup
      real(kind=8):: a(0:nx,0:ny), w(0:nx,0:ny)
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
      integer(ISZ):: prevnaux1,prevnaux2
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      save naux1,naux2,aux1,aux2

c Setup and make call to the ESSL routine for sine-transforms
 
      if (isetup == 1) then
c       --- Calculate how much space is needed. Only allocate if it is
c       --- different from last time (deallocating if needed).
        prevnaux1 = naux1
        prevnaux2 = naux2
        if (2*ny <= 16384) then
          naux1 = 50000
        else
          naux1 = 20000+.60*ny*2
        endif
        if (2*ny <= 16384) then
          naux2 = 20000
        else
          naux2 = 20000+.64*ny*2
        endif
        if (ny >= 252) naux2 = naux2 + (ny+257)*(min(128, nx))
        if (prevnaux1 > 0 .and. prevnaux1 .ne. naux1) deallocate(aux1)
        if (prevnaux2 > 0 .and. prevnaux2 .ne. naux2) deallocate(aux2)
        if (prevnaux1 .ne. naux1) allocate(aux1(naux1))
        if (prevnaux2 .ne. naux2) allocate(aux2(naux2))

        call dsinf(1,a,1+nx,1,a,1+nx,1,2*ny,nx,1.,aux1,naux1,aux2,naux2)
        return
      endif
 
      call dsinf(0,a,1+nx,1,a,1+nx,1,2*ny,nx,1.,aux1,naux1,aux2,naux2)

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
      real(kind=8):: a(0:nx,0:ny),w(2,0:nx,ny/2),c(nx-1)

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
      integer(ISZ):: prevnaux1,prevnaux2
      integer(ISZ):: prevnx,prevny,previsign
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      data prevnx/0/,prevny/0/,previsign/0/
      save naux1,naux2,aux1,aux2,prevnx,prevny,previsign
 
c     --- If needed, allocate and initialize auxiliary arrays.
      if (nx .ne. prevnx .or. ny .ne. prevny .or. isign .ne. previsign) then
c       --- Calculate how much space is needed.
        naux1 = 20000
        if (nx > 2048) naux1 = naux1 + 2.28*nx
        naux2 = 20000
        if (nx > 2048) naux2 = naux2 + 2.28*nx

c       --- Only allocate if it is larger than last time (deallocating
c       --- if needed).
        if (naux1 > prevnaux1) then
          if (allocated(aux1)) deallocate(aux1)
          allocate(aux1(naux1))
        endif
        if (naux2 > prevnaux2) then
          if (allocated(aux2)) deallocate(aux2)
          allocate(aux2(naux2))
        endif

c       --- Initialize the auxiliary arrays
        call dcft(1,w,1,nx+1,w,1,nx+1,nx,ny/2,-isign,1.,aux1,naux1,aux2,naux2)

c       --- Save the values for this intialization
        prevnx = nx
        prevny = ny
        previsign = isign
        prevnaux1 = naux1
        prevnaux2 = naux2
      endif

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
        call dcft(0,w,1,nx+1,w,1,nx+1,nx,ny/2,+1,1.,aux1,naux1,aux2,naux2)
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
        call dcft(0,w,1,nx+1,w,1,nx+1,nx,ny/2,-1,1.,aux1,naux1,aux2,naux2)
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

      return
      end
c=======================================================================
      subroutine cosqy(a,w,c,nx,ny,isetup,isign)
      integer(ISZ):: nx,ny,isetup,isign
      real(kind=8):: a(0:nx,0:ny),w(2,0:ny,nx/2),c(ny-1)

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
      integer(ISZ):: prevnaux1,prevnaux2
      integer(ISZ):: prevnx,prevny,previsign
      real(kind=8),allocatable:: aux1(:),aux2(:)
      data naux1/0/,naux2/0/
      data prevnx/0/,prevny/0/,previsign/0/
      save naux1,naux2,aux1,aux2,prevnx,prevny,previsign
 
c Setup and make call to the ESSL routine for Fourier-transforms

c     --- If needed, allocate and initialize auxiliary arrays.
      if (nx .ne. prevnx .or. ny .ne. prevny .or. isign .ne. previsign) then
c       --- Calculate how much space is needed.
        naux1 = 20000
        if (ny > 2048) naux1 = naux1 + 2.28*ny
        naux2 = 20000
        if (ny > 2048) naux2 = naux2 + 2.28*ny

c       --- Only allocate if it is larger than last time (deallocating
c       --- if needed).
        if (naux1 > prevnaux1) then
          if (allocated(aux1)) deallocate(aux1)
          allocate(aux1(naux1))
        endif
        if (naux2 > prevnaux2) then
          if (allocated(aux2)) deallocate(aux2)
          allocate(aux2(naux2))
        endif

c       --- Initialize the auxiliary arrays
        call dcft(1,w,1,ny+1,w,1,ny+1,ny,nx/2,-isign,1.,aux1,naux1,aux2,naux2)

c       --- Save the values for this intialization
        prevnx = nx
        prevny = ny
        previsign = isign
        prevnaux1 = naux1
        prevnaux2 = naux2
      endif

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
        call dcft(0,w,1,ny+1,w,1,ny+1,ny,nx/2,+1,1.,aux1,naux1,aux2,naux2)
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
        call dcft(0,w,1,ny+1,w,1,ny+1,ny,nx/2,-1,1.,aux1,naux1,aux2,naux2)
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
      real(kind=8):: a(0:nx,0:ny),w(0:ny-1,0:nx),c(nx-1)

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
      real(kind=8):: a(0:nx,0:ny),w(0:nx-1,0:ny),c(ny-1)

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
     &                   xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                   scrtch,xywork,zwork,l2symtry,l4symtry)
      use CapMat3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c Sets up capacity matrix for 3d fieldsolver. Note that the capacity matrix
c is calculated using field solves which do not do any filtering.
c (Experience has shown that when filtering is used, the inverse capacity
c matrix is very ill-conditioned and so the capacity matrix becomes garbled.)


      real(kind=8):: wx,wy,wz,cx,cy
c     real(kind=8):: tt1,tt2
      integer(ISZ):: ix,iy,iz,i,j,info

c find capacity matrix

c     --- find inverse capacity matrix
      do i=1,nc3d

c       --- start with zero phi
c       call zeroarry(phi(0,0,0),(nx+1)*(ny+1)*(1+nz))
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

c       --- Solve for fields, piece-by-piece so that the filtering can be
c       --- skipped.
        call vpois3d(2,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)
        call vpois3d(4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)
        call vpois3d(5,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)

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
c     call ssifa(cmat3d,nc3d,nc3d,kpvt3d,info)
c     if (info .ne. 0) then
c       call remark("ERROR: Capacity matrix is singular - it is likely that")
c       call remark("two or more conductor points are within the same grid cell.")
c       call kaboom(0)
c     endif
c     call ssidi(cmat3d,nc3d,nc3d,kpvt3d,tt1,tt2,scrtch,001)

c     --- These routines are from LAPACK.
#if WORDSIZE == 64
      call ssytrf("U",nc3d,cmat3d,nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),info)
#else
      call dsytrf("U",nc3d,cmat3d,nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),info)
#endif
      if (info .ne. 0) then
        call remark("ERROR: Capacity matrix is singular - it is likely that")
        call remark("two or more conductor points are within the same grid cell.")
c       print*,info
        call kaboom(0)
      endif
c     --- pcond3d is passed as workspace.
#if WORDSIZE == 64
      call ssytri("U",nc3d,cmat3d,nc3d,kpvt3d,pcond3d,info)
#else
      call dsytri("U",nc3d,cmat3d,nc3d,kpvt3d,pcond3d,info)
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
     &                    filt,xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,
     &                    xmmin,ymmin,zmmin,scrtch,xywork,zwork,
     &                    l2symtry,l4symtry)
      use CapMat3d
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

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

        call vpois3d(1,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)

        call cmset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &               xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &               scrtch,xywork,zwork,l2symtry,l4symtry)
        if (iwhich == 1) then
          call copyarry(rho(0,0,0),phi(0,0,0),(nx+1)*(ny+1)*(1+nz))
        endif

      elseif (iwhich > 1) then

c       --- Make call to do the specialized action requested.
        call vpois3d(iwhich,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry)

      endif

c     --- If a full field solve was not desired, return.
      if (iwhich > 0) return

c     --- Do the first field solve (note that rho must have been copied into
c     --- phi for this to work)
      call vpois3d(-1,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
      call copyarry(rho(0,0,0),phi(0,0,0),(nx+1)*(ny+1)*(1+nz))

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
      call vpois3d(-1,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
     &                     xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                     scrtch,xywork,zwork,l2symtry,l4symtry)
      use CapMat3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c Sets up capacity matrix for 3d fieldsolver. Note that the capacity matrix
c is calculated using field solves which do not do any filtering.
c (Experience has shown that when filtering is used, the inverse capacity
c matrix is very ill-conditioned and so the capacity matrix becomes garbled.)


      real(kind=8):: wx,wy,cx,cy
c     real(kind=8):: tt1,tt2
      integer(ISZ):: ix,iy,iz,i,j,info

c find capacity matrix

c     --- Make sure capacity matrix is fully allocated.
      nc3dz = nz
#ifdef MPIPARALLEL
      nc3dz2 = nz - 1
#else
      nc3dz2 = nz/2
#endif
      call gchange("CapMat3d",0)

c     --- find inverse capacity matrix
      do i=1,nc3d

c       --- start with zero phi
c       call zeroarry(phi(0,0,0),(nx+1)*(ny+1)*(1+nz))
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
        call vpois3d(12,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                  filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                  l2symtry,l4symtry)
        call vpois3d( 4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                  filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                  l2symtry,l4symtry)
        call vpois3d(13,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                  filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                  l2symtry,l4symtry)

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
c       call ssifa(cmat3d(1,1,iz),nc3d,nc3d,kpvt3d,info)
c       if (info .ne. 0) then
c         call remark("ERROR: Capacity matrix is singular - it is likely that")
c         call remark("two or more conductor points are within the same grid cell.")
c         call kaboom(0)
c       endif
c       call ssidi(cmat3d(1,1,iz),nc3d,nc3d,kpvt3d,tt1,tt2,scrtch,001)

c       --- These routines are from LAPACK.
#if WORDSIZE == 64
        call ssytrf("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),
     &              info)
#else
        call dsytrf("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,scrtch,(1+nx)*(1+ny),
     &              info)
#endif
        if (info .ne. 0) then
          call remark("ERROR: Capacity matrix is singular - it is likely that")
          call remark("two or more conductor points are within the same grid cell.")
c         print*,info
          call kaboom(0)
        endif
c       --- pcond3d is passed as workspace.
#if WORDSIZE == 64
        call ssytri("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,pcond3d,info)
#else
        call dsytri("U",nc3d,cmat3d(1,1,iz),nc3d,kpvt3d,pcond3d,info)
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
     &                      filt,xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,
     &                      xmmin,ymmin,zmmin,scrtch,xywork,zwork,
     &                      l2symtry,l4symtry)
      use CapMat3d
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(*)
      real(kind=8):: xlen,ylen,zlen,dx,dy,dz,xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c Applies the capacity matrix for 3d code electrostatic field solver.
c Multiplies phi by cap to get induced rho.
c Assumes that beam's rho has been set, and copied into phi.
c When filtering is being done, field solve is inexact since filtering is
c applied to the induced charge.
c For the serial version, advantage is taken by realizing the the capacity
c matrix for the kz plane iz is the same as that for kz plane nz-iz. The
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

        call vpois3d(1,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)

        call cmkzset3d(phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                 xlen,ylen,zlen,nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                 scrtch,xywork,zwork,l2symtry,l4symtry)
        if (iwhich == 1) then
          call copyarry(rho(0,0,0),phi(0,0,0),(nx+1)*(ny+1)*(1+nz))
        endif

      elseif (iwhich > 1) then

c       --- Make call to do the specialized action requested.
        call vpois3d(iwhich,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry)

      endif

c     --- If a full field solve was not desired, return.
      if (iwhich > 0) return

c     --- Do the first field solve (note that rho must have been copied into
c     --- phi for this to work). The inverse transform in z is skipped.
      call vpois3d( 2,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 3,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d(13,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
c     --- the correct thing, since the capacity matrix for plane nz-iz is the
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
      call copyarry(rho(0,0,0),phi(0,0,0),(nx+1)*(ny+1)*(1+nz))

c     --- Transform in z
      call vpois3d(10,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
      call vpois3d(12,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 5,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
     &                  filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,
     &                  l2symtry,l4symtry)
      use Pipe3d
      integer(ISZ):: iwhich,nx,ny,nz
      character(*):: pipeshpe
      real(kind=8):: rho(0:nx,0:ny,0:nz), phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen
      logical(ISZ):: l2symtry,l4symtry
 
c  External interface to capacity matrix field solver.
c  Passes external and local variables to field solver routines.
c  If iwhich is '0' or '1', also allocates local variables


c  Allocate local variables
      pipenz = nz
#ifdef MPIPARALLEL
      pipenz2 = nz - 1
#else
      pipenz2 = nz/2
#endif
      if (iwhich == 0 .or. iwhich == 1) then
        if (pipen == 0) pipen = 2*nx
        pipe8th = pipen/8 + 1
        call gallot("Pipe3d",0)
      endif

c  Call field solver
      call vpipe3d(iwhich,pipeshpe,rho,phi,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &             xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,
     &             piperadi,pipen,pipenz,pipenz2,pipe8th,pipex,pipey,
     &             pipephi,pipeq,cap3d,kpvt,l2symtry,l4symtry)

      return
      end
c=======================================================================
      subroutine pipest3d(pipeshpe,cap,phi,kxsq,kysq,kzsq,attx,atty,attz, 
     &                    filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,
     &                    piperadi,pipen,pipenz2,pipe8th,pipex,pipey,cap3d,kpvt,
     &                    l2symtry,l4symtry)
      use Constant
      integer(ISZ):: nx,ny,nz
      character(*):: pipeshpe
      real(kind=8):: cap(2*(nx+ny),*)
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen, piperadi
      integer(ISZ):: pipen,pipenz2,pipe8th
      real(kind=8):: pipex(pipen+1), pipey(pipen+1), cap3d(pipen,pipe8th,0:pipenz2)
      real(kind=8):: kpvt(pipen)
      logical(ISZ):: l2symtry,l4symtry

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
c       call zeroarry(phi(0,0,0),(nx+1)*(ny+1)*(nz+1))
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
        call vpois3d (12, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &   attx, atty, attz, filt, xlen, ylen, zlen, nx, ny, nz, scrtch,
     &   xywork,zwork, 0, l2symtry,l4symtry)
c       --- divide by k squared
        call vpois3d (4, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &   attx, atty, attz, filt, xlen, ylen, zlen, nx, ny, nz, scrtch,
     &   xywork,zwork, 0, l2symtry,l4symtry)
c       --- inverse transform in x and y to find phi(kz)
        call vpois3d (13, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &   attx, atty, attz, filt, xlen, ylen, zlen, nx, ny, nz, scrtch,
     &   xywork,zwork, 0, l2symtry,l4symtry)
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
        call ssifa(cap,2*(nx+ny),pipen,kpvt,info)
        call ssidi(cap,2*(nx+ny),pipen,kpvt,tt1,tt2,scrtch,001)

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
     &                   filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,
     &                   piperadi,pipen,pipenz,pipenz2,pipe8th,pipex,pipey,
     &                   pipephi,pipeq,cap3d,kpvt,l2symtry,l4symtry)
      integer(ISZ):: iwhich,nx,ny,nz
      character(*):: pipeshpe
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz)
      real(kind=8):: filt(5,3),scrtch(*),xywork(*),zwork(2,0:nx,0:nz)
      real(kind=8):: xlen,ylen,zlen,piperadi
      integer(ISZ):: pipen,pipenz,pipenz2,pipe8th
      real(kind=8):: pipex(pipen+1),pipey(pipen+1),cap3d(pipen,pipe8th,0:pipenz2)
      real(kind=8):: kpvt(pipen)
      real(kind=8):: pipephi(pipen,0:pipenz-1),pipeq(pipen,0:pipenz-1)
      logical(ISZ):: l2symtry,l4symtry

c applies capacity matrix, for round pipe electrostatic solver.
c multiplies phi by cap to get rho.
c assumes that rho has been set, and copied into phi.

      real(kind=8):: xx,yy
      integer(ISZ):: ix,iy,iz,iiz,i,j

c initialize capacity matrix and arrays for poisson solve
      if (iwhich == 1 .or. iwhich == 0) then
        call vpois3d(1,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                 l2symtry,l4symtry)
        call pipest3d(pipeshpe,phi,phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,
     &                piperadi,pipen,pipenz2,pipe8th,pipex,pipey,cap3d,kpvt,
     &                l2symtry,l4symtry)
        if (iwhich == 0) then
c         --- Copy rho back into phi since phi was scrambled during the
c         --- calculation of the matrix.
          call copyarry(rho,phi(0,0,0),(nx+1)*(ny+1)*(nz+1))
        endif
      endif

      if (iwhich == 1) return

c do first partial field solve (without z transform at end)
      call vpois3d( 2,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 3,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d(13,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
      call copyarry(rho,phi(0,0,0),(nx+1)*(ny+1)*(nz+1))

c     --- transform in z
      call vpois3d(10,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
      call vpois3d(12,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 4,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)
      call vpois3d( 5,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry)

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
      subroutine setcap3d(ncndpts,quadx,quady,quadz,quadcap,kpvt,nx,ny,nz,
     &                    phi,work,xywork,zwork,kxsq,kysq,kzsq,attx,atty,attz,
     &                    xlen,ylen,zlen,filt,l2symtry,l4symtry)
      integer(ISZ):: ncndpts,nx,ny,nz
      real(kind=8):: quadx(ncndpts,4), quady(ncndpts,4), quadz(ncndpts)
      real(kind=8):: quadcap(4*ncndpts,4*ncndpts)
      integer(ISZ):: kpvt(4*ncndpts)
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: work(nx+ny-4),xywork(*),zwork(2,0:nx,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3)
      real(kind=8):: xlen, ylen, zlen
      logical(ISZ):: l2symtry,l4symtry

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

      real(kind=8):: xx,yy,zz,dx,dy,dz
      integer(ISZ):: ix,iy,iz,i,j,ic,info

      dx = xlen/nx
      dy = ylen/ny
      dz = zlen/nz

c  calculate inverse cap matrix

c  first calculate capacities for one conductor
      do i=1,ncndpts
c       call zeroarry(phi,(nx+1)*(ny+1)*(nz+1))
        phi = 0.
c  deposit unit charge on grid
        ix = quadx(i,1)/dx + nx/2
        xx = quadx(i,1)/dx + nx/2 - ix
        iy = quady(i,1)/dy + ny/2
        yy = quady(i,1)/dy + ny/2 - iy
        iz = quadz(i)/dz + nz/2 - 4
        zz = quadz(i)/dz + nz/2 - 4 - iz
        phi(ix  ,iy  ,iz  ) = (1. - xx)*(1. - yy)*(1. - zz)
        phi(ix+1,iy  ,iz  ) = (     xx)*(1. - yy)*(1. - zz)
        phi(ix  ,iy+1,iz  ) = (1. - xx)*(     yy)*(1. - zz)
        phi(ix+1,iy+1,iz  ) = (     xx)*(     yy)*(1. - zz)
        phi(ix  ,iy  ,iz+1) = (1. - xx)*(1. - yy)*(     zz)
        phi(ix+1,iy  ,iz+1) = (     xx)*(1. - yy)*(     zz)
        phi(ix  ,iy+1,iz+1) = (1. - xx)*(     yy)*(     zz)
        phi(ix+1,iy+1,iz+1) = (     xx)*(     yy)*(     zz)
c  do 3-d field solve
        call vpois3d (-1, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &   attx, atty, attz, filt, xlen, ylen, zlen, nx, ny, nz, work,
     &   xywork,zwork, 0, l2symtry,l4symtry)
c  fill i'th column of inverse cap matrix
        do ic = 1,4
          do j=1,ncndpts
            ix = quadx(j,ic)/dx + nx/2
            xx = quadx(j,ic)/dx + nx/2 - ix
            iy = quady(j,ic)/dy + ny/2
            yy = quady(j,ic)/dy + ny/2 - iy
            iz = quadz(j)/dz + nz/2 - 4
            zz = quadz(j)/dz + nz/2 - 4 - iz
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
      call ssifa(quadcap,4*ncndpts,4*ncndpts,kpvt,info)

      return
      end
c==========================================================================
      subroutine capfs(iwhich,nx,ny,nz,phi,rho,
     &                 xlen,ylen,zlen,kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                 work,workz,xywork,zwork,xmmax,pipeshpe,periinz,
     &                 l2symtry,l4symtry)
      use Capvars
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz), rho(0:nx,0:ny,0:nz)
      real(kind=8):: work(nx+ny-4), workz(*),xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: xlen, ylen, zlen, filt(5,3), xmmax
      character(*):: pipeshpe
      logical(ISZ):: periinz
      logical(ISZ):: l2symtry,l4symtry

c  Capacity matrix field solve routine
c  input:
c  pipeshpe   shape of conductors, default is circlular, "hyp" is hyperbolic
c  rest of variables have the standard meanings

c  using common block since needs too many arguments to be handled by basis

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
        call vpois3d (1, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &                attx, atty, attz, filt, xlen, ylen, zlen,
     &                nx, ny, nz, work, xywork, zwork, 0, l2symtry,l4symtry)
c    set up conductor point arrays and gchange array sizes
        call cndctr3d(nx,ny,xmmax,dx,dy,quadradi,quadcent,
     &                ncndpts, quadx, quady, quadz, quadv, nzpts, nendquad,
     &                nzquad, qdslclen, pipeshpe, loadquad)
        call gchange("Capvars",0)
c    set up capacity matrix
        call setcap3d(ncndpts,quadx,quady,quadz,quadcap,kkkk,nx,ny,nz,
     &                phi,work,xywork,zwork,kxsq,kysq,kzsq,
     &                attx,atty,attz,xlen,ylen,
     &                zlen,filt,l2symtry,l4symtry)
        if (iwhich == 0) then
c         --- Copy rho back into phi since phi was scrambled during the
c         --- calculation of the matrix.
          call copyarry(rho,phi,(nx+1)*(ny+1)*(nz+1))
        endif
      endif

c  return if only being initialized
      if (iwhich == 1) return

c  do 3-d field solve
      call vpois3d(0,phi(0,0,0),phi(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &    filt,xlen,ylen,zlen,nx,ny,nz,work,xywork,zwork, 0,l2symtry,l4symtry)

c  zero induced charge array
c     call zeroarry(quadrho,numends*4*ncndpts)
      quadrho = 0.

c loop over the quadrupole focusing lenses
      do iend=1,numends

c    find z's of conducting points relative to grid
        do iquad=1,ncndpts
          workz(iquad) = (quadz(iquad) + quadend(iend))*dzi
          if (workz(iquad) < 0) workz(iquad) = 0.
          if (workz(iquad) >= nz) workz(iquad) = nz - 1.
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
        call ssisl(quadcap,4*ncndpts,4*ncndpts,kkkk,quadrho(1,1,iend))

      enddo

c  reload rho back into phi
      call copyarry(rho,phi,(nx+1)*(ny+1)*(nz+1))

c  if periodic in z, zero out iz=nz plane so that plane will only have
c  image charges (plane 0 will have both image charges and beam charge)
      if (periinz) then
        do ix=0,nx
          do iy=0,ny
            phi(ix,iy,nz) = 0.
          enddo
        enddo
      endif

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
          if (workz(iquad) >= nz) then
            workz(iquad) = nz - 1.
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

c Make phi periodic by adding planes iz=0 and nz and placing result back into
c both planes.  (Note that plane 0 has beam charge and additional image charge 
c while plane nz has only image charge.)
      if (periinz) then
        do ix=0,nx
          do iy=0,ny
            phi(ix,iy,0) = phi(ix,iy,0) + phi(ix,iy,nz)
            phi(ix,iy,nz) = phi(ix,iy,0)
          enddo
        enddo
      endif

c  finished with conductors, so...
c  now do full 3-d field solve, with induced charges
      call vpois3d (0, phi(0,0,0), phi(0,0,0), kxsq, kysq, kzsq,
     &              attx, atty, attz, filt, xlen, ylen, zlen,
     &              nx, ny, nz, work, xywork, zwork, 0,l2symtry,l4symtry)

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
     &                    zmmin,zgrid,nz,dz,
     &                    numends,nendmax,quadend,quadvlt,
     &                    vshift,quadcent,quadradi)
      integer(ISZ):: nquad,nz,numends,nendmax
      real(kind=8):: quadzs(0:nquad), quadze(0:nquad), quadde(0:nquad)
      real(kind=8):: quadvx(0:nquad), quadvy(0:nquad)
      real(kind=8):: quadvlt(nendmax), quadend(nendmax), vshift(nendmax)
      real(kind=8):: zmmin,zgrid,dz,quadcent,quadradi

c  Find starts of quads in grid, and put quad voltage into quadvlt
c  Also sets array vshift
c  dependent upon quad variables

      real(kind=8):: offset
      integer(ISZ):: iz,iq,iqprev

c  get rest of quads within mesh
      iq = 0
      iqprev = -1
      numends = 0
      do iz=1,nz-1
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
      subroutine vcap3d (iwhich, rho, phi, kxsq, kysq, kzsq, 
     &   attx, atty, attz, filt, xlen, ylen, zlen, nx, ny, nz, 
     &   scrtch, xywork, zwork, xmmax, zmmin, zgrid, pipeshpe, periinz,
     &   l2symtry,l4symtry)
      use Capvars
      use Lattice
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: kxsq(0:nx-1), kysq(0:ny-1), kzsq(0:nz)
      real(kind=8):: attx(0:nx-1), atty(0:ny-1), attz(0:nz)
      real(kind=8):: filt(5,3), scrtch(*), xywork(*), zwork(2,0:nx,0:nz)
      real(kind=8):: xlen, ylen, zlen
      real(kind=8):: xmmax,zmmin,zgrid
      character(*):: pipeshpe
      logical(ISZ):: periinz
      logical(ISZ):: l2symtry,l4symtry
 
c  External interface to capacity matrix field solver.
c  Passes external and local variables to field solver routines.
c  If iwhich is '0' or '1', also allocates local variables


      real(kind=8):: dz,offset
      integer(ISZ):: iq

      dz = zlen/nz

c  find quad starts
      call findqdnd(nquad,quadzs,quadze,quadde,quadvx,quadvy,
     &              zmmin,zgrid,nz,dz,
     &              numends,nendmax,quadend,quadvlt,vshift,quadcent,quadradi)

c  Allocate local variables
      if (iwhich == 0 .or. iwhich == 1) then
        if (loadquad) then
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
      call capfs(iwhich,nx,ny,nz,phi,rho,xlen,ylen,zlen,kxsq,kysq,kzsq,
     &           attx,atty,attz,filt,scrtch,scrtch,xywork,zwork,
     &           xmmax,pipeshpe,periinz,l2symtry,l4symtry)

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c
c PSOR3D field solver
c
c===========================================================================
      subroutine psor3df(iwhich,nx,ny,nz,phi,rho,phi1d,rho1d,rstar,dx,dy,dz,
     &                  xmmin,ymmin,zmmin,zbeam,zgrid,linbends,
     &                  l2symtry,l4symtry,scrtch,izfsmin,izfsmax)
      use PSOR3d
      use Constant
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: iwhich,nx,ny,nz,izfsmin,izfsmax
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
      real(kind=8):: phi1d(*),rho1d(*),scrtch(*)
      real(kind=8):: dx,dy,dz,xmmin,ymmin,zmmin,zbeam,zgrid
      logical(ISZ):: linbends
      logical(ISZ):: l2symtry,l4symtry

c This routine does a 3D field solve using point SOR with even-odd (red-black)
c ordering.  It makes call to the routine cond_pot which sets voltages at
c conductor grid points to the desired value.  The routine also allows for a
c bent beam-pipe.
c
c The tranverse boundaries can either be held constant, have zero normal
c derivative, or be periodic.  When BOUNDXY is zero, the boundaries are held
c constant, when 1, they have zero normal derivative, and when 2, the
c boundaries are periodic.
c
c The longitudinal boundaries can either be held constant, have zero normal
c derivative, or be periodic.  When BOUND0 or BOUNDNZ is zero, the boundaries
c are held constant, when 1, they have zero normal derivative, and when 2, the
c boundaries are periodic.
c
c Note that loops over all directions assume that nx, ny, and nz are even.
c
c The arrangement of the loops was done to increase performance.  The entire
c grid is looped over as if it were a 1D array, ignoring boundaries.
c The boundaries are then reset, the previous value was destroyed.
c
c The routine CONDBNDY is called to add sub-grid level boundary placement. The
c potential of the points near the surface is saved since they are written
c over in the 1D loop.
c
c izfsmin and izfsmax are added so that a field solve can be done with only
c part of the array.  izfsmin is the minimum value of iz for which the field
c solve is done.  The plane izfsmin acts as bound0.  izfsmax is the maximum
c value of iz.  The plane izfsmax acts as boundnz.
c
c rstar(-1) and rstar(nz+1) are set based on the axial boundary conditions.


      real(kind=8):: boundarrx(0:nx,2,0:nz)
      real(kind=8):: boundarry(0:ny,2,0:nz)
      real(kind=8):: phisave(int((nx+1)*(ny+1)*(nz+1)/(.9*nx-1)))
      real(kind=8):: rdx2,rdy2,rdz2,rdel,reps0,const
      real(kind=8):: rdx2c,rdy2c,rdz2c,reps0c,spm1
      real(kind=8):: ecdelmxi(necndbdy),ecdelmyi(necndbdy),ecdelmzi(necndbdy)
      real(kind=8):: ecdelpxi(necndbdy),ecdelpyi(necndbdy),ecdelpzi(necndbdy)
      real(kind=8):: ocdelmxi(nocndbdy),ocdelmyi(nocndbdy),ocdelmzi(nocndbdy)
      real(kind=8):: ocdelpxi(nocndbdy),ocdelpyi(nocndbdy),ocdelpzi(nocndbdy)
      integer(ISZ):: nsorerr
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,iimz,iipz,izl,izr,parity
      integer(ISZ):: i,k,ii,ie,ix,iy,iz,ic
      real(kind=8):: rs,x,r,nsorerri
      integer(ISZ):: s_parity,e_parity
      character(72):: errline

c Macro for common arguments to the psor_loop calls
#define LOOP_ARGS iimz,iipz,phi1d,rho1d,nxy,rdx2c,rdy2c,rdz2c,spm1,linbends, \
         rstar,dx,rdx2,rdy2,rdz2,rdel,const,scrtch

c Make sure that the transverse grid dimensions are even.
      if (nx .ne. 2*int(nx/2) .or. ny .ne. 2*int(ny/2)) then
        call remark("ERROR in psor3d: nx and ny must be even numbers")
        return
      endif

c If iwhich equals zero or one, do initialization
      if (iwhich == 0 .or. iwhich == 1)
     &  call psorinit(nx,ny,nz,dx,dy,dz,l2symtry,l4symtry)

c If doing initialization only, then exit
      if (iwhich == 1) return

c Set the array islctype and determine the points that make up the conductor.
c This takes extra time so should not be done if the grid is not moving
c in the lab frame.  Set gridmode to 1 to avoid this call.
      if (gridmode == 0 .or. iwhich == -2) then
c       --- Zero the number of conductor points
        ncond = 0
        necndbdy = 0
        nocndbdy = 0
c       --- Generate new conductors.
        call setcndtr(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                l2symtry,l4symtry)
      endif

c Set up for error calculations.  Set nsorerr so the sample points used to
c check for convergence are spread throughout the grid.  The expression
c (.9*nx-1) was chosen somewhat randomly. Note that the same expression
c must be used as the one above for the shape of phisave.
      nsorerr = (nx+1)*(ny+1)*(nz+1)/(.9*nx-1)
      nsorerri = 1./nsorerr

c Set temporary variables (these are used to increase performance)
      rdx2 = dx**(-2)
      rdy2 = dy**(-2)
      rdz2 = dz**(-2)
      rdel = rdz2/(rdx2 + rdy2 + rdz2)
      reps0 = 1./eps0
      const = sorparam*0.5/(rdx2 + rdy2 + rdz2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdz2c = rdz2*const
      reps0c = reps0*const
      spm1 = 1. - sorparam

c Set indices for 1d arrays used in the seven point finite difference
c form of Poisson's equation.
      nxy = (nx+1)*(ny+1)
      iimx = -1
      iipx = +1
      iimy = -nx-1
      iipy = +nx+1
      iimz = -nxy
      iipz = +nxy

c Set longitudinal indices for 1d phi array.  '1' is added to izl and izr
c since phi1d is passed as phi(0,0,-1).  If using Dirichlet boundary
c conditions, do not solve for the potential on the end planes.
      izl = izfsmin + 1
      izr = izfsmax + 1
      if (bound0 == 0) izl = izl + 1
      if (boundnz== 0) izr = izr - 1

c Do work that is bend dependent.
c Preset rho to increase performance.  (It is then reset at end of routine).
c For bends, also set up scrtch arrays.
      if (.not. linbends) then
c       --- no bends  
        do ii=nxy*izfsmin+1,nxy*izfsmax+nxy
          rho1d(ii) = rho1d(ii)*reps0c
        enddo
      else 
c       --- bends, also include curvature corrections.  Comment: Timing tests 
c       --- show that the use of 1d array is slightly faster than a 3d array.  
        do ii=nxy*izfsmin+1,nxy*izfsmax+nxy
          i = mod(ii-1,nx+1) 
          k = (ii-1)/nxy
          rs = rstar(k) 
          x  = xmmin + i*dx 
          r  = rs + x
c         rearranged to reduce divides
c         rho1d(ii) = rho1d(ii)*(rs/r)*reps0c/
c    &                ( 1. + (x/r)*((x/r)-2.)*rdel )
          rho1d(ii) = rho1d(ii)*rs*r*reps0c/
     &                ( r*r + x*(x-2.*r)*rdel )
        enddo
c       --- Fill scratch array with x values so it can be looked up
c       --- in the bent beam loop instead of calculated.
        do ii = 1, nxy
          scrtch(ii) = xmmin + mod(ii-1,nx+1)*dx
        enddo
c       --- Change rstar if using Nuemann boundary conditions
        if (bound0 == 1) rstar(-1) = rstar(1)
        if (boundnz == 1) rstar(nz+1) = rstar(nz-1)
      endif 

c Precalculate the reciprocals of the subgrid data.
      ecdelpxi = 1./ecdelpx(1:necndbdy)
      ecdelpyi = 1./ecdelpy(1:necndbdy)
      ecdelpzi = 1./ecdelpz(1:necndbdy)
      ecdelmxi = 1./ecdelmx(1:necndbdy)
      ecdelmyi = 1./ecdelmy(1:necndbdy)
      ecdelmzi = 1./ecdelmz(1:necndbdy)
      ocdelpxi = 1./ocdelpx(1:nocndbdy)
      ocdelpyi = 1./ocdelpy(1:nocndbdy)
      ocdelpzi = 1./ocdelpz(1:nocndbdy)
      ocdelmxi = 1./ocdelmx(1:nocndbdy)
      ocdelmyi = 1./ocdelmy(1:nocndbdy)
      ocdelmzi = 1./ocdelmz(1:nocndbdy)

c If doing Chebyshev accel, set sorparam to one for first iteration.
      if (lchebshv) sorparam = 1.

c Zero phi out at the edge of the mesh. This clears out any remaining values
c of phi from conductor points at the edge of the mesh. This is only done
c when gridmode==0 (i.e. when new conductors are generated each step),
c when lzerophiedge is true, when using Dirichlet boundary conditions, and
c if there are conductors. The logical allows the user to turn this off for
c special cases.
      if (lzerophiedge .and. gridmode == 0 .and.
     &    (ncond > 0 .or. necndbdy > 0 .or. nocndbdy > 0)) then
        if (boundxy == 0) then
          do iz=0,nz
            do ix=0,nx
              phi(ix,0,iz) = 0.
              phi(ix,ny,iz) = 0.
            enddo
            do iy=0,ny
              phi(0,iy,iz) = 0.
              phi(nx,iy,iz) = 0.
            enddo
          enddo
        endif
        if (bound0 == 0) then
          do iy=0,ny
            do ix=0,nx
              phi(ix,iy,0) = 0.
            enddo
          enddo
        endif
        if (boundnz == 0) then
          do iy=0,ny
            do ix=0,nx
              phi(ix,iy,nz) = 0.
            enddo
          enddo
        endif
      endif

c Put desired potential onto conductors in phi array.
      call cond_pot(nx,ny,nz,phi)

c main loop for sor iterations,  Errors are checked every isorerr iterations.

      sorerror = 1.
      soriter = 0
      do while (sorerror > sortol .and. soriter < sormaxit)

c     --- Do isorerr iterations, checking error only on the last iteration.
c     --- This is done so the error isn't checked on every iteration.
        do ie = 1,isorerr
          soriter = soriter + 1

c         --- If calculating error, save some of current phi into error array.
c         --- This takes a sampling of phi to check for the change in this 
c         --- iteration.  The explicit integer and an addition of .5 in 
c         --- the 1-d array index calculation is for a bug fix on the CRAY.
c         --- The factor of 1. in the integer causes the calculcation to be done
c         --- with real(kind=8):: numbers so that there is not an integer overflow.

          if (ie == isorerr) then
            do i=1,nsorerr
              phisave(i) = phi1d(nxy*izl +
     &                int(1.*i*nxy*(izr-izl+1)*nsorerri+.5) )
            enddo
          endif

c         --- Save values on the transverse boundaries.
c         --- Both even and odd points are saved for all transverse boundaries.
          do iz=izl-1,izr-1
            do ix=0,nx
              boundarrx(ix,1,iz) = phi(ix,0,iz)
              boundarrx(ix,2,iz) = phi(ix,ny,iz)
            enddo
            do iy=0,ny
              boundarry(iy,1,iz) = phi(0,iy,iz)
              boundarry(iy,2,iz) = phi(nx,iy,iz)
            enddo
          enddo

c         --- Save values just outside conductor surfaces.
c         --- Again, both even and odd points are saved.
          if (lcndbndy) then
            do ic=1,necndbdy
              ecndpvph(ic) = phi(iecndx(ic),iecndy(ic),iecndz(ic))
            enddo
            do ic=1,nocndbdy
              ocndpvph(ic) = phi(iocndx(ic),iocndy(ic),iocndz(ic))
            enddo
          endif

c do loop to cover even and odd points
#ifdef MPIPARALLEL
          s_parity = mod(izslave(my_index),2)
          e_parity = mod(s_parity+1,2)
#else
          s_parity = 0
          e_parity = 1
#endif
          do parity=s_parity,e_parity,e_parity-s_parity

#ifdef MPIPARALLEL
c*********************Parallel Routine*********************
c           --- exchange phi on the boundary
            call exchange_phi(nx,ny,izfsmax-izfsmin,phi(0,0,izfsmin-1),
     &                        bound0,boundnz,0)
c**********************************************************
#endif

c           --- Set guard planes in z appropriately for the boundary conditions
c           --- For periodic boundaries, this need not be done for the
c           --- parallel version since it automatically happens in exchange_phi

            if (bound0 == 1) then
              do ii = (izl-1)*nxy+1+mod(parity+izl-1,2),(izl-1)*nxy+nxy
                phi1d(ii) = phi1d(ii + 2*nxy)
              enddo
#ifndef MPIPARALLEL
            elseif (bound0 == 2) then
              do ii = (izl-1)*nxy+1+mod(parity+izl-1,2),(izl-1)*nxy+nxy
                phi1d(ii) = phi1d(ii + nz*nxy)
              enddo
#endif
            endif
            if (boundnz == 1) then
              do ii = (izr+1)*nxy+1+mod(parity+izr+1,2),(izr+1)*nxy+nxy
                phi1d(ii) = phi1d(ii - 2*nxy)
              enddo
#ifndef MPIPARALLEL
            elseif (boundnz == 2) then
              do ii = (izr+1)*nxy+1+mod(parity+izr+1,2),(izr+1)*nxy+nxy
                phi1d(ii) = phi1d(ii - nz*nxy)
              enddo
#endif
            endif

c           --- Loop over the rest of the array. Boundary points are calculated
c           --- too, even though the equation is wrong.  They are recalculated
c           --- later.  The case with and without bends are seperated since
c           --- the loops are done differently.  Without bends, the whole array
c           --- can be processed in one loop as a 1d array.  With bends, an
c           --- outer loop over z is needed.
            if (.not. linbends) then 
c             --- no bends 
              call psor_loop(parity,1,0,0,1,izl*nxy+2,(izr+1)*nxy,2,
     &                       iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
            else 
c             --- bends
              call psor_loop(parity,1,izl,izr,1,2,nxy,2,
     &                       iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
            endif 

c           --- Transverse boundaries     
c           --- Restore only even or odd boundary points to previous value
c           --- since only want to restore the values changed from the 1d
c           --- loop.  This automatically takes care of Dirichlet boundaries.
            do iz=izl-1,izr-1
              do ix=mod(iz+parity,2),nx,2
                phi(ix,0,iz) = boundarrx(ix,1,iz)
                phi(ix,ny,iz) = boundarrx(ix,2,iz)
              enddo
              do iy=mod(iz+parity,2),ny,2
                phi(0,iy,iz) = boundarry(iy,1,iz)
                phi(nx,iy,iz) = boundarry(iy,2,iz)
              enddo
            enddo

            if (boundxy == 1 .or. l2symtry .or. l4symtry) then
c             --- surfaces at ix=0, ix=nx, iy=0, and iy=ny
c             --- if only 2-fold, apply to surfaces at iy=0
c             --- if only 4-fold, apply to surfaces at ix=0 and iy=0
c             --- if also boundxy=1, then apply to all transverse surfaces
              call psor_loop(parity,1,izl,izr,1,
     &                       2,nx,2,iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
              if (boundxy == 1 .or. l4symtry) then
                call psor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,
     &                         2*(nx+1),-iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
c               --- line at transverse edge (ix=0, iy=0)
                call psor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,
     &                         2*nxy,-iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
              endif
c             --- now do planes at ix=nx and iy=ny
              if (boundxy == 1) then
                call psor_loop(parity,1,izl,izr,1,(nx+1)*ny+2,(nx+1)*(ny+1)-1,2,
     &                         iimx,iipx,iimy,-iipy,
     &LOOP_ARGS)
                call psor_loop(parity,nx+1,izl,izr,1,2*(nx+1),(nx+1)*ny,
     &                         2*(nx+1),iimx,-iipx,iimy,iipy,
     &LOOP_ARGS)
c               --- lines at other transverse edges
                call psor_loop(parity,nxy,izl,izl,1,nx+1,nxy*(izr-izl)+nx+1,
     &                         2*nxy,iimx,-iipx,-iimy,iipy,
     &LOOP_ARGS)
                call psor_loop(parity,nxy,izl,izl,1,nxy-nx,nxy*(izr-izl)+nxy-nx,
     &                         2*nxy,-iimx,iipx,iimy,-iipy,
     &LOOP_ARGS)
                call psor_loop(parity,nxy,izl,izl,1,nxy,nxy*(izr-izl)+nxy,
     &                         2*nxy,iimx,-iipx,iimy,-iipy,
     &LOOP_ARGS)
              endif
            else if (boundxy == 2) then
c             --- surfaces at ix=0, ix=nx, iy=0, and iy=ny
              call psor_loop(parity,1,izl,izr,1,
     &                       2,nx,2,iimx,iipx,(nx+1)*(ny-1),iipy,
     &LOOP_ARGS)
              call psor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,
     &                       2*(nx+1),nx-1,iipx,iimy,iipy,
     &LOOP_ARGS)
              do iz=izl-1,izr-1
                do ix=1+mod(iz+parity+1,2),nx-1,2
                  phi(ix,ny,iz) = phi(ix,0,iz)
                enddo
                do iy=1+mod(iz+parity+1,2),ny-1,2
                  phi(nx,iy,iz) = phi(0,iy,iz)
                enddo
              enddo
c             --- lines at transverse edges
              call psor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,
     &                       2*nxy,nx-1,iipx,(nx+1)*(ny-1),iipy,
     &LOOP_ARGS)
              do iz=izl-1+mod(izl-1+parity,2),izr-1,2
                phi(nx,0,iz) = phi(0,0,iz)
                phi(0,ny,iz) = phi(0,0,iz)
                phi(nx,ny,iz) = phi(0,0,iz)
              enddo
            endif
c           --- end of transverse boundaries

c           --- Apply altered difference equation to the points near the
c           --- surface of the conductor boundaries.
            if (lcndbndy) then
              if (parity == s_parity .and. necndbdy > 0) then
                call condbndy(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,sorparam,
     &                     ncndmax,necndbdy,
     &                     ecdelmx(1),ecdelmy(1),ecdelmz(1),
     &                     ecdelpx(1),ecdelpy(1),ecdelpz(1),
     &                     ecdelmxi(1),ecdelmyi(1),ecdelmzi(1),
     &                     ecdelpxi(1),ecdelpyi(1),ecdelpzi(1),
     &                     iecndx(1),iecndy(1),iecndz(1),ecvolt(1),ecndpvph(1),
     &                     bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                     izfsmin,izfsmax,dx,dy,dz,icndbndy)
              else if (nocndbdy > 0) then
                call condbndy(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,sorparam,
     &                     ncndmax,nocndbdy,
     &                     ocdelmx(1),ocdelmy(1),ocdelmz(1),
     &                     ocdelpx(1),ocdelpy(1),ocdelpz(1),
     &                     ocdelmxi(1),ocdelmyi(1),ocdelmzi(1),
     &                     ocdelpxi(1),ocdelpyi(1),ocdelpzi(1),
     &                     iocndx(1),iocndy(1),iocndz(1),ocvolt(1),ocndpvph(1),
     &                     bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                     izfsmin,izfsmax,dx,dy,dz,icndbndy)
              endif
            endif

c           --- Put desired potential onto conductors in phi array.
            call cond_pot(nx,ny,nz,phi)

c           --- do Chebyshev acceleration
            if (lchebshv) then
              sorparam = 1./(1.-pjacobi**2*sorparam*.25)
              const = sorparam*0.5/(rdx2 + rdy2 + rdz2)
              rdx2c = rdx2*const
              rdy2c = rdy2*const
              rdz2c = rdz2*const
              reps0c = reps0*const
              spm1 = 1. - sorparam
            endif

c         --- end of loop over even and odd points
          enddo

c         --- Calculate error and terminate iterations if error is within
c         --- tolerance.
c         --- The explicit integer and an addition of .5 in 
c         --- the 1-d array index calculation is for a bug fix on the CRAY.
c         --- The factor of 1. in the integer causes the calculcation to be
c         --- done with real(kind=8):: numbers so that there is not an integer overflow.

          if (ie == isorerr) then
            sorerror = 0.
            do i=1,nsorerr
              phisave(i)=abs(phi1d(nxy*izl +
     &              int(1.*i*nxy*(izr-izl+1)*nsorerri+.5)) - phisave(i))
              if (phisave(i) > sorerror) sorerror = phisave(i)
            enddo
#ifdef MPIPARALLEL
c*********************Parallel Routine*********************
c           --- calculate global sorerror
            call parallelmaxrealarray(sorerror,1)
c           write(STDERR,*) "global error = ",sorerror
c**********************************************************
#endif
          endif

c       --- end of group of iterations
        enddo 

c end of all iterations
      end do

c printout warning that maximum iterations have been exceeded
      if (soriter > sormaxit) then
        call remark("PSOR: *** Maximum sor iterations exceeded ***")
      endif

c write out number of iterations taken
      if (soriter > 0) then
        write(errline,20) soriter,sorerror
  20    format ( "PSOR: After ",i4," iterations, the error was ",1pe11.3)
        call remark(errline)
      endif

c Exchange phi at the boundaries, 0 and nz, one last time to ensure
c consistency across processors. This has the side benefit of setting
c planes -1 and nz+1 for the periodic case. Also, exchange phi in the
c guard planes, -1 and nz+1 for all processors.
#ifdef MPIPARALLEL
c*********************Parallel Routine*********************
c     --- exchange phi on the boundary
      call exchange_phi(nx,ny,izfsmax-izfsmin,phi(0,0,izfsmin-1),
     &                  bound0,boundnz,0)
      call exchange_phi(nx,ny,izfsmax-izfsmin,phi(0,0,izfsmin-1),
     &                  bound0,boundnz,-1)
c**********************************************************
#endif

c Set planes nz=-1 and nz=nz+1 based on type of boundary condition.
c This is not done in the case of a partial field solve.
c Separate code is needed for the serial and parallel versions, since
c in parallel, only the first and last processors set the guard cell
c and only for Dirichlet and Nuemann boundary conditions. The periodic
c cases are skipped by the parallel version since that exchange_phi call
c above handles it.
#ifndef MPIPARALLEL
      if (izfsmin == 0) then
        if (bound0 == 0) then
          phi(:,:,-1) = phi(:,:,0)
        else if (bound0 == 1) then
          phi(:,:,-1) = phi(:,:,1)
        else if (bound0 == 2) then
          phi(:,:,-1) = phi(:,:,nz-1)
        endif
      endif
      if (izfsmax == nz) then
        if (boundnz == 0) then
          phi(:,:,nz+1) = phi(:,:,nz)
        else if (boundnz == 1) then
          phi(:,:,nz+1) = phi(:,:,nz-1)
        else if (boundnz == 2) then
          phi(:,:,nz+1) = phi(:,:,1)
        endif
      endif
#else
      if (izfsmin == 0) then
        if (bound0 == 0 .and. my_index == 0) then
          phi(:,:,-1) = phi(:,:,0)
        else if (bound0 == 1 .and. my_index == 0) then
          phi(:,:,-1) = phi(:,:,1)
        endif
      endif
      if (izfsmax == nz) then
        if (boundnz == 0 .and. my_index == nslaves-1) then
          phi(:,:,nz+1) = phi(:,:,nz)
        else if (boundnz == 1 .and. my_index == nslaves-1) then
          phi(:,:,nz+1) = phi(:,:,nz-1)
        endif
      endif
#endif

c unset rho
      reps0c = 1./reps0c
      if (.not. linbends) then 
c       --- no bends 
        do ii=nxy*(izfsmin)+1,nxy*(izfsmax+1)
          rho1d(ii) = rho1d(ii)*reps0c
        enddo
      else 
c       --- bends.  Comment: Timing tests show that the use of a 1d array 
c       --- is slightly faster than a 3d array. 
        do ii=nxy*(izfsmin)+1,nxy*(izfsmax+1)
          i = mod(ii-1,nx+1)
          k = (ii-1)/nxy
          rs = rstar(k) 
          x  = xmmin + i*dx 
          r  = rs + x
          rho1d(ii) = rho1d(ii)/rs*reps0c*( r + x*((x/r)-2.)*rdel )
        enddo
      endif 

      return
      end
c===========================================================================
      subroutine psor_loop(i_parity,s_parity,izl,izr,izstep,
     &                     ii1,ii2,iistep,iimx,iipx,iimy,iipy,iimz,iipz,
     &                     phi1d,rho1d,nxy,rdx2c,rdy2c,rdz2c,spm1,linbends,
     &                     rstar,dx,rdx2,rdy2,rdz2,rdel,const,xx)
      real(kind=8):: phi1d(*),rho1d(*),rstar(*),xx(*)
      integer(ISZ):: i_parity,s_parity,izl,izr,izstep,nxy,ii1,ii2,iistep
      integer(ISZ):: iimx,iipx,iimy,iipy,iimz,iipz
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,dx,rdx2,rdy2,rdz2,rdel,const
      logical(ISZ):: linbends

c This routine provides the loops over the phi array which solves the 
c iterative equation.  The loops are set up to do an outer loop over
c the axial dimension (the third dimension, or z) and an inner loop
c over the transverse dimensions.  By passing the proper loop limits when
c possible, the entire array can be looped over in the inner loop, increasing
c performance.  The structure of having two seperate loops makes this routine
c usable in the calculation of boundary conditions on the grid edges.  The
c structure is also amenable to the bent field solver.
c
c The variable parity is used to ensure that the loops maintain the parity as
c passed by i_parity.  The value of parity is not necessarily the same as
c i_parity.  When parity is calculated, the expression (iz-1) is used since
c phi1d is passed starting at iz=-1 and the array location phi(0,0,0) is
c considered even.  The expression (ii1-1) is used since the 1d array starts
c at an index of one, but the first element is even.
c
c The index of rho1d is (ii-nxy) since rho has dimensions rho(,,0:nz), whereas
c the value of ii is based on an array size of (,,-1:nz+1).
c
c Inputs:
c    i_parity  :input parity
c    s_parity  :parity multiplier
c    izl       :left, or starting z index
c    izr       :right, or ending z index
c    izstep    :step for z loop
c    ii1       :starting index for inner loop
c    ii2       :ending index for inner loop
c    iistep    :step for inner loop
c    iimx      :change in 1d index to array element phi(ix-1,iy  ,iz  )
c    iipx      :change in 1d index to array element phi(ix+1,iy  ,iz  )
c    iimy      :change in 1d index to array element phi(ix  ,iy-1,iz  )
c    iipy      :change in 1d index to array element phi(ix  ,iy+1,iz  )
c    iimz      :change in 1d index to array element phi(ix  ,iy  ,iz-1)
c    iipz      :change in 1d index to array element phi(ix  ,iy  ,iz+1)
c    xx        :value of x used with bends (see scrtch in psor3d)
c The rest of the inputs are the same as in psor3d.
c
c The if test (iistep == 2) is needed for the workstation version.  Without
c the loop stride explicitly given, there is about a 50% performance reduction
c on HIF, an HP 735 workstation.

      integer(ISZ):: iz,ii,parity
      real(kind=8):: x,rs,r,xfact,del2bndc,dh,dxi

      if (.not. linbends) then 
        do iz=izl,izr,izstep
          parity = mod(i_parity*s_parity + (iz-1)*nxy + ii1-1,2)
          if (iistep == 2) then
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, 2
              phi1d(ii) = rho1d(ii-nxy) +
     &                  (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                  (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                  (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                  spm1*phi1d(ii)
            enddo
          else
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
              phi1d(ii) = rho1d(ii-nxy) +
     &                  (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                  (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                  (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                  spm1*phi1d(ii)
            enddo
          endif
        enddo
      else
c       --- bends, loop over z slices.  Comment: timing tests indicate that 
c       --- use of a 1d array is faster than a 3d array.   
        dxi = 1./dx
        do iz=izl,izr,izstep
          parity = mod(i_parity*s_parity + (iz-1)*nxy + ii1-1,2)
          rs  = rstar(iz+1)
          if ( abs(rstar(iz  )) > LARGEPOS*1.e-6 .and. 
     &         abs(rstar(iz+2)) > LARGEPOS*1.e-6 ) then 
c           --- current z-slice not in bend or bordering a bend, do not 
c           --- include curvature terms from Poisson's equation    
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
              phi1d(ii) = rho1d(ii-nxy) +
     &            (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2c +
     &            (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2c +
     &            (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2c +
     &            spm1*phi1d(ii)
            enddo  
          else  
c           --- current z-slice in a bend and near border of a bend, include 
c           --- curvature terms from Poisson's equation.  Note: if bordering 
c           --- a bend, the dh/dz "jump term" should be included.  The 
c           --- current method of treating this jump term may be inaccurate. 
c           --- r is set to 1 over (rs+x) to reduce the number of divides.
            dh = 1./rstar(iz+2) - 1./rstar(iz)
            if (abs(dh) > SMALLPOS .and. abs(rs) < LARGEPOS*1.e-6) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                r = 1./(rs + x)
                xfact    = (x*r)*((x*r)-2.) 
                del2bndc = const/( 1.+xfact*rdel )
                phi1d(ii) = rho1d(ii-nxy) + del2bndc*( 
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2+
     &              .5*(phi1d(ii+iipx)-phi1d(ii+iimx))*r*dxi + 
     &              xfact*(phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2 - 
     &              .25*(rs*r)**3*x*(phi1d(ii+iipz) - phi1d(ii+iimz))*dh*rdz2) +
     &              spm1*phi1d(ii)
              enddo 

c           --- current z-slice is near border of a bend only
            elseif (abs(dh) > SMALLPOS) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                phi1d(ii) = rho1d(ii-nxy) + const*(
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2 -
     &              .25*x*(phi1d(ii+iipz) - phi1d(ii+iimz))*dh*rdz2) +
     &              spm1*phi1d(ii)
              enddo

c           --- current z-slice is in a bend only
            elseif (abs(rs) < LARGEPOS*1.e-6) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                r = 1./(rs + x)
                xfact    = (x*r)*((x*r)-2.)
                del2bndc = const/( 1.+xfact*rdel )
                phi1d(ii) = rho1d(ii-nxy) + del2bndc*(
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2+
     &              .5*(phi1d(ii+iipx) - phi1d(ii+iimx))*r*dxi +
     &              xfact*(phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2) +
     &              spm1*phi1d(ii)
              enddo
            endif
          endif 
        enddo 
      endif 

      return
      end
c===========================================================================
      subroutine cond_pot(nx,ny,nz,phi)
      use PSOR3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
c Set conductor points to the desired potential.
      integer(ISZ):: ic

c sets phi to the desired potential on conductors
      do ic = 1,ncond
        phi(ixcond(ic),iycond(ic),izcond(ic)) = condvolt(ic)
      enddo

      return
      end
c===========================================================================
      subroutine condbndy(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,srp,
     &                    ncndmax,ncndbdy,
     &                    cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,cdelpz,
     &                    cdelmxi,cdelmyi,cdelmzi,cdelpxi,cdelpyi,cdelpzi,
     &                    icndx,icndy,icndz,cvolt,cndpvph,
     &                    bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                    izfsmin,izfsmax,dx,dy,dz,icndbndy)
      integer(ISZ):: nx,ny,nz,ncndmax,ncndbdy
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz)
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,srp
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax), cdelmz(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax), cdelpz(ncndmax)
      real(kind=8):: cdelmxi(ncndbdy), cdelmyi(ncndbdy), cdelmzi(ncndbdy)
      real(kind=8):: cdelpxi(ncndbdy), cdelpyi(ncndbdy), cdelpzi(ncndbdy)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax), icndz(ncndmax)
      real(kind=8):: cvolt(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
      integer(ISZ):: bound0,boundnz,boundxy,izfsmin,izfsmax
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: dx,dy,dz
      integer(ISZ):: icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c NOTE that rdx2cos and rdy2cos are rdx2c and rdy2c over sorparam.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: rdx2cos,rdy2cos,rdz2cos,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: denomi
      real(kind=8):: dxm,dxp,dym,dyp,dzm,dzp
      real(kind=8):: rhofix,rhovolume,dxyzi
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1,icndzp1,icndzm1

c -------------------------------------------------------------------------
c EGUN style interpolation at conductor boundaries
      if (icndbndy == 1) then

        rdx2cos = rdx2c/srp
        rdy2cos = rdy2c/srp
        rdz2cos = rdz2c/srp

c       --- loop over points near surface of rods
        do ic = 1,ncndbdy
c         --- set temporaries for boundaries
          icndxp1 = icndx(ic) + 1
          icndxm1 = icndx(ic) - 1
          icndyp1 = icndy(ic) + 1
          icndym1 = icndy(ic) - 1
          icndzp1 = icndz(ic) + 1
          icndzm1 = icndz(ic) - 1

          if (icndxm1 == -1 .and. l4symtry) icndxm1 = 1
          if (icndxm1 == -1 .and. boundxy == 0) cycle
          if (icndxm1 == -1 .and. boundxy == 1) icndxm1 = 1
          if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
          if (icndxp1 == nx+1 .and. boundxy == 0) cycle
          if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
          if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1

          if (icndym1 == -1 .and. l2symtry) icndym1 = 1
          if (icndym1 == -1 .and. l4symtry) icndym1 = 1
          if (icndym1 == -1 .and. boundxy == 0) cycle
          if (icndym1 == -1 .and. boundxy == 1) icndym1 = 1
          if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
          if (icndyp1 == ny+1 .and. boundxy == 0) cycle
          if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
          if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1

          if (icndzm1 == izfsmin-1 .and. bound0  == 0) cycle
          if (icndzm1 == izfsmin-1 .and. bound0  == 1) icndzm1 = izfsmin+1
          if (icndzm1 == izfsmin-1 .and. bound0  == 2) icndzm1 = izfsmax-1
          if (icndzp1 == izfsmax+1 .and. boundnz == 0) cycle
          if (icndzp1 == izfsmax+1 .and. boundnz == 1) icndzp1 = izfsmax-1
          if (icndzp1 == izfsmax+1 .and. boundnz == 2) icndzp1 = izfsmin+1

c         --- set temporaries with initial values
          pxm = phi(icndxm1  ,icndy(ic),icndz(ic))
          pxp = phi(icndxp1  ,icndy(ic),icndz(ic))
          pym = phi(icndx(ic),icndym1  ,icndz(ic))
          pyp = phi(icndx(ic),icndyp1  ,icndz(ic))
          pzm = phi(icndx(ic),icndy(ic),icndzm1  )
          pzp = phi(icndx(ic),icndy(ic),icndzp1  )
          denom = 1.
          dxm = 1.
          dxp = 1.
          dym = 1.
          dyp = 1.
          dzm = 1.
          dzp = 1.

c         --- the point lower in x is inside the conductor
          if (cdelmx(ic) < 1.) then
            pxm = cvolt(ic)*cdelmxi(ic)
            denom = denom + (1.-cdelmx(ic))*cdelmxi(ic)*rdx2cos
            dxm = cdelmx(ic)
          endif
c         --- the point higher in x is inside the conductor
          if (cdelpx(ic) < 1.) then
            pxp = cvolt(ic)*cdelpxi(ic)
            denom = denom + (1.-cdelpx(ic))*cdelpxi(ic)*rdx2cos
            dxp = cdelpx(ic)
          endif
c         --- the point lower in y is inside the conductor
          if (cdelmy(ic) < 1.) then
            pym = cvolt(ic)*cdelmyi(ic)
            denom = denom + (1.-cdelmy(ic))*cdelmyi(ic)*rdy2cos
            dym = cdelmy(ic)
          endif
c         --- the point higher in y is inside the conductor
          if (cdelpy(ic) < 1.) then
            pyp = cvolt(ic)*cdelpyi(ic)
            denom = denom + (1.-cdelpy(ic))*cdelpyi(ic)*rdy2cos
            dyp = cdelpy(ic)
          endif
c         --- the point lower in z is inside the conductor
          if (cdelmz(ic) < 1.) then
            pzm = cvolt(ic)*cdelmzi(ic)
            denom = denom + (1.-cdelmz(ic))*cdelmzi(ic)*rdz2cos
            dzm = cdelmz(ic)
          endif
c         --- the point lower in z is inside the conductor
          if (cdelpz(ic) < 1.) then
            pzp = cvolt(ic)*cdelpzi(ic)
            denom = denom + (1.-cdelpz(ic))*cdelpzi(ic)*rdz2cos
            dzp = cdelpz(ic)
          endif
c         --- calculate the new phi based on the boundary conditions
c         rhovolume = (dxp+dxm)*(dyp+dym)*(dzp+dzm)*0.125
          rhovolume = 1.
          phi(icndx(ic),icndy(ic),icndz(ic)) =
     &        (rho(icndx(ic),icndy(ic),icndz(ic))*rhovolume +
     &         (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c + (pzm+pzp)*rdz2c)/denom +
     &         spm1*cndpvph(ic)
        enddo

c -------------------------------------------------------------------------
      else if (icndbndy == 2) then
c EBC style interpolation, uses non-centered finite difference
        rhofix = srp*(2./(dx*dx) + 2./(dy*dy) + 2./(dz*dz))
        dxyzi = 1./(dx*dy*dz)

c       --- loop over points near surface of rods
        do ic = 1,ncndbdy
c         --- set temporaries for boundaries
          icndxp1 = icndx(ic) + 1
          icndxm1 = icndx(ic) - 1
          icndyp1 = icndy(ic) + 1
          icndym1 = icndy(ic) - 1
          icndzp1 = icndz(ic) + 1
          icndzm1 = icndz(ic) - 1

          if (icndxm1 == -1 .and. l4symtry) icndxm1 = 1
          if (icndxm1 == -1 .and. boundxy == 0) cycle
          if (icndxm1 == -1 .and. boundxy == 1) icndxm1 = 1
          if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
          if (icndxp1 == nx+1 .and. boundxy == 0) cycle
          if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
          if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1

          if (icndym1 == -1 .and. l2symtry) icndym1 = 1
          if (icndym1 == -1 .and. l4symtry) icndym1 = 1
          if (icndym1 == -1 .and. boundxy == 0) cycle
          if (icndym1 == -1 .and. boundxy == 1) icndym1 = 1
          if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
          if (icndyp1 == ny+1 .and. boundxy == 0) cycle
          if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
          if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1

          if (icndzm1 == izfsmin-1 .and. bound0  == 0) cycle
          if (icndzm1 == izfsmin-1 .and. bound0  == 1) icndzm1 = izfsmin+1
          if (icndzm1 == izfsmin-1 .and. bound0  == 2) icndzm1 = izfsmax-1
          if (icndzp1 == izfsmax+1 .and. boundnz == 0) cycle
          if (icndzp1 == izfsmax+1 .and. boundnz == 1) icndzp1 = izfsmax-1
          if (icndzp1 == izfsmax+1 .and. boundnz == 2) icndzp1 = izfsmin+1

c         --- set temporaries with initial values
          pxm = phi(icndxm1  ,icndy(ic),icndz(ic))
          pxp = phi(icndxp1  ,icndy(ic),icndz(ic))
          pym = phi(icndx(ic),icndym1  ,icndz(ic))
          pyp = phi(icndx(ic),icndyp1  ,icndz(ic))
          pzm = phi(icndx(ic),icndy(ic),icndzm1  )
          pzp = phi(icndx(ic),icndy(ic),icndzp1  )
          dxm = dx
          dxp = dx
          dym = dy
          dyp = dy
          dzm = dz
          dzp = dz

c         --- the point lower in x is inside the conductor
          if (cdelmx(ic) < 1.) then
            pxm = cvolt(ic)
            dxm = cdelmx(ic)*dx
          endif
c         --- the point higher in x is inside the conductor
          if (cdelpx(ic) < 1.) then
            pxp = cvolt(ic)
            dxp = cdelpx(ic)*dx
          endif
c         --- the point lower in y is inside the conductor
          if (cdelmy(ic) < 1.) then
            pym = cvolt(ic)
            dym = cdelmy(ic)*dy
          endif
c         --- the point higher in y is inside the conductor
          if (cdelpy(ic) < 1.) then
            pyp = cvolt(ic)
            dyp = cdelpy(ic)*dy
          endif
c         --- the point lower in z is inside the conductor
          if (cdelmz(ic) < 1.) then
            pzm = cvolt(ic)
            dzm = cdelmz(ic)*dz
          endif
c         --- the point lower in z is inside the conductor
          if (cdelpz(ic) < 1.) then
            pzp = cvolt(ic)
            dzp = cdelpz(ic)*dz
          endif
c         --- calculate the new phi based on the boundary conditions
          denomi = 1./(2./(dxm*dxp) + 2./(dym*dyp) + 2./(dzm*dzp))
c         rhovolume = (dxp+dxm)*(dyp+dym)*(dzp+dzm)*0.125*dxyzi
          rhovolume = 1.
          phi(icndx(ic),icndy(ic),icndz(ic))  =
     &        (rho(icndx(ic),icndy(ic),icndz(ic))*rhofix*rhovolume +
     &         pxm*2./(dxm*(dxm+dxp)) + pxp*2./(dxp*(dxm+dxp)) +
     &         pym*2./(dym*(dym+dyp)) + pyp*2./(dyp*(dym+dyp)) +
     &         pzm*2./(dzm*(dzm+dzp)) + pzp*2./(dzp*(dzm+dzp)))*srp*denomi +
     &        spm1*cndpvph(ic)
        enddo

      endif
c -------------------------------------------------------------------------

      return
      end
c===========================================================================
c===========================================================================
      subroutine psorinit(nx,ny,nz,dx,dy,dz,l2symtry,l4symtry)
      use PSOR3d
      use InGen
      use Lattice
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry
c Allocate various arrays


      integer(ISZ):: nc,iq

#ifdef MPIPARALLEL
c Fix longitudinal boundary conditions for the parallel version.
      if (my_index > 0) bound0 = 0
      if (my_index < nslaves-1) boundnz = 0

c     --- zparity is needed to keep track of even and odd subgrid points.
c     --- The parity is relative to the full grid. In the serial case,
c     --- it defaults to 0.
      zparity = mod(izslave(my_index),2)
#endif

c Now allocate all of the needed space.
      call gchange("PSOR3d",0)

      return
      end
c===========================================================================
c===========================================================================
      subroutine setcndtr(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                    l2symtry,l4symtry)
      use PSOR3d
      use InGen, ONLY:fstype
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      integer(ISZ):: nx,ny,nz
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c  Find type of slice, and put quad voltage into slicvx,slicvy.
c  This routine is dependent upon quad variables.
c  If quadrr is negative for that electrode a hyperbolic shape is asssumed
c  rather than a circular cylinder. (modified by I. Haber)
c  Fuzz is used to ensure proper evaluation of numeric comparisons in if's.

c     --- Call the appropriate routine.
      if(solvergeom==RZgeom) then
        call setcndtrrz(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                  l2symtry,l4symtry)
      else if (fstype == 7 .or. fstype == 11) then
        call setcndtr3dmg(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                    l2symtry,l4symtry)
c     else if(fstype==11) then
c       call setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
c    &                    l2symtry,l4symtry)
      else
        call setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                  l2symtry,l4symtry)
      end if

      return
      end
c===========================================================================
      subroutine setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                      l2symtry,l4symtry)
      use Constant
      use Surface_of_Rev
      use PSOR3d
      use Lattice
      use InjectVars
      integer(ISZ):: nx,ny,nz
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c  Find type of slice, and put quad voltage into slicvx,slicvy.
c  This routine is dependent upon quad variables.
c  If quadrr is negative for that electrode a hyperbolic shape is asssumed
c  rather than a circular cylinder. (modified by I. Haber)
c  Fuzz is used to ensure proper evaluation of numeric comparisons in if's.

      integer(ISZ):: islctype(-1:nz+1)
      real(kind=8):: slicvx(-1:nz+1),slicvy(-1:nz+1)
      integer(ISZ):: ix,iz,iq,ix_axis,iy_axis,nc,nceo
      integer(ISZ):: izp1,izm1
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax,ichange,id,ij
      integer(ISZ):: ix_near,nxinj,nyinj
      real(kind=8):: zz,offset,qdcenter,qs,qe,rs,re,xmax0,ymax0
      real(kind=8):: rodrrsq,rodmxsq,rminsq,rrr,delz,rodx,rody,rodap,delz_p
      real(kind=8):: vvv,aper,rmax,prevcent,dipocent,delx
      real(kind=8):: xcent,ycent
      real(kind=8):: fuzz,zfuzz
      logical(ISZ):: lnotice,lz_in_rod,lz_in_plate

c     --- fuzz is used to prevent the distance from a point to the surface
c     --- from being too small or zero.  Any point within fuzz of the surface
c     --- is counted to be within the conductor
      fuzz = 1.e-5
      zfuzz = dz*fuzz

c     --- loop over z planes to fill islctype array
c     --- Note that the loop ranges from -1 to nz+1, since any quads just
c     --- outside of the grid need to be found in order to apply the subgrid
c     --- scheme near the z edges of the mesh.
c     --- This refers directly to the quad arrays rather than the cquad arrays
c     --- since cquad arrays are only intended to be accessed by particles.
c     --- Key for islctype:
c     ---   0: drift
c     ---   1: quads on x axis
c     ---   2: quads on both
c     ---   3: quads on y axis
c     ---   4: plate left of quad
c     ---   5: plate right of quad
      iq = 0
      nc = 0
      nceo = 0
      lnotice = .false.
      do iz=-1,nz+1
        islctype(iz) = 0

c       --- find location in lattice array
        zz = zgrid + zmmin + iz*dz

c       --- Get id of quadrupole element
        call getquadid(zz,offset,iq,1)

c       --- set up temp variables for locating quadrupoles
        qdcenter = 0.5*(quadzs(iq) + quadze(iq)) + offset
        slicvx(iz) = quadvx(iq)
        slicvy(iz) = quadvy(iq)

c       --- calculate (in lab frame) start and end of quad
        qs = qdcenter - 0.5*(quadrl(iq) + quadgl(iq))
        qe = qdcenter + 0.5*(quadrl(iq) + quadgl(iq))

c       --- calculate (in lab frame) start and end rods
        if (quadgp(iq) < 0.) then
          rs = qs + quadgl(iq) + qdelglx(iq)
          re = qe - quadgl(iq) - qdelgly(iq)
        else
          rs = qs + quadgl(iq) + qdelgly(iq)
          re = qe - quadgl(iq) - qdelglx(iq)
        endif

c       --- Only use this quad is the aperture is > 0, and the rod radius
c       --- and length are also > 0.
        if (quadap(iq) > 0 .and. quadrr(iq) > 0. .and. quadrl(iq) > 0.) then
          if (qs <= zz+zfuzz .and. zz+zfuzz < rs .and. zz < re) then
c           --- rod on one axis, gap in the other
            if (quadgp(iq) > 0.) then
              islctype(iz) = 1
            else
              islctype(iz) = 3
            endif
          elseif (rs <= zz+zfuzz .and. zz-zfuzz <= re) then
c           --- rods on both axis
            islctype(iz) = 2
          elseif (re < zz-zfuzz .and. zz-zfuzz <= qe .and. zz > rs)then
c           --- rod on one axis, gap in the other
            if (quadgp(iq) > 0.) then
              islctype(iz) = 3
            else
              islctype(iz) = 1
            endif
          endif
        endif

c       --- find end plates
        if (lplates .and. quadpa(iq) > 0. .and.
     &      quadrr(iq) > 0. .and. quadrl(iq) > 0.) then
          if (qs-quadpw(iq)-qdelpwl(iq) <= zz+zfuzz .and.
     &        zz-zfuzz <= qs) then
            islctype(iz) = 4
          elseif (qe <= zz+zfuzz .and.
     &             zz-zfuzz <= qe+quadpw(iq)+qdelpwr(iq)) then
            islctype(iz) = 5
          endif
        endif

c       --- Check for case where both hard edged and self-consistent
c       --- models are being used for the same quadrupole.  Print notice if
c       --- this is the case.
        if (islctype(iz) == 2 .and. quadde(iq) .ne. 0.) lnotice = .true.

c       --- Estimate amount of space needed.
        if (islctype(iz) == 1 .or. islctype(iz) == 3) then
          nc = nc + 2*pi*(quadrr(iq)**2 - (quadrr(iq)-2.*max(dx,dy))**2)/(dx*dy)
          nceo = nceo + pi*quadrr(iq)/dx
        else if (islctype(iz) == 2) then
          nc = nc + 4*pi*(quadrr(iq)**2 - (quadrr(iq)-2.*max(dx,dy))**2)/(dx*dy)
          nceo = nceo + 2*pi*quadrr(iq)/dx
        else if (islctype(iz) == 4 .or. islctype(iz) == 5) then
          nc = nc + (1+nx)*(1+ny) - pi*quadpa(iq)**2/(dx*dy)
          nceo = nceo + pi*quadpa(iq)/dx/2
        endif
        if (iz > -1) then
          if ((islctype(iz) == 4 .and. islctype(iz-1) /= 4) .or.
     &        (islctype(iz) == 5 .and. islctype(iz-1) /= 5) .or.
     &        (islctype(iz-1) == 4 .and. islctype(iz) /= 4) .or.
     &        (islctype(iz-1) == 5 .and. islctype(iz) /= 5)) then
            nceo = nceo + (1+nx)*(1+ny) - pi*quadpa(iq)**2/(dx*dy)
          endif
        endif

c     --- end loop over z which fills islctype
      enddo

c     --- Change the memory allocated based on the estimates made above.
      ncondmax = max(ncondmax,ncond + nc)
      ncndmax = max(ncndmax,max(necndbdy,nocndbdy) + nceo)
      call gchange("PSOR3d",0)

c Set slicv in planes at ends of quads to allow for sub-grid placement of
c boundaries. This also takes care of the boundary planes.
      do iz=0,nz
        if (islctype(iz)==0 .and. islctype(iz-1).ne.0) then
          slicvx(iz)=slicvx(iz-1)
          slicvy(iz)=slicvy(iz-1)
        elseif (islctype(iz)==0 .and. islctype(iz+1).ne.0) then
          slicvx(iz)=slicvx(iz+1)
          slicvy(iz)=slicvy(iz+1)
        endif
      enddo

c     --- Print notice that both hard edged and self-consistent
c     --- models are being used for the same quadrupole, if
c     --- this is the case. This is done outside the loop so it is only done
c     --- once.
      if (lnotice) then
        call remark("Notice: both hard edged and self-consistent models are")
        call remark("        being used for the same quad.")
      endif

c +++++++++++++++++++++++++++
c For each rod, find points that make up the circle of the rods and points
c that are near the rod surfaces.
c The subroutine rodpoint is used for all four of the rods.
c For subgrid scale placement of boundaries,
c names with 'e' are even in x, y, and z, names with 'o' are odd.
c With no symmetries, points on and near all four rods are saved.

c     --- find the location of the center of the grid (x=0 and y=0)
      ix_axis = - xmmin/dx + fuzz
      iy_axis = - ymmin/dy + fuzz

c     --- set max so rod loops cover one quarter of grid
      xmax0 = nx - ix_axis
      ymax0 = ny - iy_axis

c     --- loop over z planes
      iq = 0
      do iz=0,nz
c       --- find location in lattice array
        zz = zgrid + zmmin + iz*dz

c       --- Get id of quadrupole element
        call getquadid(zz,offset,iq,1)

c       --- set up temp variables for locating quadrupoles
        qdcenter = 0.5*(quadzs(iq) + quadze(iq)) + offset

c       --- calculate (in lab frame) start and end of quad
        qs = qdcenter - 0.5*(quadrl(iq) + quadgl(iq))
        qe = qdcenter + 0.5*(quadrl(iq) + quadgl(iq))

c       --- set temps ...
        izp1 = iz + 1
        izm1 = iz - 1
c       --- ... to get the boundaries correct
        if (iz == 0  .and. bound0  == 1) izm1 = 1
        if (iz == 0  .and. bound0  == 2) izm1 = nz-1
        if (iz == nz .and. boundnz == 1) izp1 = nz-1
        if (iz == nz .and. boundnz == 2) izp1 = 1

c       --- Check for presence of rods or ends of rods on x axis.
c       --- Must be careful with boundaries; if either axial boundary is a
c       --- conducting wall, don't put conductor there.
        if (((islctype(iz) == 1 .or. islctype(iz) == 2) .or.
     &      ((islctype(iz) == 3 .or. islctype(iz) == 0) .and.
     &       (islctype(izm1) == 2 .or. islctype(izp1) == 2)))) then

c         --- set logical, true when z plane inside rod, false when just
c         --- outside the end of the rod
          if (islctype(iz) == 1 .or. islctype(iz) == 2) then
            lz_in_rod = .true.
          else
            lz_in_rod = .false.
          endif

c         +++ do rod on positive x axis
c         --- rod radius squared, will be negative if quadrr is negative
c         --- signifying a hyperbolic electrode.
          rodrrsq = sign((quadrr(iq) + qdelrxp(iq) - (zz-qs)*quadsl(iq))**2,
     &                   quadrr(iq))

c         --- rod radius plus one grid cell squared
          rodmxsq=(quadrr(iq) + qdelrxp(iq) - (zz-qs)*quadsl(iq)+max(dx,dy))**2

c         --- inner radius of hollow cylinder (rod is cylinder with walls
c         --- 2 grid cells thick)
          rminsq =(quadrr(iq)+qdelrxp(iq)-(zz-qs)*quadsl(iq)-2.*max(dx,dy))**2

c         --- if at end of rod, fill in all of the rod
          if (islctype(izm1) == 3 .or. islctype(izp1) == 3 .or.
     &        islctype(izm1) == 0 .or. islctype(izp1) == 0) rminsq = 0.

c         --- Set temp which can be set to zero so all of rod is included.
c         --- rrr is used as minimum radius to check for points near surface
          rrr = rodrrsq

c         --- if plane near end of rod, set rrr=0 to get all points in disc
c         --- and set delz since same for whole plane
          delz = 2.
          if (islctype(iz) == 3 .or. islctype(iz) == 0) then
            rrr = 0.
            if (islctype(izm1) == 2) delz = (zz-qs-quadrl(iq)+qdelglx(iq))/dz
            if (islctype(izp1) == 1) delz = (zz-qs)/dz
            if (islctype(izp1) == 2) delz = (zz-qe+quadrl(iq)-qdelglx(iq))/dz
            if (islctype(izm1) == 1) delz = (zz-qe)/dz
          endif

c         --- find rectangle which circumscribes the circle of the rod
c         --- if hyperbolic electrode go out to wall
          if (rodrrsq >= 0.) then
            ixmax = min(xmax0,(rodfract*2.*(quadrr(iq)+qdelrxp(iq)) +
     &                 quadap(iq) + qdelaxp(iq) + qoffx(iq))/dx + 1.)
            ixmin = max(-xmax0,(quadap(iq) + qdelaxp(iq) + qoffx(iq) +
     &                        (zz - qs)*quadsl(iq))/dx - 1.)
            iymax = min(ymax0,(quadrr(iq) + qdelrxp(iq) + qdeloxp(iq) -
     &                 (zz - qs)*quadsl(iq) + qoffy(iq))/dy + 1.)
            iymin = max(-ymax0,(- quadrr(iq) - qdelrxp(iq) +
     &                 qdeloxp(iq) + (zz - qs)*quadsl(iq) +
     &                  qoffy(iq))/dy - 1.)
          else
            ixmax = xmax0
            ixmin = max(-xmax0,(quadap(iq) + qdelaxp(iq) + qoffx(iq))/dx - 1.)
            iymax = min(ymax0, sqrt ((xmax0*dx)**2 - (quadap(iq) +
     &               qoffx(iq) + qdelaxp(iq))**2)/dy + 1 )
            iymin = max(-ymax0,-sqrt ((xmax0*dx)**2 - (quadap(iq) +
     &               qoffx(iq) + qdelaxp(iq))**2)/dy - 1 )
          endif
          if (l2symtry .or. l4symtry) iymin = 0

c         --- distance of rod center from centerline or distance from axis of
c         --- hyperbola intercept
          if (rodrrsq >= 0. ) then
            rodx = quadap(iq) +qdelaxp(iq) +quadrr(iq) + qdelrxp(iq)+qoffx(iq)
            rody = qoffy(iq) + qdeloxp(iq)
            rodap = quadap(iq) + qdelaxp(iq) + qoffx(iq) + (zz-qs)*quadsl(iq)
          else
            rodx = quadap(iq) +qdelaxp(iq) + qoffx(iq)
            rody = 0.
          endif

          call rodpoint(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
     &               lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz,
     &               quadsl(iq),rodap,fuzz,qs-zgrid-zmmin,dz,iz,
     &               slicvx(iz)+qdelvxp(iq))

c         +++ do rod on negative x axis
          if (.not. l4symtry) then
c           --- rod radius squared
            rodrrsq = sign((quadrr(iq) + qdelrxp(iq) - (zz-qs)*quadsl(iq))**2,
     &                     quadrr(iq))

c           --- rod radius plus one grid cell squared
            rodmxsq=(quadrr(iq) + qdelrxm(iq)-(zz-qs)*quadsl(iq)+max(dx,dy))**2

c           --- inner radius of hollow cylinder (rod is cylinder with walls
c           --- 2 grid cells thick)
            rminsq =(quadrr(iq)+qdelrxm(iq)-(zz-qs)*quadsl(iq)-2.*max(dx,dy))**2

c           --- if at end of rod, fill in all of the rod
            if (islctype(izm1) == 3 .or. islctype(izp1) == 3 .or.
     &          islctype(izm1) == 0 .or. islctype(izp1) == 0) rminsq = 0.

c           --- Set temp which can be set to zero so all of rod is included.
c           --- rrr is used as minimum radius to check for points near surface
            rrr = rodrrsq

c           --- if plane near end of rod, set rrr=0 to get all points in disc
c           --- and set delz since same for whole plane
            delz = 2.
            if (islctype(iz) == 3 .or. islctype(iz) == 0) then
              rrr = 0.
              if (islctype(izm1) == 2) delz=(zz-qs-quadrl(iq)+qdelglx(iq))/dz
              if (islctype(izp1) == 1) delz=(zz-qs)/dz
              if (islctype(izp1) == 2) delz=(zz-qe+quadrl(iq)-qdelglx(iq))/dz
              if (islctype(izm1) == 1) delz=(zz-qe)/dz
            endif

c           --- find rectangle which circumscribes the circle of the rod
c           --- if hyperbolic electrode go out to wall
            if (rodrrsq >= 0.) then
              ixmax = min(xmax0,(-quadap(iq) - qdelaxm(iq) + qoffx(iq) -
     &                          (zz - qs)*quadsl(iq))/dx + 1.)
              ixmin = max(-xmax0,(-rodfract*2.*(quadrr(iq)+qdelrxm(iq)) -
     &                    quadap(iq) + qoffx(iq))/dx - 1.)
              iymax = min(ymax0,(quadrr(iq) + qdelrxm(iq) + qdeloxm(iq) -
     &                    (zz - qs)*quadsl(iq) + qoffy(iq))/dy + 1.)
              iymin = max(-ymax0,(- quadrr(iq) - qdelrxm(iq) + qdeloxm(iq) +
     &                    (zz - qs)*quadsl(iq) + qoffy(iq))/dy - 1.)
            else
              ixmax = min(xmax0,(-quadap(iq) - qdelaxm(iq) + qoffx(iq))/dx+1.)
              ixmin = -xmax0
              iymax = min(ymax0, sqrt ((xmax0*dx)**2 - (quadap(iq) +
     &                 qoffx(iq) + qdelaxm(iq))**2)/dy + 1 )
              iymin = max(-ymax0, -sqrt ((xmax0*dx)**2 - (-quadap(iq) +
     &               qoffx(iq) - qdelaxm(iq))**2)/dy - 1 )
            endif
            if (l2symtry) iymin = 0

c           --- distance of rod center from centerline
            if (rodrrsq >= 0. ) then
              rodx = -quadap(iq)-qdelaxm(iq)-quadrr(iq)-qdelrxm(iq)+qoffx(iq)
              rody = qoffy(iq) + qdeloxm(iq)
              rodap = -quadap(iq) -qdelaxm(iq) +qoffx(iq) - (zz-qs)*quadsl(iq)
            else
              rodx = -quadap(iq) -qdelaxm(iq) + qoffx(iq)
              rody = 0.
            endif

            call rodpoint(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
     &               lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz,
     &               quadsl(iq),rodap,fuzz,qs-zgrid-zmmin,dz,iz,
     &               slicvx(iz)+qdelvxm(iq))

          endif
        endif


c       --- Check for presence of rods or ends of rods on y axis.
c       --- Must be careful with boundaries; if either axial boundary is a
c       --- conducting wall, don't put conductor there.
        if (((islctype(iz) == 2 .or. islctype(iz) == 3) .or.
     &      ((islctype(iz) == 1 .or. islctype(iz) == 0) .and.
     &       (islctype(izm1) == 2 .or. islctype(izp1) == 2)))) then

c         --- set logical, true when z plane inside rod, false when just
c         --- outside the end of the rod
          if (islctype(iz) == 2 .or. islctype(iz) == 3) then
            lz_in_rod = .true.
          else
            lz_in_rod = .false.
          endif

c         +++ do rod on positive y axis
c         --- rod radius squared
          rodrrsq = sign((quadrr(iq) + qdelrxp(iq) - (zz-qs)*quadsl(iq))**2,
     &                   quadrr(iq))

c         --- rod radius plus one grid cell squared
          rodmxsq=(quadrr(iq)+qdelryp(iq)-(zz-qs)*quadsl(iq)+max(dx,dy))**2

c         --- inner radius of hollow cylinder (rod is cylinder with walls
c         --- 2 grid cells thick)
          rminsq =(quadrr(iq)+qdelryp(iq)-(zz-qs)*quadsl(iq)-2.*max(dx,dy))**2

c         --- if at end of rod, fill in all of the rod
          if (islctype(izm1) == 1 .or. islctype(izp1) == 1 .or.
     &        islctype(izm1) == 0 .or. islctype(izp1) == 0) rminsq = 0.

c         --- Set temp which can be set to zero so all of rod is included.
c         --- rrr is used as minimum radius to check for points near surface
          rrr = rodrrsq

c         --- if plane near end of rod, set rrr=0 to get all points in disc
c         --- and set delz since same for whole plane
          delz = 2.
          if (islctype(iz) == 1 .or. islctype(iz) == 0) then
            rrr = 0.
            if (islctype(izm1) == 2) delz = (zz-qs-quadrl(iq)+qdelgly(iq))/dz
            if (islctype(izp1) == 3) delz = (zz-qs)/dz
            if (islctype(izp1) == 2) delz = (zz-qe+quadrl(iq)-qdelgly(iq))/dz
            if (islctype(izm1) == 3) delz = (zz-qe)/dz
          endif

c         --- find rectangle which circumscribes the circle of the rod
          if (rodrrsq >= 0.) then
           ixmax = min(xmax0,(quadrr(iq) + qdelryp(iq) - (zz - qs)*quadsl(iq) +
     &                       qoffx(iq) + qdeloyp(iq))/dx + 1.)
           ixmin = min(-xmax0,(-quadrr(iq) - qdelryp(iq) +
     &                (zz - qs)*quadsl(iq) + qoffx(iq) + qdeloyp(iq))/dx - 1.)
           iymax = min(ymax0,(rodfract*2.*(quadrr(iq)+qdelryp(iq)) +quadap(iq)+
     &                       qdelayp(iq) + qoffy(iq))/dy + 1.)
           iymin = max(-ymax0,(quadap(iq) + qdelayp(iq) + qoffy(iq) +
     &                          (zz - qs)*quadsl(iq))/dy - 1.)
          else
           ixmax = min(xmax0, sqrt((ymax0*dy)**2 - (quadap(iq) +
     &                 qoffy(iq) + qdelayp(iq))**2)/dx + 1.)
           ixmin = max(-xmax0, -sqrt((ymax0*dy)**2 - (quadap(iq) +
     &                 qoffy(iq) + qdelayp(iq))**2)/dx - 1.)
           iymax = ymax0
           iymin = max(-ymax0,(quadap(iq) + qdelayp(iq) + qoffy(iq))/dy - 1.)
          endif
          if (l4symtry) ixmin = 0

c         --- distance of rod center from centerline
          if (rodrrsq >= 0.) then
            rodx = qoffx(iq) + qdeloyp(iq)
            rody = quadap(iq) + qdelayp(iq) + quadrr(iq) +qdelryp(iq)+qoffy(iq)
            rodap = quadap(iq) + qdelayp(iq) + qoffy(iq) + (zz-qs)*quadsl(iq)
          else
            rodx = 0.
            rody = quadap(iq) + qdelayp(iq) + qoffy(iq)
          endif

          call rodpoint(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
     &               lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz,
     &               quadsl(iq),rodap,fuzz,qs-zgrid-zmmin,dz,iz,
     &               slicvy(iz)+qdelvyp(iq))

c         +++ do rod on negative y axis
          if (.not. (l2symtry .or. l4symtry)) then
c           --- rod radius squared
            rodrrsq = sign((quadrr(iq) + qdelrxp(iq) - (zz-qs)*quadsl(iq))**2,
     &                     quadrr(iq))

c           --- rod radius plus one grid cell squared
            rodmxsq=(quadrr(iq)+qdelrym(iq)-(zz-qs)*quadsl(iq)+max(dx,dy))**2

c           --- inner radius of hollow cylinder (rod is cylinder with walls
c           --- 2 grid cells thick)
            rminsq =(quadrr(iq)+qdelrym(iq)-(zz-qs)*quadsl(iq)-2.*max(dx,dy))**2

c           --- if at end of rod, fill in all of the rod
            if (islctype(izm1) == 1 .or. islctype(izp1) == 1 .or.
     &          islctype(izm1) == 0 .or. islctype(izp1) == 0) rminsq = 0.

c           --- Set temp which can be set to zero so all of rod is included.
c           --- rrr is used as minimum radius to check for points near surface
            rrr = rodrrsq

c           --- if plane near end of rod, set rrr=0 to get all points in disc
c           --- and set delz since same for whole plane
            delz = 2.
            if (islctype(iz) == 1 .or. islctype(iz) == 0) then
              rrr = 0.
              if (islctype(izm1) == 2) delz=(zz-qs-quadrl(iq)+qdelgly(iq))/dz
              if (islctype(izp1) == 3) delz=(zz-qs)/dz
              if (islctype(izp1) == 2) delz=(zz-qe+quadrl(iq)-qdelgly(iq))/dz
              if (islctype(izm1) == 3) delz=(zz-qe)/dz
            endif

c           --- find rectangle which circumscribes the circle of the rod
            if (rodrrsq >= 0.) then
              iymax = min(ymax0,(-quadap(iq) - qdelaym(iq) + qoffy(iq) -
     &                         (zz - qs)*quadsl(iq))/dy + 1.)
              iymin = max(-ymax0,(-rodfract*2.*(quadrr(iq)+qdelrym(iq)) -
     &                   quadap(iq) - qdelaym(iq) + qoffy(iq))/dy - 1.)
              ixmax = min(xmax0,(quadrr(iq) + qdelrym(iq) - (zz-qs)*quadsl(iq)+
     &                          qoffx(iq) + qdeloym(iq))/dx + 1.)
              ixmin = min(-xmax0,(-quadrr(iq) - qdelrym(iq) +
     &                 (zz - qs)*quadsl(iq) + qoffx(iq) + qdeloym(iq))/dx - 1.)
            else
              ixmax = min(xmax0, sqrt((ymax0*dy)**2 - (-quadap(iq) -
     &                   qdelaym(iq) + qoffy(iq))**2)/dx + 1.)
              ixmin = max(-xmax0, -sqrt((ymax0*dy)**2 - (-quadap(iq) -
     &                   qdelaym(iq) + qoffy(iq))**2)/dx - 1.) 
              iymax = min(ymax0, ( -quadap(iq) - qdelaym(iq) +
     &                    qoffy(iq))/dy + 1.)
              iymin = -ymax0
            endif 

c           --- distance of rod center from centerline
            if (rodrrsq >= 0. ) then
              rodx = qoffx(iq) + qdeloym(iq)
              rody = -quadap(iq)-qdelaym(iq)-quadrr(iq)-qdelrym(iq)+qoffy(iq)
              rodap = -quadap(iq) -qdelaym(iq) +qoffy(iq) - (zz-qs)*quadsl(iq)
            else
              rodx = 0.
              rody = -quadap(iq) -qdelaym(iq) +qoffy(iq)
            endif

            call rodpoint(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
     &               lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz,
     &               quadsl(iq),rodap,fuzz,qs-zgrid-zmmin,dz,iz,
     &               slicvy(iz)+qdelvym(iq))

          endif
        endif

c +++++++++++++++++++++++++++
c       --- calculates the points on the endplates including slices before
c       --- and after for purposes of subgrid placement.  Note that the
c       --- subgrid placement around the abutting rod electrode is
c       --- inconsistent with thte calculated interpolation prameters from
c       --- the rod so that an anamolous "weld-like" structure is created.
c       --- Since this structure is at the same potential as the rod and
c       --- the plate it is neglected in this coding.
 

c       --- plate on left hand side of quad
        if (islctype(iz) == 4 .and. islctype(iz-1) /= 4) then
          if (quadgp(iq) < 0.) then
            vvv = slicvy(iz) + qdelvyp(iq)
          else
            vvv = slicvx(iz) + qdelvxp(iq)
          endif

          lsrlinr = .true.
          npnts_sr = 2
          call gchange("Surface_of_Rev",0)
          z_sr(1) = qs-quadpw(iq)-qdelpwl(iq)
          z_sr(2) = qs
          r_sr = quadpa(iq) + qdelpal(iq)
          rmax = (quadpr(iq) + qdelprl(iq))
          xcent = qoffx(iq) + ix_axis*dx + xmmin
          ycent = qoffy(iq) + iy_axis*dy + ymmin
          call srfrvout3d(" ",vvv,z_sr(1),z_sr(2),qoffx(iq),qoffy(iq),rmax,
     &                    .false.,xmmin,xmmin+nx*dx,ymmin,ymmin+ny*dy,.true.,
     &                    zmmin,zmmin+nz*dz,zgrid,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,0,0,l2symtry,l4symtry)
          lsrlinr = .false.

        endif

c       --- plate on right hand side of quad
        if (islctype(iz) == 5 .and. islctype(iz-1) /= 5) then
          if (quadgp(iq) > 0.) then
            vvv = slicvy(iz) + qdelvyp(iq)
          else
            vvv = slicvx(iz) + qdelvxp(iq)
          endif

          lsrlinr = .true.
          npnts_sr = 2
          call gchange("Surface_of_Rev",0)
          z_sr(1) = qe
          z_sr(2) = qe+quadpw(iq)+qdelpwl(iq)
          r_sr = quadpa(iq) + qdelpar(iq)
          rmax = quadpr(iq) + qdelprr(iq)
          call srfrvout3d(" ",vvv,z_sr(1),z_sr(2),qoffx(iq),qoffy(iq),rmax,
     &                    .false.,xmmin,xmmin+nx*dx,ymmin,ymmin+ny*dy,.true.,
     &                    zmmin,zmmin+nz*dz,zgrid,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,0,0,l2symtry,l4symtry)
          lsrlinr = .false.

        endif

c +++++++++++++++++++++++++++
c       --- Make arrays bigger if they are about to be filled up.
c       --- If ncond is within (1+nx)*(1+ny) of ncondmax, increase ncondmax.
c       --- That is the maximum number of points needed to fill an endplate.
c       --- If the max of necndbdy and nocndbdy are within 3*max(nx,ny) if
c       --- ncndmax, increase ncndmax.  That is roughly the number of points
c       --- needed for one axial slice of a rod.
c       --- Set flag ichange so that if anything needs to be changed, all
c       --- of the changes can be made at once.
        ichange = 0
        if (ncond + (1+nx)*(1+ny) > ncondmax) then
          ncondmax = ncond + (1+nx)*(1+ny)
          ichange = 1
        endif
        if (max(necndbdy,nocndbdy) + 3*max(nx,ny) > ncndmax) then
          ncndmax = max(necndbdy,nocndbdy) + 3*max(nx,ny)
          ichange = 1
        endif
        if (ichange == 1) call gchange("PSOR3d",0)

c     --- end loop over z
      enddo

c +++++++++++++++++++++++++++
c     --- add points for electrostatic dipoles
      if (ndipo > 0) then

c       --- set prevcent (center of previous dipole) to large negative
c       --- value to ensure that first dipole is included.
        prevcent = -LARGEPOS

c       --- loop over z in the lattice array looking for dipoles
        id = 0
        do iz=0,nz
          zz = zgrid + zmmin + iz*dz
          call getdipoid(zz,offset,id,1)

c         --- calculate center of dipole
          dipocent = .5*(dipozs(id) + dipoze(id)) - zmmin - zgrid

c         --- If the center is greater than center of previous dipole,
c         --- then at a new dipole, so load in points.
          if (dipocent > prevcent) then
            prevcent = dipocent + zfuzz

c           --- calculate x position of first plate and its distance from the
c           --- grid cell below
            ix_near = dipox1(id)/dx + ix_axis + 0.5
            ix = dipox1(id)/dx + ix_axis
            delx = dipox1(id)/dx+ix_axis - int(dipox1(id)/dx+ix_axis)

c           --- check if plate is within field grid
            if (0 <= ix .and. ix <= nx) then 

c             --- calculate bounds of first plate (these are bounded to be
c             --- within the field mesh)
              iymax = min(ny,int( 0.5*dipow1(id)/dy+0.5+iy_axis))
              iymin = max(0,int(-0.5*dipow1(id)/dy-0.5+iy_axis))
              izmax = min(nz,int((dipocent + dipol1(id)*0.5)/dz + zfuzz))
              izmin = max(0,int((dipocent - dipol1(id)*0.5)/dz + zfuzz))

c             --- add more space to the arrays if necessary
              if (ncond + (iymax-iymin+1)*(izmax-izmin+1) > ncondmax) then
                ncondmax = ncond + (iymax-iymin+1)*(izmax-izmin+1)
                call gchange("PSOR3d",0)
              endif

              call dipopnt(iymin,iymax,izmin,izmax,dipov1(id),ix_near,ix,delx)

            endif

c           --- now do the same for the second dipole plate
            ix_near = dipox2(id)/dx + ix_axis + 0.5
            ix = dipox2(id)/dx + ix_axis
            delx = dipox2(id)/dx+ix_axis - int(dipox2(id)/dx+ix_axis)
            if (0 <= ix .and. ix <= nx) then 
              iymax = min(ny,int( 0.5*dipow2(id)/dy+0.5+iy_axis))
              iymin = max(0,int(-0.5*dipow2(id)/dy-0.5+iy_axis))
              izmax = min(nz,int((dipocent + dipol2(id)*0.5)/dz + 0.5))
              izmin = max(0,int((dipocent - dipol2(id)*0.5)/dz + 0.5))
              if (ncond + (iymax-iymin+1)*(izmax-izmin+1) > ncondmax) then
                ncondmax = ncond + (iymax-iymin+1)*(izmax-izmin+1)
                call gchange("PSOR3d",0)
              endif
              call dipopnt(iymin,iymax,izmin,izmax,dipov2(id),ix_near,ix,delx)
            endif

c         --- end of section that fills arrays
          endif

c       --- end of loop over z
        enddo

c     --- end of dipoles
      endif

c +++++++++++++++++++++++++++
c Include injection source
      call injctpnt(dx,dy,dz,ix_axis,iy_axis,nx,ny,zmmin)

      return
      end
c===========================================================================
      subroutine rodpoint(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
     &      lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz_in,
     &      quadsl,rodap,fuzz,qzs,dz,iz,voltage)
      use PSOR3d
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      integer(ISZ):: ixmin,ixmax,iymin,iymax,ix_axis,iy_axis,iz
      real(kind=8):: dx,dy,rodx,rody,rminsq,rodrrsq,rrr,rodmxsq,delz_in,quadsl
      real(kind=8):: rodap,fuzz,qzs,dz,voltage
      logical(ISZ):: lz_in_rod

c Subroutine that calculates conductor points for rods.  It actually only
c calculates the points for a circle on the plane iz.  To produce a rod, call
c with each value of iz within rod, and also the two values just outside the
c rod if sub-gridscale boundary placements are being used.
c Routine can be used for any general rod.
c Stores points near all rods, including points within dz of ends of rods.
c Temporarily stores points in even arrays before checking if even or odd.
c None of the input variables are affected, only the arrays in PSOR3d are
c changed.
c Note that if rodrrsq is negative, hyperbolic shaped electrodes, going out
c to the endge of the system is calculated instead of rods.
c INPUT:
c   ixmin,ixmax,iymin,iymax: rectangle which surrounds the circle (coordinates
c                            relative to zero on axis)
c   dx,dy:                   grid cell sizes
c   rodx,rody:               distance of rod center from axis
c   lz_in_rod:               logical, true if iz within rod
c   rminsq:                  inner radius of hollow cylinder squared
c   rodrrsq:                 rod radius squared
c   ix_axis,iy_axis:         grid cell indices of beam center
c   rrr:                     minimum rod radius squared: when iz within rod,
c                            same as rodrrsq, when iz outside rod, is zero
c   rodmxsq:                 rod radius plus one grid cell, all squared
c   delz_in:                 when iz is outside rod, distance from end of rod
c                            when iz is inside, must be > 1.
c   quadsl:                  rod slant (for conical rod only)
c   rodap:                   pole tip distance, (for conical rod only)
c   fuzz:                    number smalled compared to grid cell size, used
c                            to prevent precision problems
c   qzs:                     location of rod left end with respect to grid
c                            (for conical rod only)
c   dz:                      grid cell size (for conical rod only)
c   iz:                      grid cell number
c   voltage:                 voltage on rod
c
c WARNING: The coding for the conical rod most likely does not work.  It was
c          put in but never fully tested.


      integer(ISZ):: ix,iy,npoints,nsgpoints
      real(kind=8):: radsq,delx,dely,delz,rodxsq,rodysq,dxi,dyi,dzi,quadsli
      real(kind=8):: x2fuzz,y2fuzz,r2fuzz

c      if(solvergeom==RZgeom) then
c        call rodpointrz(ixmin,ixmax,iymin,iymax,dx,dy,rodx,rody,
c     &      lz_in_rod,rminsq,rodrrsq,ix_axis,iy_axis,rrr,rodmxsq,delz_in,
c     &      quadsl,rodap,fuzz,qzs,dz,iz,voltage)
c        return
c      end if
c     --- Make sure there is enough space. npoints is the maximum
c     --- number of points that could be needed. nsgpoints is 
c     --- big enough to always be larger than the number of subgrid
c     --- points.
      npoints = (ixmax-ixmin+1)*(iymax-iymin+1)
      nsgpoints = 3*(ixmax-ixmin + iymax-iymin)
      if (lz_in_rod) then
        nsgpoints = 3*(ixmax-ixmin + iymax-iymin)
      else
        nsgpoints = npoints
      endif
      if (ncond+npoints > ncondmax .or.
     &    necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax) then
        if (ncond+npoints > ncondmax) ncondmax = ncond + npoints
        if (necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax)
     &    ncndmax = max(necndbdy,nocndbdy) + nsgpoints
          call gchange("PSOR3d",0)
      endif

c     --- set temporaries
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      quadsli = 1./dvnz(quadsl)
      rodxsq = rodx**2
      rodysq = rody**2
      x2fuzz = dx*dx*fuzz
      y2fuzz = dy*dy*fuzz
      r2fuzz = min(x2fuzz,y2fuzz)

c     --- loop over rectangle
      do iy = iymin,iymax
        do ix = ixmin,ixmax

          if (rodrrsq >= 0.) then

c           --- If point in wall of cylinder, or, when whole cylinder is not
c           --- in grid, if point is on back side, and if arrays not filled
c           --- yet, select point.
            radsq = (ix*dx - rodx)**2 + (iy*dy - rody)**2
            if (lz_in_rod .and.
     &        ((rminsq <= radsq+r2fuzz .and. radsq-r2fuzz <= rodrrsq) .or.
     &        ((ix==ixmax .or. iy==iymax .or. ix==ixmin .or. iy==iymin)
     &         .and. radsq-r2fuzz <= rodrrsq)) .and. ncond < ncondmax) then
              ncond = ncond + 1
              ixcond(ncond) = ix + ix_axis
              iycond(ncond) = iy + iy_axis
              izcond(ncond) = iz
              condvolt(ncond) = voltage
            endif

          else

c           --- hyperbolic electrodes
            if (rody == 0.) then
c              --- horizontal electrode
               if (lz_in_rod .and. ((iy*dy)**2-y2fuzz < (ix*dx)**2-rodxsq)
     &            .and. ncond < ncondmax) then
                 ncond = ncond + 1
                 ixcond(ncond) = ix + ix_axis
                 iycond(ncond) = iy + iy_axis
                 izcond(ncond) = iz
                 condvolt(ncond) = voltage
               endif
            else
c              --- vertical electrode
               if (lz_in_rod .and. ((ix*dx)**2-x2fuzz < (iy*dy)**2-rodysq)
     &             .and. ncond < ncondmax) then
                 ncond = ncond + 1
                 ixcond(ncond) = ix + ix_axis
                 iycond(ncond) = iy + iy_axis
                 izcond(ncond) = iz
                 condvolt(ncond) = voltage
               endif
            endif

          endif

c         --- For circular rods, if point is within one grid cell outside
c         --- of the conductor and arrays not full, select point.   This
c         --- is only done to avoid extra calculation.  For hyperbolic
c         --- rods (rodrrsq < 0), since it is less trivial to estimate
c         --- the distance between the point and the surface, the exact
c         --- distance is calculated for all points and only the check
c         --- on whether the arrays are full is used.
          if (((rrr < radsq .and. radsq < rodmxsq) .or.
     &         rodrrsq < 0.) .and.
     &        necndbdy < ncndmax .and. nocndbdy < ncndmax) then

c           --- set to default values
            delx = 2.
            dely = 2.
            delz = 2.

c           --- find transverse distances depending on electrode shape
            if (rodrrsq >= 0.) then
c             --- do circular rod
c             --- find distance in x from boundary
              if (rodrrsq-(iy*dy-rody)**2 > 0. .and. lz_in_rod) then
                delx = ix - rodx*dxi -
     &            sign(sqrt(rodrrsq-(iy*dy-rody)**2)*dxi,ix - rodx*dxi)
              endif

c             --- find distance in y from boundary
              if (rodrrsq-(ix*dx-rodx)**2 > 0. .and. lz_in_rod) then
                dely = iy - rody*dyi -
     &            sign(sqrt(rodrrsq-(ix*dx-rodx)**2)*dyi,iy - rody*dyi)
              endif

            else
c             --- do hyperbolic electrodes

c             --- horizontal electrode
              if (rody==0) then
c               --- test that point is outside electrode
                if ((iy*dy)**2 > (ix*dx)**2 - rodxsq .and. lz_in_rod) then
                  delx = ix - sign(sqrt((iy*dy)**2 + rodxsq)*dxi, 1.*ix)
                  if ((ix*dx)**2 > rodxsq)
     &              dely = iy - sign(sqrt((ix*dx)**2 - rodxsq)*dyi, 1.*iy)
                endif
              endif

c             --- vertical electrode
              if (rodx==0) then
c               --- test that point is outside electrode
                if ((ix*dx)**2 > (iy*dy)**2 - rodysq .and. lz_in_rod) then
                  dely = iy - sign(sqrt((ix*dx)**2 + rodysq)*dyi, 1.*iy)
                  if ((iy*dy)**2 > rodysq)
     &              delx = ix - sign(sqrt((iy*dy)**2 - rodysq)*dxi, 1.*ix)
                endif
              endif

            endif
 
c           --- find distance in z from boundary
            if (quadsl .ne. 0.) delz = (qzs+(ix*dx-rodap)*quadsli)*dzi-iz
            if (.not. lz_in_rod) then

c             --- If within radius of circular rod, then set
              if (rodrrsq > 0.) then
                if (rminsq <= radsq+r2fuzz .and. radsq-r2fuzz <= rodrrsq) then 
                  delz = delz_in
                endif

c             --- If within hyperbola, then set
              else
                if (rody == 0.) then
c                 --- horizontal electrode
                  if ((iy*dy)**2-y2fuzz <= (ix*dx)**2-rodxsq) delz = delz_in
                else
c                 --- vertical electrode
                  if ((ix*dx)**2-x2fuzz <= (iy*dy)**2-rodysq) delz = delz_in
                endif
              endif
            endif

c           --- if within one cell in any direction, save the point
            if (fuzz < abs(delx) .and. abs(delx) < 1.-fuzz .or.
     &          fuzz < abs(dely) .and. abs(dely) < 1.-fuzz .or.
     &          fuzz < abs(delz) .and. abs(delz) < 1.-fuzz) then
              if (mod(ix_axis + ix + iy_axis + iy + iz + zparity,2) == 0) then
c               --- even
                necndbdy = necndbdy + 1
                ecdelmx(necndbdy) = 2.
                ecdelpx(necndbdy) = 2.
                ecdelmy(necndbdy) = 2.
                ecdelpy(necndbdy) = 2.
                ecdelmz(necndbdy) = 2.
                ecdelpz(necndbdy) = 2.
                if (fuzz <  delx .and.  delx < 1.-fuzz) ecdelmx(necndbdy)= delx
                if (fuzz <  dely .and.  dely < 1.-fuzz) ecdelmy(necndbdy)= dely
                if (fuzz <  delz .and.  delz < 1.-fuzz) ecdelmz(necndbdy)= delz
                if (fuzz < -delx .and. -delx < 1.-fuzz) ecdelpx(necndbdy)=-delx
                if (fuzz < -dely .and. -dely < 1.-fuzz) ecdelpy(necndbdy)=-dely
                if (fuzz < -delz .and. -delz < 1.-fuzz) ecdelpz(necndbdy)=-delz
                iecndx(necndbdy) = ix + ix_axis
                iecndy(necndbdy) = iy + iy_axis
                iecndz(necndbdy) = iz
                ecvolt(necndbdy) = voltage
              else
c               --- odd
                nocndbdy = nocndbdy + 1
                ocdelmx(nocndbdy) = 2.
                ocdelpx(nocndbdy) = 2.
                ocdelmy(nocndbdy) = 2.
                ocdelpy(nocndbdy) = 2.
                ocdelmz(nocndbdy) = 2.
                ocdelpz(nocndbdy) = 2.
                if (fuzz <  delx.and. delx < 1.-fuzz) ocdelmx(nocndbdy) =  delx
                if (fuzz <  dely.and. dely < 1.-fuzz) ocdelmy(nocndbdy) =  dely
                if (fuzz <  delz.and. delz < 1.-fuzz) ocdelmz(nocndbdy) =  delz
                if (fuzz < -delx.and.-delx < 1.-fuzz) ocdelpx(nocndbdy) = -delx
                if (fuzz < -dely.and.-dely < 1.-fuzz) ocdelpy(nocndbdy) = -dely
                if (fuzz < -delz.and.-delz < 1.-fuzz) ocdelpz(nocndbdy) = -delz
                iocndx(nocndbdy) = ix + ix_axis
                iocndy(nocndbdy) = iy + iy_axis
                iocndz(nocndbdy) = iz
                ocvolt(nocndbdy) = voltage
              endif
            endif
          endif
        enddo
      enddo

c     --- Print warning if conductor arrays filled up.
      if (necndbdy == ncndmax .or. nocndbdy == ncndmax .or.
     &    ncond == ncondmax) then
        call remark("WARNING: The arrays holding the conductor data have filled")
        call remark("         up.  Increase size of arrays to ensure all of")
        call remark("         the conductors have been generated.")
        call remark("         In function rodpoint.")
      endif

      return
      end
c===========================================================================
      subroutine platepnt(ixmin,ixmax,iymin,iymax,ix_axis,iy_axis,dx,dy,
     &                    aper,rmax,vvv,xoff,yoff,delz_in,iz,lz_in_plate,fuzz)
      use PSOR3d
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      integer(ISZ):: ixmin,ixmax,iymin,iymax,ix_axis,iy_axis,iz
      real(kind=8):: dx,dy,aper,rmax,vvv,xoff,yoff,delz_in,fuzz
      logical(ISZ):: lz_in_plate

c Set points which are on a tranverse plate.
c
c     ixmin    minimum value of ix
c     ixmax    maximum value of ix
c     iymin    minimum value of iy
c     iymax    maximum value of iy
c     ix_axis  x grid location of beam center
c     dx       grid cell size in x
c     iy_axis  y grid location of beam center
c     dy       grid cell size in y
c     aper     inner aperture of the plate
c     rmax     maximum radius of the plate
c     vvv      voltage on the plate
c     xoff     x offset of aperture
c     yoff     y offset of aperture
c     delz_in  fraction of cell, when outside of plate, to edge of plate
c     iz       axial grid location of the plate
c     fuzz     number smalled compared to grid cell size, used
c              to prevent precision problems


      integer(ISZ):: ix,iy
      real(kind=8):: rad,pradsq,prmxsq,paprsq,papmsq,delx,dely,delz,dxi,dyi
      real(kind=8):: r2fuzz
      integer(ISZ):: npoints,nsgpoints
      logical(ISZ):: lout

c      if(solvergeom==RZgeom) then
c        call platepntrz(ixmin,ixmax,iymin,iymax,ix_axis,iy_axis,dx,dy,
c     &                    aper,rmax,vvv,xoff,yoff,delz_in,iz,lz_in_plate,fuzz)
c        return
c      end if

      dxi = 1./dx
      dyi = 1./dy
      pradsq = rmax**2
      prmxsq = (abs(rmax) + max(dx,dy))**2
      paprsq = aper**2
      papmsq = (aper - max(dx,dy))**2
      if (aper - max(dx,dy) < 0.) papmsq = -papmsq
      r2fuzz = min(dx,dy)**2*fuzz

c     --- The sign of rmax is used to flag whether or not the include subgrid
c     --- data at the outer edge of the aperture, at rmax.
      if (rmax < 0.) then
        lout = .false.
      else
        lout = .true.
      endif

c     --- Make sure there is enough space. npoints is the maximum
c     --- number of points that could be needed. nsgpoints is 
c     --- big enough to always be larger than the number of subgrid
c     --- points.
      npoints = (ixmax-ixmin+1)*(iymax-iymin+1)
      if (lz_in_plate) then
        nsgpoints = 3*(ixmax-ixmin + iymax-iymin)
      else
        nsgpoints = npoints
      endif
      if (ncond+npoints > ncondmax .or.
     &    necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax) then
        if (ncond+npoints > ncondmax) ncondmax = ncond + npoints
        if (necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax)
     &    ncndmax = max(necndbdy,nocndbdy) + nsgpoints
          call gchange("PSOR3d",0)
      endif

c     --- loop over transverse plane
      do ix = ixmin,ixmax
        do iy = iymin,iymax
          rad = ((ix-ix_axis)*dx-xoff)**2 + ((iy-iy_axis)*dy-yoff)**2
          if (lz_in_plate .and.
     &        paprsq <= rad+r2fuzz .and. rad-r2fuzz <= pradsq .and.
     &        ncond < ncondmax) then
            ncond = ncond + 1
            ixcond(ncond) = ix
            iycond(ncond) = iy
            izcond(ncond) = iz
            condvolt(ncond) = vvv
          endif

c         --- calculate interpolation values
c         --- set to default values
          delx = 2.
          dely = 2.
          delz = 2.

c         --- if point within one grid cell outside of the plate
c         --- and arrays not full, select point.
          if ((pradsq < rad-r2fuzz .and. rad+r2fuzz < prmxsq) .and.
     &        lz_in_plate .and. lout .and.
     &        necndbdy < ncndmax .and.
     &        nocndbdy < ncndmax) then

c           --- find distance in x from boundary
            if (pradsq-((iy-iy_axis)*dy-yoff)**2 > 0.) then
              delx = ix - ix_axis - xoff*dxi -
     &               sign(sqrt(pradsq-((iy-iy_axis)*dy-yoff)**2)*dxi,
     &                    ix - ix_axis - xoff*dxi)
            endif

c           --- find distance in y from boundary
            if (pradsq-((ix-ix_axis)*dx-xoff)**2 > 0.) then
              dely = iy - iy_axis - yoff*dyi -
     &               sign(sqrt(pradsq-((ix-ix_axis)*dx-xoff)**2)*dyi,
     &                    iy - iy_axis - yoff*dyi)
            endif
          endif
 
c         --- if point within one grid cell inside of the plate
c         --- and arrays not full, select point.
          if ((papmsq < rad-r2fuzz .and. rad+r2fuzz < paprsq) .and.
     &        lz_in_plate .and.
     &        necndbdy < ncndmax .and.
     &        nocndbdy < ncndmax) then

c           --- find distance in x from boundary
            if (paprsq-((iy-iy_axis)*dy-yoff)**2 > 0.) then
              delx = ix - ix_axis - xoff*dxi -
     &               sign(sqrt(paprsq-((iy-iy_axis)*dy-yoff)**2)*dxi,
     &                    ix - ix_axis - xoff*dxi)
            endif

c           --- find distance in y from boundary
            if (paprsq-((ix-ix_axis)*dx-xoff)**2 > 0.) then
              dely = iy - iy_axis - yoff*dyi -
     &               sign(sqrt(paprsq-((ix-ix_axis)*dx-xoff)**2)*dyi,
     &                    iy - iy_axis - yoff*dyi)
            endif
          endif

c         --- find distance in z from boundary
          if (.not. lz_in_plate .and.
     &        (paprsq < rad+r2fuzz .and. rad-r2fuzz < pradsq) .and.
     &        necndbdy < ncndmax .and.
     &        nocndbdy < ncndmax) then
            delz = delz_in
          endif

c         --- if within one cell in any direction, save the point
          if (fuzz < abs(delx) .and. abs(delx) < 1.-fuzz .or.
     &        fuzz < abs(dely) .and. abs(dely) < 1.-fuzz .or.
     &        fuzz < abs(delz) .and. abs(delz) < 1.-fuzz) then
            if (mod(ix + iy + iz + zparity, 2) == 0) then
c             --- even
              necndbdy = necndbdy + 1
              ecdelmx(necndbdy) = 2.
              ecdelpx(necndbdy) = 2.
              ecdelmy(necndbdy) = 2.
              ecdelpy(necndbdy) = 2.
              ecdelmz(necndbdy) = 2.
              ecdelpz(necndbdy) = 2.
              if (fuzz <  delx.and. delx < 1.-fuzz) ecdelmx(necndbdy) =  delx
              if (fuzz <  dely.and. dely < 1.-fuzz) ecdelmy(necndbdy) =  dely
              if (fuzz <  delz.and. delz < 1.-fuzz) ecdelmz(necndbdy) =  delz
              if (fuzz < -delx.and.-delx < 1.-fuzz) ecdelpx(necndbdy) = -delx
              if (fuzz < -dely.and.-dely < 1.-fuzz) ecdelpy(necndbdy) = -dely
              if (fuzz < -delz.and.-delz < 1.-fuzz) ecdelpz(necndbdy) = -delz
              iecndx(necndbdy) = ix
              iecndy(necndbdy) = iy
              iecndz(necndbdy) = iz
              ecvolt(necndbdy) = vvv
            else
c             --- odd
              nocndbdy = nocndbdy + 1
              ocdelmx(nocndbdy) = 2.
              ocdelpx(nocndbdy) = 2.
              ocdelmy(nocndbdy) = 2.
              ocdelpy(nocndbdy) = 2.
              ocdelmz(nocndbdy) = 2.
              ocdelpz(nocndbdy) = 2.
              if (fuzz <  delx.and. delx < 1.-fuzz) ocdelmx(nocndbdy) =  delx
              if (fuzz <  dely.and. dely < 1.-fuzz) ocdelmy(nocndbdy) =  dely
              if (fuzz <  delz.and. delz < 1.-fuzz) ocdelmz(nocndbdy) =  delz
              if (fuzz < -delx.and.-delx < 1.-fuzz) ocdelpx(nocndbdy) = -delx
              if (fuzz < -dely.and.-dely < 1.-fuzz) ocdelpy(nocndbdy) = -dely
              if (fuzz < -delz.and.-delz < 1.-fuzz) ocdelpz(nocndbdy) = -delz
              iocndx(nocndbdy) = ix
              iocndy(nocndbdy) = iy
              iocndz(nocndbdy) = iz
              ocvolt(nocndbdy) = vvv
            endif
          endif
        enddo
      enddo

c     --- Print warning if conductor arrays filled up.
      if (necndbdy == ncndmax .or. nocndbdy == ncndmax .or.
     &    ncond == ncondmax) then
        call remark("WARNING: The arrays holding the conductor data have filled")
        call remark("         up.  Increase size of arrays to ensure all of")
        call remark("         the conductors have been generated.")
        call remark("         In function platepnt.")
      endif

      return
      end
c===========================================================================
      subroutine dipopnt(iymin,iymax,izmin,izmax,vvv,ix_near,ix,delx)
      use PSOR3d
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      integer(ISZ):: iymin,iymax,izmin,izmax,ix_near,ix
      real(kind=8):: vvv,delx

c Set points which are on a dipole plate.  This needs to be in a subroutine
c since the calling routine changes the PSOR3d arrays.
c
c     iymin  minimum value of iy
c     iymax  maximum value of iy
c     izmin  minimum value of iz
c     izmax  maximum value of iz
c     vvv     voltage on the dipole plate
c     ix      radial grid location of the plate
c
c If not using sub-grid boundaries, that places conductor points, otherwise
c saves data for sub-grid boundaries.


      integer(ISZ):: iy,iz,npoints,nsgpoints

c      if(solvergeom==RZgeom) then
c        call dipopntrz(iymin,iymax,izmin,izmax,vvv,ix_near,ix,delx)
c        return
c      end if

c     --- Make sure there is enough space. npoints is the maximum
c     --- number of points that could be needed. nsgpoints is 
c     --- big enough to always be larger than the number of subgrid
c     --- points.
      npoints = (izmax-izmin+1)*(iymax-iymin+1)
      nsgpoints = npoints
      if (ncond+npoints > ncondmax .or.
     &    necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax) then
        if (ncond+npoints > ncondmax) ncondmax = ncond + npoints
        if (necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax)
     &    ncndmax = max(necndbdy,nocndbdy) + nsgpoints
          call gchange("PSOR3d",0)
      endif

c     --- loop over plate to fill in data
      do iy=iymin,iymax
        do iz=izmin,izmax
          if (.not. lcndbndy .or. delx == 0) then
            ncond = ncond + 1
            ixcond(ncond) = ix_near
            iycond(ncond) = iy
            izcond(ncond) = iz
            condvolt(ncond) = vvv
          else

            if (mod(ix + iy + iz + zparity,2) == 0) then
c             --- ix is even
              necndbdy = necndbdy + 1
              ecdelmx(necndbdy) = 2.
              ecdelpx(necndbdy) = delx
              ecdelmy(necndbdy) = 2.
              ecdelpy(necndbdy) = 2.
              ecdelmz(necndbdy) = 2.
              ecdelpz(necndbdy) = 2.
              iecndx(necndbdy) = ix
              iecndy(necndbdy) = iy
              iecndz(necndbdy) = iz
              ecvolt(necndbdy) = vvv
c             --- ix+1 is odd
              nocndbdy = nocndbdy + 1
              ocdelmx(nocndbdy) = 1. - delx
              ocdelpx(nocndbdy) = 2.
              ocdelmy(nocndbdy) = 2.
              ocdelpy(nocndbdy) = 2.
              ocdelmz(nocndbdy) = 2.
              ocdelpz(nocndbdy) = 2.
              iocndx(nocndbdy) = ix + 1
              iocndy(nocndbdy) = iy
              iocndz(nocndbdy) = iz
              ocvolt(nocndbdy) = vvv
            else
c             --- ix is odd
              nocndbdy = nocndbdy + 1
              ocdelmx(nocndbdy) = 2.
              ocdelpx(nocndbdy) = delx
              ocdelmy(nocndbdy) = 2.
              ocdelpy(nocndbdy) = 2.
              ocdelmz(nocndbdy) = 2.
              ocdelpz(nocndbdy) = 2.
              iocndx(nocndbdy) = ix
              iocndy(nocndbdy) = iy
              iocndz(nocndbdy) = iz
              ocvolt(nocndbdy) = vvv
c             --- ix+1 is even
              necndbdy = necndbdy + 1
              ecdelmx(necndbdy) = 1. - delx
              ecdelpx(necndbdy) = 2.
              ecdelmy(necndbdy) = 2.
              ecdelpy(necndbdy) = 2.
              ecdelmz(necndbdy) = 2.
              ecdelpz(necndbdy) = 2.
              iecndx(necndbdy) = ix + 1
              iecndy(necndbdy) = iy
              iecndz(necndbdy) = iz
              ecvolt(necndbdy) = vvv
            endif
          endif
        enddo
      enddo

c     --- Print warning if conductor arrays filled up.
      if (necndbdy == ncndmax .or. nocndbdy == ncndmax .or.
     &    ncond == ncondmax) then
        call remark("WARNING: The arrays holding the conductor data have filled")
        call remark("         up.  Increase size of arrays to ensure all of")
        call remark("         the conductors have been generated.")
        call remark("         In function dipopnt.")
      endif

      return
      end
c===========================================================================
      subroutine injctpnt(dx,dy,dz,ix_axis,iy_axis,nx,ny,zmmin)
      use Constant
      use PSOR3d
      use InjectVars
      use InjectVars3d
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      use multigridRZ
      real(kind=8):: dx,dy,dz,zmmin
      integer(ISZ):: ix_axis,iy_axis,nx,ny

c Sets points on injection source.

      integer(ISZ):: ij,inj_ij,ix,iy,iz,npoints,nsgpoints
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      real(kind=8):: xx,yy,zz,rmaxsq,ao_dx,bo_dy,ai_dx,bi_dy
      real(kind=8):: delmx,delpx,delmy,delpy,delmz,delpz,rrsq
      real(kind=8):: fuzz

c     --- Just return if injection is turned off or the source is not to be
c     --- setup.
      if (inject == 0 .or. .not. lvinject) return

c     --- Set size of fuzz
      fuzz = 1.e-5

c     --- calculate points for injection source
      do ij=1,ninject
        inj_ij = min(ij,inj_ninj)
c       --- find bounds of loop
        ixmin = max(1,       int((-ainject(ij)+xinject(ij))/dx+ix_axis-1))
        ixmax = min(inj_nx-1,int((+ainject(ij)+xinject(ij))/dx+ix_axis+1))
        iymin = max(1,       int((-binject(ij)+yinject(ij))/dy+iy_axis-1))
        iymax = min(inj_ny-1,int((+binject(ij)+yinject(ij))/dy+iy_axis+1))

c       --- Define maximum and minimum transverse extents of source
        ao_dx = ainject(ij) + dx
        bo_dy = binject(ij) + dy
        ai_dx = ainjmin(ij)
        bi_dy = binjmin(ij)

c       --- Make sure there is enough space. npoints is the maximum
c       --- number of points that could be needed. nsgpoints is 
c       --- big enough to always be larger than the number of subgrid
c       --- points.
        npoints = pi*(ao_dx*bo_dy - ai_dx*bi_dy)/(dx*dy) + 1
        nsgpoints = npoints
        if (ncond+npoints > ncondmax .or.
     &      necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax) then
          if (ncond+npoints > ncondmax) ncondmax = ncond + npoints
          if (necndbdy+nsgpoints > ncndmax .or. nocndbdy+nsgpoints > ncndmax)
     &      ncndmax = max(necndbdy,nocndbdy) + nsgpoints
          call gchange("PSOR3d",0)
        endif

c       --- Loop over transverse points
        do ix = ixmin,ixmax
          do iy = iymin,iymax
            xx = (ix-ix_axis)*dx - xinject(ij)
            yy = (iy-iy_axis)*dy - yinject(ij)
c           --- skip points not on the surface
            if (((xx*bo_dy)**2 + (yy*ao_dx)**2 > (ao_dx*bo_dy)**2) .or.
     &          ((xx*bi_dy)**2 + (yy*ai_dx)**2 < (ai_dx*bi_dy)**2) .or.
     &          ((xx**2 + yy**2) > rmaxsq)) cycle

c           --- If rinject is not set by the user, rely on inj_grid (surface
c           --- will be piecewise linear), otherwise make the surface spherical.
            if (rinject(ij) == LARGEPOS) then
              iz = int((inj_grid(ix,iy,inj_ij) - zmmin)/dz + 1)
              zz = zmmin + iz*dz
              delmx = 2.
              delpx = 2.
              delmy = 2.
              delpy = 2.
              delmz = iz - (inj_grid(ix,iy,inj_ij) - zmmin)/dz
              delpz = 2.
              if (inj_grid(ix-1,iy,inj_ij) > zz)
     &          delmx = delmz/(inj_grid(ix-1,iy,inj_ij)-inj_grid(ix,iy,inj_ij))*dz
              if (inj_grid(ix+1,iy,inj_ij) > zz)
     &          delpx = delmz/(inj_grid(ix+1,iy,inj_ij)-inj_grid(ix,iy,inj_ij))*dz
              if (inj_grid(ix,iy-1,inj_ij) > zz)
     &          delmy = delmz/(inj_grid(ix,iy-1,inj_ij)-inj_grid(ix,iy,inj_ij))*dz
              if (inj_grid(ix,iy+1,inj_ij) > zz)
     &          delpy = delmz/(inj_grid(ix,iy+1,inj_ij)-inj_grid(ix,iy,inj_ij))*dz
            else
              zz = (zinject(ij) - zmmin) + (xx**2 + yy**2)/
     &             (rinject(ij) + sqrt(rinject(ij)**2 - (xx**2 + yy**2)))
              iz = int(zz/dz + 1)
              delmx = 2.
              delpx = 2.
              delmy = 2.
              delpy = 2.
              delmz = iz - zz/dz
              delpz = 2.
              if (rinject(ij) > 0.) then
                rrsq = 2.*rinject(ij)*zz - zz**2
                delpx = (sqrt(rrsq - yy**2) - xx)/dx
                delpy = (sqrt(rrsq - xx**2) - yy)/dy
              endif
            endif
c           --- Now add conductor point or subgrid point
            if (delmx < fuzz .or. delpx < fuzz .or.
     &          delmy < fuzz .or. delpy < fuzz .or.
     &          delmz < fuzz) then
c             --- If any del's are small, then all del's are small, so add
c             --- a conductor point.
              ncond = ncond + 1
              ixcond(ncond) = ix
              iycond(ncond) = iy
              izcond(ncond) = iz
              condvolt(ncond) = vinject(ij)
            else if (delmx < 1.-fuzz .or. delpx < 1.-fuzz .or.
     &               delmy < 1.-fuzz .or. delpy < 1.-fuzz .or.
     &               delmz < 1.-fuzz) then
c             --- Add subgrid point
              if (mod(ix+iy+iz,2) == 0) then
c               --- even
                necndbdy = necndbdy + 1
                iecndx(necndbdy) = ix
                iecndy(necndbdy) = iy
                iecndz(necndbdy) = iz
                ecdelmx(necndbdy) = delmx
                ecdelmy(necndbdy) = delmy
                ecdelmz(necndbdy) = delmz
                ecdelpx(necndbdy) = delpx
                ecdelpy(necndbdy) = delpy
                ecdelpz(necndbdy) = delpz
                ecvolt(necndbdy) = vinject(ij)
              else
c               --- odd
                nocndbdy = nocndbdy + 1
                iocndx(nocndbdy) = ix
                iocndy(nocndbdy) = iy
                iocndz(nocndbdy) = iz
                ocdelmx(nocndbdy) = delmx
                ocdelmy(nocndbdy) = delmy
                ocdelmz(nocndbdy) = delmz
                ocdelpx(nocndbdy) = delpx
                ocdelpy(nocndbdy) = delpy
                ocdelpz(nocndbdy) = delpz
                ocvolt(nocndbdy) = vinject(ij)
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
c===========================================================================
c===========================================================================
      subroutine srfrvout(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                    xmin,xmax,ymin,ymax,lshell,
     &                    zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      use Surface_of_Rev
      use InGen, ONLY:fstype
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the outside of a surface of revolution.
c The Basis routine rofzfunc should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvout("rofz",...)
c
c Input:
c   rofzfunc  is the name of the basis routine which calculates the radius as a
c             function of z.
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   rmax      is the maximum radius of the conductor.
c   lfill     is the logical requesting that the whole conductor be filled
c             with points.
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.

      real(kind=8):: zz
      integer(ISZ):: i
      real(kind=8):: srfrv_f

c     --- Table-ize the data if srfrv_pernz > 0.
      if (srfrv_pernz > 0 .and. .not. lsrlinr) then
        npnts_sr = (zmax - zmin)/dz*srfrv_pernz + 1
        zz = (zmax - zmin)/(npnts_sr - 1)
        call gchange("Surface_of_Rev",0)
        do i=1,npnts_sr
          z_sr(i) = zmin + (i-1)*zz
          r_sr(i) = srfrv_f(z_sr(i),rofzfunc,0)
        enddo
        lsrlinr = .true.
      endif

c     --- Call the appropriate version of the routine.
      if(solvergeom==RZgeom) then
        call srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                    xmin,xmax,ymin,ymax,lshell,
     &                    zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      else if (fstype == 7 .or. fstype == 11) then
        call srfrvout3dmg(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                    xmin,xmax,ymin,ymax,lshell,
     &                    zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      else
        call srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                  xmin,xmax,ymin,ymax,lshell,
     &                  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                  ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      end if

c     --- Clean up the table-ized data.
      if (srfrv_pernz > 0) then
        lsrlinr = .false.
      endif


      return
      end
c===========================================================================
      subroutine srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                      xmin,xmax,ymin,ymax,lshell,
     &                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the outside of a surface of revolution.
c The Basis routine rofzfunc should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvout("rofz",...)
c
c Input:
c   rofzfunc  is the name of the basis routine which calculates the radius as a
c             function of z.
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   rmax      is the maximum radius of the conductor.
c   lfill     is the logical requesting that the whole conductor be filled
c             with points.
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.

      integer(ISZ):: ix,iy,iz,i
      real(kind=8):: rrr,rad,rrrsq,rrrsq1,rrrmax,rrrmin
      real(kind=8):: rmaxsq,zgmin,zgmax,r_zmin,r_zmax
      real(kind=8):: rrrlower,rrrhigher,zlower,zhigher
      real(kind=8):: xxmin,xxmax,yymin,yymax
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: jxmin,jxmax,jymin,jymax
      integer(ISZ):: izbmin,izbmax,iterrf
      real(kind=8):: fuzz,xfuzz,yfuzz,rfuzz,zfuzz,xx,yy,zz
      real(kind=8):: delmz0,delpz0
      real(kind=8):: delmx,delmy,delmz,delpx,delpy,delpz
      real(kind=8):: max_dxdy,dxi,dyi,dzi
      real(kind=8):: ddz,znew,rnew
      real(kind=8):: srfrv_f
      real(kind=8):: x(2),f(2)

      rmaxsq = rmax*rmax
      max_dxdy = max(dx,dy)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c Set fuzz which is used to avoid problems with roundoff.
      fuzz = 1.e-10
      xfuzz = dx*fuzz
      yfuzz = dy*fuzz
      rfuzz = min(xfuzz,yfuzz)
      zfuzz = dz*fuzz

c Find maximal transverse extent of the conductors with an additional grid
c cell added to allow subgrid-scale resolution at the edge.
      jxmin = int(max(xmin*dxi,(-rmax+xcent)*dxi)) + ix_axis - 1
      jymin = int(max(ymin*dyi,(-rmax+ycent)*dyi)) + iy_axis - 1
      jxmax = int(min(xmax*dxi,(+rmax+xcent)*dxi)) + ix_axis + 1
      jymax = int(min(ymax*dyi,(+rmax+ycent)*dyi)) + iy_axis + 1

c Set limits for transverse loops based on the transverse boundary conditions
c and the symmetry.
      ixmin = max(0   ,jxmin)
      iymin = max(0   ,jymin)
      ixmax = min(nx  ,jxmax)
      iymax = min(ny  ,jymax)

c Calculate bounding box relative to center of surface of revolution.
      xxmin = xmin - xcent
      xxmax = xmax - xcent
      yymin = ymin - ycent
      yymax = ymax - ycent

c Set limits for longitudinal loops based on boundary conditions.  If
c using Dirichlet boundaries, don't put conductor at planes iz=0 and nz.
      izbmin = 0
      izbmax = nz

c Set limits for longitudinal loops.  They are set so that the minimum
c end is rounded up and maximum end is rounded down, ensuring that iz is
c within the conductor.
      zgmin = zmmin + zbeam
      zgmax = zmmax + zbeam
      izmin = nz+int((zmin-zgmax)*dzi-zfuzz)
      izmax =    int((zmax-zgmin)*dzi+zfuzz)

c Perform error checking, to make sure that izmin and izmax are within the
c grid.
      izmin = max(izbmin,izmin)
      izmax = min(izbmax,izmax)
     
c If izmin > izbmax or izmax < izbmin then the conductor is completely off
c of the grid, so return.
      if (izmin > izbmax .or. izmax < izbmin) then
        call remark("WARNING:srfrvout: conductor is completely off the grid.")
        return
      endif

c Set conductor points.  If lfill or at end, set all points outside radius.
      do iz=izmin,izmax
        zz = iz*dz+zgmin
        zlower = max(zz-dz,zmin)
        zhigher = min(zz+dz,zmax)
c       --- Make sure there is enough room in the arrays.
        if (ncond + (1+nx)*(1+ny) > ncondmax) then
          ncondmax = ncond + (1+nx)*(1+ny)
          call gchange("PSOR3d",0)
        endif
c       --- Calculate the surface radius at iz.
        rrr = srfrv_f(zz,rofzfunc,0)
        rrrsq = rrr*rrr
c       --- Calculate the maximum radius to be filled.
c       --- If not filling the whole conductor, set rrrmax large enough so
c       --- that enough conductor is filled to match neighboring points in iz.
        if (lfill .or. (lshell .and. (iz == izmin .or. iz == izmax))) then
          rrrmax = rmax
        else
          rrrmax = max((rrr+max_dxdy),
     &              (srfrv_f(zlower,rofzfunc,0)+max_dxdy),
     &              (srfrv_f(zhigher,rofzfunc,0)+max_dxdy))
        endif
c       --- Loop over transverse plane.
        do ix=ixmin,ixmax
          do iy=iymin,iymax
            xx = (ix-ix_axis)*dx-xcent
            yy = (iy-iy_axis)*dy-ycent
            rad = sqrt(xx**2 + yy**2)

c           --- The points selected must meet a number of conditions.
c           --- > The radius must be between the surface radius and rmax.
c           --- > The position must be within the min/max transverse extent
c           --- > of the conductor.
c           --- The points selected must also meet one of the following.
c           --- > A radius less than rrrmax.
c           --- > A radius within one grid cell of rmax.
c           --- > A position within a grid cell of the transverse bounds.
c           --- The last two conditions put an optional shell around an
c           --- unfilled conductor.
            if ((rrr-rfuzz <= rad .and. rad <= rmax+rfuzz .and.
     &           xxmin <= xx .and. xx <= xxmax .and.
     &           yymin <= yy .and. yy <= yymax) .and.
     &          (rad <= rrrmax .or.
     &           (lshell .and. (rmax-rad < max_dxdy .or.
     &                          xx-xxmin < max_dxdy .or.
     &                          yy-yymin < max_dxdy .or.
     &                          xxmax-xx < max_dxdy .or.
     &                          yymax-yy < max_dxdy      )))) then
              ncond = ncond + 1
              ixcond(ncond) = ix
              iycond(ncond) = iy
              izcond(ncond) = iz
              condvolt(ncond) = volt
            endif
          enddo
        enddo
      enddo
c End of conductor points section

c Set subgrid scale arrays
      if (lcndbndy) then

c       --- Calculate r at zmin and zmax.
        r_zmin = srfrv_f(zmin,rofzfunc,0)
        r_zmax = srfrv_f(zmax,rofzfunc,0)

c       --- Recalculate limits for longitudinal loops.  Now, they are rounded
c       --- to include the plane just outside the conductor.
        izmin = max(izbmin,int((zmin-zgmin)*dzi-zfuzz))
        izmax = min(izbmax,nz+int((zmax-zgmax)*dzi+zfuzz))

c       -- Set all points along and just outside conductor.
        do iz=izmin,izmax
          zz = iz*dz+zgmin
          zlower = max(zz-dz,zmin)
          zhigher = min(zz+dz,zmax)

          if (necndbdy + (1+nx)*(1+ny) > ncndmax .or.
     &        nocndbdy + (1+nx)*(1+ny) > ncndmax) then
            ncndmax = max(necndbdy,nocndbdy) + (1+nx)*(1+ny)
            call gchange("PSOR3d",0)
          endif

c         --- Calculate radius of the surface.
c         --- When zz is outside the conductor, get radius at end of
c         --- conductor.  rrrsq is used to chose which grid points may
c         --- be within the subgrid-scale region.  rrrsq1 is used in the
c         --- calculation of the transverse distances.  When zz is
c         --- outside the conductor, only points within the radius rmax
c         --- (when lshell is true) or r_zmin or r_zmax (otherwise) are
c         --- desired, but the radius used in the calculation is set to
c         --- be large so that transverse subgrid resolution will not be done.
          if (zz < zmin-zfuzz) then
            rrr = r_zmin
            if (lshell) then
              rrrsq = rmax*rmax
            else
              rrrsq = rrr*rrr
            endif
            rrrsq1 = LARGEPOS
          elseif (zz > zmax+zfuzz) then
            rrr = r_zmax
            if (lshell) then
              rrrsq = rmax*rmax
            else
              rrrsq = rrr*rrr
            endif
            rrrsq1 = LARGEPOS
          else
            rrr = srfrv_f(zz,rofzfunc,0)
            rrrsq = rrr*rrr
            rrrsq1 = rrrsq
          endif

c         --- Find radius one grid cell lower in z.  delmz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz > zmax, the delmz0 is set to the distance
c         --- between the z plane and the end of the conductor.
c         --- When zz < zmin, rrrlower is set to rrr so the
c         --- iteration is not done.
          rrrlower = srfrv_f(zlower,rofzfunc,0)
          delmz0 = 2.
          if (zz < zmin) then
            rrrlower = rrr
          elseif (zz > zmax) then
            delmz0 = (zz - zmax)*dzi
          endif

c         --- Find radius one grid cell higher in z.  delpz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz < zmin, the delmp0 is set to the distance
c         --- between the z plane and the end of the conductor.
c         --- When zz > zmax, rrrhigher is set to rrr so the
c         --- iteration is not done.
          rrrhigher = srfrv_f(zhigher,rofzfunc,0)
          delpz0 = 2.
          if (zz > zmax) then
            rrrhigher = rrr
          elseif (zz < zmin) then
            delpz0 = (zmin - zz)*dzi
          endif

c         --- Find minimum of rrr and radius of surface at nearby z points.
          rrrmin = min(rrr-max_dxdy,rrrlower,rrrhigher)

          do ix=ixmin,ixmax
            do iy=iymin,iymax
              xx = (ix-ix_axis)*dx-xcent
              yy = (iy-iy_axis)*dy-ycent
              rad = sqrt(xx**2 + yy**2)

c             --- Calculate distances in z first.  Select points which
c             --- are transversely near the surface of revolution and
c             --- are within the bounding box.

c             --- Calculate distance to conductor at lower z.
              delmz = delmz0

c             --- If grid location is between rrr and rrrlower and between
c             --- rmax and rrrlower, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrlower <= rad .and. rad < min(rrr,rmax)-rfuzz
     &            .and. delmz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = min(zz,zmax)
                f(1) = rrr - rad
                x(2) = zlower
                f(2) = rrrlower - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = min((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmax)
                  f(i) = srfrv_f(x(i),rofzfunc,0) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zlower <= x(i) .and. x(i) < zz) then
                  delmz = (zz - x(i))*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rofzfunc,0)
                    if (rnew > rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew <= rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz = (zz - znew)*dzi
                endif
              elseif (rad < rrrlower .or. rad > rmax) then
                delmz = 2.
              endif

c             --- Calculate distance to conductor at higher z.
              delpz = delpz0

c             --- If grid location is between rrr and rrrhigher and between
c             --- rmax and rrrhigher, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrhigher <= rad .and. rad < min(rrr,rmax)-rfuzz
     &            .and. delpz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = max(zz,zmin)
                f(1) = rrr - rad
                x(2) = zhigher
                f(2) = rrrhigher - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = max((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmin)
                  f(i) = srfrv_f(x(i),rofzfunc,0) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zz < x(i) .and. x(i) <= zhigher) then
                  delpz = (x(i) - zz)*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rofzfunc,0)
                    if (rnew <= rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew > rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz = (znew - zz)*dzi
                endif
              elseif (rad < rrrhigher .or. rad > rmax) then
                delpz = 2.
              endif

c             --- Now calculate transverse distances.
              delmx = 2.
              delpx = 2.
              delmy = 2.
              delpy = 2.

c             --- Calculate distances to surface of revolution.
              if (rrrmin+rfuzz < rad .and. rad <= min(rrr,rmax)+rfuzz) then
c               --- Note that for z slices outside the conductor, rrrsq1 is a
c               --- large number, so all transverse distance all also large.
                delmx = (sqrt(max(0.,rrrsq1 - yy**2)) + xx)*dxi
                delpx = (sqrt(max(0.,rrrsq1 - yy**2)) - xx)*dxi
                delmy = (sqrt(max(0.,rrrsq1 - xx**2)) + yy)*dyi
                delpy = (sqrt(max(0.,rrrsq1 - xx**2)) - yy)*dyi
              endif

c             --- If requested, select points on shell around rmax and
c             --- bounding box.
              if (lshell .and. zmin-zfuzz <= zz .and. zz <= zmax+zfuzz) then
c               --- Select grid points within 1 grid cell of rmax.
                if (rmax+rfuzz < rad .and. rad < rmax+max_dxdy-rfuzz) then
                  if (yy**2 < rmaxsq) then
                    delmx = ( xx - sqrt(max(0.,rmaxsq - yy**2)))*dxi
                    delpx = (-xx - sqrt(max(0.,rmaxsq - yy**2)))*dxi
                  endif
                  if (xx**2 < rmaxsq) then
                    delmy = ( yy - sqrt(max(0.,rmaxsq - xx**2)))*dyi
                    delpy = (-yy - sqrt(max(0.,rmaxsq - xx**2)))*dyi
                  endif
                endif
c               --- Select points within 1 grid cell outside of bounding box
c               --- which are also within rrr and rmax.
                if (rrr-rfuzz <= rad .and. rad <= rmax+rfuzz) then
                  if (xx < xxmin) delpx = (xxmin-xx)*dxi
                  if (yy < yymin) delpy = (yymin-yy)*dyi
                  if (xx > xxmax) delmx = (xx-xxmax)*dxi
                  if (yy > yymax) delmy = (yy-yymax)*dyi
                endif
              endif

c             --- Now, ensure that the point on the surface from which the
c             --- distance was calculated is within the bounding box.
              if (xx-delmx*dx < xxmin) delmx = 2.
              if (xx+delpx*dx > xxmax) delpx = 2.
              if (yy-delmy*dy < yymin) delmy = 2.
              if (yy+delpy*dy > yymax) delpy = 2.

c             --- Make sure that distances are greater than zero.
              if (delmx < fuzz) delmx = 2.
              if (delpx < fuzz) delpx = 2.
              if (delmy < fuzz) delmy = 2.
              if (delpy < fuzz) delpy = 2.
              if (delmz < fuzz) delmz = 2.
              if (delpz < fuzz) delpz = 2.

c             -- Check if any of the distances are less than 1.
              if (delmx < 1.-fuzz .or. delpx < 1.-fuzz .or.
     &            delmy < 1.-fuzz .or. delpy < 1.-fuzz .or.
     &            delmz < 1.-fuzz .or. delpz < 1.-fuzz) then
                if (mod(ix+iy+iz+zparity,2) == 0) then
c                 --- even
                  necndbdy = necndbdy + 1
                  iecndx(necndbdy) = ix
                  iecndy(necndbdy) = iy
                  iecndz(necndbdy) = iz
                  ecdelmx(necndbdy) = delmx
                  ecdelmy(necndbdy) = delmy
                  ecdelmz(necndbdy) = delmz
                  ecdelpx(necndbdy) = delpx
                  ecdelpy(necndbdy) = delpy
                  ecdelpz(necndbdy) = delpz
                  ecvolt(necndbdy) = volt
                else
c                 --- odd
                  nocndbdy = nocndbdy + 1
                  iocndx(nocndbdy) = ix
                  iocndy(nocndbdy) = iy
                  iocndz(nocndbdy) = iz
                  ocdelmx(nocndbdy) = delmx
                  ocdelmy(nocndbdy) = delmy
                  ocdelmz(nocndbdy) = delmz
                  ocdelpx(nocndbdy) = delpx
                  ocdelpy(nocndbdy) = delpy
                  ocdelpz(nocndbdy) = delpz
                  ocvolt(nocndbdy) = volt
                endif
              endif
            enddo
          enddo
        enddo
      endif

      return
      end
c========================================================================
c========================================================================
      subroutine srfrvin(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                   xmin,xmax,ymin,ymax,lshell,
     &                   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                   ix_axis,iy_axis,
     &                   xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      use Surface_of_Rev
      use InGen, ONLY:fstype
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmin
      real(kind=8):: xmin,xmax,ymin,ymax
      logical(ISZ):: lfill,lshell
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the inside of a surface of revolution.
c The Basis routine rofzfunc should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvin("rofz",...)
c
c Input:
c   rofzfunc  is the name of the basis routine which calculates the radius as a
c             function of z.
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   rmin      is the minimum radius of the conductor.
c   lfill     is the logical requesting that the whole conductor be filled
c             with points.
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.
c See comments in srfrvout.

      real(kind=8):: zz,srfrv_f
      integer(ISZ):: i

c     --- Table-ize the data if srfrv_pernz > 0.
      if (srfrv_pernz > 0 .and. .not. lsrlinr) then
        npnts_sr = (zmax - zmin)/dz*srfrv_pernz + 1
        zz = (zmax - zmin)/(npnts_sr - 1)
        call gchange("Surface_of_Rev",0)
        do i=1,npnts_sr
          z_sr(i) = zmin + (i-1)*zz
          r_sr(i) = srfrv_f(z_sr(i),rofzfunc,0)
        enddo
        lsrlinr = .true.
      endif

      if (fstype == 7 .or. fstype == 11) then
        call srfrvin3dmg(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                 xmin,xmax,ymin,ymax,lshell,
     &                 zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                 ix_axis,iy_axis,
     &                 xmesh,ymesh,l2symtry,l4symtry)
      else
c     --- Call the appropriate routine
        call srfrvin3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                 xmin,xmax,ymin,ymax,lshell,
     &                 zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                 ix_axis,iy_axis,
     &                 xmesh,ymesh,l2symtry,l4symtry)
      endif

c     --- Clean up the table-ized data.
      if (srfrv_pernz > 0) then
        lsrlinr = .false.
      endif

      return
      end
c========================================================================
      subroutine srfrvin3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                     xmin,xmax,ymin,ymax,lshell,
     &                     zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                     ix_axis,iy_axis,
     &                     xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmin
      real(kind=8):: xmin,xmax,ymin,ymax
      logical(ISZ):: lfill,lshell
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the inside of a surface of revolution.
c The Basis routine rofzfunc should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvin("rofz",...)
c
c Input:
c   rofzfunc  is the name of the basis routine which calculates the radius as a
c             function of z.
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   rmin      is the minimum radius of the conductor.
c   lfill     is the logical requesting that the whole conductor be filled
c             with points.
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.
c See comments in srfrvout.


      integer(ISZ):: ix,iy,iz,i
      real(kind=8):: rrr,rad,rrrsq,rrrsq1,rrrmax,rrrmin
      real(kind=8):: rminsq,zgmin,zgmax,r_zmin,r_zmax
      real(kind=8):: rrrlower,rrrhigher,zlower,zhigher
      real(kind=8):: xxmin,xxmax,yymin,yymax
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: jxmin,jxmax,jymin,jymax
      integer(ISZ):: izbmin,izbmax,iterrf
      real(kind=8):: fuzz,xfuzz,yfuzz,rfuzz,zfuzz,xx,yy,zz
      real(kind=8):: delmz0,delpz0
      real(kind=8):: delmx,delmy,delmz,delpx,delpy,delpz
      real(kind=8):: max_dxdy,dxi,dyi,dzi
      real(kind=8):: ddz,znew,rnew
      real(kind=8):: srfrv_f
      real(kind=8):: x(2),f(2)

      rminsq = rmin*rmin
      max_dxdy = max(dx,dy)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c Set fuzz which is used to avoid problems with roundoff.
      fuzz = 1.e-10
      xfuzz = dx*fuzz
      yfuzz = dy*fuzz
      rfuzz = min(xfuzz,yfuzz)
      zfuzz = dz*fuzz

c Find maximal transverse extent of the conductors with an additional grid
c cell added to allow subgrid-scale resolution at the edge.
      jxmin = int(xmin*dxi) + ix_axis - 1
      jymin = int(ymin*dyi) + iy_axis - 1
      jxmax = int(xmax*dxi) + ix_axis + 1
      jymax = int(ymax*dyi) + iy_axis + 1

c Set limits for transverse loops based on the transverse boundary conditions
c and the symmetry.
      ixmin = max(0 ,jxmin)
      iymin = max(0 ,jymin)
      ixmax = min(nx,jxmax)
      iymax = min(ny,jymax)

c Calculate bounding box relative to center of surface of revolution.
      xxmin = xmin - xcent
      xxmax = xmax - xcent
      yymin = ymin - ycent
      yymax = ymax - ycent

c Set limits for longitudinal loops based on boundary conditions.  If
c using Dirichlet boundaries, don't put conductor at planes iz=0 and nz.
      izbmin = 0
      izbmax = nz

c Set limits for longitudinal loops.  They are set so that the minimum
c end is rounded up and maximum end is rounded down, ensuring that iz is
c within the conductor.
      zgmin = zmmin + zbeam
      zgmax = zmmax + zbeam
      izmin = nz+int((zmin-zgmax)*dzi-zfuzz)
      izmax = int((zmax-zgmin)*dzi+zfuzz)

c Perform error checking, to make sure that izmin and izmax are within the
c grid.
      izmin = max(izbmin,izmin)
      izmax = min(izbmax,izmax)
     
c If izmin > izbmax or izmax < izbmin then the conductor is completely off
c of the grid, so return.
      if (izmin > izbmax .or. izmax < izbmin) then
        call remark("WARNING:srfrvout: conductor is completely off the grid.")
        return
      endif

c Set conductor points.  If lfill or at end, set all points inside radius.
      do iz=izmin,izmax
        zz = iz*dz+zgmin
        zlower = max(zz-dz,zmin)
        zhigher = min(zz+dz,zmax)
c       --- Make sure there is enough room in the arrays.
        if (ncond + (1+nx)*(1+ny) > ncondmax) then
          ncondmax = ncond + (1+nx)*(1+ny)
          call gchange("PSOR3d",0)
        endif
c       --- Calculate the surface radius at iz.
        rrr = srfrv_f(zz,rofzfunc,0)
        if (lfill .or. (lshell .and. (iz == izmin .or. iz == izmax))) then
          rrrmin = rmin
        else
          rrrmin = min((rrr-max_dxdy),
     &              (srfrv_f(zlower,rofzfunc,0)-max_dxdy),
     &              (srfrv_f(zhigher,rofzfunc,0)-max_dxdy))
        endif
c       --- Calculate transverse bounds for this iz plane.
        jxmin = max(ixmin, int((-rrr+xcent)*dxi)+ix_axis)
        jymin = max(iymin, int((-rrr+ycent)*dyi)+iy_axis)
        jxmax = min(ixmax, int(( rrr+xcent)*dxi)+ix_axis)
        jymax = min(iymax, int(( rrr+ycent)*dyi)+iy_axis)
c       --- Loop over transverse plane.  Select all points within the
c       --- range from rrrmin to rrr.
        do ix=jxmin,jxmax
          do iy=jymin,jymax
            xx = (ix-ix_axis)*dx-xcent
            yy = (iy-iy_axis)*dy-ycent
            rad = sqrt(xx**2 + yy**2)

c           --- The points selected must meet a number of conditions.
c           --- > The radius must be between the surface radius and rmin.
c           --- > The position must be within the min/max transverse extent
c           --- > of the conductor.
c           --- The points selected must also meet one of the following.
c           --- > A radius greater than rrrmin.
c           --- > A radius within one grid cell of rmin.
c           --- > A position within a grid cell of the transverse bounds.
c           --- The last two conditions put an optional shell around an
c           --- unfilled conductor.
            if ((rmin-rfuzz <= rad .and. rad <= rrr+rfuzz .and.
     &           xxmin <= xx .and. xx <= xxmax .and.
     &           yymin <= yy .and. yy <= yymax) .and.
     &          (rad >= rrrmin .or.
     &           (lshell .and. (rad-rmin < max_dxdy .or.
     &                          xx-xxmin < max_dxdy .or.
     &                          yy-yymin < max_dxdy .or.
     &                          xxmax-xx < max_dxdy .or.
     &                          yymax-yy < max_dxdy      )))) then
              ncond = ncond + 1
              ixcond(ncond) = ix
              iycond(ncond) = iy
              izcond(ncond) = iz
              condvolt(ncond) = volt
            endif
          enddo
        enddo
      enddo
c End of conductor points section

c Set subgrid scale arrays
      if (lcndbndy) then

c       --- Calculate r at zmin and zmax.
        r_zmin = srfrv_f(zmin,rofzfunc,0)
        r_zmax = srfrv_f(zmax,rofzfunc,0)

c       --- Recalculate limits for longitudinal loops.  Now, they are rounded
c       --- to include the plane just outside the conductor.
        izmin = max(izbmin,int((zmin-zgmin)*dzi-zfuzz))
        izmax = min(izbmax,nz+int((zmax-zgmax)*dzi+zfuzz))

c       -- Loop over z planes
        do iz=izmin,izmax
          zz = iz*dz+zgmin
          zlower = max(zz-dz,zmin)
          zhigher = min(zz+dz,zmax)

          if (necndbdy + (1+nx)*(1+ny) > ncndmax .or.
     &        nocndbdy + (1+nx)*(1+ny) > ncndmax) then
            ncndmax = max(necndbdy,nocndbdy) + (1+nx)*(1+ny)
            call gchange("PSOR3d",0)
          endif

c         --- Calculate radius of the surface.
c         --- When zz is outside the conductor, get radius at end of
c         --- conductor.  rrrsq is used to chose which grid points may
c         --- be within the subgrid-scale region.  rrrsq1 is used in the
c         --- calculation of the transverse distances.  When zz is
c         --- outside the conductor, only points within the radius rmax
c         --- (when lshell is true) or r_zmin or r_zmax (otherwise) are
c         --- desired, but the radius used in the calculation is set to
c         --- be small so that transverse subgrid resolution will not be done.
          if (zz < zmin-zfuzz) then
            rrr = r_zmin
            if (lshell) then
              rrrsq = rmin*rmin
            else
              rrrsq = rrr*rrr
            endif
            rrrsq1 = 0.
          elseif (zz > zmax+zfuzz) then
            rrr = r_zmax
            if (lshell) then
              rrrsq = rmin*rmin
            else
              rrrsq = rrr*rrr
            endif
            rrrsq1 = 0.
          else
            rrr = srfrv_f(zz,rofzfunc,0)
            rrrsq = rrr*rrr
            rrrsq1 = rrrsq
          endif

c         --- Find radius one grid cell lower in z.  delmz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz > zmax, the delmz0 is set to the distance
c         --- between the z plane and the end of the conductor.
c         --- When zz < zmin, rrrlower is set to rrr so the
c         --- iteration is not done.
          rrrlower = srfrv_f(zlower,rofzfunc,0)
          delmz0 = 2.
          if (zz < zmin) then
            rrrlower = rrr
          elseif (zz > zmax) then
            delmz0 = (zz - zmax)*dzi
          endif

c         --- Find radius one grid cell higher in z.  delpz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz < zmin, the delpz0 is set to the distance
c         --- between the z plane and the end of the conductor.
c         --- When zz > zmax, rrrhigher is set to rrr so the
c         --- iteration is not done.
          rrrhigher = srfrv_f(zhigher,rofzfunc,0)
          delpz0 = 2.
          if (zz > zmax) then
            rrrhigher = rrr
          elseif (zz < zmin) then
            delpz0 = (zmin - zz)*dzi
          endif

c         --- Find maximum of rrr and of radius at neighboring z points.
          rrrmax = max(rrr+max_dxdy,rrrlower,rrrhigher)

c         --- Limits of transverse loop.  Covers full radius of current
c         --- z slice and neighboring slices.
          jxmin = max(ixmin, int((-rrrmax+xcent-dx)*dxi)+ix_axis)
          jymin = max(iymin, int((-rrrmax+ycent-dy)*dyi)+iy_axis)
          jxmax = min(ixmax, int(( rrrmax+xcent+dx)*dxi)+ix_axis)
          jymax = min(iymax, int(( rrrmax+ycent+dy)*dyi)+iy_axis)

c         --- Transverse loop
          do ix=jxmin,jxmax
            do iy=jymin,jymax
              xx = (ix-ix_axis)*dx-xcent
              yy = (iy-iy_axis)*dy-ycent
              rad = sqrt(xx**2 + yy**2)

c             --- Calculate distances in z first.  Select points which
c             --- are transversely near the surface of revolution and
c             --- are within the bounding box.

c             --- Calculate distance to conductor at lower z.
              delmz = delmz0

c             --- If grid location is between rrr and rrrlower and between
c             --- rmin and rrrlower, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrlower >= rad .and. rad > max(rrr,rmin)+rfuzz
     &            .and. delmz > zfuzz*dzi .and.
     &            xxmin <= xx .and. xx <= xxmax .and.
     &            yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = min(zz,zmax)
                f(1) = rrr - rad
                x(2) = zlower
                f(2) = rrrlower - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = min((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmax)
                  f(i) = srfrv_f(x(i),rofzfunc,0) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zlower <= x(i) .and. x(i) < zz) then
                  delmz = (zz - x(i))*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rofzfunc,0)
                    if (rnew < rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew >= rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz = (zz - znew)*dzi
                endif
              elseif (rad > rrrlower .or. rad < rmin) then
                delmz = 2.
              endif

c             --- Calculate distance to conductor at higher z.
              delpz = delpz0

c             --- If grid location is between rrr and rrrhigher and between
c             --- rmin and rrrhigher, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrhigher >= rad .and. rad > max(rrr,rmin)+rfuzz
     &            .and. delpz > zfuzz*dzi .and.
     &            xxmin <= xx .and. xx <= xxmax .and.
     &            yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = max(zz,zmin)
                f(1) = rrr - rad
                x(2) = zhigher
                f(2) = rrrhigher - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = max((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmin)
                  f(i) = srfrv_f(x(i),rofzfunc,0) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zz < x(i) .and. x(i) <= zhigher) then
                  delpz = (x(i) - zz)*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rofzfunc,0)
                    if (rnew >= rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew < rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz = (znew - zz)*dzi
                endif
              elseif (rad > rrrhigher .or. rad < rmin) then
                delpz = 2.
              endif

c             --- Now calculate transverse distances.
              delmx = 2.
              delpx = 2.
              delmy = 2.
              delpy = 2.

c             --- Calculate distances to surface of revolution.
              if (max(rrr,rmin)-rfuzz <= rad .and. rad < rrrmax-rfuzz) then
c               --- Note that for z slices outside the conductor, rrrsq1 is
c               --- zero, so distances are set to 2.
                if (yy**2 < rrrsq1) then
                  delmx = ( xx - sqrt(max(0.,rrrsq - yy**2)))*dxi
                  delpx = (-xx - sqrt(max(0.,rrrsq - yy**2)))*dxi
                endif
                if (xx**2 < rrrsq1) then
                  delmy = ( yy - sqrt(max(0.,rrrsq - xx**2)))*dyi
                  delpy = (-yy - sqrt(max(0.,rrrsq - xx**2)))*dyi
                endif
              endif

c             --- If requested, select points on shell around rmin and
c             --- bounding box.
              if (lshell .and. zmin-zfuzz <= zz .and. zz <= zmax+zfuzz) then
c               --- Select points within 1 grid cell of rmin
                if (rmin-max_dxdy+rfuzz < rad .and. rad < rmin-rfuzz) then
                  delmx = (sqrt(max(0.,rminsq - yy**2)) + xx)*dxi
                  delpx = (sqrt(max(0.,rminsq - yy**2)) - xx)*dxi
                  delmy = (sqrt(max(0.,rminsq - xx**2)) + yy)*dyi
                  delpy = (sqrt(max(0.,rminsq - xx**2)) - yy)*dyi
                endif
c               --- Select points within 1 grid cell outside of bounding box
c               --- which are also within rrr and rmin.
                if (rmin-rfuzz <= rad .and. rad <= rrr+rfuzz) then
                  if (xx < xxmin) delpx = (xxmin-xx)*dxi
                  if (yy < yymin) delpy = (yymin-yy)*dyi
                  if (xx > xxmax) delmx = (xx-xxmax)*dxi
                  if (yy > yymax) delmy = (yy-yymax)*dyi
                endif
              endif

c             --- Now, ensure that the point on the surface from which the
c             --- distance was calculated is within the bounding box.
              if (xx-delmx*dx < xxmin) delmx = 2.
              if (xx+delpx*dx > xxmax) delpx = 2.
              if (yy-delmy*dy < yymin) delmy = 2.
              if (yy+delpy*dy > yymax) delpy = 2.

c             --- Make sure that distances are greater than zero.
              if (delmx < fuzz) delmx = 2.
              if (delpx < fuzz) delpx = 2.
              if (delmy < fuzz) delmy = 2.
              if (delpy < fuzz) delpy = 2.
              if (delmz < fuzz) delmz = 2.
              if (delpz < fuzz) delpz = 2.

c             -- Check if any of the distances are less than 1.
              if (delmx < 1.-fuzz .or. delpx < 1.-fuzz .or.
     &            delmy < 1.-fuzz .or. delpy < 1.-fuzz .or.
     &            delmz < 1.-fuzz .or. delpz < 1.-fuzz) then
                if (mod(ix+iy+iz+zparity,2) == 0) then
                  necndbdy = necndbdy + 1
                  iecndx(necndbdy) = ix
                  iecndy(necndbdy) = iy
                  iecndz(necndbdy) = iz
                  ecdelmx(necndbdy) = delmx
                  ecdelmy(necndbdy) = delmy
                  ecdelmz(necndbdy) = delmz
                  ecdelpx(necndbdy) = delpx
                  ecdelpy(necndbdy) = delpy
                  ecdelpz(necndbdy) = delpz
                  ecvolt(necndbdy) = volt
                else
                  nocndbdy = nocndbdy + 1
                  iocndx(nocndbdy) = ix
                  iocndy(nocndbdy) = iy
                  iocndz(nocndbdy) = iz
                  ocdelmx(nocndbdy) = delmx
                  ocdelmy(nocndbdy) = delmy
                  ocdelmz(nocndbdy) = delmz
                  ocdelpx(nocndbdy) = delpx
                  ocdelpy(nocndbdy) = delpy
                  ocdelpz(nocndbdy) = delpz
                  ocvolt(nocndbdy) = volt
                endif
              endif
            enddo
          enddo
        enddo
      endif

      return
      end
c===========================================================================
      subroutine srfrvinout(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                      lzend,xmin,xmax,ymin,ymax,lshell,
     &                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      use Surface_of_Rev
      use InGen, ONLY:fstype
      use InGen3d, ONLY:solvergeom,RZgeom,XYZgeom
      character(*):: rminofz,rmaxofz
      real(kind=8):: volt,zmin,zmax,xcent,ycent
      logical(ISZ):: lzend,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the outside of a surface of revolution.
c The Basis routines rminofz and rmaxofz should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvout("rofz",...)
c
c Input:
c   rminofz   is the name of the basis routine which calculates the radius as a
c             function of z.
c   rmaxofz   (same)
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   lzend     is the logical requesting that the end of the conductor be
c             included
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.

      real(kind=8):: zz,srfrv_f
      integer(ISZ):: i,nc,nceo

c     --- Table-ize the data if srfrv_pernz > 0.
      if (srfrv_pernz > 0 .and. .not. lsrminlinr) then
        npnts_srmin = (zmax - zmin)/dz*srfrv_pernz + 1
        zz = (zmax - zmin)/(npnts_srmin - 1)
        call gchange("Surface_of_Rev",0)
        do i=1,npnts_srmin
          z_srmin(i) = zmin + (i-1)*zz
          r_srmin(i) = srfrv_f(z_srmin(i),rminofz,1)
        enddo
        lsrminlinr = .true.
      endif
      if (srfrv_pernz > 0 .and. .not. lsrmaxlinr) then
        npnts_srmax = (zmax - zmin)/dz*srfrv_pernz + 1
        zz = (zmax - zmin)/(npnts_srmax - 1)
        call gchange("Surface_of_Rev",0)
        do i=1,npnts_srmax
          z_srmax(i) = zmin + (i-1)*zz
          r_srmax(i) = srfrv_f(z_srmax(i),rmaxofz,2)
        enddo
        lsrmaxlinr = .true.
      endif

c     --- Call the appropriate routine
      if(solvergeom==RZgeom) then
         call srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                     lzend,xmin,xmax,ymin,ymax,lshell,
     &                     zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                     ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      else if (fstype == 7 .or. fstype == 11) then
        call srfrvinout3dmg(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                      lzend,xmin,xmax,ymin,ymax,lshell,
     &                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
c     else if(fstype==11) then
c        call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
c    &                   lzend,xmin,xmax,ymin,ymax,lshell,
c    &                   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
c    &                   ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      else
        call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                    lzend,xmin,xmax,ymin,ymax,lshell,
     &                    zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      end if

c     --- Clean up the table-ized data.
      if (srfrv_pernz > 0) then
        lsrminlinr = .false.
        lsrmaxlinr = .false.
      endif

      return 
      end
c===========================================================================
      subroutine srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                        lzend,xmin,xmax,ymin,ymax,lshell,
     &                        zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                        ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      use PSOR3d
      character(*):: rminofz,rmaxofz
      real(kind=8):: volt,zmin,zmax,xcent,ycent
      logical(ISZ):: lzend,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

c Sets up a conductor represented by the outside of a surface of revolution.
c The Basis routines rminofz and rmaxofz should be of the form
c    function rofz
c      $z = srfrv_z
c      $r = f($z)
c      srfrv_r = $z
c    endf
c where f() is the radius as a function of z and srfrv_z and srfrv_r are
c compiled variables which pass data into and out of rofz.
c The begining of the calling sequence would then be
c call srfrvout("rofz",...)
c
c Input:
c   rminofz   is the name of the basis routine which calculates the radius as a
c             function of z.
c   rmaxofz   (same)
c   volt      is the voltage on the conductor.
c   zmin      is the minimum z of the conductor.
c   zmax      is the maximum z of the conductor.
c   xcent     is the x center of the conductor.
c   ycent     is the y center of the conductor.
c   lzend     is the logical requesting that the end of the conductor be
c             included
c   lshell    is the logical requesting that the shell be subgrid resolved
c   xmin,xmax,ymin,ymax
c             min and max transverse extent of conductor.
c   zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
c   xmesh,ymesh,l2symtry,l4symtry
c            are all variables describing the grid.
c Output is put directly into the conductor arrays of PSOR3d.

      integer(ISZ):: ix,iy,iz,i
      real(kind=8):: rad,rrrmax,rrrmin,rmin,rmax
      real(kind=8):: zgmin,zgmax
      real(kind=8):: rmin_zmin,rmin_zmax
      real(kind=8):: rmax_zmin,rmax_zmax
      real(kind=8):: rrrminlower,rrrminhigher,rrrmaxlower,rrrmaxhigher
      real(kind=8):: zlower,zhigher
      real(kind=8):: xxmin,xxmax,yymin,yymax
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: jxmin,jxmax,jymin,jymax
      integer(ISZ):: izbmin,izbmax,iterrf
      real(kind=8):: fuzz,xfuzz,yfuzz,rfuzz,zfuzz,xx,yy,zz
      real(kind=8):: delmz0,delpz0
      real(kind=8):: delmx,delmy,delmz,delpx,delpy,delpz
      real(kind=8):: max_dxdy,dxi,dyi,dzi
      real(kind=8):: ddz,znew,rnew
      real(kind=8):: srfrv_f
      real(kind=8):: x(2),f(2)

      max_dxdy = max(dx,dy)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c Set fuzz which is used to avoid problems with roundoff.
      fuzz = 1.e-10
      xfuzz = dx*fuzz
      yfuzz = dy*fuzz
      rfuzz = min(xfuzz,yfuzz)
      zfuzz = dz*fuzz

c Find maximal transverse extent of the conductors with an additional grid
c cell added to allow subgrid-scale resolution at the edge.
      jxmin = 0
      jymin = 0
      jxmax = nx
      jymax = ny

c Set limits for transverse loops based on the transverse boundary conditions
c and the symmetry.
      ixmin = max(0   ,jxmin)
      iymin = max(0   ,jymin)
      ixmax = min(nx  ,jxmax)
      iymax = min(ny  ,jymax)

c Calculate bounding box relative to center of surface of revolution.
      xxmin = xmin - xcent
      xxmax = xmax - xcent
      yymin = ymin - ycent
      yymax = ymax - ycent

c Set limits for longitudinal loops based on boundary conditions.  If
c using Dirichlet boundaries, don't put conductor at planes iz=0 and nz.
      izbmin = 0
      izbmax = nz

c Set limits for longitudinal loops.  They are set so that the minimum
c end is rounded up and maximum end is rounded down, ensuring that iz is
c within the conductor.
      zgmin = zmmin + zbeam
      zgmax = zmmax + zbeam
      izmin = nz+int((zmin-zgmax)*dzi-zfuzz)
      izmax =    int((zmax-zgmin)*dzi+zfuzz)

c Perform error checking, to make sure that izmin and izmax are within the
c grid.
      izmin = max(izbmin,izmin)
      izmax = min(izbmax,izmax)
     
c If izmin > izbmax or izmax < izbmin then the conductor is completely off
c of the grid, so return.
      if (izmin > izbmax .or. izmax < izbmin) then
        call remark("WARNING:srfrvout: conductor is completely off the grid.")
        return
      endif

c Set conductor points.
      do iz=izmin,izmax
        zz = iz*dz+zgmin
        zlower = max(zz-dz,zmin)
        zhigher = min(zz+dz,zmax)
c       --- Make sure there is enough room in the arrays.
        if (ncond + (1+nx)*(1+ny) > ncondmax) then
          ncondmax = ncond + (1+nx)*(1+ny)
          call gchange("PSOR3d",0)
        endif
c       --- Calculate the surface radius at iz.
        rrrmin = srfrv_f(zz,rminofz,1)
        rrrmax = srfrv_f(zz,rmaxofz,2)
c       --- Loop over transverse plane.
        do ix=ixmin,ixmax
          do iy=iymin,iymax
            xx = (ix-ix_axis)*dx-xcent
            yy = (iy-iy_axis)*dy-ycent
            rad = sqrt(xx**2 + yy**2)

c           --- The position must be within the min/max transverse extent
c           --- of the conductor.
            if ((rrrmin-rfuzz <= rad .and. rad <= rrrmax+rfuzz .and.
     &           xxmin <= xx .and. xx <= xxmax .and.
     &           yymin <= yy .and. yy <= yymax)) then
              ncond = ncond + 1
              ixcond(ncond) = ix
              iycond(ncond) = iy
              izcond(ncond) = iz
              condvolt(ncond) = volt
            endif
          enddo
        enddo
      enddo
c End of conductor points section

c Set subgrid scale arrays
      if (lcndbndy) then

c       --- Calculate r at zmin and zmax.
        rmin_zmin = srfrv_f(zmin,rminofz,1)
        rmin_zmax = srfrv_f(zmax,rminofz,1)
        rmax_zmin = srfrv_f(zmin,rmaxofz,2)
        rmax_zmax = srfrv_f(zmax,rmaxofz,2)

c       --- Recalculate limits for longitudinal loops.  Now, they are rounded
c       --- to include the plane just outside the conductor.
        izmin = max(izbmin,int((zmin-zgmin)*dzi-zfuzz))
        izmax = min(izbmax,nz+int((zmax-zgmax)*dzi+zfuzz))

c       -- Set all points along and just outside conductor.
        do iz=izmin,izmax
          zz = iz*dz+zgmin
          zlower = max(zz-dz,zmin)
          zhigher = min(zz+dz,zmax)

          if (necndbdy + (1+nx)*(1+ny) > ncndmax .or.
     &        nocndbdy + (1+nx)*(1+ny) > ncndmax) then
            ncndmax = max(necndbdy,nocndbdy) + (1+nx)*(1+ny)
            call gchange("PSOR3d",0)
          endif

c         --- Calculate radius of the surface.
c         --- When zz is outside the conductor, get radius at end of
c         --- conductor.
          if (zz < zmin-zfuzz) then
            rmin = rmin_zmin
            rmax = rmax_zmin
          elseif (zz > zmax+zfuzz) then
            rmin = rmin_zmax
            rmax = rmax_zmax
          else
            rmin = srfrv_f(zz,rminofz,1)
            rmax = srfrv_f(zz,rmaxofz,2)
          endif

c         --- Find radius one grid cell lower in z.  delmz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz > zmax, the delmz0 is set to the distance
c         --- between the z plane and the end of the conductor.
          rrrminlower = srfrv_f(zlower,rminofz,1)
          rrrmaxlower = srfrv_f(zlower,rmaxofz,2)
          delmz0 = 2.
          if (zz > zmax .and. lzend) then
            delmz0 = (zz - zmax)*dzi
          endif

c         --- Find radius one grid cell higher in z.  delpz0 is set to 2
c         --- in most cases, assuming that the conductor is far away in z.
c         --- When zz < zmin, the delmp0 is set to the distance
c         --- between the z plane and the end of the conductor.
          rrrminhigher = srfrv_f(zhigher,rminofz,1)
          rrrmaxhigher = srfrv_f(zhigher,rmaxofz,2)
          delpz0 = 2.
          if (zz < zmin .and. lzend) then
            delpz0 = (zmin - zz)*dzi
          endif

c         --- Find minimum of rrr and radius of surface at nearby z points.
          rrrmin = min(rmin-max_dxdy,rrrminlower,rrrminhigher)
          rrrmax = max(rmax+max_dxdy,rrrmaxlower,rrrmaxhigher)

          do ix=ixmin,ixmax
            do iy=iymin,iymax
              xx = (ix-ix_axis)*dx-xcent
              yy = (iy-iy_axis)*dy-ycent
              rad = sqrt(xx**2 + yy**2)

c             --- Calculate distances in z first.  Select points which
c             --- are transversely near the surface of revolution and
c             --- are within the bounding box.

c             --- Calculate distance to conductor at lower z.
              delmz = delmz0

c             --- If grid location is between rrr and rrrlower and between
c             --- rmax and rrrlower, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrminlower <= rad .and. rad < rmin-rfuzz
     &            .and. delmz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = min(zz,zmax)
                f(1) = rmin - rad
                x(2) = zlower
                f(2) = rrrminlower - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = min((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmax)
                  f(i) = srfrv_f(x(i),rminofz,1) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zlower <= x(i) .and. x(i) < zz) then
                  delmz = (zz - x(i))*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rminofz,1)
                    if (rnew > rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew <= rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz = (zz - znew)*dzi
                endif
              elseif (rrrmaxlower >= rad .and. rad > rmax+rfuzz
     &            .and. delmz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = min(zz,zmax)
                f(1) = rad - rmax
                x(2) = zlower
                f(2) = rad - rrrmaxlower
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = min((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmax)
                  f(i) = rad - srfrv_f(x(i),rmaxofz,2)
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zlower <= x(i) .and. x(i) < zz) then
                  delmz = (zz - x(i))*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rmaxofz,2)
                    if (rnew < rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew >= rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz = (zz - znew)*dzi
                endif
              elseif (rad < rrrminlower .or.
     &                rad > rrrmaxlower) then
                delmz = 2.
              endif

c             --- Calculate distance to conductor at higher z.
              delpz = delpz0

c             --- If grid location is between rrr and rrrhigher and between
c             --- rmax and rrrhigher, find intersection point.  The zfuzz
c             --- ensures that a distance of zero is not calculated.
              if (rrrminhigher <= rad .and. rad < rmin-rfuzz
     &            .and. delpz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = max(zz,zmin)
                f(1) = rmin - rad
                x(2) = zhigher
                f(2) = rrrminhigher - rad
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = max((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmin)
                  f(i) = srfrv_f(x(i),rminofz,1) - rad
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zz < x(i) .and. x(i) <= zhigher) then
                  delpz = (x(i) - zz)*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rminofz,1)
                    if (rnew <= rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew > rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz = (znew - zz)*dzi
                endif
              elseif (rrrmaxhigher >= rad .and. rad > rmax+rfuzz
     &            .and. delpz > zfuzz*dzi .and.
     &          xxmin <= xx .and. xx <= xxmax .and.
     &          yymin <= yy .and. yy <= yymax) then
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
                x(1) = max(zz,zmin)
                f(1) = rad - rmax
                x(2) = zhigher
                f(2) = rad - rrrmaxhigher
                iterrf = 0
                do while (abs(f(2))+abs(f(1)) > rfuzz .and. iterrf < 12)
                  i = 1
                  if (abs(f(2)) > abs(f(1))) i = 2
                  x(i) = max((x(1)*f(2)-x(2)*f(1))/(f(2)-f(1)),zmin)
                  f(i) = rad - srfrv_f(x(i),rmaxofz,2)
                  iterrf = iterrf + 1
                end do
                if (abs(f(2))+abs(f(1)) < rfuzz .and.
     &              zz < x(i) .and. x(i) <= zhigher) then
                  delpz = (x(i) - zz)*dzi
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,rmaxofz,2)
                    if (rnew >= rad) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew < rad) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz = (znew - zz)*dzi
                endif
              elseif (rad < rrrminhigher .or.
     &                rad > rrrmaxhigher) then
                delpz = 2.
              endif

c             --- Now calculate transverse distances.
              delmx = 2.
              delpx = 2.
              delmy = 2.
              delpy = 2.

c             --- Calculate distances to surface of revolution.
              if (rrrmin+rfuzz < rad .and. rad <= rmin+rfuzz .and.
     &            zmin <= zz .and. zz <= zmax) then
                delmx = (sqrt(max(0.,rmin**2 - yy**2)) + xx)*dxi
                delpx = (sqrt(max(0.,rmin**2 - yy**2)) - xx)*dxi
                delmy = (sqrt(max(0.,rmin**2 - xx**2)) + yy)*dyi
                delpy = (sqrt(max(0.,rmin**2 - xx**2)) - yy)*dyi
              elseif (rrrmax-rfuzz > rad .and. rad >= rmax-rfuzz .and.
     &            zmin <= zz .and. zz <= zmax) then
                if (yy**2 < rmax**2) then
                  delmx = ( xx - sqrt(max(0.,rmax**2 - yy**2)))*dxi
                  delpx = (-xx - sqrt(max(0.,rmax**2 - yy**2)))*dxi
                endif
                if (xx**2 < rmax**2) then
                  delmy = ( yy - sqrt(max(0.,rmax**2 - xx**2)))*dyi
                  delpy = (-yy - sqrt(max(0.,rmax**2 - xx**2)))*dyi
                endif
              endif

c             --- If requested, select points on shell around bounding box.
              if (lshell .and. zmin-zfuzz <= zz .and. zz <= zmax+zfuzz) then
c               --- Select points within 1 grid cell outside of bounding box
c               --- which are also within rrr and rmax.
                if (rmin-rfuzz <= rad .and. rad <= rmax+rfuzz) then
                  if (xx < xxmin) delpx = (xxmin-xx)*dxi
                  if (yy < yymin) delpy = (yymin-yy)*dyi
                  if (xx > xxmax) delmx = (xx-xxmax)*dxi
                  if (yy > yymax) delmy = (yy-yymax)*dyi
                endif
              endif

c             --- Now, ensure that the point on the surface from which the
c             --- distance was calculated is within the bounding box.
              if (xx-delmx*dx < xxmin) delmx = 2.
              if (xx+delpx*dx > xxmax) delpx = 2.
              if (yy-delmy*dy < yymin) delmy = 2.
              if (yy+delpy*dy > yymax) delpy = 2.

c             --- Make sure that distances are greater than zero.
              if (delmx < fuzz) delmx = 2.
              if (delpx < fuzz) delpx = 2.
              if (delmy < fuzz) delmy = 2.
              if (delpy < fuzz) delpy = 2.
              if (delmz < fuzz) delmz = 2.
              if (delpz < fuzz) delpz = 2.

c             -- Check if any of the distances are less than 1.
              if (delmx < 1.-fuzz .or. delpx < 1.-fuzz .or.
     &            delmy < 1.-fuzz .or. delpy < 1.-fuzz .or.
     &            delmz < 1.-fuzz .or. delpz < 1.-fuzz) then
                if (mod(ix+iy+iz+zparity,2) == 0) then
c                 --- even
                  necndbdy = necndbdy + 1
                  iecndx(necndbdy) = ix
                  iecndy(necndbdy) = iy
                  iecndz(necndbdy) = iz
                  ecdelmx(necndbdy) = delmx
                  ecdelmy(necndbdy) = delmy
                  ecdelmz(necndbdy) = delmz
                  ecdelpx(necndbdy) = delpx
                  ecdelpy(necndbdy) = delpy
                  ecdelpz(necndbdy) = delpz
                  ecvolt(necndbdy) = volt
                else
c                 --- odd
                  nocndbdy = nocndbdy + 1
                  iocndx(nocndbdy) = ix
                  iocndy(nocndbdy) = iy
                  iocndz(nocndbdy) = iz
                  ocdelmx(nocndbdy) = delmx
                  ocdelmy(nocndbdy) = delmy
                  ocdelmz(nocndbdy) = delmz
                  ocdelpx(nocndbdy) = delpx
                  ocdelpy(nocndbdy) = delpy
                  ocdelpz(nocndbdy) = delpz
                  ocvolt(nocndbdy) = volt
                endif
              endif
            enddo
          enddo
        enddo
      endif

      return
      end
c========================================================================
      real(kind=8) function srfrv_f(zz,rofzfunc,icase)
      use PSOR3d
      use Surface_of_Rev
      real(kind=8):: zz
      integer(ISZ):: icase
      character(*):: rofzfunc
c     This function returns the value of radius as a function of z.
c     That can be calculated either from a piecewise-linear curve, or from
c     the basis routine rofzfunc.
c     The radius returned is always greater than or equal to zero.

      integer(ISZ):: i,ii

      if (lsrlinr .and. icase == 0) then

c       --- Calculate r(z) from a piecewise-linear curve.

c       --- Find which piece z is in.
        ii = 0
        do i=1,npnts_sr
          if (zz > z_sr(i)) ii = i
        enddo

c       --- Calculate r from z.
        if (ii == 0) then
          srfrv_r = r_sr(1)
        elseif (ii == npnts_sr) then
          srfrv_r = r_sr(npnts_sr)
        else
          srfrv_r = r_sr(ii) + (zz - z_sr(ii))*
     &                   (r_sr(ii+1) - r_sr(ii))/(z_sr(ii+1) - z_sr(ii))
        endif

      elseif (lsrminlinr .and. icase == 1) then

c       --- Calculate rmin(z) from a piecewise-linear curve.

c       --- Find which piece z is in.
        ii = 0
        do i=1,npnts_srmin
          if (zz > z_srmin(i)) ii = i
        enddo

c       --- Calculate r from z.
        if (ii == 0) then
          srfrv_r = r_srmin(1)
        elseif (ii == npnts_srmin) then
          srfrv_r = r_srmin(npnts_srmin)
        else
          srfrv_r = r_srmin(ii) + (zz - z_srmin(ii))*
     &              (r_srmin(ii+1) - r_srmin(ii))/(z_srmin(ii+1) - z_srmin(ii))
        endif

      elseif (lsrmaxlinr .and. icase == 2) then

c       --- Calculate rmax(z) from a piecewise-linear curve.

c       --- Find which piece z is in.
        ii = 0
        do i=1,npnts_srmax
          if (zz > z_srmax(i)) ii = i
        enddo

c       --- Calculate r from z.
        if (ii == 0) then
          srfrv_r = r_srmax(1)
        elseif (ii == npnts_srmax) then
          srfrv_r = r_srmax(npnts_srmax)
        else
          srfrv_r = r_srmax(ii) + (zz - z_srmax(ii))*
     &              (r_srmax(ii+1) - r_srmax(ii))/(z_srmax(ii+1) - z_srmax(ii))
        endif

      else

c       --- Get r(z) from the basis function.
        srfrv_z = zz
        call execuser(rofzfunc)

      endif

      srfrv_f = max(0.,srfrv_r)

      return
      end
c===========================================================================
      subroutine setcndtr_rz(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                       l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

      l4symtry = .true.
      call setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                l2symtry,l4symtry)

      return
      end
c===========================================================================
      subroutine srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                       xmin,xmax,ymin,ymax,lshell,
     &                       zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                       ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

      call srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                xmin,xmax,ymin,ymax,lshell,
     &                zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

      return
      end
c===========================================================================
      subroutine srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                      lzend,xmin,xmax,ymin,ymax,lshell,
     &                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      character(*):: rminofz,rmaxofz
      real(kind=8):: volt,zmin,zmax,xcent,ycent
      logical(ISZ):: lzend,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

      call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                  lzend,xmin,xmax,ymin,ymax,lshell,
     &                  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                  ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

      return
      end
c===========================================================================
      subroutine setcndtr_3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                       l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

      l4symtry = .true.
      call setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                l2symtry,l4symtry)

      return
      end
c===========================================================================
      subroutine srfrvout_3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                    xmin,xmax,ymin,ymax,lshell,
     &                    zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                    ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

      call srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                xmin,xmax,ymin,ymax,lshell,
     &                zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

      return
      end
c===========================================================================
      subroutine srfrvinout_3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                      lzend,xmin,xmax,ymin,ymax,lshell,
     &                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
      character(*):: rminofz,rmaxofz
      real(kind=8):: volt,zmin,zmax,xcent,ycent
      logical(ISZ):: lzend,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry

      call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,
     &                  lzend,xmin,xmax,ymin,ymax,lshell,
     &                  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                  ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

      return
      end
c===========================================================================

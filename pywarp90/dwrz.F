#include "top.h"
c=============================================================================
c@(#) File DWRZ.M, version $Revision: 3.5 $, $Date: 2002/04/24 21:36:25 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is part of the package WRZ of code WARP
c  RZ electrostatic PIC code, cylindrical geometry, for beam problems
c  This file contains the diagnostic and plotting routines.
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c  Debbie Callahan, LLNL, (510)294-4066
c=============================================================================
      subroutine prntparz(lprntpara)
      use InMeshrz
      use Io
      use Ch_var
      logical(ISZ):: lprntpara
 
c  Prints out various parameters to a plot frame and an output file or tty
 

c Exit if paramters are not to be printed
      if (.not. lprntpara) return

c     --- Call script version of this routine
      call execuser("printparametersrz")
 
c  Write to tty
   30 format(1x,a,i8,a)
      write (STDOUT,30) "Number of grid points in r = ",nr," "
      write (STDOUT,30) "Number of grid points in z = ",nz," "
 
c  Write to text output file
      if (warpout > -1) then
        call edit (warpout, "nr")
        call edit (warpout, "nz")
      endif
 
      return
      end
c=============================================================================
      subroutine pltfldrz(fld,freqflag)
      use InPltCtlrz
      character(3):: fld
      integer(ISZ):: freqflag
 
c  Master control for doing field plots
 
      return
      end
c=============================================================================
      subroutine srhoaxrz
      use Z_arrays
      use InMeshrz
      use Picglbrz
      use Fieldsrz
c  Sets 1d array for the charge density on the axis
      integer(ISZ):: i
      do i=0,nz
        rhoax(i) = rho(0,i)
      enddo
      return
      end
c=============================================================================
      subroutine sphiaxrz
      use Z_arrays
      use InMeshrz
      use Picglbrz
      use Fieldsrz
c  Sets 1d array for the potential on the axis
      integer(ISZ):: i
      do i=0,nz
        phiax(i) = phi(0,i) + eoffrz*zmesh(i) - phiresist(i)
      enddo
      return
      end
c=============================================================================
      subroutine sezaxrz
      use Z_arrays
      use InMeshrz
      use Fieldsrz
      use Picglbrz
c  Sets 1d array for the space charge E field on the axis
      integer(ISZ):: i
      do i=0,nz
        ezax(i) = ezfld(0,i)
      enddo
      return
      end
c=============================================================================
      subroutine rhodiarz
      use InGenrz
      use InDiag
      use InMeshrz
      use Fieldsrz
      use Picglb
      use Picglbrz
      use Z_Moments
      use Win_Moments

      integer(ISZ):: iwin,iz,ir,ifwindow
      real(kind=8):: z1,z0,zwin

c  Sets rho window diagnostics (they need data from the rz mesh)
      ifwindow = 1
      if (ifwindow .eq. 1) then

        do iwin = 0,NWINDOWS
           if (zwindows(1,iwin).ne.zwindows(2,iwin)) then
              zwin = 0.5 * (zwindows(1,iwin) + zwindows(2,iwin))
              iz = (zwin - zmmin) / dz
              z1 = (zwin - zmmin) / dz - iz
              z0 = 1. - z1
              rhomid(iwin) = z0 * rho(0,iz)
     &                     + z1 * rho(0,iz+1)
              rhomax(iwin) = 0.
              do ir = 0,nr
                 rhomax(iwin) = max(rhomax(iwin),
     &                                z0*rho(ir,iz) + z1*rho(ir,iz+1))
              enddo
           endif
        enddo
      endif

c  Sets rho grid diagnostics (they need data from the rz mesh)
      if (ifzmmnt .gt. 0) then

c       --- rho on axis
        do iz = 0,nzmmnt
           rhomidz(iz) = rho(0,iz)
        enddo

c       --- max rho at each z (vectorized)
        do iz = 0,nzmmnt
           rhomaxz(iz) = -LARGEPOS
        enddo
        do ir = 0,nr
          do iz = 0,nzmmnt
            rhomaxz(iz) = max(rhomaxz(iz),rho(ir,iz))
          enddo
        enddo
      endif

      return
      end

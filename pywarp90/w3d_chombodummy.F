#include "top.h"
c=============================================================================
c dummy version of these routines which should get overridden by the linker
c if the libraries are put in the right order on the link line
c
      subroutine amrfieldsolve(iwhich)
      integer(ISZ) iwhich
      print*,  iwhich 
      print*,"*** W3DEXE: stopping.  dummy amrfieldsolve getting called"
      return
      end

      subroutine setAMRGrids(nx, ny, nz, dx, xmmin, ymmin, zmmin,
     &    numLevels, refratio, i, j, k, level, numtags)
      integer(ISZ) nx, ny, nz, numlevels
      integer(ISZ) i(*), j(*), k(*), refratio(*)
      double precision dx, xmmin, ymmin, zmmin
      integer(ISZ) level(*), numtags
      print*,"***: stopping.  dummy tagGrid called "
      return
      end

      subroutine returnphi(phiout)
      double precision phiout(*)
      print*,"***: stopping.  dummy returnphi called"
      return
      end

      subroutine returnphic(phiout)
      double precision phiout(*)
      print*,"***: stopping.  dummy returnphic called"
      return
      end

      subroutine returnrho(rhoout)
      double precision rhoout(*)
      print*,"***: stopping.  dummy returnrho called"
      return
      end

c=============================================================================
c dummy version of these routines which should get overridden by the linker
c if the libraries are put in the right order on the link line
c
      subroutine amrfieldsolve(iwhich)
      integer iwhich
      print*,  iwhich 
      print*,"*** W3DEXE: stopping.  dummy amrfieldsolve getting called"
      return
      end

      subroutine setAMRGrid(nx, ny, nz, dx, xmmin, ymmin, zmmin,
     &    numLevels, refratio, i, j, k, level, numtags)
      integer nx, ny, nz, numlevels
      integer i(*), j(*), k(*), refratio(*)
      double precision dx, xmmin, ymmin, zmmin
      integer level(*), numtags
      print*,"***: stopping.  dummy tagGrid called "
      return
      end


#include "top.h"
c=============================================================================
c@(#) File CIR.M, version $Revision: 3.7 $, $Date: 2001/07/18 23:54:05 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  It is the source file for the package CIR of the PIC code WARP,
c  but it may be useful by itself.  It consists of the envelope/fluid
c  model which originated in CIRCE.
c
c=============================================================================
      subroutine cirinit
      use CIRversion

c  Called at first reference to package (not nec. a "run" etc.).

      call cirvers (STDOUT)

      return
      end
c=============================================================================
      subroutine cirvers (iout)
      use CIRversion
      integer(ISZ):: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"CIRCE",verscir)
      return
      end
c=============================================================================
      subroutine cirgen()
      use Constant
      use Beam_acc
      use InGen
      use InPart
      use Picglb
      use Lattice
      use LatticeInternal
      use CIRvars
      use CIRflags
      use CIRtmp
      use CIRfield
      use CIRhist

c  Invoked by the GENERATE command, it allocates storage
c  for CIRCE.

c  Send a startup message to the standard output

      call remark(" ***  CIRCE generating")

c  Calculate derived quantities

      call derivqty

c     --- allocate dynamic arrays ("group change")
      if (nit == 0) nit = 101
      call gchange("CIRvars",0)
      nittmp = nit
      call gchange("CIRtmp",0)
      nitfield = nit
      call gchange("CIRfield",0)
      lhcir = 100
      if (nithist == 0) nithist = nit
      call gchange("CIRhist",0)

c     --- reset Lattice arrays to final size
      call resetlat

c     --- Initialize the LatticeInternal arrays
      if (ltdependent) then
        nzl = nit + 2
        nzlmax = nzl
        dzl = (zimax - zimin)/nit
        zlmin = zimin - dzl
        zlmax = zimax - dzl
        dzli = 1./dzl
        call gchange("LatticeInternal",0)
      else
        nzl = 0
        zlmin = 0.
        zlmax = 0.
      endif
      call setlattzt(zlcir,time,0)

c     --- Create an initial beam
      call initbeam

      return
      end
c=============================================================================
      subroutine cirexe()
      use CIRvars
      use CIRgeneral
      use CIRhist
      use CIRflags

c  Invoked by the STEP command, it runs the CIRCE code.
c  passes sigma and sigma0 back to TOP through OutParams

      real(kind=8):: wtimeoff

c  Send a startup message to the standard output

      if (lcirout) then
        call remark(" ***  CIRCE running")
      endif

c  Execute envelope solve, and time it
      call wtimeon
      call cirx
      cirtime = wtimeoff()
 
c     --- Resize history arrays
      if (lsavehist .and. lhcir > jhcir) then
        lhcir = jhcir
        call gchange("CIRhist",0)
      endif

      return
      end
c=============================================================================
      subroutine cirfin()

c  Deallocates dynamic storage used by test driver

      call gfree("CIRvars")
      call gfree("CIRtmp")
      call gfree("CIRfield")
      call gfree("CIRhist")

      return
      end
c=============================================================================
      subroutine initbeam()
      use Constant
      use Beam_acc
      use Picglb
      use InGen
      use InPart
      use CIRgeneral
      use CIRvars
      use CIRvartmp
      use CIRflags
      use CIRtmp
      use CIRfield
      use CIRhist

c  For convenience, this routine accesses the WARP database and is not
c  linked to the universe solely by its calling sequence.

      real(kind=8):: currise,curfall,deltbeam

c     --- allocate dynamic arrays ("group change")
      if (nit == 0) nit = 101
      call gchange("CIRvars",0)
      nittmp = nit
      call gchange("CIRtmp",0)
      nitfield = nit
      call gchange("CIRfield",0)
      if (nithist == 0) nithist = nit
      call gchange("CIRhist",0)

c     --- Set initial conditions
      currise = (1. - straight)/2.
      curfall = (1. - straight)/2.
      deltbeam = (zimax - zimin)/dvnz(vbeam)
      call setbeam(nit,var,a0,b0,ap0,bp0,x0,y0,xp0,yp0,
     &             ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot,
     &             vbeam/clight,gammabar,vtilt,tiltmid,tiltlen,
     &             emitx,emity,time,zbeam,aion,zion,
     &             icurload,iprofile,ivelload,igrid,iemload,lperveance,
     &             lendzero)

      return
      end
c=============================================================================
      subroutine cirx()
      use Constant
      use Beam_acc
      use Picglb
      use InGen
      use InPart
      use Lattice
      use LatticeInternal
      use Mult_data
      use CIRgeneral
      use CIRvars
      use CIRvartmp
      use CIRflags

c  For convenience, this routine accesses the lattice database and is not
c  linked to the universe solely by its calling sequence.

      integer(ISZ):: ii

c     --- Main loop of envelope code
c     --- Run the CIRCE kernel
      call cirrun(nit,var,zlcir,zucir,dscir,nscir,aion,zion,
     &            icharge,lezbeam,lperveance,lemittance,lallez,llinear,limage,
     &            lendzero,lsavehist,lfixed)

c     --- Copy the data into separated arrays.
      do ii=1,nit
        acir(ii)   = var( 1,ii)
        apcir(ii)  = var( 2,ii)
        bcir(ii)   = var( 3,ii)
        bpcir(ii)  = var( 4,ii)
        xcir(ii)   = var( 5,ii)
        xpcir(ii)  = var( 6,ii)
        ycir(ii)   = var( 7,ii)
        ypcir(ii)  = var( 8,ii)
        tcir(ii)   = var( 9,ii)
        vzcir(ii)  = var(10,ii)
        enxcir(ii) = var(11,ii)
        enycir(ii) = var(12,ii)
        cur(ii)    = var(13,ii)
        dq(ii)     = var(14,ii)
        den(ii)    = var(15,ii)
      enddo

      return
      end
c=============================================================================
      subroutine cirrun(nit,y,zlcir,zucir,dscir,nscir,aion,zion,
     &            icharge,lezbeam,lperveance,lemittance,lallez,llinear,limage,
     &            lendzero,lsavehist,lfixed)
      use Constant
      use Picglb
      use LatticeInternal
      use CIRvartmp
      use CIRfield
      use CIRhist
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: zlcir,zucir,dscir,aion,zion
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,llinear,limage
      logical(ISZ):: lendzero,lsavehist,lfixed

c Run the CIRCE kernel, using the inputted array, y and the starting data.

      real(kind=8):: s,ds
      integer(ISZ):: nscir,il
      logical(ISZ):: ltdependent = .false.

c     --- Make sure that the sign of the step size and the relation
c     --- between zlcir and zucir are consistent.
      if (dscir > 0. .and. zucir < zlcir) then
        call remark("Error: Either make dscir < 0 or make zlcir < zucir")
        return
      endif
      if (dscir < 0. .and. zlcir < zucir) then
        call remark("Error: Either make dscir > 0 or make zlcir > zucir")
        return
      endif

c     --- Allocate temp space
      nitvartmp = nit
      call gchange("CIRvartmp",0)
      nitfield = nit
      call gchange("CIRfield",0)

c     --- Make sure that nzl is zero (some other package may have changed it).
      if (.not. ltdependent) nzl = 0

c  Set fundamental constants, derived quantities, initial conditions

c     --- Make history arrays big enough
      if (lsavehist) then
        il = max(100,100*int((zucir - zlcir)/dscir/max(1,nhcir)/100. + 1))
        il = max(lhcir,il)
        if (nithist == 0) nithist = nit
        if (il > lhcir) then
          lhcir = il
          call gchange("CIRhist",0)
        endif
c       jhcir = -1
      endif

c     --- Make sure that the input is reasonable
      call checkenvelope(y,nit)

c     --- Save the initial history
      call savehistcir(zlcir,y,nit,nscir,lsavehist)

c     --- Set the starting position in the lab frame.
      s = zlcir

c     --- Check that the distance traveled is less than the final distance
c     --- to travel.
      do while (abs(s-zlcir) < abs(zucir-zlcir))

        nscir = nscir + 1

c       --- If doing fast method, set field and pipe quantities now.
c       if(ifast.ne.0)then
c         call getpipe
c         call circegetimage()
c         call getezbeam()
c       endif

c       --- Advance the envelope and orbits;
        call circesetlatt(s,time)
        call setstep(s,ds,zlcir,zucir,dscir)
        call integrt(y,yt,dy1,dy2,dy3,nit,s,ds,aion,zion,
     &            icharge,lezbeam,lperveance,lemittance,lallez,llinear,limage,
     &            lendzero,lfixed)

        s = s + ds

        call savehistcir(s,y,nit,nscir,lsavehist)

      end do

      return
      end
c=============================================================================
      subroutine circesetlatt(z,time)
      use CIRvars
      use CIRflags
      use LatticeInternal
      real(kind=8):: z,time

c Adjusts the internal lattice frame to cover the whole beam before setting
c up the lattice.

      if (ltdependent) then
        dzl = (var(9,nit) - var(9,1))/nit
        zlmin = var(9,1) - dzl
        zlmax = var(9,nit) - dzl
        dzli = 1./dzl
      else
        nzl = 0
        zlmin = 0.
        zlmax = 0.
      endif

      call setlattzt(z,time,0)

      return
      end
c=============================================================================
      real(kind=8) function dstoelement(s,zs,ze,ds)
      real(kind=8):: s,zs,ze,ds
c Compare ds to distance to end or start of element

c     --- Default value
      dstoelement = ds

      if (ds > 0.) then
c       --- Forward simulation
        if (zs <= s .and. s < ze) then
          dstoelement = min(ds,ze - s)
        elseif (s < zs) then
          dstoelement = min(ds,zs - s)
        endif

      else
c       --- Backward simulation
        if (zs < s .and. s <= ze) then
          dstoelement = max(ds,zs - s)
        elseif (s > ze) then
          dstoelement = max(ds,ze - s)
        endif

      endif

      return
      end
c=============================================================================
      subroutine setstep(s,ds,zlcir,zucir,dscir)
      use Lattice
      use LatticeInternal
      real(kind=8):: s,ds,zlcir,zucir,dscir

c Calculate the step size to the start (or end) of the next (current)
c lattice element.

      real(kind=8):: dstoelement
      integer(ISZ):: io

c     --- Next step size can not go beyond of range of calculation.
      if (zlcir <= zucir) then
c       --- Forward simulation
        ds = min(dscir,zucir - s)
      else
c       --- Backward simulation
        ds = max(dscir,zucir - s)
      endif

c     --- Hard edge quadrupoles
      if (quads) then
        do io=1,nquadol
          ds = dstoelement(s,cquadzs(0,io),cquadze(0,io),ds)
        enddo
      endif

c     --- Hard edge multipoles, can include linearly focusing
c     --- quadrupole or azimuthally symmetric field
      if (heles) then
        do io=1,nheleol
          ds = dstoelement(s,chelezs(0,io),cheleze(0,io),ds)
        enddo
      endif

c     --- Extract any linear focusing field component from the multipole
c     --- (emlt and mmlt) elements.
      if (emlts) then
        do io=1,nemltol
          ds = dstoelement(s,cemltzs(0,io),cemltze(0,io),ds)
        enddo
      endif
      if (mmlts) then
        do io=1,nmmltol
          ds = dstoelement(s,cmmltzs(0,io),cmmltze(0,io),ds)
        enddo
      endif

c     --- Get accelerating gradient from accl gaps.
      if (accls) then
        do io=1,nacclol
          ds = dstoelement(s,cacclzs(0,io),cacclze(0,io),ds)
        enddo
      endif

c     --- Get inverse radius of curvature from bend elements
      if (bends) ds = dstoelement(s,cbendzs(0),cbendze(0),ds)

c     --- Get bending field from dipole elements
      if (dipos) then
        do io=1,ndipool
          ds = dstoelement(s,cdipozs(0,io),cdipoze(0,io),ds)
        enddo
      endif

      return
      end
c=============================================================================
      subroutine extebcir(z,dz,dzfrac,y,nit)
      use Constant
      use Beam_acc
      use Picglb
      use Lattice
      use LatticeInternal
      use Mult_data
      use CIRfield
      integer(ISZ):: nit
      real(kind=8):: z,dz,dzfrac,y(16,nit)

c  Computes the forces needed for the general-lattice envelope calculation.
c  Includes both hard edged quadrupole elements 
c  and hard-edge multipole elements, and elements defined by their 
c  multipoles.  Only linear forces consistent with a truncated second-order 
c  envelope formulation are calculated.  These include a linear applied 
c  forces (electric and magnetic quadrupole fields, and a radial 
c  E-field proportional to the radial coodrinate r), and 
c  linear space-charge forces.  
c  Input:
c     z    is the starting location of the current time step
c     dz   is the distance to the intermediate step in the Runge-Kutta step
c     y    is the array hold the envelope information
c     nit  is the number of slices
c  Note that z, dz, and dzfrac are passed in seperately, instead of just
c  using z+dz. z is used to check whether the beam is in an element.
c  z+dz*dzfrac is used to as the beam position to fetch the fields which
c  depend on z location (emlt and mmlt).
c  Note that the check on whether z is within an elemnent needs to also check
c  on the sign of dz. The logic elsewhere in the code assumes that when
c  z == start of element, it is in, and when z == end of element, it is out.
c  But which end of the element is the start and which is the end depends on
c  the sign of dz, whether the beam is moving forward or backward.

      real(kind=8):: dedrr,dbdxq,dedxq
      integer(ISZ):: io,iele,ii,iem,imm,iim,iie,ia
      real(kind=8):: rotfuzz,wie,wim,wi

c     --- Get the current lattice parameters
      call circesetlatt(z,time)

c     --- Default value of rpipe
      if (rwall > 0.) rpipe = rwall

c     --- Type of transverse boundaries - used in calculation of
c     --- image forces. In all cases, rpipe is the aperture.
c     ---   0 round pipe (default value)
c     ---   1 electric quadrupole rods
c     ---   2 dipole plates
      ipipetype = 0

c     --- Zero out the fields
      do ii=1,nit
        ex(ii) = 0.
        ey(ii) = 0.
        ez(ii) = 0.
        ears(ii) = 0.
        bx(ii) = 0.
        by(ii) = 0.
        bz(ii) = 0.
        dedx(ii) = 0.
        dbdx(ii) = 0.
      enddo

c---------------------------------------------------------------------------
c     --- Hard edge quadrupoles

      if (quads) then

        do io=1,nquadol
c         --- If in a quad, get field, pipe, offset, and images.
          if (((cquadzs(0,io) <= z .and. z <  cquadze(0,io)) .and. dz>0.) .or.
     &        ((cquadzs(0,io) <  z .and. z <= cquadze(0,io)) .and. dz<0.)) then
            iele = cquadid(0,io)
            do ii=1,nit
              dedx(ii) = dedx(ii) + quadde(iele)
              dbdx(ii) = dbdx(ii) + quaddb(iele)
            enddo
            xpipeerr = cqoffx(0,io)
            ypipeerr = cqoffy(0,io)
            if (cquadde(0,io) .ne. 0) ipipetype = 1
            if (quadap(iele) > 0.) rpipe = quadap(iele)

c           --- Get time varying accelerating field.
            if (ntquad > 0) then
              do ii=1,nit
                ia = int((y(9,ii) - quadts(iele))/quaddt(iele))
                if (0 <= ia .and. ia < ntquad) then
                  wi = (y(9,ii) - quadts(iele))/quaddt(iele) - ia
                  dedx(ii) = dedx(ii) + quadet(ia  ,iele)*(1. - wi) +
     &                                  quadet(ia+1,iele)*      wi
                  dbdx(ii) = dbdx(ii) + quadbt(ia  ,iele)*(1. - wi) +
     &                                  quadbt(ia+1,iele)*      wi
                endif
              enddo
            endif

          endif
        enddo

      endif 

c---------------------------------------------------------------------------
c     --- Hard edge multipoles, can include linearly focusing 
c     --- quadrupole or azimuthally symmetric field 
 
      if (heles) then

        do io=1,nheleol
          if (((chelezs(0,io) <= z .and. z <  cheleze(0,io)) .and. dz>0.) .or.
     &        ((chelezs(0,io) <  z .and. z <= cheleze(0,io)) .and. dz<0.)) then
            iele = cheleid(0,io)
            if (heleap(iele) > 0.) rpipe = heleap(iele)
            xpipeerr = heleox(iele)
            ypipeerr = heleoy(iele)

            rotfuzz = 1.e-3
            dbdxq = 0.
            dedxq = 0.
            dedrr = 0.
            do ii = 1, max(helenm(iele),helene(iele)) 
              if ((nint(hele_n(ii,iele)) == 2) .and. 
     &            (nint(hele_v(ii,iele)) == 0)) then
c               --- quadrupole field  
                if ( abs(cos(helepm(ii,iele))*sin(helepm(ii,iele))) > rotfuzz )
     &         call remark("Warning: magnetic multipole structure incompatible")
                if ( abs(cos(helepe(ii,iele))*sin(helepe(ii,iele))) > rotfuzz )
     &         call remark("Warning: electric multipole structure incompatible")
                dbdxq = dbdxq + heleam(ii,iele)*cos(helepm(ii,iele))
                dedxq = dedxq + heleae(ii,iele)*cos(helepe(ii,iele))
                if (heleae(ii,iele) .ne. 0) ipipetype = 1
              endif 
              if ((nint(hele_n(ii,iele)) == 0) .and. 
     &            (nint(hele_v(ii,iele)) == 0)) then 
c               --- uniform field 
                dedrr = dedrr - 0.5*heleep(ii,iele)
              endif  
            enddo             

            do ii=1,nit
              dedx(ii) = dedx(ii) + dedxq
              dbdx(ii) = dbdx(ii) + dbdxq
            enddo

          endif
        enddo
      endif 

c---------------------------------------------------------------------------
c     --- Extract any linear focusing field component from the multipole 
c     --- (emlt and mmlt) elements.

      if (emlts) then

        do io=1,nemltol
c         --- If z is within an element, fetch the field.
          if (((cemltzs(0,io) <= z .and. z <  cemltze(0,io)) .and. dz>0.) .or.
     &        ((cemltzs(0,io) <  z .and. z <= cemltze(0,io)) .and. dz<0.)) then
            iele = cemltid(0,io)

c           --- Make sure the element was defined.
            iem = emltid(iele)
            if (iem > 0) then
              if (emltap(iele) > 0.) rpipe = emltap(iele)
              xpipeerr = cemltox(0,io)
              ypipeerr = cemltoy(0,io)
              ipipetype = 1

c             --- Find location of z in the electric multipole grid.
              iie = int((z + dz*dzfrac - cemltzs(0,io))/dzemlt(iem))
              wie =     (z + dz*dzfrac - cemltzs(0,io))/dzemlt(iem) - iie

c             --- Lookup field strengths with linear interpolation.
              dedxq = 0.
              dedrr = 0.
              do ii=1,nesmult
                if ((nint(emlt_n(ii)) == 2) .and. 
     &              (nint(emlt_v(ii)) == 0)) then
                  dedxq = dedxq + (esemlt(iie  ,ii,iem)*(1.-wie)+
     &                             esemlt(iie+1,ii,iem)*wie)*
     &                            (emltsc(iele) + emltsf(iele))
                endif
                if ((nint(emlt_n(ii)) == 0) .and. 
     &              (nint(emlt_v(ii)) == 0)) then
                  dedrr = dedrr - 0.5*(esemltp(iie  ,ii,iem)*(1.-wie)+
     &                                 esemltp(iie+1,ii,iem)*wie)*
     &                                (emltsc(iele) + emltsf(iele))
                endif
              enddo

              do ii=1,nit
                dedx(ii) = dedx(ii) + dedxq
              enddo

            endif
          endif
        enddo
      endif

      if (mmlts) then

        do io=1,nmmltol
c         --- If z is within an element, fetch the field.
          if (((cmmltzs(0,io) <= z .and. z <  cmmltze(0,io)) .and. dz>0.) .or.
     &        ((cmmltzs(0,io) <  z .and. z <= cmmltze(0,io)) .and. dz<0.)) then
            iele = cmmltid(0,io)

c           --- Make sure the element was defined.
            imm = mmltid(iele)
            if (imm > 0) then
              if (mmltap(iele) > 0.) rpipe = mmltap(iele)
              xpipeerr = cmmltox(0,io)
              ypipeerr = cmmltoy(0,io)

c             --- Find location of z in the magnetic multipole grid.
              iim = int((z + dz*dzfrac - cmmltzs(0,io))/dzmmlt(imm))
              wim =     (z + dz*dzfrac - cmmltzs(0,io))/dzmmlt(imm) - iim

c             --- Lookup field strength with linear interpolation.
              dbdxq = 0.
              do ii=1,nmsmult
                if ((nint(mmlt_n(ii)) == 2) .and. 
     &              (nint(mmlt_v(ii)) == 0)) then
                  dbdxq = dbdxq +  (msmmlt(iim  ,ii,imm)*(1.-wim)+
     &                              msmmlt(iim+1,ii,imm)*wim)*
     &                             (mmltsc(iele) + mmltsf(iele))
                endif
              enddo

              do ii=1,nit
                dbdx(ii) = dbdx(ii) + dbdxq
              enddo

            endif
          endif
        enddo
      endif

c---------------------------------------------------------------------------
c     --- Get accelerating gradient from accl gaps.
      if (accls) then

        do io=1,nacclol
          if (((cacclzs(0,io) <= z .and. z <  cacclze(0,io)) .and. dz>0.) .or.
     &        ((cacclzs(0,io) <  z .and. z <= cacclze(0,io)) .and. dz<0.)) then
            iele = cacclid(0,io)
            do ii=1,nit
              ez(ii) = ez(ii) + acclez(iele)
            enddo
            if (acclap(iele) > 0.) rpipe = acclap(iele)

c           --- Get time varying accelerating field.
            if (ntaccl > 0) then
              do ii=1,nit
                ia = int((y(9,ii) - acclts(iele))/accldt(iele))
                if (0 <= ia .and. ia < ntaccl) then
                  wi = (y(9,ii) - acclts(iele))/accldt(iele) - ia
                  ez(ii) = ez(ii) + acclet(ia  ,iele)*(1. - wi) +
     &                              acclet(ia+1,iele)*      wi
                endif
              enddo
            endif

          endif

        enddo
      endif

c---------------------------------------------------------------------------
c     --- Get inverse radius of curvature from bend elements
      if (bends) then

        if (((cbendzs(0) <= z .and. z <  cbendze(0)) .and. dz > 0.) .or.
     &      ((cbendzs(0) <  z .and. z <= cbendze(0)) .and. dz < 0.)) then
          iele = cbendid(0)
          bendcurv = 1./bendrc(iele)
          if (bendap(iele) > 0.) rpipe = bendap(iele)
          sinrot = 0.
          cosrot = 1.
        endif

      endif

c---------------------------------------------------------------------------
c     --- Get inverse radius of curvature from bend elements
      if (drfts) then

        do io=1,ndrftol
          if (((cdrftzs(0,io) <= z .and. z <  cdrftze(0,io)) .and. dz>0.) .or.
     &        ((cdrftzs(0,io) <  z .and. z <= cdrftze(0,io)) .and. dz<0.)) then
            iele = cdrftid(0,io)
            if (drftap(iele) > 0.) rpipe = drftap(iele)
          endif
        enddo

      endif

c---------------------------------------------------------------------------
c     --- Get bending field from dipole elements
      if (dipos) then

        do io=1,ndipool
          if (((cdipozs(0,io) <= z .and. z <  cdipoze(0,io)) .and. dz>0.) .or.
     &        ((cdipozs(0,io) <  z .and. z <= cdipoze(0,io)) .and. dz<0.)) then
            iele = cdipoid(0,io)
            do ii=1,nit
              ex(ii) = dipoex(iele)
              ey(ii) = dipoey(iele)
              bx(ii) = dipobx(iele)
              by(ii) = dipoby(iele)
            enddo
            if (dipoap(iele) > 0.) rpipe = dipoap(iele)
            if (dipoex(iele) .ne. 0) ipipetype = 2
          endif
        enddo

      endif

c---------------------------------------------------------------------------

c --- Print warning if rpipe is zero
      if (rpipe == 0.) then
        call remark("WARNING: rpipe is zero, make sure that either top.rwall")
        call remark("and/or the lattice element apertures are set.")
      endif

      return
      end
c=============================================================================
      subroutine derivs(s,ds,dsfrac,y,dy,nit,aion,zion,icharge,
     &                  lezbeam,lperveance,lemittance,lallez,llinear,
     &                  limage,lendzero,lfixed)
      use Constant
      use CIRfield
      use CIRtmp
      integer(ISZ):: nit
      real(kind=8):: s,ds,dsfrac,y(16,nit),dy(16,nit)
      real(kind=8):: aion,zion
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero
      logical(ISZ):: lfixed

c  Computes the forces needed for the general-lattice envelope calculation.
c  Includes both residence corrected hard edged quadrupole elements 
c  and hard-edge multipole elements, and elements defined by their 
c  multipoles.  Only linear forces consistent with a truncated second-order 
c  envelope formulation are calculated.  These include a linear applied 
c  forces (electric and magnetic quadrupole fields, and a radial 
c  E-field proportional to the radial coodrinate r), and 
c  linear space-charge forces.  

c  Note: CIRtmp is needed for the fixed field solve only

      integer(ISZ):: ii
      real(kind=8):: alfcuri, denom
      real(kind=8):: y1i,y3i,y1py3i,vval,vvali,xval,yval,gam2,gam2i,gambetai
      real(kind=8):: rigidityi,bendfldx,bendfldy,bendex,bendey,quadfld,perv
      real(kind=8):: xemit,yemit,xbend,dtds,bend,dlnbeta,dbetads
      real(kind=8):: rpipei2
      logical(ISZ):: lfailed

      lfailed = .false.
      
c     --- Get the applied fields.
      call extebcir(s,ds,dsfrac,y,nit)

c     --- Get image forces based on pipe type set above.
c     --- This subroutine call puts the image forces into the temporary
c     --- arrays in the group CIRfield.
c     ---  fx,fy,gxx,gxy,gyx,gyy
      call circegetimage(y,nit,llinear,limage,lendzero)

c     --- Get self axial fields of beam.
      call getezbeam(nit,y,ezbeam,rpipe,icharge,lezbeam,lendzero,lfixed,lfailed)
      if (lfailed) then
        print *, "Slice overtaking has occurred"
      endif

c     --- Set some temporaries.
      rpipei2 = 1./rpipe**2
      alfcuri = zion*echarge/(4.*pi*eps0*amu*aion*clight**3)

c---------------------------------------------------------------------------
c     --- Loop over time slices of beam
      do ii=1,nit

c       --- Beam size reciprocals
        if (lendzero .and. (ii==1 .or. ii==nit)) then
          y1i = 0.
          y3i = 0.
          y1py3i = 0.
        else
          y1i = 1./y(1,ii)
          y3i = 1./y(3,ii)
          y1py3i = 1./(y(1,ii) + y(3,ii))
        endif

c       --- longitudinal velocity
        vval = clight*y(10,ii)
        vvali = 1./vval

c       --- centroid coordinates in displaced elements
        xval = y(5,ii) - xpipeerr
        yval = y(7,ii) - ypipeerr

c       --- relativistic factor gam2 = gamma**2 and gambeta = gamma*beta
        gam2 = 1.0/(1.0 - y(10,ii)**2)
        gam2i = (1.0 - y(10,ii)**2)
        gambetai = sqrt(gam2i)/y(10,ii)

c       --- magnetic rigidity reciprocal
        rigidityi = gambetai*zion*echarge/(aion*amu*clight)

c       --- bend-force terms
        bendfldx = (ex(ii)*vvali - by(ii))*rigidityi
        bendfldy = (ey(ii)*vvali + bx(ii))*rigidityi

        bendex = ex(ii)*vvali*rigidityi
        bendey = ey(ii)*vvali*rigidityi

c       --- quadrupole-force coefficient
        quadfld = ( - dedx(ii)*vvali + dbdx(ii))*rigidityi

c       --- perveance perv (labeled K in the envelope equation)
        if (lperveance) then
          perv = 2.0*y(13,ii)*(alfcuri*gambetai**3)
        else
          perv = 0.
        endif

c       --- normalized emittance to unnormalized emittance
        if (lemittance) then
          xemit = gambetai*y(11,ii)
          yemit = gambetai*y(12,ii)
        else
          xemit = 0.
          yemit = 0.
        endif

c       --- d ln(beta)/ds factor
        xbend = bendcurv*(cosrot*y(5,ii) + sinrot*y(7,ii))
        dtds =  (1.0 + xbend)*vvali

        bend = bendcurv*(cosrot*y(6,ii) + sinrot*y(8,ii))*gam2i/(1.0 + xbend)
        if (lallez) then
          dlnbeta = dtds*(ez(ii) + ezbeam(ii))*(gam2i*rigidityi)
        else
          dlnbeta = 0.
        endif

c       --- Fix dlnbeta
        if (lfixed) then
          denom = 1. - gtemp(ii) * alfcuri * y(13,ii) * gambetai**3
          if (denom < 0.) then
            print *, "Longitudinal current limit exceeded"
            lfailed = .true.
          endif
          dlnbeta = dlnbeta / denom
        endif

        dbetads = y(10,ii)*(dlnbeta - bend)
        if (llinear) then
          dlnbeta = gam2*dlnbeta
        else
          dlnbeta = gam2*(dlnbeta - bend)
        endif

c       --- check if slice overtaking has occurred
c       --- or the current limit was exceeded;
c       --- if so, print the current position s
        if (lfailed) then
          print *, " at s =", s
          stop
        endif

c       --- Now, set derivatives
c       --- da/ds
        dy(1,ii) = y(2,ii)

c       --- d(da/ds)/ds
        dy(2,ii) = - quadfld*y(1,ii)
     .                  + 2.0*perv*y1py3i
     .                  + xemit**2*y1i**3
     .                  + fx(ii)*perv*y(1,ii)*rpipei2
     .                  + (cosrot*bendcurv*y(1,ii))
     .                   *(cosrot*bendcurv + 2.0*bendfldx)
     .                  - dlnbeta*y(2,ii)

c       --- db/ds
        dy(3,ii) = y(4,ii)

c       --- d(db/ds)/ds
        dy(4,ii) = quadfld*y(3,ii)
     .                 + 2.0*perv*y1py3i
     .                 + yemit**2*y3i**3
     .                 + fy(ii)*perv*y(3,ii)*rpipei2
     .                 + (sinrot*bendcurv*y(3,ii))
     .                  *(sinrot*bendcurv + 2.0*bendfldy)
     .                 - dlnbeta*y(4,ii)

c       --- dx/ds
        dy(5,ii) = y(6,ii)

c       --- d(dx/ds)/ds
        dy(6,ii) = -quadfld*(y(5,ii) - xpipeerr)
     .                 + (perv*rpipei2)*(gxx(ii)*xval + gxy(ii)*yval)
     .                 + (1.0 + xbend)*(cosrot*bendcurv + bendex)
     .                 - (1.0 + 2.0*xbend)*by(ii)*rigidityi
     .                 - dlnbeta*y(6,ii)

c       --- dy/ds
        dy(7,ii) = y(8,ii)

c       --- d(dy/ds)/ds
        dy(8,ii) = quadfld*(y(7,ii) - ypipeerr)
     .                  + (perv*rpipei2)*(gyx(ii)*xval + gyy(ii)*yval)
     .                  + (1.0 + xbend)*(sinrot*bendcurv + bendey)
     .                  + (1.0 + 2.0*xbend)*bx(ii)*rigidityi
     .                  - dlnbeta*y(8,ii)

c       --- dt/ds
        dy(9,ii) = dtds

c       --- dbeta/ds
        dy(10,ii) = dbetads

      enddo

      return
      end
c=============================================================================
      subroutine integrt(y,yt,dy1,dy2,dy3,nit,s,ds,aion,zion,
     &            icharge,lezbeam,lperveance,lemittance,lallez,llinear,limage,
     &            lendzero,lfixed)
      integer(ISZ):: nit
      real(kind=8):: y(16,nit),yt(16,nit),dy1(16,nit),dy2(16,nit),dy3(16,nit)
      real(kind=8):: s,ds,aion,zion
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,lfixed

c INTEGRT advances the set of the envelope equations using a fourth-order
c Runge-Kutta method. The routine is adapted from the subroutine RK4 from
c Numerical Recipes by Press et al. and differs mainly in the use
c of two dimensional arrays rather than linear arrays.

      real(kind=8):: halfds,sixthds
      integer(ISZ):: i,j

      halfds = 0.5*ds
      sixthds = ds/6.0

c     call checkenvelope(y,nit)
      call derivs(s,ds,0.,y,dy1,nit,aion,zion,icharge,
     &            lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,
     &            lfixed)

      do i=1,10
        do j=1,nit
          yt(i,j) = y(i,j) + halfds*dy1(i,j)
        enddo
      enddo
      do i=11,15
        do j=1,nit
          yt(i,j) = y(i,j)
        enddo
      enddo
      call checkenvelope(yt,nit)

      call derivs(s,ds,0.5,yt,dy2,nit,aion,zion,icharge,
     &            lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,
     &            lfixed)

      do i=1,10
        do j=1,nit
          yt(i,j) = y(i,j) + halfds*dy2(i,j)
        enddo
      enddo
      call checkenvelope(yt,nit)

      call derivs(s,ds,0.5,yt,dy3,nit,aion,zion,icharge,
     &            lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,
     &            lfixed)

      do i=1,10
        do j=1,nit
          yt(i,j) = y(i,j) + ds*dy3(i,j)
          dy3(i,j) = dy2(i,j) + dy3(i,j)
        enddo
      enddo
      call checkenvelope(yt,nit)

      call derivs(s,ds,1.0,yt,dy2,nit,aion,zion,icharge,
     &            lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,
     &            lfixed)

      do i=1,10
        do j=1,nit
          y(i,j) = y(i,j) + sixthds*(dy1(i,j) + dy2(i,j) + 2.0*dy3(i,j))
        enddo
      enddo
      call checkenvelope(y,nit)
c     call checkenvelope(yt,nit)

      return
      end
c=============================================================================
      subroutine checkenvelope(y,nit)
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
c     real(kind=8):: temp
c Checks to make sure that the envelope is positive
      integer(ISZ):: i
      do i=1,nit
        if (y(1,i) <= 0.) then
          y(1,i) = -y(1,i)
          y(2,i) = -y(2,i)
        endif
        if (y(3,i) <= 0.) then
          y(3,i) = -y(3,i)
          y(4,i) = -y(4,i)
        endif
      enddo
c Cheap fix to make sure the time is increasing monotonically
c     do i=2,nit
c       if (y(9,i) < y(9,i-1)) then
c         temp = y(9,i)
c         y(9,i) = y(9,i-1)
c         y(9,i-1) = temp
c       endif
c     enddo
      return
      end
c=============================================================================
      subroutine circegetimage(y,nit,llinear,limage,lendzero)
      use Constant
      use CIRfield
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      logical(ISZ):: llinear,limage,lendzero

c This subroutine was renamed circegetimage, instead of getimage, since a
c subroutine getimage exists already in the Hermes package.
c CIRCEGETIMAGE calculates image-field coefficients as functions of
c time for the pipe type ipipetype.
c   0 round pipe (default)
c   1 electric quadrupole rods (hyperbolic sufaces)
c   2 dipole plates (parallel plates)
c Current values of a, b, X, and Y from 'y' are used.
c The assumed forms for the image-field components scaled by
c rigidity*beta*clight*gamma**2 are
c   Ex = (K/R**2)*(Fxx*(x-X) + Fxy*(y-Y) + Gxx*X + Gxy*Y)
c   Ey = (K/R**2)*(Fyx*(x-X) + Fyy*(y-Y) + Gyx*X + Gyy*Y).
c Fxy and Fyx do not to contribute to the averaged equations, and
c the remaining coefficients are stored in common block CIRtmp as...
c   fx(ii) = Fxx
c   fy(ii) = Fyy
c   gxx(ii) = Gxx
c   gxy(ii) = Gxy
c   gyx(ii) = Gyx
c   gyy(ii) = Gyy
c
c If linear == .true., terms of the order of (a,b,X,Y)**2/rpipe**2
c in the image coefficients are discarded.
c
c Image coefficients are zero if the centroid displacement
c X**2 + Y**2 exceeds rpipe**2.

      integer(ISZ):: ii
      real(kind=8):: logcon
      real(kind=8):: pi2,cosrot2,sinrot2,xval,yval,xyrpipe,abrpipe,displace,geomfac
      real(kind=8):: f1,f2,g1,g2
      real(kind=8):: dispx,dispy
      real(kind=8):: c0,c1,c2,d0,d1
      real(kind=8):: arg0,arg1,arg2
      logical(ISZ):: lzeroenv = .false.
      save lzeroenv

c     --- Check if images are to be applied. If not, zero arrays and return.
      if (.not. limage) then
        do ii=1,nit
          fx(ii) = 0.0
          fy(ii) = 0.0
          gxx(ii) = 0.0
          gxy(ii) = 0.0
          gyx(ii) = 0.0
          gyy(ii) = 0.0
        enddo
        return
      endif

c     --- Constants
      pi2 = pi**2
      cosrot2 = cosrot**2
      sinrot2 = sinrot**2

c     --- Image forces for dipoles (parallel plates)
      if (ipipetype == 2) then
        do ii=1,nit
          xval = y(5,ii) - xpipeerr
          yval = y(7,ii) - ypipeerr
          if (llinear) then
            xyrpipe = 0.
            abrpipe = 0.
          else
            xyrpipe = ((xval*cosrot + yval*sinrot)/rpipe)**2
            abrpipe = (y(1,ii)**2 - y(3,ii)**2)/rpipe**2
          endif
          f1 = 1.0 - 0.375*pi2*xyrpipe
          f2 = 0.025*pi2*abrpipe
          g1 = 1.0 - (pi2/12.0)*(xyrpipe + 0.1875*abrpipe
     .              *(cosrot2 - sinrot2))

          displace = (xval**2 + yval**2)/rpipe**2
          if (displace < 1.0) then
            geomfac = pi2/24.0
          else
            geomfac = 0.0
          endif
          fx(ii) = geomfac*cosrot2*(f1 - f2*(1.0 - 0.25*sinrot2))
          fy(ii) = geomfac*sinrot2*(f1 + f2*(1.0 - 0.25*cosrot2))
          gxx(ii) = 3.0*geomfac*cosrot2*g1
          gxy(ii) = 3.0*geomfac*sinrot*cosrot*g1
          gyx(ii) = 3.0*geomfac*sinrot*cosrot*g1
          gyy(ii) = 3.0*geomfac*sinrot2*g1
        enddo

      elseif (ipipetype == 1) then
c     --- Set image coefficients for electric quadrupoles
c     --- Laslett formulation is used here

        do ii=1,nit
          xval = y(5,ii) - xpipeerr
          yval = y(7,ii) - ypipeerr
          dispx = xval**2/rpipe**2
          dispy = yval**2/rpipe**2

          displace = dispx + dispy
          if (displace < 1.0) then
            geomfac = 1.0
          else
            geomfac = 0.0
          endif

          if ((y(1,ii) <= 0.0).or.(y(3,ii) <= 0.0)) then
            if (.not. (lzeroenv .or. lendzero)) then
              print*," fatal error: beam radius less than zero at slice ",ii
              call kaboom(-1)
              lzeroenv = .true.
            endif
          endif
          logcon = y(1,ii)*y(3,ii)*log10(y(1,ii)/y(3,ii))/rpipe**2
          c0 = 0.785 + 3.0*logcon
          d0 = 0.785 - 3.0*logcon
          c1 = 0.834 + 3.0*logcon
          d1 = 0.834 - 3.0*logcon
          c2 = - (2.486 + (6.5*logcon)**2)
          fx(ii) = 0.0
          fy(ii) = 0.0
          gxx(ii) = geomfac*(c0 + 2.0*c1*dispx + c2*dispy)
          gxy(ii) = 0.0
          gyx(ii) = 0.0
          gyy(ii) = geomfac*(d0 + 2.0*d1*dispy + c2*dispx)
        enddo

      else

c     --- set image coefficients for circular pipe
c     --- image-force geometric factors f1, g1, and g2 for circular pipe

        do ii=1,nit
          xval = y(5,ii) - xpipeerr
          yval = y(7,ii) - ypipeerr
          arg0 = (y(1,ii)**2-y(3,ii)**2)/(4.0*rpipe**2)
          arg1 = (xval**2 - yval**2)/rpipe**2
          arg2 = (xval**2 + yval**2)/rpipe**2
          if (llinear) then
            f1 = 0.
            g1 = 0.
            g2 = 0.
          else
            f1 = arg1 + arg0*(1.0 + 6.0*arg2 + 6.0*arg0*arg1)
            g1 = arg1 + arg0*(1.0 + 3.0*arg2 + 2.0*arg0*arg1)
            g2 = 2.0*(1.0-2.0*arg0**2)*xval*yval/rpipe**2
          endif

c   turn off image force if X**2+Y**2 >= rpipe**2

          displace = arg2
          if (displace < 1.0) then
            geomfac = 1.0
          else
            geomfac = 0.0
          endif
          fx(ii)  = geomfac*f1
          fy(ii)  = -geomfac*f1
          gxx(ii) = geomfac*(1.0 + g1)
          gxy(ii) = geomfac*g2
          gyx(ii) = geomfac*g2
          gyy(ii) = geomfac*(1.0 - g1)
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getezbeam(nit,y,eval,rpipe,icharge,lezbeam,lendzero,lfixed,
     &                     lfailed)
      use Constant
      use CIRtmp
      use CIRbessdata
      integer(ISZ):: nit
      real(kind=8):: y(16,nit),eval(nit),rpipe
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lendzero,lfixed,lfailed

c GETEZBEAM calculates the axial space-charge field of a beam
c using various approximations.
c beam data is passed in array y, and the space-charge field
c at slice boundaries is returned in array eval.
c The field is set to zero if lezbeam if false.
c The space-charge model is determined by icharge
c   icharge = 1 uses simple g-factor space-charge model
c   icharge = 2 includes envelope variation in space-charge model
c   icharge = 3 includes end effects in space-charge model
c   icharge = 4 uses a more accurate model of end effects
c   icharge = 5 uses a Bessel-series model of end effects
c The beam current at slice boundaries cur is also calculated here.
c For icharge = 3, a simple generalization of the g-factor model
c is used, rather than the more accurate version obtained by
c expansion about the integrand maximum.
c note: modified expressions should be derived for parallel-plate
c boundaries and for electrostatic quadrupoles.

      real(kind=8):: frpieps
      real(kind=8):: rmin,abmin,aval,bval,del1,del2,fact1,fact2,fact3,part,epsfact,zetaval
      real(kind=8):: endval,exfact,zetapeak,gval,zeron,denomn,zetafact
      real(kind=8):: part0,rad0,bet0,den0,dden0,denv0,arg0
      real(kind=8):: part1,rad1,bet1,den1,dden1,denv1,arg1
      integer(ISZ):: ii,istart,istop,ii0,ii1,iimid,n
      real(kind=8):: bessp0,bessp1,bessp2,con0,con1,bessfun0,bessfun1,bessfun2,bessfun3
      real(kind=8):: bess0,bess1,bess2
      real(kind=8):: bessj0,bessj1
      real(kind=8):: didt,dvdt
      logical(ISZ):: lzeroenv = .false.
      save lzeroenv

      iimid = max(1,(nit+1)/2)
      frpieps = 4.0*pi*eps0

c calculate space-charge geometric factor g
c   note: in this version, a and b are allowed to be negative

      if (icharge .eq. 0) then
        do ii=1,nit
          gtemp(ii) = 2*log(1.6)
        enddo
      else
        rmin = 0.01*rpipe
        abmin = rpipe
        do 10 ii=1,nit
          aval = y(1,ii)
          bval = y(3,ii)
          abmin = min(abmin,aval,bval)
          aval = min(max(rmin,aval),rpipe)
          bval = min(max(rmin,bval),rpipe)
          rad(ii) = sqrt(aval*bval)
          if(icharge < 2)then
            gtemp(ii) = log(rpipe**2/(aval*bval))
          else
            gtemp(ii) = log(rpipe**2/(aval*bval)) + 0.5
          endif
   10   continue

c exit if a or b is less than zero

        if(abmin < 0.0)then
          if (.not. lzeroenv) then
            print*," fatal error: beam radius is less than zero"
            call kaboom(-1)
            lzeroenv = .true.
          endif
        endif
      endif

c calculate current and line-charge density times c at slice midpoints

      curmid(nit)=0.0
      do 20 ii=1,nit-1
        curmid(ii) = y(14,ii)/(y(9,ii+1) - y(9,ii))
   20 continue

c calculate beam quantities on slice boundaries
c   cur = beam current cur
c   dlna = logarithmic time derivative of a
c   dlnb = logarithmic time derivative of b
c   den = line-charge density times c
c   denv = envelope-variation term 0.5*den*(d(ab)/dt)/(ab*clight)
c   dden = time-derivative of the line-charge density d(den)/clight*dt
c   fixed model:
c   dden = beta * z-derivative of the line charge density beta*clight*d(lambda)/dz
c          calculated from (beta/clight)*d(current/(beta**2))/dt

c   values at beam head

      del1 = y(9,2) - y(9,1)
      del2 = y(9,3) - y(9,2)
      fact1 = (2.0 + del2/del1)/(clight*(del1 + del2))
      fact2 = del1/(clight*del2*(del1 + del2))
c   note: dlna and dlnb are used as diagnostics here
      if (.not. lendzero) then
        dlna(1) = (fact1*(y(1,2) - y(1,1))
     .           - fact2*(y(1,3) - y(1,2)))/y(1,1)
        dlnb(1) = (fact1*(y(3,2) - y(3,1))
     .           - fact2*(y(3,3) - y(3,2)))/y(3,1)
      else
        dlna(1) = 0.
        dlnb(1) = 0.
      endif

      y(13,1) = 0.0
      y(15,1) = 0.0
      denv(1) = 0.0

c   values at beam tail

      del1 = y(9,nit) - y(9,nit-1)
      del2 = y(9,nit-1) - y(9,nit-2)
      fact1 = (2.0 + del2/del1)/(clight*(del1 + del2))
      fact2 = del1/(clight*del2*(del1 + del2))
c   note: dlna and dlnb are used as diagnostics here
      if (.not. lendzero) then
        dlna(nit) = (fact1*(y(1,nit) - y(1,nit-1))
     .             - fact2*(y(1,nit-1) - y(1,nit-2)))/y(1,nit)
        dlnb(nit) = (fact1*(y(3,nit) - y(3,nit-1))
     .             - fact2*(y(3,nit-1) - y(3,nit-2)))/y(3,nit)
      else
        dlna(nit) = 0.
        dlnb(nit) = 0.
      endif

      y(13,nit) = 0.0
      y(15,nit) = 0.0
      denv(nit) = 0.0

c   values in beam body

      do 30 ii=2,nit-1

        del1 = y(9,ii) - y(9,ii-1)
        del2 = y(9,ii+1) - y(9,ii)
        part = del1/(del1 + del2)
        fact1 = del1/(clight*del2*(del1 + del2))
        fact2 = del2/(clight*del1*(del1 + del2))

        if (del1 < 0.) then
          lfailed = .true.
        endif

        dlna(ii) = (fact1*(y(1,ii+1)-y(1,ii))
     .         + fact2*(y(1,ii)-y(1,ii-1)))/y(1,ii)
        dlnb(ii) = (fact1*(y(3,ii+1)-y(3,ii))
     .         + fact2*(y(3,ii)-y(3,ii-1)))/y(3,ii)

        y(13,ii) = (1.0-part)*curmid(ii-1) + part*curmid(ii)
        y(15,ii) = y(13,ii) / y(10,ii)
        denv(ii) = 0.5*y(15,ii)*(dlna(ii) + dlnb(ii))
   30 continue
      if (del2 < 0.) then
        lfailed = .true.
      endif
   
      if (lfixed) then
        do ii=2,nit-1
          del1 = y(9,ii) - y(9,ii-1)
          del2 = y(9,ii+1) - y(9,ii)
          fact1 = 2./(del1*(del1+del2))
          fact2 = 2./(del2*(del1+del2))
          didt = fact2*y(14,ii) - fact1*y(14,ii-1)
          fact1 = del1/(del2*(del1 + del2))
          fact2 = del2/(del1*(del1 + del2))
          fact3 = (del2 - del1)/(del1*del2)
          dvdt = fact1*y(10,ii+1) - fact2*y(10,ii-1) + fact3 * y(10,ii)
          dden(ii) = (didt - 2*y(13,ii)*dvdt/y(10,ii))/(y(10,ii)*clight)
        enddo
        del1 = y(9,2) - y(9,1)
        didt = 2*y(14,1)/del1**2
        dden(1) = didt/(y(10,1)*clight)
        del2 = y(9,nit) - y(9,nit-1)
        didt = -2*y(14,nit-1)/del2**2
        dden(nit) = didt/(y(10,nit)*clight)
      else
        denmid(nit)=0.0
        do ii=1,nit-1
          denmid(ii) = 2.0*curmid(ii)/(y(10,ii+1) + y(10,ii))
        enddo
        dden(1) = 2.0*denmid(1)/((y(9,2) - y(9,1))*clight)
        dden(nit) = -2.0*denmid(nit-1)/((y(9,nit) - y(9,nit-1))*clight)
        do ii=2,nit-1
          dden(ii) = 2.0*(denmid(ii) - denmid(ii-1)) /((y(9,ii+1) - y(9,ii-1))*clight)
        enddo
      endif

      if (lfixed) then
c ---   Add da/ds term in denv
        do ii=2,nit-1
          denv(ii) = denv(ii) - 0.5*y(10,ii)*y(2,ii)/y(1,ii)
          denv(ii) = denv(ii) - 0.5*y(10,ii)*y(4,ii)/y(3,ii)
        enddo
      endif

c calculate space-charge field at slice boundaries
c   eval = longitudinal space-charge field

      if((icharge < 0).or.(icharge>5).or.(.not. lezbeam))then
c     set space-charge field to zero

        do 50 ii=1,nit
          eval(ii) = 0.0
   50   continue

      elseif(icharge == 0 .or. icharge == 1)then
c       simple g-factor model

        do 60 ii=1,nit
          eval(ii) = gtemp(ii)*dden(ii)/(clight*frpieps*y(10,ii))
   60   continue

      elseif(icharge == 2)then
c       general g-factor model

        do 70 ii=1,nit
          epsfact = 1.0/(clight*frpieps*y(10,ii))
          eval(ii) = epsfact*(gtemp(ii)*dden(ii) - denv(ii))
   70   continue

      elseif(icharge >= 3)then
c       space-charge models with end effects

c       calculate distance zeta from nearest end

        zetaval = 0.0
        zeta(1) = zetaval
        do 80 ii=2,nit
          zetaval = zetaval + 0.5*clight*(y(9,ii) - y(9,ii-1)) *
     .                            (y(10,ii) + y(10,ii-1))
          zeta(ii) = zetaval
   80   continue

        do 85 ii = iimid,nit
          zeta(ii) = zeta(nit) - zeta(ii)
   85   continue

c       find boundaries on ends
c       --- Initially set start and stop to middle of the beam, since
c       --- in some cases, zeta < endval for all beam slices.

        endval = 8.0*rpipe/besszeros(1)
        istart = iimid
        istop  = iimid
        do 90 ii=2,nit
          if((zeta(ii) > endval) .and. (zeta(ii-1) < endval))then
            istart = ii
          elseif((zeta(ii) < endval) .and. (zeta(ii-1) > endval))then
            istop = ii-1
          endif
   90   continue

c       use simple g-factor model in beam midsection
        do 100 ii=istart,istop
          epsfact = 1.0/(clight*frpieps*y(10,ii))
          eval(ii) = epsfact*(gtemp(ii)*dden(ii) - denv(ii))
  100   continue

c calculate space-charge for models with end effects

        if(icharge == 3)then
c         simple end-effect model expanding about end

          arg0 = besszeros(1)/rpipe
          do 110 ii=1,istart-1
            epsfact = 1.0/(clight*frpieps*y(10,ii))
            zetafact = arg0*zeta(ii)
            exfact = exp(-zetafact)
            etemp(ii) = 1.0 - 0.5*(1.0 + zetafact)*exfact
            ftemp(ii) = 0.5*arg0*exfact
            eval(ii) = epsfact*(gtemp(ii)*(etemp(ii)*dden(ii)
     .                - ftemp(ii)*y(15,ii)*y(10,ii))  - etemp(ii)*denv(ii))
  110     continue
          do 120 ii=istop+1,nit
            epsfact = 1.0/(clight*frpieps*y(10,ii))
            zetafact = arg0*zeta(ii)
            exfact = exp(-zetafact)
            etemp(ii) = 1.0 - 0.5*(1.0 + zetafact)*exfact
            ftemp(ii) = 0.5*arg0*exfact
            eval(ii) = epsfact*(gtemp(ii)*(etemp(ii)*dden(ii)
     .                + ftemp(ii)*y(15,ii)*y(10,ii))  - etemp(ii)*denv(ii))
  120     continue

        elseif(icharge == 4)then
c         refined end-effect model expanding about integrand maximum

            rad0 = 0.0
            rad1 = 0.0
            arg0 = besszeros(1)/rpipe
            zetapeak = 1.0/arg0
            do 130 ii=1,iimid
              ii0 = ii
              ii1 = nit-ii+1
              if(zeta(ii0+1) > zetapeak)then
                if(rad0 == 0.0)then
                  part0 = (zetapeak-zeta(ii0))/(zeta(ii0+1)-zeta(ii0))
                  part1 = 1.0-part0
                  bet0 = part0*y(10,ii0+1) + part1*y(10,ii0)
                  rad0 = part0*rad(ii0+1) + part1*rad(ii0)
                  den0 = part0*y(15,ii0+1) + part1*y(15,ii0)
                  dden0 = part0*dden(ii0+1) + part1*dden(ii0)
                  denv0 = part0*denv(ii0+1) + part1*denv(ii0)
                endif
              endif
              if(zeta(ii1-1) > zetapeak)then
                if(rad1 == 0.0)then
                  part0 = (zetapeak-zeta(ii1))/(zeta(ii1-1)-zeta(ii1))
                  part1 = 1.0-part0
                  bet1 = part0*y(10,ii1-1) + part1*y(10,ii1)
                  rad1 = part0*rad(ii1-1) + part1*rad(ii1)
                  den1 = part0*y(15,ii1-1) + part1*y(15,ii1)
                  dden1 = part0*dden(ii1-1) + part1*dden(ii1)
                  denv1 = part0*denv(ii1-1) + part1*denv(ii1)
                endif
              endif
              if(rad0*rad1.ne.0.0)go to 140
  130       continue

  140     epsfact = 1.0/(clight*frpieps)
          do 150 ii=1,istart-1
            zetafact = arg0*zeta(ii)
            exfact = exp(-zetafact)
            etemp(ii) = 1.0 - 0.5*(1.0 + max(1.0,zetafact))*exfact
            ftemp(ii) = 0.5*arg0*exfact
            if(zeta(ii) < zetapeak)then
              gval = 0.5 + 2.0*log(rpipe/rad0)
              eval(ii) = epsfact*(gval*(etemp(ii)*dden0 + ftemp(ii)*den0*bet0)
     .                      - etemp(ii)*denv0)/bet0
            else
              eval(ii) = epsfact*(gtemp(ii)*(etemp(ii)*dden(ii)
     .            + ftemp(ii)*y(15,ii)*y(10,ii)) - etemp(ii)*denv(ii))/y(10,ii)
            endif
  150     continue
          do 160 ii=istop+1,nit
            zetafact = arg0*zeta(ii)
            exfact = exp(-zetafact)
            etemp(ii) = 1.0 - 0.5*(1.0 + max(1.0,zetafact))*exfact
            ftemp(ii) = 0.5*arg0*exfact
            if(zeta(ii) < zetapeak)then
              gval = 0.5 + 2.0*log(rpipe/rad0)
              eval(ii) = epsfact*(gval*(etemp(ii)*dden1 - ftemp(ii)*den1*bet1)
     .                      - etemp(ii)*denv1)/bet1
            else
              eval(ii) = epsfact*(gtemp(ii)*(etemp(ii)*dden(ii)
     .            - ftemp(ii)*y(15,ii)*y(10,ii)) - etemp(ii)*denv(ii))/y(10,ii)
            endif
  160     continue

        elseif(icharge >= 5)then
c         approximate Bessel-series model
         
          do 170 ii=1,nit
            sum0(ii) = 0.0
            sum1(ii) = 0.0
  170     continue

          do 210 n=1,nmax
            rad0 = 0.0
            rad1 = 0.0

            if(n < 10)then
              zeron = besszeros(n)
              denomn = bessdenom(n)
            else
              zeron = pi*(n - 0.25)
              denomn = 2.0*(n - 0.25)
            endif

            arg0 = zeron/rpipe
            zetapeak = 1.0/arg0
            do 180 ii=1,iimid
              ii0 = ii
              ii1 = nit-ii+1
              if(zeta(ii0+1) > zetapeak)then
                if(rad0 == 0.0)then
                  part0 = (zetapeak-zeta(ii0))/(zeta(ii0+1)-zeta(ii0))
                  part1 = 1.0-part0
                  bet0 = part0*y(10,ii0+1) + part1*y(10,ii0)
                  rad0 = part0*rad(ii0+1) + part1*rad(ii0)
                  den0 = part0*y(15,ii0+1) + part1*y(15,ii0)
                  dden0 = part0*dden(ii0+1) + part1*dden(ii0)
                  denv0 = part0*denv(ii0+1) + part1*denv(ii0)
                endif
              endif
              if(zeta(ii1-1) > zetapeak)then
                if(rad1 == 0.0)then
                  part0 = (zetapeak-zeta(ii1))/(zeta(ii1-1)-zeta(ii1))
                  part1 = 1.0-part0
                  bet1 = part0*y(10,ii1-1) + part1*y(10,ii1)
                  rad1 = part0*rad(ii1-1) + part1*rad(ii1)
                  den1 = part0*y(15,ii1-1) + part1*y(15,ii1)
                  dden1 = part0*dden(ii1-1) + part1*dden(ii1)
                  denv1 = part0*denv(ii1-1) + part1*denv(ii1)
                endif
              endif
  180       continue

            bessp0 = bessj0(arg0*rad0)
            bessp1 = bessj1(arg0*rad0)
            bessp2 = (2.0/(arg0*rad0))*bessp1 - bessp0
            do 190 ii=1,istart-1
              arg1 = arg0*zeta(ii)
              con1 = exp(-arg1)
              con0 = 1.0 - 0.5*(1.0 + max(1.0,arg1))*con1
              if(zeta(ii) <= zetapeak)then
                bessfun0 = 2.0*bessp1/(arg0*rad0*denomn)
                bessfun1 = 2.0*dden0*bessp1/(arg0*bet0*rad0)
                bessfun2 = 2.0*denv0*bessp2/bet0
                bessfun3 = den0*bessp1/rad0
              else
                bess0 = bessj0(arg0*rad(ii))
                bess1 = bessj1(arg0*rad(ii))
                bess2 = (2.0/(arg0*rad(ii)))*bess1 - bess0
                bessfun0 = 2.0*bess1/(arg0*rad(ii)*denomn)
                bessfun1 = 2.0*dden(ii)*bess1/
     .                             (arg0*y(10,ii)*rad(ii))
                bessfun2 = 2.0*denv(ii)*bess2/y(10,ii)
                bessfun3 = y(15,ii)*bess1/rad(ii)
              endif
              sum0(ii) = sum0(ii) + con0*bessfun0*(bessfun1 - bessfun2)
              sum1(ii) = sum1(ii) + con1*bessfun0*bessfun3

  190       continue

            bessp0 = bessj0(arg0*rad1)
            bessp1 = bessj1(arg0*rad1)
            bessp2 = (2.0/(arg0*rad1))*bessp1 - bessp0
            do 200 ii=istop+1,nit
              arg1 = arg0*zeta(ii)
              con1 = exp(-arg1)
              con0 = 1.0 - 0.5*(1.0 + max(1.0,arg1))*con1
              if(zeta(ii) <= zetapeak)then
                bessfun0 = 2.0*bessp1/(arg0*rad1*denomn)
                bessfun1 = 2.0*dden1*bessp1/(arg0*bet1*rad1)
                bessfun2 = 2.0*denv1*bessp2/bet1
                bessfun3 = den1*bessp1/rad1
              else
                bess0 = bessj0(arg0*rad(ii))
                bess1 = bessj1(arg0*rad(ii))
                bess2 = (2.0/(arg0*rad(ii)))*bess1 - bess0
                bessfun0 = 2.0*bess1/(arg0*rad(ii)*denomn)
                bessfun1 = 2.0*dden(ii)*bess1/
     .                             (arg0*y(10,ii)*rad(ii))
                bessfun2 = 2.0*denv(ii)*bess2/y(10,ii)
                bessfun3 = y(15,ii)*bess1/rad(ii)
              endif
              sum0(ii) = sum0(ii) + con0*bessfun0*(bessfun1 - bessfun2)
              sum1(ii) = sum1(ii) - con1*bessfun0*bessfun3
  200       continue
  210     continue

          con0 = 4.0/(clight*frpieps)
          do 220 ii=1,istart-1
            eval(ii) = con0*(sum0(ii) + sum1(ii))
  220     continue
          do 240 ii=istop+1,nit
            eval(ii) = con0*(sum0(ii) + sum1(ii))
  240     continue

        endif
      endif

      return
      end
c=============================================================================
      real(kind=8) function bessj0(x)
      real(kind=8):: x

c     BESSJ0 calculates the Bessel function J0 using matched analytic
c     approximations.  The routine is taken from Numerical Recipes.

      real(kind=8):: p1 =  1.0
      real(kind=8):: p2 = -0.1098628627e-02
      real(kind=8):: p3 =  0.2734510407e-04
      real(kind=8):: p4 = -0.2073370639e-05
      real(kind=8):: p5 =  0.2093887211e-06
      real(kind=8):: q1 = -0.1562499995e-01
      real(kind=8):: q2 =  0.1430488765e-03
      real(kind=8):: q3 = -0.6911147651e-05
      real(kind=8):: q4 =  0.7621095161e-06
      real(kind=8):: q5 = -0.934945152e-07
      real(kind=8):: r1 =  57568490574.0
      real(kind=8):: r2 = -13362590354.0
      real(kind=8):: r3 =  651619640.7
      real(kind=8):: r4 = -11214424.18
      real(kind=8):: r5 =  77392.33017
      real(kind=8):: r6 = -184.9052456
      real(kind=8):: s1 =  57568490411.
      real(kind=8):: s2 =  1029532985.
      real(kind=8):: s3 =  9494680.718
      real(kind=8):: s4 =  59272.64853
      real(kind=8):: s5 =  267.8532712
      real(kind=8):: s6 =  1.0

      real(kind=8):: ax,y,z,xx

      ax = abs(x)
      if(ax < 8.0)then
        y = x**2
        bessj0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     .          /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        z = 8.0/ax
        y = z**2
        xx = ax - 0.785398164
        bessj0 = sqrt(0.636619772/ax)*
     .             (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))
     .           - z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif

      return
      end
c=============================================================================
      real(kind=8) function bessj1(x)
      real(kind=8):: x

c     BESSJ1 calculates the Bessel function J0 using matched analytic
c     approximations.  The routine is taken from Numerical Recipes.

      real(kind=8):: r1 = 72362614232.0
      real(kind=8):: r2 = -7895059235.0
      real(kind=8):: r3 = 242396853.1
      real(kind=8):: r4 = -2972611.439
      real(kind=8):: r5 = 15704.48260
      real(kind=8):: r6 = -30.16036606
      real(kind=8):: s1 = 144725228442.0
      real(kind=8):: s2 = 2300535178.0
      real(kind=8):: s3 = 18583304.74
      real(kind=8):: s4 = 99447.43394
      real(kind=8):: s5 = 376.9991397
      real(kind=8):: s6 = 1.0
      real(kind=8):: p1 = 1.0
      real(kind=8):: p2 = 0.183105e-02
      real(kind=8):: p3 = -0.3516396496e-04
      real(kind=8):: p4 = 0.2457520174e-05
      real(kind=8):: p5 = -.240337019e-06
      real(kind=8):: q1 = 0.04687499995
      real(kind=8):: q2 = -0.2002690873e-03
      real(kind=8):: q3 = 0.8449199096e-05
      real(kind=8):: q4 = -0.88228987e-06
      real(kind=8):: q5 = 0.105787412e-06

      real(kind=8):: ax,y,z,xx

      ax = abs(x)
      if(ax < 8.0)then
        y = x**2
        bessj1 = x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     .            /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else

        z = 8./ax
        y = z**2
        xx = ax - 2.356194491
        bessj1 = sqrt(0.636619772/ax)*
     .             (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))
     .             -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     .             *sign(1.0,x)
      endif

      return
      end
c=============================================================================
      subroutine savehistcir(s,y,nit,nscir,lsavehist)
      use CIRhist
      integer(ISZ):: nit,nscir
      real(kind=8):: s,y(16,nit)
      logical(ISZ):: lsavehist

c     --- Check if saving history
      if (lsavehist) then
c       --- Check if saving history this step. Requires (nhcir > 0).
        if (mod(nscir,nhcir) == 0) then
          if (jhcir >= lhcir) then
            lhcir = lhcir + 100
            if (nithist == 0) nithist = nit
            call gchange("CIRhist",0)
          endif
          call savehistcir1(s,y,nit)
        endif
      endif

      return
      end
c=============================================================================
      subroutine savehistcir1(s,y,nit)
      use CIRhist
      integer(ISZ):: nit
      real(kind=8):: s,y(16,nit)

      integer(ISZ):: i,ii
      real(kind=8):: ri,ristep

      jhcir = jhcir + 1
      hscir(jhcir) = s

c     --- ri and ristep are intentionally real(kind=8):: numbersa. This gives an
c     --- evenly spaced list of slices to save the data of, including the
c     --- end points, as well as the middle if nithist is odd.
      ristep = (nit-1.)/(nithist-1.)
      ri = 1.
      do i=1,nithist
        ii = int(ri + 0.5)
        ri = ri + ristep
        hacir(i,jhcir)   = y( 1,ii)
        hapcir(i,jhcir)  = y( 2,ii)
        hbcir(i,jhcir)   = y( 3,ii)
        hbpcir(i,jhcir)  = y( 4,ii)
        hxcir(i,jhcir)   = y( 5,ii)
        hxpcir(i,jhcir)  = y( 6,ii)
        hycir(i,jhcir)   = y( 7,ii)
        hypcir(i,jhcir)  = y( 8,ii)
        htcir(i,jhcir)   = y( 9,ii)
        hvzcir(i,jhcir)  = y(10,ii)
        henxcir(i,jhcir) = y(11,ii)
        henycir(i,jhcir) = y(12,ii)
        hcur(i,jhcir)    = y(13,ii)
        hdq(i,jhcir)     = y(14,ii)
        hden(i,jhcir)    = y(15,ii)
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine getphase(sig0,sig1,hlp,betaval,curval,tval,sval,emitx,emity,
     &                    aion,zion,lperveance)
      real(kind=8):: sig0,sig1,hlp,betaval,curval,tval,sval,emitx,emity,aion,zion
      logical(ISZ):: lperveance

c Calls getphaset for one slice.

      real(kind=8):: tsig0(1),tsig1(1),ty(16)

      ty( 9) = tval
      ty(10) = betaval
      ty(11) = emitx
      ty(12) = emity
      ty(13) = curval
      call getphaset(1,tsig0,tsig1,hlp,sval,ty,aion,zion,lperveance)
      sig0 = tsig0(1)
      sig1 = tsig1(1)

      return
      end
c=============================================================================
      subroutine getphaset(nit,sig0,sig1,hlp,sval,y,aion,zion,lperveance)
      use Constant
      use Lattice
      use LatticeInternal
      integer(ISZ):: nit
      real(kind=8):: sig0(nit),sig1(nit),hlp,y(16,nit),sval,aion,zion
      logical(ISZ):: lperveance

c     GETPHASE estimates the undepressed phase advance sig0 in radians
c     per lattice period and the depressed phase advance sig1 using
c     the approximate expressions by Lee, Fessenden, and Laslett.
c     betaval is the value of beta at lattice position sval.
c     curval is the beam current used to calculate the perveance.
c     tval is the time the middle beam slice arrives at sval.
c     The average radius is either calculated from the geometric mean
c     of the midpoint beam radii or estimated from the perveance and
c     normalized emittances emitnx0 and emitny0.
c     The principle omissions from this formulation are the absence
c     of the effects of bend magnets and the presumption of simple
c     periodic focusing.

      integer(ISZ):: iq1,iq2,ii
      real(kind=8):: offset
      real(kind=8):: gambeta,rigidity,eta,fstr,focus,arg
      real(kind=8):: betafreq,emitavg,radavg,perv
      real(kind=8):: alfcuri
c     real(kind=8):: sig2

c Print warning if there are not enough quads to make a lattice period

      if (nquad == 0) then
        call remark("warning: full period not found, phase advance unset")
        return
      endif

c --- Use slice at middle of beam
      ii = (nit+1)/2

c Set internal lattice parameters
      call circesetlatt(sval,0.)

c calculate quadrupole occupancy in nearest lattice period

      if (nquad > 0) then
        offset = 0.
        iq1 = minval(cquadid(0,:))
        iq2 = iq1 + 1
        if (iq2 > nquad) then
          if (zlatperi > 0.) then
            iq2 = 0
            offset = zlatperi
          else
            iq2 = iq2 - 1
            call remark("Warning: full lattice period not found for quads")
          endif
        endif

        hlp = (quadzs(iq2) + offset - quadzs(iq1))
        eta = (quadze(iq1) - quadzs(iq1) + quadze(iq2) - quadzs(iq2))/(2.*hlp)
        fstr = (abs(quadde(iq1))/(y(10,ii)*clight) + abs(quaddb(iq1)))*
     &          (quadze(iq1) - quadzs(iq1)) +
     &         (abs(quadde(iq2))/(y(10,ii)*clight) + abs(quaddb(iq2)))*
     &          (quadze(iq2) - quadzs(iq2))

      elseif (nhele > 0) then
        offset = 0.
        iq1 = minval(cheleid(0,:))
        iq2 = iq1 + 1
        if (iq2 > nhele) then
          if (zlatperi > 0.) then
            iq2 = 0
            offset = zlatperi
          else
            iq2 = iq2 - 1
            call remark("Warning: full lattice period not found for heles")
          endif
        endif
        do ii = 1, max(helenm(iq1),helene(iq1)) 
          if ((nint(hele_n(ii,iq1)) == 2) .and. 
     &        (nint(hele_v(ii,iq1)) == 0)) then
            hlp = (helezs(iq2) + offset - helezs(iq1))
            eta = (heleze(iq1)-helezs(iq1) + heleze(iq2)-helezs(iq2))/(2.*hlp)
            fstr = (abs(heleae(ii,iq1))/(y(10,ii)*clight)+abs(heleam(ii,iq1)))*
     &              (heleze(iq1) - helezs(iq1)) +
     &             (abs(heleae(ii,iq2))/(y(10,ii)*clight)+abs(heleam(ii,iq2)))*
     &              (heleze(iq2) - helezs(iq2))
          endif
        enddo

      endif

      alfcuri = zion*echarge/(4.*pi*eps0*amu*aion*clight**3)
      do ii=1,nit

        gambeta = y(10,ii)/sqrt(1.0-y(10,ii)**2)
        rigidity = gambeta/(zion*echarge/(aion*amu*clight))

        focus = fstr/(eta*rigidity)

c       --- calculate undepressed phase advance sig0 in radians
        if(eta.ne.0.0)then
          arg = ((1.5-eta)/3.0)*(eta*focus*hlp**2)**2
          arg = max(0.0,min(2.0,arg))
          sig0(ii) = max(SMALLPOS,acos(1.0 - arg))
        else
          sig0(ii) = SMALLPOS
        endif

        if(sig0(ii)<=SMALLPOS) call remark("warning: unphysical phase-advance")

c       --- calculate depressed phase advance

        if(y(13,ii) == 0.0)then
c         --- set value to sig0 if curval = 0
          sig1(ii) = sig0(ii)

        else
c         --- calculate value from average equilibrium radius
c         --- note: an envelope-equation approximation is used here rather
c         ---       than the dubious expression from Lee, Fessenden & Laslett

c         --- perveance perv (labeled K in the envelope equation)
          if (lperveance) then
            perv = 2.0*y(13,ii)*(alfcuri/gambeta**3)
          else
            perv = 0.
          endif

          betafreq = 0.5*sig0(ii)/hlp
          emitavg = 0.5*(y(11,ii) + y(12,ii))
          radavg = sqrt(0.5*(perv + sqrt(perv**2
     .                    + 4.0*(betafreq*emitavg)**2))/betafreq**2)

          sig1(ii) = 2.0*hlp*emitavg/radavg**2
c         sig2 = sqrt(max(0.0,sig0(ii)**2 - 4.0*perv*(hlp/radavg)**2))

        endif

      enddo

      return
      end
c=============================================================================
      subroutine setbeam(nit,y,a0,b0,ap0,bp0,x0,y0,xp0,yp0,
     &                   ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot,
     &                   beta,gammabar,tilt,tiltmid,tiltlen,emitx,emity,
     &                   time,zbeam,aion,zion,
     &                   icurload,iprofile,ivelload,igrid,iemload,lperveance,
     &                   lendzero)
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: a0,b0,ap0,bp0,x0,y0,xp0,yp0
      real(kind=8):: ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot
      real(kind=8):: beta,gammabar,tilt,tiltmid,tiltlen
      real(kind=8):: emitx,emity,time,zbeam,aion,zion
      integer(ISZ):: icurload,iprofile,ivelload,igrid,iemload
      logical(ISZ):: lperveance,lendzero
c Sets up beam data based on input arguments
      call settval(nit,y,deltbeam,time,igrid)
      call setenvelope(nit,y,a0,b0,ap0,bp0,x0,y0,xp0,yp0,lendzero)
      call setcur(nit,y,ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot,
     &            icurload,iprofile,beta,tilt,tiltmid,tiltlen,ivelload)
      call setemit(nit,y,emitx,emity,beta,gammabar,ibeam,time,zbeam,aion,zion,
     &             iemload,lperveance)
      call checkenvelope(y,nit)

      return
      end
c=============================================================================
      subroutine setenvelope(nit,y,a0,b0,ap0,bp0,x0,y0,xp0,yp0,lendzero)
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: a0,b0,ap0,bp0,x0,y0,xp0,yp0
      logical(ISZ):: lendzero

c Sets envelope size and centroid

      integer(ISZ):: ii

      do ii=1,nit
        y(1,ii) = a0
        y(2,ii) = ap0
        y(3,ii) = b0
        y(4,ii) = bp0
        y(5,ii) = x0
        y(6,ii) = xp0
        y(7,ii) = y0
        y(8,ii) = yp0
      enddo
      if (lendzero) then
        y(1,1) = 0.
        y(2,1) = 0.
        y(3,1) = 0.
        y(4,1) = 0.
        y(5,1) = 0.
        y(6,1) = 0.
        y(7,1) = 0.
        y(8,1) = 0.
        y(1,nit) = 0.
        y(2,nit) = 0.
        y(3,nit) = 0.
        y(4,nit) = 0.
        y(5,nit) = 0.
        y(6,nit) = 0.
        y(7,nit) = 0.
        y(8,nit) = 0.
      endif


      return
      end
c=============================================================================
      subroutine setcur(nit,y,
     &                  ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot,
     &                  icurload,iprofile,
     &                  beta,tilt,tiltmid,tiltlen,
     &                  ivelload)
      use Constant
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: ibeam,currise,curfall,deltbeam,footrise,tfoot,curfoot
      integer(ISZ):: icurload,iprofile
      real(kind=8):: beta,tilt,tiltmid,tiltlen
      integer(ISZ):: ivelload

c     SETCUR calculates the initial velocity and current profiles as 
c     functions of slice arrival time t.  
c     Beta and current values are calculated at nit time values
c     given in array y(9,) and stored in arrays y(10,) (beta) and y(13,) (cur)
c     The line-charge profile is controlled by flag iprofile,
c     The flag icurload controls whether the line charge
c     is set up as a function of t or axial distance s.

      integer(ISZ):: ii
      real(kind=8):: valmax,part,arg1,arg2,trise,tfall,expfact,offset
      real(kind=8):: head,tail,tfootrise,tflat,sbeam,srise,sfall,sval

c calculate slice velocities

      call setbeta(beta,nit,y,tilt,tiltmid,tiltlen,ivelload)

c start current calculation

      valmax = 0.0

      if(icurload == 0)then

c set current as a function of t

        if(iprofile <= 0)then
          do ii=1,nit
            part = (y(9,ii) - y(9,1))/(y(9,nit) - y(9,1))
            arg1 = part/(currise+0.00001)
            arg2 = (1.0001-part)/(curfall+0.00001)
            y(13,ii) = y(10,ii)*tanh(arg1)*tanh(arg2)
            valmax = max(valmax,y(13,ii))
          enddo
        elseif(iprofile < 4)then
          trise = currise*deltbeam
          tfall = curfall*deltbeam
          do 24 ii=1,nit
            if(y(9,ii)-y(9,1) < trise)then
              y(13,ii) = 1.0 - ((trise-y(9,ii)+y(9,1))/trise)**iprofile
            elseif(y(9,ii)-y(9,1) > deltbeam-tfall)then
              y(13,ii) = 1.0-((y(9,ii)-y(9,1)-deltbeam+tfall)/tfall)**iprofile
            else
              y(13,ii) = 1.0
            endif
            y(13,ii) = y(10,ii)*y(13,ii)
            valmax = max(valmax,y(13,ii))
   24     continue
        elseif(iprofile == 4)then
c       Gaussian fall-off at ends
          trise = currise*deltbeam
          tfall = curfall*deltbeam
          expfact = 2.277
          offset = exp(-expfact)
          do 26 ii=1,nit
            if(y(9,ii)-y(9,1) < trise)then
              y(13,ii) = exp(-expfact*((trise-y(9,ii)+y(9,1))/trise)**2)-offset
            elseif(y(9,ii)-y(9,1) > deltbeam-tfall)then
              y(13,ii) = exp(-expfact*
     .                   ((y(9,ii)-y(9,1)-deltbeam+tfall)/tfall)**2) - offset
            else
              y(13,ii) = 1.0 - offset
            endif
            y(13,ii) = y(10,ii)*y(13,ii)
            valmax = max(valmax,y(13,ii))
   26     continue            
        elseif (iprofile == 5)then
c       Fermi-function ends
          trise = currise*deltbeam
          tfall = curfall*deltbeam
          expfact = 9.2
          do 28 ii=1,nit
            head = exp(expfact*((0.5*trise-y(9,ii)+y(9,1))/trise))
            tail = exp(expfact*((y(9,ii)-y(9,1)-deltbeam+0.5*tfall)/tfall))
            y(13,ii) = y(10,ii)/((1.0 + head)*(1.0 + tail))
            valmax = max(valmax,y(13,ii))
   28     continue
        elseif(iprofile == 6)then
c       main pulse with low-current foot
          trise = currise*deltbeam
          tfall = curfall*deltbeam
          tfootrise = footrise*deltbeam
          tflat = deltbeam - tfootrise - tfoot - tfall
          expfact = 3.0
          offset = exp(-expfact) 
          do 30 ii=1,nit
            if(y(9,ii)-y(9,1) < tfootrise)then
              y(13,ii) = curfoot*(1.0-((tfootrise-y(9,ii)+y(9,1))/tfootrise)**2)
            elseif(y(9,ii)-y(9,1) < (tfootrise + tfoot))then
              y(13,ii) = curfoot + (ibeam - curfoot)*
     .          exp(-expfact*((tfootrise + tfoot - y(9,ii)+y(9,1))/trise)**2)
            elseif(y(9,ii)-y(9,1) < (deltbeam - tfall))then
              y(13,ii) = ibeam
            else
              y(13,ii) = ibeam*
     .           exp(-expfact*((y(9,ii)-y(9,1) - (deltbeam-tfall))/tfall)**2)
     .                       - offset 
            endif
            valmax = max(valmax,y(13,ii))
   30     continue
        endif

      else 
        
c set current as a function of s
c note: other iprofile options should be added 

        sbeam = clight*y(10,nit)*(y(9,nit)-y(9,1))

        if(iprofile <= 0)then
c       hyperbollic-tangent profile

          do 50 ii=1,nit
            part = clight*y(10,ii)*(y(9,ii)-y(9,1))/sbeam
            arg1 = part/(currise+0.00001)
            arg2 = (1.0001-part)/(curfall+0.00001)
            y(13,ii) = (y(10,ii)**2)*tanh(arg1)*tanh(arg2)
            valmax = max(valmax,y(13,ii))
   50     continue

        else
c       power-law fall off

          srise = currise*sbeam
          sfall = curfall*sbeam
          do 60 ii=1,nit
            sval = clight*y(10,ii)*(y(9,ii)-y(9,1))
            if(sval < srise)then
              y(13,ii) = 1.0 - ((srise-sval)/srise)**iprofile
            elseif(sval > sbeam-sfall)then
              y(13,ii) = 1.0 - ((sval-sbeam+sfall)/sfall)**iprofile
            else
              y(13,ii) = 1.0
            endif
            y(13,ii) = (y(10,ii)**2)*y(13,ii)
            valmax = max(valmax,y(13,ii))
   60     continue

        endif
      endif

      do ii=1,nit
        y(13,ii) = ibeam*y(13,ii)/valmax
c       --- Set density from current and beta
        y(15,ii) = y(13,ii)/y(10,ii)
      enddo

      do ii=1,nit-1
c       --- set charge per slice dq
        y(14,ii) = 0.5*(y(13,ii+1)+y(13,ii))*(y(9,ii+1)-y(9,ii))
      enddo
      y(14,ii) = 0.

      return
      end
c=========================================================================
      subroutine setbeta(beta,nit,y,tilt,tiltmid,tiltlen,ivelload)
      real(kind=8):: beta
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: tilt,tiltmid,tiltlen
      integer(ISZ):: ivelload

c     SETBETA sets initial beta values at nit arrival times given in
c     array tval and stores them in array betaval.  
c     The average initial beta and head-to-tail relative velocity tilt
c     are given by beta and tile in common block sliceset.
c     The type of velocity tilt is determined by the flag ivelload
c     in common block cntrl.
c     This version uses a hyperbolic-tangent form for the velocity
c     variation, which reverts to a linear variation when the scale
c     length tiltlen is much larger than unity and the tilt midpoint
c     tiltmid is 0.5.  The 

      integer(ISZ):: ii
      real(kind=8):: duration,tmid,tscale,c0

c calculate beam midpoint

      duration = y(9,nit) - y(9,1)
      tmid = tiltmid*duration
      tscale = tiltlen*duration

c select type of velocity tilt

      if (ivelload == 0) then

c   impose linear tile in v-t space

        c0 = tilt/(tanh((y(9,nit)-tmid)/tscale) - tanh((y(9,1)-tmid)/tscale))
     
        do ii=1,nit
          y(10,ii) = beta*(1.0 + c0*tanh((y(9,ii) - tmid)/tscale))
        enddo

      else

c   impose linear tilt in v-s space
c   note: this algorithm still uses a linear tilt

        c0 =(exp(-tilt)-1.0)/duration
     
        do ii=1,nit
          y(10,ii) = beta/(1.0 + c0*(y(9,ii) - tmid))
        enddo

      endif

      return
      end
c===========================================================================
      subroutine settval(nit,y,deltbeam,time0,igrid)
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: deltbeam,time0
      integer(ISZ):: igrid

c Set initial arrival times

      integer(ISZ):: ii
      real(kind=8):: con0

      if(igrid <= 0)then
c       --- uniform spacing
        con0 = deltbeam/(nit-1)
        do ii=1,nit
          y(9,ii) = con0*(ii-1) + time0
        enddo
      else
c       --- cell size decreasing near ends
        con0 = deltbeam/(nit - 9)
        y(9,1) = 0.0 + time0
        y(9,2) = 0.1*con0 + time0
        y(9,3) = 0.3*con0 + time0
        y(9,4) = 0.6*con0 + time0
        y(9,5) = con0 + time0
        y(9,6) = 1.5*con0 + time0
        y(9,7) = 2.0*con0 + time0
        do ii= 8,nit-6
          y(9,ii) = y(9,ii-1) + con0
        enddo
        y(9,nit-5) = y(9,nit-6) + 0.5*con0
        y(9,nit-4) = y(9,nit-5) + 0.5*con0
        y(9,nit-3) = y(9,nit-4) + 0.4*con0
        y(9,nit-2) = y(9,nit-3) + 0.3*con0
        y(9,nit-1) = y(9,nit-2) + 0.2*con0
        y(9,nit  ) = y(9,nit-1) + 0.1*con0
      endif

      return
      end
c===========================================================================
      subroutine setemit(nit,y,emitx,emity,beta,gammabar,ibeam,time,zbeam,
     &                   aion,zion,iemload,lperveance)
      use Constant
      use Lattice
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: emitx,emity,beta,gammabar,ibeam,time,zbeam,aion,zion
      integer(ISZ):: iemload
      logical(ISZ):: lperveance

c     SETEMIT sets the initial beam emittance depending on the flag
c     iemload. Allowed options are
c       iemload = 0 constant normalized emittance
c       iemload = 1 constant transverse temperature
c       iemload = 2 constant beam charge-dentsity
c       iemload = 3 constant beam radius
c     Expressions used here are based on a smooth-focusing model. 

      integer(ISZ):: k,ii
      real(kind=8):: radcon,ratio,dencon,emitcon,pervcon,perv0,tempcon,gamcon
      real(kind=8):: sig0,sig1,hlp,s
      real(kind=8):: emitnx0,emitny0,alfcuri

c determine initial type of focusing

      if (nquad > 0) then
        if (quadde(0) .ne. 0) then
          k = 1
        else
          k = 2
        endif
      else
        k=0
      endif

c     --- Calculate normalized emittance
      emitnx0 = emitx*beta*gammabar
      emitny0 = emity*beta*gammabar

c set initial beam emittance

      if (iemload == 0) then

c       constant normalized emittance
        do ii=1,nit
          y(11,ii) = emitnx0
          y(12,ii) = emitny0
        enddo

      else

        call getphase(sig0,sig1,hlp,beta,ibeam,time,zbeam,emitx,emity,
     &                aion,zion,lperveance)

        emitcon = 0.5*(emitnx0**2 + emitny0**2)/(gammabar*beta)**2
        if (lperveance) then
          alfcuri = zion*echarge/(4.*pi*eps0*amu*aion*clight**3)
          pervcon = 2.0*alfcuri/gammabar**3
        else
          pervcon = 0.
        endif
        perv0 = pervcon*ibeam/beta**3
        radcon = 2.0*((hlp/sig0)**2)*(perv0 + 
     .             sqrt(perv0**2 + emitcon*(sig0/hlp)**2)) 
     
        if (iemload == 1) then

c       constant transverse temperature
          tempcon = emitcon/radcon
          gamcon = 2.0*gammabar*hlp/beta
          do ii=1,nit
            if (k == 1) then
              s = sig0*(beta/y(10,ii))**2
            else
              s = sig0*(beta/y(10,ii))
            endif
            ratio = (gamcon/s)*sqrt((tempcon*y(10,ii)**2
     .              + (pervcon*y(13,ii)/y(10,ii)))/radcon)
            y(11,ii) = ratio*emitnx0
            y(12,ii) = ratio*emitny0
          enddo

        elseif (iemload == 2) then

c       constant beam density
          dencon = ibeam/(radcon*beta*clight)
          do ii=1,nit
            if (k == 1) then
              s = sig0*(beta/y(10,ii))**2
            else
              s = sig0*(beta/y(10,ii))
            endif
            ratio = (y(13,ii)/(beta*clight*y(10,ii)))*
     .              sqrt((((0.5*y(10,ii)*s/hlp)**2/dencon)
     .               - pervcon*clight)/(emitcon*dencon))
            ratio = max(0.005,ratio)
            y(11,ii) = ratio*emitnx0
            y(12,ii) = ratio*emitny0
          enddo

        else

c       constant beam radius
          do ii=1,nit
            if (k == 1) then
              s = sig0*(beta/y(10,ii))**2
            else
              s = sig0*(beta/y(10,ii))
            endif
            ratio = sqrt((radcon*(0.5*y(10,ii)*s/hlp)**2
     .              - pervcon*y(13,ii)/y(10,ii))*radcon/emitcon)/beta
            y(11,ii) = ratio*emitnx0
            y(12,ii) = ratio*emitny0
          enddo
        endif
      endif
   
      return
      end
c===========================================================================
      subroutine setears(nit,y,ntval,earval,tval,z,iearset,iearform,icharge,
     &                   ieargap,lendzero,lfixed)
      use Constant
      use Beam_acc
      use Lattice
      use LatticeInternal
      use CIRtmp
      use CIRfield
      integer(ISZ):: nit,ntval
      real(kind=8):: y(16,nit),earval(ntval), tval(ntval),z
      integer(ISZ):: iearset,iearform,icharge,ieargap
      logical(ISZ):: lendzero, lfixed

c     SETEARS calculates ear-voltage values at nit times specified 
c     in array tval are stored in array earval.  Times are specified 
c     relative to the beam midpoint.
c     Two methods are used
c     If iearset = 0, a single set of ear-field values is calculated
c        based on the initial current profile and velocity.  
c     If iearset = 1, a set of ear-field values appropriate for the
c        ear cell nearest to the value z is calculated.
c     If iearset = 2, 

c     Modeling choices are set by iearform and icharge
c     If iearform = 0, a g-factor model is used 
c        If icharge = 1, ear fields are modeled by a simple expression
c           with the field proportional to the time derivative of the 
c           line-charge density current/beta.
c        If icharge > 1, a more general expression is used that 
c           includes the variation of the beam envelope.
c     if iearform = 1, the ear voltage is calculated from the beam
c        space-charge field, using the model specified by icharge.

c     The characteristic length in ear calculations is set by ieargap
c     If ieargap = 0, an appropriate average field is calculated for
c        the given number of ear cells, and the residence factor
c        L/Lgap is not changed during compression.
c     If ieargap = 1, the ear magnitude is based on the distance to 
c        next ear cell to account for varying spacing of ear gaps.
c     If ieargap = 2, the ear magnitude is based on the 
c        average of the distance from the last ear cell and that to 
c        next to account for varying spacing of ear cells.

      integer(ISZ):: inextear,ii,i,iear
      real(kind=8):: sear,slastear,snextear,amin,bmin,del0,del1,del2
      real(kind=8):: fact1,fact2,epsfact,part0,part1,aval,bval,betaval
      real(kind=8):: ddenval,denval,gval
      real(kind=8):: offset
      logical(ISZ):: lfailed
      real(kind=8):: eargap = -1.0
      logical(ISZ):: printnow = .true.
      save printnow
  
c     --- Update internal lattice and get a value for rpipe.
      if (iearset > 0 .or. ieargap > 0 .or. iearform > 0) then
        call extebcir(z,SMALLPOS,0.,y,nit)
      endif

c set spacing eargap between ear cells 

      if(eargap < 0.0)then
        if (naccl == 0) then
          eargap = 0.0
          call remark(" warning: no ear cells. ears cannot be set") 
          return
        endif
      elseif(eargap == 0.0)then
        return
      endif 
 
      if((ieargap == 0).or.(iearset == 0))then
c     calculate average value
        if (zlatperi > 0.) then
          eargap = zlatperi/naccl
        else
          eargap = (acclze(naccl) - acclzs(0))/naccl
        endif

      elseif(ieargap >= 1)then 
c     calculate local value based on ear cell spacing

        iear = minval(cacclid(0,:))
        sear = 0.5*(acclzs(iear) + acclze(iear)) + zlatstrt
        if (iear == 0) then
          slastear = zlatstrt
        else
          slastear = 0.5*(acclzs(iear-1) + acclze(iear-1)) + zlatstrt
        endif
       
c       find midpoint of next ear cell   
        inextear = iear + 1
        offset = zlatstrt
        if (inextear > naccl) then
          if (zlatperi > 0.) then
            inextear = 0
            offset = zlatperi + zlatstrt
          else
            inextear = iear
          endif
        endif
        snextear = 0.5*(acclzs(inextear) + acclze(inextear)) + offset

c       calculate distance eargap for ear calculation
c       --- When current gap is last gap, use distance between last two
c       --- gaps as a guess.
        if(ieargap == 1 .and. inextear > iear)then
          eargap = snextear - sear
        elseif (inextear > iear) then
          eargap = 0.5*(snextear - slastear)
        else
          eargap = snextear - slastear
        endif
      endif

c calculate tabular ear-voltage values

      if(iearform == 0)then
c     generate ear values using g-factor model 
c       load line-charge density times c = current/beta into array den

        amin = y(1,1)
        bmin = y(3,1)
        do ii=1,nit
          amin = min(amin,y(1,ii))
          bmin = min(bmin,y(3,ii))
        enddo

c     calculate current derivative for ii = 1

        del0 = y(9,2) - y(9,1)
        del1 = y(9,3) - y(9,2)
        del2 = clight*(del0 + del1)
        fact1 = (2.0 + del1/del0)/del2
        fact2 = del0/(del1*del2)

        dlna(1) = (fact1*(y(1,2) - y(1,1))
     .          - fact2*(y(1,3) - y(1,2)))/y(1,1)
        dlnb(1) = (fact1*(y(3,2) - y(3,1))
     .          - fact2*(y(3,3) - y(3,2)))/y(3,1)   
        denv(1) = 0.5*y(15,1)*(dlna(1) + dlnb(2))      
        dden(1) = fact1*(y(15,2) - y(15,1))
     .          - fact2*(y(15,3) - y(15,2))

c     calculate current derivative for ii = nit

        del0 = y(9,nit) - y(9,nit-1)
        del1 = y(9,nit-1) - y(9,nit-2)
        del2 = clight*(del0 + del1)
        fact1 = (2.0 + del1/del0)/del2
        fact2 = del0/(del1*del2)

        dlna(nit) = (fact1*(y(1,nit) - y(1,nit-1))
     .      - fact2*(y(1,nit-1) - y(1,nit-2)))/y(1,nit)
        dlnb(nit) = (fact1*(y(3,nit) - y(3,nit-1))
     .      - fact2*(y(3,nit-1) - y(3,nit-2)))/y(3,nit)
        denv(nit) = 0.5*y(15,nit)*(dlna(nit) + dlnb(nit))      
        dden(nit) = fact1*(y(15,nit)-y(15,nit-1))
     .              - fact2*(y(15,nit-1)-y(15,nit-2))

c     calculate current derivatives in beam interior

        do ii=2,nit-1
          del0 = y(9,ii) - y(9,ii-1)
          del1 = y(9,ii+1) - y(9,ii)
          del2 = clight*(del0 + del1)
          fact1 = del1/(del0*del2)
          fact2 = del0/(del1*del2)

          dlna(ii) = (fact1*(y(1,ii  ) - y(1,ii-1))
     .              + fact2*(y(1,ii+1) - y(1,ii  )))/y(1,ii)
          dlnb(ii) = (fact1*(y(3,ii  ) - y(3,ii-1))
     .              + fact2*(y(3,ii+1) - y(3,ii  )))/y(3,ii)
          denv(ii) = 0.5*y(15,ii)*(dlna(ii) + dlnb(ii))      
          dden(ii) = fact1*(y(15,ii  ) - y(15,ii-1))
     .             + fact2*(y(15,ii+1) - y(15,ii  ))

        enddo

c     set ear voltages by linear interpolation of dden and denv

        ii = 2
        epsfact = eargap/(clight*4.*pi*eps0)
        do i=1,ntval
          do while (tval(i) > y(9,ii) .and. ii < nit)
            ii = min(nit,ii + 1)
          end do

          part0 = (tval(i) - y(9,ii-1))/(y(9,ii) - y(9,ii-1)) 
          part1 = 1.0 - part0
          
          aval = max(amin,part0*y(1,ii) + part1*y(1,ii-1))
          bval = max(bmin,part0*y(3,ii) + part1*y(3,ii-1))
          betaval = part0*y(10,ii) + part1*y(10,ii-1)
          ddenval = part0*dden(ii) + part1*dden(ii-1)
          denval = part0*denv(ii) + part1*denv(ii-1)

          if(icharge < 2)then
            gval = log(rwall**2/(aval*bval))
          else
            gval = log(rwall**2/(aval*bval)) + 0.5
          endif

          if(icharge < 2)then
            earval(i) = - (epsfact/betaval)* gval*ddenval
          else
            earval(i) = - (epsfact/betaval)* (gval*ddenval - denval)
          endif

        enddo

      elseif(iearform >= 1)then 
c     generate ear-voltage values from local space-charge field ezbeam

        call getezbeam(nit,y,ezval,rpipe,icharge,.true.,lendzero,lfixed,lfailed)

        do ii=1,nit
          ezval(ii) = - eargap*ezval(ii)
        enddo

c interpolate ear voltage onto earval 

        ii = 2
        do i=1,ntval
          do while (tval(i) > y(9,ii) .and. ii < nit)
            ii = min(nit,ii + 1)
          end do

          part0 = (tval(i) - y(9,ii-1))/(y(9,ii) - y(9,ii-1))
          part1 = 1.0 - part0
          earval(i) = part0*ezval(ii) + part1*ezval(ii-1)
        enddo

c       if (printnow) then
c         printnow = .false.
c         print*,"from warp import *"
c         print*,"setearsdata = array(["
c         do i=1,nit
c           print*,"[",ezval(i),",",earval(i),",",y(9,i),",",tval(i),"],"
c         enddo
c         print*,"])"
c       endif
      
      endif

      return
      end
c===========================================================================
      subroutine getlen(nit,y,beamlen)
      use Constant
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: beamlen

c     GETLEN calculates the beam length beamlen in meters

      integer(ISZ):: ii
      real(kind=8):: sum0

c calculate beam length in meters

      sum0 = 0.0
      do ii=1,nit-1
        sum0 = sum0 + (clight*y(9,ii+1) - clight*y(9,ii))
     .                      *(y(10,ii) + y(10,ii+1))
      enddo
      beamlen = 0.5*sum0

      return
      end
c===========================================================================
      subroutine gettemp(nit,y,s,beamtemp,delv,aion,zion,
     &                   icharge,lezbeam,lperveance,lemittance,lallez,
     &                   llinear,limage,lendzero,lfixed)
      use Constant
      use CIRvartmp
      integer(ISZ):: nit
      real(kind=8):: y(16,nit),s
      real(kind=8):: beamtemp,delv,aion,zion
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero
      logical(ISZ):: lfixed

c     GETTEMP estimates the beam longitudinal "temperature" in eV
c     beamtemp from the mean-squared deviation of beta from that 
c     of a drifting distribution with beta-avg(beta) proportional to
c     z-avg(z).

      integer(ISZ):: ii
      real(kind=8):: sum0,sum1,sum2,sum3,c0,c1,c2,weight
      real(kind=8):: qtot
      real(kind=8):: boltzcon = 1.3807e-23
      real(kind=8):: degcon = 8.6173e-05

c calculate local s-derivatives of the beam variables var

      call derivs(s,SMALLPOS,0.,y,dy1,nit,aion,zion,icharge,
     &            lezbeam,lperveance,lemittance,lallez,llinear,limage,lendzero,lfixed)

c calculate required beam integrals and constant s

      sum0 = 0.0
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
      qtot = 0.0

      do ii=1,nit-1
c       --- avg(t)
        sum0 = sum0 + y(14,ii)*(y(9,ii) + y(9,ii+1))
c       --- avg(t**2)
        sum1 = sum1 + y(14,ii)*(y(9,ii)**2 + y(9,ii+1)**2)
c       --- avg(t')
        sum2 = sum2 + y(14,ii)*(dy1(9,ii) + dy1(9,ii+1))
c       --- avg(tt')
        sum3 = sum3 + y(14,ii)*(y(9,ii)*dy1(9,ii) + y(9,ii+1)*dy1(9,ii+1))
c       --- total charge
        qtot = qtot + y(14,ii)
      enddo

      sum0 = 0.5*sum0/qtot
      sum1 = 0.5*sum1/qtot
      sum2 = 0.5*sum2/qtot
      sum3 = 0.5*sum3/qtot

c calculate the velocity deviation

      c0 = sum1 - sum0**2
      c1 = sum1*sum2 - sum0*sum3
      c2 = sum3 - sum0*sum2
      
      do ii=1,nit
        dy2(1,ii) = c0/(c1+c2*y(9,ii))
        dy3(1,ii) = (dy2(1,ii) - clight*y(10,ii))**2
      enddo

c calculate mean-squared deviation

      sum0 = 0.0
      sum1 = 0.0

      do 30 ii=1,nit-1
        weight = 0.5*y(14,ii)/qtot
        sum0 = sum0 + weight*(dy2(1,ii) + dy2(1,ii+1))
        sum1 = sum1 + weight*(dy3(1,ii) + dy3(1,ii+1))
   30 continue

c calculate approximate equivalent temperature in eV
c     degcon converts K to eV
c     boltzcon is Boltzmann constant in JK**(-1)
     
      delv = sqrt(sum1)/sum0
      beamtemp = aion*amu*degcon*delv*sum0**2/boltzcon
       
      return
      end
c===========================================================================
      subroutine resetbeta(nit,y,aion,zion,inout)
      use Constant
      use Lattice
      use LatticeInternal
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: aion,zion
      character(*):: inout

c     RESETBETA modifies beta for each beam slice to account for the
c     change of electric potential entering and exiting electric 
c     dipoles.  In this version, hard-edged dipoles are asumed.
c     The element number in bend array bendlist is given by nbend.
c     When inout="in", the beam is entering the dipole, and it is 
c     exiting when inout = "out".
c     We note that time derivatives of a,b,X, and Y
c     are constant entering and exiting hard-edge dipole, so the
c     z derivatives are adjusted to account for the change in dz/dt.
c     Making changes in array var rather than var0 is appropriate
c     with the present integrator.
c     The present version is incorrect for dipole with fringe fields.

      integer(ISZ):: idipo,ii
      real(kind=8):: dipov,vslope,factor,betaold,betanew,ratio

c check bend class

      idipo = minval(cdipoid(0,:))
      if (dipoex(idipo) == 0. .and. dipoey(idipo) == 0. .and.
     &    dipov1(idipo) == 0. .and. dipov2(idipo) == 0.) return

c set constants

      if (dipox2(idipo) .ne. dipox1(idipo)) then
        vslope = (dipov2(idipo)-dipov1(idipo))/
     &           (dipox2(idipo)-dipox1(idipo))
      else
        vslope = 0.
      endif
      factor = 2.*zion*echarge/(aion*amu*clight**2)
      if(inout == "out")then
        factor = -factor
      endif
   
c adjust beta values

      do ii=1,nit
        dipov = y(5,ii)*dipoex(idipo) + y(7,ii)*dipoey(idipo) +
     &          vslope*(y(5,ii) - dipox1(idipo)) + dipov1(idipo)

        betaold = y(10,ii)
        betanew = sqrt(max(0.0,y(10,ii)**2 + factor*dipov))
        if(betanew == 0.0)then
          call remark("fatal: beam energy reduced to zero at dipole edge")
        endif
        ratio = betaold/betanew

        y(2,ii) = ratio*y(2,ii)
        y(4,ii) = ratio*y(4,ii)
        y(6,ii) = ratio*y(6,ii)
        y(8,ii) = ratio*y(8,ii)
        y(10,ii) = betanew
      enddo

      return
      end
c===========================================================================
      subroutine setacc(schedule,nit,y,ekinstart,ekinend,zlcir,zucir,deltbeam,
     &                  curmax,currise,curfall,footrise,tfoot,curfoot,
     &                  icurload,iprofile,beta0,tilt,tiltmid,tiltlen,ivelload,
     &                  lperveance,aion,zion)
      use Constant
      use Lattice
      use CIRsetacc
      character(*):: schedule
      integer(ISZ):: nit
      real(kind=8):: y(16,nit)
      real(kind=8):: ekinstart,ekinend,zlcir,zucir,deltbeam
      real(kind=8):: curmax,currise,curfall,footrise,tfoot,curfoot
      integer(ISZ):: icurload,iprofile,ivelload
      real(kind=8):: beta0,tilt,tiltmid,tiltlen,aion,zion
      logical(ISZ):: lperveance

c     SETACC uses a modified version of Kim-Smith theory to set up
c     acceleration voltages in all cells.
c     Tabulated voltage values are written to an external file and 
c     used during the subsequent runs.  
c     In this version, the target beam duration and mid-point 
c     beta are specified as functions of the longitudinal
c     position s in the interpreter routine named in schedule. Instantaneous 
c     acceleration is assumed in the cells. 
c     note: We assume sector bends and ignore any initial tilt.
c Assumes that space for accl elements has been fully allocated.

      real(kind=8):: curvfact,delstart,stotal,con0,cdelt0,snextacc,start,stop
      real(kind=8):: smidacc,dels,bendangsum,bendlensum,radfact
      real(kind=8):: quadfact,betamidnext,deltnext
      real(kind=8):: ctmidnext,deltratio,sumhead,sumtail,curv,curvavg
      real(kind=8):: sbend,sig0,sig1,sigfact,bendfact,sumfact
      real(kind=8):: wavenum,cdelt,ratio,xheadapp,xtailapp,cdeltold
      real(kind=8):: delbetaold,delsi,delbetahead,delbetatail,cdelti,pervfact
      real(kind=8):: beti,gambeti,pervi,gammold,delbeta
      real(kind=8):: alfcuri,snextbend,hlp,rpipe
      integer(ISZ):: i,ibend0,ibend,iacc,inextacc,iezmid,ib

c check parameters
      if(naccl == 0)then
        call remark("warning: no accelerating cells. Acc fields cannot be set")
        return
      elseif(ekinend == ekinstart)then
        call remark("warning: no energy gain specified. Acc fields set to zero")
      endif

      radfact = pi/180.0

c Allocate space
      nitacc = nit
      call gchange("CIRsetacc",0)

c set index for beam mid-point
      iezmid = (nit+1)/2

c set switch for inclusion of curvature term in displcement expression
      curvfact = 1.0

c set starting position and total path length
      delstart = 0.5*(acclzs(0) + acclze(0)) + zlatstrt
      stotal = zucir - delstart

c set initial time for beam slices
      con0 = deltbeam/(nit-1)
      do i=1,nit
        y(9,i) = con0*(i-1)
      enddo

c set initial beta and current values
      call setcur(nit,y,curmax,currise,curfall,deltbeam,footrise,tfoot,curfoot,
     &            icurload,iprofile,beta0,tilt,tiltmid,tiltlen,ivelload)

c initialize iterated quantities at first cell
      do i=1,nit
        ct0(i) = clight*y(9,i)
        ctold(i) = ct0(i) + delstart/y(10,i)
        betaold(i) = y(10,i)
      enddo
      cdelt0 = ct0(nit) - ct0(1)

c set number of bends before first acceleration cell
      ibend0 = -1
      do while (ibend0+1 <= nbend .and. bendzs(ibend0+1) < acclzs(0))
        ibend0 = ibend0 + 1
      end do

c begin stepping through acceleration cells

      ibend = ibend0

      do iacc=0,naccl

c find location of next acceleration cell mid-point
        inextacc = iacc+1
        if (inextacc <= naccl) then
          snextacc = 0.5*(acclzs(inextacc) + acclze(inextacc)) + zlatstrt
        else
          snextacc = zucir
        endif

c calculate distance dels between cells and gap length gaplen
        start = acclzs(iacc) + zlatstrt
        stop = acclze(iacc) + zlatstrt
        smidacc = 0.5*(start + stop)
        dels = snextacc - smidacc

c count bends before next acceleration cell and find total angle
        if(nbend > 0)then
          ibend = ibend0
          bendangsum = 0.0
          bendlensum = 0.0
          snextbend = bendzs(ibend) + zlatstrt
          do while (snextbend < snextacc .and. ibend < nbend)
            ibend = ibend + 1
            bendangsum=bendangsum+(bendze(ibend)-bendzs(ibend))/bendrc(ibend)
            bendlensum = bendlensum + (bendze(ibend)-bendzs(ibend))
            snextbend = bendzs(ibend) + zlatstrt
          end do
        endif

c Assume magnetic quads
        quadfact = 2.0

c       --- Call interpreter routine which returns the midpoint velocity and
c       --- duration of the beam at the current accl element.
        iaccl = iacc
        call execuser(schedule)
        betamidnext = betamid_set
        deltnext = deltbeam_set

c calculate midpoint arrival time at next acceleration cell
        ctmidnext = ctold(iezmid) + dels/betamidnext

c calculate arrival times of all beam slices
        deltratio = clight*deltnext/cdelt0
        do i=1,nit
          ctnext(i) = ctmidnext + deltratio*(ct0(i)-ct0(iezmid))
        enddo

c estimate beta values between acceleration cells, ignoring bends
        do i=1,nit
          betanext(i) = dels/(ctnext(i) - ctold(i))
        enddo

        if(ibend > ibend0)then

c   estimate next beta for head and tail
          sumhead = 0.0
          sumtail = 0.0 
          curv = curvfact*bendangsum/bendlensum            
          curvavg = bendangsum/dels

          do ib=ibend0+1,ibend
            sbend = bendzs(ib) + zlatstrt

c          get undepressed phase advance at beam mid-point
c          note: the time here is wrong but unimportant for fixed quads
            call getphase(sig0,sig1,hlp,betamidnext,0.0,ctmidnext/clight,
     &                      sbend,0.,0.,0.,0.,.false.)
            sigfact = quadfact*(1.0 - cos(sig0))/sin(sig0)

c          calculate sums at head and tail
c             --- Assume magnetic bend
            bendfact = 1.0
            sumfact = bendfact*radfact*curvavg*
     &                        (bendze(ibend)-bendzs(ibend))/bendrc(ibend)

            delbeta = (betaold(1)/betaold(iezmid)) - 1.0
            wavenum = (sig0 - sigfact*delbeta)/(2.0*hlp)
            sumhead = sumhead + sumfact/(wavenum**2 + bendfact*curv*curvavg)
            delbeta = (betaold(nit)/betaold(iezmid)) - 1.0
            wavenum = (sig0 - sigfact*delbeta)/(2.0*hlp)
            sumtail = sumtail + sumfact/(wavenum**2 + bendfact*curv*curvavg)
          enddo

c         estimate next beta at head and tail
          cdelt = ctnext(1) - ctold(1)
          ratio = 0.5*(sumhead + dels)/cdelt
          betanext(1) = ratio + sqrt(ratio**2 - betamidnext*sumhead/cdelt)
          cdelt = ctnext(nit) - ctold(nit)
          ratio = 0.5*(sumtail + dels)/cdelt
          betanext(nit) = ratio + sqrt(ratio**2 - betamidnext*sumtail/cdelt)

c diagnostic: estimate displacment at ends for mag bends

          delbeta = (betaold(1)/betaold(iezmid)) - 1.0
          wavenum = (sig0 - sigfact*delbeta)/(2.0*hlp)
          delbeta = 1.0 - (betaold(iezmid)/betaold(1))
          xheadapp = bendfact*curvavg*delbeta/
     &                      (wavenum**2 + bendfact*curv*curvavg)
          delbeta = (betaold(nit)/betaold(iezmid)) - 1.0
          wavenum = (sig0 - sigfact*delbeta)/(2.0*hlp)
          delbeta = 1.0 - (betanext(iezmid)/betanext(nit))
          xtailapp = bendfact*curvavg*delbeta/
     &                      (wavenum**2 + bendfact*curv*curvavg)
c end diagnostic

c   estimate path-length change for remaining slices
          sumhead = 0.0
          sumtail = 0.0
          do i=1,nit
            sumslice(i) = 0.0
          enddo
          cdeltold = ctold(nit) - ctold(1)
          do ib=ibend0,ibend
            sbend = bendzs(ib) + zlatstrt
            rpipe = bendap(ib)

c           get undepressed phase advance at beam midpoint
c           note: the time here is wrong but unimportant for fixed quads
            call getphase(sig0,sig1,hlp,betamidnext,0.0,ctmidnext/clight,
     &                      sbend,0.,0.,0.,0.,.false.)
            sigfact = quadfact*(1.0 - cos(sig0))/sin(sig0)

c           calculate sums at heat and tail
c             --- Assume magnetic bend
            bendfact = 1.0
            sumfact = bendfact*radfact*curvavg*
     &                        (bendze(ibend)-bendzs(ibend))/bendrc(ibend)

            delbetaold = 1.0 - betaold(iezmid)/betaold(1)
            wavenum = (sig0 - sigfact*delbetaold)/(2.0*hlp)
            sumhead = sumhead + sumfact/(wavenum**2 + bendfact*curv*curvavg) 
            delbetaold = 1.0 - betaold(iezmid)/betaold(nit)
            wavenum = (sig0 - sigfact*delbetaold)/(2.0*hlp)
            sumtail = sumtail + sumfact/(wavenum**2 + bendfact*curv*curvavg)

c           estimate cdelt at the bend            
            delsi = sbend - zucir
            delbeta = (1.0/betanext(nit)) - (1.0/betanext(1))
            delbetahead =  1.0 - (betamidnext/betanext(1)**2)
            delbetatail =  1.0 - (betamidnext/betanext(nit)**2)
            cdelti = cdeltold + delsi*delbeta
     .                    + (sumtail*delbetatail - sumhead*delbetahead)

            if (lperveance) then
              alfcuri = zion*echarge/(4.*pi*eps0*amu*aion*clight**3)
              pervfact = 2.0*alfcuri
            else
              pervfact = 0.
            endif

            do i = 2,nit-1

c             estimate sigma0 at bend
              delbetaold = 1.0 - betaold(iezmid)/betaold(i)
              wavenum = (sig0 - sigfact*delbetaold)/(2.0*hlp)

c             estimate perveance term at bend
              beti = betamidnext*(1.0 + delbetaold)
              gambeti = beti/sqrt(1.0 - beti**2)
              pervi = pervfact*cdelt0*y(13,i)/(cdelti*gambeti**3)

c           evaluate contributions to pathlength sums
              sumslice(i) = sumslice(i) + sumfact/
     .                 (wavenum**2 - (pervi/rpipe**2) + bendfact*curv*curvavg)
            enddo
          enddo  

c estimate beta for all remaining slices

          do i = 2,nit-1
            cdelt = ctnext(i) - ctold(i)
            ratio = 0.5*(sumslice(i) + dels)/cdelt
            betanext(i) = ratio+sqrt(ratio**2-betamidnext*sumslice(i)/cdelt)
          enddo            
        endif

c calculate voltage values at iacc

        do i=1,nit
           gammold = sqrt(1.0/(1.0 - betaold(i)**2))
           tout(i) = ctold(i)/clight
           accout(i) = 0.5*(clight*gammold**3/
     &      (zion*echarge/(aion*amu*clight)))*(betanext(i)**2 - betaold(i)**2)
           earout(i) = 0.0
           ctold(i) = ctnext(i) 
           betaold(i) = betanext(i)
        enddo

c     write acceleration voltages 
        acclts(iacc) = tout(1)
        accldt(iacc) = tout(2) - tout(1)
        do i=1,nit
          acclet(i-1,iacc) = (accout(i) + earout(i))/(acclze(iacc)-acclzs(iacc))
        enddo

c     update iterated quantities

        if(nbend > 0)then
          ibend0 = ibend
        endif

c     end loop over acceleration cells
      enddo

c diagnostic: plot estimated displacement, beta, and current

c     if(nbend == 0) return

c     iturn = itrn
c     itsav = itmax
c     itmax = nit

c     s = snext
c     tmid = ctnext(iezmid)/c
c     call getphase(sig0,sig1,betamidnext,0.0,tmid,snext)
c     call getphase(sighead,sig1,betanext(1),0.0,tmid,snext)
c     call getphase(sigtail,sig1,betanext(nit),0.0,tmid,snext)

c     pervfact = 2.0*iperv*qovamu/alfcur
c     sigfact = quadfact*(1.0 - cos(sig0))/sin(sig0)
c     sighchk = sig0 - sigfact*((betanext(1)/betamidnext) - 1.0)
c     sigtchk = sig0 - sigfact*((betanext(nit)/betamidnext) - 1.0)

c     do i=1,nit
c       t(i) = ctnext(i)/c
c       beta(i) = betanext(i)
c       cur(i) = y(13,i)*cdelt0/(ctnext(nit) - ctnext(1)) 

c       sigi = sig0 - sigfact*((beta(i)/betamidnext) - 1.0)
c       gambeti = beta(i)/sqrt(1.0 - beta(i)**2)
c       pervi = pervfact*cur(i)/gambeti**3
c       wavenumi = sigi/(2.0*hlp)
c       delbeti = 1.0 - betamidnext/beta(i)
c       x(i) = bendfact*curvavg*delbeti 
c    .                   /(wavenumi**2 - (pervi/piperad**2)
c    .                                    + bendfact*curv*curvavg)
c       var(1,i) = 0.0
c       var(2,i) = 0.0
c       var(3,i) = 0.0
c       var(4,i) = 0.0
c       var(5,i) = x(i)
c       var(6,i) = 0.0
c       var(7,i) = 0.0
c       var(8,i) = 0.0
c       var(9,i) = t(i)
c       var(10,i) = beta(i)
c     enddo
c     call plot0(var)
c  
c     iturn = 0
c     itmax = itsav
c     s = 0.0

c end diagnostic

      return
      end

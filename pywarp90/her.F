#include "top.h"
c=============================================================================
c@(#) File HER.M, version $Revision: 3.7 $, $Date: 2001/07/19 17:21:45 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  It is the source file for the package HER of the PIC code WARP,
c  based on the envelope/fluid model which originated in CIRCE.
c
c=============================================================================
c "Moneo, fuge litora Circes"
c "I warn you, keep away from Circe's shores!"
c Macareus, a companion of Odysseus, to Aeneas, in
c Metamorphoses 14.245 by P. Ovidius Naso.
c=============================================================================
      subroutine herinit
      use HERversion

c  Called at first reference to package (not nec. a "run" etc.).

      call hervers (STDOUT)

      return
      end
c=============================================================================
      subroutine hervers (iout)
      use HERversion
      integer(ISZ):: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"HERMES",versher)
      return
      end
c=============================================================================
      subroutine hergen()
      use HERvars
      use HERflags
      use HERtmp
      use HERfield
      use HERhist
      use HERvartmp
      use LatticeInternal
      use Picglb
      use Picglbrz
      use InMeshrz
      use Fieldsrz
      use InPart

c  Invoked by the GENERATE command, it allocates storage
c  for HERMES.

      integer(ISZ):: ii

c  Send a startup message to the standard output

      call remark(" ***  HERMES generating")

c  Calculate derived quantities

      call derivqty
      
c     --- Reset counter
      it = 0

c     --- allocate dynamic arrays ("group change")
      if (niz == 0) niz = 101
      call gchange("HERvars",0)
      niztmp = niz
      call gchange("HERtmp",0)
      nizfield = niz
      call gchange("HERfield",0)
      nizhist = niz
      lhher = 1000
      call gchange("HERhist",0)
      nizvartmp = niz
      call gchange("HERvartmp",0)

c     --- reset Lattice arrays to final size
      call resetlat
      
c     --- Initialize the LatticeInternal arrays
c     --- The choice for nzl is arbitrary; better choices may exist.
      nzl = niz + 2
      nzlmax = nzl
      dzl = (zimax - zimin)/niz
      zlmin = zimin - dzl
      zlmax = zimax + dzl
      dzli = 1./dzl
      call gchange("LatticeInternal",0)
      call setlattzt(0.,0.,0)
      
      if (icharge == 0) then
c --- Set the g-factor for the constant g-factor model
        do ii=1,niz
          gtemp(ii) = 2*log(1.6)
        enddo
      endif

      if (icharge == 5 .or. icharge == 6) then
c     Read in Bessel data
        call execuser ("readbessel")
      endif

      if (icharge == 7) then
c ---   Allocate space for the RZ field solver
        call remark (" ***  field solver package RZ generating for Hermes  ***")

c ---   Create the dynamic arrays for fields
        nmrz  = max(nr, nz)
        nrz = (nr+1) * (nz+1)
        call gallot ("Fieldsrz", 0)
 
c ---   Calculate mesh dimensioning quantities
        dr = (rmmax - rmmin) / nr
        dz = (zmmax - zmmin) / nz
        do ii = 0, nr
          rmesh(ii) = ii * dr + rmmin
        enddo

        do ii = 0, nz
          zmesh(ii) = ii * dz + zmmin
        enddo

        call zeroarry(rho,nrz)

c ---   Initial call to fieldsolver in order to initialize attz, kzsq, etc.
        call vprz (1)

      endif

      return
      end
c=============================================================================
      subroutine herexe()

      use InGen
      use Picglb
      use Beam_acc
      use HERvars
      use HERvartmp
      use HERgeneral
      use HERflags

c  Invoked by the STEP command, it runs the HERMES code.
      
      lfail = .false.
      
      if (it == 0) then
c       --- Start at time = 0.
        time = 0.
c       --- Save the initial history
        call savehermesvars(time,var,niz,it,lsavehist)
c       --- Find initial forces
        call getderivs(time,dt,var,dy,niz,aion,zion,fviscous,icharge,iimage,
     &            lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail)
        if (lfail) return
        if (lherout) then
c         --- Send a startup message to the standard output
          call remark(" ***  HERMES running")
        endif
      endif
      
c     --- Run the HERMES kernel
      it = it+1
      time = time + dt
c     --- Advance the envelope and orbits (one step)
      call stephr (var,dy,niz,time,dt,aion,zion,fviscous,
     &              icharge,iimage,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &              lviscous,lcurgrid,lfail)
      if (lfail) return

      call savehermesvars(time,var,niz,it,lsavehist)

      return
      end

c=============================================================================
      subroutine herfin()

c  Deallocates dynamic storage used by test driver

      call gfree("HERvars")
      call gfree("HERfield")
      call gfree("HERtmp")
      call gfree("HERvartmp")
      call gfree("HERbessel")
      call gfree("HERhist")

      return
      end

c=============================================================================
      subroutine extebher(t,dt,y,niz)
      use Constant
      use Beam_acc
      use Picglb
      use Lattice
      use LatticeInternal
      use Mult_data
      use HERfield
      integer(ISZ):: niz
      real(kind=8):: t,dt,y(16,niz)

c  Computes the forces needed for the general-lattice envelope calculation.
c  Includes both hard edged quadrupole elements
c  and hard-edge multipole elements, and elements defined by their
c  multipoles.  Only dipole and quadrupolar fields are included. Uniform 
c  focusing fields are ignored for now, but can be implemented easily.
c  Input:
c     t    is the starting time of the current time step
c     dt   is the size of the current time step
c     y    is the array hold the envelope information
c     niz  is the number of slices

      real(kind=8):: vz,vzi,zl,zr,fl,fr,frac,dta,dzi,mltsf,angle,rotfuzz,term
      integer(ISZ):: io,iele,iz,ii,izm,j
      
      dta = abs(dt)
      
      zlmin = y(9,1)
      zlmax = y(9,niz)
      do iz=1,niz
        zlmin = min (zlmin,y(9,iz))
        zlmax = max (zlmax,y(9,iz))
      enddo
      dzl = (zlmax - zlmin)/niz
      zlmin = zlmin - dzl
      zlmax = zlmax + dzl
      dzli = 1./dzl
      
      call setlattzt(0.,t,0)

      if (rwall > 0.) then
        do iz=1,niz
c         --- Default value of rpipe
           rpipe(iz) = rwall
        enddo
      endif

      do iz=1,niz
c       --- Zero out the fields
        ex(iz) = 0.
        ey(iz) = 0.
        ez(iz) = 0.
        bx(iz) = 0.
        by(iz) = 0.
        bz(iz) = 0.
        dedx(iz) = 0.
        dbdx(iz) = 0.
      enddo

c---------------------------------------------------------------------------
c       --- Hard edge quadrupoles

      if (quads) then
c       --- If in a quad, get field and pipe.

        do io=1,nquadol
          do iz=1,niz
            vz = y(10,iz)*clight
            vzi = 1./dvnz(vz)
c           --- find z-cell in which particle lies
            j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = y(9,iz) - vz * dta / 2.
            fl = 0.
            if (zl.ge.cquadzs(j,io) .and. zl.lt.cquadze(j,io)) then
              fl = 1.
              iele = cquadid(j,io)
            endif
c           --- "right" end of velocity advance step
            zr = y(9,iz) + vz * dta / 2.
            fr = 0.
            if (zr.ge.cquadzs(j,io) .and. zr.lt.cquadze(j,io)) fr = 1.
c            --- residence fraction
            frac = fl
            if (fl .gt. fr) frac = (cquadze(j,io)-zl)*vzi/dta
            if (fr .gt. fl) frac = (zr-cquadzs(j,io))*vzi/dta
c           --- set the field
            if (frac .gt. 0.) then
              dedx(iz) = dedx(iz) + cquadde(j,io) * frac
              dbdx(iz) = dbdx(iz) + cquaddb(j,io) * frac
            endif
            if (ntquad > 0) then
              ii = int((t - quadts(iele))/quaddt(iele))
              if (0 <= ii .and. ii < ntquad) then
                fr = (t - quadts(iele))/quaddt(iele) - ii
                term = quadet(ii,iele)*frac*(1.-fr) + quadet(ii+1,iele)*frac*fr
                dedx(iz) = dedx(iz) + term
                term = quadbt(ii,iele)*frac*(1.-fr) + quadbt(ii+1,iele)*frac*fr
                dbdx(iz) = dbdx(iz) + term
              endif
            endif
c           --- use current position to find the aperture, offset and rotation angle
            if (y(9,iz).ge.cquadzs(j,io) .and. y(9,iz).lt.cquadze(j,io)
     &         .and. quadap(iele) > 0.) then
              rpipe(iz) = quadap(iele)
              pipeox(iz) = qoffx(iele)
              pipeoy(iz) = qoffy(iele)
              pipeph(iz) = 0.
            endif
          enddo
        enddo
      endif

c---------------------------------------------------------------------------
c     --- Hard edge multipoles, can include linearly focusing 
c     --- quadrupole or azimuthally symmetric field 
 
      if (heles) then

        do io=1,nheleol
          do iz=1,niz
            vz = y(10,iz)*clight
            vzi = 1./dvnz(vz)
c           --- find z-cell in which particle lies
            j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = y(9,iz) - vz * dta / 2.
            fl = 0.
            if (zl.ge.chelezs(j,io) .and. zl.lt.cheleze(j,io)) then
              fl = 1.
              iele = cheleid(j,io)
            endif
c           --- "right" end of velocity advance step
            zr = y(9,iz) + vz * dta / 2.
            fr = 0.
            if (zr.ge.chelezs(j,io) .and. zr.lt.cheleze(j,io)) fr = 1.
c            --- residence fraction
            frac = fl
            if (fl .gt. fr) frac = (cheleze(j,io)-zl)*vzi/dta
            if (fr .gt. fl) frac = (zr-chelezs(j,io))*vzi/dta
c           --- set the field
            if (frac .gt. 0.) then
              do ii = 1, max(helenm(iele),helene(iele)) 
                if ((nint(hele_n(ii,iele)) == 0) .and.
     &                   (nint(hele_v(ii,iele)) == 0)) then
c                 --- longitudinal field
c                 --- radial field is ignored for now
                  ez(iz) = ez(iz) + heleae(ii,iele)
                  bz(iz) = bz(iz) + heleam(ii,iele)
                else if ((nint(hele_n(ii,iele)) == 1) .and. 
     &              (nint(hele_v(ii,iele)) == 0)) then
c                 --- dipole field
                  ex(iz) =   heleae(ii,iele) * cos(helepe(ii,iele))
                  ey(iz) = - heleae(ii,iele) * sin(helepe(ii,iele))
                  bx(iz) =   heleam(ii,iele) * sin(helepm(ii,iele))
                  by(iz) =   heleam(ii,iele) * cos(helepm(ii,iele))
                else if ((nint(hele_n(ii,iele)) == 2) .and. 
     &              (nint(hele_v(ii,iele)) == 0)) then
c                 --- quadrupole field  
                  rotfuzz = 1.e-3
                  if ( abs(cos(helepm(ii,iele))*sin(helepm(ii,iele))) > rotfuzz )
     &              call remark("Warning: magnetic multipole structure incompatible")
                  if ( abs(cos(helepe(ii,iele))*sin(helepe(ii,iele))) > rotfuzz )
     &              call remark("Warning: electric multipole structure incompatible")
                  term = heleae(ii,iele)*cos(helepe(ii,iele)) * frac
                  dedx(iz) = dedx(iz) + term
                  term = heleam(ii,iele)*cos(helepm(ii,iele)) * frac
                  dbdx(iz) = dbdx(iz) + term
c                 --- using the rotation of the electrostatic multipole for the pipe rotation
                  angle = helepe(ii,iele)
                
                endif
              enddo
            endif           
            if (y(9,iz).ge.chelezs(j,io) .and. y(9,iz).lt.cheleze(j,io)
     &          .and. heleap(iele) > 0.) then
c                  --- use current position to find the aperture, offset and rotation angle
              rpipe(iz) = heleap(iele)
              pipeox(iz) = heleox(iele)
              pipeoy(iz) = heleoy(iele)
              pipeph(iz) = angle
            endif
          enddo
        enddo
      endif
      
c---------------------------------------------------------------------------
c     --- Extract any linear focusing field component from the multipole
c     --- (emlt and mmlt) elements.

      if (emlts) then

        do io=1,nemltol
          do iz=1,niz
c           --- find z-cell in which slice boundary lies
            j = max(0., (y(9,iz) - zlmin)*dzli + 0.5)
            iele = cemltim(j,io)
            dzi = 1./dzemlt(iele)
            mltsf = cemltsc(j,io) + cemltsf(j,io)
c           --- find z location relative to multipole data
            izm = int((y(9,iz)-cemltzs(j,io))*dzi)
            fr  =     (y(9,iz)-cemltzs(j,io))*dzi - izm
            if (izm .lt. 1) then
              izm = 1
              fr = 0.
              iele = 0
            elseif (izm .gt. nzemlt(iele)-1) then
              izm = nzemlt(iele)-1
              fr = 0.
            iele = 0
            endif
            fl = 1. - fr
            do ii=1,nesmult
              if (iele .gt. 0) then
                if (emlt_n(ii) == 0 .and. emlt_v(ii) == 0) then
c                 --- longitudinal field
c                 --- radial field ignored for now
                  ez(iz) = ez(iz)
     &                   + mltsf*(esemlt (izm,ii,iele)*fl + esemlt (izm+1,ii,iele)*fr)
                elseif (emlt_n(ii) == 1 .and. emlt_v(ii) == 0) then
                  ex(iz) = ex(iz) + mltsf*esemlt(izm  ,ii,iele)*fl*cos(emltph(iele)+esemltph(izm  ,ii,iele))
                  ex(iz) = ex(iz) + mltsf*esemlt(izm+1,ii,iele)*fl*cos(emltph(iele)+esemltph(izm+1,ii,iele))
                  ey(iz) = ey(iz) - mltsf*esemlt(izm  ,ii,iele)*fl*sin(emltph(iele)+esemltph(izm  ,ii,iele))
                  ey(iz) = ey(iz) - mltsf*esemlt(izm+1,ii,iele)*fl*sin(emltph(iele)+esemltph(izm+1,ii,iele))
                elseif (emlt_n(ii) == 2 .and. emlt_v(ii) == 0) then
c                 --- Apply the rest of the multipoles.
                  dedx(iz) = dedx(iz)
     &                     + mltsf*(esemlt(izm,ii,iele)*fl + esemlt(izm+1,ii,iele)*fr)
                endif
              endif
            enddo
c           --- use current position to find the aperture, offset and rotation angle

            if (y(9,iz).ge.cemltzs(j,io) .and. y(9,iz).lt.cemltze(j,io)
     &         .and. emltap(iele) > 0.) then
              rpipe(iz) = emltap(iele)
              pipeox(iz) = emltox(iele)
              pipeoy(iz) = emltoy(iele)
              pipeph(iz) = emltph(iele)
            endif
          enddo
        enddo
      endif

      if (mmlts) then

        do io=1,nemltol
          do iz=1,niz
c           --- find z-cell in which slice boundary lies
            j = max(0., (y(9,iz) - zlmin)*dzli + 0.5)
            iele = cmmltim(j,io)
            dzi = 1./dzmmlt(iele)
            mltsf = cmmltsc(j,io) + cmmltsf(j,io)
c           --- find z location relative to multipole data
            izm = int((y(9,iz)-cmmltzs(j,io))*dzi)
            fr  =     (y(9,iz)-cmmltzs(j,io))*dzi - izm
            if (izm .lt. 1) then
              izm = 1
              fr = 0.
              iele = 0
            elseif (izm .gt. nzmmlt(iele)-1) then
              izm = nzmmlt(iele)-1
              fr = 0.
              iele = 0
            endif
            fl = 1. - fr
            do ii=1,nmsmult
              if (iele .gt. 0) then
                if (mmlt_n(ii) == 0 .and. mmlt_v(ii) == 0) then
c                 --- Apply accelerating field
c                 --- Radial field is ignored for now
                  bz(iz) = bz(iz)
     &                   +     mltsf*(msmmlt (izm,ii,iele)*fl + msmmlt (izm+1,ii,iele)*fr)
                elseif (mmlt_n(ii) == 1 .and. mmlt_v(ii) == 0) then
                  bx(iz) = bx(iz) + mltsf*msmmlt(izm  ,ii,iele)*fl*sin(mmltph(iele)+msmmltph(izm  ,ii,iele))
                  bx(iz) = bx(iz) + mltsf*msmmlt(izm+1,ii,iele)*fl*sin(mmltph(iele)+msmmltph(izm+1,ii,iele))
                  by(iz) = by(iz) + mltsf*msmmlt(izm  ,ii,iele)*fl*cos(mmltph(iele)+msmmltph(izm  ,ii,iele))
                  by(iz) = by(iz) + mltsf*msmmlt(izm+1,ii,iele)*fl*cos(mmltph(iele)+msmmltph(izm+1,ii,iele))
               elseif (mmlt_n(ii) == 2 .and. mmlt_v(ii) == 0) then
c                 --- Apply the rest of the multipoles.
                  term = mltsf*(msmmlt(izm,ii,iele)*fl + msmmlt(izm+1,ii,iele)*fr)
                  dbdx(iz) = dbdx(iz) + term
                endif
              endif
            enddo
c           --- use current position to find the aperture, offset and rotation angle
            if (y(9,iz).ge.cmmltzs(j,io) .and. y(9,iz).lt.cmmltze(j,io)
     &          .and. mmltap(iele) > 0.) then
              rpipe(iz) = mmltap(iele)
              pipeox(iz) = mmltox(iele)
              pipeoy(iz) = mmltoy(iele)
              pipeph(iz) = mmltph(iele)
            endif
          enddo
        enddo
      endif

c---------------------------------------------------------------------------
c     --- Get accelerating gradient from accl gaps.
      if (accls) then
c       --- If in a quad, get field and pipe.

        do io=1,nacclol
          do iz=1,niz
            vz = y(10,iz)*clight
            vzi = 1./dvnz(vz)
c           --- find z-cell in which particle lies
            j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = y(9,iz) - vz * dta / 2.
            fl = 0.
            if (zl.ge.cacclzs(j,io) .and. zl.lt.cacclze(j,io)) then
              fl = 1.
              iele = cacclid(j,io)
            endif
c           --- "right" end of velocity advance step
            zr = y(9,iz) + vz * dta / 2.
            fr = 0.
            if (zr.ge.cacclzs(j,io) .and. zr.lt.cacclze(j,io)) fr = 1.
c            --- residence fraction
            frac = fl
            if (fl .gt. fr) frac = (cacclze(j,io)-zl)*vzi/dta
            if (fr .gt. fl) frac = (zr-cacclzs(j,io))*vzi/dta
c           --- set the field
            if (frac .gt. 0.) then
              ez(iz) = ez(iz) + cacclez(j,io) * frac
            endif
c           --- Get time varying accelerating field.
            if (ntaccl > 0) then
              ii = int((t - acclts(iele))/accldt(iele))
              if (0 <= ii .and. ii < ntaccl) then
                fr = (t - acclts(iele))/accldt(iele) - ii
                ez(iz) = ez(iz) + acclet(ii  ,iele)*frac*(1.-fr) +
     &                            acclet(ii+1,iele)*frac*    fr
              endif
            endif
c           --- use current position to find the aperture, offset and rotation angle
            if (y(9,iz).ge.cacclzs(j,io) .and. y(9,iz).lt.cacclze(j,io)
     &          .and. acclap(iele) > 0.) then
              rpipe(iz) = acclap(iele)
              pipeox(iz) = acclox(iele)
              pipeoy(iz) = accloy(iele)
              pipeph(iz) = 0.
            endif
          enddo
        enddo
      endif
c---------------------------------------------------------------------------
c     --- Get inverse radius of curvature from bend elements
      if (bends) then

        do iz=1,niz
c         --- find z-cell in which particle lies
          j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c         --- use current position to find the inverse radius of curvature, aperture, and offset
          if (y(9,iz).ge.cbendzs(j) .and. y(9,iz).lt.cbendze(j)) then
            iele = cbendid(0)
            bendcurv(iz) = 1./bendrc(iele)
            if (bendap(iele) > 0.) then
              rpipe(iz) = bendap(iele)
              pipeox(iz) = bendox(iele)
              pipeoy(iz) = bendoy(iele)
              pipeph(iz) = 0.
            endif
          endif
        enddo
      endif

c---------------------------------------------------------------------------

c     --- Get bending field from dipole elements
      if (dipos) then

        do io=1,ndipool
          do iz=1,niz
c           --- find z-cell in which particle lies
            j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = y(9,iz) - vz * dta / 2.
            fl = 0.
            if (zl.ge.cdipozs(j,io) .and. zl.lt.cdipoze(j,io)) then
              fl = 1.
              iele = cdipoid(j,io)
            endif
c           --- "right" end of velocity advance step
            zr = y(9,iz) + vz * dta / 2.
            fr = 0.
            if (zr.ge.cdipozs(j,io) .and. zr.lt.cdipoze(j,io)) fr = 1.
c            --- residence fraction
            frac = fl
            if (fl .gt. fr) frac = (cdipoze(j,io)-zl)*vzi/dta
            if (fr .gt. fl) frac = (zr-cdipozs(j,io))*vzi/dta
c           --- set the field
            if (frac .gt. 0.) then
              ex(iz) = ex(iz) + cdipoex(j,io) * frac
              ey(iz) = ey(iz) + cdipoey(j,io) * frac
              bx(iz) = bx(iz) + cdipobx(j,io) * frac
              by(iz) = by(iz) + cdipoby(j,io) * frac
            endif

c            --- use current position to find the aperture, offset and rotation angle
            if (y(9,iz).ge.cdipozs(j,io) .and. y(9,iz).lt.cdipoze(j,io)) then
              if (dipoap(iele) > 0.) then
                rpipe(iz) = dipoap(iele)
                pipeox(iz) = 0.
                pipeoy(iz) = 0.
                pipeph(iz) = 0.
              endif
            endif
          enddo
        enddo

      endif

c---------------------------------------------------------------------------
c       --- Get aperture from drift elements
      if (drfts) then
        
        do io=1,ndrftol
          do iz=1,niz
c           --- find z-cell in which particle lies
            j = max(0., (y(9,iz) - zlmin) * dzli + 0.5)
c            --- use current position to find the aperture, offset and rotation angle
            if (y(9,iz).ge.cdrftzs(j,io) .and. y(9,iz).lt.cdrftze(j,io)) then
              iele = cdrftid(j,io)
              if (drftap(iele) > 0.) then
                rpipe(iz) = drftap(iele)
                pipeox(iz) = drftox(iele)
                pipeoy(iz) = drftoy(iele)
                pipeph(iz) = 0.
              endif
            endif
          enddo
        enddo
        
      endif

c---------------------------------------------------------------------------

      do iz=1,niz
c ---   Print warning if rpipe is zero
        if (rpipe(iz) == 0.) then
          call remark("WARNING: rpipe is zero, make sure that either top.rwall and/or the lattice element apertures are set.")
        endif
      enddo

      return
      end
c=============================================================================
      subroutine getderivs(t,dt,y,dy,niz,aion,zion,fviscous,
     &                  icharge,iimage,lezbeam,lperveance,lemittance,lallez,
     &                  lezcenter,lviscous,lcurgrid,lfail)
      use Constant
      use HERfield
      integer(ISZ):: niz
      real(kind=8):: t,dt,y(16,niz),dy(16,niz)
      real(kind=8):: aion,zion
      real(kind=8):: fviscous
      integer(ISZ):: icharge,iimage
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail

c  Computes the forces needed for the general-lattice envelope calculation.
c  Includes both residence corrected hard edged quadrupole elements
c  and hard-edge multipole elements, and elements defined by their
c  multipoles.  Only electric and magnetic quadrupole fields, dipoles, as well as
c  longitudinal electric and magnetic fields are taken into account, together
c  with the linear space-charge forces.

      integer(ISZ):: iz
      real(kind=8):: alfcur
      real(kind=8):: y1i,y3i,y1py3i,gamma,gambeta
      real(kind=8):: rigidityi,quadfld,perv
      real(kind=8):: xemit,yemit,dsdt,dlnbeta,dbetadt
      real(kind=8):: bendex,bendey,xbend
      real(kind=8):: xval,yval,rpipei2
      real(kind=8):: aviscous(niz)
      

c     --- Get the applied fields.
      call extebher(t,dt,y,niz)

c     --- Image forces are not yet implemented in HERMES

c     --- Calculate the beam radius if needed
      if (icharge < 7) call getradius(y,niz,rpipe,icharge)
      if (icharge < 5) call getgfactor(niz,rpipe,icharge)
      
c     --- Calculate auxiliary variables needed for fieldsolhr and getcurrent
      call setrhohr (niz,y,rpipe,icharge,lfail)
      if (lfail) return

c     --- Calculate the current (needed for the envelope equation)
      call getcurrent(y,niz,rpipe,icharge,lcurgrid,lfail)
      if (lfail) return
      
c     --- Get self axial fields of beam.
      call fieldsolhr (niz,y,ezbeam,rpipe,icharge,lezbeam,lezcenter,lfail)
      if (lfail) return
      
c     --- Get the image effect on the beam envelope
      if (iimage > 0) then
        call getimage(y,niz,iimage,lfail)
        if (lfail) return
      endif
      
      if (lviscous .and. lezbeam) then
        call getviscous(y,niz,aviscous,fviscous)
      endif

c     --- Set some temporaries.
      alfcur = (4.*pi*eps0*amu*aion*clight**3)/(zion*echarge)
      
c---------------------------------------------------------------------------
c     --- Loop over time slices of beam
      do iz=1,niz

c       --- centroid coordinates in displaced elements
        xval = y(5,iz) - pipeox(iz)
        yval = y(7,iz) - pipeoy(iz)
        rpipei2 = 1./rpipe(iz)**2
      
c       --- Beam size reciprocals
        if (y(1,iz) == 0. .or. y(3,iz) == 0.) then
          y1i = 0.
          y3i = 0.
          y1py3i = 0.
        else
          y1i = 1./y(1,iz)
          y3i = 1./y(3,iz)
          y1py3i = 1./(y(1,iz) + y(3,iz))
        endif

c       --- pathlength correction
        xbend = bendcurv(iz)*y(5,iz)

c       --- longitudinal velocity
        vz = clight*y(10,iz)
        dsdt = vz/(1.+xbend)

c       --- relativistic factor gamma and gambeta = gamma*beta
        gamma = 1/sqrt(1.0 - y(10,iz)**2)
        gambeta = y(10,iz)*gamma

c       --- magnetic rigidity reciprocal
        rigidityi = zion*echarge/(aion*amu*clight*gambeta)

c       --- quadrupole-force coefficient
        quadfld = (dedx(iz) - dbdx(iz) * vz)*rigidityi
c       --- bz ignored for now

c       --- perveance perv (labeled K in the envelope equation)
        if (lperveance) then
          perv = 2.0*y(13,iz)/(alfcur*gambeta**3)
        else
          perv = 0.
        endif

c       --- normalized emittance to unnormalized emittance
        if (lemittance) then
          xemit = y(11,iz)/gambeta
          yemit = y(12,iz)/gambeta
        else
          xemit = 0.
          yemit = 0.
        endif

        if (lallez) then
          dbetadt = (ezbeam(iz)+ez(iz))*zion*echarge/(aion*amu*clight*gamma**3)
        else
          dbetadt = 0.
        endif

        if (lviscous .and. lallez) then
          dbetadt = dbetadt + aviscous(iz)
        endif
        
c       --- bend-force terms
        bendfldx = (ex(iz)/vz - by(ii))*rigidityi
        bendfldy = (ey(iz)/vz + bx(ii))*rigidityi

        bendex = ex(iz)*rigidityi/vz
        bendey = ey(iz)*rigidityi/vz

        dlnbeta = (dbetadt/y(10,iz))*gamma**2

c       --- Now, set derivatives
c       --- da/dt
        dy(1,iz) = y(2,iz)*dsdt

c       --- d(da/ds)/dt
        dy(2,iz) = quadfld*y(1,iz)
     &           + ( 2.0*perv*y1py3i
     &           +   xemit**2*y1i**3
     &           +   fx(iz)*perv*y(1,iz)*rpipei2
     &           +   bendcurv(iz)*y(1,iz)*(bendcurv(iz) + 2.0*bendfldx)) * dsdt
     &           - dlnbeta*y(2,iz)

c       --- db/dt
        dy(3,iz) = y(4,iz)*dsdt

c       --- d(db/ds)/dt
        dy(4,iz) = -quadfld*y(3,iz)
     &           + ( 2.0*perv*y1py3i
     &           +   yemit**2*y3i**3
     &           +   fy(iz)*perv*y(3,iz)*rpipei2) * dsdt
     &           - dlnbeta*y(4,iz)

c       --- dx/ds
        dy(5,iz) = y(6,iz)*dsdt

c       --- d(dx/ds)/dt
        dy(6,iz) = quadfld*xval
     &           + ( perv*(gxx(iz)*xval + gxy(iz)*yval)*rpipei2
     &             + (1.0 + xbend)*(bendcurv(iz) + bendex)
     &             - (1.0 + 2.0*xbend)*by(iz)*rigidityi) * dsdt
     &           - dlnbeta*y(6,iz)

c       --- dy/dt
        dy(7,iz) = y(8,iz)*dsdt

c       --- d(dy/ds)/ds
        dy(8,iz) = -quadfld*yval
     &           + ( perv*(gyx(iz)*xval + gyy(iz)*yval)*rpipei2
     &             + (1.0 + xbend)*bendey
     &             + (1.0 + 2.0*xbend)*bx(iz)*rigidityi) * dsdt
     &           - dlnbeta*y(8,iz)

c       --- ds/dt
        dy(9,iz) = dsdt

c       --- dbeta/dt
        dy(10,iz) = dbetadt

      enddo
      
      return
      end
c=============================================================================
      subroutine stephr(y,dy,niz,t,dt,aion,zion,fviscous,
     &            icharge,iimage,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &            lviscous,lcurgrid,lfail)
      use Constant
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),dy(16,niz)
      real(kind=8):: t,dt,aion,zion
      real(kind=8):: fviscous
      integer(ISZ):: icharge,iimage
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail

c STEPHR  advances the set of the envelope equations using the
c isochronous leapfrog algorithm.

      real(kind=8):: halfdt, ds, dslight
      integer(ISZ):: j
      
      halfdt = 0.5*dt
      dslight = clight * dt
      
      do j=1,niz
        y( 2,j) = y( 2,j) + halfdt*dy( 2,j)
        y( 4,j) = y( 4,j) + halfdt*dy( 4,j)
        y( 6,j) = y( 6,j) + halfdt*dy( 6,j)
        y( 8,j) = y( 8,j) + halfdt*dy (8,j)
        y(10,j) = y(10,j) + halfdt*dy(10,j)
      enddo
      
      do j=1,niz
        ds = y(10,j) * dslight
        y(1,j) = y(1,j) + y(2,j) * ds
        y(3,j) = y(3,j) + y(4,j) * ds
        y(5,j) = y(5,j) + y(6,j) * ds
        y(7,j) = y(7,j) + y(8,j) * ds
        y(9,j) = y(9,j) + ds
      enddo
      
      do j=1,niz
        if (y(1,j) < 0. .or. y(3,j) < 0.) then
          call remark ("Error: Beam radius became negative")
          lfail = .true.
          return
        endif
        if ((y(1,j) .ne. y(1,j)) .or. (y(3,j) .ne. y(3,j))) then
          call remark ("Error: Beam radius became NaN")
          lfail = .true.
          return
        endif
      enddo
      
      call getderivs(t,dt,y,dy,niz,aion,zion,fviscous,icharge,iimage,
     &            lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail)
      if (lfail) return

      do j=1,niz
        y( 2,j) = y( 2,j) + halfdt*dy( 2,j)
        y( 4,j) = y( 4,j) + halfdt*dy( 4,j)
        y( 6,j) = y( 6,j) + halfdt*dy( 6,j)
        y( 8,j) = y( 8,j) + halfdt*dy (8,j)
        y(10,j) = y(10,j) + halfdt*dy(10,j)
      enddo
     
      return
      end
c=============================================================================
      subroutine setrhohr(niz,y,rpipe,icharge,lfail)
      use HERtmp
      use Constant
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),rpipe(niz)
      integer(ISZ):: icharge
      logical(ISZ):: lfail

c SETRHOHR calculates the line-charge density within the slices.
c If icharge == 7, it also loads the charge onto the RZ grid

      real(kind=8):: fact1,fact2,part1,part2
      real(kind=8):: rm,dtemp,ddtemp
      real(kind=8):: zslmin, zslmax
      integer(ISZ):: ii,ir,irm,iz,izl,izr

      if (icharge == 5 .or. icharge == 6) return
c --- denmid not needed in these cases
      
c --- first calculate line-charge density inside the slices
      denmid(niz)=0.0
      do ii=1,niz-1
        denmid(ii) = y(14,ii)/abs(y(9,ii+1) - y(9,ii))
      enddo
      
c --- load charge if icharge == 7 only
      if (icharge < 7) return

      rmmax = rpipe(1)
      zmmin = y(9,1)
      zmmax = y(9,niz)
      do ii=1,niz
        rmmax = max(rmmax,rpipe(ii))
        zmmin = min(zmmin,y(9,ii))
        zmmax = max(zmmax,y(9,ii))
      enddo
      dr = rmmax / nr
      zmmin = zmmin - 4. * rmmax
      zmmax = zmmax + 4. * rmmax
      dz = (zmmax - zmmin) / nz
      do ii = 0, nr
        rmesh(ii) = ii * dr
      enddo
      do ii = 0, nz
        zmesh(ii) = ii * dz + zmmin
      enddo
      call zeroarry(rho,nrz)
c --- Load charge
      do ii=1,niz-1
        zslmin = min (y(9,ii),y(9,ii+1))
        zslmax = max (y(9,ii),y(9,ii+1))
        izl = (zslmin-zmmin)/dz
        izr = (zslmax-zmmin)/dz
        rm = (y(1,ii)*y(3,ii) + y(1,ii+1)*y(3,ii+1)) / 2.
        dtemp = denmid(ii) / (pi * rm)
        rm = sqrt(rm)/dr
        irm = rm
        rm = rm - irm
        if (irm + 1> nr) then
          call remark ("Error in setrhohr: Beam radius is larger than the aperture.")
          lfail = .true.
          return
        endif
        fact1 = (zslmin-zmmin)/dz - izl
        fact2 = (zslmax-zmmin)/dz - izr
        if (irm > 0) then
          part1 = (irm/2. - 0.25 - rm*(irm+rm/2.)*((rm**2)/2.+rm*irm-1-2*irm)/(2.*irm+1.))/irm
          part2 = ((rm*(irm+rm/2.))**2)/((irm+1.)*(2.*irm+1.))
        else
          part1 = 4.*(rm**2)*(1-(rm**2)/2.) - 2.
          part2 = (rm**2)/4.
        endif
        do iz = izl,izr+1
          rho(0,iz) = rho(0,iz) + 2. * dtemp
          do ir=1,irm-1  
            rho(ir,iz) = rho(ir,iz) + dtemp
          enddo
          rho(irm,  iz) = rho(irm,  iz) + part1 * dtemp
          rho(irm+1,iz) = rho(irm+1,iz) + part2 * dtemp
        enddo
c ---   Subtract missing contribution to izl
        ddtemp = dtemp * (0.5+fact1-0.5*(fact1**2))
        rho(0,izl) = rho(0,izl) - 2. * ddtemp
        do ir=1,irm-1  
          rho(ir,izl) = rho(ir,izl) - ddtemp
        enddo
        rho(irm,  izl) = rho(irm,  izl) - part1 * ddtemp
        rho(irm+1,izl) = rho(irm+1,izl) - part2 * ddtemp
c ---   Subtract missing contribution to izl+1
        ddtemp = dtemp * 0.5 * (fact1**2)
        rho(0,izl+1) = rho(0,izl+1) - 2. * ddtemp
        do ir=1,irm-1  
          rho(ir,izl+1) = rho(ir,izl+1) - ddtemp
        enddo
        rho(irm,  izl+1) = rho(irm,  izl+1) - part1 * ddtemp
        rho(irm+1,izl+1) = rho(irm+1,izl+1) - part2 * ddtemp
c ---   Subtract missing contribution to izr   
        ddtemp = dtemp * 0.5 * ((1-fact2)**2)
        rho(0,izr) = rho(0,izr) - 2. * ddtemp
        do ir=1,irm-1  
          rho(ir,izr) = rho(ir,izr) - ddtemp
        enddo
        rho(irm,  izr) = rho(irm,  izr) - part1 * ddtemp
        rho(irm+1,izr) = rho(irm+1,izr) - part2 * ddtemp
c ---   Subtract missing contribution to izr+1
        ddtemp = dtemp * 0.5 * (2.-(fact2**2))       
        rho(0,izr+1) = rho(0,izr+1) - 2. * ddtemp
        do ir=1,irm-1  
          rho(ir,izr+1) = rho(ir,izr+1) - ddtemp
        enddo
        rho(irm,  izr+1) = rho(irm,  izr+1) - part1 * ddtemp
        rho(irm+1,izr+1) = rho(irm+1,izr+1) - part2 * ddtemp
      enddo

      return
      end

c=============================================================================
      subroutine fieldsolhr(niz,y,eval,rpipe,icharge,lezbeam,lezcenter,lfail)
      use Constant
      use HERtmp
      use HERbessel
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),eval(niz),rpipe(niz)
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lezcenter,lfail

c FIELDSOLHR calculates the axial space-charge field of a beam
c using various approximations.
c beam data is passed in array y, and the space-charge field
c at slice boundaries is returned in array eval.
c The field is set to zero if lezbeam if false.
c The space-charge model is determined by icharge
c   icharge = 0 uses constant g-factor space-charge model
c   icharge = 1 uses simple g-factor space-charge model
c   icharge = 2 includes envelope variation in space-charge model
c   icharge = 3 not implemented yet
c   icharge = 4 not implemented yet
c   icharge = 5 uses a Bessel-series expansion with a flat line charge density within each slice
c   icharge = 6 uses a Bessel-series expansion with a smoothly varying line charge density
c   icharge = 7 uses WARP's RZ solver

      real(kind=8):: frpieps
      real(kind=8):: del1,del2,fact1,fact2,fact3,part1,part2
      real(kind=8):: rm,phi1,phi2,phi3
      real(kind=8):: kr, arg, besselterm(niz),left,right,ddtemp
      real(kind=8):: sl,su

      integer(ISZ):: ii,ir,irm,iz

      frpieps = 4.0*pi*eps0

c calculate space-charge field at slice boundaries
c   eval = longitudinal space-charge field

      if((icharge < 0).or.(icharge>7).or.(.not. lezbeam))then
c     set space-charge field to zero

        do ii=1,niz
          eval(ii) = 0.0
        enddo

      elseif(icharge == 0 .or. icharge == 1)then
c       simple g-factor model

        do ii=1,niz
          eval(ii) = - gtemp(ii) * dden(ii) / frpieps
        enddo

      elseif(icharge == 2)then
c       general g-factor model

        do ii=1,niz
          eval(ii) = (- gtemp(ii)*dden(ii) + y(15,ii)*denv(ii)) / frpieps
        enddo
      
      elseif(icharge == 5)then
c       Bessel series expansion with a flat line charge density within each slice

        rm = 0.
        do iz=1,niz
          rm = max(rm,rpipe(iz))
        enddo
        fact1 = 2. / (eps0 * pi *(rm**2))

        do iz=1,niz
          eval(iz) = 0.
        enddo

        do ir=1,nbessel
          kr = besselzero(ir)/rm
          call getbesselterms(niz,kr,rad,besselterm)
          fact2 = fact1 * besselfactor(ir)
          deval(1) = 0.
          sl = 0.
          do iz=1,niz-1
            left = besselterm(iz)
            right = besselterm(iz+1)
            arg = kr*(y(9,iz+1)-y(9,iz))
            fact3=exp(-arg)
            ddtemp = right-left*fact3-(right-left)*(1.-fact3)/arg
            ddtemp = ddtemp*y(14,iz)/arg
            sl=fact3*sl+ddtemp
            deval(iz+1)=sl
          enddo
          su = 0.
          do iz=niz-1,1,-1
            right = besselterm(iz+1)
            left = besselterm(iz)
            arg = kr*(y(9,iz+1)-y(9,iz))
            fact3=exp(-arg)
            ddtemp = left-right*fact3+(right-left)*(1.-fact3)/arg
            ddtemp = ddtemp*y(14,iz)/arg
            su=fact3*su+ddtemp
            deval(iz)=deval(iz)-su
          enddo
          do iz=1,niz
            fact3 = fact2*besselterm(iz)
            eval(iz)=eval(iz)+fact3*deval(iz)
          enddo
        enddo

      elseif(icharge == 6)then
c       Bessel series expansion with a smoothly varying line charge density

        rm = 0.
        do iz=1,niz
          rm = max(rm,rpipe(iz))
        enddo
        fact1 = 2. / (eps0 * pi *(rm**2))

        do iz=1,niz
          eval(iz) = 0.
        enddo

        do ir=1,nbessel
          kr = besselzero(ir)/rm
          call getbesselterms(niz,kr,rad,besselterm)
          fact2 = fact1 * besselfactor(ir)
          deval(1) = 0.
          sl = 0.
          do iz=1,niz-1
            left = besselterm(iz)*y(15,iz)
            right = besselterm(iz+1)*y(15,iz+1)
            arg = kr*(y(9,iz+1)-y(9,iz))
            fact3=exp(-arg)
            ddtemp = right-left*fact3-(right-left)*(1.-fact3)/arg
            ddtemp = ddtemp/kr
            sl=fact3*sl+ddtemp
            deval(iz+1)=sl
          enddo
          su = 0.
          do iz=niz-1,1,-1
            right = besselterm(iz+1)*y(15,iz+1)
            left = besselterm(iz)*y(15,iz)
            arg = kr*(y(9,iz+1)-y(9,iz))
            fact3=exp(-arg)
            ddtemp = left-right*fact3+(right-left)*(1.-fact3)/arg
            ddtemp = ddtemp/kr
            su=fact3*su+ddtemp
            deval(iz)=deval(iz)-su
          enddo
          do iz=1,niz
            fact3 = fact2*besselterm(iz)
            eval(iz)=eval(iz)+fact3*deval(iz)
          enddo
        enddo

      elseif(icharge == 7)then
c       use Warp's RZ package to find the field

c ---   The charge has been loaded already using setrhohr
c ---   Calculate potential
        call vprz (1)
        call fieldsolrz (0)
c ---   Find Ez on slices
        if (lezcenter) then
          do ii=1,niz
            iz = 0.5 + (y(9,ii)-zmmin)/dz
            del2 = (iz+1)*dz-y(9,ii)+zmmin
            del1 = 2*dz-del2
            fact1 =  del2-dz/2.
            fact2 =  del1-del2
            fact3 = -del1+dz/2.
            eval(ii) = (fact1*phi(0,iz-1) + fact2*phi(0,iz) + fact3*phi(0,iz+1)) / dz**2
          enddo
        else
          do ii=1,niz
            iz = 0.5 + (y(9,ii)-zmmin)/dz
            del2 = (iz+1)*dz-y(9,ii)+zmmin
            del1 = 2*dz-del2
            fact1 =  del2-dz/2.
            fact2 =  del1-del2
            fact3 = -del1+dz/2.
            rm = sqrt(y(1,ii)*y(3,ii))/dr
            if (rm > 0.) then
              irm = rm
              rm = rm - irm
              part1 = irm - 0.5 - 2.*rm*(irm+rm/2.)*((rm**2)/2.+rm*irm-1-2*irm)/(2.*irm+1.)
              part2 = 2.*((rm*(irm+rm/2.))**2)/(2.*irm+1.)
              phi1 = 0.5*phi(0,iz-1)
              phi2 = 0.5*phi(0,iz)
              phi3 = 0.5*phi(0,iz+1)
              do ir = 1,irm-1
                phi1 = phi1 + 2.*ir*phi(ir,iz-1)
                phi2 = phi2 + 2.*ir*phi(ir,iz)
                phi3 = phi3 + 2.*ir*phi(ir,iz+1)
              enddo
              phi1 = phi1 + part1*phi(irm,iz-1) + part2*phi(irm+1,iz-1)
              phi2 = phi2 + part1*phi(irm,iz)   + part2*phi(irm+1,iz)
              phi3 = phi3 + part1*phi(irm,iz+1) + part2*phi(irm+1,iz+1)
              eval(ii) = (fact1*phi1 + fact2*phi2 + fact3*phi3) / (((irm+rm)*dz)**2)
            else
              phi1 = phi(0,iz-1)
              phi2 = phi(0,iz)
              phi3 = phi(0,iz+1)
              eval(ii) = (fact1*phi1 + fact2*phi2 + fact3*phi3) / (dz**2)
            endif
          enddo
        endif
      endif
      
      return
      end
c=============================================================================

      subroutine getbesselterms(niz,kr,rad,besselterm)
      integer(ISZ):: niz
      real(kind=8):: kr
      real(kind=8):: rad(niz)
      real(kind=8):: besselterm(niz)

c     getbesselterms calculates the Bessel function J1 at x,
c     (using matched analytic approximations),
c     divided by x:
c     besselterm(x) = J1(x)/x
c     in which x = kr*rad
c     The routine is based on bessj1 from Numerical Recipes.

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

      real(kind=8):: y,z,x,xx
      integer(ISZ):: ii
      
      do ii=1,niz
        x = kr*rad(ii)
        if(x < 8.0)then
          y = x**2
          besselterm(ii) = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     .                    /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
        else
          z = 8./x
          y = z**2
          xx = x - 2.356194491
          besselterm(ii) = sqrt(0.636619772/x)*
     .                      (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))
     .                      -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))/x
        endif
      enddo

      return
      end
c=============================================================================
      subroutine getradius(y,niz,rpipe,icharge)
      use HERtmp
      use Constant
      integer(ISZ):: niz,icharge
      real(kind=8):: y(16,niz), rpipe(niz)

      integer(ISZ):: ii
      real(kind=8):: rmin,aval,bval

c --- calculate the beam radius as needed to find Ez
      if (icharge < 5) then
        do ii=1,niz
          rmin = 0.01*rpipe(ii)
          aval = y(1,ii)
          bval = y(3,ii)
          aval = min(max(rmin,aval),rpipe(ii))
          bval = min(max(rmin,bval),rpipe(ii))
          rad(ii) = sqrt(aval*bval)
        enddo
      else
        do ii=1,niz
          aval = y(1,ii)
          bval = y(3,ii)
          rad(ii) = sqrt(aval*bval)
        enddo
      endif
      
      return
      end
c=============================================================================
      subroutine getgfactor(niz,rpipe,icharge)
      use HERtmp
      integer(ISZ):: niz,icharge
      real(kind=8):: rpipe(niz)

      integer(ISZ):: ii
      
c --- calculate space-charge geometric factor g
c --- For icharge == 0, the g-factor has been set by the generate
      if (icharge > 1) then
        do ii=1,niz
          gtemp(ii) = 2.*log(rpipe(ii)/rad(ii)) + 0.5
        enddo
      else
        do ii=1,niz
          gtemp(ii) = 2.*log(rpipe(ii)/rad(ii))
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getcurrent(y,niz,rpipe,icharge,lcurgrid,lfail)
      use HERtmp
      use Constant
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      integer(ISZ):: niz,icharge
      real(kind=8):: y(16,niz), rpipe(niz)
      logical(ISZ):: lcurgrid,lfail

      integer(ISZ):: ii, iz, ir
      real(kind=8):: del1, del2, part1, part2, fact1, fact2
      real(kind=8):: dlna, dlnb
      real(kind=8):: zb, denleft, denright

c --- calculate beam quantities on slice boundaries
c ---  cur = beam current cur
c ---  den = line-charge density
c ---  denv = envelope-variation term 0.5*den*(d(ab)/dt)/(ab*clight)

      if (icharge == 7 .and. lcurgrid) then
      
c ---   Interpolate the line charge density from the grid
        fact1 = 2.*pi*dr*dr
        do ii=1, niz
          zb = (y(9,ii)-zmmin)/dz
          iz = zb
          denleft  = 0.125*rho(0,iz)
          denright = 0.125*rho(0,iz+1)
          do ir = 1,nr
            denleft  = denleft  + ir*rho(ir,iz)
            denright = denright + ir*rho(ir,iz+1)
          enddo
          part2 = zb - iz
          part1 = 1. - part2
          y(15,ii) = fact1 * (part1 * denleft + part2 * denright)
          y(13,ii) = y(15,ii) * y(10,ii) * clight
        enddo
      
      else

c ---   values at beam head
        y(13,1) = 0.0
        y(15,1) = 0.0
        denv(1) = 0.0
c ---   values at beam tail
        y(13,niz) = 0.0
        y(15,niz) = 0.0
        denv(niz) = 0.0
c ---   values in beam body
        do ii=2,niz-1
          del1 = y(9,ii) - y(9,ii-1)
          del2 = y(9,ii+1) - y(9,ii)
          part1 = del1/(del1 + del2)
          fact1 = del1/(del2*(del1 + del2))
          fact2 = del2/(del1*(del1 + del2))
          if (del1 < 0.) then
            print *, "Slice overtaking has occurred at s = ", y(9,ii-1)
            print *, "Slice ",ii-2," overtook slice ",ii-1
c ---       following Python-convention for array indexing
            lfail = .true.
            return
          endif
          if (del1 .ne. del1) then
            print *, "At s = ", y(9,ii-1)
            print *, "Position of slice", ii-1, " became NaN"
            lfail = .true.
            return
          endif
          dlna = (fact1*(y(1,ii+1)-y(1,ii)) + fact2*(y(1,ii)-y(1,ii-1)))/y(1,ii)
          dlnb = (fact1*(y(3,ii+1)-y(3,ii)) + fact2*(y(3,ii)-y(3,ii-1)))/y(3,ii)
          y(15,ii) = (1.0-part1)*denmid(ii-1) + part1*denmid(ii)
          y(13,ii) = y(15,ii) * y(10,ii) *clight
          denv(ii) = 0.5*y(15,ii)*(dlna + dlnb)
        enddo
        if (del2 < 0.) then
          print *, "Slice overtaking has occurred at s = ", y(9,niz-1)
          print *, "Slice ",niz-2," overtook slice ",niz-1
c ---     following Python-convention for array indexing
          lfail = .true.
          return
        endif
        if (del2 .ne. del2) then
          print *, "At s = ", y(9,niz-1)
          print *, "Position of slice ", niz-1, " became NaN"
          lfail = .true.
          return
        endif
      endif

c --- calculate partial z-derivative of the line-charge density d(den)/dz 
      if (icharge < 5) then
c ---   don't bother otherwise
        del1 = y(9,2) - y(9,1)
        dden(1) = 2*y(14,1)/del1**2
        del2 = y(9,niz) - y(9,niz-1)
        dden(niz) = -2*y(14,niz-1)/del2**2
        do ii=2,niz-1
          dden(ii) = 2.0*(denmid(ii) - denmid(ii-1)) /( y(9,ii+1) - y(9,ii-1))
        enddo
      endif
      
      return
      end
c=============================================================================
      subroutine getviscous(y,niz,aviscous,fviscous)
      use Constant
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),aviscous(niz),fviscous

c add artificial viscosity to the Ez field

      integer(ISZ):: iz
      real(kind=8):: del1, del2, fact
      real(kind=8):: dzviscous
      real(kind=8):: dbetadz(niz)
      real(kind=8):: temp(niz)

      dzviscous = fviscous * (y(9,niz)-y(9,1)) 
      del1 = y(9,2) - y(9,1)
      del2 = y(9,niz) - y(9,niz-1)
      dbetadz(1)   = (y(10,2)  -y(10,1))    /del1
      dbetadz(niz) = (y(10,niz)-y(10,niz-1))/del2
      do iz=2,niz-1
        del1 = y(9,iz) - y(9,iz-1)
        del2 = y(9,iz+1) - y(9,iz)
        dbetadz(iz) = ((del1/del2)*y(10,iz+1)-(del2/del1)*y(10,iz-1))/(del1+del2)
     &              + (1/del1-1/del2)*y(10,iz)
      enddo
      temp(1)   = 0.
      temp(niz) = 0.
      do iz=2,niz-1
        temp(iz) = (dbetadz(iz)**2)*y(1,iz)*y(3,iz)/y(15,iz)
c       because this is what we actually need
      enddo
      fact = - clight*(dzviscous*y(15,niz/2)/(y(1,niz/2)*y(3,niz/2)))**2
      aviscous(1)   = 0.
      aviscous(niz) = 0.
      do iz=2,niz-1
        del1 = y(9,iz) - y(9,iz-1)
        del2 = y(9,iz+1) - y(9,iz)
        if (dbetadz(iz)*y(10,iz) < 0.) then
          aviscous(iz) = (((del1/del2)*temp(iz+1)
     &                 - (del2/del1)*temp(iz-1))/(del1+del2)
     &                + (1/del1-1/del2)*temp(iz)) * fact * y(1,iz) * y(3,iz) / y(15,iz)
        else
          aviscous(iz) = 0.
        endif
      enddo
      
      return
      end
c=============================================================================
      subroutine savehermesvars(t,y,niz,it,lsavehist)
      use HERhist
      integer(ISZ):: niz,it
      real(kind=8):: t,y(16,niz)
      logical(ISZ):: lsavehist

      integer(ISZ):: iz

c     --- Check if saving history
      if (lsavehist) then
c       --- Check if saving history this step. Requires (nhher > 0).
        if (mod(it,nhher) == 0) then
          if (jhher >= lhher) then
            lhher = lhher + 1000
            call gchange("HERhist",0)
          endif
          jhher = jhher + 1
          hther(jhher) = t
          do iz=1,niz
            haher(iz,jhher)   = y( 1,iz)
            hapher(iz,jhher)  = y( 2,iz)
            hbher(iz,jhher)   = y( 3,iz)
            hbpher(iz,jhher)  = y( 4,iz)
            hxher(iz,jhher)   = y( 5,iz)
            hxpher(iz,jhher)  = y( 6,iz)
            hyher(iz,jhher)   = y( 7,iz)
            hypher(iz,jhher)  = y( 8,iz)
            hsher(iz,jhher)   = y( 9,iz)
            hvzher(iz,jhher)  = y(10,iz)
            henxher(iz,jhher) = y(11,iz)
            henyher(iz,jhher) = y(12,iz)
            hcur(iz,jhher)    = y(13,iz)
            hdq(iz,jhher)     = y(14,iz)
            hden(iz,jhher)    = y(15,iz)
          enddo
        endif
      endif

      return
      end

c=============================================================================

      subroutine resizehermeshist()
      use HERhist
      use HERflags

c  At the end of a run, resizes the history arrays

c     --- Resize history arrays
      if (lsavehist .and. lhher > jhher) then
        lhher = jhher
        call gchange("HERhist",0)
      endif

      return
      end

c=============================================================================

      subroutine getimage(y,niz,iimage,lfail)
      use Constant
      use HERfield
      integer(ISZ):: niz
      real(kind=8):: y(16,niz)
      integer(ISZ):: iimage
      logical(ISZ):: lfail

c GETIMAGE calculates image-field coefficients for the different
c slice boundaries for the pipe type ipipetype.
c   0 round pipe (default)
c   1 electric quadrupole rods (hyperbolic sufaces)
c   2 dipole plates (parallel plates)
c Current values of a, b, X, and Y from 'y' are used.
c The assumed forms for the image-field components scaled by
c rigidity*beta*clight*gamma**2 are
c   Ex = (K/R**2)*(Fxx*(x-X) + Fxy*(y-Y) + Gxx*X + Gxy*Y)
c   Ey = (K/R**2)*(Fyx*(x-X) + Fyy*(y-Y) + Gyx*X + Gyy*Y).
c Fxy and Fyx do not to contribute to the averaged equations, and
c the remaining coefficients are stored in common block HERtmp as...
c   fx(ii) = Fxx
c   fy(ii) = Fyy
c   gxx(ii) = Gxx
c   gxy(ii) = Gxy
c   gyx(ii) = Gyx
c   gyy(ii) = Gyy
c
c The integer flag iimage determines to what order image effects are included
c If iimage == 0, then image effects are neglected
c If iimage == 1, only linear terms are included. Terms of the order of 
c                 (a,b,X,Y)**2/rpipe**2 in the image coefficients are discarded.
c If iimage == 2, these terms are retained.
c
c Image coefficients are zero if the centroid displacement
c X**2 + Y**2 exceeds rpipe**2.
c
c This subroutine was copied from Circe and modified for Hermes.
c It has not been checked yet.

      integer(ISZ):: iz
      real(kind=8):: logcon
      real(kind=8):: pi2,cosrot,sinrot,cosrot2,sinrot2
      real(kind=8):: xval,yval,xyrpipe,abrpipe,geomfac
      real(kind=8):: f1,f2,g1,g2
      real(kind=8):: dispx,dispy
      real(kind=8):: c0,c1,c2,d0,d1
      real(kind=8):: arg0,arg1,arg2

c     --- Check if images are to be applied. If not, zero arrays and return.
      if ((iimage <= 0) .or. (iimage > 2) 
     &  .or. (ipipetype < 0) .or. (ipipetype > 2)) then
        do iz=1,niz
          fx(iz) = 0.0
          fy(iz) = 0.0
          gxx(iz) = 0.0
          gxy(iz) = 0.0
          gyx(iz) = 0.0
          gyy(iz) = 0.0
        enddo
        return
      endif

c     --- Constants
      pi2 = pi**2

c     --- Setup image coefficients

      if (ipipetype == 0) then
c       --- set image coefficients for circular pipe
        do iz=1,niz
          xval = y(5,iz) - pipeox(iz)
          yval = y(7,iz) - pipeoy(iz)
          arg0 = (y(1,iz)**2-y(3,iz)**2)/(4.0*rpipe(iz)**2)
          arg1 = (xval**2 - yval**2)/rpipe(iz)**2
          arg2 = (xval**2 + yval**2)/rpipe(iz)**2
c         --- turn off image force if X**2+Y**2 >= rpipe**2
          displace = arg2
          if (displace > 1.0) then
            fx(iz)  = 0.0
            fy(iz)  = 0.0
            gxx(iz) = 0.0
            gxy(iz) = 0.0
            gyx(iz) = 0.0
            gyy(iz) = 0.0
          else
            if (iimage == 1) then
              fx(iz)  = 0.
              fy(iz)  = 0.
              gxx(iz) = 1.
              gxy(iz) = 0.
              gyx(iz) = 0.
              gyy(iz) = 1.
            elseif (iimage == 2) then
              f1 = arg1 + arg0*(1.0 + 6.0*arg2 + 6.0*arg0*arg1)
              g1 = arg1 + arg0*(1.0 + 3.0*arg2 + 2.0*arg0*arg1)
              g2 = 2.0*(1.0-2.0*arg0**2)*xval*yval/rpipe(iz)**2
              fx(iz)  = f1
              fy(iz)  = -f1
              gxx(iz) = 1.0 + g1
              gxy(iz) = g2
              gyx(iz) = g2
              gyy(iz) = 1.0 - g1
            endif
          endif
        enddo
       
      elseif (ipipetype == 1) then
c     --- Set image coefficients for electric quadrupoles
c     --- Laslett formulation is used here

        do iz=1,niz
          if ((y(1,iz) <= 0.0).or.(y(3,iz) <= 0.0)) then
            call remark ("ERROR: In getimage: beam radius less than zero")
            lfail = .true.
            return
          endif
          
          xval = y(5,iz) - pipeox(iz)
          yval = y(7,iz) - pipeox(iz)
          dispx = xval**2/rpipe(iz)**2
          dispy = yval**2/rpipe(iz)**2

          displace = dispx + dispy
          if (displace > 1.0) then
            fx(iz) = 0.0
            fy(iz) = 0.0
            gxx(iz) = 0.0
            gxy(iz) = 0.0
            gyx(iz) = 0.0
            gyy(iz) = 0.0
          else
            logcon = y(1,iz)*y(3,iz)*log10(y(1,iz)/y(3,iz))/rpipe(iz)**2
            c0 = 0.785 + 3.0*logcon
            d0 = 0.785 - 3.0*logcon
            c1 = 0.834 + 3.0*logcon
            d1 = 0.834 - 3.0*logcon
            c2 = - (2.486 + (6.5*logcon)**2)
            fx(iz) = 0.0
            fy(iz) = 0.0
            gxx(iz) = c0 + 2.0*c1*dispx + c2*dispy
            gxy(iz) = 0.0
            gyx(iz) = 0.0
            gyy(iz) = d0 + 2.0*d1*dispy + c2*dispx
          endif
        enddo

      elseif (ipipetype == 2) then
c       --- Image forces for dipoles (parallel plates)
        do iz=1,niz
          xval = y(5,iz) - pipeox(iz)
          yval = y(7,iz) - pipeoy(iz)
          if (xval**2 + yval**2 > rpipe(iz)**2) then
            fx(iz) = 0.
            fy(iz) = 0.
            gxx(iz) = 0.
            gxy(iz) = 0.
            gyx(iz) = 0.
            gyy(iz) = 0.
          else
            geomfac = pi2/24.0
            cosrot = cos(pipeph(iz))
            sinrot = sin(pipeph(iz))
            cosrot2 = cosrot**2
            sinrot2 = sinrot**2
            if (iimage == 1) then
              fx(iz) = geomfac*cosrot2
              fy(iz) = geomfac*sinrot2
              gxx(iz) = 3.0*geomfac*cosrot2
              gxy(iz) = 3.0*geomfac*sinrot*cosrot
              gyx(iz) = 3.0*geomfac*sinrot*cosrot
              gyy(iz) = 3.0*geomfac*sinrot2
            elseif (iimage == 2) then
              xyrpipe = ((xval*cosrot + yval*sinrot)/rpipe(iz))**2
              abrpipe = (y(1,iz)**2 - y(3,iz)**2)/rpipe(iz)**2
              g1 = 1.0 - (pi2/12.0)*(xyrpipe + 0.1875*abrpipe*(cosrot2 - sinrot2))
              f1 = 1.0 - 0.375*pi2*xyrpipe
              f2 = 0.025*pi2*abrpipe
              fx(iz) = geomfac*cosrot2*(f1 - f2*(1.0 - 0.25*sinrot2))
              fy(iz) = geomfac*sinrot2*(f1 + f2*(1.0 - 0.25*cosrot2))
              gxx(iz) = 3.0*geomfac*cosrot2*g1
              gxy(iz) = 3.0*geomfac*sinrot*cosrot*g1
              gyx(iz) = 3.0*geomfac*sinrot*cosrot*g1
              gyy(iz) = 3.0*geomfac*sinrot2*g1
            endif
          endif
        enddo

      endif

      return
      end
c=============================================================================

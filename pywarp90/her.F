#include "top.h"
c=============================================================================
c@(#) File HER.M, version $Revision: 3.4 $, $Date: 2001/06/26 17:29:04 $
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

c     --- allocate dynamic arrays ("group change")
      if (niz == 0) niz = 101
      call gchange("HERvars",0)
      niztmp = niz
      call gchange("HERtmp",0)
      nizfield = niz
      call gchange("HERfield",0)
      lhher = 100
      if (nizhist == 0) nizhist = niz
      call gchange("HERhist",0)

c     --- Allocate temp space
      nizvartmp = niz
      call gchange("HERvartmp",0)
      nizfield = niz
      call gchange("HERfield",0)
c     --- reset Lattice arrays to final size
      call resetlat

c     --- Initialize the LatticeInternal arrays
      nzl = niz + 2
      nzlmax = nzl
      dzl = (zimax - zimin)/niz
      zlmin = zimin - dzl
      zlmax = zimax + dzl
      dzli = 1./dzl
      call gchange("LatticeInternal",0)
      call setlattzt(0.,0.,0)

c     if (icharge == 0) then
c --- Set the g-factor for the constant g-factor model
c       do ii=1,niz
c         gtemp(ii) = 2*log(1.6)
c       enddo
c     endif
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
      use HERvars
      use HERgeneral
      use HERhist
      use HERflags

c  Invoked by the STEP command, it runs the HERMES code.

      real(kind=8):: wtimeoff
      
c  Send a startup message to the standard output

      if (lherout) then
        call remark(" ***  HERMES running")
      endif

c  Read in Bessel data if needed
      if (icharge == 5 .or. icharge == 6) then
        call readbessel
      endif
c  Execute envelope solve, and time it
      call wtimeon
      call herx
      hertime = wtimeoff()
 
c     --- Resize history arrays
      if (lsavehist .and. lhher > jhher) then
        lhher = jhher
        call gchange("HERhist",0)
      endif

      return
      end
c=============================================================================
      subroutine herfin()

c  Deallocates dynamic storage used by test driver

      call gfree("HERvars")
      call gfree("HERtmp")
      call gfree("HERfield")
      call gfree("HERhist")

      return
      end
c=============================================================================
      subroutine herx()
      use Constant
      use Beam_acc
      use Picglb
      use InGen
      use InPart
      use Lattice
      use LatticeInternal
      use Mult_data
      use HERgeneral
      use HERvars
      use HERvartmp
      use HERflags
      use HERhist

c  For convenience, this routine accesses the lattice database and is not
c  linked to the universe solely by its calling sequence.

      integer(ISZ):: ii,il
      logical(ISZ):: lfail

c     --- Make sure that the sign of the step size and the relation
c     --- between zlher and zuher are consistent.
      if (dther > 0. .and. zuher < zlher) then
        call remark("Error: Either change the sign of dther, or make zlher < zuher")
        return
      endif
      if (dther < 0. .and. zlher < zuher) then
        call remark("Error: Either change the sign of dther, or make zlher > zuher")
        return
      endif
      
c     --- Make history arrays big enough
      if (lsavehist) then
        il = max(100,100*int((zuher - zlher)/(vbeam*dther)/max(1,nhher)/100. + 1))
        il = max(lhher,il)
        if (nizhist == 0) nizhist = niz
        if (il > lhher) then
          lhher = il
          call gchange("HERhist",0)
        endif
c       jhher = -1
      endif

c     --- Main loop of envelope code
c     --- Run the HERMES kernel
      call herrun(niz,var,dther,nther,aion,zion,zlher,zuher,fviscous,islice,
     &            nsteps,icharge,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &            lviscous,lsavehist,lcurgrid,lfail)
     
c     --- Copy the data into separated arrays.
      do ii=1,niz
        aher(ii)   = var( 1,ii)
        apher(ii)  = var( 2,ii)
        bher(ii)   = var( 3,ii)
        bpher(ii)  = var( 4,ii)
        xher(ii)   = var( 5,ii)
        xpher(ii)  = var( 6,ii)
        yher(ii)   = var( 7,ii)
        ypher(ii)  = var( 8,ii)
        sher(ii)   = var( 9,ii)
        vzher(ii)  = var(10,ii)
        enxher(ii) = var(11,ii)
        enyher(ii) = var(12,ii)
        cur(ii)    = var(13,ii)
        dq(ii)     = var(14,ii)
        den(ii)    = var(15,ii)
      enddo

      return
      end
c=============================================================================
      subroutine herrun(niz,y,dther,nther,aion,zion,zlher,zuher,fviscous,
     &            islice,nsteps,icharge,lezbeam,lperveance,lemittance,lallez,
     &            lezcenter,lviscous,lsavehist,lcurgrid,lfail)
      use Constant
      use Picglb
      use LatticeInternal
      use HERvartmp
      use HERfield
      use HERhist
      integer(ISZ):: niz,nther
      real(kind=8):: y(16,niz)
      real(kind=8):: dther,aion,zion,zlher,zuher
      real(kind=8):: fviscous
      integer(ISZ):: icharge,islice,nsteps
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid
      logical(ISZ):: lsavehist,lfail

c Run the HERMES kernel, using the inputted array, y and the starting data.

      real(kind=8):: t
      integer(ISZ):: i
      
      t = 0.
c     --- Start at time = 0
      lfail = .false.

c  Set fundamental constants, derived quantities, initial conditions

c     --- Save the initial history
      call savehisther(t,y,niz,nther,lsavehist)
      
c     --- Find initial forces
      call getderivs(t,dther,y,dy,niz,aion,zion,fviscous,icharge,
     &            lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail)
     
c    --- Check that the distance traveled is less than the final distance
c    --- to travel.

      if (nsteps==0) then

        do while ((y(9,islice+1)<=zuher .and. y(10,islice+1)*dther > 0.) 
     &    .or. (y(9,islice+1)>=zuher .and. y(10,islice+1)*dther < 0.))

          nther = nther+1
c         --- Advance the envelope and orbits;
          call onestep (y,dy,niz,t,dther,aion,zion,fviscous,
     &              icharge,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &              lviscous,lcurgrid,lfail)
          if (lfail) exit
          t = t + dther
          call savehisther(t,y,niz,nther,lsavehist)
        end do

      else
      
        do i=1,nsteps

          nther = nther+1
c           --- Advance the envelope and orbits;
          call onestep (y,dy,niz,t,dther,aion,zion,fviscous,
     &              icharge,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &              lviscous,lcurgrid,lfail)
          if (lfail) exit
          t = t + dther
          call savehisther(t,y,niz,nther,lsavehist)
        enddo
         
      endif
      
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
c  multipoles.  Only linear forces consistent with a truncated second-order 
c  envelope formulation are calculated.  These include a linear applied 
c  forces (electric and magnetic quadrupole fields, and a radial 
c  E-field proportional to the radial coordinate r), and 
c  linear space-charge forces.  
c  Input:
c     t    is the starting time of the current time step
c     dt   is the size of the current time step
c     y    is the array hold the envelope information
c     niz  is the number of slices

      real(kind=8):: vz,vzi,zl,zr,fl,fr,frac,dta
      integer(ISZ):: io,iele,ii,j
      
      dta = abs(dt)
      
      zlmin = y(9,1)
      zlmax = y(9,niz)
      do ii=1,niz
        zlmin = min (zlmin,y(9,ii))
        zlmax = max (zlmax,y(9,ii))
      enddo
      dzl = (zlmax - zlmin)/niz
      zlmin = zlmin - dzl
      zlmax = zlmax + dzl
      dzli = 1./dzl
      
      call setlattzt(0.,t,0)

      if (rwall > 0.) then
        do ii=1,niz
c         --- Default value of rpipe
           rpipe(ii) = rwall
        enddo
      endif

      do ii=1,niz
c       --- Zero out the fields
        dedx(ii) = 0.
        dbdx(ii) = 0.
      enddo

c---------------------------------------------------------------------------
c       --- Hard edge quadrupoles

      if (quads) then
c       --- If in a quad, get field and pipe.

        do io=1,nquadol
          do ii=1,niz
            vz = y(10,ii)*clight
            vzi = 1./dvnz(vz)
c           --- find z-cell in which particle lies
            j = max(0., (y(9,ii) - zlmin) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = y(9,ii) - vz * dta / 2.
            fl = 0.
            if (zl.ge.cquadzs(j,io) .and. zl.lt.cquadze(j,io)) then
              fl = 1.
              iele = cquadid(j,io)
              if (quadap(iele) > 0.) rpipe(ii) = quadap(iele)
            endif
c           --- "right" end of velocity advance step
            zr = y(9,ii) + vz * dta / 2.
            fr = 0.
            if (zr.ge.cquadzs(j,io) .and. zr.lt.cquadze(j,io)) fr = 1.
c            --- residence fraction
            frac = fl
            if (fl .gt. fr) frac = (cquadze(j,io)-zl)*vzi/dta
            if (fr .gt. fl) frac = (zr-cquadzs(j,io))*vzi/dta
c           --- set the field
            if (frac .gt. 0.) then
              dbdx(ii) = dbdx(ii) + cquaddb(j,io) * frac
              dedx(ii) = dedx(ii) + cquadde(j,io) * frac
            endif
          enddo
        enddo
      endif

c---------------------------------------------------------------------------
c       --- Get aperture from drift elements
      if (drfts) then
        
        do io=1,ndrftol
          do ii=1,niz
c           --- find z-cell in which particle lies
            j = max(0., (y(9,ii) - zlmin) * dzli + 0.5)
c           --- use current position to find the aperture
            zl = y(9,ii)
            if (zl.ge.cdrftzs(j,io) .and. zl.lt.cdrftze(j,io)) then
              iele = cdrftid(j,io)
              if (drftap(iele) > 0.) rpipe(ii) = drftap(iele)
            endif
          enddo
        enddo
        
      endif

c---------------------------------------------------------------------------

      do ii=1,niz
c ---   Print warning if rpipe is zero
        if (rpipe(ii) == 0.) then
          call remark("WARNING: rpipe is zero, make sure that either top.rwall")
          call remark("and/or the lattice element apertures are set.")
        endif
        
      enddo

      return
      end
c=============================================================================
      subroutine getderivs(t,dt,y,dy,niz,aion,zion,fviscous,
     &                  icharge,lezbeam,lperveance,lemittance,lallez,
     &                  lezcenter,lviscous,lcurgrid,lfail)
      use Constant
      use HERfield
      integer(ISZ):: niz
      real(kind=8):: t,dt,y(16,niz),dy(16,niz)
      real(kind=8):: aion,zion
      real(kind=8):: fviscous
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail

c  Computes the forces needed for the general-lattice envelope calculation.
c  Includes both residence corrected hard edged quadrupole elements 
c  and hard-edge multipole elements, and elements defined by their 
c  multipoles.  Only linear forces consistent with a truncated second-order 
c  envelope formulation are calculated.  These include a linear applied 
c  forces (electric and magnetic quadrupole fields, and a radial 
c  E-field proportional to the radial coordinate r), and 
c  linear space-charge forces.  

      integer(ISZ):: ii
      real(kind=8):: alfcur
      real(kind=8):: y1i,y3i,y1py3i,gamma,gambeta
      real(kind=8):: rigidityi,quadfld,perv
      real(kind=8):: xemit,yemit,dsdt,dlnbeta,dbetadt
      real(kind=8):: rpipei2
      real(kind=8):: aviscous(niz)

c     --- Get the applied fields.
      call extebher(t,dt,y,niz)
      
c     --- Image forces are not yet implemented in HERMES

c     --- Calculate the beam radius if needed
      if (icharge < 7) call getradius(y,niz,rpipe,icharge)
      if (icharge < 5) call getgfactor(niz,rpipe,icharge)
      
c     --- Calculate auxiliary variables needed for getselffield and getcurrent
      call loadcharge (niz,y,rpipe,icharge,lfail)
      if (lfail) return

c     --- Calculate the current (needed for the envelope equation)
      call getcurrent(y,niz,rpipe,icharge,lcurgrid,lfail)
      if (lfail) return
      
c     --- Get self axial fields of beam.
      call getselffield (niz,y,ezbeam,rpipe,icharge,lezbeam,lezcenter,lfail)
      if (lfail) return
      
      if (lviscous .and. lezbeam) then
        call getviscous(y,niz,aviscous,fviscous)
      endif

c     --- Set some temporaries.
      alfcur = (4.*pi*eps0*amu*aion*clight**3)/(zion*echarge)
      
c---------------------------------------------------------------------------
c     --- Loop over time slices of beam
      do ii=1,niz

        rpipei2 = 1./rpipe(ii)**2
      
c       --- Beam size reciprocals
        if (y(1,ii) == 0. .or. y(3,ii) == 0.) then
          y1i = 0.
          y3i = 0.
          y1py3i = 0.
        else
          y1i = 1./y(1,ii)
          y3i = 1./y(3,ii)
          y1py3i = 1./(y(1,ii) + y(3,ii))
        endif

c       --- longitudinal velocity
        dsdt = clight*y(10,ii)

c       --- relativistic factor gamma and gambeta = gamma*beta
        gamma = 1/sqrt(1.0 - y(10,ii)**2)
        gambeta = y(10,ii)*gamma

c       --- magnetic rigidity reciprocal
        rigidityi = zion*echarge/(aion*amu*clight*gambeta)

c       --- quadrupole-force coefficient
        quadfld = ( - dedx(ii)/dsdt + dbdx(ii))*rigidityi

c       --- perveance perv (labeled K in the envelope equation)
        if (lperveance) then
          perv = 2.0*y(13,ii)/(alfcur*gambeta**3)
        else
          perv = 0.
        endif

c       --- normalized emittance to unnormalized emittance
        if (lemittance) then
          xemit = y(11,ii)/gambeta
          yemit = y(12,ii)/gambeta
        else
          xemit = 0.
          yemit = 0.
        endif

        if (lallez) then
          dbetadt = ezbeam(ii)*zion*echarge/(aion*amu*clight*gamma**3)
        else
          dbetadt = 0.
        endif

        if (lviscous .and. lallez) then
          dbetadt = dbetadt + aviscous(ii)
        endif
        
        dlnbeta = (dbetadt/y(10,ii))*gamma**2

c       --- Now, set derivatives
c       --- da/dt
        dy(1,ii) = y(2,ii)*dsdt

c       --- d(da/ds)/dt
        dy(2,ii) = (- quadfld*y(1,ii) + 2.0*perv*y1py3i + xemit**2*y1i**3) * dsdt - dlnbeta*y(2,ii)

c       --- db/dt
        dy(3,ii) = y(4,ii)*dsdt

c       --- d(db/ds)/dt
        dy(4,ii) = (quadfld*y(3,ii) + 2.0*perv*y1py3i + yemit**2*y3i**3) * dsdt - dlnbeta*y(4,ii)

c       --- dx/ds
        dy(5,ii) = y(6,ii)*dsdt

c       --- d(dx/ds)/dt
        dy(6,ii) = 0.
c       --- Centroid motion is not implemented in HERMES

c       --- dy/dt
        dy(7,ii) = y(8,ii)*dsdt

c       --- d(dy/ds)/ds
        dy(8,ii) = 0.

c       --- ds/dt
        dy(9,ii) = dsdt

c       --- dbeta/dt
        dy(10,ii) = dbetadt

      enddo
      
      return
      end
c=============================================================================
      subroutine onestep(y,dy,niz,t,dt,aion,zion,fviscous,
     &            icharge,lezbeam,lperveance,lemittance,lallez,lezcenter,
     &            lviscous,lcurgrid,lfail)
      use Constant
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),dy(16,niz)
      real(kind=8):: t,dt,aion,zion
      real(kind=8):: fviscous
      integer(ISZ):: icharge
      logical(ISZ):: lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail

c ONESTEP  advances the set of the envelope equations using the
c isochronous leapfrog algorithm.

      real(kind=8):: halfdt, ds, dslight
      integer(ISZ):: j
      
     
      halfdt = 0.5*dt
      dslight = clight * dt
      
      do j=1,niz
        y( 2,j) = y( 2,j) + halfdt*dy( 2,j)
        y( 4,j) = y( 4,j) + halfdt*dy( 4,j)
        y(10,j) = y(10,j) + halfdt*dy(10,j)
      enddo
      
      do j=1,niz
        ds = y(10,j) * dslight
        y(1,j) = y(1,j) + y(2,j) * ds
        y(3,j) = y(3,j) + y(4,j) * ds
        y(9,j) = y(9,j) + ds
      enddo
      
      do j=1,niz
        if (y(1,j) < 0. .or. y(3,j) < 0.) then
          print *, "Error: Beam radius became negative"
          lfail = .true.
          return
        endif
        if ((y(1,j) .ne. y(1,j)) .or. (y(3,j) .ne. y(3,j))) then
          print *, "Error: Beam radius became NaN"
          lfail = .true.
          return
        endif
      enddo
      
      call getderivs(t,dt,y,dy,niz,aion,zion,fviscous,icharge,
     &            lezbeam,lperveance,lemittance,lallez,lezcenter,lviscous,lcurgrid,lfail)
      if (lfail) return

      do j=1,niz
        y( 2,j) = y( 2,j) + halfdt*dy( 2,j)
        y( 4,j) = y( 4,j) + halfdt*dy( 4,j)
        y(10,j) = y(10,j) + halfdt*dy(10,j)
      enddo
     
      return
      end
c=============================================================================
      subroutine loadcharge(niz,y,rpipe,icharge,lfail)
      use HERtmp
      use Constant
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      integer(ISZ):: niz
      real(kind=8):: y(16,niz),rpipe(niz)
      integer(ISZ):: icharge
      logical(ISZ):: lfail

c LOADCHARGE calculates the line-charge density within the slices.
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
          print *, "Error in getselffield:"
          print *, "Beam radius is larger than the aperture."
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
      subroutine getselffield(niz,y,eval,rpipe,icharge,lezbeam,lezcenter,lfail)
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

c GETSELFFIELD calculates the axial space-charge field of a beam
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

c ---   The charge has been loaded already using loadcharge
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
          iz = 0.5 + (y(9,1)-zmmin)/dz
          del2 = (iz+1)*dz-y(9,1)+zmmin
          del1 = 2*dz-del2
          fact1 =  del2-dz/2.
          fact2 =  del1-del2
          fact3 = -del1+dz/2.
          eval(1)   = (fact1*phi(1,iz-1) + fact2*phi(1,iz) + fact3*phi(1,iz+1)) / dz**2
          iz = 0.5 + (y(9,niz)-zmmin)/dz
          del2 = (iz+1)*dz-y(9,niz)+zmmin
          del1 = 2*dz-del2
          fact1 =  del2-dz/2.
          fact2 =  del1-del2
          fact3 = -del1+dz/2.
          eval(niz) = (fact1*phi(1,iz-1) + fact2*phi(1,iz) + fact3*phi(1,iz+1)) / dz**2
          do ii=2,niz-1
            iz = 0.5 + (y(9,ii)-zmmin)/dz
            del2 = (iz+1)*dz-y(9,ii)+zmmin
            del1 = 2*dz-del2
            fact1 =  del2-dz/2.
            fact2 =  del1-del2
            fact3 = -del1+dz/2.
            rm = sqrt(y(1,ii)*y(3,ii))/dr
            irm = rm
            rm = rm - irm
            if (irm > 0) then
              part1 = (irm - 0.5 - 2.*rm*(irm+rm/2.)*((rm**2)/2.+rm*irm-1-2*irm)/(2.*irm+1.))/irm
              part2 = 2.*((rm*(irm+rm/2.))**2)/((irm+1.)*(2.*irm+1.))
            else
              part1 = (1-(rm**2)/2.) - 0.5
              part2 = (rm**2)/2.
            endif
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
      if (icharge == 0) then
        do ii=1,niz
          gtemp(ii) = 2.*log(1.6)
        enddo
      elseif (icharge > 1) then
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

      integer(ISZ):: ii
      real(kind=8):: del1, del2, fact
      real(kind=8):: dzviscous
      real(kind=8):: dbetadz(niz)
      real(kind=8):: temp(niz)

      dzviscous = fviscous * (y(9,niz)-y(9,1)) 
      del1 = y(9,2) - y(9,1)
      del2 = y(9,niz) - y(9,niz-1)
      dbetadz(1)   = (y(10,2)  -y(10,1))    /del1
      dbetadz(niz) = (y(10,niz)-y(10,niz-1))/del2
      do ii=2,niz-1
        del1 = y(9,ii) - y(9,ii-1)
        del2 = y(9,ii+1) - y(9,ii)
        dbetadz(ii) = ((del1/del2)*y(10,ii+1)-(del2/del1)*y(10,ii-1))/(del1+del2)
     &              + (1/del1-1/del2)*y(10,ii)
      enddo
      temp(1)   = 0.
      temp(niz) = 0.
      do ii=2,niz-1
        temp(ii) = (dbetadz(ii)**2)*y(1,ii)*y(3,ii)/y(15,ii)
c       because this is what we actually need
      enddo
      fact = - clight*(dzviscous*y(15,niz/2)/(y(1,niz/2)*y(3,niz/2)))**2
      aviscous(1)   = 0.
      aviscous(niz) = 0.
      do ii=2,niz-1
        del1 = y(9,ii) - y(9,ii-1)
        del2 = y(9,ii+1) - y(9,ii)
        if (dbetadz(ii)*y(10,ii) < 0.) then
          aviscous(ii) = (((del1/del2)*temp(ii+1)
     &                 - (del2/del1)*temp(ii-1))/(del1+del2)
     &                + (1/del1-1/del2)*temp(ii)) * fact * y(1,ii) * y(3,ii) / y(15,ii)
        else
          aviscous(ii) = 0.
        endif
      enddo
      
      return
      end
c=============================================================================
      subroutine savehisther(t,y,niz,nther,lsavehist)
      use HERhist
      integer(ISZ):: niz,nther
      real(kind=8):: t,y(16,niz)
      logical(ISZ):: lsavehist

c     --- Check if saving history
      if (lsavehist) then
c       --- Check if saving history this step. Requires (nhher > 0).
        if (mod(nther,nhher) == 0) then
          if (jhher >= lhher) then
            lhher = lhher + 100
            if (nizhist == 0) nizhist = niz
            call gchange("HERhist",0)
          endif
          call savehisther1(t,y,niz)
        endif
      endif

      return
      end
c=============================================================================
      subroutine savehisther1(t,y,niz)
      use HERhist
      integer(ISZ):: niz
      real(kind=8):: t,y(16,niz)

      integer(ISZ):: ii

      jhher = jhher + 1
      hther(jhher) = t
      
      do ii=1,nizhist
        haher(ii,jhher)   = y( 1,ii)
        hapher(ii,jhher)  = y( 2,ii)
        hbher(ii,jhher)   = y( 3,ii)
        hbpher(ii,jhher)  = y( 4,ii)
        hxher(ii,jhher)   = y( 5,ii)
        hxpher(ii,jhher)  = y( 6,ii)
        hyher(ii,jhher)   = y( 7,ii)
        hypher(ii,jhher)  = y( 8,ii)
        hsher(ii,jhher)   = y( 9,ii)
        hvzher(ii,jhher)  = y(10,ii)
        henxher(ii,jhher) = y(11,ii)
        henyher(ii,jhher) = y(12,ii)
        hcur(ii,jhher)    = y(13,ii)
        hdq(ii,jhher)     = y(14,ii)
        hden(ii,jhher)    = y(15,ii)
      enddo

      return
      end
c=============================================================================
      subroutine readbessel
      use HERbessel
      
      logical(ISZ):: lread = .false.
      save lread

c --- This subroutine reads in the zeros of the Bessel function J0, as well as
c --- besselfactor which is (J1(x))^(-2) in which J1 is the first Bessel
c --- function and x is a zero of the Bessel function J0.   
c --- These values were taken from Mathematica

      if (lread) return

      call gallot("HERbessel",0)
      
      besselzero(1) = 2.40482555769574668
      besselzero(2) = 5.52007811028631323
      besselzero(3) = 8.65372791311373212
      besselzero(4) = 11.7915344390570653
      besselzero(5) = 14.9309177085008819
      besselzero(6) = 18.0710639679159541
      besselzero(7) = 21.2116366298815162
      besselzero(8) = 24.3524715307504324
      besselzero(9) = 27.4934791320408722
      besselzero(10) = 30.6346064684323327
      besselzero(11) = 33.7758202135737883
      besselzero(12) = 36.9170983536641861
      besselzero(13) = 40.0584257646283337
      besselzero(14) = 43.1997917131767916
      besselzero(15) = 46.3411883716618611
      besselzero(16) = 49.4826098973978467
      besselzero(17) = 52.6240518411150226
      besselzero(18) = 55.7655107550199957
      besselzero(19) = 58.906983926080958
      besselzero(20) = 62.0484691902271823
      besselzero(21) = 65.1899648002068676
      besselzero(22) = 68.3314693298567998
      besselzero(23) = 71.4729816035937393
      besselzero(24) = 74.6145006437018487
      besselzero(25) = 77.7560256303880503
      besselzero(26) = 80.897555871137623
      besselzero(27) = 84.0390907769381989
      besselzero(28) = 87.180629843641153
      besselzero(29) = 90.3221726372104782
      besselzero(30) = 93.4637187819447845
      besselzero(31) = 96.6052679509962608
      besselzero(32) = 99.7468198586805954
      besselzero(33) = 102.888374254194792
      besselzero(34) = 106.029930917647763
      besselzero(35) = 109.171489650901198
      besselzero(36) = 112.313050281501314
      besselzero(37) = 115.454612654593402
      besselzero(38) = 118.596176631727301
      besselzero(39) = 121.73774208874124
      besselzero(40) = 124.879308913965081
      besselzero(41) = 128.020877006687872
      besselzero(42) = 131.162446275845789
      besselzero(43) = 134.304016638894041
      besselzero(44) = 137.445588020833375
      besselzero(45) = 140.587160353367401
      besselzero(46) = 143.728733574169931
      besselzero(47) = 146.870307626246692
      besselzero(48) = 150.01188245737711
      besselzero(49) = 153.153458019624792
      besselzero(50) = 156.295034268906967
      besselzero(51) = 159.436611164614934
      besselzero(52) = 162.578188669278445
      besselzero(53) = 165.719766748268293
      besselzero(54) = 168.861345369531923
      besselzero(55) = 172.002924503358372
      besselzero(56) = 175.144504122168109
      besselzero(57) = 178.286084200325367
      besselzero(58) = 181.427664713969782
      besselzero(59) = 184.569245640865489
      besselzero(60) = 187.710826960264931
      besselzero(61) = 190.8524086527866
      besselzero(62) = 193.993990700304408
      besselzero(63) = 197.135573085847504
      besselzero(64) = 200.277155793509908
      besselzero(65) = 203.418738808368049
      besselzero(66) = 206.560322116406247
      besselzero(67) = 209.701905704448687
      besselzero(68) = 212.843489560097332
      besselzero(69) = 215.985073671675475
      besselzero(70) = 219.126658028176057
      besselzero(71) = 222.26824261921414
      besselzero(72) = 225.409827434983833
      besselzero(73) = 228.55141246621824
      besselzero(74) = 231.692997704153125
      besselzero(75) = 234.834583140493347
      besselzero(76) = 237.976168767381413
      besselzero(77) = 241.117754577369769
      besselzero(78) = 244.259340563393489
      besselzero(79) = 247.400926718747005
      besselzero(80) = 250.54251303706061
      besselzero(81) = 253.684099512280392
      besselzero(82) = 256.825686138648645
      besselzero(83) = 259.967272910685665
      besselzero(84) = 263.108859823173757
      besselzero(85) = 266.250446871141432
      besselzero(86) = 269.392034049849016
      besselzero(87) = 272.533621354775368
      besselzero(88) = 275.67520878160554
      besselzero(89) = 278.816796326218874
      besselzero(90) = 281.958383984678562
      besselzero(91) = 285.09997175322109
      besselzero(92) = 288.241559628247259
      besselzero(93) = 291.383147606312888
      besselzero(94) = 294.52473568412083
      besselzero(95) = 297.666323858513024
      besselzero(96) = 300.807912126463516
      besselzero(97) = 303.949500485071366
      besselzero(98) = 307.091088931554256
      besselzero(99) = 310.232677463242723
      besselzero(100) = 313.374266077574192
      besselzero(101) = 316.515854772087923
      besselzero(102) = 319.65744354441985
      besselzero(103) = 322.799032392297968
      besselzero(104) = 325.940621313537892
      besselzero(105) = 329.082210306038547
      besselzero(106) = 332.223799367778571
      besselzero(107) = 335.365388496812011
      besselzero(108) = 338.506977691265253
      besselzero(109) = 341.64856694933392
      besselzero(110) = 344.790156269278824
      besselzero(111) = 347.931745649424106
      besselzero(112) = 351.073335088153593
      besselzero(113) = 354.21492458390853
      besselzero(114) = 357.356514135185054
      besselzero(115) = 360.498103740531528
      besselzero(116) = 363.639693398546759
      besselzero(117) = 366.781283107877253
      besselzero(118) = 369.922872867215701
      besselzero(119) = 373.064462675298713
      besselzero(120) = 376.206052530905221
      besselzero(121) = 379.347642432854659
      besselzero(122) = 382.489232380004793
      besselzero(123) = 385.63082237125128
      besselzero(124) = 388.772412405524914
      besselzero(125) = 391.91400248179109
      besselzero(126) = 395.055592599047997
      besselzero(127) = 398.19718275632546
      besselzero(128) = 401.338772952683608
      besselzero(129) = 404.480363187211988
      besselzero(130) = 407.621953459027874
      besselzero(131) = 410.763543767275951
      besselzero(132) = 413.905134111126482
      besselzero(133) = 417.046724489774867
      besselzero(134) = 420.188314902440929
      besselzero(135) = 423.329905348366875
      besselzero(136) = 426.471495826818003
      besselzero(137) = 429.613086337080751
      besselzero(138) = 432.754676878462163
      besselzero(139) = 435.89626745028962
      besselzero(140) = 439.037858051909335
      besselzero(141) = 442.179448682686437
      besselzero(142) = 445.32103934200391
      besselzero(143) = 448.462630029261877
      besselzero(144) = 451.604220743877071
      besselzero(145) = 454.745811485282835
      besselzero(146) = 457.887402252927611
      besselzero(147) = 461.028993046275026
      besselzero(148) = 464.170583864803099
      besselzero(149) = 467.312174708004057
      besselzero(150) = 470.453765575383542
      besselzero(151) = 473.595356466460515
      besselzero(152) = 476.736947380766551
      besselzero(153) = 479.878538317845304
      besselzero(154) = 483.020129277252419
      besselzero(155) = 486.161720258555085
      besselzero(156) = 489.303311261331597
      besselzero(157) = 492.444902285170904
      besselzero(158) = 495.586493329672528
      besselzero(159) = 498.728084394446025
      besselzero(160) = 501.869675479110632
      besselzero(161) = 505.011266583295181
      besselzero(162) = 508.152857706637561
      besselzero(163) = 511.294448848784721
      besselzero(164) = 514.436040009392048
      besselzero(165) = 517.577631188123454
      besselzero(166) = 520.719222384651025
      besselzero(167) = 523.860813598654395
      besselzero(168) = 527.00240482982128
      besselzero(169) = 530.143996077846146
      besselzero(170) = 533.285587342431544
      besselzero(171) = 536.427178623286149
      besselzero(172) = 539.568769920125923
      besselzero(173) = 542.710361232673399
      besselzero(174) = 545.851952560657061
      besselzero(175) = 548.993543903812142
      besselzero(176) = 552.13513526187965
      besselzero(177) = 555.276726634606454
      besselzero(178) = 558.418318021745108
      besselzero(179) = 561.559909423053671
      besselzero(180) = 564.7015008382958
      besselzero(181) = 567.843092267240212
      besselzero(182) = 570.984683709660689
      besselzero(183) = 574.126275165335986
      besselzero(184) = 577.267866634049831
      besselzero(185) = 580.409458115590304
      besselzero(186) = 583.551049609750283
      besselzero(187) = 586.692641116327262
      besselzero(188) = 589.834232635122557
      besselzero(189) = 592.975824165942189
      besselzero(190) = 596.117415708595999
      besselzero(191) = 599.259007262897736
      besselzero(192) = 602.400598828665323
      besselzero(193) = 605.542190405720237
      besselzero(194) = 608.683781993887684
      besselzero(195) = 611.825373592996424
      besselzero(196) = 614.966965202878945
      besselzero(197) = 618.108556823370935
      besselzero(198) = 621.250148454311279
      besselzero(199) = 624.391740095542591
      besselzero(200) = 627.533331746909972
      besselzero(201) = 630.674923408262078
      besselzero(202) = 633.816515079450493
      besselzero(203) = 636.958106760329556
      besselzero(204) = 640.099698450756715
      besselzero(205) = 643.241290150591993
      besselzero(206) = 646.382881859698255
      besselzero(207) = 649.52447357794103
      besselzero(208) = 652.666065305188158
      besselzero(209) = 655.80765704131041
      besselzero(210) = 658.949248786180863
      besselzero(211) = 662.090840539674907
      besselzero(212) = 665.232432301670417
      besselzero(213) = 668.374024072047667
      besselzero(214) = 671.515615850688618
      besselzero(215) = 674.657207637478251
      besselzero(216) = 677.798799432302967
      besselzero(217) = 680.940391235051656
      besselzero(218) = 684.081983045615249
      besselzero(219) = 687.223574863886454
      besselzero(220) = 690.365166689760023
      besselzero(221) = 693.506758523132749
      besselzero(222) = 696.648350363903112
      besselzero(223) = 699.789942211971371
      besselzero(224) = 702.931534067239915
      besselzero(225) = 706.073125929612466
      besselzero(226) = 709.214717798994609
      besselzero(227) = 712.356309675293886
      besselzero(228) = 715.497901558418725
      besselzero(229) = 718.639493448280042
      besselzero(230) = 721.781085344789552
      besselzero(231) = 724.922677247861013
      besselzero(232) = 728.064269157409249
      besselzero(233) = 731.205861073351126
      besselzero(234) = 734.34745299560431
      besselzero(235) = 737.489044924088421
      besselzero(236) = 740.630636858723789
      besselzero(237) = 743.772228799432788
      besselzero(238) = 746.91382074613859
      besselzero(239) = 750.055412698765877
      besselzero(240) = 753.197004657240665
      besselzero(241) = 756.338596621489944
      besselzero(242) = 759.480188591442129
      besselzero(243) = 762.621780567026785
      besselzero(244) = 765.763372548174459
      besselzero(245) = 768.904964534817203
      besselzero(246) = 772.046556526887695
      besselzero(247) = 775.18814852432003
      besselzero(248) = 778.329740527049285
      besselzero(249) = 781.471332535011509
      besselzero(250) = 784.612924548143997
      besselzero(251) = 787.754516566384666
      besselzero(252) = 790.896108589672941
      besselzero(253) = 794.037700617948694
      besselzero(254) = 797.179292651153126
      besselzero(255) = 800.320884689227973
      besselzero(256) = 803.462476732116215
      besselzero(257) = 806.60406877976154
      besselzero(258) = 809.745660832108704
      besselzero(259) = 812.887252889103351
      besselzero(260) = 816.028844950691123
      besselzero(261) = 819.170437016819974
      besselzero(262) = 822.312029087437679
      besselzero(263) = 825.453621162492723
      besselzero(264) = 828.59521324193448
      besselzero(265) = 831.736805325713746
      besselzero(266) = 834.878397413781314
      besselzero(267) = 838.019989506088869
      besselzero(268) = 841.161581602588981
      besselzero(269) = 844.303173703234932
      besselzero(270) = 847.444765807980715
      besselzero(271) = 850.586357916780499
      besselzero(272) = 853.727950029589877
      besselzero(273) = 856.869542146364615
      besselzero(274) = 860.011134267061194
      besselzero(275) = 863.152726391636981
      besselzero(276) = 866.294318520049522
      besselzero(277) = 869.435910652257604
      besselzero(278) = 872.577502788219838
      besselzero(279) = 875.71909492789608
      besselzero(280) = 878.860687071246182
      besselzero(281) = 882.002279218231421
      besselzero(282) = 885.14387136881254
      besselzero(283) = 888.285463522951701
      besselzero(284) = 891.427055680611246
      besselzero(285) = 894.568647841754049
      besselzero(286) = 897.710240006343518
      besselzero(287) = 900.851832174343769
      besselzero(288) = 903.993424345718743
      besselzero(289) = 907.135016520433978
      besselzero(290) = 910.276608698454658
      besselzero(291) = 913.4182008797465
      besselzero(292) = 916.559793064276107
      besselzero(293) = 919.701385252010439
      besselzero(294) = 922.842977442916279
      besselzero(295) = 925.984569636961651
      besselzero(296) = 929.126161834114761
      besselzero(297) = 932.267754034343987
      besselzero(298) = 935.409346237618422
      besselzero(299) = 938.550938443907689
      besselzero(300) = 941.692530653181236
      besselzero(301) = 944.834122865409575
      besselzero(302) = 947.975715080563396
      besselzero(303) = 951.11730729861339
      besselzero(304) = 954.258899519531134
      besselzero(305) = 957.400491743288384
      besselzero(306) = 960.542083969857252
      besselzero(307) = 963.683676199210381
      besselzero(308) = 966.825268431320594
      besselzero(309) = 969.96686066616089
      besselzero(310) = 973.108452903705156
      besselzero(311) = 976.250045143926925
      besselzero(312) = 979.391637386800795
      besselzero(313) = 982.53322963230083
      besselzero(314) = 985.674821880402518
      besselzero(315) = 988.816414131080634
      besselzero(316) = 991.958006384311019
      besselzero(317) = 995.099598640068983
      besselzero(318) = 998.241190898331254
      besselzero(319) = 1001.38278315907403
      besselzero(320) = 1004.5243754222739
      besselzero(321) = 1007.66596768790783
      besselzero(322) = 1010.80755995595361
      besselzero(323) = 1013.94915222638859
      besselzero(324) = 1017.09074449919012
      besselzero(325) = 1020.23233677433688
      besselzero(326) = 1023.37392905180735
      besselzero(327) = 1026.51552133157997
      besselzero(328) = 1029.65711361363365
      besselzero(329) = 1032.79870589794775
      besselzero(330) = 1035.94029818450117
      besselzero(331) = 1039.08189047327437
      besselzero(332) = 1042.22348276424692
      besselzero(333) = 1045.36507505739884
      besselzero(334) = 1048.50666735271058
      besselzero(335) = 1051.64825965016306
      besselzero(336) = 1054.78985194973695
      besselzero(337) = 1057.93144425141316
      besselzero(338) = 1061.07303655517348
      besselzero(339) = 1064.21462886099926
      besselzero(340) = 1067.35622116887229
      besselzero(341) = 1070.49781347877415
      besselzero(342) = 1073.6394057906875
      besselzero(343) = 1076.78099810459459
      besselzero(344) = 1079.92259042047811
      besselzero(345) = 1083.06418273832072
      besselzero(346) = 1086.20577505810533
      besselzero(347) = 1089.34736737981507
      besselzero(348) = 1092.48895970343395
      besselzero(349) = 1095.6305520289451
      besselzero(350) = 1098.77214435633231
      besselzero(351) = 1101.91373668557934
      besselzero(352) = 1105.05532901667069
      besselzero(353) = 1108.19692134959035
      besselzero(354) = 1111.33851368432323
      besselzero(355) = 1114.48010602085356
      besselzero(356) = 1117.62169835916624
      besselzero(357) = 1120.76329069924685
      besselzero(358) = 1123.90488304108004
      besselzero(359) = 1127.04647538465119
      besselzero(360) = 1130.18806772994562
      besselzero(361) = 1133.32966007694913
      besselzero(362) = 1136.47125242564817
      besselzero(363) = 1139.61284477602786
      besselzero(364) = 1142.75443712807489
      besselzero(365) = 1145.89602948177526
      besselzero(366) = 1149.03762183711544
      besselzero(367) = 1152.17921419408209
      besselzero(368) = 1155.32080655266189
      besselzero(369) = 1158.46239891284175
      besselzero(370) = 1161.60399127460856
      besselzero(371) = 1164.74558363794966
      besselzero(372) = 1167.88717600285241
      besselzero(373) = 1171.02876836930414
      besselzero(374) = 1174.17036073729219
      besselzero(375) = 1177.31195310680458
      besselzero(376) = 1180.4535454778291
      besselzero(377) = 1183.5951378503533
      besselzero(378) = 1186.73673022436609
      besselzero(379) = 1189.87832259985504
      besselzero(380) = 1193.01991497680837
      besselzero(381) = 1196.16150735521497
      besselzero(382) = 1199.30309973506332
      besselzero(383) = 1202.44469211634208
      besselzero(384) = 1205.58628449903992
      besselzero(385) = 1208.7278768831462
      besselzero(386) = 1211.86946926864935
      besselzero(387) = 1215.01106165553895
      besselzero(388) = 1218.15265404380435
      besselzero(389) = 1221.29424643343442
      besselzero(390) = 1224.43583882441897
      besselzero(391) = 1227.577431216748
      besselzero(392) = 1230.71902361041063
      besselzero(393) = 1233.86061600539687
      besselzero(394) = 1237.00220840169694
      besselzero(395) = 1240.14380079930019
      besselzero(396) = 1243.2853931981973
      besselzero(397) = 1246.42698559837827
      besselzero(398) = 1249.56857799983334
      besselzero(399) = 1252.71017040255316
      besselzero(400) = 1255.85176280652799
      besselzero(401) = 1258.99335521174893
      besselzero(402) = 1262.13494761820577
      besselzero(403) = 1265.27654002589029
      besselzero(404) = 1268.41813243479273
      besselzero(405) = 1271.5597248449042
      besselzero(406) = 1274.70131725621582
      besselzero(407) = 1277.84290966871871
      besselzero(408) = 1280.98450208240421
      besselzero(409) = 1284.12609449726344
      besselzero(410) = 1287.26768691328796
      besselzero(411) = 1290.40927933046933
      besselzero(412) = 1293.55087174879868
      besselzero(413) = 1296.69246416826844
      besselzero(414) = 1299.83405658886974
      besselzero(415) = 1302.97564901059437
      besselzero(416) = 1306.11724143343455
      besselzero(417) = 1309.25883385738228
      besselzero(418) = 1312.40042628242914
      besselzero(419) = 1315.54201870856757
      besselzero(420) = 1318.68361113579002
      besselzero(421) = 1321.82520356408828
      besselzero(422) = 1324.96679599345457
      besselzero(423) = 1328.10838842388223
      besselzero(424) = 1331.24998085536305
      besselzero(425) = 1334.39157328788948
      besselzero(426) = 1337.53316572145442
      besselzero(427) = 1340.67475815605052
      besselzero(428) = 1343.81635059167092
      besselzero(429) = 1346.95794302830762
      besselzero(430) = 1350.09953546595418
      besselzero(431) = 1353.24112790460349
      besselzero(432) = 1356.382720344248
      besselzero(433) = 1359.52431278488172
      besselzero(434) = 1362.66590522649733
      besselzero(435) = 1365.80749766908793
      besselzero(436) = 1368.94909011264709
      besselzero(437) = 1372.09068255716792
      besselzero(438) = 1375.232275002644
      besselzero(439) = 1378.37386744906886
      besselzero(440) = 1381.51545989643587
      besselzero(441) = 1384.65705234473857
      besselzero(442) = 1387.79864479397074
      besselzero(443) = 1390.94023724412597
      besselzero(444) = 1394.08182969519822
      besselzero(445) = 1397.22342214718109
      besselzero(446) = 1400.36501460006856
      besselzero(447) = 1403.50660705385443
      besselzero(448) = 1406.6481995085331
      besselzero(449) = 1409.78979196409799
      besselzero(450) = 1412.9313844205436
      besselzero(451) = 1416.0729768778637
      besselzero(452) = 1419.21456933605298
      besselzero(453) = 1422.35616179510548
      besselzero(454) = 1425.49775425501534
      besselzero(455) = 1428.63934671577696
      besselzero(456) = 1431.78093917738502
      besselzero(457) = 1434.92253163983329
      besselzero(458) = 1438.06412410311672
      besselzero(459) = 1441.20571656723033
      besselzero(460) = 1444.3473090321679
      besselzero(461) = 1447.48890149792437
      besselzero(462) = 1450.63049396449432
      besselzero(463) = 1453.77208643187261
      besselzero(464) = 1456.91367890005407
      besselzero(465) = 1460.05527136903339
      besselzero(466) = 1463.19686383880531
      besselzero(467) = 1466.33845630936487
      besselzero(468) = 1469.48004878070719
      besselzero(469) = 1472.62164125282701
      besselzero(470) = 1475.76323372571955
      besselzero(471) = 1478.90482619937984
      besselzero(472) = 1482.04641867380315
      besselzero(473) = 1485.18801114898391
      besselzero(474) = 1488.3296036249181
      besselzero(475) = 1491.4711961016005
      besselzero(476) = 1494.61278857902657
      besselzero(477) = 1497.75438105719205
      besselzero(478) = 1500.89597353609178
      besselzero(479) = 1504.03756601572098
      besselzero(480) = 1507.17915849607582
      besselzero(481) = 1510.32075097715123
      besselzero(482) = 1513.46234345894288
      besselzero(483) = 1516.60393594144604
      besselzero(484) = 1519.74552842465682
      besselzero(485) = 1522.88712090857014
      besselzero(486) = 1526.0287133931822
      besselzero(487) = 1529.17030587848863
      besselzero(488) = 1532.31189836448544
      besselzero(489) = 1535.45349085116766
      besselzero(490) = 1538.59508333853174
      besselzero(491) = 1541.73667582657287
      besselzero(492) = 1544.87826831528751
      besselzero(493) = 1548.01986080467157
      besselzero(494) = 1551.16145329472061
      besselzero(495) = 1554.30304578543054
      besselzero(496) = 1557.44463827679791
      besselzero(497) = 1560.58623076881808
      besselzero(498) = 1563.72782326148795
      besselzero(499) = 1566.86941575480265
      besselzero(500) = 1570.01100824875887
      besselzero(501) = 1573.15260074335299
      besselzero(502) = 1576.29419323858037
      besselzero(503) = 1579.4357857344381
      besselzero(504) = 1582.57737823092181
      besselzero(505) = 1585.71897072802814
      besselzero(506) = 1588.86056322575335
      besselzero(507) = 1592.0021557240938
      besselzero(508) = 1595.14374822304549
      besselzero(509) = 1598.28534072260541
      besselzero(510) = 1601.42693322276912
      besselzero(511) = 1604.56852572353403
      besselzero(512) = 1607.71011822489606
      besselzero(513) = 1610.85171072685229
      besselzero(514) = 1613.99330322939853
      besselzero(515) = 1617.13489573253133
      besselzero(516) = 1620.27648823624792
      besselzero(517) = 1623.4180807405445
      besselzero(518) = 1626.55967324541776
      besselzero(519) = 1629.70126575086453
      besselzero(520) = 1632.84285825688124
      besselzero(521) = 1635.98445076346443
      besselzero(522) = 1639.12604327061162
      besselzero(523) = 1642.26763577831853
      besselzero(524) = 1645.40922828658331
      besselzero(525) = 1648.55082079540116
      besselzero(526) = 1651.69241330477012
      besselzero(527) = 1654.83400581468664
      besselzero(528) = 1657.97559832514762
      besselzero(529) = 1661.11719083615021
      besselzero(530) = 1664.2587833476905
      besselzero(531) = 1667.40037585976663
      besselzero(532) = 1670.54196837237505
      besselzero(533) = 1673.68356088551228
      besselzero(534) = 1676.82515339917621
      besselzero(535) = 1679.96674591336336
      besselzero(536) = 1683.10833842807099
      besselzero(537) = 1686.24993094329544
      besselzero(538) = 1689.39152345903558
      besselzero(539) = 1692.53311597528677
      besselzero(540) = 1695.67470849204688
      besselzero(541) = 1698.81630100931317
      besselzero(542) = 1701.95789352708288
      besselzero(543) = 1705.09948604535281
      besselzero(544) = 1708.24107856412084
      besselzero(545) = 1711.38267108338348
      besselzero(546) = 1714.52426360313836
      besselzero(547) = 1717.66585612338325
      besselzero(548) = 1720.8074486441145
      besselzero(549) = 1723.94904116533043
      besselzero(550) = 1727.09063368702757
      besselzero(551) = 1730.23222620920407
      besselzero(552) = 1733.37381873185663
      besselzero(553) = 1736.51541125498276
      besselzero(554) = 1739.65700377858034
      besselzero(555) = 1742.79859630264644
      besselzero(556) = 1745.94018882717927
      besselzero(557) = 1749.08178135217511
      besselzero(558) = 1752.22337387763182
      besselzero(559) = 1755.36496640354755
      besselzero(560) = 1758.5065589299198
      besselzero(561) = 1761.64815145674547
      besselzero(562) = 1764.78974398402215
      besselzero(563) = 1767.93133651174798
      besselzero(564) = 1771.0729290399203
      besselzero(565) = 1774.21452156853725
      besselzero(566) = 1777.35611409759535
      besselzero(567) = 1780.49770662709283
      besselzero(568) = 1783.63929915702819
      besselzero(569) = 1786.78089168739751
      besselzero(570) = 1789.92248421819999
      besselzero(571) = 1793.0640767494328
      besselzero(572) = 1796.20566928109326
      besselzero(573) = 1799.34726181317961
      besselzero(574) = 1802.48885434568961
      besselzero(575) = 1805.63044687862053
      besselzero(576) = 1808.77203941197102
      besselzero(577) = 1811.9136319457386
      besselzero(578) = 1815.05522447992078
      besselzero(579) = 1818.19681701451507
      besselzero(580) = 1821.33840954952077
      besselzero(581) = 1824.48000208493397
      besselzero(582) = 1827.62159462075395
      besselzero(583) = 1830.76318715697823
      besselzero(584) = 1833.90477969360432
      besselzero(585) = 1837.04637223063045
      besselzero(586) = 1840.18796476805448
      besselzero(587) = 1843.32955730587499
      besselzero(588) = 1846.47114984408915
      besselzero(589) = 1849.61274238269517
      besselzero(590) = 1852.75433492169128
      besselzero(591) = 1855.89592746107499
      besselzero(592) = 1859.03752000084487
      besselzero(593) = 1862.17911254099917
      besselzero(594) = 1865.32070508153502
      besselzero(595) = 1868.46229762245137
      besselzero(596) = 1871.60389016374608
      besselzero(597) = 1874.74548270541667
      besselzero(598) = 1877.88707524746243
      besselzero(599) = 1881.02866778988015
      besselzero(600) = 1884.17026033266878
      besselzero(601) = 1887.31185287582619
      besselzero(602) = 1890.45344541935094
      besselzero(603) = 1893.59503796324091
      besselzero(604) = 1896.73663050749397
      besselzero(605) = 1899.87822305210904
      besselzero(606) = 1903.0198155970833
      besselzero(607) = 1906.16140814241603
      besselzero(608) = 1909.30300068810474
      besselzero(609) = 1912.44459323414801
      besselzero(610) = 1915.58618578054442
      besselzero(611) = 1918.72777832729148
      besselzero(612) = 1921.86937087438813
      besselzero(613) = 1925.01096342183188
      besselzero(614) = 1928.15255596962202
      besselzero(615) = 1931.29414851775607
      besselzero(616) = 1934.43574106623294
      besselzero(617) = 1937.57733361505017
      besselzero(618) = 1940.71892616420704
      besselzero(619) = 1943.86051871370142
      besselzero(620) = 1947.00211126353153
      besselzero(621) = 1950.1437038136963
      besselzero(622) = 1953.28529636419325
      besselzero(623) = 1956.42688891502167
      besselzero(624) = 1959.56848146617943
      besselzero(625) = 1962.7100740176651
      besselzero(626) = 1965.85166656947691
      besselzero(627) = 1968.99325912161415
      besselzero(628) = 1972.13485167407398
      besselzero(629) = 1975.27644422685568
      besselzero(630) = 1978.41803677995748
      besselzero(631) = 1981.55962933337832
      besselzero(632) = 1984.70122188711606
      besselzero(633) = 1987.84281444116928
      besselzero(634) = 1990.98440699553655
      besselzero(635) = 1994.12599955021683
      besselzero(636) = 1997.26759210520796
      besselzero(637) = 2000.40918466050854
      besselzero(638) = 2003.5507772161174
      besselzero(639) = 2006.69236977203357
      besselzero(640) = 2009.83396232825484
      besselzero(641) = 2012.97555488477986
      besselzero(642) = 2016.11714744160774
      besselzero(643) = 2019.25873999873628
      besselzero(644) = 2022.40033255616458
      besselzero(645) = 2025.5419251138913
      besselzero(646) = 2028.68351767191512
      besselzero(647) = 2031.82511023023426
      besselzero(648) = 2034.96670278884783
      besselzero(649) = 2038.10829534775361
      besselzero(650) = 2041.24988790695116
      besselzero(651) = 2044.3914804664387
      besselzero(652) = 2047.53307302104481
      besselzero(653) = 2050.67466559137612
      besselzero(654) = 2053.81625814162571
      besselzero(655) = 2056.95785071219594
      besselzero(656) = 2060.09944326334082
      besselzero(657) = 2063.24103583415485
      besselzero(658) = 2066.38262838618036
      besselzero(659) = 2069.52422095724308
      besselzero(660) = 2072.6658135101319
      besselzero(661) = 2075.80740608144953
      besselzero(662) = 2078.94899863518789
      besselzero(663) = 2082.09059120676398
      besselzero(664) = 2085.23218376133723
      besselzero(665) = 2088.37377633317711
      besselzero(666) = 2091.51536888857014
      besselzero(667) = 2094.65696146067825
      besselzero(668) = 2097.7985540168782
      besselzero(669) = 2100.94014658925809
      besselzero(670) = 2104.08173914625029
      besselzero(671) = 2107.22333171890552
      besselzero(672) = 2110.36492427667798
      besselzero(673) = 2113.50651684961299
      besselzero(674) = 2116.64810940815373
      besselzero(675) = 2119.78970198137073
      besselzero(676) = 2122.93129454065709
      besselzero(677) = 2126.07288711414366
      besselzero(678) = 2129.21447967418142
      besselzero(679) = 2132.35607224805034
      besselzero(680) = 2135.49766480893854
      besselzero(681) = 2138.63925738282256
      besselzero(682) = 2141.78084994312856
      besselzero(683) = 2144.92244251559061
      besselzero(684) = 2148.06403507816234
      besselzero(685) = 2151.20562765825917
      besselzero(686) = 2154.34722021581137
      besselzero(687) = 2157.48881279603121
      besselzero(688) = 2160.63040535444139
      besselzero(689) = 2163.77199793479269
      besselzero(690) = 2166.91359049404397
      besselzero(691) = 2170.05518307453337
      besselzero(692) = 2173.19677563461111
      besselzero(693) = 2176.33836821524573
      besselzero(694) = 2179.47996077613392
      besselzero(695) = 2182.62155335692131
      besselzero(696) = 2185.76314591860576
      besselzero(697) = 2188.90473849955125
      besselzero(698) = 2192.04633106201685
      besselzero(699) = 2195.18792364312664
      besselzero(700) = 2198.32951620636007
      besselzero(701) = 2201.47110878764085
      besselzero(702) = 2204.61270135162745
      besselzero(703) = 2207.75429393308453
      besselzero(704) = 2210.89588649781143
      besselzero(705) = 2214.03747907945014
      besselzero(706) = 2217.17907164490402
      besselzero(707) = 2220.32066422672924
      besselzero(708) = 2223.46225679289766
      besselzero(709) = 2226.60384937491429
      besselzero(710) = 2229.74544194178481
      besselzero(711) = 2232.88703452399861
      besselzero(712) = 2236.02862709156102
      besselzero(713) = 2239.17021967397378
      besselzero(714) = 2242.31181224220588
      besselzero(715) = 2245.45340482480693
      besselzero(716) = 2248.59499739371626
      besselzero(717) = 2251.73658997662062
      besselzero(718) = 2254.87818254630711
      besselzero(719) = 2258.01977512914309
      besselzero(720) = 2261.16136769815545
      besselzero(721) = 2264.30296027946464
      besselzero(722) = 2267.4445528506939
      besselzero(723) = 2270.58614543967074
      besselzero(724) = 2273.7277380057228
      besselzero(725) = 2276.86933059474405
      besselzero(726) = 2280.01092316158526
      besselzero(727) = 2283.15251575065936
      besselzero(728) = 2286.29410831827417
      besselzero(729) = 2289.43570090741044
      besselzero(730) = 2292.57729347578376
      besselzero(731) = 2295.71888606498841
      besselzero(732) = 2298.86047863410647
      besselzero(733) = 2302.00207122338884
      besselzero(734) = 2305.14366379323654
      besselzero(735) = 2308.28525638260284
      besselzero(736) = 2311.42684895316685
      besselzero(737) = 2314.56844154262464
      besselzero(738) = 2317.71003411389164
      besselzero(739) = 2320.85162670344757
      besselzero(740) = 2323.9932192754038
      besselzero(741) = 2327.1348118650641
      besselzero(742) = 2330.27640443769801
      besselzero(743) = 2333.41799702746888
      besselzero(744) = 2336.55958960076759
      besselzero(745) = 2339.70118219065525
      besselzero(746) = 2342.84277476460633
      besselzero(747) = 2345.98436735461657
      besselzero(748) = 2349.12595992920802
      besselzero(749) = 2352.26755251934438
      besselzero(750) = 2355.40914509456555
      besselzero(751) = 2358.55073768484402
      besselzero(752) = 2361.69233026070068
      besselzero(753) = 2364.83392285110616
      besselzero(754) = 2367.97551542747308
      besselzero(755) = 2371.11710801789632
      besselzero(756) = 2374.25870059514831
      besselzero(757) = 2377.40029318707968
      besselzero(758) = 2380.54188576671421
      besselzero(759) = 2383.68347835716676
      besselzero(760) = 2386.82507092970919
      besselzero(761) = 2389.96666352628661
      besselzero(762) = 2393.10825609957955
      besselzero(763) = 2396.24984869613744
      besselzero(764) = 2399.39144127016357
      besselzero(765) = 2402.53303386671035
      besselzero(766) = 2405.6746264414559
      besselzero(767) = 2408.81621903800091
      besselzero(768) = 2411.95781161345168
      besselzero(769) = 2415.09940421000424
      besselzero(770) = 2418.24099678614512
      besselzero(771) = 2421.38258938271277
      besselzero(772) = 2424.52418195953134
      besselzero(773) = 2427.66577455612253
      besselzero(774) = 2430.80736713360412
      besselzero(775) = 2433.94895973022684
      besselzero(776) = 2437.09055230835991
      besselzero(777) = 2440.23214490502082
      besselzero(778) = 2443.37373748379205
      besselzero(779) = 2446.51533008049915
      besselzero(780) = 2449.65692265989565
      besselzero(781) = 2452.79851525665559
      besselzero(782) = 2455.94010783666671
      besselzero(783) = 2459.08170043348572
      besselzero(784) = 2462.22329301409947
      besselzero(785) = 2465.36488561098465
      besselzero(786) = 2468.50647819219082
      besselzero(787) = 2471.64807078914749
      besselzero(788) = 2474.78966337092476
      besselzero(789) = 2477.93125596794094
      besselzero(790) = 2481.07284855029863
      besselzero(791) = 2484.21444114749956
      besselzero(792) = 2487.3560337305447
      besselzero(793) = 2490.49762632753513
      besselzero(794) = 2493.63921890972406
      besselzero(795) = 2496.78081150495057
      besselzero(796) = 2499.92240408936527
      besselzero(797) = 2503.06399669262891
      besselzero(798) = 2506.20558927138947
      besselzero(799) = 2509.34718187456784
      besselzero(800) = 2512.4887744540354
      besselzero(801) = 2515.63036705713738
      besselzero(802) = 2518.77195963729638
      besselzero(803) = 2521.91355224033308
      besselzero(804) = 2525.05514482116888
      besselzero(805) = 2528.19673742415051
      besselzero(806) = 2531.33833000564845
      besselzero(807) = 2534.47992260858434
      besselzero(808) = 2537.62151519073153
      besselzero(809) = 2540.76310779362968
      besselzero(810) = 2543.90470037641236
      besselzero(811) = 2547.04629297928253
      besselzero(812) = 2550.18788556268738
      besselzero(813) = 2553.32947816553801
      besselzero(814) = 2556.47107074955304
      besselzero(815) = 2559.61266335239168
      besselzero(816) = 2562.75425593700356
      besselzero(817) = 2565.89584853983776
      besselzero(818) = 2569.03744112503584
      besselzero(819) = 2572.17903372787404
      besselzero(820) = 2575.32062631364544
      besselzero(821) = 2578.46221891649429
      besselzero(822) = 2581.6038115028288
      besselzero(823) = 2584.74540410569542
      besselzero(824) = 2587.88699669258326
      besselzero(825) = 2591.02858929547341
      besselzero(826) = 2594.17018188289238
      besselzero(827) = 2597.3117744857963
      besselzero(828) = 2600.45336707375662
      besselzero(829) = 2603.59495967679865
      besselzero(830) = 2606.73655226541223
      besselzero(831) = 2609.87814486819269
      besselzero(832) = 2613.01973745588611
      besselzero(833) = 2616.16133005682449
      besselzero(834) = 2619.30292264673525
      besselzero(835) = 2622.44451525580668
      besselzero(836) = 2625.58610783991719
      besselzero(837) = 2628.7277004488545
      besselzero(838) = 2631.86929303363692
      besselzero(839) = 2635.01088564245167
      besselzero(840) = 2638.15247822789222
      besselzero(841) = 2641.29407083659462
      besselzero(842) = 2644.43566342267866
      besselzero(843) = 2647.57725603128069
      besselzero(844) = 2650.71884861799267
      besselzero(845) = 2653.86044122650321
      besselzero(846) = 2657.00203381383071
      besselzero(847) = 2660.14362642225998
      besselzero(848) = 2663.28521901018921
      besselzero(849) = 2666.42681161854611
      besselzero(850) = 2669.56840420706509
      besselzero(851) = 2672.70999681535847
      besselzero(852) = 2675.85158940445477
      besselzero(853) = 2678.9931820126922
      besselzero(854) = 2682.13477460235383
      besselzero(855) = 2685.27636721054552
      besselzero(856) = 2688.41795980076003
      besselzero(857) = 2691.55955240891264
      besselzero(858) = 2694.70114499966939
      besselzero(859) = 2697.84273760779047
      besselzero(860) = 2700.9843301990779
      besselzero(861) = 2704.12592280717589
      besselzero(862) = 2707.26751539898602
      besselzero(863) = 2710.40910800706669
      besselzero(864) = 2713.55070059937641
      besselzero(865) = 2716.6922932074291
      besselzero(866) = 2719.83388580024953
      besselzero(867) = 2722.97547840840348
      besselzero(868) = 2726.11707100184741
      besselzero(869) = 2729.25866360969582
      besselzero(870) = 2732.40025620216231
      besselzero(871) = 2735.54184880809536
      besselzero(872) = 2738.68344140278008
      besselzero(873) = 2741.82503401695498
      besselzero(874) = 2744.9666266056929
      besselzero(875) = 2748.10821921969283
      besselzero(876) = 2751.24981180907646
      besselzero(877) = 2754.3914044229143
      besselzero(878) = 2757.53299701292808
      besselzero(879) = 2760.6745896266136
      besselzero(880) = 2763.81618221724334
      besselzero(881) = 2766.95777483078808
      besselzero(882) = 2770.09936742202089
      besselzero(883) = 2773.24096003543507
      besselzero(884) = 2776.38255262725586
      besselzero(885) = 2779.52414524055013
      besselzero(886) = 2782.6657378329469
      besselzero(887) = 2785.8073304461306
      besselzero(888) = 2788.94892303909047
      besselzero(889) = 2792.09051565217292
      besselzero(890) = 2795.23210824568391
      besselzero(891) = 2798.37370085867398
      besselzero(892) = 2801.51529345272365
      besselzero(893) = 2804.65688606563023
      besselzero(894) = 2807.79847866020748
      besselzero(895) = 2810.94007127303946
      besselzero(896) = 2814.08166386813221
      besselzero(897) = 2817.22325648089686
      besselzero(898) = 2820.36484907649498
      besselzero(899) = 2823.50644168920172
      besselzero(900) = 2826.64803428529598
      besselzero(901) = 2829.78962689795032
      besselzero(902) = 2832.93121949451798
      besselzero(903) = 2836.07281210710979
      besselzero(904) = 2839.21440470416364
      besselzero(905) = 2842.35599731682508
      besselzero(906) = 2845.49758991447898
      besselzero(907) = 2848.63918252679493
      besselzero(908) = 2851.78077512341943
      besselzero(909) = 2854.92236773374941
      besselzero(910) = 2858.06396033260057
      besselzero(911) = 2861.20555295129452
      besselzero(912) = 2864.34714554405331
      besselzero(913) = 2867.48873816253713
      besselzero(914) = 2870.63033075591977
      besselzero(915) = 2873.77192337420517
      besselzero(916) = 2876.91351596819711
      besselzero(917) = 2880.05510858629509
      besselzero(918) = 2883.19670118088247
      besselzero(919) = 2886.33829379880513
      besselzero(920) = 2889.47988639397301
      besselzero(921) = 2892.62147901173136
      besselzero(922) = 2895.76307160746679
      besselzero(923) = 2898.90466422507131
      besselzero(924) = 2902.04625682136185
      besselzero(925) = 2905.18784943882213
      besselzero(926) = 2908.32944203565446
      besselzero(927) = 2911.47103465298062
      besselzero(928) = 2914.6126272503432
      besselzero(929) = 2917.75421986754501
      besselzero(930) = 2920.89581246542486
      besselzero(931) = 2924.03740508251175
      besselzero(932) = 2927.17899768089839
      besselzero(933) = 2930.32059029787781
      besselzero(934) = 2933.46218289675952
      besselzero(935) = 2936.6037755136416
      besselzero(936) = 2939.74536811300701
      besselzero(937) = 2942.88696072979902
      besselzero(938) = 2946.02855332963908
      besselzero(939) = 2949.1701459463492
      besselzero(940) = 2952.31173854665218
      besselzero(941) = 2955.4533311632877
      besselzero(942) = 2958.59492376404489
      besselzero(943) = 2961.73651638061308
      besselzero(944) = 2964.87810898181471
      besselzero(945) = 2968.01970159832162
      besselzero(946) = 2971.16129419995811
      besselzero(947) = 2974.30288681641297
      besselzero(948) = 2977.44447941847525
      besselzero(949) = 2980.58607203488268
      besselzero(950) = 2983.72766463736169
      besselzero(951) = 2986.86925725372898
      besselzero(952) = 2990.01084985661602
      besselzero(953) = 2993.1524424729492
      besselzero(954) = 2996.29403507623592
      besselzero(955) = 2999.43562769254068
      besselzero(956) = 3002.57722029621954
      besselzero(957) = 3005.71881291250164
      besselzero(958) = 3008.86040551656508
      besselzero(959) = 3012.00199813282898
      besselzero(960) = 3015.143590737269
      besselzero(961) = 3018.28518335352091
      besselzero(962) = 3021.42677595832997
      besselzero(963) = 3024.56836857457478
      besselzero(964) = 3027.70996117974577
      besselzero(965) = 3030.85155379598836
      besselzero(966) = 3033.99314640151462
      besselzero(967) = 3037.13473901775987
      besselzero(968) = 3040.27633162363297
      besselzero(969) = 3043.41792423988533
      besselzero(970) = 3046.55951684610126
      besselzero(971) = 3049.70110946236427
      besselzero(972) = 3052.84270206891417
      besselzero(973) = 3055.98429468519361
      besselzero(974) = 3059.12588729207302
      besselzero(975) = 3062.267479908372
      besselzero(976) = 3065.40907251557337
      besselzero(977) = 3068.55066513189544
      besselzero(978) = 3071.69225773941434
      besselzero(979) = 3074.8338503557635
      besselzero(980) = 3077.97544296359326
      besselzero(981) = 3081.11703557997307
      besselzero(982) = 3084.25862818810792
      besselzero(983) = 3087.40022080452192
      besselzero(984) = 3090.54181341295741
      besselzero(985) = 3093.68340602941138
      besselzero(986) = 3096.82499863814308
      besselzero(987) = 3099.96659125462858
      besselzero(988) = 3103.10818386362363
      besselzero(989) = 3106.24977648013977
      besselzero(990) = 3109.39136908951629
      besselzero(991) = 3112.53296170630733
      besselzero(992) = 3115.67455431584905
      besselzero(993) = 3118.81614693106268
      besselzero(994) = 3121.95773953739408
      besselzero(995) = 3125.09933215131185
      besselzero(996) = 3128.24092476527537
      besselzero(997) = 3131.38251737928507
      besselzero(998) = 3134.52410999334096
      besselzero(999) = 3137.66570260744526
      besselzero(1000) = 3140.80729522159751
      besselzero(1001) = 3143.94888783579995
      besselzero(1002) = 3147.09048045005257
      besselzero(1003) = 3150.2320730643567
      besselzero(1004) = 3153.37366567871324
      besselzero(1005) = 3156.51525829312218
      besselzero(1006) = 3159.65685090758574
      besselzero(1007) = 3162.79844352210481
      besselzero(1008) = 3165.9400361366794
      besselzero(1009) = 3169.08162875131083
      besselzero(1010) = 3172.22322136599998
      besselzero(1011) = 3175.36481398074776
      besselzero(1012) = 3178.5064065955555
      besselzero(1013) = 3181.64799921042362
      besselzero(1014) = 3184.78959182535259
      besselzero(1015) = 3187.93118444034373
      besselzero(1016) = 3191.07277705539882
      besselzero(1017) = 3194.21436967051786
      besselzero(1018) = 3197.35596228570218
      besselzero(1019) = 3200.49755490095222
      besselzero(1020) = 3203.63914751626888
      besselzero(1021) = 3206.78074013165304
      besselzero(1022) = 3209.92233274710647
      besselzero(1023) = 3213.06392536262919
      besselzero(1024) = 3216.20551797822296
      besselfactor(1) = 3.71038068570940149
      besselfactor(2) = 8.63705996767315831
      besselfactor(3) = 13.5710343180866699
      besselfactor(4) = 18.5056422933907036
      besselfactor(5) = 23.4403774800316755
      besselfactor(6) = 28.3751507697616389
      besselfactor(7) = 33.3099385387485558
      besselfactor(8) = 38.2447327558289487
      besselfactor(9) = 43.17953018986592
      besselfactor(10) = 48.1143293719934028
      besselfactor(11) = 53.0491295694953457
      besselfactor(12) = 57.9839303891639268
      besselfactor(13) = 62.918731607151157
      besselfactor(14) = 67.8535330896593791
      besselfactor(15) = 72.7883347533742419
      besselfactor(16) = 77.7231365445768496
      besselfactor(17) = 82.6579384275732031
      besselfactor(18) = 87.5927403780174707
      besselfactor(19) = 92.5275423789166495
      besselfactor(20) = 97.4623444181659515
      besselfactor(21) = 102.397146486983792
      besselfactor(22) = 107.331948578894964
      besselfactor(23) = 112.266750689053695
      besselfactor(24) = 117.201552813784504
      besselfactor(25) = 122.136354950265824
      besselfactor(26) = 127.071157096306741
      besselfactor(27) = 132.00595925018892
      besselfactor(28) = 136.940761410551625
      besselfactor(29) = 141.875563576308056
      besselfactor(30) = 146.810365746582611
      besselfactor(31) = 151.745167920665302
      besselfactor(32) = 156.679970097975518
      besselfactor(33) = 161.614772278035756
      besselfactor(34) = 166.549574464209086
      besselfactor(35) = 171.484376648335379
      besselfactor(36) = 176.419178834246892
      besselfactor(37) = 181.353981021707255
      besselfactor(38) = 186.288783210517223
      besselfactor(39) = 191.223585400506657
      besselfactor(40) = 196.158387591531564
      besselfactor(41) = 201.093189783467574
      besselfactor(42) = 206.027991976208735
      besselfactor(43) = 210.962794169662704
      besselfactor(44) = 215.897596363749677
      besselfactor(45) = 220.832398558401266
      besselfactor(46) = 225.767200753556718
      besselfactor(47) = 230.702002949163098
      besselfactor(48) = 235.63680514517471
      besselfactor(49) = 240.571607341550386
      besselfactor(50) = 245.506409538254732
      besselfactor(51) = 250.441211735255375
      besselfactor(52) = 255.376013932524958
      besselfactor(53) = 260.310816130037858
      besselfactor(54) = 265.24561832777227
      besselfactor(55) = 270.180420525708253
      besselfactor(56) = 275.115222723828267
      besselfactor(57) = 280.050024922115969
      besselfactor(58) = 284.98482712055706
      besselfactor(59) = 289.919629319139193
      besselfactor(60) = 294.854431517850414
      besselfactor(61) = 299.789233716680225
      besselfactor(62) = 304.724035915619495
      besselfactor(63) = 309.65883811465904
      besselfactor(64) = 314.593640313792377
      besselfactor(65) = 319.528442513010757
      besselfactor(66) = 324.463244712309029
      besselfactor(67) = 329.398046911680619
      besselfactor(68) = 334.332849111120555
      besselfactor(69) = 339.26765131062373
      besselfactor(70) = 344.202453510186012
      besselfactor(71) = 349.13725570980274
      besselfactor(72) = 354.072057909470805
      besselfactor(73) = 359.006860109186209
      besselfactor(74) = 363.941662308945712
      besselfactor(75) = 368.876464508746959
      besselfactor(76) = 373.811266708586975
      besselfactor(77) = 378.746068908463495
      besselfactor(78) = 383.680871108373366
      besselfactor(79) = 388.615673308315621
      besselfactor(80) = 393.550475508286989
      besselfactor(81) = 398.485277708287188
      besselfactor(82) = 403.420079908313589
      besselfactor(83) = 408.354882108364147
      besselfactor(84) = 413.28968430843842
      besselfactor(85) = 418.224486508534454
      besselfactor(86) = 423.159288708651715
      besselfactor(87) = 428.094090908787983
      besselfactor(88) = 433.028893108942902
      besselfactor(89) = 437.963695309114608
      besselfactor(90) = 442.898497509303457
      besselfactor(91) = 447.833299709507137
      besselfactor(92) = 452.768101909726095
      besselfactor(93) = 457.702904109958819
      besselfactor(94) = 462.637706310204688
      besselfactor(95) = 467.57250851046308
      besselfactor(96) = 472.507310710733285
      besselfactor(97) = 477.442112911014326
      besselfactor(98) = 482.376915111306292
      besselfactor(99) = 487.311717311608738
      besselfactor(100) = 492.246519511920688
      besselfactor(101) = 497.181321712241253
      besselfactor(102) = 502.116123912570788
      besselfactor(103) = 507.050926112909028
      besselfactor(104) = 511.985728313254462
      besselfactor(105) = 516.920530513608422
      besselfactor(106) = 521.855332713969133
      besselfactor(107) = 526.790134914336327
      besselfactor(108) = 531.72493711471036
      besselfactor(109) = 536.659739315090167
      besselfactor(110) = 541.594541515476457
      besselfactor(111) = 546.529343715867899
      besselfactor(112) = 551.464145916265024
      besselfactor(113) = 556.398948116667035
      besselfactor(114) = 561.333750317074997
      besselfactor(115) = 566.268552517486778
      besselfactor(116) = 571.203354717903799
      besselfactor(117) = 576.138156918324462
      besselfactor(118) = 581.072959118749743
      besselfactor(119) = 586.007761319178843
      besselfactor(120) = 590.942563519611674
      besselfactor(121) = 595.877365720048324
      besselfactor(122) = 600.812167920488615
      besselfactor(123) = 605.74697012093237
      besselfactor(124) = 610.681772321378879
      besselfactor(125) = 615.616574521828674
      besselfactor(126) = 620.551376722281577
      besselfactor(127) = 625.486178922737146
      besselfactor(128) = 630.420981123196089
      besselfactor(129) = 635.355783323656897
      besselfactor(130) = 640.290585524120281
      besselfactor(131) = 645.225387724586862
      besselfactor(132) = 650.160189925055043
      besselfactor(133) = 655.094992125525799
      besselfactor(134) = 660.02979432599913
      besselfactor(135) = 664.964596526474327
      besselfactor(136) = 669.899398726951211
      besselfactor(137) = 674.834200927430138
      besselfactor(138) = 679.76900312791173
      besselfactor(139) = 684.703805328394743
      besselfactor(140) = 689.638607528878555
      besselfactor(141) = 694.573409729365032
      besselfactor(142) = 699.508211929853285
      besselfactor(143) = 704.443014130342515
      besselfactor(144) = 709.377816330833432
      besselfactor(145) = 714.312618531325504
      besselfactor(146) = 719.247420731819886
      besselfactor(147) = 724.182222932315333
      besselfactor(148) = 729.117025132812024
      besselfactor(149) = 734.051827333309603
      besselfactor(150) = 738.986629533808958
      besselfactor(151) = 743.921431734309646
      besselfactor(152) = 748.856233934810955
      besselfactor(153) = 753.791036135313419
      besselfactor(154) = 758.725838335816682
      besselfactor(155) = 763.660640536321277
      besselfactor(156) = 768.595442736827116
      besselfactor(157) = 773.530244937334288
      besselfactor(158) = 778.465047137841814
      besselfactor(159) = 783.399849338350673
      besselfactor(160) = 788.334651538860509
      besselfactor(161) = 793.269453739370344
      besselfactor(162) = 798.204255939881335
      besselfactor(163) = 803.139058140393302
      besselfactor(164) = 808.073860340906513
      besselfactor(165) = 813.008662541420434
      besselfactor(166) = 817.943464741934711
      besselfactor(167) = 822.878266942449521
      besselfactor(168) = 827.813069142965396
      besselfactor(169) = 832.747871343481982
      besselfactor(170) = 837.682673543999101
      besselfactor(171) = 842.617475744517285
      besselfactor(172) = 847.552277945035115
      besselfactor(173) = 852.487080145554365
      besselfactor(174) = 857.421882346073616
      besselfactor(175) = 862.356684546594109
      besselfactor(176) = 867.291486747114781
      besselfactor(177) = 872.226288947635808
      besselfactor(178) = 877.161091148157545
      besselfactor(179) = 882.095893348679638
      besselfactor(180) = 887.030695549203151
      besselfactor(181) = 891.965497749726665
      besselfactor(182) = 896.900299950250001
      besselfactor(183) = 901.835102150774403
      besselfactor(184) = 906.769904351298983
      besselfactor(185) = 911.704706551823917
      besselfactor(186) = 916.639508752349385
      besselfactor(187) = 921.574310952875563
      besselfactor(188) = 926.509113153401564
      besselfactor(189) = 931.44391535392792
      besselfactor(190) = 936.378717554455874
      besselfactor(191) = 941.313519754982941
      besselfactor(192) = 946.24832195550983
      besselfactor(193) = 951.183124156038672
      besselfactor(194) = 956.11792635656716
      besselfactor(195) = 961.05272855709618
      besselfactor(196) = 965.987530757625378
      besselfactor(197) = 970.922332958154221
      besselfactor(198) = 975.85713515868413
      besselfactor(199) = 980.791937359214571
      besselfactor(200) = 985.726739559744657
      besselfactor(201) = 990.661541760275987
      besselfactor(202) = 995.596343960806606
      besselfactor(203) = 1000.53114616133731
      besselfactor(204) = 1005.4659483618682
      besselfactor(205) = 1010.40075056240064
      besselfactor(206) = 1015.33555276293308
      besselfactor(207) = 1020.27035496346485
      besselfactor(208) = 1025.20515716399707
      besselfactor(209) = 1030.13995936452906
      besselfactor(210) = 1035.07476156506283
      besselfactor(211) = 1040.00956376559505
      besselfactor(212) = 1044.94436596612927
      besselfactor(213) = 1049.87916816666282
      besselfactor(214) = 1054.81397036719615
      besselfactor(215) = 1059.74877256773036
      besselfactor(216) = 1064.6835747682648
      besselfactor(217) = 1069.6183769687988
      besselfactor(218) = 1074.55317916933346
      besselfactor(219) = 1079.48798136986812
      besselfactor(220) = 1084.42278357040366
      besselfactor(221) = 1089.35758577093766
      besselfactor(222) = 1094.29238797147299
      besselfactor(223) = 1099.22719017200876
      besselfactor(224) = 1104.16199237254386
      besselfactor(225) = 1109.09679457308052
      besselfactor(226) = 1114.03159677361562
      besselfactor(227) = 1118.96639897415228
      besselfactor(228) = 1123.90120117468828
      besselfactor(229) = 1128.83600337522471
      besselfactor(230) = 1133.7708055757607
      besselfactor(231) = 1138.70560777629714
      besselfactor(232) = 1143.64040997683402
      besselfactor(233) = 1148.57521217737135
      besselfactor(234) = 1153.51001437790801
      besselfactor(235) = 1158.44481657844555
      besselfactor(236) = 1163.37961877898199
      besselfactor(237) = 1168.31442097951976
      besselfactor(238) = 1173.24922318005797
      besselfactor(239) = 1178.18402538059596
      besselfactor(240) = 1183.11882758113394
      besselfactor(241) = 1188.05362978167083
      besselfactor(242) = 1192.9884319822095
      besselfactor(243) = 1197.9232341827477
      besselfactor(244) = 1202.85803638328525
      besselfactor(245) = 1207.79283858382346
      besselfactor(246) = 1212.72764078436168
      besselfactor(247) = 1217.66244298490078
      besselfactor(248) = 1222.59724518543944
      besselfactor(249) = 1227.53204738597831
      besselfactor(250) = 1232.46684958651675
      besselfactor(251) = 1237.40165178705607
      besselfactor(252) = 1242.33645398759517
      besselfactor(253) = 1247.27125618813383
      besselfactor(254) = 1252.20605838867382
      besselfactor(255) = 1257.14086058921248
      besselfactor(256) = 1262.07566278975225
      besselfactor(257) = 1267.01046499029135
      besselfactor(258) = 1271.94526719083045
      besselfactor(259) = 1276.8800693913711
      besselfactor(260) = 1281.81487159191021
      besselfactor(261) = 1286.74967379244931
      besselfactor(262) = 1291.68447599298929
      besselfactor(263) = 1296.61927819352973
      besselfactor(264) = 1301.55408039406994
      besselfactor(265) = 1306.48888259460904
      besselfactor(266) = 1311.42368479514992
      besselfactor(267) = 1316.35848699569013
      besselfactor(268) = 1321.2932891962299
      besselfactor(269) = 1326.228091396771
      besselfactor(270) = 1331.16289359731077
      besselfactor(271) = 1336.09769579785075
      besselfactor(272) = 1341.03249799839163
      besselfactor(273) = 1345.96730019893295
      besselfactor(274) = 1350.90210239947339
      besselfactor(275) = 1355.83690460001382
      besselfactor(276) = 1360.77170680055448
      besselfactor(277) = 1365.70650900109447
      besselfactor(278) = 1370.64131120163579
      besselfactor(279) = 1375.57611340217755
      besselfactor(280) = 1380.51091560271821
      besselfactor(281) = 1385.44571780325842
      besselfactor(282) = 1390.38052000380041
      besselfactor(283) = 1395.31532220434172
      besselfactor(284) = 1400.25012440488216
      besselfactor(285) = 1405.184926605423
      besselfactor(286) = 1410.11972880596463
      besselfactor(287) = 1415.05453100650546
      besselfactor(288) = 1419.98933320704656
      besselfactor(289) = 1424.92413540758944
      besselfactor(290) = 1429.85893760813054
      besselfactor(291) = 1434.79373980867137
      besselfactor(292) = 1439.72854200921318
      besselfactor(293) = 1444.66334420975499
      besselfactor(294) = 1449.59814641029716
      besselfactor(295) = 1454.53294861083808
      besselfactor(296) = 1459.46775081137953
      besselfactor(297) = 1464.40255301192134
      besselfactor(298) = 1469.33735521246262
      besselfactor(299) = 1474.27215741300621
      besselfactor(300) = 1479.20695961354731
      besselfactor(301) = 1484.14176181408858
      besselfactor(302) = 1489.07656401463128
      besselfactor(303) = 1494.01136621517309
      besselfactor(304) = 1498.94616841571455
      besselfactor(305) = 1503.88097061625618
      besselfactor(306) = 1508.81577281679906
      besselfactor(307) = 1513.75057501734078
      besselfactor(308) = 1518.68537721788303
      besselfactor(309) = 1523.62017941842537
      besselfactor(310) = 1528.55498161896719
      besselfactor(311) = 1533.48978381950988
      besselfactor(312) = 1538.42458602005249
      besselfactor(313) = 1543.3593882205951
      besselfactor(314) = 1548.29419042113674
      besselfactor(315) = 1553.22899262167934
      besselfactor(316) = 1558.16379482222195
      besselfactor(317) = 1563.09859702276376
      besselfactor(318) = 1568.03339922330708
      besselfactor(319) = 1572.96820142384952
      besselfactor(320) = 1577.90300362439195
      besselfactor(321) = 1582.8378058249342
      besselfactor(322) = 1587.77260802547681
      besselfactor(323) = 1592.70741022601818
      besselfactor(324) = 1597.64221242656168
      besselfactor(325) = 1602.57701462710366
      besselfactor(326) = 1607.51181682764681
      besselfactor(327) = 1612.44661902819004
      besselfactor(328) = 1617.38142122873362
      besselfactor(329) = 1622.3162234292765
      besselfactor(330) = 1627.25102562981831
      besselfactor(331) = 1632.18582783036181
      besselfactor(332) = 1637.1206300309038
      besselfactor(333) = 1642.05543223144667
      besselfactor(334) = 1646.99023443198955
      besselfactor(335) = 1651.92503663253269
      besselfactor(336) = 1656.85983883307486
      besselfactor(337) = 1661.79464103361818
      besselfactor(338) = 1666.72944323416061
      besselfactor(339) = 1671.66424543470419
      besselfactor(340) = 1676.59904763524707
      besselfactor(341) = 1681.53384983579101
      besselfactor(342) = 1686.46865203633433
      besselfactor(343) = 1691.40345423687499
      besselfactor(344) = 1696.33825643742035
      besselfactor(345) = 1701.27305863796181
      besselfactor(346) = 1706.20786083850495
      besselfactor(347) = 1711.14266303904827
      besselfactor(348) = 1716.07746523959185
      besselfactor(349) = 1721.01226744013509
      besselfactor(350) = 1725.94706964067823
      besselfactor(351) = 1730.88187184122084
      besselfactor(352) = 1735.81667404176336
      besselfactor(353) = 1740.75147624230686
      besselfactor(354) = 1745.68627844285124
      besselfactor(355) = 1750.62108064339305
      besselfactor(356) = 1755.55588284393522
      besselfactor(357) = 1760.49068504447987
      besselfactor(358) = 1765.42548724502417
      besselfactor(359) = 1770.36028944556589
      besselfactor(360) = 1775.29509164611054
      besselfactor(361) = 1780.22989384665253
      besselfactor(362) = 1785.16469604719745
      besselfactor(363) = 1790.09949824774068
      besselfactor(364) = 1795.03430044828249
      besselfactor(365) = 1799.9691026488267
      besselfactor(366) = 1804.90390484937108
      besselfactor(367) = 1809.83870704991396
      besselfactor(368) = 1814.77350925045826
      besselfactor(369) = 1819.70831145100149
      besselfactor(370) = 1824.64311365154472
      besselfactor(371) = 1829.57791585208759
      besselfactor(372) = 1834.51271805263118
      besselfactor(373) = 1839.44752025317441
      besselfactor(374) = 1844.382322453718
      besselfactor(375) = 1849.31712465426194
      besselfactor(376) = 1854.25192685480553
      besselfactor(377) = 1859.18672905534982
      besselfactor(378) = 1864.12153125589341
      besselfactor(379) = 1869.05633345643664
      besselfactor(380) = 1873.99113565697987
      besselfactor(381) = 1878.92593785752346
      besselfactor(382) = 1883.86074005806776
      besselfactor(383) = 1888.79554225861099
      besselfactor(384) = 1893.73034445915315
      besselfactor(385) = 1898.66514665969816
      besselfactor(386) = 1903.59994886024139
      besselfactor(387) = 1908.53475106078604
      besselfactor(388) = 1913.46955326132963
      besselfactor(389) = 1918.40435546187215
      besselfactor(390) = 1923.33915766241609
      besselfactor(391) = 1928.27395986295897
      besselfactor(392) = 1933.20876206350291
      besselfactor(393) = 1938.14356426404614
      besselfactor(394) = 1943.07836646459151
      besselfactor(395) = 1948.01316866513509
      besselfactor(396) = 1952.94797086567832
      besselfactor(397) = 1957.88277306622298
      besselfactor(398) = 1962.8175752667655
      besselfactor(399) = 1967.75237746731015
      besselfactor(400) = 1972.68717966785338
      besselfactor(401) = 1977.62198186839732
      besselfactor(402) = 1982.55678406894091
      besselfactor(403) = 1987.49158626948521
      besselfactor(404) = 1992.42638847002879
      besselfactor(405) = 1997.36119067057309
      besselfactor(406) = 2002.29599287111792
      besselfactor(407) = 2007.2307950716608
      besselfactor(408) = 2012.165597272205
      besselfactor(409) = 2017.10039947274966
      besselfactor(410) = 2022.03520167329253
      besselfactor(411) = 2026.97000387383496
      besselfactor(412) = 2031.90480607437917
      besselfactor(413) = 2036.83960827492382
      besselfactor(414) = 2041.77441047546714
      besselfactor(415) = 2046.70921267601091
      besselfactor(416) = 2051.64401487655512
      besselfactor(417) = 2056.57881707709888
      besselfactor(418) = 2061.51361927764309
      besselfactor(419) = 2066.44842147818685
      besselfactor(420) = 2071.38322367873106
      besselfactor(421) = 2076.31802587927527
      besselfactor(422) = 2081.25282807981815
      besselfactor(423) = 2086.18763028036457
      besselfactor(424) = 2091.12243248090745
      besselfactor(425) = 2096.05723468145166
      besselfactor(426) = 2100.99203688199498
      besselfactor(427) = 2105.92683908253919
      besselfactor(428) = 2110.86164128308295
      besselfactor(429) = 2115.79644348362716
      besselfactor(430) = 2120.73124568417048
      besselfactor(431) = 2125.66604788471469
      besselfactor(432) = 2130.60085008525934
      besselfactor(433) = 2135.53565228580133
      besselfactor(434) = 2140.47045448634687
      besselfactor(435) = 2145.40525668689152
      besselfactor(436) = 2150.34005888743484
      besselfactor(437) = 2155.27486108798038
      besselfactor(438) = 2160.20966328852326
      besselfactor(439) = 2165.14446548906747
      besselfactor(440) = 2170.07926768961123
      besselfactor(441) = 2175.01406989015544
      besselfactor(442) = 2179.9488720906992
      besselfactor(443) = 2184.88367429124341
      besselfactor(444) = 2189.81847649178806
      besselfactor(445) = 2194.75327869233227
      besselfactor(446) = 2199.68808089287604
      besselfactor(447) = 2204.62288309342113
      besselfactor(448) = 2209.55768529396578
      besselfactor(449) = 2214.4924874945091
      besselfactor(450) = 2219.42728969505198
      besselfactor(451) = 2224.36209189559619
      besselfactor(452) = 2229.29689409614129
      besselfactor(453) = 2234.23169629668638
      besselfactor(454) = 2239.16649849722748
      besselfactor(455) = 2244.10130069777213
      besselfactor(456) = 2249.03610289831856
      besselfactor(457) = 2253.97090509886233
      besselfactor(458) = 2258.90570729940609
      besselfactor(459) = 2263.8405094999503
      besselfactor(460) = 2268.77531170049318
      besselfactor(461) = 2273.71011390103827
      besselfactor(462) = 2278.64491610158204
      besselfactor(463) = 2283.57971830212624
      besselfactor(464) = 2288.5145205026709
      besselfactor(465) = 2293.44932270321555
      besselfactor(466) = 2298.3841249037602
      besselfactor(467) = 2303.31892710430397
      besselfactor(468) = 2308.2537293048464
      besselfactor(469) = 2313.18853150539283
      besselfactor(470) = 2318.12333370593526
      besselfactor(471) = 2323.05813590648124
      besselfactor(472) = 2327.99293810702545
      besselfactor(473) = 2332.92774030756788
      besselfactor(474) = 2337.86254250811165
      besselfactor(475) = 2342.7973447086563
      besselfactor(476) = 2347.73214690920228
      besselfactor(477) = 2352.6669491097456
      besselfactor(478) = 2357.60175131029026
      besselfactor(479) = 2362.53655351083313
      besselfactor(480) = 2367.47135571137823
      besselfactor(481) = 2372.40615791192421
      besselfactor(482) = 2377.34096011246887
      besselfactor(483) = 2382.27576231301174
      besselfactor(484) = 2387.21056451355639
      besselfactor(485) = 2392.14536671409927
      besselfactor(486) = 2397.08016891464392
      besselfactor(487) = 2402.01497111518813
      besselfactor(488) = 2406.94977331573323
      besselfactor(489) = 2411.88457551627699
      besselfactor(490) = 2416.81937771682076
      besselfactor(491) = 2421.75417991736541
      besselfactor(492) = 2426.68898211791006
      besselfactor(493) = 2431.62378431845338
      besselfactor(494) = 2436.55858651899803
      besselfactor(495) = 2441.49338871954402
      besselfactor(496) = 2446.42819092008556
      besselfactor(497) = 2451.36299312063021
      besselfactor(498) = 2456.2977953211762
      besselfactor(499) = 2461.23259752172129
      besselfactor(500) = 2466.16739972226373
      besselfactor(501) = 2471.10220192280838
      besselfactor(502) = 2476.03700412335214
      besselfactor(503) = 2480.97180632389768
      besselfactor(504) = 2485.90660852444056
      besselfactor(505) = 2490.84141072498477
      besselfactor(506) = 2495.77621292553075
      besselfactor(507) = 2500.71101512607496
      besselfactor(508) = 2505.64581732661961
      besselfactor(509) = 2510.58061952716382
      besselfactor(510) = 2515.51542172770714
      besselfactor(511) = 2520.45022392825135
      besselfactor(512) = 2525.38502612879467
      besselfactor(513) = 2530.31982832934065
      besselfactor(514) = 2535.25463052988531
      besselfactor(515) = 2540.18943273042996
      besselfactor(516) = 2545.12423493097461
      besselfactor(517) = 2550.05903713151749
      besselfactor(518) = 2554.99383933206258
      besselfactor(519) = 2559.92864153260768
      besselfactor(520) = 2564.86344373315011
      besselfactor(521) = 2569.7982459336961
      besselfactor(522) = 2574.73304813424075
      besselfactor(523) = 2579.66785033478407
      besselfactor(524) = 2584.60265253533005
      besselfactor(525) = 2589.53745473587382
      besselfactor(526) = 2594.47225693641669
      besselfactor(527) = 2599.40705913696312
      besselfactor(528) = 2604.34186133750689
      besselfactor(529) = 2609.27666353804976
      besselfactor(530) = 2614.21146573859575
      besselfactor(531) = 2619.14626793913907
      besselfactor(532) = 2624.08107013968372
      besselfactor(533) = 2629.01587234022704
      besselfactor(534) = 2633.95067454077302
      besselfactor(535) = 2638.88547674131723
      besselfactor(536) = 2643.82027894186188
      besselfactor(537) = 2648.75508114240654
      besselfactor(538) = 2653.68988334295075
      besselfactor(539) = 2658.62468554349407
      besselfactor(540) = 2663.55948774403827
      besselfactor(541) = 2668.49428994458337
      besselfactor(542) = 2673.42909214512669
      besselfactor(543) = 2678.3638943456709
      besselfactor(544) = 2683.298696546216
      besselfactor(545) = 2688.2334987467602
      besselfactor(546) = 2693.16830094730397
      besselfactor(547) = 2698.1031031478504
      besselfactor(548) = 2703.0379053483946
      besselfactor(549) = 2707.97270754893926
      besselfactor(550) = 2712.90750974948391
      besselfactor(551) = 2717.84231195002634
      besselfactor(552) = 2722.77711415057233
      besselfactor(553) = 2727.71191635111431
      besselfactor(554) = 2732.64671855165985
      besselfactor(555) = 2737.58152075220362
      besselfactor(556) = 2742.51632295274916
      besselfactor(557) = 2747.45112515329426
      besselfactor(558) = 2752.38592735383802
      besselfactor(559) = 2757.32072955438178
      besselfactor(560) = 2762.25553175492732
      besselfactor(561) = 2767.19033395546976
      besselfactor(562) = 2772.12513615601663
      besselfactor(563) = 2777.05993835656084
      besselfactor(564) = 2781.99474055710505
      besselfactor(565) = 2786.92954275765147
      besselfactor(566) = 2791.86434495819435
      besselfactor(567) = 2796.799147158739
      besselfactor(568) = 2801.73394935928179
      besselfactor(569) = 2806.66875155982697
      besselfactor(570) = 2811.6035537603734
      besselfactor(571) = 2816.53835596091628
      besselfactor(572) = 2821.4731581614604
      besselfactor(573) = 2826.40796036200381
      besselfactor(574) = 2831.34276256255077
      besselfactor(575) = 2836.277564763094
      besselfactor(576) = 2841.21236696363901
      besselfactor(577) = 2846.14716916418419
      besselfactor(578) = 2851.08197136472654
      besselfactor(579) = 2856.01677356527261
      besselfactor(580) = 2860.95157576581691
      besselfactor(581) = 2865.88637796635978
      besselfactor(582) = 2870.8211801669039
      besselfactor(583) = 2875.75598236744944
      besselfactor(584) = 2880.69078456799445
      besselfactor(585) = 2885.62558676853697
      besselfactor(586) = 2890.56038896908376
      besselfactor(587) = 2895.4951911696277
      besselfactor(588) = 2900.42999337017182
      besselfactor(589) = 2905.36479557071647
      besselfactor(590) = 2910.29959777125935
      besselfactor(591) = 2915.23439997180631
      besselfactor(592) = 2920.16920217234989
      besselfactor(593) = 2925.10400437289561
      besselfactor(594) = 2930.03880657344026
      besselfactor(595) = 2934.97360877398207
      besselfactor(596) = 2939.90841097452815
      besselfactor(597) = 2944.8432131750728
      besselfactor(598) = 2949.77801537561479
      besselfactor(599) = 2954.71281757616122
      besselfactor(600) = 2959.64761977670711
      besselfactor(601) = 2964.58242197724874
      besselfactor(602) = 2969.51722417779464
      besselfactor(603) = 2974.45202637833929
      besselfactor(604) = 2979.38682857888448
      besselfactor(605) = 2984.3216307794286
      besselfactor(606) = 2989.25643297997378
      besselfactor(607) = 2994.19123518051844
      besselfactor(608) = 2999.12603738106309
      besselfactor(609) = 3004.0608395816073
      besselfactor(610) = 3008.99564178215195
      besselfactor(611) = 3013.93044398269394
      besselfactor(612) = 3018.8652461832417
      besselfactor(613) = 3023.80004838378324
      besselfactor(614) = 3028.73485058433056
      besselfactor(615) = 3033.66965278487343
      besselfactor(616) = 3038.60445498541764
      besselfactor(617) = 3043.53925718596318
      besselfactor(618) = 3048.47405938650783
      besselfactor(619) = 3053.40886158705071
      besselfactor(620) = 3058.34366378759714
      besselfactor(621) = 3063.27846598814135
      besselfactor(622) = 3068.21326818868467
      besselfactor(623) = 3073.1480703892311
      besselfactor(624) = 3078.08287258977486
      besselfactor(625) = 3083.01767479031996
      besselfactor(626) = 3087.95247699086417
      besselfactor(627) = 3092.88727919140749
      besselfactor(628) = 3097.82208139195169
      besselfactor(629) = 3102.75688359249723
      besselfactor(630) = 3107.69168579304011
      besselfactor(631) = 3112.62648799358609
      besselfactor(632) = 3117.56129019413119
      besselfactor(633) = 3122.49609239467762
      besselfactor(634) = 3127.43089459521961
      besselfactor(635) = 3132.3656967957647
      besselfactor(636) = 3137.30049899631069
      besselfactor(637) = 3142.23530119685401
      besselfactor(638) = 3147.17010339739644
      besselfactor(639) = 3152.1049055979442
      besselfactor(640) = 3157.03970779848663
      besselfactor(641) = 3161.97450999903306
      besselfactor(642) = 3166.90931219957594
      besselfactor(643) = 3171.84411440012059
      besselfactor(644) = 3176.77891660066613
      besselfactor(645) = 3181.71371880121034
      besselfactor(646) = 3186.64852100175366
      besselfactor(647) = 3191.58332320229964
      besselfactor(648) = 3196.51812540284563
      besselfactor(649) = 3201.45292760339028
      besselfactor(650) = 3206.3877298039336
      besselfactor(651) = 3211.32253200447781
      besselfactor(652) = 3216.25733418877768
      besselfactor(653) = 3221.19213642157919
      besselfactor(654) = 3226.12693859039234
      besselfactor(655) = 3231.06174082215025
      besselfactor(656) = 3235.99654299199102
      besselfactor(657) = 3240.93134522273862
      besselfactor(658) = 3245.86614739357282
      besselfactor(659) = 3250.8009496233412
      besselfactor(660) = 3255.73575179513863
      besselfactor(661) = 3260.67055402395888
      besselfactor(662) = 3265.60535619668934
      besselfactor(663) = 3270.54015842459611
      besselfactor(664) = 3275.47496059822318
      besselfactor(665) = 3280.40976282524443
      besselfactor(666) = 3285.34456499974414
      besselfactor(667) = 3290.27936722590874
      besselfactor(668) = 3295.21416940125266
      besselfactor(669) = 3300.14897162658504
      besselfactor(670) = 3305.08377380274254
      besselfactor(671) = 3310.01857602727733
      besselfactor(672) = 3314.95337820422397
      besselfactor(673) = 3319.88818042798206
      besselfactor(674) = 3324.82298260569786
      besselfactor(675) = 3329.75778482870055
      besselfactor(676) = 3334.69258700712734
      besselfactor(677) = 3339.62738922935864
      besselfactor(678) = 3344.56219140852883
      besselfactor(679) = 3349.49699363034847
      besselfactor(680) = 3354.43179581058581
      besselfactor(681) = 3359.36659803085913
      besselfactor(682) = 3364.30140020768231
      besselfactor(683) = 3369.23620242190446
      besselfactor(684) = 3374.17100460427121
      besselfactor(685) = 3379.10580684091144
      besselfactor(686) = 3384.04060900594583
      besselfactor(687) = 3388.97541124142165
      besselfactor(688) = 3393.91021340760579
      besselfactor(689) = 3398.84501564195184
      besselfactor(690) = 3403.77981780924443
      besselfactor(691) = 3408.71462004249536
      besselfactor(692) = 3413.64942221086398
      besselfactor(693) = 3418.5842244430602
      besselfactor(694) = 3423.51902661246887
      besselfactor(695) = 3428.45382884364147
      besselfactor(696) = 3433.38863101405555
      besselfactor(697) = 3438.32343324424005
      besselfactor(698) = 3443.25823541562803
      besselfactor(699) = 3448.19303764485463
      besselfactor(700) = 3453.12783981718496
      besselfactor(701) = 3458.06264204548474
      besselfactor(702) = 3462.99744421872413
      besselfactor(703) = 3467.93224644612729
      besselfactor(704) = 3472.86704862025041
      besselfactor(705) = 3477.80185084678672
      besselfactor(706) = 3482.73665302176471
      besselfactor(707) = 3487.67145524745947
      besselfactor(708) = 3492.60625742326169
      besselfactor(709) = 3497.54105964814598
      besselfactor(710) = 3502.47586182474535
      besselfactor(711) = 3507.41066404884671
      besselfactor(712) = 3512.34546622622634
      besselfactor(713) = 3517.28026844956298
      besselfactor(714) = 3522.21507062765937
      besselfactor(715) = 3527.14987285021219
      besselfactor(716) = 3532.08467502906486
      besselfactor(717) = 3537.01947725120202
      besselfactor(718) = 3541.95427943113605
      besselfactor(719) = 3546.88908165170513
      besselfactor(720) = 3551.82388382816904
      besselfactor(721) = 3556.7586860425977
      besselfactor(722) = 3561.69348822469
      besselfactor(723) = 3566.62829046187201
      besselfactor(724) = 3571.56309262637439
      besselfactor(725) = 3576.49789486237779
      besselfactor(726) = 3581.43269702803968
      besselfactor(727) = 3586.36749926289777
      besselfactor(728) = 3591.30230142968631
      besselfactor(729) = 3596.23710366343729
      besselfactor(730) = 3601.17190583131741
      besselfactor(731) = 3606.10670806399502
      besselfactor(732) = 3611.04151023292718
      besselfactor(733) = 3615.97631246457007
      besselfactor(734) = 3620.91111463452275
      besselfactor(735) = 3625.84591686516111
      besselfactor(736) = 3630.78071903610055
      besselfactor(737) = 3635.71552126576813
      besselfactor(738) = 3640.65032343766414
      besselfactor(739) = 3645.58512566639203
      besselfactor(740) = 3650.51992783920686
      besselfactor(741) = 3655.45473006703014
      besselfactor(742) = 3660.38953224074115
      besselfactor(743) = 3665.32433446768113
      besselfactor(744) = 3670.259136642259
      besselfactor(745) = 3675.19393886834766
      besselfactor(746) = 3680.1287410437653
      besselfactor(747) = 3685.0635432690284
      besselfactor(748) = 3689.99834544525473
      besselfactor(749) = 3694.93314766971625
      besselfactor(750) = 3699.86794984672595
      besselfactor(751) = 3704.80275207044673
      besselfactor(752) = 3709.73755424826734
      besselfactor(753) = 3714.67235647120741
      besselfactor(754) = 3719.60715864945745
      besselfactor(755) = 3724.54196087128508
      besselfactor(756) = 3729.47676305115161
      besselfactor(757) = 3734.4115652765546
      besselfactor(758) = 3739.34636746274851
      besselfactor(759) = 3744.28116968235169
      besselfactor(760) = 3749.21597184512256
      besselfactor(761) = 3754.1507740828286
      besselfactor(762) = 3759.08557624681671
      besselfactor(763) = 3764.0203784833246
      besselfactor(764) = 3768.95518064848911
      besselfactor(765) = 3773.88998288383703
      besselfactor(766) = 3778.82478505013841
      besselfactor(767) = 3783.75958728437211
      besselfactor(768) = 3788.69438945177483
      besselfactor(769) = 3793.6291916849214
      besselfactor(770) = 3798.56399385339482
      besselfactor(771) = 3803.49879608548846
      besselfactor(772) = 3808.43359825499447
      besselfactor(773) = 3813.36840048607417
      besselfactor(774) = 3818.30320265657832
      besselfactor(775) = 3823.23800488667675
      besselfactor(776) = 3828.17280705814724
      besselfactor(777) = 3833.10760928729408
      besselfactor(778) = 3838.0424114597007
      besselfactor(779) = 3842.97721368792455
      besselfactor(780) = 3847.91201586123587
      besselfactor(781) = 3852.8468180885703
      besselfactor(782) = 3857.7816202627595
      besselfactor(783) = 3862.71642248923186
      besselfactor(784) = 3867.651224664273
      besselfactor(785) = 3872.58602688990691
      besselfactor(786) = 3877.52082906577157
      besselfactor(787) = 3882.45563129060187
      besselfactor(788) = 3887.39043346723267
      besselfactor(789) = 3892.32523569122434
      besselfactor(790) = 3897.26003786865593
      besselfactor(791) = 3902.19484009220973
      besselfactor(792) = 3907.12964227079027
      besselfactor(793) = 3912.06444449267465
      besselfactor(794) = 3916.99924666755894
      besselfactor(795) = 3921.93404888289833
      besselfactor(796) = 3926.86885106377747
      besselfactor(797) = 3931.80365330337266
      besselfactor(798) = 3936.73845546549614
      besselfactor(799) = 3941.67325770383882
      besselfactor(800) = 3946.60805986719864
      besselfactor(801) = 3951.54286210432346
      besselfactor(802) = 3956.47766426888125
      besselfactor(803) = 3961.41246650482656
      besselfactor(804) = 3966.34726867054521
      besselfactor(805) = 3971.28207090535312
      besselfactor(806) = 3976.21687307218873
      besselfactor(807) = 3981.1516753058953
      besselfactor(808) = 3986.08647747381805
      besselfactor(809) = 3991.02127970645526
      besselfactor(810) = 3995.95608187542481
      besselfactor(811) = 4000.89088410703119
      besselfactor(812) = 4005.8256862770154
      besselfactor(813) = 4010.76048850762489
      besselfactor(814) = 4015.69529067859143
      besselfactor(815) = 4020.63009290823636
      besselfactor(816) = 4025.56489508015324
      besselfactor(817) = 4030.49969730885937
      besselfactor(818) = 4035.43449948169996
      besselfactor(819) = 4040.36930170949837
      besselfactor(820) = 4045.30410388322714
      besselfactor(821) = 4050.23890611015247
      besselfactor(822) = 4055.17370828474543
      besselfactor(823) = 4060.10851051082433
      besselfactor(824) = 4065.04331268625396
      besselfactor(825) = 4069.97811491150862
      besselfactor(826) = 4074.91291708771719
      besselfactor(827) = 4079.8477193121263
      besselfactor(828) = 4084.782521489144
      besselfactor(829) = 4089.71732371310814
      besselfactor(830) = 4094.6521258913
      besselfactor(831) = 4099.58692811356151
      besselfactor(832) = 4104.52173028799105
      besselfactor(833) = 4109.45653250359654
      besselfactor(834) = 4114.39133468412077
      besselfactor(835) = 4119.32613692441318
      besselfactor(836) = 4124.26093908585844
      besselfactor(837) = 4129.19574132486744
      besselfactor(838) = 4134.13054348756769
      besselfactor(839) = 4139.06534572534212
      besselfactor(840) = 4144.0001478892583
      besselfactor(841) = 4148.93495012583457
      besselfactor(842) = 4153.86975229093202
      besselfactor(843) = 4158.80455452635011
      besselfactor(844) = 4163.73935669258621
      besselfactor(845) = 4168.67415892688608
      besselfactor(846) = 4173.60896109421997
      besselfactor(847) = 4178.54376332743537
      besselfactor(848) = 4183.47856549583774
      besselfactor(849) = 4188.41336772800509
      besselfactor(850) = 4193.34816989744307
      besselfactor(851) = 4198.28297212859169
      besselfactor(852) = 4203.21777429902443
      besselfactor(853) = 4208.15257652919072
      besselfactor(854) = 4213.08737870059424
      besselfactor(855) = 4218.02218092980574
      besselfactor(856) = 4222.9569831021454
      besselfactor(857) = 4227.89178533043941
      besselfactor(858) = 4232.82658750368412
      besselfactor(859) = 4237.76138973108373
      besselfactor(860) = 4242.69619190520597
      besselfactor(861) = 4247.63099413174494
      besselfactor(862) = 4252.56579630672248
      besselfactor(863) = 4257.50059853242923
      besselfactor(864) = 4262.43540070819459
      besselfactor(865) = 4267.37020293303448
      besselfactor(866) = 4272.30500510963029
      besselfactor(867) = 4277.23980733401631
      besselfactor(868) = 4282.17460951180139
      besselfactor(869) = 4287.1094117344608
      besselfactor(870) = 4292.04421390841339
      besselfactor(871) = 4296.97901612429067
      besselfactor(872) = 4301.91381830445429
      besselfactor(873) = 4306.84862054546435
      besselfactor(874) = 4311.78342270619463
      besselfactor(875) = 4316.71822494591062
      besselfactor(876) = 4321.65302710791898
      besselfactor(877) = 4326.5878293463782
      besselfactor(878) = 4331.52263150962113
      besselfactor(879) = 4336.45743374685999
      besselfactor(880) = 4341.39223591130196
      besselfactor(881) = 4346.32703814736665
      besselfactor(882) = 4351.26184031296503
      besselfactor(883) = 4356.19664254789018
      besselfactor(884) = 4361.13144471461211
      besselfactor(885) = 4366.06624694843237
      besselfactor(886) = 4371.00104911623787
      besselfactor(887) = 4375.93585134899143
      besselfactor(888) = 4380.87065351785032
      besselfactor(889) = 4385.80545574956914
      besselfactor(890) = 4390.74025791944411
      besselfactor(891) = 4395.67506015015929
      besselfactor(892) = 4400.60986232102102
      besselfactor(893) = 4405.5446645507672
      besselfactor(894) = 4410.47946672258195
      besselfactor(895) = 4415.41426895139288
      besselfactor(896) = 4420.34907112412689
      besselfactor(897) = 4425.28387335202921
      besselfactor(898) = 4430.21867552565407
      besselfactor(899) = 4435.1534777526842
      besselfactor(900) = 4440.08827992717947
      besselfactor(901) = 4445.0230821533605
      besselfactor(902) = 4449.95788432865868
      besselfactor(903) = 4454.89268655395687
      besselfactor(904) = 4459.82748873010237
      besselfactor(905) = 4464.7622909549387
      besselfactor(906) = 4469.69709313229479
      besselfactor(907) = 4474.63189535536809
      besselfactor(908) = 4479.56669752881886
      besselfactor(909) = 4484.50149974498657
      besselfactor(910) = 4489.43630192475939
      besselfactor(911) = 4494.37110416654483
      besselfactor(912) = 4499.30590632651927
      besselfactor(913) = 4504.24070856698044
      besselfactor(914) = 4509.17551072825251
      besselfactor(915) = 4514.11031296743026
      besselfactor(916) = 4519.04511512996265
      besselfactor(917) = 4523.97991736790405
      besselfactor(918) = 4528.9147195316568
      besselfactor(919) = 4533.84952176839917
      besselfactor(920) = 4538.78432393333319
      besselfactor(921) = 4543.71912616891382
      besselfactor(922) = 4548.65392833498471
      besselfactor(923) = 4553.58873056944624
      besselfactor(924) = 4558.52353273662291
      besselfactor(925) = 4563.45833496999286
      besselfactor(926) = 4568.39313713824246
      besselfactor(927) = 4573.32793937055903
      besselfactor(928) = 4578.26274153984424
      besselfactor(929) = 4583.19754377114385
      besselfactor(930) = 4588.13234594143182
      besselfactor(931) = 4593.06714817174555
      besselfactor(932) = 4598.00195034300251
      besselfactor(933) = 4602.93675257236323
      besselfactor(934) = 4607.87155474455634
      besselfactor(935) = 4612.80635697299334
      besselfactor(936) = 4617.74115914609595
      besselfactor(937) = 4622.67596137363678
      besselfactor(938) = 4627.61076354762046
      besselfactor(939) = 4632.54556577429532
      besselfactor(940) = 4637.48036794913254
      besselfactor(941) = 4642.41517017497251
      besselfactor(942) = 4647.34997235062685
      besselfactor(943) = 4652.28477457565681
      besselfactor(944) = 4657.21957675211406
      besselfactor(945) = 4662.15437897635531
      besselfactor(946) = 4667.08918115358351
      besselfactor(947) = 4672.02398337706803
      besselfactor(948) = 4676.95878555504229
      besselfactor(949) = 4681.89358777779496
      besselfactor(950) = 4686.82838995649397
      besselfactor(951) = 4691.76319217853077
      besselfactor(952) = 4696.69799435792612
      besselfactor(953) = 4701.63279657927635
      besselfactor(954) = 4706.56759875935382
      besselfactor(955) = 4711.50240098003703
      besselfactor(956) = 4716.43720316076731
      besselfactor(957) = 4721.37200538080747
      besselfactor(958) = 4726.30680756217191
      besselfactor(959) = 4731.24160978158414
      besselfactor(960) = 4736.17641196356675
      besselfactor(961) = 4741.11121418237324
      besselfactor(962) = 4746.04601636495004
      besselfactor(963) = 4750.98081858317389
      besselfactor(964) = 4755.91562076632445
      besselfactor(965) = 4760.85042298397986
      besselfactor(966) = 4765.78522516769353
      besselfactor(967) = 4770.72002738479739
      besselfactor(968) = 4775.65482956904663
      besselfactor(969) = 4780.58963178562379
      besselfactor(970) = 4785.52443397039617
      besselfactor(971) = 4790.45923618646174
      besselfactor(972) = 4795.39403837173239
      besselfactor(973) = 4800.32884058730324
      besselfactor(974) = 4805.2636427730615
      besselfactor(975) = 4810.19844498814919
      besselfactor(976) = 4815.13324717438706
      besselfactor(977) = 4820.06804938901201
      besselfactor(978) = 4825.00285157570285
      besselfactor(979) = 4829.93765378987394
      besselfactor(980) = 4834.87245597701509
      besselfactor(981) = 4839.80725819074919
      besselfactor(982) = 4844.74206037831489
      besselfactor(983) = 4849.67686259162622
      besselfactor(984) = 4854.61166477961025
      besselfactor(985) = 4859.54646699252101
      besselfactor(986) = 4864.48126918091006
      besselfactor(987) = 4869.41607139340071
      besselfactor(988) = 4874.3508735820944
      besselfactor(989) = 4879.28567579415783
      besselfactor(990) = 4884.2204779835372
      besselfactor(991) = 4889.15528019594614
      besselfactor(992) = 4894.09008238532639
      besselfactor(993) = 4899.02488459227125
      besselfactor(994) = 4903.95968677105287
      besselfactor(995) = 4908.89448897341473
      besselfactor(996) = 4913.82929117566646
      besselfactor(997) = 4918.76409337780806
      besselfactor(998) = 4923.69889557984308
      besselfactor(999) = 4928.63369778177862
      besselfactor(1000) = 4933.56849998361379
      besselfactor(1001) = 4938.50330218536015
      besselfactor(1002) = 4943.43810438700792
      besselfactor(1003) = 4948.37290658857309
      besselfactor(1004) = 4953.30770879005299
      besselfactor(1005) = 4958.24251099145651
      besselfactor(1006) = 4963.17731319277477
      besselfactor(1007) = 4968.1121153940273
      besselfactor(1008) = 4973.04691759521055
      besselfactor(1009) = 4977.98171979632453
      besselfactor(1010) = 4982.91652199737988
      besselfactor(1011) = 4987.85132419837484
      besselfactor(1012) = 4992.78612639931651
      besselfactor(1013) = 4997.72092860020755
      besselfactor(1014) = 5002.65573080104797
      besselfactor(1015) = 5007.59053300184753
      besselfactor(1016) = 5012.5253352026009
      besselfactor(1017) = 5017.46013740332053
      besselfactor(1018) = 5022.39493960400551
      besselfactor(1019) = 5027.32974180466207
      besselfactor(1020) = 5032.26454400529554
      besselfactor(1021) = 5037.19934620589704
      besselfactor(1022) = 5042.13414840648699
      besselfactor(1023) = 5047.06895060705829
      besselfactor(1024) = 5052.00375280761715
      
      return
      end
c=============================================================================

#include "top.h"
c=============================================================================
c@(#) File top_lattice.F, version $Revision: 3.21 $, $Date: 2008/06/23 23:14:36 $
c# Copyright (c) 1990-2006, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This file contains routines dealing with the accelerator lattice
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)495-2961
c============================================================================ 
      subroutine countdrfts()
      use Lattice

c     --- count drfts 

      integer(ISZ):: idrft,nndrft

      nndrft = -1
      do idrft = 0, ndrft
         if (drftzs(idrft) .ne. drftze(idrft)) nndrft = idrft
      enddo
      ndrft = nndrft

      return 
      end 
c============================================================================ 
      subroutine setupdrfts()
      use Lattice

      integer(ISZ):: idrft

      if (zlatperi > 0.) then
        if (ndrft >= 0) then
          if (drftzs(0) == drftze(0)) then
             drftzs(0) = drftzs(ndrft) - zlatperi
             drftze(0) = drftze(ndrft) - zlatperi
             drftap(0) = drftap(ndrft)
             drftox(0) = drftox(ndrft)
             drftoy(0) = drftoy(ndrft)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      drfts = .false.
      do idrft = 0, ndrft
         if (drftzs(idrft) .ne. drftze(idrft)) drfts = .true.
      enddo

      return
      end
c============================================================================ 
      subroutine countbends()
      use Lattice

c     --- count bends 

      integer(ISZ):: ib,nnbend

      nnbend = -1
      do ib = 0, nbend
         if (bendzs(ib) .ne. bendze(ib)) nnbend = ib
      enddo
      nbend = nnbend

      return 
      end 
c============================================================================ 
      subroutine setupbends()
      use Lattice

      integer(ISZ):: ib

      if (zlatperi > 0.) then
        if (nbend >= 0) then
          if (bendzs(0) == bendze(0)) then
             bendzs(0) = bendzs(nbend) - zlatperi
             bendze(0) = bendze(nbend) - zlatperi
             bendrc(0) = bendrc(nbend)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      bends = .false.
      do ib = 0, nbend
         if (bendrc(ib) .ne. 0.) bends = .true.
      enddo

      return
      end
c============================================================================ 
      subroutine countdipos()
      use Lattice

c     --- count dipoles 

      integer(ISZ):: id,nndipo

      nndipo = -1
      do id = 0, ndipo
         if (dipozs(id) .ne. dipoze(id)) nndipo = id
      enddo
      ndipo = nndipo

      return 
      end 
c============================================================================ 
      subroutine setupdipos()
      use Constant
      use Lattice
      use Beam_acc

c     --- if autoset dipoles are specified, then derive dipole locations 
c     --- consistently from bend locations and radii 

      integer(ISZ):: id

      if (diposet) then 
        do id = 0, ndipo
           if (id <= nbend) then
              if (bendrc(id) .ne. 0.) then
c                --- auto-set dipole field to match bend radius of curvature
                 if (dipoby(id) == 0. .and. dipoex(id) == 0.) then 
                    dipoby(id) = aion * amu / (zion * echarge)
     &                            * gammabar * vbeam / dvnz(bendrc(id))
                 endif
c                --- auto-set dipole starts and ends to match those of bends
                 if (dipozs(id) == dipoze(id)) then
                    dipozs(id) = bendzs(id)
                    dipoze(id) = bendze(id)
                 endif
c                --- auto-set dipole entrance & exit angles for "box dipoles"
                 if (dipotype == "box") then
                  dipota(id) = tan(.5*(dipoze(id)-dipozs(id))/dvnz(bendrc(id)))
                  dipotb(id) = tan(.5*(dipoze(id)-dipozs(id))/dvnz(bendrc(id)))
                 endif
              endif
           endif
        enddo
      endif 

      if (zlatperi > 0.) then
        if (ndipo >= 0) then
          if (dipozs(0) == dipoze(0)) then
             dipozs(0) = dipozs(ndipo) - zlatperi
             dipoze(0) = dipoze(ndipo) - zlatperi
             dipoby(0) = dipoby(ndipo)
             dipobx(0) = dipobx(ndipo)
             dipoex(0) = dipoex(ndipo)
             dipoey(0) = dipoey(ndipo)
             dipox1(0) = dipox1(ndipo)
             dipox2(0) = dipox2(ndipo)
             dipov1(0) = dipov1(ndipo)
             dipov2(0) = dipov2(ndipo)
             dipol1(0) = dipol1(ndipo)
             dipol2(0) = dipol2(ndipo)
             dipow1(0) = dipow1(ndipo)
             dipow2(0) = dipow2(ndipo)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      dipos = .false.
      do id = 0, ndipo
         if (dipoby(id).ne.0. .or. dipoex(id).ne.0.) dipos = .true.
         if (dipobx(id).ne.0. .or. dipoey(id).ne.0.) dipos = .true.
      enddo

      return
      end
c============================================================================ 
      subroutine countquads()
      use Lattice

c     --- count quads 

      integer(ISZ):: iq,nnquad

      nnquad = -1
      do iq = 0, nquad
         if (quadzs(iq) .ne. quadze(iq)) nnquad = iq
      enddo
      nquad = nnquad

      return 
      end 
c============================================================================ 
      subroutine setupquads()
      use Lattice

      integer(ISZ):: iq,l
      real(kind=8):: rnorm

      if (zlatperi > 0.) then
        if (nquad >= 0) then
          if (quadzs(0) == quadze(0)) then
             quadzs(0) = quadzs(nquad) - zlatperi
             quadze(0) = quadze(nquad) - zlatperi
             quaddb(0) = quaddb(nquad)
             quadde(0) = quadde(nquad)
             quadph(0) = quadph(nquad)
             quadts(0) = quadts(nquad)
             quaddt(0) = quaddt(nquad)
             do l = 0,ntquad
               quadet(l,0) = quadet(l,nquad)
               quadbt(l,0) = quadbt(l,nquad)
             enddo
             quadvx(0) = quadvx(nquad)
             quadvy(0) = quadvy(nquad)
             quadap(0) = quadap(nquad)
             quadrr(0) = quadrr(nquad)
             quadrl(0) = quadrl(nquad)
             quadgl(0) = quadgl(nquad)
             quadgp(0) = quadgp(nquad)
             quadpw(0) = quadpw(nquad)
             quadpa(0) = quadpa(nquad)
             quadpr(0) = quadpr(nquad)
             quadsl(0) = quadsl(nquad)
             qdelglx(0) = qdelglx(nquad)
             qdelgly(0) = qdelgly(nquad)
             qdelaxp(0) = qdelaxp(nquad)
             qdelaxm(0) = qdelaxm(nquad)
             qdelayp(0) = qdelayp(nquad)
             qdelaym(0) = qdelaym(nquad)
             qdelrxp(0) = qdelrxp(nquad)
             qdelrxm(0) = qdelrxm(nquad)
             qdelryp(0) = qdelryp(nquad)
             qdelrym(0) = qdelrym(nquad)
             qdelvxp(0) = qdelvxp(nquad)
             qdelvxm(0) = qdelvxm(nquad)
             qdelvyp(0) = qdelvyp(nquad)
             qdelvym(0) = qdelvym(nquad)
             qdeloxp(0) = qdeloxp(nquad)
             qdeloxm(0) = qdeloxm(nquad)
             qdeloyp(0) = qdeloyp(nquad)
             qdeloym(0) = qdeloym(nquad)
             qdelpwl(0) = qdelpwl(nquad)
             qdelpwr(0) = qdelpwr(nquad)
             qdelpal(0) = qdelpal(nquad)
             qdelpar(0) = qdelpar(nquad)
             qdelprl(0) = qdelprl(nquad)
             qdelprr(0) = qdelprr(nquad)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      quads = .false.
      if (ntquad > 0) quads = .true.
      do iq = 0, nquad
         if (quaddb(iq).ne.0. .or. quadde(iq).ne.0.) quads = .true.
      enddo

c     --- If qerrxrms or qerryrms not zero, then set qoff arrays using cutoff
c     --- Gaussian distribution of offsets
      if (qerrxrms .ne. 0. .or. qerryrms .ne. 0.) then
        do l = 0,iqerr-1
          qoffx(l) = 0.
          qoffy(l) = 0.
        enddo
        do l = iqerr, nqerr
          qoffx(l) = qerrxrms*rnorm()
          qoffy(l) = qerryrms*rnorm()
        enddo
      endif

      return
      end
c============================================================================ 
      subroutine countsexts()
      use Lattice

c     --- count sextupoles 

      integer(ISZ):: is,nnsext

      nnsext = -1
      do is = 0, nsext
         if (sextzs(is) .ne. sextze(is)) nnsext = is
      enddo
      nsext = nnsext

      return
      end
c============================================================================ 
      subroutine setupsexts()
      use Lattice

      integer(ISZ):: is

      if (zlatperi > 0.) then
        if (nsext >= 0) then
          if (sextzs(0) == sextze(0)) then
             sextzs(0) = sextzs(nsext) - zlatperi
             sextze(0) = sextze(nsext) - zlatperi
             sextdb(0) = sextdb(nsext)
             sextde(0) = sextde(nsext)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      sexts = .false.
      do is = 0, nsext
         if (sextdb(is).ne.0. .or. sextde(is).ne.0.) sexts = .true.
      enddo

      return
      end
c============================================================================ 
      subroutine countheles()
      use Lattice

c     --- count hard-edge multipole elements and number of multipole 
c     --- components in each element

      integer(ISZ):: ih,ii,nnhele,nnhmlt

      nnhele = -1
      nnhmlt = -1
      do ih = 0, nhele
         if (helezs(ih) .ne. heleze(ih)) nnhele = ih
         do ii = 1, nhmlt 
            if (heleae(ii,ih) .ne. 0. .or. heleep(ii,ih) .ne. 0.) helene(ih)=ii
            if (heleam(ii,ih) .ne. 0. .or. helemp(ii,ih) .ne. 0.) helenm(ih)=ii
            if (nnhmlt < helene(ih)) nnhmlt = helene(ih)
            if (nnhmlt < helenm(ih)) nnhmlt = helenm(ih)           
         enddo
      enddo
      nhele = nnhele
      nhmlt = nnhmlt  

      return 
      end 
c============================================================================ 
      subroutine setupheles()
      use Lattice

      integer(ISZ):: ii,ih

      if (zlatperi > 0.) then
        if (nhele >= 0) then
          if (helezs(0) == heleze(0)) then
             helezs(0) = helezs(nhele) - zlatperi
             heleze(0) = heleze(nhele) - zlatperi
             do ii = 1, nhmlt 
                heleae(ii,0) = heleae(ii,nhele)
                heleep(ii,0) = heleep(ii,nhele)
                heleam(ii,0) = heleam(ii,nhele)
                helemp(ii,0) = helemp(ii,nhele)
                hele_n(ii,0) = hele_n(ii,nhele)
                hele_v(ii,0) = hele_v(ii,nhele)
                helepe(ii,0) = helepe(ii,nhele)
                helepm(ii,0) = helepm(ii,nhele)
             enddo 
             helene(0) = helene(nhele)  
             helenm(0) = helenm(nhele)  
             heleox(0) = heleox(nhele)   
             heleoy(0) = heleoy(nhele) 
             helerr(0) = helerr(nhele)
             helerl(0) = helerl(nhele)
             helegl(0) = helegl(nhele)
             helegp(0) = helegp(nhele)
             helepw(0) = helepw(nhele)
             helepa(0) = helepa(nhele)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      heles = .false.
      do ih = 0, nhele
         do ii = 1,nhmlt 
            if (heleae(ii,ih).ne.0. .or. heleam(ii,ih).ne.0. .or.
     &          heleep(ii,ih).ne.0. .or. helemp(ii,ih).ne.0.) heles = .true.
         enddo 
      enddo

      return
      end
c============================================================================ 
      subroutine countaccls()
      use Lattice

c     --- count acceleration gaps 

      integer(ISZ):: ia,nnaccl

      nnaccl = -1
      do ia = 0, naccl
         if (acclzs(ia) .ne. acclze(ia)) nnaccl = ia
      enddo
      naccl = nnaccl

      return 
      end 
c============================================================================ 
      subroutine setupaccls()
      use Lattice

      integer(ISZ):: ia,l

      if (zlatperi > 0.) then
        if (naccl >= 0) then
          if (acclzs(0) == acclze(0)) then
             acclzs(0) = acclzs(naccl) - zlatperi
             acclze(0) = acclze(naccl) - zlatperi
             acclez(0) = acclez(naccl)
             acclxw(0) = acclxw(naccl)
             acclsw(0) = acclsw(naccl)
             acclts(0) = acclts(naccl)
             accldt(0) = accldt(naccl)
             do l = 0,ntaccl
               acclet(l,0) = acclet(l,naccl)
             enddo
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      accls = .false.
      if (ntaccl > 0) accls = .true.
      do ia = 0, naccl
         if (acclez(ia).ne.0) accls = .true.
      enddo

      return
      end
c============================================================================ 
      subroutine countemlts()
      use Lattice

c     --- count electric multipole data sets 

      integer(ISZ):: ie,nnemlt

      nnemlt = -1
      do ie = 0, nemlt
         if (emltzs(ie) .ne. emltze(ie) .or. emltid(ie) > 0) nnemlt = ie
      enddo
      nemlt = nnemlt

      return
      end
c============================================================================ 
      subroutine setupemlts()
      use Lattice
      use Mult_data

      integer(ISZ):: ie,i,iz
      real(kind=8):: dzio2,emax,ephmax

      if (zlatperi > 0.) then
        if (nemlt >= 0) then
c       --- Note that these elemente are different since only the centers are
c       --- specified.
          if (emltid(0) == 0) then
             emltzs(0) = emltzs(nemlt) - zlatperi
             emltze(0) = emltze(nemlt) - zlatperi
             emltph(0) = emltph(nemlt)
             emltsf(0) = emltsf(nemlt)
             emltsc(0) = emltsc(nemlt)
             emltid(0) = emltid(nemlt)
             emltap(0) = emltap(nemlt) 
             emltrr(0) = emltrr(nemlt)
             emltrl(0) = emltrl(nemlt)
             emltgl(0) = emltgl(nemlt)
             emltgp(0) = emltgp(nemlt)
             emltpw(0) = emltpw(nemlt)
             emltpa(0) = emltpa(nemlt)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.  emlts = .false.
      do ie = 0, nemlt
         if (emltid(ie) > nemltsets) then
           call kaboom("setupemlts: emltid must be between 1 and nemltsets, or <= 0")
           return
         endif
         if (emltid(ie) > 0) emlts = .true.
      enddo

c     --- Precalculate the axial derivatives of the multipole moments
c     --- if they were not supplied by the user.
      do i=1,nemltsets
        dzio2 = .5/dzemlt(i)
        do ie=1,nesmult
c         --- Axial derivative of the moment strengths
c         --- Find the maximum value.
          emax = 0.
          do iz=0,nzemltmax
            if (abs(esemltp(iz,ie,i)) > emax) emax = abs(esemltp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (emax == 0.) then
            esemltp(0,ie,i) = 0.
            do iz=1,nzemltmax-1
              esemltp(iz,ie,i) = (esemlt(iz+1,ie,i) - esemlt(iz-1,ie,i))*dzio2
            enddo
            esemltp(nzemltmax,ie,i) = 0.
          endif
c         --- Axial derivative of the phase angle
c         --- Find the maximum value.
          ephmax = 0.
          do iz=0,nzemltmax
            if (abs(esemltphp(iz,ie,i)) > ephmax)
     &        ephmax = abs(esemltphp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (ephmax == 0.) then
            esemltphp(0,ie,i) = 0.
            do iz=1,nzemltmax-1
              esemltphp(iz,ie,i)=(esemltph(iz+1,ie,i)-esemltph(iz-1,ie,i))*dzio2
            enddo
            esemltphp(nzemltmax,ie,i) = 0.
          endif
        enddo
      enddo

      return
      end
c============================================================================ 
      subroutine countmmlts()
      use Lattice

c     --- count magnetic multipole data sets 

      integer(ISZ):: im,nnmmlt

      nnmmlt = -1
      do im = 0, nmmlt
         if (mmltzs(im) .ne. mmltze(im) .or. mmltid(im) > 0) nnmmlt = im
      enddo
      nmmlt = nnmmlt

      return 
      end 
c============================================================================ 
      subroutine setupmmlts()
      use Lattice
      use Mult_data

      integer(ISZ):: im,i,iz
      real(kind=8):: dzio2,mmax,mphmax

      if (zlatperi > 0.) then
        if (nmmlt >= 0) then
          if (mmltid(0) == 0) then
             mmltzs(0) = mmltzs(nmmlt) - zlatperi
             mmltze(0) = mmltze(nmmlt) - zlatperi
             mmltph(0) = mmltph(nmmlt)
             mmltsf(0) = mmltsf(nmmlt)
             mmltsc(0) = mmltsc(nmmlt)
             mmltid(0) = mmltid(nmmlt)
             mmltap(0) = mmltap(nmmlt) 
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      mmlts = .false.
      do im = 0, nmmlt
         if (mmltid(im) > nmmltsets) then
           call kaboom("setupmmlts: mmltid must be between 1 and nmmltsets, or <= 0")
           return
         endif
         if (mmltid(im) > 0) mmlts = .true.
      enddo

c     --- Precalculate the axial derivatives of the multipole moments
c     --- if they were not supplied by the user.
      do i=1,nmmltsets
        dzio2 = .5/dzmmlt(i)
        do im=1,nmsmult
c         --- Axial derivative of the moment strengths
c         --- Find the maximum value.
          mmax = 0.
          do iz=0,nzmmltmax
            if (abs(msmmltp(iz,im,i)) > mmax) mmax = abs(msmmltp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mmax == 0.) then
            msmmltp(0,im,i) = 0.
            do iz=1,nzmmltmax-1
              msmmltp(iz,im,i) = (msmmlt(iz+1,im,i) - msmmlt(iz-1,im,i))*dzio2
            enddo
            msmmltp(nzmmltmax,im,i) = 0.
          endif
c         --- Axial derivative of the phase angle
c         --- Find the maximum value.
          mphmax = 0.
          do iz=0,nzmmltmax
            if (abs(msmmltphp(iz,im,i)) > mphmax)
     &        mphmax = abs(msmmltphp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mphmax == 0.) then
            msmmltphp(0,im,i) = 0.
            do iz=1,nzmmltmax-1
              msmmltphp(iz,im,i)=(msmmltph(iz+1,im,i)-msmmltph(iz-1,im,i))*dzio2
            enddo
            msmmltphp(nzmmltmax,im,i) = 0.
          endif
        enddo
      enddo

      return
      end
c============================================================================ 
      subroutine countegrds()
      use Lattice

c     --- count egridded field data sets s

      integer(ISZ):: ie,nnegrd

      nnegrd = -1
      do ie = 0, negrd
         if (egrdzs(ie) .ne. egrdze(ie) .or. egrdid(ie) > 0) nnegrd = ie
      enddo
      negrd = nnegrd

      return 
      end 
c============================================================================ 
      subroutine setupegrds()
      use Lattice
      use EGRDdata

      integer(ISZ):: im,i

      if (zlatperi > 0.) then
        if (negrd >= 0) then
          if (egrdid(0) == 0) then
             egrdzs(0) = egrdzs(negrd) - zlatperi
             egrdze(0) = egrdze(negrd) - zlatperi
             egrdxs(0) = egrdxs(negrd)
             egrdys(0) = egrdys(negrd)
             egrdid(0) = egrdid(negrd)
             egrdsf(0) = egrdsf(negrd)
             egrdsc(0) = egrdsc(negrd)
             egrdsy(0) = egrdsy(negrd)
             egrdox(0) = egrdox(negrd)
             egrdoy(0) = egrdoy(negrd)
             egrdph(0) = egrdph(negrd) 
             egrdap(0) = egrdap(negrd)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      egrds = .false.
      do im = 0, negrd
         if (egrdid(im) > egrdns) then
           call kaboom("setupegrds: egrdid must be between 1 and egrdns, or <= 0")
           return
         endif
         if (egrdid(im) > 0) egrds = .true.
      enddo

c     --- Calculate 1/egrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if egrdns=0.
      do i=1,egrdns
        if (egrddxi(i) == 0.) egrddxi(i) = 1./egrddx(i)
        if (egrddyi(i) == 0.) egrddyi(i) = 1./egrddy(i)
        if (egrddzi(i) == 0.) egrddzi(i) = 1./egrddz(i)
      enddo

c     --- Calculate sines and cosines of egrdph for efficiency
      do i=0,negrd
        egrdsp(i) = sin(egrdph(i)) 
        egrdcp(i) = cos(egrdph(i)) 
      enddo 

      return
      end
c============================================================================ 
      subroutine countbgrds()
      use Lattice

c     --- count bgridded field data sets s

      integer(ISZ):: ib,nnbgrd

      nnbgrd = -1
      do ib = 0, nbgrd
         if (bgrdzs(ib) .ne. bgrdze(ib) .or. bgrdid(ib) > 0) nnbgrd = ib
      enddo
      nbgrd = nnbgrd

      return 
      end 
c============================================================================ 
      subroutine setupbgrds()
      use Lattice
      use BGRDdata

      integer(ISZ):: im,i

      if (zlatperi > 0.) then
        if (nbgrd >= 0) then
          if (bgrdid(0) == 0) then
             bgrdzs(0) = bgrdzs(nbgrd) - zlatperi
             bgrdze(0) = bgrdze(nbgrd) - zlatperi
             bgrdxs(0) = bgrdxs(nbgrd)
             bgrdys(0) = bgrdys(nbgrd)
             bgrdid(0) = bgrdid(nbgrd)
             bgrdsf(0) = bgrdsf(nbgrd)
             bgrdsc(0) = bgrdsc(nbgrd)
             bgrdsy(0) = bgrdsy(nbgrd)
             bgrdox(0) = bgrdox(nbgrd)
             bgrdoy(0) = bgrdoy(nbgrd)
             bgrdph(0) = bgrdph(nbgrd) 
             bgrdot(0) = bgrdot(nbgrd)
             bgrdop(0) = bgrdop(nbgrd)
             bgrdap(0) = bgrdap(nbgrd)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      bgrds = .false.
      do im = 0, nbgrd
         if (bgrdid(im) > bgrdns) then
           call kaboom("setupbgrds: bgrdid must be between 1 and bgrdns, or <= 0")
           return
         endif
         if (bgrdid(im) > 0) bgrds = .true.
      enddo

c     --- Calculate 1/bgrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if bgrdns=0.
      do i=1,bgrdns
        if (bgrddxi(i) == 0.) bgrddxi(i) = 1./bgrddx(i)
        if (bgrddyi(i) == 0. .and. bgrddy(i) .ne. 0.) bgrddyi(i) = 1./bgrddy(i)
        if (bgrddzi(i) == 0.) bgrddzi(i) = 1./bgrddz(i)
      enddo

      return
      end
c============================================================================ 
      subroutine countpgrds()
      use Lattice

      integer(ISZ):: ip,nnpgrd

c     --- count potential gridded field data sets 

      nnpgrd = -1
      do ip = 0, npgrd
         if (pgrdzs(ip) .ne. pgrdze(ip) .or. pgrdid(ip) > 0) nnpgrd = ip
      enddo
      npgrd = nnpgrd

      return
      end
c============================================================================ 
      subroutine setuppgrds()
      use Lattice
      use PGRDdata

      integer(ISZ):: im,i

      if (zlatperi > 0.) then
        if (npgrd >= 0) then
          if (pgrdid(0) == 0) then
             pgrdzs(0) = pgrdzs(npgrd) - zlatperi
             pgrdze(0) = pgrdze(npgrd) - zlatperi
             pgrdxs(0) = pgrdxs(npgrd)
             pgrdys(0) = pgrdys(npgrd)
             pgrdid(0) = pgrdid(npgrd)
             pgrdsf(0) = pgrdsf(npgrd)
             pgrdsc(0) = pgrdsc(npgrd)
             pgrdox(0) = pgrdox(npgrd)
             pgrdoy(0) = pgrdoy(npgrd)
             pgrdph(0) = pgrdph(npgrd) 
             pgrdap(0) = pgrdap(npgrd)
             pgrdrr(0) = pgrdrr(npgrd)
             pgrdrl(0) = pgrdrl(npgrd)
             pgrdgl(0) = pgrdgl(npgrd)
             pgrdgp(0) = pgrdgp(npgrd)
             pgrdpw(0) = pgrdpw(npgrd)
             pgrdpa(0) = pgrdpa(npgrd)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      pgrds = .false.
      do im = 0, npgrd
         if (pgrdid(im) > pgrdns) then
           call kaboom("setuppgrds: pgrdid must be between 1 and pgrdns, or <= 0")
           return
         endif
         if (pgrdid(im) > 0) pgrds = .true.
      enddo

c     --- Calculate 1/pgrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if pgrdns=0.
      do i=1,pgrdns
        if (pgrddxi(i) == 0.) pgrddxi(i) = 1./pgrddx(i)
        if (pgrddyi(i) == 0.) pgrddyi(i) = 1./pgrddy(i)
        if (pgrddzi(i) == 0.) pgrddzi(i) = 1./pgrddz(i)
      enddo

c     --- Calculate sines and cosines of pgrdph for efficiency
      do i=0,npgrd
        pgrdsp(i) = sin(pgrdph(i)) 
        pgrdcp(i) = cos(pgrdph(i)) 
      enddo 

      return
      end
c============================================================================ 
      subroutine countbsqgrads()
      use Lattice

      integer(ISZ):: ib,nnbsqgrad

c     --- count B squared gridded field data sets 

      nnbsqgrad = -1
      do ib = 0, nbsqgrad
         if (bsqgradzs(ib) .ne. bsqgradze(ib) .or. bsqgradid(ib) > 0) then
           nnbsqgrad = ib
         endif
      enddo
      nbsqgrad = nnbsqgrad

      return 
      end 
c============================================================================ 
      subroutine setupbsqgrads()
      use Lattice
      use BSQGRADdata

      integer(ISZ):: ib,id,i
      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,temp

      if (zlatperi > 0.) then
        if (nbsqgrad >= 0) then
          if (bsqgradid(0) == 0) then
             bsqgradzs(0) = bsqgradzs(nbsqgrad) - zlatperi
             bsqgradze(0) = bsqgradze(nbsqgrad) - zlatperi
             bsqgradxs(0) = bsqgradxs(nbsqgrad)
             bsqgradys(0) = bsqgradys(nbsqgrad)
             bsqgradid(0) = bsqgradid(nbsqgrad)
             bsqgradsf(0) = bsqgradsf(nbsqgrad)
             bsqgradsc(0) = bsqgradsc(nbsqgrad)
             bsqgradsy(0) = bsqgradsy(nbsqgrad)
             bsqgradox(0) = bsqgradox(nbsqgrad)
             bsqgradoy(0) = bsqgradoy(nbsqgrad)
             bsqgradph(0) = bsqgradph(nbsqgrad) 
             bsqgradap(0) = bsqgradap(nbsqgrad)
          endif
        endif
      endif

c     --- Set flags for existence of the elements - checking if the fields
c     --- are specified.
      bsqgrads = .false.
      do ib = 0, nbsqgrad
         if (bsqgradid(ib) > bsqgradns) then
           call kaboom("setupbsqgrads: bsqgradid must be between 1 and bsqgradns, or <= 0")
           return
         endif
         if (bsqgradid(ib) > 0) bsqgrads = .true.
      enddo

c     --- Calculate 1/bsqgraddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if bsqgradns=0.
      do i=1,bsqgradns
        if (bsqgraddxi(i) == 0.) bsqgraddxi(i) = 1./bsqgraddx(i)
        if (bsqgraddyi(i) == 0.) bsqgraddyi(i) = 1./bsqgraddy(i)
        if (bsqgraddzi(i) == 0.) bsqgraddzi(i) = 1./bsqgraddz(i)
      enddo

c     --- Calculate sines and cosines of bsqgradph for efficiency
      do i=0,nbsqgrad
        bsqgradsp(i) = sin(bsqgradph(i)) 
        bsqgradcp(i) = cos(bsqgradph(i)) 
      enddo 

      return
      end
c============================================================================ 
      subroutine calculatebsqgrad()
      use Subtimerstop
      use Lattice
      use BSQGRADdata
      use InGen

      integer(ISZ):: nn,ib,id
      integer(ISZ):: ix,iy,iz,nx,ny,nz
      real(kind=8):: xx,yy,zz,temp,dxi,dyi,dzi
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()
      
c     --- If bsqgrad is calculated, then setup the position
c     --- arrays. XXX This assumes that there is a one to one correspondence
c     --- betweeen bsqgrad elements and datasets, i.e. that the same id is
c     --- not used for multiple element instances.
c     --- The data is put into the bsqgradtemp array to avoid the
c     --- proliforation of multiple arrays.
c     --- The 11 is for 3 position, vz, gaminv, 3 B and 3 E.
c     --- The E is just wasted space since some of the field gathering
c     --- routines also gather E.
      if (bsqgradntemp == 0) then
        bsqgradntemp = 11
        call gchange("BSQGRADdata",0)

        do ib = 0, nbsqgrad
          id = bsqgradid(ib)
          do iz=0,bsqgradnz
            do iy=0,bsqgradny
              do ix=0,bsqgradnx

                xx = bsqgradxs(ib) + bsqgraddx(id)*ix
                yy = bsqgradys(ib) + bsqgraddy(id)*iy
                zz = bsqgradzs(ib) + bsqgraddz(id)*iz

c               --- Transverse rotation to take into account for a rotation 
c               --- of the field element.
                if ( bsqgradph(ib) .ne. 0. ) then
                  temp = xx  
                  xx =  temp*bsqgradcp(ib) - yy*bsqgradsp(ib) 
                  yy = +temp*bsqgradsp(ib) + yy*bsqgradcp(ib) 
                endif 

c               --- Add in any offsets
                xx = xx - bsqgradox(ib)
                yy = yy - bsqgradoy(ib)

                bsqgradtemp(ix,iy,iz,0,id) = xx
                bsqgradtemp(ix,iy,iz,1,id) = yy
                bsqgradtemp(ix,iy,iz,2,id) = zz
                bsqgradtemp(ix,iy,iz,3,id) = 0.
                bsqgradtemp(ix,iy,iz,4,id) = 1.

              enddo
            enddo
          enddo
        enddo
      endif

c     --- Now fetch the B fields
      nz = bsqgradnz
      nx = bsqgradnx
      ny = bsqgradny
      nn = (1+bsqgradnx)*(1+bsqgradny)*(1+bsqgradnz)
      do ib = 0, nbsqgrad
        id = bsqgradid(ib)

c       bsqgradtemp(0:nx,0:ny,0:nz, 5,id) = 0.
c       bsqgradtemp(0:nx,0:ny,0:nz, 6,id) = 0.
c       bsqgradtemp(0:nx,0:ny,0:nz, 7,id) = 0.
        bsqgradtemp(0:nx,0:ny,0:nz, 8,id) = bx0
        bsqgradtemp(0:nx,0:ny,0:nz, 9,id) = by0
        bsqgradtemp(0:nx,0:ny,0:nz,10,id) = bz0

c       --- handle uniform fields
c       call applyuniformfields(nn,bsqgradtemp(0:nx,0:ny,0:nz, 7,id),
c    &                               bsqgradtemp(0:nx,0:ny,0:nz,10,id))

c       --- Note that the uz and the time step sizes are all zero. This
c       --- has the effect for the residence corrected elements to find
c       --- the field at the given location without any correction.

c       --- handle quads
        call applyquad(nn,bsqgradtemp(0,0,0, 0,id),
     &                 bsqgradtemp(0,0,0, 1,id),nn,
     &                 bsqgradtemp(0,0,0, 2,id),
     &                 bsqgradtemp(0,0,0, 3,id),
     &                 bsqgradtemp(0,0,0, 4,id),
     &                 0.,0.,0.,.false.,
     &                 bsqgradtemp(0,0,0, 5,id),
     &                 bsqgradtemp(0,0,0, 6,id),
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id))

c       --- handle dipos 
        call applydipo(nn,nn,bsqgradtemp(0,0,0, 2,id),
     &                 bsqgradtemp(0,0,0, 3,id),
     &                 bsqgradtemp(0,0,0, 4,id),
     &                 0.,0.,0.,.false.,bsqgradtemp(0,0,0, 5,id),
     &                 bsqgradtemp(0,0,0, 6,id),
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id))

c       --- handle sexts
        call applysext(nn,bsqgradtemp(0,0,0, 0,id),
     &                 bsqgradtemp(0,0,0, 1,id),nn,
     &                 bsqgradtemp(0,0,0, 2,id),
     &                 bsqgradtemp(0,0,0, 3,id),
     &                 bsqgradtemp(0,0,0, 4,id),0.,0.,0.,.false.,
     &                 bsqgradtemp(0,0,0, 5,id),
     &                 bsqgradtemp(0,0,0, 6,id),
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id))

c       --- handle hard-edge electric and magnetic multipoles
        call applyhele(nn,bsqgradtemp(0,0,0, 0,id),
     &                 bsqgradtemp(0,0,0, 1,id),nn,
     &                 bsqgradtemp(0,0,0, 2,id),
     &                 bsqgradtemp(0,0,0, 3,id),
     &                 bsqgradtemp(0,0,0, 4,id),0.,0.,0.,.false.,
     &                 bsqgradtemp(0,0,0, 5,id),
     &                 bsqgradtemp(0,0,0, 6,id),
     &                 bsqgradtemp(0,0,0, 7,id),
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id),
     &                 bsqgradtemp(0,0,0,10,id))

c       --- handle magnetostatic multipole components
        call applymmlt(nn,bsqgradtemp(0,0,0, 0,id),
     &                 bsqgradtemp(0,0,0, 1,id),nn,
     &                 bsqgradtemp(0,0,0, 2,id),0.,0.,0.,.false.,
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id),
     &                 bsqgradtemp(0,0,0,10,id))

c       --- handle magnetic fields from 3-D grid
        call applybgrd(nn,bsqgradtemp(0,0,0, 0,id),
     &                 bsqgradtemp(0,0,0, 1,id),nn,
     &                 bsqgradtemp(0,0,0, 2,id),.false.,
     &                 bsqgradtemp(0,0,0, 8,id),
     &                 bsqgradtemp(0,0,0, 9,id),
     &                 bsqgradtemp(0,0,0,10,id))


c       --- Calculate B dot B
        do iz=0,bsqgradnz
          do iy=0,bsqgradny
            do ix=0,bsqgradnx
              bsqgradtemp(ix,iy,iz,5,id) =
     &                        bsqgradtemp(ix,iy,iz,8,id)**2 +
     &                        bsqgradtemp(ix,iy,iz,9,id)**2 +
     &                        bsqgradtemp(ix,iy,iz,10,id)**2
            enddo
          enddo
        enddo

c       --- Calculate the gradient
c       --- First, do the interior points
        dxi = bsqgraddxi(id)
        dyi = bsqgraddyi(id)
        dzi = bsqgraddzi(id)
        bsqgrad(1,:,:,1:nz-1,id) = (bsqgradtemp(:,:,2:nz,5,id) -
     &                              bsqgradtemp(:,:,0:nz-2,5,id))*0.5*dzi
        if (bsqgradnc == 3) then
          bsqgrad(2,1:nx-1,:,:,id) = (bsqgradtemp(2:nx,:,:,5,id) -
     &                                bsqgradtemp(0:nx-2,:,:,5,id))*0.5*dxi
          bsqgrad(3,:,1:ny-1,:,id) = (bsqgradtemp(:,2:ny,:,5,id) -
     &                                bsqgradtemp(:,0:ny-2,:,5,id))*0.5*dyi
        endif

c       --- Now do the boundaries
c       --- First upper boundaries, one-sided differences
        bsqgrad(1,:,:,0,id) = (bsqgradtemp(:,:,1,5,id) -
     &                         bsqgradtemp(:,:,0,5,id))*dzi
        bsqgrad(1,:,:,nz,id) = (bsqgradtemp(:,:,nz,5,id) -
     &                          bsqgradtemp(:,:,nz-1,5,id))*dzi
        if (bsqgradnc == 3) then
          bsqgrad(2,nx,:,:,id) = (bsqgradtemp(nx,:,:,5,id) -
     &                            bsqgradtemp(nx-1,:,:,5,id))*dxi
          bsqgrad(3,:,ny,:,id) = (bsqgradtemp(:,ny,:,5,id) -
     &                            bsqgradtemp(:,ny-1,:,5,id))*dyi
c         --- now lower boundaries. Treat differently depending on symmetry.
          if (bsqgradsy(ib) == 0) then
c           --- no symmetry
            bsqgrad(2,0,:,:,id) = (bsqgradtemp(1,:,:,5,id) -
     &                             bsqgradtemp(0,:,:,5,id))*dxi
            bsqgrad(3,:,0,:,id) = (bsqgradtemp(:,1,:,5,id) -
     &                             bsqgradtemp(:,0,:,5,id))*dyi
          else if (bsqgradsy(ib) == 2) then
c           --- quadrupole symmetry
            bsqgrad(2,0,:,:,id) = 0.
            bsqgrad(3,:,0,:,id) = 0.
          else
            call kaboom("bsqgrad symmetry value not supported")
            return
          endif
        endif
      enddo

!$OMP MASTER
      if (ltoptimesubs) timecalculatebsqgrad = timecalculatebsqgrad +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c============================================================================ 
c============================================================================ 
c============================================================================ 
      subroutine resetlatdrft()
      use Lattice
      use LatticeInternal

      call countdrfts
      call gchange("Lattice", 0)
      call setupdrfts

      ndrftol = 0
      if (ndrft >= 0)
     &  call getelemoverlaps(ndrft,drftzs,drftze,drftol,zlatbuffer,ndrftol)
      call gchange("LatticeInternal",0)

      if (ndrft >= 0)
     &  call getelemoverlapindices(ndrft,drftol,ndrftol,odrftoi,odrftio,odrftii,
     &                             odrftnn)

      return
      end
c============================================================================ 
      subroutine resetlatbend()
      use Lattice
      use LatticeInternal

      call countbends
      call gchange("Lattice", 0)
      call setupbends

      nbendol = 0
      if (nbend >= 0)
     &  call getelemoverlaps(nbend,bendzs,bendze,bendol,zlatbuffer,nbendol)
      call gchange("LatticeInternal",0)
      if (nbend >= 0)
     &  call getelemoverlapindices(nbend,bendol,nbendol,obendoi,obendio,obendii,
     &                             obendnn)

      return
      end
c============================================================================ 
      subroutine resetlatdipo()
      use Lattice
      use LatticeInternal

      call countdipos
      call gchange("Lattice", 0)
      call setupdipos

      ndipool = 0
      if (ndipo >= 0)
     &  call getelemoverlaps(ndipo,dipozs,dipoze,dipool,zlatbuffer,ndipool)
      call gchange("LatticeInternal",0)
      if (ndipo >= 0)
     &  call getelemoverlapindices(ndipo,dipool,ndipool,odipooi,odipoio,odipoii,
     &                             odiponn)

      return
      end
c============================================================================ 
      subroutine resetlatquad()
      use Lattice
      use LatticeInternal

      call countquads
      call gchange("Lattice", 0)
      call setupquads

      nquadol = 0
      if (nquad >= 0)
     &  call getelemoverlaps(nquad,quadzs,quadze,quadol,zlatbuffer,nquadol)
      call gchange("LatticeInternal",0)
      if (nquad >= 0)
     &  call getelemoverlapindices(nquad,quadol,nquadol,oquadoi,oquadio,oquadii,
     &                             oquadnn)

      return
      end
c============================================================================ 
      subroutine resetlatsext()
      use Lattice
      use LatticeInternal

      call countsexts
      call gchange("Lattice", 0)
      call setupsexts

      nsextol = 0
      if (nsext >= 0)
     &  call getelemoverlaps(nsext,sextzs,sextze,sextol,zlatbuffer,nsextol)
      call gchange("LatticeInternal",0)
      if (nsext >= 0)
     &  call getelemoverlapindices(nsext,sextol,nsextol,osextoi,osextio,osextii,
     &                             osextnn)

      return
      end
c============================================================================ 
      subroutine resetlathele()
      use Lattice
      use LatticeInternal

      call countheles
      call gchange("Lattice", 0)
      call setupheles

      nheleol = 0
      if (nhele >= 0)
     &  call getelemoverlaps(nhele,helezs,heleze,heleol,zlatbuffer,nheleol)
      call gchange("LatticeInternal",0)
      if (nhele >= 0)
     &  call getelemoverlapindices(nhele,heleol,nheleol,oheleoi,oheleio,oheleii,
     &                             ohelenn)

      return
      end
c============================================================================ 
      subroutine resetlataccl()
      use Lattice
      use LatticeInternal

      call countaccls
      call gchange("Lattice", 0)
      call setupaccls

      nacclol = 0
      if (naccl >= 0)
     &  call getelemoverlaps(naccl,acclzs,acclze,acclol,zlatbuffer,nacclol)
      call gchange("LatticeInternal",0)
      if (naccl >= 0)
     &  call getelemoverlapindices(naccl,acclol,nacclol,oaccloi,oacclio,oacclii,
     &                             oacclnn)

      return
      end
c============================================================================ 
      subroutine resetlatemlt()
      use Lattice
      use LatticeInternal

      call countemlts
      call gchange("Lattice", 0)
      call setupemlts

      nemltol = 0
      if (nemlt >= 0)
     &  call getelemoverlaps(nemlt,emltzs,emltze,emltol,zlatbuffer,nemltol)
      call gchange("LatticeInternal",0)
      if (nemlt >= 0)
     &  call getelemoverlapindices(nemlt,emltol,nemltol,oemltoi,oemltio,oemltii,
     &                             oemltnn)

      return
      end
c============================================================================ 
      subroutine resetlatmmlt()
      use Lattice
      use LatticeInternal

      call countmmlts
      call gchange("Lattice", 0)
      call setupmmlts

      nmmltol = 0
      if (nmmlt >= 0)
     &  call getelemoverlaps(nmmlt,mmltzs,mmltze,mmltol,zlatbuffer,nmmltol)
      call gchange("LatticeInternal",0)
      if (nmmlt >= 0)
     &  call getelemoverlapindices(nmmlt,mmltol,nmmltol,ommltoi,ommltio,ommltii,
     &                             ommltnn)

      return
      end
c============================================================================ 
      subroutine resetlategrd()
      use Lattice
      use LatticeInternal

      call countegrds
      call gchange("Lattice", 0)
      call setupegrds

      negrdol = 0
      if (negrd >= 0)
     &  call getelemoverlaps(negrd,egrdzs,egrdze,egrdol,zlatbuffer,negrdol)
      call gchange("LatticeInternal",0)
      if (negrd >= 0)
     &  call getelemoverlapindices(negrd,egrdol,negrdol,oegrdoi,oegrdio,oegrdii,
     &                             oegrdnn)

      return
      end
c============================================================================ 
      subroutine resetlatbgrd()
      use Lattice
      use LatticeInternal

      call countbgrds
      call gchange("Lattice", 0)
      call setupbgrds

      nbgrdol = 0
      
      if (nbgrd < 0) return

      call bgrdgetfulllength()
      call getelemoverlaps(nbgrd,bgrdfs,bgrdfe,bgrdol,zlatbuffer,nbgrdol)
      call gchange("LatticeInternal",0)
      call getelemoverlapindices(nbgrd,bgrdol,nbgrdol,obgrdoi,obgrdio,obgrdii,
     &                           obgrdnn)

      return
      end
c============================================================================ 
      subroutine resetlatpgrd()
      use Lattice
      use LatticeInternal

      call countpgrds
      call gchange("Lattice", 0)
      call setuppgrds

      npgrdol = 0
      if (npgrd >= 0)
     &  call getelemoverlaps(npgrd,pgrdzs,pgrdze,pgrdol,zlatbuffer,npgrdol)
      call gchange("LatticeInternal",0)
      if (npgrd >= 0)
     &  call getelemoverlapindices(npgrd,pgrdol,npgrdol,opgrdoi,opgrdio,opgrdii,
     &                             opgrdnn)

      return
      end
c============================================================================ 
      subroutine resetlatbsqgrad()
      use Lattice
      use LatticeInternal

      call countbsqgrads
      call gchange("Lattice", 0)
      call setupbsqgrads

      nbsqgradol = 0
      if (nbsqgrad >= 0)
     &  call getelemoverlaps(nbsqgrad,bsqgradzs,bsqgradze,bsqgradol,
     &                       zlatbuffer,nbsqgradol)
      call gchange("LatticeInternal",0)
      if (nbsqgrad >= 0)
     &  call getelemoverlapindices(nbsqgrad,bsqgradol,nbsqgradol,obsqgradoi,
     &                             obsqgradio,obsqgradii,obsqgradnn)

      return
      end
c============================================================================ 
c============================================================================ 
      subroutine resetlat()
      use Lattice

c  Resizes the lattice arrays to the final lengths needed by the actual data.
c  Sets dipole field, lattice element 0, flags for bends, etc.
c  This is somewhat more efficient since the gchange routines are only
c  called after all of the element have been updated.

      call countdrfts
      call countbends
      call countdipos
      call countquads
      call countsexts
      call countheles
      call countaccls
      call countemlts
      call countmmlts
      call countegrds
      call countpgrds
      call countbsqgrads

c     --- resize lattice arrays (use corrected ndipo size if autoset dipoles, 
c     --- bends imply dipoles for autoset dipoles, see following code)
      if (diposet) ndipo = max(ndipo, nbend)    
      call gchange("Lattice", 0)

c     --- the rest of the lattice reset is carried out in a work routine 
c     --- to avoid possibilites of a bug associated with resizing 
c     --- multidimensional dynamics arrays (used for hard-edge multipole 
c     --- elements) and then using them in the same routine  
      call setupdrfts
      call setupbends
      call setupdipos
      call setupquads
      call setupsexts
      call setupheles
      call setupaccls
      call setupemlts
      call setupmmlts
      call setupegrds
      call setuppgrds
      call setupbsqgrads

c     --- Get the overlap level of the lattice elements
      call getoverlaps()

c     --- Note that bgrd elements must be handled specially.
      call resetlatbgrd()

c     --- resetlat no longer needs to be called, as long as no changes are
c     --- made to the lattice.
      lresetlat = .false.

      return 
      end 
c============================================================================ 
      subroutine getoverlaps()
      use Lattice
      use LatticeInternal

c For each of the elements types, get the level of overlaps
c zlatbuffer is the amount of space (on meters) around each element. It is
c used so that elements right next to each other but not strictly overlapping
c are considered to be overlapping.

c     --- Note that bgrd elements must be handled specially.

      ndrftol = 0
      nbendol = 0
      ndipool = 0
      nquadol = 0
      nsextol = 0
      nheleol = 0
      nacclol = 0
      nemltol = 0
      nmmltol = 0
      negrdol = 0
      npgrdol = 0
      nbsqgradol = 0

      if (ndrft >= 0)
     &  call getelemoverlaps(ndrft,drftzs,drftze,drftol,zlatbuffer,ndrftol)
      if (nbend >= 0)
     &  call getelemoverlaps(nbend,bendzs,bendze,bendol,zlatbuffer,nbendol)
      if (ndipo >= 0)
     &  call getelemoverlaps(ndipo,dipozs,dipoze,dipool,zlatbuffer,ndipool)
      if (nquad >= 0)
     &  call getelemoverlaps(nquad,quadzs,quadze,quadol,zlatbuffer,nquadol)
      if (nsext >= 0)
     &  call getelemoverlaps(nsext,sextzs,sextze,sextol,zlatbuffer,nsextol)
      if (nhele >= 0)
     &  call getelemoverlaps(nhele,helezs,heleze,heleol,zlatbuffer,nheleol)
      if (naccl >= 0)
     &  call getelemoverlaps(naccl,acclzs,acclze,acclol,zlatbuffer,nacclol)
      if (nemlt >= 0)
     &  call getelemoverlaps(nemlt,emltzs,emltze,emltol,zlatbuffer,nemltol)
      if (nmmlt >= 0)
     &  call getelemoverlaps(nmmlt,mmltzs,mmltze,mmltol,zlatbuffer,nmmltol)
      if (negrd >= 0)
     &  call getelemoverlaps(negrd,egrdzs,egrdze,egrdol,zlatbuffer,negrdol)
      if (npgrd >= 0)
     &  call getelemoverlaps(npgrd,pgrdzs,pgrdze,pgrdol,zlatbuffer,npgrdol)
      if (nbsqgrad >= 0)
     &  call getelemoverlaps(nbsqgrad,bsqgradzs,bsqgradze,bsqgradol,zlatbuffer,nbsqgradol)

c Now, setup indicing arrays in LatticeInternal
      call gchange("LatticeInternal",0)

      if (ndrft >= 0)
     &  call getelemoverlapindices(ndrft,drftol,ndrftol,odrftoi,odrftio,odrftii,
     &                             odrftnn)
      if (nbend >= 0)
     &  call getelemoverlapindices(nbend,bendol,nbendol,obendoi,obendio,obendii,
     &                             obendnn)
      if (ndipo >= 0)
     &  call getelemoverlapindices(ndipo,dipool,ndipool,odipooi,odipoio,odipoii,
     &                             odiponn)
      if (nquad >= 0)
     &  call getelemoverlapindices(nquad,quadol,nquadol,oquadoi,oquadio,oquadii,
     &                             oquadnn)
      if (nsext >= 0)
     &  call getelemoverlapindices(nsext,sextol,nsextol,osextoi,osextio,osextii,
     &                             osextnn)
      if (nhele >= 0)
     &  call getelemoverlapindices(nhele,heleol,nheleol,oheleoi,oheleio,oheleii,
     &                             ohelenn)
      if (naccl >= 0)
     &  call getelemoverlapindices(naccl,acclol,nacclol,oaccloi,oacclio,oacclii,
     &                             oacclnn)
      if (nemlt >= 0)
     &  call getelemoverlapindices(nemlt,emltol,nemltol,oemltoi,oemltio,oemltii,
     &                             oemltnn)
      if (nmmlt >= 0)
     &  call getelemoverlapindices(nmmlt,mmltol,nmmltol,ommltoi,ommltio,ommltii,
     &                             ommltnn)
      if (negrd >= 0)
     &  call getelemoverlapindices(negrd,egrdol,negrdol,oegrdoi,oegrdio,oegrdii,
     &                             oegrdnn)
      if (npgrd >= 0)
     &  call getelemoverlapindices(npgrd,pgrdol,npgrdol,opgrdoi,opgrdio,opgrdii,
     &                             opgrdnn)
      if (nbsqgrad >= 0)
     &  call getelemoverlapindices(nbsqgrad,bsqgradol,nbsqgradol,obsqgradoi,obsqgradio,obsqgradii,
     &                             obsqgradnn)

      return
      end
c============================================================================ 
c============================================================================ 
      subroutine getelemoverlaps(nelem,elemzs,elemze,elemol,zlatbuffer,nelemol)
      integer(ISZ):: nelem,nelemol
      real(kind=8):: elemzs(0:nelem),elemze(0:nelem)
      integer(ISZ):: elemol(0:nelem)
      real(kind=8):: zlatbuffer

c Loops over lattice elements and for elements of the same type that overlap,
c keep track of the level of overlaps.

      integer(ISZ):: ie,ii
      logical(ISZ):: overlaps(1+nelem)

c     --- The first element is always at the lowest overlap level.
c     --- Don't set it if it is flagged to ignore overlaps.
      if (elemol(0) .ne. -1) elemol(0) = 1

c     --- Loop over the rest of the elements...
      do ie=1,nelem

c       --- Skip element if it is flagged to ignore overlaps.
        if (elemol(ie) == -1) cycle
      
c       --- First assume that this element doesn't overlap any others.
        overlaps = .false.

c       --- Loop over the previous elements and check which ones overlap the
c       --- current one.
        do ii=0,ie-1
c         --- Skip element if it is flagged to ignore overlaps.
          if (elemol(ii) == -1) cycle
c         --- For each overlap, set the corresponding overlap level to true.
          if (elemze(ie)+zlatbuffer > elemzs(ii) .and.
     &        elemze(ii)+zlatbuffer > elemzs(ie)) then
            overlaps(elemol(ii)) = .true.
          endif
        enddo

c       --- Loop over the list of overlap levels, looking for the lowest
c       --- level which has no overlap
        ii = 1
        do while (overlaps(ii) .and. ii <= 1+nelem)
          ii = ii + 1
        enddo
        elemol(ie) = ii

      enddo

c     --- Get the number of overlap levels. The max(1,...) is needed in case
c     --- all of the elements had the flag set to ignore overlaps, then
c     --- maxval(elemol) would be -1.
      nelemol = max(1,maxval(elemol))

      return
      end
c============================================================================ 
      subroutine getelemoverlapindices(nelem,elemol,
     &                                 nelemol,oelemoi,oelemio,oelemii,oelemnn)
      use Lattice
      use LatticeInternal
      integer(ISZ),intent(in):: nelem,nelemol
      integer(ISZ),intent(in):: elemol(0:nelem)
      integer(ISZ),intent(out):: oelemoi(0:nelem),oelemio(0:nelem)
      integer(ISZ),intent(out):: oelemii(nelemol),oelemnn(nelemol)

c Generate indicing information which is used to make finding elements
c efficient.
c    nelem    the number of elements
c    elemol   the overlap level of each element
c    nelemol  the number of overlapping levels
c The output is
c    oelemoi   a list of the element indices group by overlap level
c    oelemio   revervse list of the element indices group by overlap level
c    oelemii   the starting index in that list for each overlap level
c    oelemnn   the number of elements in each overlap level

      integer(ISZ):: io,ie,ii

c     --- Loop over the number of levels of overlap. For each level, get
c     --- the element indices that are at that level and increment
c     --- the number of elements at that level.
c     --- The absolute value of elemol is used since it can be set to -1,
c     --- signifying that overlaps are to be ignored. In that case, the
c     --- element will by default be in overlap level 1.
      ii = 0
      do io=1,nelemol
        oelemnn(io) = 0
        oelemii(io) = ii
        do ie=0,nelem
          if (abs(elemol(ie)) == io) then
            oelemoi(ii) = ie
            oelemio(ie) = ii - oelemii(io) + 1
            oelemnn(io) = oelemnn(io) + 1
            ii = ii + 1
          endif
        enddo
      enddo

      return
      end
c============================================================================ 
      subroutine setlatt
      use Subtimerstop
      use InGen
      use Picglb
      use Lattice,Only: lresetlat
c  Sets the lattice data for the current beam location.
c  These ride with the beam, just as the self-field does, and
c  should be loaded at the same time.
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (lresetlat) call resetlat()
      call setlattzt(zbeam,time,fstype)

!$OMP MASTER
      if (ltoptimesubs) timesetlatt = timesetlatt + wtime() - substarttime
!$OMP END MASTER

      return
      end
c============================================================================ 
      subroutine setlattzt(zbeam,time,fstype)
      use Lattice
      use LatticeInternal
      real(kind=8):: zbeam,time
      integer(ISZ):: fstype

c  Sets the lattice data at zbeam and time into 1d arrays.


c  If the transition from one element to the next falls to left of cell 
c  center, associate the "right" element with the cell; here, idipo(j)=i+1.

c                                                            (dipole elements
c                                                             begin and end
c            |<- cell j-1 ->|<-- cell j -->|<- cell j+1 ->|   at midpoints of 
c                                                             dipo-drifts)
c        |                      |                                |
c        |<----- element i ---->|<-------- element i+1 --------->|
c        |                      |                                |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c | dipo-drift  |dipo i |  dipo-drift   |   dipo i+1    |   dipo-drift    |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c              /         \             /                 \              z --->
c          dipozs(i)  dipoze(i)    dipozs(i+1)         dipoze(i+1)     
c       cdipozs(j-1)  cdipoze(j-1) cdipozs(j)          cdipoze(j)
c                                  cdipozs(j+1)        cdipoze(j+1)

      integer(ISZ):: io,iz
      integer(ISZ):: iquadert,iquaderu,iquaderr,ii,imultert,imulteru,imulterr
      integer(ISZ):: id
      real(kind=8):: zz,offset,wi

c     --- Save the z location and time at which the internal lattice frame
c     --- is set.
      zlframe = zbeam
      zltime = time

c     --- Assume at the start that there are not any elements in the mesh.
      lindrft = .false.
      linbend = .false.
      lindipo = .false.
      linquad = .false.
      linsext = .false.
      linhele = .false.
      linemlt = .false.
      linmmlt = .false.
      linaccl = .false.
      linegrd = .false.
      linbgrd = .false.
      linpgrd = .false.
      linbsqgrad = .false.

      do iz = 0, nzl
c        --- compute z of center of cell (may be off lattice)
c        --- If nzl is zero, then find the elements at zbeam. Otherwise,
c        --- find the elements at each of the grid points.
         if (nzl > 0) then
           zz = iz*dzl + zlmin + zbeam + 0.5*dzl
         else
           zz = zbeam
         endif

c        --- find index of nearest drft element
         do io=1,ndrftol
           call getelemid(zz,offset,ndrft,drftzs,drftze,drftol,
     &                    odrftnn(io),odrftoi(odrftii(io)),odrftio,
     &                    cdrftid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cdrftzs(iz,io) = drftzs(cdrftid(iz,io)) + offset
           cdrftze(iz,io) = drftze(cdrftid(iz,io)) + offset
           if (cdrftzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cdrftze(iz,io))
     &         lindrft(io) = .true.
         enddo

c        --- find index of nearest bend element
c        --- Note that while some setup is done for overlapping bends, they
c        --- are actually not supported since the idea of overlapping bends
c        --- doesn't make much sense. The extra indexing is done so that the
c        --- getelemid routine can still be used.
         if (nbend >= 0) then
           call getelemid(zz,offset,nbend,bendzs,bendze,bendol,
     &                    obendnn(1),obendoi(obendii(1)),obendio,
     &                    cbendid(iz),1)
c          --- load lattice info into comoving 1d array
           cbendzs(iz) = bendzs(cbendid(iz)) + offset
           cbendze(iz) = bendze(cbendid(iz)) + offset
           cbendrc(iz) = bendrc(cbendid(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbendzs(iz) <= zlmax+zbeam+zlatbuffer .and.
     &                        zlmin+zbeam-zlatbuffer <= cbendze(iz))
     &       linbend = .true.
         endif

c        --- find index of nearest dipole element
c        --- The cdipo arrays calculated here are needed for residence
c        --- corrections (when dipos is true), and internal conductors.
         do io=1,ndipool
           call getelemid(zz,offset,ndipo,dipozs,dipoze,dipool,
     &                    odiponn(io),odipooi(odipoii(io)),odipoio,
     &                    cdipoid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cdipozs(iz,io) = dipozs(cdipoid(iz,io)) + offset
           cdipoze(iz,io) = dipoze(cdipoid(iz,io)) + offset
           cdipoby(iz,io) = dipoby(cdipoid(iz,io))
           cdipobx(iz,io) = dipobx(cdipoid(iz,io))
           cdipota(iz,io) = dipota(cdipoid(iz,io))
           cdipotb(iz,io) = dipotb(cdipoid(iz,io))
           cdipoex(iz,io) = dipoex(cdipoid(iz,io))
           cdipoey(iz,io) = dipoey(cdipoid(iz,io))
           if (cdipozs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cdipoze(iz,io))
     &       lindipo(io) = .true.
         enddo

c        --- find index of nearest quad element
c        --- The cquad arrays calculated here are needed for residence
c        --- corrections (when quads is true) and the internal conductors.
         do io=1,nquadol
           call getelemid(zz,offset,nquad,quadzs,quadze,quadol,
     &                    oquadnn(io),oquadoi(oquadii(io)),oquadio,
     &                    cquadid(iz,io),io)
           id = cquadid(iz,io)
c          --- load lattice info into comoving 1d array
           cquadzs(iz,io) = quadzs(id) + offset
           cquadze(iz,io) = quadze(id) + offset
           cquaddb(iz,io) = quaddb(id)
           cquadde(iz,io) = quadde(id)
           cquadvx(iz,io) = quadvx(id)
           cquadvy(iz,io) = quadvy(id)
           if (ntquad > 0 .and. 
     &         quadts(id)<=time .and. time<quadts(id)+ntquad*quaddt(id)) then
             ii = int((time - quadts(id))/quaddt(id))
             wi = (time - quadts(id))/quaddt(id) - ii
             cquadde(iz,io) = cquadde(iz,io) + quadet(ii  ,id)*(1. - wi) +
     &                                         quadet(ii+1,id)*      wi 
             cquaddb(iz,io) = cquaddb(iz,io) + quadbt(ii  ,id)*(1. - wi) +
     &                                         quadbt(ii+1,id)*      wi 
           endif
c          --- set index for offset quads.  NOTE that this uses "nqerr"
c          --- and not "nquad".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via nqerr=nquad, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING, especially when iqerr is nonzero.
           if (zlatperi > 0.) then
             iquadert = id + (int(10.+(zz-zlatstrt)/zlatperi)-10)*nquad
           else
             iquadert = id
           endif
           iquaderu = max(iquadert, 0)
           if (nqerr == 0) then
             iquaderr = 0
           else
             iquaderr = mod(iquaderu, nqerr) 
           endif
           cqoffx(iz,io) = qoffx(iquaderr)
           cqoffy(iz,io) = qoffy(iquaderr)
           if (cquadzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cquadze(iz,io))
     &       linquad(io) = .true.
         enddo

c        --- find index of nearest sext element
c        --- The csext arrays calculated here are needed for residence
c        --- corrections (when sexts is true).
         do io=1,nsextol
           call getelemid(zz,offset,nsext,sextzs,sextze,sextol,
     &                    osextnn(io),osextoi(osextii(io)),osextio,
     &                    csextid(iz,io),io)
c          --- load lattice info into comoving 1d array
           csextzs(iz,io) = sextzs(csextid(iz,io)) + offset
           csextze(iz,io) = sextze(csextid(iz,io)) + offset
           csextdb(iz,io) = sextdb(csextid(iz,io))
           csextde(iz,io) = sextde(csextid(iz,io))
           if (csextzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= csextze(iz,io))
     &       linsext(io) = .true.
         enddo

c        --- find index of nearest hard-edge multipole element
         do io=1,nheleol
           call getelemid(zz,offset,nhele,helezs,heleze,heleol,
     &                    ohelenn(io),oheleoi(oheleii(io)),oheleio,
     &                    cheleid(iz,io),io)
c          --- load lattice info into comoving 1d array
           chelezs(iz,io) = helezs(cheleid(iz,io)) + offset
           cheleze(iz,io) = heleze(cheleid(iz,io)) + offset
           if (chelezs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cheleze(iz,io))
     &       linhele(io) = .true.
         enddo

c        --- find index of nearest acceleration element
         do io=1,nacclol
           call getelemid(zz,offset,naccl,acclzs,acclze,acclol,
     &                    oacclnn(io),oaccloi(oacclii(io)),oacclio,
     &                    cacclid(iz,io),io)
           id = cacclid(iz,io)
           if (acclze(id) + offset > acclzstt) then
             cacclzs(iz,io) = acclzs(id) + offset
             cacclze(iz,io) = acclze(id) + offset
             cacclez(iz,io) = acclez(id)
             cacclxw(iz,io) = acclxw(id) 
             cacclsw(iz,io) = acclsw(id) 
             if (acclts(id) <= time .and.
     &                         time < acclts(id)+ntaccl*accldt(id)) then
               ii = int((time - acclts(id))/accldt(id))
               wi = (time - acclts(id))/accldt(id) - ii
               cacclez(iz,io) = acclez(id) + acclet(ii  ,id)*(1. - wi) +
     &                                       acclet(ii+1,id)*      wi 
             endif
             if (cacclzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cacclze(iz,io))
     &         linaccl(io) = .true.
           else
             cacclzs(iz,io) = 0.
             cacclze(iz,io) = 0.
             cacclez(iz,io) = 0.
             cacclxw(iz,io) = 0.
             cacclsw(iz,io) = 0.  
           endif
         enddo

c        --- find index of nearest emlt element
c        --- The cemlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         do io=1,nemltol
           call getelemid(zz,offset,nemlt,emltzs,emltze,emltol,
     &                    oemltnn(io),oemltoi(oemltii(io)),oemltio,
     &                    cemltid(iz,io),io)
c          --- load lattice info into comoving 1d array
           id = cemltid(iz,io)
           cemltzs(iz,io) = emltzs(id) + offset
           cemltze(iz,io) = emltze(id) + offset
           cemltph(iz,io) = emltph(id)
           cemltsf(iz,io) = emltsf(id)
           cemltsc(iz,io) = emltsc(id)
           cemltim(iz,io) = emltid(id)
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cemltzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cemltze(iz,io))
     &       linemlt(io) = .true.
c          --- set index for offset emlts.  NOTE that this uses "neerr"
c          --- and not "nemlt".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via neerr=nemlt, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING
           if (zlatperi > 0.) then
             imultert = id+(int(10.+(zz-zlatstrt)/zlatperi)-10)*nemlt
           else
             imultert = id
           endif
           imulteru = max(imultert, 0)
           if (neerr == 0) then
             imulterr = 0
           else
             imulterr = mod(imulteru, neerr)
           endif
           cemltox(iz,io) = emltox(imulterr)
           cemltoy(iz,io) = emltoy(imulterr)
         enddo

c        --- find index of nearest mmlt element
c        --- The cmmlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         do io=1,nmmltol
           call getelemid(zz,offset,nmmlt,mmltzs,mmltze,mmltol,
     &                    ommltnn(io),ommltoi(ommltii(io)),ommltio,
     &                    cmmltid(iz,io),io)
c          --- load lattice info into comoving 1d array
           id = cmmltid(iz,io)
           cmmltzs(iz,io) = mmltzs(id) + offset
           cmmltze(iz,io) = mmltze(id) + offset
           cmmltph(iz,io) = mmltph(id)
           cmmltsf(iz,io) = mmltsf(id)
           cmmltsc(iz,io) = mmltsc(id)
           cmmltim(iz,io) = mmltid(id)
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cmmltzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cmmltze(iz,io))
     &       linmmlt(io) = .true.
c          --- set index for offset mmlts.  NOTE that this uses "nmerr"
c          --- and not "nmmlt".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via nmerr=nmmlt, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING
           if (zlatperi > 0.) then
             imultert = id+(int(10.+(zz-zlatstrt)/zlatperi)-10)*nmmlt
           else
             imultert = id
           endif
           imulteru = max(imultert, 0)
           if (nmerr == 0) then
             imulterr = 0
           else
             imulterr = mod(imulteru, nmerr)
           endif
           cmmltox(iz,io) = mmltox(imulterr)
           cmmltoy(iz,io) = mmltoy(imulterr)
         enddo

c        --- find index of nearest egrd element
c        --- The cegrdid array calculated here is needed for the magnetic
c        --- elements described by E field data on a 3-D grid.
         do io=1,negrdol
           call getelemid(zz,offset,negrd,egrdzs,egrdze,egrdol,
     &                    oegrdnn(io),oegrdoi(oegrdii(io)),oegrdio,
     &                    cegrdid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cegrdzs(iz,io) = egrdzs(cegrdid(iz,io)) + offset
           cegrdze(iz,io) = egrdze(cegrdid(iz,io)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cegrdzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cegrdze(iz,io))
     &       linegrd(io) = .true.
         enddo

c        --- find index of nearest bgrd element
c        --- The cbgrdid array calculated here is needed for the magnetic
c        --- elements described by B field data on a 3-D grid.
c        --- Note that the full lengths are used here, in case there is
c        --- a rotation off the z-axis.
         do io=1,nbgrdol
           call getelemid(zz,offset,nbgrd,bgrdfs,bgrdfe,bgrdol,
     &                    obgrdnn(io),obgrdoi(obgrdii(io)),obgrdio,
     &                    cbgrdid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cbgrdzs(iz,io) = bgrdfs(cbgrdid(iz,io)) + offset
           cbgrdze(iz,io) = bgrdfe(cbgrdid(iz,io)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbgrdzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cbgrdze(iz,io))
     &       linbgrd(io) = .true.
         enddo

c        --- find index of nearest pgrd element
c        --- The cpgrdid array calculated here is needed for the electrostatic
c        --- elements described by potential data on a 3-D grid.
         do io=1,npgrdol
           call getelemid(zz,offset,npgrd,pgrdzs,pgrdze,pgrdol,
     &                    opgrdnn(io),opgrdoi(opgrdii(io)),opgrdio,
     &                    cpgrdid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cpgrdzs(iz,io) = pgrdzs(cpgrdid(iz,io)) + offset
           cpgrdze(iz,io) = pgrdze(cpgrdid(iz,io)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cpgrdzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cpgrdze(iz,io))
     &       linpgrd(io) = .true.
         enddo

c        --- find index of nearest bsqgrad element
c        --- The cbsqgradid array calculated here is needed for the magnetic
c        --- elements described by grad B^2 field data on a 3-D grid.
         do io=1,nbsqgradol
           call getelemid(zz,offset,nbsqgrad,bsqgradzs,bsqgradze,bsqgradol,
     &                    obsqgradnn(io),obsqgradoi(obsqgradii(io)),obsqgradio,
     &                    cbsqgradid(iz,io),io)
c          --- load lattice info into comoving 1d array
           cbsqgradzs(iz,io) = bsqgradzs(cbsqgradid(iz,io)) + offset
           cbsqgradze(iz,io) = bsqgradze(cbsqgradid(iz,io)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbsqgradzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                           zlmin+zbeam-zlatbuffer <= cbsqgradze(iz,io))
     &       linbsqgrad(io) = .true.
         enddo

      enddo

c     --- For each element type, check to see if any elements are in the mesh.
      lindrft(0)    = any(lindrft)
      lindipo(0)    = any(lindipo)
      linquad(0)    = any(linquad)
      linsext(0)    = any(linsext)
      linhele(0)    = any(linhele)
      linemlt(0)    = any(linemlt)
      linmmlt(0)    = any(linmmlt)
      linaccl(0)    = any(linaccl)
      linegrd(0)    = any(linegrd)
      linbgrd(0)    = any(linbgrd)
      linpgrd(0)    = any(linpgrd)
      linbsqgrad(0) = any(linbsqgrad)

c     --- All parallel slave processors must know if any processors
c     --- are in a bend.  There is probably a better way of doing this.
#ifdef MPIPARALLEL
      if (bends) then
        wi = 0.
        if (linbend) wi = 1.
        call parallelmaxrealarray(wi,1)
        if (wi > 0.5) linbend = .true.
      endif
#endif

      return
      end
c=============================================================================
      subroutine getelemid(z,offset,nelem,elemzs,elemze,elemol,
     &                     oelemnn,oelemoi,oelemio,id,io)
      use Lattice,Only: zlatstrt,zlatperi
      integer(ISZ):: nelem,oelemnn,id,io
      real(kind=8):: z,offset,elemzs(0:nelem),elemze(0:nelem)
      integer(ISZ):: elemol(0:nelem),oelemoi(oelemnn),oelemio(0:nelem)

c Given information about the starts and stops of lattice elements, this
c routine finds the number of the element closest to the specified z value.
c It assumes that some previous value of the element number is passed in
c through id. This is used as a starting point for the search.

c  If the transition from one element to the next falls to left of z,
c  associate the "right" element with z; here, ielem(j)=i+1.

c                                                            (elements
c                                                             begin and end
c            |<- cell j-1 ->|<-- cell j -->|<- cell j+1 ->|   at midpoints of 
c                                                             elem-drifts)
c        |                      |                                |
c        |<----- element i ---->|<-------- element i+1 --------->|
c        |                      |                                |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c | elem-drift  |elem i |  elem-drift   |   elem i+1    |   elem-drift    |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c              /         \             /                 \              z --->
c          elemzs(i)  elemze(i)    elemzs(i+1)         elemze(i+1)     
c       celemzs(j-1)  celemze(j-1) celemzs(j)          celemze(j)
c                                  celemzs(j+1)        celemze(j+1)

      integer(ISZ):: ii
      real(kind=8):: zc

c The Lattice is only needed for zlatstrt and zlatperi.

c --- Offset of the lattice, for periodic repeat.
      if (zlatperi > 0.) then
        offset = (int(10.+(z-zlatstrt)/zlatperi)-10)*zlatperi+zlatstrt
      else
        offset = zlatstrt
      endif

c --- Check if there is only one element. If so, return id is that element.
      if (oelemnn == 1) then
        id = oelemoi(1)
        return
      endif

c --- Make sure id actually refers to an element of the overlap level.
c --- The absolute value of elemol is used since it can be set to -1,
c --- signifying that overlaps are to be ignored. In that case, the
c --- element will by default be in overlap level 1.
      if (abs(elemol(id)) .ne. io) id = oelemoi(1)

c --- Check if z is somehow greater than the current element (id) location.
c --- If so, reset it to 0 and start from the beginning. This happens if
c --- during a periodic repeat of the lattice or if setlatt was called out
c --- of order. Also, make sure that id is within the proper bounds.
      if (id > oelemoi(oelemnn)) then
        id = oelemoi(1)
      else if (id > oelemoi(1)) then
        if (id > oelemoi(oelemnn)) then
          id = oelemoi(1)
        else if (z < .5*(elemze(oelemoi(oelemio(id)-1))+elemzs(id))+offset) then
          id = oelemoi(1)
        endif
      else if (id < oelemoi(1)) then
        id = oelemoi(1)
      endif
      id = oelemoi(1)

c --- Make ii consistent with id
      ii = oelemio(id)

c --- Scan sequentially through the list until the element is found.
c --- Note that loop is written in this way so that if ii reaches oelemnn
c --- there is not an out of bounds reference (oelemoi(ii+1)). The initial
c --- statement, zc = 0., is just to insure the zc has value.
      zc = 0.
      if (ii < oelemnn) zc = 0.5*(elemze(id) + elemzs(oelemoi(ii+1))) + offset
      do while (ii < oelemnn .and. z > zc)
        ii = ii + 1
        id = oelemoi(ii)
        if (ii < oelemnn) zc = 0.5*(elemze(id) + elemzs(oelemoi(ii+1))) + offset
      end do

      return
      end
c=============================================================================
      subroutine getdipoid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the dipole at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (ndipo < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,ndipo,dipozs,dipoze,dipool,
     &                 odiponn(io),odipooi(odipoii(io):),odipoio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,ndipo,dipozs,dipoze,dipool,
     &                 odiponn(1),odipooi(odipoii(1):),odipoio,id,1)
c       --- Loop over rest of overlap levels looking for a dipo with lower
c       --- z start.
        id_temp = id
        do ii=2,ndipool
          call getelemid(zz,offset_temp,ndipo,dipozs,dipoze,dipool,
     &                   odiponn(ii),odipooi(odipoii(ii):),odipoio,id_temp,ii)
          if (dipozs(id_temp)+offset_temp < dipozs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getquadid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the quadrupole at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (nquad < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nquad,quadzs,quadze,quadol,
     &                 oquadnn(io),oquadoi(oquadii(io):),oquadio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nquad,quadzs,quadze,quadol,
     &                 oquadnn(1),oquadoi(oquadii(1):),oquadio,id,1)
c       --- Loop over rest of overlap levels looking for a quad with lower
c       --- z start.
        id_temp = id
        do ii=2,nquadol
          call getelemid(zz,offset_temp,nquad,quadzs,quadze,quadol,
     &                   oquadnn(ii),oquadoi(oquadii(ii):),oquadio,id_temp,ii)
          if (quadzs(id_temp)+offset_temp < quadzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getheleid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the hele at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (nhele < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nhele,helezs,heleze,heleol,
     &                 ohelenn(io),oheleoi(oheleii(io):),oheleio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nhele,helezs,heleze,heleol,
     &                 ohelenn(1),oheleoi(oheleii(1):),oheleio,id,1)
c       --- Loop over rest of overlap levels looking for a hele with lower
c       --- z start.
        id_temp = id
        do ii=2,nheleol
          call getelemid(zz,offset_temp,nhele,helezs,heleze,heleol,
     &                   ohelenn(ii),oheleoi(oheleii(ii):),oheleio,id_temp,ii)
          if (helezs(id_temp)+offset_temp < helezs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getemltid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the emlt at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (nemlt < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nemlt,emltzs,emltze,emltol,
     &                 oemltnn(io),oemltoi(oemltii(io):),oemltio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nemlt,emltzs,emltze,emltol,
     &                 oemltnn(1),oemltoi(oemltii(1):),oemltio,id,1)
c       --- Loop over rest of overlap levels looking for a emlt with lower
c       --- z start.
        id_temp = id
        do ii=2,nemltol
          call getelemid(zz,offset_temp,nemlt,emltzs,emltze,emltol,
     &                   oemltnn(ii),oemltoi(oemltii(ii):),oemltio,id_temp,ii)
          if (emltzs(id_temp)+offset_temp < emltzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getmmltid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the mmlt at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (nmmlt < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nmmlt,mmltzs,mmltze,mmltol,
     &                 ommltnn(io),ommltoi(ommltii(io):),ommltio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nmmlt,mmltzs,mmltze,mmltol,
     &                 ommltnn(1),ommltoi(ommltii(1):),ommltio,id,1)
c       --- Loop over rest of overlap levels looking for a mmlt with lower
c       --- z start.
        id_temp = id
        do ii=2,nmmltol
          call getelemid(zz,offset_temp,nmmlt,mmltzs,mmltze,mmltol,
     &                   ommltnn(ii),ommltoi(ommltii(ii):),ommltio,id_temp,ii)
          if (mmltzs(id_temp)+offset_temp < mmltzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getacclid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the accl at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (naccl < 0) return

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,naccl,acclzs,acclze,acclol,
     &                 oacclnn(io),oaccloi(oacclii(io):),oacclio,id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,naccl,acclzs,acclze,acclol,
     &                 oacclnn(1),oaccloi(oacclii(1):),oacclio,id,1)
c       --- Loop over rest of overlap levels looking for a accl with lower
c       --- z start.
        id_temp = id
        do ii=2,nacclol
          call getelemid(zz,offset_temp,naccl,acclzs,acclze,acclol,
     &                   oacclnn(ii),oaccloi(oacclii(ii):),oacclio,id_temp,ii)
          if (acclzs(id_temp)+offset_temp < acclzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine applyuniformfields(np,ex,ey,ez,bx,by,bz)
      use InGen
      integer(ISZ):: np
      real(kind=8):: ex(np),ey(np),ez(np),bx(np),by(np),bz(np)

c Apply uniform fields.

      if (ex0/=0.) ex = ex + ex0
      if (ey0/=0.) ey = ey + ey0
      if (ez0/=0.) ez = ez + ez0
      if (bx0/=0.) bx = bx + bx0
      if (by0/=0.) by = by + by0
      if (bz0/=0.) bz = bz + bz0

      return
      end
c============================================================================ 
      subroutine getbend(np,npz,zp,uzp,gaminv,bendres,bendradi,dtl,dtr,lslice)
      use InGen
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: dtl,dtr
      real(kind=8):: zp(npz),uzp(npz),gaminv(npz)
      logical(ISZ):: lslice
      real(kind=8):: bendres(np),bendradi(np)

c  Gets residence factor and radius for bends.
c  For periodic runs, assumes mesh period length = lattice period length

      integer(ISZ):: ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8),allocatable:: vz(:),dzi(:)

      if (.not. bends .or. .not. linbend) then
        bendres = 0.
        return
      endif

      allocate(vz(npz),dzi(npz))

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

c  Extract the local bend properties for each particle off the general lattice

      if (lslice) then
c       --- find z-cell in which particle lies
        iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c       --- Precalculate these quantities. zl is the min of the two and
c       --- zr is the max. This generalizes the routine, allowing left
c       --- moving particles, vz < 0.
        z1 = zp(1) + vz(1)*dtl
        z2 = zp(1) + vz(1)*dtr
c       --- "left" end of velocity advance step
        zl = min(z1,z2)
        fl = 0.
        if (zl >= cbendzs(iz) .and. zl < cbendze(iz)) fl = 1.
c       --- "right" end of velocity advance step
        zr = max(z1,z2)
        fr = 0.
        if (zr >= cbendzs(iz) .and. zr < cbendze(iz)) fr = 1.
c       --- residence fraction
        frac = fl
        if (fl > fr) frac = (cbendze(iz)-zl)*dzi(1)
        if (fr > fl) frac = (zr-cbendzs(iz))*dzi(1)
        bendres(1) = frac
        bendradi(1) = cbendrc(iz)
        return
      endif

      do ip = 1, np

c       --- find z-cell in which particle lies
        iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c       --- Precalculate these quantities. zl is the min of the two and
c       --- zr is the max. This generalizes the routine, allowing left
c       --- moving particles, vz < 0.
        z1 = zp(ip) + vz(ip)*dtl
        z2 = zp(ip) + vz(ip)*dtr
c       --- "left" end of velocity advance step
        zl = min(z1,z2)
        fl = 0.
        if (zl >= cbendzs(iz) .and. zl < cbendze(iz)) fl = 1.
c       --- "right" end of velocity advance step
        zr = max(z1,z2)
        fr = 0.
        if (zr >= cbendzs(iz) .and. zr < cbendze(iz)) fr = 1.
c       --- residence fraction
        frac = fl
        if (fl > fr) frac = (cbendze(iz)-zl)*dzi(ip)
        if (fr > fl) frac = (zr-cbendzs(iz))*dzi(ip)

        bendres(ip) = frac
        bendradi(ip) = cbendrc(iz)
      enddo

      deallocate(vz,dzi)

      return
      end
c=============================================================================
      subroutine applybend(np,xp,uzp,npz,bendres,bendradi,m,q,lslice,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),uzp(np),bendres(npz),bendradi(npz)
      real(kind=8):: m,q
      logical(ISZ):: lslice
      real(kind=8):: by(np)

c Apply transformation for bending field.

      integer(ISZ):: ip
      real(kind=8):: qi,bres,brad

      if (.not. bends .or. .not. linbend) return

      if (lslice) then
        bres = bendres(1)
        brad = bendradi(1)
        if (bres == 0.) return
      endif
      qi = 1./q
      do ip = 1,np
        if (.not. lslice) then
          bres = bendres(ip)
          brad = bendradi(ip)
        endif
        if (bres > 0.) then
          by(ip) = by(ip) - bres*(m*qi)*uzp(ip)/(brad + xp(ip))
        endif
      enddo

      return
      end
c=============================================================================
      subroutine applyquad(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse focusing element quad.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   quadrupole electric field
c   bx, by   quadrupole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac,xpmqoff,ypmqoff,xph,yph
      real(kind=8):: cosph,sinph
      real(kind=8):: dodec
      real(kind=8), allocatable:: vz(:),dzi(:)
      integer(ISZ):: allocerror

      if (.not. quads .or. .not. linquad(0)) return
      allocate(vz(npz),dzi(npz),stat=allocerror)
      if (allocerror /= 0) then
        print*,"applyquad: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",npz
        call kaboom("applyquad: allocation error")
        return
      endif

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nquadol
        if (.not. linquad(io)) cycle

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= cquadzs(iz,io) .and. zl < cquadze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= cquadzs(iz,io) .and. zr < cquadze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cquadze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-cquadzs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= cquadzs(iz,io) .and. zl < cquadze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= cquadzs(iz,io) .and. zr < cquadze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cquadze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-cquadzs(iz,io))*dzi(ip)
          endif
      
c         --- set the field
          if (frac > 0.) then
            xpmqoff = xp(ip) - cqoffx(iz,io)
            ypmqoff = yp(ip) - cqoffy(iz,io)
            if (quadph(cquadid(iz,io)) .ne. 0.) then
              cosph = cos(quadph(cquadid(iz,io)))
              sinph = sin(quadph(cquadid(iz,io)))
              xph = +xpmqoff*cosph + ypmqoff*sinph
              yph = -xpmqoff*sinph + ypmqoff*cosph
              xpmqoff = xph
              ypmqoff = yph
            endif
            bx(ip) = bx(ip) + cquaddb(iz,io)*frac*ypmqoff
            by(ip) = by(ip) + cquaddb(iz,io)*frac*xpmqoff
            ex(ip) = ex(ip) + cquadde(iz,io)*frac*xpmqoff
            ey(ip) = ey(ip) - cquadde(iz,io)*frac*ypmqoff
            dodec = quaddo(cquadid(iz,io))
            if (dodec .ne. 0.) then
              ex(ip) = ex(ip) + dodec*cquadde(iz,io)*frac*
     &                 (6.*xpmqoff**5 - 60.*xpmqoff**3*ypmqoff**2 +
     &                 30.*xpmqoff*ypmqoff**4)
              ey(ip) = ey(ip) - dodec*cquadde(iz,io)*frac*
     &                 (6.*ypmqoff**5 - 60.*ypmqoff**3*xpmqoff**2 +
     &                 30.*ypmqoff*xpmqoff**4)
            endif
            if (quadph(cquadid(iz,io)) .ne. 0.) then
              xph = bx(ip)
              yph = by(ip)
              bx(ip) = +xph*cosph - yph*sinph
              by(ip) = +xph*sinph + yph*cosph
              xph = ex(ip)
              yph = ey(ip)
              ex(ip) = +xph*cosph - yph*sinph
              ey(ip) = +xph*sinph + yph*cosph
            endif
          endif
        enddo
      enddo
      deallocate(vz,dzi)

      return
      end
c=============================================================================
      subroutine applydipo(np,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse bending element dipo.
c Input:
c   np       number of particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   dipole electric field
c   bx, by   dipole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8), allocatable:: vz(:),dzi(:)
      integer(ISZ):: allocerror

      if (.not. dipos .or. .not. lindipo(0)) return
      allocate(vz(npz),dzi(npz),stat=allocerror)
      if (allocerror /= 0) then
        print*,"applydipo: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",npz
        call kaboom("applydipo: allocation error")
        return
      endif

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,ndipool
        if (.not. lindipo(io)) cycle

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= cdipozs(iz,io) .and. zl < cdipoze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= cdipozs(iz,io) .and. zr < cdipoze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cdipoze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-cdipozs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= cdipozs(iz,io) .and. zl < cdipoze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= cdipozs(iz,io) .and. zr < cdipoze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cdipoze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-cdipozs(iz,io))*dzi(ip)
          endif

c         --- set the field
          if (frac > 0.) then
            by(ip) = by(ip) + cdipoby(iz,io)*frac
            bx(ip) = bx(ip) + cdipobx(iz,io)*frac
            ex(ip) = ex(ip) + cdipoex(iz,io)*frac
            ey(ip) = ey(ip) + cdipoey(iz,io)*frac
          endif
        enddo
      enddo
      deallocate(vz,dzi)

      return
      end
c=========================================================================== 
      subroutine zbendcor(pgroup,np,ipmin,ddt,bendres,bendradi)
      use ParticleGroupmodule
      use Lattice,Only:bends
      use LatticeInternal,Only:linbend
      use Subtimerstop
      type(ParticleGroup):: pgroup
      integer(ISZ):: np,ipmin
      real(kind=8):: ddt
      real(kind=8):: bendres(ipmin:ipmin+np-1), bendradi(ipmin:ipmin+np-1)

c  Corrects particle axial positions for bending; zp is position
c  "along the centerline."

      integer(ISZ):: ip
      real(kind=8):: xprv,xc
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (.not. bends .or. .not. linbend) then
!$OMP MASTER
        if (ltoptimesubs) timezbendcor = timezbendcor + wtime() - substarttime
!$OMP END MASTER
        return
      endif

      do ip=ipmin,ipmin+np-1
        xprv = pgroup%xp(ip) - ddt*pgroup%uxp(ip)*pgroup%gaminv(ip)
        xc = 0.5*(pgroup%xp(ip) + xprv)
        pgroup%zp(ip) = pgroup%zp(ip) +
     &                  ddt*pgroup%uzp(ip)*pgroup%gaminv(ip)*bendres(ip)*
     &                  (bendradi(ip)/(bendradi(ip) + xc) - 1.)
      enddo

!$OMP MASTER
      if (ltoptimesubs) timezbendcor = timezbendcor + wtime() - substarttime
!$OMP END MASTER
      return
      end      
c=========================================================================== 
c=========================================================================== 
      subroutine sledgcor(pgroup,np,ipmin,zpo,zbeam,zbeamo,vbeam,m,q,lslice)
      use ParticleGroupmodule
      use Lattice
      use LatticeInternal
      type(ParticleGroup):: pgroup
      integer(ISZ):: np,ipmin
      real(kind=8):: zpo(ipmin:ipmin+np-1)
      real(kind=8):: zbeam,zbeamo,vbeam
      real(kind=8):: q, m
      logical(ISZ):: lslice

c Applies position and velocity jump corrections on entry/exit of dipoles 
c with slanted faces.  (Must be called after all z advancement is done this 
c timestep, so that entry/exit of the bend is detected exactly once.)    
c Accounts for the fact that the jumped velocities act on the particles for 
c less than the entire andvance step on entry/exit of the dipole.  

      integer(ISZ):: io,ip,iz
      real(kind=8),allocatable:: vzi(:)
      real(kind=8):: qoverm,factor

      if (.not. dipos .or. .not. lindipo(0)) return

      allocate(vzi(ipmin:ipmin+np-1))

      qoverm = q/m
      if (lslice) then
        vzi = 1./dvnz(vbeam)
      else
        vzi = pgroup%uzp(ipmin:ipmin+np-1)*pgroup%gaminv(ipmin:ipmin+np-1)
        vzi = 1./dvnz(vzi)
      endif

      do io=1,ndipool
        if (.not. lindipo(io)) cycle

        if (lslice) then
          iz = 0
          if (.not.((zbeamo < cdipozs(iz,io) .and. zbeam >= cdipozs(iz,io)).or.
     &              (zbeamo < cdipoze(iz,io) .and. zbeam >= cdipoze(iz,io))))
     &      cycle
        endif
        do ip = ipmin,ipmin+np-1
          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = min(nzl,int(max(0.,(pgroup%zp(ip)-zlmin-zlframe)*dzli+0.5)))
          endif
          if ((zpo(ip) < cdipozs(iz,io)) .and.
     &        (pgroup%zp(ip) >= cdipozs(iz,io))) then 
c           --- entered dipole
            factor = cdipoby(iz,io)*cdipota(iz,io)*qoverm
            if (factor == 0) cycle
            pgroup%uxp(ip) = pgroup%uxp(ip) + pgroup%xp(ip)*factor
            pgroup%uyp(ip) = pgroup%uyp(ip) - pgroup%yp(ip)*factor
            pgroup%xp(ip)  = pgroup%xp(ip) +
     &            (pgroup%zp(ip)-cdipozs(iz,io))*vzi(ip)*pgroup%xp(ip)*factor
            pgroup%yp(ip)  = pgroup%yp(ip) -
     &            (pgroup%zp(ip)-cdipozs(iz,io))*vzi(ip)*pgroup%yp(ip)*factor
          elseif ((zpo(ip) <  cdipoze(iz,io)) .and.
     &             (pgroup%zp(ip) >= cdipoze(iz,io))) then  
c           --- exited dipole
            factor = cdipoby(iz,io)*cdipotb(iz,io)*qoverm
            if (factor == 0) cycle
            pgroup%uxp(ip) = pgroup%uxp(ip) + pgroup%xp(ip)*factor
            pgroup%uyp(ip) = pgroup%uyp(ip) - pgroup%yp(ip)*factor
            pgroup%xp(ip)  = pgroup%xp(ip) +
     &            (pgroup%zp(ip)-cdipoze(iz,io))*vzi(ip)*pgroup%xp(ip)*factor
            pgroup%yp(ip)  = pgroup%yp(ip) -
     &            (pgroup%zp(ip)-cdipoze(iz,io))*vzi(ip)*pgroup%yp(ip)*factor
          endif 
        enddo
      enddo

      deallocate(vzi)

      return
      end
c=============================================================================
      subroutine applysext(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse element sext.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   sextupole electric field
c   bx, by   sextupole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: xpmqoff,ypmqoff
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8), allocatable:: vz(:),dzi(:)
      integer(ISZ):: allocerror

      if (.not. sexts .or. .not. linsext(0)) return
      allocate(vz(npz),dzi(npz),stat=allocerror)
      if (allocerror /= 0) then
        print*,"applysext: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",npz
        call kaboom("applysext: allocation error")
        return
      endif

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nsextol
        if (.not. linsext(io)) cycle

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= csextzs(iz,io) .and. zl < csextze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= csextzs(iz,io) .and. zr < csextze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (csextze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-csextzs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= csextzs(iz,io) .and. zl < csextze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= csextzs(iz,io) .and. zr < csextze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (csextze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-csextzs(iz,io))*dzi(ip)
          endif

c         --- set the field
          if (frac > 0.) then
            xpmqoff = xp(ip)
            ypmqoff = yp(ip)
            bx(ip) = bx(ip) + csextdb(iz,io)*frac*3.*(xpmqoff**2 - ypmqoff**2)
            by(ip) = by(ip) + csextdb(iz,io)*frac*(-6.)*xpmqoff*ypmqoff
            ex(ip) = ex(ip) + csextde(iz,io)*frac*3.*(xpmqoff**2 - ypmqoff**2)
            ey(ip) = ey(ip) + csextde(iz,io)*frac*(-6.)*xpmqoff*ypmqoff
          endif
        enddo
      enddo
      deallocate(vz,dzi)

      return
      end
c=============================================================================
      subroutine applyhele(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,ez,bx,by,bz)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np),bx(np),by(np),bz(np)

c Apply the transverse element hele.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey, ez   electric field
c   bx, by, bz   magnetic field

      integer(ISZ):: io,ip,iz,iele,ii
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8):: n,v,f,rpow,cosnt,sinnt
      real(kind=8):: xpmqoff,ypmqoff
      integer(ISZ), allocatable:: tiele(:)
      real(kind=8), allocatable, dimension(:):: vz,dzi,tfrac,txpmqoff,typmqoff,
     &                                          trpmqoff,ttpmqoff
      integer(ISZ):: allocerror

      if (.not. heles .or. .not. linhele(0)) return

      allocate(tiele(np),vz(npz),dzi(npz),tfrac(np),txpmqoff(np),typmqoff(np),
     &         trpmqoff(np),ttpmqoff(np),stat=allocerror)
      if (allocerror /= 0) then
        print*,"applyhele: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",np," and ",npz
        call kaboom("applyhele: allocation error")
        return
      endif

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nheleol
        if (.not. linhele(io)) cycle

c       --- calculate coordinates, element indices, and residence
c       --- fractions over particle block and store in temp arrays
c       --- These factors will be used repeatedly when the fields
c       --- are accumulated at each particle resulting from
        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= chelezs(iz,io) .and. zl < cheleze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= chelezs(iz,io) .and. zr < cheleze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cheleze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-chelezs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= chelezs(iz,io) .and. zl < cheleze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= chelezs(iz,io) .and. zr < cheleze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cheleze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-chelezs(iz,io))*dzi(ip)
          endif

c         --- find element index
          iele = cheleid(iz,io)
          tiele(ip) = iele
          tfrac(ip) = frac
c         --- set the field
          if (frac > 0.) then
c           --- x,y coordinates about multipole center
            xpmqoff = xp(ip) - heleox(iele)
            ypmqoff = yp(ip) - heleoy(iele)
            txpmqoff(ip) = xpmqoff
            typmqoff(ip) = ypmqoff
c           --- r,theta coordinates about multipole center
            trpmqoff(ip) = sqrt(xpmqoff**2 + ypmqoff**2)
            ttpmqoff(ip) = atan2(ypmqoff,dvnz(xpmqoff))
          endif
        enddo
c       --- accumulate field contributions at each particle in
c       --- block looping (outermost) over each multipole component
c       --- for vectorization
        do ii = 1, nhmlt
c         --- loop over particle block using stored coordinates, element
c         --- indices, and residence fractions
          do ip = 1,np
            if (tfrac(ip) == 0.) cycle
            iele = tiele(ip)
            n = hele_n(ii,iele)
            v = hele_v(ii,iele)
c           --- electric multipoles
            if (ii <= helene(iele)) then
              if (n == 0) then
c               --- accelerating field
c               --- Note that this really should not be used since the
c               --- accl element does a better job.  The transverse
c               --- fields should only be kicks at the entrance and exit.
                f = 1./(2.*v+2.)
                rpow = trpmqoff(ip)**(2*v)
                ex(ip) = ex(ip)-tfrac(ip)*heleep(ii,iele)*rpow*f*txpmqoff(ip)
                ey(ip) = ey(ip)-tfrac(ip)*heleep(ii,iele)*rpow*f*typmqoff(ip)
                ez(ip) = ez(ip)+tfrac(ip)*heleae(ii,iele)*rpow
              else if (n == 1 .and. v == 0) then
c               --- The pure dipole term needs special treatment since it
c               --- breaks down at r = 0.
                cosnt = cos(helepe(ii,iele))
                sinnt = sin(helepe(ii,iele))
                ex(ip) = ex(ip) + tfrac(ip)*heleae(ii,iele)*cosnt
                ey(ip) = ey(ip) + tfrac(ip)*heleae(ii,iele)*sinnt
              else
c               --- rest of components
                f = 1. + 2.*v/n
                cosnt = cos(n*(ttpmqoff(ip)-helepe(ii,iele)))
                sinnt = sin(n*(ttpmqoff(ip)-helepe(ii,iele)))
                rpow = trpmqoff(ip)**(n-2+2*v)
                ex(ip) = ex(ip) + tfrac(ip)*heleae(ii,iele)*rpow*
     &                         (f*txpmqoff(ip)*cosnt + typmqoff(ip)*sinnt)
                ey(ip) = ey(ip) + tfrac(ip)*heleae(ii,iele)*rpow*
     &                         (f*typmqoff(ip)*cosnt - txpmqoff(ip)*sinnt)
              endif
            endif
c           --- magnetic multipoles
            if (ii <= helenm(iele)) then
              if (n == 0) then
c               --- Solenoid field
c               --- The transverse fields should only be kicks at the entrance
c               --- and exit.
                f = 1./(2.*v+2.)
                rpow = trpmqoff(ip)**(2*v)
                bx(ip) = bx(ip)-tfrac(ip)*helemp(ii,iele)*rpow*f*txpmqoff(ip)
                by(ip) = by(ip)-tfrac(ip)*helemp(ii,iele)*rpow*f*typmqoff(ip)
                bz(ip) = bz(ip)+tfrac(ip)*heleam(ii,iele)*rpow
              else if (n == 1 .and. v == 0) then
c               --- The pure dipole term needs special treatment since it
c               --- breaks down at r = 0.
                cosnt = cos(helepm(ii,iele))
                sinnt = sin(helepm(ii,iele))
                bx(ip) = bx(ip) - tfrac(ip)*heleam(ii,iele)*sinnt
                by(ip) = by(ip) + tfrac(ip)*heleam(ii,iele)*cosnt
              else
c               --- rest of components
                f = 1. + 2.*v/n
                cosnt = cos(n*(ttpmqoff(ip)-helepm(ii,iele)))
                sinnt = sin(n*(ttpmqoff(ip)-helepm(ii,iele)))
                rpow = trpmqoff(ip)**(n-2+2*v)
                bx(ip) = bx(ip) + tfrac(ip)*heleam(ii,iele)*rpow*
     &                         (f*txpmqoff(ip)*sinnt - typmqoff(ip)*cosnt)
                by(ip) = by(ip) + tfrac(ip)*heleam(ii,iele)*rpow*
     &                         (f*typmqoff(ip)*sinnt + txpmqoff(ip)*cosnt)
              endif
            endif
          enddo
        enddo
      enddo
      deallocate(tiele,vz,dzi,tfrac,txpmqoff,typmqoff,trpmqoff,ttpmqoff)

      return
      end
c=============================================================================
      subroutine applyemlt(np,xp,yp,npz,zp,dtl,dtr,dt,lslice,ex,ey,ez)
      use Lattice
      use LatticeInternal
      use Mult_data
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np)

c Apply the element emlt.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey, ez   electric field

c  See HIF-note 96-10, equation D17, with E replacing B and
c  -n*psi replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Er = - sum_v  [ E'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Env*(1+2*v/n)*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Et = - sum_nv [ Env*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Ez = sum_v  [ E0v*r**(2*v) ] +
c         sum_nv [ E'nv/n*r**(2*v+n)*cos(n*(theta-psi)) ]
c
c    Ex = Er*cos(theta) - Et*sin(theta)
c    Ey = Er*sin(theta) + Et*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  Additional scale factors are included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factors are added together (in mltlocat) and multiply the applied field.
c  The addition is done in mltlocat so that it is only done once per particle
c  and not once per particle per multipole component.

      integer(ISZ):: io,in,ip
      real(kind=8):: n,v,f,fz,coeffp,coeff,rpow,alpha,alphap,cosnt,sinnt
      integer(ISZ), allocatable:: izm(:), imm(:)
      real(kind=8), allocatable, dimension(:):: xxx, yyy, zzz, ttt, wzm, rr, tt,
     &                                          mltph, mltsf, tex, tey, tez
      integer(ISZ):: allocerror

      if (.not. emlts .or. .not. linemlt(0)) return

      allocate(izm(np), imm(np), xxx(np), yyy(np), zzz(np), ttt(np), wzm(np), 
     &         rr(np), tt(np), mltph(np), mltsf(np), tex(np), tey(np), tez(np),
     &         stat=allocerror)
      if (allocerror /= 0) then
        print*,"applyemlt: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",np
        call kaboom("applyemlt: allocation error")
        return
      endif

      do io=1,nemltol
        if (.not. linemlt(io)) cycle

c       --- Get location of particle relative to moment data.
        call mltlocat(np,xp,yp,npz,zp,nemltsets,dzemlt,nzemlt,nzl,nzlmax,
     &                nemltol,io,cemltzs,cemltze,cemltph,cemltsf,cemltsc,
     &                cemltim,cemltox,cemltoy,imm,izm,wzm,xxx,yyy,zzz,ttt,
     &                rr,tt,mltph,mltsf,zlmin,zlframe,dzli,
     &                bends,linbend,cbendzs,cbendze,cbendrc,lslice)

c       --- zero temporary field arrays
        tex = 0.
        tey = 0.
        tez = 0.

c       --- accumulate the E field of the electrostatic multipoles
        do in=1,nesmult
          n = emlt_n(in)
          v = emlt_v(in)
          if (n == 0) then
c           --- Apply accelerating field and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                     (esemltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = rr(ip)**(2*v)
                tex(ip) = tex(ip) - coeffp*rpow*f*xxx(ip)
                tey(ip) = tey(ip) - coeffp*rpow*f*yyy(ip)
                tez(ip) = tez(ip) + coeff*rpow
              endif
            enddo
          else if (n==1 .and. v == 0) then
c           --- The pure dipole term needs special treatment since it
c           --- breaks down at r = 0.
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = mltph(ip) +
     &                      esemltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      esemltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (esemltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (esemltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(alpha)
                sinnt = sin(alpha)
                tex(ip) = tex(ip) + coeff*cosnt
                tey(ip) = tey(ip) + coeff*sinnt
                alpha = tt(ip) - alpha
                cosnt = cos(alpha)
                sinnt = sin(alpha)
                tez(ip) = tez(ip) + rr(ip)*(coeffp*cosnt - alphap*coeff*sinnt)
              endif
            enddo
          else if (n > 0 .and. v >= 0) then
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      esemltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      esemltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (esemltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (esemltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = rr(ip)**(n-2+2*v)
                tex(ip)=tex(ip) + coeff*rpow*(f*xxx(ip)*cosnt + yyy(ip)*sinnt)
                tey(ip)=tey(ip) + coeff*rpow*(f*yyy(ip)*cosnt - xxx(ip)*sinnt)
                tez(ip)=tez(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*cosnt - alphap*coeff*sinnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)*cos(ttt(ip)) - tez(ip)*sin(ttt(ip))
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tex(ip)*sin(ttt(ip)) + tez(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tez(ip)
          enddo
        endif

      enddo
      deallocate(izm, imm, xxx, yyy, zzz, ttt, wzm, 
     &           rr, tt, mltph, mltsf, tex, tey, tez)

      return
      end
c=============================================================================
      subroutine applymmlt(np,xp,yp,npz,zp,dtl,dtr,dt,lslice,bx,by,bz)
      use Lattice
      use LatticeInternal
      use Mult_data
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: bx(np),by(np),bz(np)

c Apply the element mmlt.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   bx, by, bz   magnetic field

c  See HIF-note 96-10, equation D17, with -n*psi-pi/2 replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Br = - sum_v  [ B'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Bnv*(1+2*v/n)*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Bt = sum_nv [ Bnv*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Bz = sum_v  [ B0v*r**(2*v) ] +
c         sum_nv [ B'nv/n*r**(2*v)*sin(n*(theta-psi)) ]
c
c    Bx = Br*cos(theta) - Bt*sin(theta)
c    By = Br*sin(theta) + Bt*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  Additional scale factors are included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factors are added together (in mltlocat) and multiply the applied field.
c  The addition is done in mltlocat so that it is only done once per particle
c  and not once per particle per multipole component.

      integer(ISZ):: io,in,ip
      real(kind=8):: n,v,f,fz,coeffp,coeff,rpow,alpha,alphap,cosnt,sinnt
      integer(ISZ), allocatable:: izm(:), imm(:)
      real(kind=8), allocatable, dimension(:):: xxx, yyy, zzz, ttt, wzm, rr, tt,
     &                                          mltph, mltsf, tbx, tby, tbz
      integer(ISZ):: allocerror

      if (.not. mmlts .or. .not. linmmlt(0)) return

      allocate(izm(np), imm(np), xxx(np), yyy(np), zzz(np), ttt(np), wzm(np), 
     &         rr(np), tt(np), mltph(np), mltsf(np), tbx(np), tby(np), tbz(np),
     &         stat=allocerror)
      if (allocerror /= 0) then
        print*,"applymmlt: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",np
        call kaboom("applymmlt: allocation error")
        return
      endif

      do io=1,nmmltol
        if (.not. linmmlt(io)) cycle

c       --- Get location of particle relative to moment data.
        call mltlocat(np,xp,yp,npz,zp,nmmltsets,dzmmlt,nzmmlt,nzl,nzlmax,
     &                nmmltol,io,cmmltzs,cmmltze,cmmltph,cmmltsf,cmmltsc,
     &                cmmltim,cmmltox,cmmltoy,imm,izm,wzm,xxx,yyy,zzz,ttt,
     &                rr,tt,mltph,mltsf,zlmin,zlframe,dzli,
     &                bends,linbend,cbendzs,cbendze,cbendrc,lslice)

c       --- zero the temporary field arrays
        tbx = 0.
        tby = 0.
        tbz = 0.

c       --- accumulate the B field of the magnetostatic multipoles
        do in=1,nmsmult
          n = mmlt_n(in)
          v = mmlt_v(in)
          if (n == 0) then
c           --- Apply solenoidal field component and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                    (msmmltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                     msmmltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                    (msmmlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                     msmmlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = rr(ip)**(2*v)
                tbx(ip) = tbx(ip) - coeffp*rpow*f*xxx(ip)
                tby(ip) = tby(ip) - coeffp*rpow*f*yyy(ip)
                tbz(ip) = tbz(ip) + coeff*rpow
              endif
            enddo
          else if (n == 1 .and. v == 0) then
c           --- The pure dipole term needs special treatment since it
c           --- breaks down at r = 0.
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = mltph(ip) +
     &                      msmmltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      msmmltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (msmmltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (msmmlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      msmmlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (msmmltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(alpha)
                sinnt = sin(alpha)
                tbx(ip) = tbx(ip) - coeff*rpow*sinnt
                tby(ip) = tby(ip) + coeff*rpow*cosnt
                alpha = tt(ip) - alpha
                cosnt = cos(alpha)
                sinnt = sin(alpha)
                tbz(ip) = tbz(ip) + rr(ip)*(coeffp*sinnt + alphap*coeff*cosnt)
              endif
            enddo
          else if (n > 0 .and. v >= 0) then
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      msmmltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      msmmltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (msmmltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (msmmlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      msmmlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (msmmltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = rr(ip)**(n-2+2*v)
                tbx(ip)=tbx(ip) + coeff*rpow*(f*xxx(ip)*sinnt - yyy(ip)*cosnt)
                tby(ip)=tby(ip) + coeff*rpow*(f*yyy(ip)*sinnt + xxx(ip)*cosnt)
                tbz(ip)=tbz(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*sinnt + alphap*coeff*cosnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)*cos(ttt(ip)) - tbz(ip)*sin(ttt(ip))
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbx(ip)*sin(ttt(ip)) + tbz(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbz(ip)
          enddo
        endif
      enddo
      deallocate(izm, imm, xxx, yyy, zzz, ttt, wzm, 
     &           rr, tt, mltph, mltsf, tbx, tby, tbz)

      return
      end
c=============================================================================
      subroutine mltlocat(np,xp,yp,npz,zp,nmltsets,dzmlt,nzmlt,nzl,nzlmax,
     &                    nol,io,cmltzs,cmltze,cmltph,cmltsf,cmltsc,
     &                    cmltim,cmltox,cmltoy,
     &                    imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                    mltph,mltsf,zlmin,zlframe,dzli,bends,linbend,
     &                    cbendzs,cbendze,cbendrc,lslice)
      integer(ISZ):: np,npz,nmltsets,nzl,nzlmax,nol,io
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dzmlt(nmltsets)
      integer(ISZ):: nzmlt(nmltsets)
      real(kind=8):: cmltzs(0:nzlmax,nol),cmltze(0:nzlmax,nol)
      real(kind=8):: cmltph(0:nzlmax,nol)
      real(kind=8):: cmltsf(0:nzlmax,nol),cmltsc(0:nzlmax,nol)
      integer(ISZ):: cmltim(0:nzlmax,nol)
      real(kind=8):: cmltox(0:nzlmax,nol),cmltoy(0:nzlmax,nol)
      integer(ISZ):: imm(np),izm(np)
      real(kind=8):: wzm(np),xxx(np),yyy(np),zzz(np)
      real(kind=8):: ttt(np),rr(np),tt(np),mltph(np),mltsf(np)
      logical(ISZ):: bends,linbend
      real(kind=8):: cbendzs(0:nzlmax),cbendze(0:nzlmax),cbendrc(0:nzlmax)
      real(kind=8):: zlmin,zlframe,dzli
      logical(ISZ):: lslice

c     --- Calculate grid location and polar coordinates of particles
c     --- for mlt elements.  Includes change of coordinates out of
c     --- warped coordinates (when a bend is between the
c     --- particle and the nearest mult).
c  Assumes that at most one bend will be between particle and the mult center.

      integer(ISZ):: iz,ip
      real(kind=8):: zz,cent,dzi,rrr,bzs,bze,bendrci

      if (lslice) then
        zz = zp(1)
        iz = 0
      endif
      if (.not. bends .or. .not. linbend) then
        do ip=1,np
          if (.not. lslice) then
            zz = zp(ip)
c           --- find z-cell in which particle lies
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
          endif
          imm(ip) = cmltim(iz,io)
          if (imm(ip) <= 0) cycle
          dzi = 1./dzmlt(imm(ip))
          mltph(ip) = cmltph(iz,io)
          mltsf(ip) = cmltsc(iz,io) + cmltsf(iz,io)
c         --- find z location relative to multipole data
          izm(ip) = int((zz-cmltzs(iz,io))*dzi)
          wzm(ip) =     (zz-cmltzs(iz,io))*dzi - izm(ip)
          if (zz < cmltzs(iz,io)) then
            izm(ip) = 0
            wzm(ip) = 0
            imm(ip) = 0
          elseif (zz > cmltze(iz,io)) then
            izm(ip) = nzmlt(imm(ip))
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Set temporary coordinates and add offset if particle within
c         --- an element.
          if (imm(ip) > 0) then
            xxx(ip) = xp(ip) - cmltox(iz,io)
            yyy(ip) = yp(ip) - cmltoy(iz,io)
            zzz(ip) = zz
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      else
        do ip=1,np
          if (.not. lslice) then
            zz = zp(ip)
c           --- find z-cell in which particle lies
            iz = max(0., (zz - zlmin - zlframe)*dzli)
          endif
c         --- set temporaries
          imm(ip) = cmltim(iz,io)
          if (imm(ip) <= 0) cycle
          dzi = 1./dzmlt(imm(ip))
          cent = 0.5*(cmltzs(iz,io) + cmltze(iz,io))
          xxx(ip) = xp(ip)
          zzz(ip) = zz
          if (cmltzs(iz,io) <= zz .and. zz <= cmltze(iz,io)) then
            mltph(ip) = cmltph(iz,io)
            mltsf(ip) = cmltsc(iz,io) + cmltsf(iz,io)
c           --- apply coordinate change left of mult, using index of (iz)
            if (zz <= cent) then
              bendrci = 1./cbendrc(iz)
c             --- find end of bend: if bend extends beyond center of mult, use
c             --- that as the end of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendze(iz) < cent) then
                bze = cbendze(iz)
              else
                bze = cent
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(iz)
              if (zz < cbendzs(iz) .and. cbendzs(iz) < cent) then
                ttt(ip) = (bze - cbendzs(iz))*bendrci
                zzz(ip) = bze - rrr*sin(ttt(ip)) -
     &                    (cbendzs(iz) - zz)*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz) -
     &                    (cbendzs(iz) - zz)*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zz < bze .and. cbendzs(iz) < cent) then
                ttt(ip) = (bze - zz)*bendrci
                zzz(ip) = bze - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz)
              endif
c           --- apply coordinate change right of mult, using index of (iz+1)
            elseif (zz >= cent) then
              bendrci = 1./cbendrc(iz+1)
c             --- find start of bend: if bend extends beyond center of mult, use
c             --- that as the start of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendzs(iz+1) < cent) then
                bzs = cent
              else
                bzs = cbendzs(iz+1)
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(iz+1)
              if (zz > cbendze(iz+1) .and. cbendze(iz+1) > cent) then
                ttt(ip) =  - (cbendze(iz+1) - bzs)*bendrci
                zzz(ip) = bzs - rrr*sin(ttt(ip)) +
     &                    (zz - cbendze(iz+1))*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz+1) +
     &                    (zz - cbendze(iz+1))*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zz > bzs .and. cbendze(iz+1) > cent) then
                ttt(ip) = - (zz - bzs)*bendrci
                zzz(ip) = bzs - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz+1)
              endif
            endif
          endif
c         --- end of coordinate transformation

          izm(ip) = int((zzz(ip)-cmltzs(iz,io))*dzi)
          wzm(ip) =     (zzz(ip)-cmltzs(iz,io))*dzi - izm(ip)
          if (zzz(ip) < cmltzs(iz,io)) then
            izm(ip) = 0
            wzm(ip) = 0
            imm(ip) = 0
          elseif (zzz(ip) > cmltze(iz,io)) then
            izm(ip) = nzmlt(imm(ip))
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Add offset
          if (imm(ip) > 0) then
            xxx(ip) = xxx(ip) - cmltox(iz,io)
            yyy(ip) = yp(ip) -  cmltoy(iz,io)
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      endif

      return
      end
c=========================================================================== 
      subroutine applyegrd(np,xp,yp,npz,zp,lslice,ex,ey,ez)
      use Lattice
      use LatticeInternal
      use EGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np)

c  Sets magnetic field for particles from data sets
c  containing Ex, Ey, and Ez on 3-D grids.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     ex,ey,ez       magnetic field

c  The field is:
c     Ex = u0 * v0 * w0 * egrdex(i  ,j  ,k  ,egrdid)
c        + u1 * v0 * w0 * egrdex(i+1,j  ,k  ,egrdid)
c        + u0 * v1 * w0 * egrdex(i  ,j+1,k  ,egrdid)
c        + ...

      integer(ISZ):: io,ip,i,j,k,iz,ii,id,jp1
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,efac,xsign,ysign,zsign,txp,typ,temp
      real(kind=8):: txp0,typ0
      real(kind=8):: tex,tey

c  Evaluation of E, vectorized over particles

      if (.not. egrds .or. .not. linegrd(0)) return

      do io=1,negrdol
        if (.not. linegrd(io)) cycle

        if (lslice) then
c         --- All particles are in the same z-cell
          zz = zp(1)
          iz = 0
          ii = cegrdid(iz,io)
          id = egrdid(ii)
          k = (zz - cegrdzs(iz,io))*egrddzi(id)
c         --- skip field accumulation if particles are outside of axial grid
          if (zz < cegrdzs(iz,io) .or. cegrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
            zz = zp(ip)
c           --- find the location of the particle in the internal lattice arrays
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cegrdid(iz,io)
            id = egrdid(ii)
            if (zz < cegrdzs(iz,io) .or. cegrdze(iz,io) < zz) cycle
          endif

c         --- Cycle if this element is turned off
          if (id <= 0) cycle

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - egrdox(ii)
          typ = yp(ip) - egrdoy(ii)
c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( egrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*egrdcp(ii) + typ*egrdsp(ii) 
            typ = -temp*egrdsp(ii) + typ*egrdcp(ii) 
          endif 

c         --- If data is in RZ, then replace x with the radius and y=ys
          if (egrdrz(id)) then
            txp0 = txp
            typ0 = typ
            txp = sqrt(txp**2 + typ**2)
            typ = egrdys(ii)
          endif

c         --- Shift coordinates to measure from the edge of the field grid
          txp = txp - egrdxs(ii)
          typ = typ - egrdys(ii)

c         --- Set default sign of E field
          xsign = 1.
          ysign = 1.
          zsign = 1.

c         --- If E is quadrupolar symmetric, make transformations.
c         --- When the particle is in one of the even quadrants (either
c         --- x<0 or y<0 but not both), the transformation is done ey
c         --- swapping x and y, and ey swapping Ex and Ey (done at the end
c         --- of the loop). The sign of Ez is also changed. In the third
c         --- quadrant, the signs of both Ex and Ey are changed.
c                 Quadrupole symmetries on field grid
c           
c             Quadrant      E_x            E_y             E_z 
c             --------------------------------------------------------
c                I          E_x( x, y,z)   E_y( x, y,z)    E_z( x, y,z) 
c                II         E_x(-x, y,z)  -E_y(-x, y,z)   -E_z(-x, y,z)  
c                III       -E_x(-x,-y,z)  -E_y(-x,-y,z)    E_z(-x,-y,z) 
c                IV        -E_x( x,-y,z)   E_x( x,-y,z)   -E_z( x,-y,z) 

          if (egrdsy(ii) == 2) then
c           --- Get quadrant that the particle is in.
            if (txp < 0.) then
              ysign = -1.
              txp = -txp
            endif
            if (typ < 0.) then
              xsign = -1.
              typ = -typ
            endif
c           --- If in even quadrant...
            if (xsign*ysign < 0.) then
c             --- Switch sign of Ez.
              zsign = -1.
c             --- Swap x and y
cc              temp = txp
cc              txp = typ
cc              typ = temp
            endif
          endif

c         --- find location of particle in E field grid
          i =  txp*egrddxi(id)
          j =  typ*egrddyi(id)
          k = (zz - cegrdzs(iz,io))*egrddzi(id)

c         --- Calculate linear weights
          u1 = txp*egrddxi(id) - i
          v1 = typ*egrddyi(id) - j
          w1 = (zz - cegrdzs(iz,io))*egrddzi(id) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the E field grid
          if ((0 <= i .and. i < egrdnx) .and.
     &        (0 <= j .and. j < egrdny .or. egrdrz(id)) .and. 
     &        (0 <= k .and. k < egrdnz)) then

            efac = egrdsc(ii) + egrdsf(ii)

            jp1 = j + 1
            if (egrdrz(id)) then
              jp1 = j
              v0 = 1.
              v1 = 0.
            endif

            tex = xsign*efac*(u0*v0*w0*egrdex(i  ,j  ,k  ,id) +
     &                        u1*v0*w0*egrdex(i+1,j  ,k  ,id) +
     &                        u0*v1*w0*egrdex(i  ,jp1,k  ,id) +
     &                        u1*v1*w0*egrdex(i+1,jp1,k  ,id) +
     &                        u0*v0*w1*egrdex(i  ,j  ,k+1,id) +
     &                        u1*v0*w1*egrdex(i+1,j  ,k+1,id) +
     &                        u0*v1*w1*egrdex(i  ,jp1,k+1,id) +
     &                        u1*v1*w1*egrdex(i+1,jp1,k+1,id))

            tey = ysign*efac*(u0*v0*w0*egrdey(i  ,j  ,k  ,id) +
     &                        u1*v0*w0*egrdey(i+1,j  ,k  ,id) +
     &                        u0*v1*w0*egrdey(i  ,jp1,k  ,id) +
     &                        u1*v1*w0*egrdey(i+1,jp1,k  ,id) +
     &                        u0*v0*w1*egrdey(i  ,j  ,k+1,id) +
     &                        u1*v0*w1*egrdey(i+1,j  ,k+1,id) +
     &                        u0*v1*w1*egrdey(i  ,jp1,k+1,id) +
     &                        u1*v1*w1*egrdey(i+1,jp1,k+1,id))

c            if (egrdsy(ii) == 2 .and. xsign*ysign < 0.) then
c              temp = tex 
c              tex = tey
c              tey = temp
c            endif 

c           --- rotate transverse field componets for correct lab frame 
c               orientation if the element is rotated.  
            if ( egrdph(ii) .ne. 0. ) then
              temp = tex 
              tex = temp*egrdcp(ii) - tey*egrdsp(ii) 
              tey = temp*egrdsp(ii) + tey*egrdcp(ii)         
            endif 

            if (egrdrz(id)) then
c             --- For the RZ case, tbx now holds Er. Convert it to Ex and Ey.
c             --- Note that this ignores Etheta (in tby) for now.
              if (txp .ne. 0) then
                tey = tex*typ0/txp
                tex = tex*txp0/txp
              else
                tey = 0.
                tex = tex
              endif
            endif

c           --- tranverse field 
            ex(ip) = ex(ip) + tex
            ey(ip) = ey(ip) + tey

            ez(ip) = ez(ip) + zsign*efac*(u0*v0*w0*egrdez(i  ,j  ,k  ,id) +
     &                                    u1*v0*w0*egrdez(i+1,j  ,k  ,id) +
     &                                    u0*v1*w0*egrdez(i  ,j+1,k  ,id) +
     &                                    u1*v1*w0*egrdez(i+1,j+1,k  ,id) +
     &                                    u0*v0*w1*egrdez(i  ,j  ,k+1,id) +
     &                                    u1*v0*w1*egrdez(i+1,j  ,k+1,id) +
     &                                    u0*v1*w1*egrdez(i  ,j+1,k+1,id) +
     &                                    u1*v1*w1*egrdez(i+1,j+1,k+1,id))

          endif

        enddo
      enddo

      return
      end
c=========================================================================== 
      subroutine applybgrd(np,xp,yp,npz,zp,lslice,bx,by,bz)
      use Lattice
      use LatticeInternal
      use BGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: bx(np),by(np),bz(np)

c  Sets magnetic field for particles from data sets
c  containing Bx, By, and Bz on 3-D grids.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     bx,by,bz       magnetic field

c  The field is:
c     Bx = u0 * v0 * w0 * bgrdbx(i  ,j  ,k  ,bgrdid)
c        + u1 * v0 * w0 * bgrdbx(i+1,j  ,k  ,bgrdid)
c        + u0 * v1 * w0 * bgrdbx(i  ,j+1,k  ,bgrdid)
c        + ...

      integer(ISZ):: io,ip,i,j,k,iz,ii,id,jp1
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,bfac,xsign,ysign,zsign,zcent,zlen
      real(kind=8):: xt,yt,zt,x1,y1,z1
      real(kind=8):: xt0,yt0
      real(kind=8):: bxt,byt,bzt,bx1,by1,bz1
      real(kind=8),allocatable:: ca(:),sa(:),ct(:),st(:),cp(:),sp(:)

c  Evaluation of B, vectorized over particles

      if (.not. bgrds .or. .not. linbgrd(0)) return

c     --- Precalculate sines and cosines for efficiency
      allocate(ca(0:nbgrd),sa(0:nbgrd))
      allocate(ct(0:nbgrd),st(0:nbgrd))
      allocate(cp(0:nbgrd),sp(0:nbgrd))
      ca = cos(bgrdph)
      sa = sin(bgrdph)
      ct = cos(bgrdot)
      st = sin(bgrdot)
      cp = cos(bgrdop)
      sp = sin(bgrdop)

      do io=1,nbgrdol
        if (.not. linbgrd(io)) cycle

        if (lslice) then
c         --- All particles are in the same z-cell
          zz = zp(1)
          iz = 0
          ii = cbgrdid(iz,io)
          id = bgrdid(ii)
c         --- skip field accumulation if particles are outside of axial grid
          if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
            zz = zp(ip)
c           --- find the location of the particle in the internal lattice arrays
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cbgrdid(iz,io)
            id = bgrdid(ii)
            if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
          endif

c         --- Cycle if this element is turned off
          if (id <= 0) cycle

          zcent = 0.5*(cbgrdzs(iz,io) + cbgrdze(iz,io))
          zlen = bgrdze(ii) - bgrdzs(ii)

c         --- find particle coordinates in frame of gridded field
          xt = xp(ip) - bgrdox(ii)
          yt = yp(ip) - bgrdoy(ii)
          zt = zz     - zcent

c         --- transverse rotation to take into account an active rotation 
c         --- of the field element. Rotation is relative to the z center.
c         --- Particles are rotated to the frame of the element.
c         --- Later, the field components accumulated are rotated back. 
c         --- First, rotate by bgrdph about the z axis.
          if (bgrdop(ii) .ne. 0.) then
            x1 = xt
            y1 = yt
            xt = +x1*cp(ii) + y1*sp(ii)
            yt = -x1*sp(ii) + y1*cp(ii)
          endif

c         --- then by bgrdot about the new y axis
          if (bgrdot(ii) .ne. 0.) then
            x1 = xt
            z1 = zt
            xt = x1*ct(ii) - z1*st(ii)
            zt = x1*st(ii) + z1*ct(ii)
          endif

c         --- then by bgrdop about the new z axis
          if (bgrdph(ii) .ne. 0.) then
            x1 = xt
            y1 = yt
            xt = +x1*ca(ii) + y1*sa(ii)
            yt = -x1*sa(ii) + y1*ca(ii)
          endif 

c         --- If data is in RZ, then replace x with the radius and y=ys
          if (bgrdrz(id)) then
            xt0 = xt
            yt0 = yt
            xt = sqrt(xt**2 + yt**2)
            yt = bgrdys(ii)
          endif

c         --- Shift coordinates to measure from the edge of the field grid
          xt = xt - bgrdxs(ii)
          yt = yt - bgrdys(ii)
          zt = zt + zlen/2.

c         --- Set default sign of B field
          xsign = 1.
          ysign = 1.
          zsign = 1.

c         --- If B is quadrupolar symmetric, make transformations.
c         --- When the particle is in one of the even quadrants (either
c         --- x<0 or y<0 but not both), the transformation is done by
c         --- swapping x and y, and by swapping Bx and By (done at the end
c         --- of the loop). The sign of Bz is also changed. In the third
c         --- quadrant, the signs of both Bx and By are changed.
c                 Quadrupole symmetries on field grid
c           
c             Quadrant      B_x            B_y             B_z 
c             --------------------------------------------------------
c                I          B_x( x, y,z)   B_y( x, y,z)    B_z( x, y,z) 
c                II         B_x(-x, y,z)  -B_y(-x, y,z)   -B_z(-x, y,z)  
c                III       -B_x(-x,-y,z)  -B_y(-x,-y,z)    B_z(-x,-y,z) 
c                IV        -B_x( x,-y,z)   B_x( x,-y,z)   -B_z( x,-y,z) 

          if (bgrdsy(ii) == 2) then
c           --- Get quadrant that the particle is in.
            if (xt < 0.) then
              ysign = -1.
              xt = -xt
            endif
            if (yt < 0.) then
              xsign = -1.
              yt = -yt
            endif
c           --- If in even quadrant...
            if (xsign*ysign < 0.) then
c             --- Switch sign of Bz.
              zsign = -1.
c             --- Swap x and y
cc              temp = xt
cc              xt = yt
cc              yt = temp
            endif
          endif

c         --- find location of particle in B field grid
          i = xt*bgrddxi(id)
          j = yt*bgrddyi(id)
          k = zt*bgrddzi(id)

c         --- Calculate linear weights
          u1 = xt*bgrddxi(id) - i
          v1 = yt*bgrddyi(id) - j
          w1 = zt*bgrddzi(id) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if ((0 <= i .and. i < bgrdnx) .and.
     &        (0 <= j .and. j < bgrdny .or. bgrdrz(id)) .and.
     &        (0 <= k .and. k < bgrdnz)) then

            bfac = bgrdsc(ii) + bgrdsf(ii)

            jp1 = j + 1
            if (bgrdrz(id)) then
              jp1 = j
              v0 = 1.
              v1 = 0.
            endif

            bxt = xsign*bfac*(u0*v0*w0*bgrdbx(i  ,j  ,k  ,id) +
     &                        u1*v0*w0*bgrdbx(i+1,j  ,k  ,id) +
     &                        u0*v1*w0*bgrdbx(i  ,jp1,k  ,id) +
     &                        u1*v1*w0*bgrdbx(i+1,jp1,k  ,id) +
     &                        u0*v0*w1*bgrdbx(i  ,j  ,k+1,id) +
     &                        u1*v0*w1*bgrdbx(i+1,j  ,k+1,id) +
     &                        u0*v1*w1*bgrdbx(i  ,jp1,k+1,id) +
     &                        u1*v1*w1*bgrdbx(i+1,jp1,k+1,id))

            byt = ysign*bfac*(u0*v0*w0*bgrdby(i  ,j  ,k  ,id) +
     &                        u1*v0*w0*bgrdby(i+1,j  ,k  ,id) +
     &                        u0*v1*w0*bgrdby(i  ,jp1,k  ,id) +
     &                        u1*v1*w0*bgrdby(i+1,jp1,k  ,id) +
     &                        u0*v0*w1*bgrdby(i  ,j  ,k+1,id) +
     &                        u1*v0*w1*bgrdby(i+1,j  ,k+1,id) +
     &                        u0*v1*w1*bgrdby(i  ,jp1,k+1,id) +
     &                        u1*v1*w1*bgrdby(i+1,jp1,k+1,id))

            bzt = zsign*bfac*(u0*v0*w0*bgrdbz(i  ,j  ,k  ,id) +
     &                        u1*v0*w0*bgrdbz(i+1,j  ,k  ,id) +
     &                        u0*v1*w0*bgrdbz(i  ,jp1,k  ,id) +
     &                        u1*v1*w0*bgrdbz(i+1,jp1,k  ,id) +
     &                        u0*v0*w1*bgrdbz(i  ,j  ,k+1,id) +
     &                        u1*v0*w1*bgrdbz(i+1,j  ,k+1,id) +
     &                        u0*v1*w1*bgrdbz(i  ,jp1,k+1,id) +
     &                        u1*v1*w1*bgrdbz(i+1,jp1,k+1,id))

c            if (bgrdsy(ii) == 2 .and. xsign*ysign < 0.) then
c              temp = bxt 
c              bxt = byt
c              byt = temp
c            endif 

            if (bgrdrz(id)) then
c             --- For the RZ case, bxt now holds Br. Convert it to Bx and By.
c             --- Note that this ignores Btheta (in byt) for now.
              if (xt .ne. 0) then
                byt = bxt*yt0/xt
                bxt = bxt*xt0/xt
              else
                byt = 0.
                bxt = bxt
              endif
            endif

c           --- Rotate field componets back to the lab frame. These are the
c           --- same rotations as above but in opposite order with the sign
c           --- of the angle reversed.
            if (bgrdph(ii) .ne. 0.) then
              bx1 = bxt
              by1 = byt
              bxt = +bx1*ca(ii) - by1*sa(ii)
              byt = +bx1*sa(ii) + by1*ca(ii)
            endif 

            if (bgrdot(ii) .ne. 0.) then
              bx1 = bxt
              bz1 = bzt
              bxt = +bx1*ct(ii) + bz1*st(ii)
              bzt = -bx1*st(ii) + bz1*ct(ii)
            endif

            if (bgrdop(ii) .ne. 0.) then
              bx1 = bxt
              by1 = byt
              bxt = +bx1*cp(ii) - by1*sp(ii)
              byt = +bx1*sp(ii) + by1*cp(ii)
            endif

c           --- add in the resulting field 
            bx(ip) = bx(ip) + bxt
            by(ip) = by(ip) + byt
            bz(ip) = bz(ip) + bzt

          endif

        enddo
      enddo

      deallocate(ca,sa)
      deallocate(ct,st)
      deallocate(cp,sp)

      return
      end
c============================================================================ 
      subroutine bgrdgetfulllength()
      use Lattice
      use LatticeInternal
      use BGRDdata

c Calculate the full axial length of the elements, including the extra length
c that may be taken up if the element has a rotation away from the z-axis.
c The results are put in bgrdfs and bgrdfe.

      integer(ISZ):: ib,id
      real(kind=8):: zcent
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax
      real(kind=8):: ca,sa,ct,st,cp,sp
      real(kind=8):: xl(8),yl(8),zl(8)
      real(kind=8):: xr(8),yr(8),zr(8)
      real(kind=8):: x1(8),y1(8),z1(8)
      real(kind=8):: x2(8),y2(8),z2(8)

      do ib = 0,nbgrd
        id = bgrdid(ib)

c       --- The transverse mins and maxs of the grid
        xmin = bgrdxs(ib)
        xmax = bgrdxs(ib) + bgrddx(id)*bgrdnx
        ymin = bgrdys(ib)
        ymax = bgrdys(ib) + bgrddy(id)*bgrdny

        if (bgrdrz(id)) then
c         --- For RZ, extend the mins and maxs to cover the full transverse
c         --- extent as if the grid were 3D. The y extent is the same as x,
c         --- but with the appropriate y center.
          xmin = bgrdxs(ib) - bgrddx(id)*bgrdnx
          ymin = xmin - bgrdxs(ib) + bgrdys(ib)
          ymax = xmax - bgrdxs(ib) + bgrdys(ib)
        endif

c       --- Get the zmin and zmax relative to the center of the grid.
        zcent = 0.5*(bgrdzs(ib) + bgrdze(ib))
        zmin = bgrdzs(ib) - zcent
        zmax = bgrdze(ib) - zcent

c       --- Generate the locations of the eight corners of the grid in
c       --- the rotated frame.
        xr = (/xmin,xmax,xmin,xmax,xmin,xmax,xmin,xmax/)
        yr = (/ymin,ymin,ymax,ymax,ymin,ymin,ymax,ymax/)
        zr = (/zmin,zmin,zmin,zmin,zmax,zmax,zmax,zmax/)

c       --- Precalculate the cosines and sines.
        ca = cos(bgrdph(ib))
        sa = sin(bgrdph(ib))
        ct = cos(bgrdot(ib))
        st = sin(bgrdot(ib))
        cp = cos(bgrdop(ib))
        sp = sin(bgrdop(ib))

c       --- Convert the corners into the lab frame by three successive
c       --- rotations, first by bgrdph about the z axis,
        x1 = +xr*ca - yr*sa
        y1 = +xr*sa + yr*ca
        z1 = zr
c       --- then by bgrdot about the new y axis,
        x2 = +x1*ct + z1*st
        y2 = y1
        z2 = -x1*st + z1*ct
c       --- then by bgrdop about the new z axis.
c       xl = +x2*cp - y2*sp
c       yl = +x2*sp + y2*cp
        zl = z2

        bgrdfs(ib) = zcent + minval(zl)
        bgrdfe(ib) = zcent + maxval(zl)

      enddo

      return
      end
c=========================================================================== 
      subroutine applybsqgrad(np,xp,yp,npz,zp,lslice,b)
      use Lattice
      use LatticeInternal
      use BSQGRADdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: b(bsqgradnc,np)

c  Sets magnetic field for particles from data sets
c  containing Bx, By, and Bz on 3-D grids. Grid can also contain
c  additional vector information. It is assumed that the additional
c  information has the same symmetry as the magnetic field.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     b           magnetic field and additional information

c  The field is:
c     b(:,ip) = u0 * v0 * w0 * bsqgrad(:,i  ,j  ,k  ,bsqgradid)
c             + u1 * v0 * w0 * bsqgrad(:,i+1,j  ,k  ,bsqgradid)
c             + u0 * v1 * w0 * bsqgrad(:,i  ,j+1,k  ,bsqgradid)
c             + ...

      integer(ISZ):: io,ip,i,j,k,iz,ii,id,ic
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,bfac,xsign,ysign,zsign,txp,typ,temp
      real(kind=8):: tb(bsqgradnc)

c  Evaluation of B, vectorized over particles

      if (.not. bsqgrads .or. .not. linbsqgrad(0)) return

      do io=1,nbsqgradol
        if (.not. linbsqgrad(io)) cycle

        if (lslice) then
c         --- All particles are in the same z-cell
          zz = zp(1)
          iz = 0
          ii = cbsqgradid(iz,io)
          id = bsqgradid(ii)
          k = (zz - cbsqgradzs(iz,io))*bsqgraddzi(id)
c         --- skip field accumulation if particles are outside of axial grid
          if (zz < cbsqgradzs(iz,io) .or. cbsqgradze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
            zz = zp(ip)
c           --- find the location of the particle in the internal lattice arrays
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cbsqgradid(iz,io)
            id = bsqgradid(ii)
            if (zz < cbsqgradzs(iz,io) .or. cbsqgradze(iz,io) < zz) cycle
          endif

c         --- Cycle if this element is turned off
          if (id <= 0) cycle

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - bsqgradox(ii)
          typ = yp(ip) - bsqgradoy(ii)
c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( bsqgradph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*bsqgradcp(ii) + typ*bsqgradsp(ii) 
            typ = -temp*bsqgradsp(ii) + typ*bsqgradcp(ii) 
          endif 

c         --- Shift coordinates to measure from the edge of the field grid
          txp = txp - bsqgradxs(ii)
          typ = typ - bsqgradys(ii)

c         --- Set default sign of B field
          xsign = 1.
          ysign = 1.
          zsign = 1.

c         --- If B is quadrupolar symmetric, make transformations.
c         --- When the particle is in one of the even quadrants (either
c         --- x<0 or y<0 but not both), the transformation is done by
c         --- swapping x and y, and by swapping Bx and By (done at the end
c         --- of the loop). The sign of Bz is also changed. In the third
c         --- quadrant, the signs of both Bx and By are changed.
c                 Quadrupole symmetries on field grid
c           
c             Quadrant      B_x            B_y             B_z 
c             --------------------------------------------------------
c                I          B_x( x, y,z)   B_y( x, y,z)    B_z( x, y,z) 
c                II         B_y( y,-x,z)  -B_x( y,-x,z)   -B_z( y,-x,z)  
c                III       -B_x(-x,-y,z)  -B_y(-x,-y,z)    B_z(-x,-y,z) 
c                IV        -B_y(-y, x,z)   B_x(-y, x,z)   -B_z(-y, x,z) 

          if (bsqgradsy(ii) == 2) then
c           --- Get quadrant that the particle is in.
            if (txp < 0.) then
              xsign = -1.
              txp = -txp
            endif
            if (typ < 0.) then
              ysign = -1.
              typ = -typ
            endif
c           --- If in even quadrant...
            if (xsign*ysign < 0.) then
c             --- Switch sign of Bz.
              zsign = -1.
c             --- Swap x and y
              temp = txp
              txp = typ
              typ = temp
            endif
          endif

c         --- find location of particle in B field grid
          i =  txp*bsqgraddxi(id)
          j =  typ*bsqgraddyi(id)
          k = (zz - cbsqgradzs(iz,io))*bsqgraddzi(id)

c         --- Calculate linear weights
          u1 = txp*bsqgraddxi(id) - i
          v1 = typ*bsqgraddyi(id) - j
          w1 = (zz - cbsqgradzs(iz,io))*bsqgraddzi(id) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if (0 <= i .and. i < bsqgradnx .and.
     &        0 <= j .and. j < bsqgradny .and.
     &        0 <= k .and. k < bsqgradnz) then

            bfac = bsqgradsc(ii) + bsqgradsf(ii)

            tb = bfac*(u0*v0*w0*bsqgrad(:,i  ,j  ,k  ,id) +
     &                 u1*v0*w0*bsqgrad(:,i+1,j  ,k  ,id) +
     &                 u0*v1*w0*bsqgrad(:,i  ,j+1,k  ,id) +
     &                 u1*v1*w0*bsqgrad(:,i+1,j+1,k  ,id) +
     &                 u0*v0*w1*bsqgrad(:,i  ,j  ,k+1,id) +
     &                 u1*v0*w1*bsqgrad(:,i+1,j  ,k+1,id) +
     &                 u0*v1*w1*bsqgrad(:,i  ,j+1,k+1,id) +
     &                 u1*v1*w1*bsqgrad(:,i+1,j+1,k+1,id))


           if(.false.) then
c --- original
            if (bsqgradsy(ii) == 2 .and. xsign*ysign < 0.) then
              do ic=0,bsqgradnc/3-1
                tb(1+ic*3) = xsign*tb(1+ic*3)
                tb(2+ic*3) = ysign*tb(2+ic*3)
                tb(3+ic*3) = zsign*tb(3+ic*3)
              enddo
              if (xsign*ysign < 0.) then
                do ic=0,bsqgradnc/3-1
                  temp = tb(1+ic*3)
                  tb(1+ic*3) = tb(2+ic*3)
                  tb(2+ic*3) = temp
                enddo
              endif 
            endif
           else
c --- modified
            if (bsqgradsy(ii) == 2) then
              do ic=0,bsqgradnc/3-1
                tb(1+ic*3) = zsign*tb(1+ic*3)
                tb(2+ic*3) = ysign*tb(2+ic*3)
                tb(3+ic*3) = xsign*tb(3+ic*3)
              enddo
              if (xsign*ysign < 0.) then
                do ic=0,bsqgradnc/3-1
                  temp = tb(3+ic*3)
                  tb(3+ic*3) = tb(2+ic*3)
                  tb(2+ic*3) = temp
                enddo
              endif 
            endif 
           end if

c           --- rotate transverse field componets for correct lab frame 
c               orientation if the element is rotated.  
            if ( bsqgradph(ii) .ne. 0. ) then
              do ic=0,bsqgradnc/3-1
                temp = tb(1+ic*3)
                tb(1+ic*3) = temp*bsqgradcp(ii) - tb(2+ic*3)*bsqgradsp(ii) 
                tb(2+ic*3) = temp*bsqgradsp(ii) + tb(2+ic*3)*bsqgradcp(ii)         
              enddo
            endif 

c           --- Accumulate the interpolated field.
            b(:,ip) = b(:,ip) + tb

          endif

        enddo
      enddo

      return
      end
c=========================================================================== 
      subroutine applypgrd(np,xp,yp,npz,zp,lslice,ex,ey,ez)
      use Lattice
      use LatticeInternal
      use PGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np)

c  Sets electric field for particles from data
c  sets containing the potential on a 3-D grid.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     ex,ey,ez    electric field

c  The field is:
c     Ex = u0 * v0 * w0 * ex(i  ,j  ,k  ,pgrdid)
c        + u1 * v0 * w0 * ex(i+1,j  ,k  ,pgrdid)
c        + u0 * v1 * w0 * ex(i  ,j+1,k  ,pgrdid)
c        + ...
c
c  Note that this routine is very similar to sete3d routine.  It might be
c  a good idea to combine them to reduce code complexity, but that would make
c  sete3d more complex and possibly slower.  Since the pgrd routine will not
c  be used very often and sete3d is used always, the author felt it was better
c  to have a seperate routine for pgrd.

      integer(ISZ):: io,ip,i,j,k,iz,ii,id,im1,jm1
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,xfac,yfac,zfac,txp,typ,temp

c  Evaluation of E, vectorized over particles

      if (.not. pgrds .or. .not. linpgrd(0)) return

      do io=1,npgrdol
        if (.not. linpgrd(io)) cycle

        if (lslice) then
c         --- find the location of the particle in the internal lattice arrays
          zz = zp(1)
          iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
          ii = cpgrdid(iz,io)
          id = pgrdid(ii)
          k = (zz - cpgrdzs(iz,io))*pgrddzi(id)
          if (zz < cpgrdzs(iz,io) .or. cpgrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
c           --- find the location of the particle in the internal lattice arrays
            zz = zp(ip)
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cpgrdid(iz,io)
            id = pgrdid(ii)
c           --- cycle if lost particle
            if (zz < cpgrdzs(iz,io) .or. cpgrdze(iz,io) < zz) cycle
          endif

c         --- Cycle if this element is turned off
          if (id <= 0) cycle

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - pgrdox(ii)
          typ = yp(ip) - pgrdoy(ii)

c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( pgrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*pgrdcp(ii) + typ*pgrdsp(ii) 
            typ = -temp*pgrdsp(ii) + typ*pgrdcp(ii) 
          endif 

c         --- find location of particle in potential grid
          i =  abs(txp - pgrdxs(ii))*pgrddxi(id)
          j =  abs(typ - pgrdys(ii))*pgrddyi(id)
          k = (zz - cpgrdzs(iz,io))*pgrddzi(id)

c         --- Calculate linear weights
          u1 =  abs(txp - pgrdxs(ii))*pgrddxi(id) - i
          v1 =  abs(typ - pgrdys(ii))*pgrddyi(id) - j
          w1 = (zz - cpgrdzs(iz,io))*pgrddzi(id) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if (0 <= i .and. i < pgrdnx .and.
     &        0 <= j .and. j < pgrdny .and.
     &        0 <= k .and. k < pgrdnz) then

c           --- make sure all indices refer to first quadrant
            im1 = i - 1
            jm1 = j - 1
            if (i == 0) im1 = 1
            if (j == 0) jm1 = 1

c           --- adjust sign of E field for approiate quadrant
            xfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddxi(id)*0.5
            yfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddyi(id)*0.5
            zfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddzi(id)*0.5
            if (xp(ip) < pgrdxs(ii)) xfac = -xfac
            if (yp(ip) < pgrdys(ii)) yfac = -yfac

            ex(ip) = ex(ip) + xfac*
     &            (u0*v0*w0*(pgrd(im1,j  ,k  ,id) - pgrd(i+1,j  ,k  ,id))
     &           + u1*v0*w0*(pgrd(i  ,j  ,k  ,id) - pgrd(i+2,j  ,k  ,id))
     &           + u0*v1*w0*(pgrd(im1,j+1,k  ,id) - pgrd(i+1,j+1,k  ,id))
     &           + u1*v1*w0*(pgrd(i  ,j+1,k  ,id) - pgrd(i+2,j+1,k  ,id))
     &           + u0*v0*w1*(pgrd(im1,j  ,k+1,id) - pgrd(i+1,j  ,k+1,id))
     &           + u1*v0*w1*(pgrd(i  ,j  ,k+1,id) - pgrd(i+2,j  ,k+1,id))
     &           + u0*v1*w1*(pgrd(im1,j+1,k+1,id) - pgrd(i+1,j+1,k+1,id))
     &           + u1*v1*w1*(pgrd(i  ,j+1,k+1,id) - pgrd(i+2,j+1,k+1,id)))

            ey(ip) = ey(ip) + yfac*
     &            (u0*v0*w0*(pgrd(i  ,jm1,k  ,id) - pgrd(i  ,j+1,k  ,id))
     &           + u1*v0*w0*(pgrd(i+1,jm1,k  ,id) - pgrd(i+1,j+1,k  ,id))
     &           + u0*v1*w0*(pgrd(i  ,j  ,k  ,id) - pgrd(i  ,j+2,k  ,id))
     &           + u1*v1*w0*(pgrd(i+1,j  ,k  ,id) - pgrd(i+1,j+2,k  ,id))
     &           + u0*v0*w1*(pgrd(i  ,jm1,k+1,id) - pgrd(i  ,j+1,k+1,id))
     &           + u1*v0*w1*(pgrd(i+1,jm1,k+1,id) - pgrd(i+1,j+1,k+1,id))
     &           + u0*v1*w1*(pgrd(i  ,j  ,k+1,id) - pgrd(i  ,j+2,k+1,id))
     &           + u1*v1*w1*(pgrd(i+1,j  ,k+1,id) - pgrd(i+1,j+2,k+1,id)))

            ez(ip) = ez(ip) + zfac*
     &            (u0*v0*w0*(pgrd(i  ,j  ,k-1,id) - pgrd(i  ,j  ,k+1,id))
     &           + u1*v0*w0*(pgrd(i+1,j  ,k-1,id) - pgrd(i+1,j  ,k+1,id))
     &           + u0*v1*w0*(pgrd(i  ,j+1,k-1,id) - pgrd(i  ,j+1,k+1,id))
     &           + u1*v1*w0*(pgrd(i+1,j+1,k-1,id) - pgrd(i+1,j+1,k+1,id))
     &           + u0*v0*w1*(pgrd(i  ,j  ,k  ,id) - pgrd(i  ,j  ,k+2,id))
     &           + u1*v0*w1*(pgrd(i+1,j  ,k  ,id) - pgrd(i+1,j  ,k+2,id))
     &           + u0*v1*w1*(pgrd(i  ,j+1,k  ,id) - pgrd(i  ,j+1,k+2,id))
     &           + u1*v1*w1*(pgrd(i+1,j+1,k  ,id) - pgrd(i+1,j+1,k+2,id)))

          endif

        enddo
      enddo

      return
      end
c=============================================================================
      subroutine fetchphi_from_pgrd(np,xp,yp,zp,p)
      use Lattice
      use LatticeInternal
      use PGRDdata
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np),p(np)

c Fetch phi from pgrd at points.

      integer(ISZ):: io,ip,i,j,k,iz,ii,id
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,txp,typ,temp

      if (.not. pgrds .or. .not. linpgrd(0)) return

      do io=1,npgrdol
        if (.not. linpgrd(io)) cycle

        do ip = 1, np

c         --- find the location of the particle in the internal lattice arrays
          zz = zp(ip)
          iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
          ii = cpgrdid(iz,io)
          id = pgrdid(ii)

c         --- Cycle if this element is turned off
          if (id <= 0) cycle

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - pgrdox(ii)
          typ = yp(ip) - pgrdoy(ii)

c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( pgrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*pgrdcp(ii) + typ*pgrdsp(ii) 
            typ = -temp*pgrdsp(ii) + typ*pgrdcp(ii) 
          endif 

c         --- find location of particle in potential grid
          i =  abs(txp - pgrdxs(ii))*pgrddxi(id)
          j =  abs(typ - pgrdys(ii))*pgrddyi(id)
          k = (zz - cpgrdzs(iz,io))*pgrddzi(id)

c         --- Calculate linear weights
          u1 =  abs(txp - pgrdxs(ii))*pgrddxi(id) - i
          v1 =  abs(typ - pgrdys(ii))*pgrddyi(id) - j
          w1 = (zz - cpgrdzs(iz,io))*pgrddzi(id) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the field grid
          if (0 <= i .and. i < pgrdnx .and.
     &        0 <= j .and. j < pgrdny .and.
     &        0 <= k .and. k < pgrdnz) then

            p(ip) = p(ip) + (pgrdsc(ii)+pgrdsf(ii)) *
     &            (u0*v0*w0*pgrd(i  ,j  ,k  ,id)
     &           + u1*v0*w0*pgrd(i+1,j  ,k  ,id)
     &           + u0*v1*w0*pgrd(i  ,j+1,k  ,id)
     &           + u1*v1*w0*pgrd(i+1,j+1,k  ,id)
     &           + u0*v0*w1*pgrd(i  ,j  ,k+1,id)
     &           + u1*v0*w1*pgrd(i+1,j  ,k+1,id)
     &           + u0*v1*w1*pgrd(i  ,j+1,k+1,id)
     &           + u1*v1*w1*pgrd(i+1,j+1,k+1,id))

          endif

        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applyaccl(np,xp,zp,uzp,gaminv,dtl,dtr,dt,q,m,lslice,ez)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: xp(np),zp(np),uzp(np),gaminv(np)
      real(kind=8):: dtl,dtr,dt,q,m
      logical(ISZ):: lslice
      real(kind=8):: ez(np)

c Apply the accelerating gap element accl.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c   qoverm   charge over mass
c Output:
c   ez       axial electric field

      real(kind=8):: qoverm
      real(kind=8):: dti,vz,gapez,frac,z1,z2,vn,zl,zr,cacclz,caccll,vni
      real(kind=8):: zs,ze,zz
      real(kind=8):: frac1,frac2
      integer(ISZ):: io,ip,iz

      if (.not. accls .or. .not. linaccl(0) .or. m == 0.) return

      qoverm = q/m

c     --- XXX This code needs to be fixed for the case of backwards moving
c     --- XXX particles. A partial fix has been made to prevent some errors.

      dti = 1./(dtr-dtl)

c     --- Finite length gaps
      if (.not. lacclzl) then

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

c         --- For velocity correction, first calculate the velocity at
c         --- time level n.
          do ip=1,np
            vz = uzp(ip)*gaminv(ip)
c           if (vz == 0.) vz = LARGEPOS
c           --- find z-cell in which particle lies
            iz = int(max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5))
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)

            if (vz >= 0.) then
              zz = zp(ip)
              zs = cacclzs(iz,io)
              ze = cacclze(iz,io)
            else
c             --- If a particle is moving backwards, flip the z orientation
c             --- around. This is done so that the comparisons in the if
c             --- statements below are correct. They assume that a particle is
c             --- moving in the direction from zs to ze.
              zz = -zp(ip)
              vz = -vz
              gapez = -gapez
              zs = -cacclze(iz,io)
              ze = -cacclzs(iz,io)
            endif

            z1 = zz + dtl*vz - 0.5*gapez*qoverm*dtl**2
            z2 = zz + dtl*vz
            if (zz <= zs) then
              vn = vz
            elseif (z1 < zs) then
              vn = sqrt(vz**2 + 2.*gapez*qoverm*(zz - zs))
            elseif (zz <= ze) then
              vn = vz - gapez*qoverm*dtl
            elseif (z2 < ze) then
              vn = 0.5*(vz + 0.5*gapez*qoverm*dt + sqrt((vz +
     &             0.5*gapez*qoverm*dt)**2 -
     &             4.*gapez*qoverm*(zz - ze)))
            else
              vn = vz
            endif

c           --- Calculate the fraction of time in the gap.  Cases inside
c           --- and outside are included implicitly in the max and min calls.
c           --- Note that the max's inside the sqrt are for idiot proofing.
            if (zz <= zs) then
              vni = 1./vn
              frac = max((+dtr + (zz - zs)*vni)*dti , 0.)
            elseif (zz <= ze) then
              frac1 = max(0.,min(1.,(dtr + 2.*(zz - zs)/
     &                    (vz + vn))*dti))
              frac2 = max(0.,min(1.,(-dtl + 2.*(ze - zz)/
     &                (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &                (ze - zz))) + vn))*dti))
              frac = min(frac1,frac2)
            else
              vni = 1./vn
              frac = max((-dtl - (zz - ze)*vni)*dti , 0.)
            endif

c           --- add acceleration field to Ez field
c           --- gapez is recalculated since its sign above is sometimes
c           --- changed (for a backward moving particle).
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
            ez(ip) = ez(ip) + gapez*frac

          enddo
        enddo

      else

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

c         --- Zero length gaps
          do ip = 1, np
            vz = uzp(ip)*gaminv(ip)
            iz = max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5)
c           --- "left" and "right" end of velocity advance step
            if (vz >= 0.) then
              zl = zp(ip) + vz*dtl
              zr = zp(ip) + vz*dtr
            else
              zl = zp(ip) + vz*dtr
              zr = zp(ip) + vz*dtl
            endif
c           --- Gap center and length
            cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
            caccll = cacclze(iz,io) - cacclzs(iz,io)
c           --- Add acceleration to velocity.
c           --- Calculate the change in velocity and convert that into
c           --- and Ez - the expression reduces to q*V/m for small V.
            if (zl <= cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              ez(ip) = ez(ip) + (sqrt(vz**2 + 2.*qoverm*gapez*caccll)-vz)/
     &                          qoverm*dti
            endif
          enddo
        enddo

      endif

      return
      end
c=========================================================================== 
      subroutine zgapcorr(np,zp,xp,uzp,gaminv,dtl,dtr,dt,m,q,time)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: zp(np),xp(np),uzp(np),gaminv(np)
      real(kind=8):: dtl,dtr,dt,m,q,time

      integer(ISZ):: io,ip,iz
      real(kind=8):: gapez,zp2,delta,qoverm,cacclz,caccll,deltav,zr
      real(kind=8):: zz,zs,ze
      real(kind=8):: vz(np),vzi(np)

      qoverm = q/m

c     Add residence correction to particle position when acceleration is done

      if (.not. accls .or. .not. linaccl(0)) return

      vz = uzp*gaminv
      vzi = 1./dvnz(vz)

c     --- Finite length gap
      if (.not. lacclzl) then

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

          do ip=1,np
c           --- find index of nearest accl element
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)

            if (vz(ip) >= 0.) then
              zz = zp(ip)
              zs = cacclzs(iz,io)
              ze = cacclze(iz,io)
            else
              zz = -zp(ip)
              zs = -cacclze(iz,io)
              ze = -cacclzs(iz,io)
              gapez = -gapez
              vz(ip) = -vz(ip)
              vzi(ip) = -vzi(ip)
            endif

c           --- Approximate position of particle on next step
c           --- with simple residence correction for gaps.
c           --- Calculate time of entrance or exit, delta.
c           --- Space charge forces are ignored.
            zp2 = zz + vz(ip)*dt
            delta = 0

            if (zz <= zs .and. zs <= zp2) then
              delta = (zs - zz)*vzi(ip)
            elseif (zs < zz .and. zz <= ze
     &          .and. ze <= zp2 + gapez*qoverm*dt*(dtr-dtl)) then
              delta = 2.*(ze - zz)/
     &                (sqrt((vz(ip) - gapez*qoverm*dtl)**2 +
     &                2.*(ze - zz)*gapez*qoverm) +
     &                vz(ip) - gapez*qoverm*dtl)
              zp2 = zp2 + gapez*qoverm*dt*(dtr-dtl)
            endif
c           --- calculate correction on particle position
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
            if (delta > 0.5*dt) delta = dt - delta
            if (zz <= zs .and. zs <= zp2) then
              zp(ip) = zp(ip) + 0.5*gapez*qoverm*delta**2
            elseif (zz <= ze .and. ze <= zp2) then
              zp(ip) = zp(ip) - 0.5*gapez*qoverm*delta**2
            endif
          enddo
        enddo

      else
c       --- Zero length gap

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

          do ip = 1, np
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
            cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
            caccll = cacclze(iz,io) - cacclzs(iz,io)
c           --- "right" end of velocity advance step
            zr = zp(ip) + vz(ip)*dtr
c           --- calculate correction on position
            if (zp(ip) > cacclz .and. cacclz > zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) - (cacclz - zp(ip))*vzi(ip)*deltav
            endif
            if (zr > cacclz .and. cacclz > zp(ip) + vz(ip)*dt) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) + (dt - (cacclz - zp(ip))*vzi(ip))*deltav
            endif
          enddo
        enddo

      endif

      return
      end
c=============================================================================
      subroutine applyacclxy(np,xp,zp,uzp,gaminv,dtp,dtl,dtr,m,q,dt,lslice,ez)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: dtl,dtr,m,q,dt
      real(kind=8):: xp(np),zp(np)
      real(kind=8):: uzp(np),gaminv(np),dtp(np)
      logical(ISZ):: lslice
      real(kind=8):: ez(np)

c Apply acceleration for finite length gaps.
c This is done seperately from the rest of the elements (which are done
c in extebxy) since the residence corrections depends on dtp which
c changes during the iteration dealing with changing Vz.  It would be much
c less efficient to put the entire extebxy within that iteration loop so
c the accl elements was seperated out.

      real(kind=8):: dti,oneodt,qoverm,moverq,gapez,z1,z2,vn,frac
      real(kind=8):: cacclz,caccll,zl,zr
      real(kind=8):: dtl_a,dtr_a,dti_a,vz_a
      integer(ISZ):: io,iz,ip

c  --- Note that the accl elements in the slice code are slightly different
c  --- than other hard edged elements since the residence correction is
c  --- based on fraction of time spent inside element rather than fraction of
c  --- distance.  So, the particles individual velocity and dt are used.

      if (.not. accls .or. .not. linaccl(0)) return

      dti = 1./(dtr-dtl)
      oneodt = 1./dt
      qoverm = q/m

c     --- All particles lie at the same z-cell
      iz = 0

      if (.not. lacclzl) then

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

c         --- For velocity correction, first calculate the velocity at
c         --- time level n.
          do ip=1,np
            dtl_a = dtl*oneodt*dtp(ip)
            dtr_a = dtr*oneodt*dtp(ip)
            dti_a = dti*dt/dtp(ip)
            vz_a = uzp(ip)*gaminv(ip)
c           if (vz_a == 0.) vz_a = LARGEPOS
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
            z1 = zlframe + dtl_a*vz_a - 0.5*gapez*qoverm*dtl_a**2
            z2 = zlframe + dtl_a*vz_a
            if (zlframe <= cacclzs(iz,io)) then
              vn = vz_a
            elseif (z1 < cacclzs(iz,io)) then
              vn = sqrt(vz_a**2 + 2.*gapez*qoverm*(zlframe - cacclzs(iz,io)))
            elseif (zlframe <= cacclze(iz,io)) then
              vn = vz_a - gapez*qoverm*dtl_a
            elseif (z2 < cacclze(iz,io)) then
              vn = 0.5*(vz_a + 0.5*gapez*qoverm*dtp(ip) + sqrt((vz_a +
     &             0.5*gapez*qoverm*dtp(ip))**2 -
     &             4.*gapez*qoverm*(zlframe - cacclze(iz,io))))
            else
              vn = vz_a
            endif
c           --- Calculate the fraction of time in the gap.  Cases inside
c           --- and outside are included implicitly in the max and min calls.
c           --- Note that the max's inside the sqrt are for idiot proofing.
            if (zlframe <= cacclzs(iz,io)) then
              frac = max((dtr_a + (zlframe - cacclzs(iz,io))/vn)*dti_a , 0.)
            elseif (zlframe <= (cacclzs(iz,io) + cacclze(iz,io))*0.5) then
              frac = min((dtr_a + 2.*(zlframe - cacclzs(iz,io))/
     &               (sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &               (zlframe - cacclzs(iz,io)))) + vn))*dti_a, 1.)
            elseif (zlframe <= cacclze(iz,io)) then
              frac = min((-dtl_a + 2.*(cacclze(iz,io) - zlframe)/
     &               (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &               (cacclze(iz,io) - zlframe))) + vn))*dti_a, 1.)
            else
              frac = max((-dtl_a - (zlframe - cacclze(iz,io))/vn)*dti_a , 0.)
            endif
c           --- add acceleration field to Ez field
            ez(ip) = ez(ip) + gapez*frac
          enddo
        enddo

      else

c       --- Zero-length gaps.
        do io=1,nacclol
          if (.not. linaccl(io)) cycle

          moverq = m/q
          cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
          caccll = cacclze(iz,io) - cacclzs(iz,io)
c         --- Add acceleration to velocity.
c         --- Calculate the change in velocity and convert that into
c         --- and Ez - the expression reduces to q*V/m for small V.
          do ip = 1, np
            dtl_a = dtl*oneodt*dtp(ip)
            dtr_a = dtr*oneodt*dtp(ip)
            dti_a = dti*dt/dtp(ip)
            vz_a = uzp(ip)*gaminv(ip)
            zl = zlframe + dtl_a*vz_a
            zr = zlframe + dtr_a*vz_a
            if (zl <= cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              ez(ip) = ez(ip) + (sqrt(vz_a**2 + 2.*qoverm*gapez*caccll) - vz_a)
     &                 *moverq*dti_a
            endif
          enddo
        enddo

      endif

c  ------------------------------------------------------------
c  End of acceleration
c  ------------------------------------------------------------

      return
      end
c=========================================================================== 
      subroutine zgapcorrxy(np,zp,xp,uzp,gaminv,dtp,dtl,dtr,dt,m,q,time)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: zp(np), xp(np), uzp(np), gaminv(np), dtp(np)
      real(kind=8):: dtl,dtr,dt,m,q,time

c Add residence correction to particle position when acceleration is done

      integer(ISZ):: io,ip,iz
      real(kind=8):: gapez,zp2,delta,qoverm,vz(np),vzi(np)
      real(kind=8):: dtl_a,dtr_a
      real(kind=8):: cacclz,caccll,zr,deltav

      if (.not. accls .or. .not. linaccl(0)) return

      qoverm = q/m
      vz = uzp*gaminv
      vzi = 1./dvnz(vz)

c     --- All particles lie at the same z-cell
      iz = 0

c     --- Finite length gap
      if (.not. lacclzl) then

        do io=1,nacclol
          if (.not. linaccl(io)) cycle

          do ip=1,np
c           --- Scale time step size by particles time step.
            dtl_a = dtl/dt*dtp(ip)
            dtr_a = dtr/dt*dtp(ip)
c           --- Fetch accelerating voltage
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
c           --- Approximate position of particle on next step
c           --- with simple residence correction for gaps.
c           --- Calculate time of entrance or exit, delta.
c           --- Space charge forces are ignored.
            zp2 = zlframe + vz(ip)*dtp(ip)
            delta = 0
            if (zlframe <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
              delta = (cacclzs(iz,io) - zlframe)*vzi(ip)
            elseif (cacclzs(iz,io) < zlframe .and. zlframe <= cacclze(iz,io)
     &   .and. cacclze(iz,io) <= zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)) then
              delta = 2.*(cacclze(iz,io) - zlframe)/
     &                (sqrt((vz(ip) - gapez*qoverm*dtl_a)**2 +
     &                2.*(cacclze(iz,io) - zlframe)*gapez*qoverm) +
     &                vz(ip) - gapez*qoverm*dtl_a)
              zp2 = zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)
            endif
c           --- calculate correction on particle position
            if (delta > 0.5*dtp(ip)) delta = dtp(ip) - delta
            if (zlframe <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
              zp(ip) = zp(ip) + 0.5*gapez*qoverm*delta**2
            elseif (zlframe <= cacclze(iz,io) .and. cacclze(iz,io) <= zp2) then
              zp(ip) = zp(ip) - 0.5*gapez*qoverm*delta**2
            endif
          enddo
        enddo

      else

c       --- Zero length gap
        do io=1,nacclol
          if (.not. linaccl(io)) cycle

          cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
          caccll = cacclze(iz,io) - cacclzs(iz,io)
          do ip = 1, np
c           --- "right" end of velocity advance step
            zr = zlframe + vz(ip)*dtr
c           --- calculate correction on position
            if (zlframe < cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) - (cacclz - zlframe)*vzi(ip)*deltav
            endif
            if (zr < cacclz .and. cacclz < zlframe + vz(ip)*dtp(ip)) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) + (dtp(ip) - (cacclz - zlframe)*vzi(ip))*deltav
            endif
          enddo
        enddo

      endif

      return
      end
c=============================================================================
      subroutine acclbfrm(zcorrection)
      use Constant
      use Beam_acc
      use InGen
      use InPart
      use Lattice
      use LatticeInternal
      use Picglb
      real(kind=8):: zcorrection
c Accelerate the beam frame using the accl lattice elements.
c Beam frame is treated as a particle of species 1 located at the center of
c the beam frame.
c
c Output argument, zcorrection, is the correction to the beam frame's
c longitudinal postion. Note that zbeam can not be directly changed here
c since that may force some particles off of the mesh. Instead, the
c zcorrection is returned so that it can be dealt with in a reasonable way
c by the calling routine.

      integer(ISZ):: io,astat,iz
      real(kind=8):: zs,ze,zc
      real(kind=8):: gapez,z1,z2,vn,frac,qoverm
      real(kind=8):: zz,zl,zr,vbeamfrmo,deltav
      real(kind=8):: zznew

c     --- This must be set to zero even if there are no elements since
c     --- it may be used in the calling code.
      zcorrection = 0.

      if (.not. accls .or. .not. linaccl(0) .or. vbeamfrm == 0.) return

      qoverm = echarge*zion/(aion*amu)

      if (acclbeamframe .ne. 0.) then
        zz = acclbeamframe + zlframe
        iz = nint((acclbeamframe-zlmin)/dzl)
        iz = max(0,min(nzl,iz))
      else
        iz = nzl/2
        zz = zlmin + iz*dzl + zlframe
      endif

      do io=1,nacclol
        if (.not. linaccl(io)) cycle

c       --- Check if this gap should accelerate the beam frame.
        if (cacclsw(iz,io) .ne. 0) cycle

c       --- Set some temporaries
        gapez = cacclez(iz,io)
        zs = cacclzs(iz,io)
        ze = cacclze(iz,io)

c       --- If not using zero length gaps then do residence corrections
c       --- across the gap.
        if (.not. lacclzl) then

c         --- Calculate grid frame velocity at time level n
          z1 = zz - 0.5*dt*vbeamfrm - 0.5*gapez*qoverm*(dt*0.5)**2
          z2 = zz - 0.5*dt*vbeamfrm
          if (zz <= zs) then
            vn = vbeamfrm
          elseif (z1 < zs) then
            vn = sqrt(vbeamfrm**2 + 2.*gapez*qoverm*(zz - zs))
          elseif (zz <= ze) then
            vn = vbeamfrm + 0.5*gapez*qoverm*dt
          elseif (z2 < ze) then
            vn = 0.5*(vbeamfrm + 0.5*gapez*qoverm*dt + sqrt((vbeamfrm +
     &           0.5*gapez*qoverm*dt)**2 -
     &           4.*gapez*qoverm*(zz - ze)))
          else
            vn = vbeamfrm
          endif
c         --- Calculate the fraction of time in the gap.  Cases inside
c         --- and outside are included implicitly in the max and min calls.
c         --- Note that the max's inside the sqrt are for idiot proofing.
          if (zz <= zs) then
            frac = max(0.5 + (zz - zs)/(vn*dt) , 0.)
          elseif (zz <= (zs + ze)*0.5) then
            frac = min(0.5 + 2.*(zz - zs)/
     &             (dt*(sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &             (zz - zs))) + vn)), 1.)
          elseif (zz <= ze) then
            frac = min(0.5 + 2.*(ze - zz)/
     &             (dt*(sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &             (ze - zz))) + vn)), 1.)
          else
            frac = max(0.5 - (zz - ze)/(vn*dt) , 0.)
          endif
c         --- do correction on position
          zznew = zz
          call zgapcorr(1,zznew,0.,vbeamfrm,1.,-dt*0.5,dt*0.5,dt,
     &                  aion*amu,zion*echarge,time)
          zcorrection = zcorrection + zznew - zz
c         --- add acceleration to beam frame velocity
          vbeamfrm = vbeamfrm + gapez*qoverm*dt*frac
c         --- also change beam frame velocity
          vbeam = vbeamfrm

        else
c         --- Zero-length gap

c         --- location of gap center
          zc = 0.5*(zs + ze)
c         --- "left" end of velocity advance step
          zl = zz - 0.5*vbeamfrm*dt
c         --- "right" end of velocity advance step
          zr = zz + 0.5*vbeamfrm*dt
c         --- Save current beam velocity
          vbeamfrmo = vbeamfrm
c         --- add acceleration to beam velocity
          if (zl .le. zc .and. zc .lt. zr) then
            vbeamfrm = sqrt(vbeamfrm**2 + 2.*qoverm*gapez*(ze - zs))
            vbeam = vbeamfrm
          endif
c         --- calculate correction on position of beam frame
          frac = 0.
          deltav = 0.
          if (zz < zc .and. zc < zr) then
            frac = - (zc - zz)/(vbeamfrmo*dt)
            deltav = vbeamfrm - vbeamfrmo
          endif
          if (zr < zc .and. zc < zz+vbeamfrmo*dt) then
            frac = 1. - (zc - zz)/(vbeamfrmo*dt)
            deltav = sqrt(vbeamfrm**2 + 2.*qoverm*gapez*(ze - zs))
     &               - vbeamfrm
          endif
c         --- Set the correction to the beam frame position
          zcorrection = zcorrection + frac*deltav*dt

        endif
      enddo

      return
      end
c=============================================================================

from warp import *
from lattice import *
drifts_version = "$Id: drifts.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

def initdrifts():
  if top.drfts: return

  # --- First, get list of all of the existing elements in one nice place.
  # --- This makes use of the classes setup in the lattice module.
  elems = []
  if top.bends:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nbend+1):
      elist.append(Bend(zs=top.bendzs[ii],ze=top.bendze[ii],ap=top.bendap[ii],
                        ox=0.,oy=0.,
                        rc=top.bendrc[ii]))
  if top.dipos:
    elist = []
    elems.append(elist)
    for ii in xrange(top.ndipo+1):
      elist.append(Dipo(zs=top.dipozs[ii],ze=top.dipoze[ii],ap=top.dipoap[ii],
                        ox=0.,oy=0.,
                        ex=top.dipoex[ii],ey=top.dipoey[ii],
                        bx=top.dipobx[ii],by=top.dipoby[ii]))
  if top.quads:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nquad+1):
      elist.append(Quad(zs=top.quadzs[ii],ze=top.quadze[ii],ap=top.quadap[ii],
                        ox=top.qoffx[ii],oy=top.qoffy[ii],
                        de=top.quadde[ii],db=top.quaddb[ii]))
  if top.sexts:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nsext+1):
      elist.append(Sext(zs=top.sextzs[ii],ze=top.sextze[ii],ap=top.sextap[ii],
                        ox=0.,oy=0.,
                        de=top.sextde[ii],db=top.sextdb[ii]))
  if top.heles:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nhele+1):
      elist.append(Hele(zs=top.helezs[ii],ze=top.heleze[ii],ap=top.heleap[ii],
                        ox=top.heleox[ii],oy=top.heleoy[ii],
                        ae=top.heleae[ii],am=top.heleam[ii]))
  if top.accls:
    elist = []
    elems.append(elist)
    for ii in xrange(top.naccl+1):
      elist.append(Accl(zs=top.acclzs[ii],ze=top.acclze[ii],ap=top.acclap[ii],
                        ox=top.acclox[ii],oy=top.accloy[ii],
                        ez=top.acclez[ii]))
  if top.emlts:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nemlt+1):
      elist.append(Emlt(zs=top.emltzs[ii],ze=top.emltze[ii],ap=top.emltap[ii],
                        ox=top.emltox[ii],oy=top.emltoy[ii],
                        id=top.emltid[ii]))
  if top.mmlts:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nmmlt+1):
      elist.append(Mmlt(zs=top.mmltzs[ii],ze=top.mmltze[ii],ap=top.mmltap[ii],
                        ox=top.mmltox[ii],oy=top.mmltoy[ii],
                        id=top.mmltid[ii]))
  if top.bgrds:
    elist = []
    elems.append(elist)
    for ii in xrange(top.nbgrd+1):
      elist.append(Bgrd(zs=top.bgrdzs[ii],ze=top.bgrdze[ii],ap=top.bgrdap[ii],
                        ox=top.bgrdox[ii],oy=top.bgrdoy[ii],
                        id=top.bgrdid[ii]))
  if top.pgrds:
    elist = []
    elems.append(elist)
    for ii in xrange(top.npgrd+1):
      elist.append(Pgrd(zs=top.pgrdzs[ii],ze=top.pgrdze[ii],ap=top.pgrdap[ii],
                        ox=top.pgrdox[ii],oy=top.pgrdoy[ii],
                        id=top.pgrdid[ii]))

  # --- Create a counter for each element type
  ielem = len(elems)*[0]

  # --- Get start and end of first element, whatever type it may be.
  elemzs = top.zlatperi
  for ie in xrange(len(elems)):
    e = elems[ie][0]
    if e.zs < elemzs:
      elemzs = e.zs
      elemze = e.ze
      celem = e
  while 1:
    for ie in xrange(len(elems)):
      if ielem[ie] < len(elems[ie]):
        e = elems[ie][ielem[ie]]
        if e.zs < elemze:
          ielem[ie] = ielem[ie] + 1
          if e.ze > elemze:
            elemze = e.ze
            celem = e
            break
    break

  # --- Create list of drifts.
  driftlist = []

  # --- Now, create fist drift if the first element starts after zero.
  if elemzs > 0:
    driftlist.append(Drft(zs=0.,ze=elemzs,ap=celem.aperture,ox=celem.offset_x,oy=celem.offset_y))

  # --- Set last z location of any element
  lastz = elemze

  # --- Find the rest of the drift spaces between elements.
  while 1:
    elemzs = top.zlatperi
    for ie in xrange(len(elems)):
      if ielem[ie] < len(elems[ie]):
        e = elems[ie][ielem[ie]]
        if e.zs < elemzs:
          elemzs = e.zs
          elemze = e.ze
          celem = e
    while 1:
      for ie in xrange(len(elems)):
        if ielem[ie] < len(elems[ie]):
          e = elems[ie][ielem[ie]]
          if e.zs < elemze:
            ielem[ie] = ielem[ie] + 1
            if e.ze > elemze:
              elemze = e.ze
              celem = e
              break
      break
    driftlist.append(Drft(zs=lastz,ze=elemzs,ap=celem.aperture,ox=celem.offset_x,oy=celem.offset_y))
    lastz = elemze
    if elemzs == top.zlatperi:
      break


  # --- Copy the data into the builtin Lattice arrays
  ndrft = len(driftlist)-1
  top.ndrft = ndrft
  gchange("Lattice")
  top.drfts = 0
  for ii in xrange(ndrft+1):
    top.drftzs[ii] = driftlist[ii].zs
    top.drftze[ii] = driftlist[ii].ze
    top.drftox[ii] = driftlist[ii].offset_x
    top.drftoy[ii] = driftlist[ii].offset_y
    top.drftap[ii] = driftlist[ii].aperture
    if top.drftze[ii] > top.drftzs[ii]: top.drfts = 1


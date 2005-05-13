"""
# Module b_fields.py
#
# by: Rami A. Kishek (based on algorithm by I. Haber)
# created: September 22, 2000
#
#	Last Modified: Jan. 5, 2005
#
# This script provides the following functions used in to generate
# bgrd data for solenoids using a formula
#
#   make_sol	... Uses Kishek's formula to generate a short solenoid
#       and derive nonlinear terms.  Parameters can be supplied.
#
#   make_long_sole  	... Uses Y. Zou's formula to generate a long
#       solenoid and derive nonlinear terms.  Parameters can be supplied.
#
#   plot_bofz	... Plot magnetic field of brgd elements vs z
#   plot_bofx	... Plot magnetic field of brgd elements vs x
#   plot_bofy	... Plot magnetic field of brgd elements vs y
#   comp_bofz   ... Plots comparisons of different bgrdids
#
# ==================
"""
b_fields_version = "$Id: b_fields.py,v 1.4 2005/05/13 06:01:16 ramiak Exp $"
from warp import *
from rami_scripts import *

def b_fieldsdoc():
  import b_fields
  print b_fields.__doc__

#=====================

def make_sol( ispec=0, bcn=3.43e-02, dcn=4.82e-02, ccn=0.0170, bz0=76.0,
              llin=no, zshift=0, hwid=None, zlen=None):
    """ make_sol( ispec=0, bcn=3.43e-02, dcn=4.82e-02, ccn=0.0170,
                  bz0=76.0, llin=no, zshift=0, hwid=None, zlen=None)

# --- Calculate actual bgrd fields for solenoid for an assumed on-axis
        z magnetic field of the form:

  Bz(0) = bz0*exp(-(z/dcn)**2)*[1./cosh(zz/bcn) + ccn*(sinh(zz/bcn))**2]

    The magnetic field for all positions is calculated by a series
expansion of the axial fields.

Defaults correspond to UMER injector solenoid
- bcn, dcn are fitting parameters to experimental data, in m
- bz0 in Gauss
- zshift is distance in m left (-) or right (+) that solenoid
      center is shifted.  Zero for symmetric truncation.
- zlen=top.bgrdnz*top.bgrddz[ispec] is the truncation length
- hwid=w3d.xmmax is the half width of the bgrd array
- llin is a logical flag for discarding nonlinear terms of the expansion
    """
    if not hwid: hwid = w3d.xmmax
    if not zlen: zlen = top.bgrdnz*top.bgrddz[ispec]

    # --- Set up vector temporaries
    xxi = -hwid + top.bgrddx[ispec]*iota(0, top.bgrdnx)

    z_ctr = top.bgrddz[ispec]*top.bgrdnz/2.0 + zshift
    bz0_ext = bz0*1.0e-4    # convert to Tesla

    in_start = nint((z_ctr-0.5*zlen)/top.bgrddz[ispec])
    in_end   = nint((z_ctr+0.5*zlen)/top.bgrddz[ispec])

    # -- Calculate Solenoid field by expanding formula
    #
    # ---First the coefficients at a particular z are calculated.
    for iiz in range( in_start, in_end+1 ):
        zz = iiz*top.bgrddz[ispec] - z_ctr
        g1 = bz0_ext*exp(-(zz/dcn)**2)
        g2 = 1./cosh(zz/bcn)
        g3 = ccn*(sinh(zz/bcn))**2
        dcn2 = dcn*dcn
        aleph = -2./dcn2
        ba = 1./bcn
        ta = 1./bcn
        #
        g1p = aleph*zz*g1
        g2p = -ba*(g2*g2*sinh(zz/bcn))
        g3p = 2.*ccn*ta*sinh(zz/bcn)*cosh(zz/bcn)
        g1pp = aleph*(g1 + zz*g1p)
        g2pp = -ba*g2*(ba + 2.*g2p*sinh(zz/bcn))
        g3pp = 2.*ta*ta*(g3 + ccn*(cosh(zz/bcn))**2)
        g1ppp = aleph*(2.*g1p + zz*g1pp)
        g2ppp = -3.*ba*ba*g2p - 2.*ba*(g2*g2pp + g2p*g2p)*sinh(zz/bcn)
        g3ppp = 4.*ta*ta*g3p
        g1pppp = aleph*(3.*g1pp + zz*g1ppp)
        g2pppp = -5*ba*ba*g2pp - 2.*ba*ba*g2p*g2p*cosh(zz/bcn) - \
                               2.*ba*(3*g2p*g2pp + g2*g2ppp)*sinh(zz/bcn)
        g3pppp = 4*ta*ta*g3pp
        #
        bcofs = [g1*(g2+g3),
                 g1p*(g2+g3) + g1*(g2p+g3p),
                 g1pp*(g2+g3) + 2.*g1p*(g2p+g3p) + g1*(g2pp+g3pp),
                 g1ppp*(g2+g3) + 3.*g1pp*(g2p+g3p) + 3.*g1p*(g2pp+g3pp) + g1*(g2ppp+g3ppp),
                 g1pppp*(g2+g3) + 4.*g1ppp*(g2p+g3p) + 6.*g1pp*(g2pp+g3pp) +
                                  4.*g1p*(g2ppp+g3ppp) + g1*(g2pppp+g3pppp)]
        bcof1 = bcofs[0]
        bcof2 = bcofs[1]*0.5
        bcof3 = bcofs[2]*0.25
        bcof4 = bcofs[3]/16.0
        bcof5 = bcofs[4]/64.
        #
        # ---After the coefficients are calculated the fields at a position x,y
        # ---are calculated from the off-axis expansion
        #
        for iiy in range( 0, top.bgrdny+1 ):
            # --- Inner loop is done with  vectors
            yy = -hwid + iiy*top.bgrddy[ispec]
            r2i = xxi[:]**2 + yy**2
            if llin:
                top.bgrdbx[:,iiy,iiz,ispec] = top.bgrdbx[:,iiy,iiz,ispec] + (-bcof2)*xxi[:]
                top.bgrdby[:,iiy,iiz,ispec] = top.bgrdby[:,iiy,iiz,ispec] + (-bcof2)*yy
                top.bgrdbz[:,iiy,iiz,ispec] = top.bgrdbz[:,iiy,iiz,ispec] + bcof1
            else:
                top.bgrdbx[:,iiy,iiz,ispec] = top.bgrdbx[:,iiy,iiz,ispec] + \
						                        (-bcof2 + bcof4*r2i[:])*xxi[:]
                top.bgrdby[:,iiy,iiz,ispec] = top.bgrdby[:,iiy,iiz,ispec] + \
						                        (-bcof2 + bcof4*r2i[:])*yy
                top.bgrdbz[:,iiy,iiz,ispec] = top.bgrdbz[:,iiy,iiz,ispec] + bcof1 - \
                                                bcof3*r2i[:] + bcof5*r2i[:]**2

#=====================

def make_long_sole( ispec=0, aaa=5.0408e-2, ccc=0.5027, lact=137.08e-2,
                    bz0=50.0, llin=no, zshift=0, hwid=None, zlen=None):
    """ make_long_sole( ispec=0, aaa=5.0408e-2, ccc=0.5027, lact=137.08e-2,
                        bz0=50.0, llin=no, zshift=0, hwid=None, zlen=None)

# --- Calculate actual bgrd fields for a long solenoid for an
        assumed on-axis z magnetic field of the form:

  Bz(0) = bz0*ccc*[(z+hlact)/sqrt(z+hlact)**2+aaa2))
                    - (z-hlact)/sqrt((z-hlact)**2+aaa2))]
     where:   hlact = lact/2.0

    The magnetic field for all positions is calculated by a series
expansion of the axial fields.

Defaults correspond to UMER injector solenoid
- aaa is a fitting parameter to experimental data, in m
- ccc is a dimensionless normalization parameter
- lact is the effective length, in m
- bz0 in Gauss
- zshift is distance in m left (-) or right (+) that solenoid
      center is shifted.  Zero for symmetric truncation.
- zlen=top.bgrdnz*top.bgrddz[ispec] is the truncation length
- hwid=w3d.xmmax is the half width of the bgrd array
- llin is a logical flag for discarding nonlinear terms of the expansion
    """
    if not hwid: hwid = w3d.xmmax
    if not zlen: zlen = top.bgrdnz*top.bgrddz[ispec]

    # --- Set up vector temporaries
    xxi = -hwid + top.bgrddx[ispec]*iota(0, top.bgrdnx)

    z_ctr = top.bgrddz[ispec]*top.bgrdnz/2.0 + zshift
    bz0_ext = bz0*1.0e-4    # convert to Tesla

    in_start = nint((z_ctr-0.5*zlen)/top.bgrddz[ispec])
    in_end   = nint((z_ctr+0.5*zlen)/top.bgrddz[ispec])

    # -- Calculate Solenoid field by expanding formula
    #
    # ---First the coefficients at a particular z are calculated.
    for iiz in range( in_start, in_end+1 ):
        zz = iiz*top.bgrddz[ispec] - z_ctr
        aaa2 = aaa*aaa
        hlact = lact/2.0
        g1 = bz0_ext*ccc
        g2 = (zz+hlact)/((zz+hlact)**2+aaa2)**0.5
        g3 = - (zz-hlact)/((zz-hlact)**2+aaa2)**0.5
        #
        g2p = 1.0/((zz+hlact)**2+aaa2)**0.5 - (zz+hlact)**2*((zz+hlact)**2+aaa2)**(-1.5)
        g2pp = -3.0*(zz+hlact)*((zz+hlact)**2+aaa2)**(-1.5) +   \
                            3.0*((zz+hlact)**3)*((zz+hlact)**2+aaa2)**(-2.5)
        g2ppp = -3.0*((zz+hlact)**2+aaa2)**(-1.5) + \
                            18.0*(zz+hlact)**2*((zz+hlact)**2+aaa2)**(-2.5) -   \
		                    15.0*((zz+hlact)**4)*((zz+hlact)**2+aaa2)**(-3.5)
        g2pppp = 45.0*(zz+hlact)*((zz+hlact)**2+aaa2)**(-2.5) - \
		                    150*(zz+hlact)**3*((zz+hlact)**2+aaa2)**(-3.5) +    \
		                    105.0*(zz+hlact)**5*((zz+hlact)**2+aaa2)**(-4.5)
        #
        g3p = -1.0/((zz-hlact)**2+aaa2)**0.5 + (zz-hlact)**2*((zz-hlact)**2+aaa2)**(-1.5)
        g3pp = 3.0*(zz-hlact)*((zz-hlact)**2+aaa2)**(-1.5) -    \
		                    3.0*((zz-hlact)**3)*((zz-hlact)**2+aaa2)**(-2.5)
        g3ppp = 3.0*((zz-hlact)**2+aaa2)**(-1.5) -  \
		                    18.0*(zz-hlact)**2*((zz-hlact)**2+aaa2)**(-2.5) +   \
 		                    15.0*((zz-hlact)**4)*((zz-hlact)**2+aaa2)**(-3.5)
        g3pppp = -45.0*(zz-hlact)*((zz-hlact)**2+aaa2)**(-2.5) +    \
		                    150*(zz-hlact)**3*((zz-hlact)**2+aaa2)**(-3.5) -    \
		                    105.0*(zz-hlact)**5*((zz-hlact)**2+aaa2)**(-4.5)
        #
        bcofs = [g1*(g2+g3), g1*(g2p+g3p), g1*(g2pp+g3pp),
                             g1*(g2ppp+g3ppp), g1*(g2pppp+g3pppp)]
        bcof1 = bcofs[0]
        bcof2 = bcofs[1]*0.5
        bcof3 = bcofs[2]*0.25
        bcof4 = bcofs[3]/16.0
        bcof5 = bcofs[4]/64.
        #
        # ---After the coefficients are calculated the fields at a position x,y
        # ---are calculated from the off-axis expansion
        #
        for iiy in range( 0, top.bgrdny+1 ):
            # --- Inner loop is done with  vectors
            yy = -hwid + iiy*top.bgrddy[ispec]
            r2i = xxi[:]**2 + yy**2
            if llin:
                top.bgrdbx[:,iiy,iiz,ispec] = top.bgrdbx[:,iiy,iiz,ispec] + (-bcof2)*xxi[:]
                top.bgrdby[:,iiy,iiz,ispec] = top.bgrdby[:,iiy,iiz,ispec] + (-bcof2)*yy
                top.bgrdbz[:,iiy,iiz,ispec] = top.bgrdbz[:,iiy,iiz,ispec] + bcof1
            else:
                top.bgrdbx[:,iiy,iiz,ispec] = top.bgrdbx[:,iiy,iiz,ispec] + \
						                        (-bcof2 + bcof4*r2i[:])*xxi[:]
                top.bgrdby[:,iiy,iiz,ispec] = top.bgrdby[:,iiy,iiz,ispec] + \
						                        (-bcof2 + bcof4*r2i[:])*yy
                top.bgrdbz[:,iiy,iiz,ispec] = top.bgrdbz[:,iiy,iiz,ispec] + bcof1 - \
                                                bcof3*r2i[:] + bcof5*r2i[:]**2


#===========

def plot_bofx(nspec=0, yplot=0, zplot=0, plot_title="", kwdict={},  **kw):
    """ plot_bofx(nspec=0, yplot=0, zplot=0, plot_title="", kwdict={},  **kw)
    Plots bgrd field of element species 'nspec' vs. axis // x
yplot, zplot define coords of desired plot axis
    - 'xscale': 100.
    - 'xoffset': 0.0
    - 'yscale': 1.0e4
    - 'yoffset': 0.0
    - 'width': 2.0
    - 'marks': 1
    - 'msize': 1.0
    - 'titleb', 'titlel', 'titlet', 'titler'
    - 'titles': 1
    """
    # --- Dictionary specifying plot defaults
    pldef = {'xscale': 100., 'xoffset': 0.0, 'yscale': 1.0e4, 'yoffset': 0.0,
             'width': 2.0, 'marks': 1, 'msize': 1.0,
             'titleb': "X (cm)", 'titlel': "B (Gauss)",
             'titlet': "B-field of Bgrd element vs. x", 'titler': "", 'titles': 1}
    if plot_title <> "": pldef['titlet'] = plot_title
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    for key in pldef.keys():    exec key+"=pldef['"+key+"']"
    #
    iiz = nint(zplot/top.bgrddz[nspec]) + nint(top.bgrdnz/2)
    iiy = nint(yplot/top.bgrddy[nspec]) + nint(top.bgrdny/2)
    #
    # --- Calculate x positions in grid to use for plotting fields
    xgrd = iota(-top.bgrdnx/2, top.bgrdnx/2)*top.bgrddx[nspec]
    #
    pldefault(width=width, marks=marks, msize=msize)
    plg( yoffset+yscale*top.bgrdbx[:,iiy,iiz,nspec],
                                xoffset+xscale*xgrd[:], marker="x", type="solid")
    plg( yoffset+yscale*top.bgrdby[:,iiy,iiz,nspec],
                                xoffset+xscale*xgrd[:], marker="y", type="dot")
    plg( yoffset+yscale*top.bgrdbz[:,iiy,iiz,nspec],
                                xoffset+xscale*xgrd[:], marker="z", type="dash")
    if titles:  ptitles(titlet, titleb, titlel, titler)
    #
    fma()

#===========

def plot_bofy(nspec=0, xplot=0, zplot=0, plot_title="", kwdict={},  **kw):
    """ plot_bofy(nspec=0, xplot=0, zplot=0, plot_title="", kwdict={},  **kw)
    Plots bgrd field of element species 'nspec' vs. axis // x
yplot, zplot define coords of desired plot axis
    - 'xscale': 100.
    - 'xoffset': 0.0
    - 'yscale': 1.0e4
    - 'yoffset': 0.0
    - 'width': 2.0
    - 'marks': 1
    - 'msize': 1.0
    - 'titleb', 'titlel', 'titlet', 'titler'
    - 'titles': 1
    """
    # --- Dictionary specifying plot defaults
    pldef = {'xscale': 100., 'xoffset': 0.0, 'yscale': 1.0e4, 'yoffset': 0.0,
             'width': 2.0, 'marks': 1, 'msize': 1.0,
             'titleb': "X (cm)", 'titlel': "B (Gauss)",
             'titlet': "B-field of Bgrd element vs. x", 'titler': "", 'titles': 1}
    if plot_title <> "": pldef['titlet'] = plot_title
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    for key in pldef.keys():    exec key+"=pldef['"+key+"']"
    #
    iiz = nint(zplot/top.bgrddz[nspec]) + nint(top.bgrdnz/2)
    iix = nint(xplot/top.bgrddx[nspec]) + nint(top.bgrdnx/2)
    #
    # --- Calculate x positions in grid to use for plotting fields
    ygrd = iota(-top.bgrdny/2, top.bgrdny/2)*top.bgrddy[nspec]
    #
    pldefault(width=width, marks=marks, msize=msize)
    plg( yoffset+yscale*top.bgrdbx[iix,:,iiz,nspec],
                                xoffset+xscale*ygrd[:], marker="x", type="solid")
    plg( yoffset+yscale*top.bgrdby[iix,:,iiz,nspec],
                                xoffset+xscale*ygrd[:], marker="y", type="dot")
    plg( yoffset+yscale*top.bgrdbz[iix,:,iiz,nspec],
                                xoffset+xscale*ygrd[:], marker="z", type="dash")
    if titles:  ptitles(titlet, titleb, titlel, titler)
    #
    fma()

#===========

def plot_bofz(nspec=0, xplot=0, yplot=0, plot_title="", kwdict={},  **kw):
    """ plot_bofz(nspec=0, xplot=0, yplot=0, plot_title="", kwdict={},  **kw)
    Plots bgrd field of element species 'nspec' vs. axis // z
xplot, yplot define coords of desired plot axis
    - 'xscale': 100.
    - 'xoffset': 0.0
    - 'yscale': 1.0e4
    - 'yoffset': 0.0
    - 'width': 2.0
    - 'marks': 1
    - 'msize': 1.0
    - 'titleb', 'titlel', 'titlet', 'titler'
    - 'titles': 1
    """
    # --- Dictionary specifying plot defaults
    pldef = {'xscale': 100., 'xoffset': 0.0, 'yscale': 1.0e4, 'yoffset': 0.0,
             'width': 2.0, 'marks': 1, 'msize': 1.0,
             'titleb': "Z (cm)", 'titlel': "B (Gauss)",
             'titlet': "B-field of Bgrd element vs. z", 'titler': "", 'titles': 1}
    if plot_title <> "": pldef['titlet'] = plot_title
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    for key in pldef.keys():    exec key+"=pldef['"+key+"']"
    #
    iix = nint(xplot/top.bgrddx[nspec]) + nint(top.bgrdnx/2)
    iiy = nint(yplot/top.bgrddy[nspec]) + nint(top.bgrdny/2)
    #
    # --- Calculate x positions in grid to use for plotting fields
    zgrd = iota(-top.bgrdnz/2, top.bgrdnz/2)*top.bgrddz[nspec]
    #
    pldefault(width=width, marks=marks, msize=msize)
    plg( yoffset+yscale*top.bgrdbx[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="x", type="solid")
    plg( yoffset+yscale*top.bgrdby[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="y", type="dot")
    plg( yoffset+yscale*top.bgrdbz[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="z", type="dash")
    if titles:  ptitles(titlet, titleb, titlel, titler)
    #
    fma()

#===========

def comp_bofz(nspec=0, xplot=0, yplot=0, plot_title="", kwdict={},  **kw):
    """ comp_bofz(nspec=0, xplot=0, yplot=0, plot_title="", kwdict={},  **kw)
    Plots bgrd field of element species 'nspec' vs. axis // z
xplot, yplot define coords of desired plot axis
    - 'xscale': 100.
    - 'xoffset': 0.0
    - 'yscale': 1.0e4
    - 'yoffset': 0.0
    - 'width': 2.0
    - 'marks': 1
    - 'msize': 1.0
    - 'titleb', 'titlel', 'titlet', 'titler'
    - 'titles': 1
    """
    # --- Dictionary specifying plot defaults
    pldef = {'xscale': 100., 'xoffset': 0.0, 'yscale': 1.0e4, 'yoffset': 0.0,
             'width': 2.0, 'marks': 1, 'msize': 1.0,
             'titleb': "Z (cm)", 'titlel': "B (Gauss)",
             'titlet': "B-field of Bgrd element vs. z", 'titler': "", 'titles': 1}
    if plot_title <> "": pldef['titlet'] = plot_title
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    for key in pldef.keys():    exec key+"=pldef['"+key+"']"
    #
    iix = nint(xplot/top.bgrddx[nspec]) + nint(top.bgrdnx/2)
    iiy = nint(yplot/top.bgrddy[nspec]) + nint(top.bgrdny/2)
    #
    # --- Calculate x positions in grid to use for plotting fields
    zgrd = iota(-top.bgrdnz/2, top.bgrdnz/2)*top.bgrddz[nspec]
    #
    pldefault(width=width, marks=marks, msize=msize)
    plg( yoffset+yscale*top.bgrdbx[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="x", type="solid")
    plg( yoffset+yscale*top.bgrdby[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="y", type="dot")
    plg( yoffset+yscale*top.bgrdbz[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="z", type="dash")
    if titles:  ptitles(titlet, titleb, titlel, titler)
	#
    plg( yoffset+yscale*solbx[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="x", type="solid")
    plg( yoffset+yscale*solby[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="y", type="dot")
    plg( yoffset+yscale*solbz[iix,iiy,:,nspec],
                                xoffset+xscale*zgrd[:], marker="z", type="dash")
    fma()

#===========


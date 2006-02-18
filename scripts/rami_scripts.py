""" Module rami_scripts.py
# by: Rami A. Kishek    (some functions based on Grote's scripts)
# created: Aug. 30, 2000
#
#	Last Modified: 5/13/2005
#
# Contains many convenience functions used in Rami's Python Input decks
# and other modules:
#
# ==================
# CONVENIENCE FUNCS:
#   calc_dedr	... Calculates dedr based on const beam size, perv, and emit
#   change_t_step 	... Changes time step in middle of run
# * change_pipe_rad ... Changes pipe radius in middle of run
#   apply_aper      ... Applies an aperture to beam
#   accelerate_beam	... applies zero-width acceleration gap
#   make_nice_output 	... Prints nice one-liners
#   rmarray         ... Creates array of Gaussian random #s same size as x
#   tableread       ... Reads a table of numeric data from formatted txt file
#   tablewrite      ... Writes an array of numeric data into formatted txt file
#
# MOMENT CALC + SAVING:
# * init_cmom   ... Initializes history z-moments needed for calc_mom()
# * calc_mom	... Calculates additional moments to be saved
# * seek_name   ... Function to assist in finding WARP variable names for saving
#   save_data	... Saves selected parameters for future plotting .pdb
#								(compatible w/ old BASIS create / restore)
#   save_data3  ... Saves z-moment mesh info in 3D Lab frame simus
#   save_long   ... Saves more variables than save_data (those calculated in calc_mom)
#
# PLOTTING:
#   gen_plot    ... generic history plotting routine with strobing, defaults to vzbar
#   plot_env	... Plot 2*rms beam enveloppe
#   plot_emit	... Plot 4*rms beam Emittance
#   plot_cent	... Plot beam centroid in x,y
#   plot_vz 	... Plot vz,thermal
#   plot_de     ... Plot ENERGY Spread
#   plot_vzbar  ... Plot ave vz
#   plot_ekin   ... Plot ave ENERGY
#   plot_np		... Plot np as function of t
#   plot_curr	... Plot beam current as function of t
#   plot_linechg	... Plot Line Charge as function of t
#   plot_rot	... Plot beam rotation angle
#   plot_remit  ... Plot generalized Emittances with rotations, Eg, Eh
#   plot_disp	... Plot dispersion function vz * < x * (uzp-vz) > / < (uzp-vz)^2 >
#   plot_demit  ... Plot generalized Emittance with dispersion, Ed, Ed2
#   plot_comp	... Compare several runs
#   comp_exp    ... Compares Current Run to Experimental Output of UMERBeam.m
#
# HISTOGRAMMING:
#   rhistogram	... Simple histogramming function using numeric package
#   plot_hist   ... Uses rhistogram to plot a histogram
#
# ==================
# OPTIONAL:
#   save_multi  ... Save multiple runs with multiple random seeds
# ==================

MODS #####################################

11/19/04 - fixed dispersion functions
5/11/05  - added relative plotting option
5/13/05  - partially fixed hepsd calc expressions
"""
#===========
# GLOBALS + INITIALIZATION

yes = 1; no = 0
aper_dist = 0.0
sfact = 1.0
lcalc_mom = no
import sys, __main__
from warp import *
from histplots import *

rami_scripts_version = "$Id: rami_scripts.py,v 1.6 2006/02/18 01:57:00 dave Exp $"
def rami_scriptsdoc():
  import rami_scripts
  print rami_scripts.__doc__

# -- Default Variables to be saved
def_vars1 = ["hpnum", "hvzbar", "hvzrms",
             "hxrms", "hyrms", "hepsx", "hepsy", "hxbar", "hybar"]
def_vars2 = ["zscale", "nx", "ny", "xmmax", "ymmax"]
def_vars =  def_vars1 + def_vars2
add_vars = ["hepsg", "hepsh"]                   # -- saved if _long
add_hmom = ["hepsd", "hdisp", "hrot"]           # -- calculated using calc_mom

# -- Moments required to calculate additional moments:
req_mom  = ["vzbar", "vzrms", "xvzbar", "vxvzbar", "xsqbar", "xpsqbar", "xxpbar", "epsnx",
	        "xybar", "xbar", "ybar", "xrms", "yrms", "vxbar", "vxrms"]  #"xvxbar", "vxsqbar"

#===========

def swhere( cond, x, y ):
    """ swhere( cond, x, y )
        Version of Numpy's 'where' that returns a scalar
    """
    return where(cond, x, y)[0]

#===========

def calc_dedr(a0=None, emit=None, ibeam=None, ekin=None):
	"""top.dedr = calc_dedr(a0=None, emit=None, ibeam=None, ekin=None)
		Parameters default to those in TOP
	Given beam parameters, calculates required dedr to keep
	beam at a constant size of a0
	"""
	if a0==None: 	a0=top.a0
	if emit==None:	emit=top.emit
	if ibeam==None:	ibeam=top.ibeam
	if ekin ==None:	ekin=top.ekin
	#
	gen_perv = 1.515e4*abs(ibeam)/(ekin**1.5)
	kappa = (gen_perv/(a0**2)) + ((emit**2)/(a0**4))
	return  -kappa*2.*ekin/(top.zion)

#===========

def accelerate_beam(accelv):
    """ accelerate_beam(accelv)
    	where accelv is the additional VOLTAGE
    # --- Apply acceleration to beam and beam frame.
    # --- This replaces the accl elements and models a zero-length gap.
    """
    accelv2 = abs(2.*top.sq[0]/top.sm[0])*accelv
    top.vbeamfrm = sqrt(top.vbeamfrm**2 + accelv2)
    top.vbeam = top.vbeamfrm
    top.uzp = sqrt(top.uzp**2 + accelv2)

#===========

def change_t_step(step_size):
    """ change_t_step(step_size)
  	# --- Change the step size
    """
    wxy.ds = step_size
    top.dt = wxy.ds/top.vbeam
    top.pid[:,wxy.dtpid-1] = wxy.ds/where( greater(top.uzp,0.), top.uzp, top.vbeam )

#===========

def change_pipe_rad(pipe_rad):
    """ change_pipe_rad(pipe_rad)
  	# --- Change radius at which particles are removed
  	# --- Change capacity matrix to new pipe_rad
    """
    top.prwallz = top.prwall = pipe_rad
    symm_fact = swhere( w3d.l4symtry, pi/2.0,
                        swhere(w3d.l2symtry, pi, 2.0*pi) )
    fxy.xcond[0:fxy.ncxy] = pipe_rad*cos(symm_fact*iota(fxy.ncxy)/fxy.ncxy)
    fxy.ycond[0:fxy.ncxy] = pipe_rad*sin(symm_fact*iota(fxy.ncxy)/fxy.ncxy)
    fxy.vcond[0:fxy.ncxy] = 0.e0            # Put pipe at ground
    fieldsol(0)                           	# recalculate matrix and fields
    w3d.phi[:,:,-1] = w3d.phi[:,:,1]		# Bug fix
    w3d.phi[:,:, 0] = w3d.phi[:,:,1]	

#===========

def apply_aper(aper_rad, aper_thick, naper_thick):
    """ apply_aper(aper_rad, aper_thick, naper_thick)
        Applies and aperture of radius 'aper_rad' and thickness 'aper_thick'
        and steps 'naper_thick' steps through it.
    """
    pr1 = top.prwall
    dzz = wxy.ds
    change_pipe_rad(aper_rad)
    change_t_step(aper_thick/naper_thick)
    step(naper_thick)
    change_t_step(dzz)
    change_pipe_rad(pr1)


#===========

def make_nice_output(nwin=0, *files):
    """ make_nice_output(nwin=0, *files)
	# --- Generate nice ouput to the terminal and any extra files """

    oneliner = "it = %5d zbeam = %9.5f 2*yrms = %7.2f xbar = %6.2f emitx = %7.2f nplive = %8d\n" %\
                (top.it,top.zbeam-aper_dist, 2.0*top.yrms[nwin]*1.e3, top.xbar[nwin]*1.e3,
                 top.epsx[nwin]*1.e6, top.nplive)
    files = (sys.stdout,)+files         # Write to terminal
    for fid in files:
        fid.write(oneliner)


#===========

def rmarray(x):
    """ rmarray(x)
    # --- Crude approximation to a Gaussian with mean 0, standard deviation 1
    # --- Advantage is, it cuts off smoothly at 3 root 2.
    # --- A similar one (used in some older AF codes) adds 12 variates,
    # --- subtracts 6.
    """
    return sqrt(2.)*(ranf(x)+ranf(x)+ranf(x)+ranf(x)+ranf(x)+ranf(x)-3.)


#===========

def init_cmom(l3d=yes, lhist=yes, lcalc_mom=yes):
    """ init_cmom(l3d=yes, lhist=yes, lcalc_mom=yes)
    Invoke at beginning of run if 3d lab frame and want z-moment history.
    Saves histories of proper moments needed for additional moment calculation
    using calc_mom().
    """
    if lhist:
        if l3d:
            top.lhvzofz = true
            for mom in def_vars1[1:]:
                exec "top.l"+mom+"z = true"
            if lcalc_mom:
                for mom in req_mom[0:]:
                    exec "top.lh"+mom+"z = true"
                for mom in add_vars:
                    exec "top.l"+mom+"z = true"
#         elif lcalc_mom:
#             for mom in req_mom[0:]:
#                 exec "top.lh"+mom+" = true"


#===========

def calc_mom(l3d=no, lhist=no):
    """ calc_mom(l3d=no, lhist=no)
    # To be invoked at end of run; l3d=yes for zmoments;
    #   lhist=yes for history of z-mom (BE SURE to INVOKE init_cmom() at start of run)
	# Calculates:
	#		hepsd = emittance conserved in dispersion
	#		hdisp  = Dispersion moment <x*delta>/<delta^2>
	#       hepsg, hepsh = Generalized emits w/ QROT
	#       hrot   = Beam rotation angle
	# Note: ignores centroids
    """
    global hrot, hdisp, hepsd, hepsd2
	# --- Select appropriate variables:
    m = {}; pre = 'h'; suf = ''
    if l3d:
        suf = 'z'
        if lhist:   req_mom[0] = 'vzof'
        else: pre = ''
    for mom in req_mom:
        m[mom] = eval("top."+pre+mom+suf)
    if lhist and l3d:   m["vzbar"] = m["vzof"]

    # --- Rotation Moments
    hdxy = m["xybar"] - (m["xbar"]*m["ybar"])
    hrot = 0.5*arctan(2*hdxy/(m["xrms"]**2-m["yrms"]**2))

    # --- Dispersion moments
    hdisp    = m["vzbar"]*(m["xvzbar"]-(m["xbar"]*m["vzbar"]))/(m["vzrms"])**2
    print "ave Dispersion", ave(hdisp[0,:])

    # --- hepsnd

    # --- Code to correct for missing moments
    m["vxsqbar"] = (m["vxrms"]**2) + (m["vxbar"]**2)
    if "hxvxbar" in dir(top):   m["xvxbar"] = top.hxvxbar   #!! NEEDS FIXING
    elif "hxvxbar" in dir(__main__): m["xvxbar"] = __main__.hxvxbar
    else:
        print "\n\n!! You are making an inaccurate assumption about 'xvxbar' !!\n\n"
        m["xvxbar"] = m["xxpbar"]*m["vzbar"]

    hDd   = (m["xvzbar"] - (m["xbar"]*m["vzbar"]))/(m["vzrms"])
    hDpd  = (m["vxvzbar"]-(m["vxbar"]*m["vzbar"]))/(m["vzrms"])

    hepsd = sqrt(((m["xrms"]**2)-(hDd**2)) * (m["vxrms"]**2-(hDpd**2)) -
                 (m["xvxbar"]-m["xbar"]*m["vxbar"]-hDd*hDpd)**2 )*4.0/m["vzbar"]


# ==================

def seek_name(var, l3d=no, lhist=no):
    """ seek_name(var, l3d=no, lhist=no)
        Given an output var name, find actual variable name
    """
    pre = 'top.h'; suf = ''
    if l3d:
        suf = 'z'
        if lhist:   pre = 'transpose(top.h'; suf = 'z)'
        else:       pre = 'top.'
    if (lhist and l3d):
        if var == "hvzbar":
            try:                return "transpose(top.hvzbarz)"
            except NameError:   return "transpose(top.hvzofz)"
        if var == "hzbeam":     return "top.hzbeam"
        if var == "hpnum":      return "top.pnumz"
        if var == "hlinechg":   return "transpose(top.hlinechg)"
        if var == "hcurr":      return "transpose(top.hcurrz)"
    if var in top.varlist("dump"):  return pre+var[1:]+suf
    if var in w3d.varlist("dump"):  return "w3d."+var
    if var in globals().keys():     return var
    if var == "zscale":
        if l3d: return "top.zmntmesh"
        else:   return "top.hzbeam"
    print var, "not saved!"
    return ""

# ==================

def save_data(crun="0", vars=def_vars, l3d=no, lhist=no, llw=no):
    """ save_data(crun="0", vars=def_vars, l3d=no, lhist=no, llw=no)
    #   Uses PW.PW routine to save vars in pdb format.
    #       default vars stored in list 'def_vars'
    #   Read output using
    #       restore("filename.pdb")
    #   or
    #       out = PR.PR("filename.pdb")
    #       out.???
    """
    runid = arraytostr(top.runid)
    outfile = PW.PW("data."+runid+crun+".pdb")
    #
    for vname in vars:
        __main__.__dict__[vname+runid] = eval(seek_name(vname, l3d, lhist))
        outfile.write( vname+runid, eval(vname+runid, __main__.__dict__) )
    if llw:     # Lab Windows used
        for vname in ['hlinechg', 'hcurr', 'htime']+def_vars1:
            __main__.__dict__[vname[1:]+'lw'+runid] = eval('transpose(top.'+vname[1:]+'lw)')
            outfile.write( vname[1:]+'lw'+runid, eval(vname[1:]+'lw'+runid, __main__.__dict__) )
#        __main__.__dict__["timelw"+runid] = eval('transpose(top.timelw)')
#        outfile.write( 'timelw'+runid, eval('timelw'+runid, __main__.__dict__) )
        __main__.__dict__["ESPREADlw"+runid] = transpose(eval(
            '(emass/echarge)*(top.vzbarlw*top.vzrmslw)' ))
        outfile.write( 'ESPREADlw'+runid, eval('ESPREADlw'+runid, __main__.__dict__) )
        __main__.__dict__["ekinlw"+runid] = transpose(eval(
            '0.5*(emass/echarge)*(top.vzbarlw**2)' ))
        outfile.write( 'ekinlw'+runid, eval('ekinlw'+runid, __main__.__dict__) )
    outfile.close()

# ==================

def save_data3(crun="0", vars=def_vars, lhist=no, llw=no):
    """ save_data3(crun="0", vars=def_vars, lhist=no, llw=no)
        Same as save_data with l3d=yes, for 3D runs with zmoment mesh
    """
    save_data(crun=crun, vars=vars, l3d=yes, lhist=lhist, llw=llw)

# ==================

def save_long(crun="0", vars=def_vars+add_vars+add_hmom, l3d=no, lhist=no, llw=no):
    """ save_long(crun="0", vars=def_vars+add_vars+add_hmom, l3d=no, lhist=no, llw=no)
	#   Uses PW.PW routine to save vars in pdb format.
	#       default vars stored in variables 'def_vars' + 'add_vars'
	#	Read output using
	#       restore("filename.pdb")
	#   or
	#		out = PR.PR("filename.pdb")
	#		out.???
    """
    save_data(crun, vars, l3d, lhist, llw)

#===========

def gen_plot(vars=(), xaxis="", runid=None, kwdict={},  **kw):
    """ gen_plot(vars=(), xaxis="", runid=None,  kwdict={}, **kw)
	    Generic plotter of any number of identical variables, with strobing
	    Understands any of Grote's arguments [hpdoc() for help] +  the following:

	    vars = tuple of names of variables plotted (defaults to "vzbar")
	    xaxis = name of x-axis variable (defaults to "zscale")
	    runid = runid of run to be plotted.  (None defaults to current run)
	            "" plots vars with no runid attached to their name
	    nwin = window to be plotted, (def None if array 1-D)
	    begin = beginning index (0)
	    strobe = # of points to be strobed (1)
    """
    # --- Dictionary specifying plot defaults
    pldef = {'nwin': None, 'begin':0, 'end': -1, 'strobe':1,
             'type': "solid", 'color': 'fg', 'width': 2.0,
             'marks': 0, 'marker': None, 'msize': 1.0,
             'logplot': no, 'lrel': 0,
             'xmin': 'e', 'xmax': 'e', 'ymin': 'e', 'ymax': 'e',
             'xscale': 1.0, 'xoffset': 0.0, 'yscale': 1.0, 'yoffset': 0.0,
             'titleb': "S (m)", 'titlel': "", 'titlet': "", 'titler': "", 'titles': 1}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    for key in pldef.keys():    exec key+"=pldef['"+key+"']"
    #
    if runid is None:
        runid = arraytostr(top.runid)
        titlet = titlet+' '+runid
    if xaxis == "": xaxis = "zscale"
    if vars == ():  vars = ("hvzbar",)
    #
    if (len(eval(xaxis+runid, __main__.__dict__).shape) > 1) and (nwin is not None) and (
        eval(xaxis+runid, __main__.__dict__).shape[0] > 1):
        absc = xoffset + xscale*eval(xaxis+runid, __main__.__dict__)[nwin,  begin:end:strobe]
    else:
        absc = xoffset + xscale*eval(xaxis+runid, __main__.__dict__)[...,  begin:end:strobe]
    for vname in vars:
        if 'x' in vname: marker = 'x'
        elif 'y' in vname: marker = 'y'
        elif 'g' in vname:  marker = 'g'
        elif 'h' in vname[1:]:  marker = 'h'
        elif '2' in vname:  marker = '2'
        if nwin is None:
            ordinate = eval(vname+runid, __main__.__dict__)[ ..., begin:end:strobe]
        else:
            ordinate = eval(vname+runid, __main__.__dict__)[nwin, begin:end:strobe]
        oord = yoffset + yscale*ordinate
        if lrel:    # Plot ordinate relative to initial value
            oord = oord/oord[..., 0]
        if logplot: oord = log(oord)
        plg( oord, absc, color=color, width=width, type=type,
             marks=marks, marker=marker, msize=msize )
    limits(xmin, xmax, ymin, ymax)
    if titles:  ptitles(titlet, titleb, titlel, titler)

#===========

def plot_env(runid=None,  kwdict={}, **kw):
    """ plot_env(runid=None, kwdict={}, **kw)
        Plots ENVS in x and y, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 2.e3,
             'titlel': "X,Y 2*RMS Envs (mm)", 'titlet': "Beam Envelope from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    try:    gen_plot(("hxrms", "hyrms"), "zscale", runid, pldef)
    except NameError:
        gen_plot(("xenv", "yenv"), "zscale", runid, pldef, yscale=1.e3)

#===========

def plot_xenv(runid=None,  kwdict={}, **kw):
    """ plot_xenv(runid=None, kwdict={}, **kw)
        Plots ENVS in x, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 2.e3,
             'titlel': "X,Y 2*RMS Envs (mm)", 'titlet': "Beam Envelope from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    try:    gen_plot(("hxrms",), "zscale", runid, pldef)
    except NameError:
        gen_plot(("xenv",), "zscale", runid, pldef, yscale=1.e3)

#===========

def plot_emit(runid=None,  kwdict={}, **kw):
    """ plot_emit(runid=None, kwdict={}, **kw)
        Plots EMITS in x and y, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 1.e6,
             'titlel': "X,Y 4*RMS Emits (mm-mr)", 'titlet': "Beam Emmittance from"}
    pldef.update(kwdict);    pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hepsx", "hepsy"), "zscale", runid, pldef)
    except NameError:   gen_plot(("xemit", "yemit"), "zscale", runid, pldef)

#===========

def plot_xemit(runid=None,  kwdict={}, **kw):
    """ plot_xemit(runid=None, kwdict={}, **kw)
        Plots EMITS in x, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 1.e6,
             'titlel': "X,Y 4*RMS Emits (mm-mr)", 'titlet': "Beam Emmittance from"}
    pldef.update(kwdict);    pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hepsx",), "zscale", runid, pldef)
    except NameError:   gen_plot(("xemit",), "zscale", runid, pldef)

#===========

def plot_nemit(runid=None,  kwdict={}, **kw):
    """ plot_nemit(runid=None, kwdict={}, **kw)
        Plots NORMALIZED EMITS in x and y, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 1.e6,
             'titlel': "X,Y 4*RMS Emits (mm-mr)", 'titlet': "Beam Emmittance from"}
    pldef.update(kwdict);    pldef.update(kw)    # Override defaults & import new params
    gen_plot(("hepsx"+runid+"*hvzbar"+runid+"/clight",
              "hepsy"+runid+"*hvzbar"+runid+"/clight"),
              "zscale"+runid, runid="", kwdict=pldef)

#===========

def plot_cent(runid=None,  kwdict={}, **kw):
    """ plot_cent(runid=None, kwdict={}, **kw):
        Plots CENTROIDS in x and y, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 1.e3,
             'titlel': "<X>, <Y> (mm)", 'titlet': "Beam Centroid from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hxbar", "hybar"), "zscale", runid, pldef)
    except NameError:   gen_plot(("xbar", "ybar"), "zscale", runid, pldef)

#===========

def plot_vz(runid=None,  kwdict={}, **kw):
    """ plot_vz(runid=None, kwdict={}, **kw)
        Plots VZ,th, using gen_plot
    """
    pldef = {'nwin': 0, 'titlel': "Vz,th (m/s)", 'titlet': "Vz,th from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hvzrms", ), "zscale", runid, pldef)
    except NameError:   gen_plot(("vzrms", ), "zscale", runid, pldef)

#===========

def plot_curr(runid=None,  kwdict={}, **kw):
    """ plot_curr(runid=None, kwdict={}, **kw)
        Plots Beam Current, using gen_plot
    """
    pldef = {'nwin': 0, 'yscale': 1.e3,
             'titlel': "I (mA)", 'titlet': "Beam Current from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    gen_plot(("hcurr", ), "zscale", runid, pldef)

#===========

def plot_linechg(runid=None,  kwdict={}, **kw):
    """ plot_linechg(runid=None, kwdict={}, **kw)
        Plots Line Charge, using gen_plot
    """
    pldef = {'nwin': 0, 'titlel': "lambda (C/m)", 'titlet': "Line charge from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    gen_plot(("hlinechg", ), "zscale", runid, pldef)

#===========

def plot_vzbar(runid=None,  kwdict={}, **kw):
    """ plot_vzbar(runid=None, kwdict={}, **kw)
        Plots ave VZ, using gen_plot
    """
    pldef = {'nwin': 0, 'titlel': "Vz-BAR (m/s)", 'titlet': "Vz-BAR from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hvzbar", ), "zscale", runid, pldef)
    except NameError:   gen_plot(("vzbar", ), "zscale", runid, pldef)

#===========

def plot_de(runid=None,  kwdict={}, **kw):
    """ plot_de(runid=None, kwdict={}, **kw)
        Plots Energy Spread Delta-E, using gen_plot
    """
    if runid is None:   runid = arraytostr(top.runid)
    try:
        vzb = eval("hvzbar"+runid, __main__.__dict__)
        vzth = eval("hvzrms"+runid, __main__.__dict__)
    except NameError:
        vzb = eval("vzbar"+runid, __main__.__dict__)
        vzth = eval("vzrms"+runid, __main__.__dict__)
    __main__.__dict__["ESPREAD"+runid] = (emass/echarge)*(vzb*vzth)
    pldef = {'nwin': 0, 'titlel': "Delta-E (eV)", 'titlet': "Delta-E from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    gen_plot(("ESPREAD", ), "zscale", runid, pldef)

#===========

def plot_ekin(runid=None,  kwdict={}, **kw):
    """ plot_ekin(runid=None, kwdict={}, **kw)
        Plots Beam Kinetic Energy using gen_plot
    """
    if runid is None:   runid = arraytostr(top.runid)
    try:                vzb = eval("hvzbar"+runid, __main__.__dict__)
    except NameError:   vzb = eval("vzbar"+runid, __main__.__dict__)
    __main__.__dict__["ekin"+runid] = 0.5*(emass/echarge)*(vzb**2)
    pldef = {'nwin': 0, 'titlel': "Ekin (eV)", 'titlet': "E_kin from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    gen_plot(("ekin", ), "zscale", runid, pldef)

#===========

def plot_np(runid=None,  kwdict={}, **kw):
    """ plot_np(runid=None, kwdict={}, **kw)
        Plots NUM Particles, using gen_plot
    """
    pldef = {'nwin': 0, 'titlel': "NP", 'titlet': "Number Particles from"}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hpnum", ), "zscale", runid, pldef)
    except NameError:   gen_plot(("npleft", ), "zscale", runid, pldef)

#===========

def plot_rot(runid=None,  kwdict={}, **kw):
    """ plot_rot(runid=None, kwdict={}, **kw)
        Plots Beam Rotation Angle, using gen_plot
    """
    pldef = {'nwin': 0, 'yscale': 180./pi,
             'titlel': "Rot. Ang. (deg.)", 'titlet': "Beam Rotation Angle from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hrot", ), "zscale", runid, pldef)
    except NameError:   gen_plot(("rot", ), "zscale", runid, pldef, nwin=None)

#===========

def plot_remit(runid=None,  kwdict={}, **kw):
    """ plot_remit(runid=None, kwdict={}, **kw)
        Plots GEN EMITS w/ Rotation, Eg, Eh, using gen_plot
    """
    pldef = {'nwin': 0, 'marks': 1, 'yscale': 1.e6,
             'titlel': "4*RMS Emits (mm-mr)", 'titlet': "E_g, E_h from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    try:                gen_plot(("hepsg", "hepsh"), "zscale", runid, pldef)
    except NameError:   gen_plot(("remit", "hemit"), "zscale", runid, pldef, nwin=None)

#===========

def plot_disp(runid=None,  kwdict={}, **kw):
    """ plot_disp(runid=None, kwdict={}, **kw)
        Plots Marco's dispersion function, using gen_plot
    """
    pldef = {'nwin': 0,
             'titlel': "<x*delta>/<delta**2>  (m)", 'titlet': "Dispersion Function from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    gen_plot(("hdisp", ), "zscale", runid, pldef)

#===========

def plot_demit(runid=None, kwdict={}, **kw):
    """ plot_demit(runid=None, kwdict={}, **kw)
        Plots GEN EMIT w/ Dispersion in x, Ed (+Ed2), using gen_plot
    """
    pldef = {'nwin': 0, 'yscale': 1.e6,
             'titlel': "4*RMS Emits (mm-mr)", 'titlet': "E_d from"}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    gen_plot(("hepsd",), "zscale", runid, pldef)

#===========

def plot_comp(runs={}, plots=("env", "emit", "cent"), kwdict={}, **kw):
    """ plot_comp(runs={}, plots=("env", "emit", "cent"), kwdict={}, **kw)
	    Generates comparison plots between any number of simulations using gen_plot
            runs  = dictionary relating strings of runids+crun to dictionary of
                    properties
                    eg., {"abcd0": {'type': "solid", 'yscale': 2.0, ...}, ...}
            plots = tuple of strings specifying plots be compared
    """
    if  len(runs) < 2: runs[arraytostr(top.runid)+"0"] = {'type': "solid"}
    pldef = {'titlet': ''}
    for run in runs.keys():
        restore("data."+run+".pdb")
        pldef['titlet'] = pldef['titlet']+ run+" ("+runs[run]['type']+"); "
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    #
    for plot in plots:
        for run in runs.keys():
            temp = pldef.copy()
            temp.update(runs[run])
            exec "plot_"+plot+"('"+run[:-1]+"', kwdict=temp)"
        fma()

#===========

def tablewrite(outfile, tdata, sep='\t', fmt="%12.8f"):
    """ tablewrite(outfile, tdata, sep='\t', fmt="%12.8f")
        Writes a table (array) of numeric data into a formatted txt file
        with file HANDLE 'outfile', placing the numbers in columns separated
        by separator 'sep'. 'fmt' is an optional formatting string.
    """
    fmtst = fmt+sep
    for line in tdata:
        for el in line:
            outfile.write(fmtst % el)
        outfile.write('\n')


#=========== This doesn't always work!!!
#def tableread(expfile, sep='\t'):
#    """ tdata = tableread(expfile, sep='\t')
#        Reads a table of numeric data in formatted txt file 'expfile',
#        which contains the numbers in columns separated by separator 'sep'.
#    """
#    import string
#    tdata = open(expfile, 'r').readlines()
#    tdata = map(lambda s : string.split(string.strip(s), sep), tdata)
#    return array(map(lambda s:map(float,s),tdata))

def tableread(expfile, sep='\t', R1=0, C1=0):
    """ tdata = tableread(expfile, sep='\t', R1=0, C1=0)
        Reads a table of numeric data in formatted txt file 'expfile',
        which contains the numbers in columns separated by separator 'sep'.
        R1 and C1 are the first row and column to begin download
        R2 and C2 are the last row and column
    """
    import string
    tdata = open(expfile, 'r').readlines()
    for ii in range(0,R1): del tdata[ii]
    for line in tdata:
        lineno = tdata.index(line)
        line = string.split(string.strip(line), sep)
        for ii in range(0,C1): del line[ii]
        tdata[lineno] = line
        for ii in range(0, len(line)):
            el = line[ii]
            if el == '': el = 0.0
            tdata[lineno][ii] = float(el)
    return array(tdata, 'd')


#===========

def comp_exp(expfile, sep='\t', run=None,  kwdict={}, **kw):
    """ comp_exp(expfile, sep='\t', run=None, kwdict={}, **kw)
        Compares Current Run to Experimental Output of UMERBeam.m
    Plots centroid, envelope, and rotation angle, using gen_plot
    Experimental Data supplied using txt file 'expfile', which contains
    the numbers in columns separated by separator 'sep'.
    'run' = 'runid'+'crun'
    """
    runid = None
    if run is not None:
        restore("data."+run+".pdb")
        runid = run[:-1]
    pldef = {'nwin': 0, 'marks': 0, 'msize': 1.25, 'xscale': 1.0}
    pldef.update(kwdict); pldef.update(kw)    # Override defaults & import new params
    if pldef['xscale'] == 100.0: pldef['titleb'] = "S (cm)"
    if pldef['xscale'] == 1000.0: pldef['titleb'] = "S (mm)"

    # --- Load Experimental Data
    expdata = tableread(expfile, sep)
    z = 1.e-2*pldef["xscale"]*expdata[:,2]
    a = expdata[:,3]
    b = expdata[:,4]
    xc = expdata[:,5]
    yc = expdata[:,6]
    rot = expdata[:,7]

    # --- Centroid
    plot_cent(runid, kwdict=pldef)
    plg( xc, z, type="none", marker='\5', msize=pldef['msize'] )
    plg( yc, z, type="none", marker='\4', msize=pldef['msize']  )
    fma()

    # --- Envelope
    plot_env(runid, kwdict=pldef)
    plg( a, z, type="none", marker='\5', msize=pldef['msize'] )
    plg( b, z, type="none", marker='\4', msize=pldef['msize'] )
    fma()

    # --- Rotation Angle
    if lcalc_mom:
        plot_rot(runid, kwdict=pldef)
        plg( rot, z, type="none", marker='\4', msize=pldef['msize'] )
        fma()


#==========

def rhistogram(a, nbins=30.):
	"""rhistogram(a, nbins=30.)
		returns  two arrays of length nbins:
		    bins: containing the x axis
		    f: containing the distribution
	"""
	bot = min(a.flat); top = max(a.flat)
	bins = arange(bot, top, (top-bot)/nbins)
	n = searchsorted(sort(a), bins)
	n = concatenate([n, [len(a)]])
	return bins, (0.0+n[1:]-n[:-1])/len(a)

#===========

def plot_hist(dat, nbins=None, kwdict={}, **kw):
    """ plot_hist(dat, nbins=None, kwdict={}, **kw)
        Plots histogram of 'dat'
    """
    # -- Caluclate Histogram
    if nbins is None:   nbins = max( 10, nint(len(dat)/10.0) )
    bins, oord = rhistogram(dat, nbins)
    # --- Dictionary specifying plot defaults
    pldef = {'xscale': 1.0, 'yscale': 1.0, 'titlel': "# occurrences",
             'titlet': "", 'titleb': "", 'titler': "", 'titles': 1}
    pldef.update(kwdict);  pldef.update(kw)    # Override defaults & import new params
    # --- Actual plotting
    pla( pldef['yscale']*oord, pldef['xscale']*bins )
    if pldef['titles']:
        ptitles(pldef['titlet'],pldef['titleb'],pldef['titlel'],pldef['titler'])

#===========

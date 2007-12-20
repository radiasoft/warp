"""
Module postrack.py
by: Rami A. Kishek
Created: Nov 27, 2001

Last Modified: Feb. 1, 2002

Contains a streamlined improved version of Lyapunov Calculations
Contains functions for post-processing saved particle trajectories:
Can be imported anytime after a simulation is over, and can work
on saved data.

top.runid needs to be defined prior to importing (or reloading).

plot_moms()     ... Plot Moments for each species
plot_delay()    ... Generate set of delay coordinate plots
findindices()   ... Searches for indices in vector corresponding to low/up limits
calc_lyap()     ... Calculates Lyapunov exponent for each clump
save_lyap()     ... Saves Lyapunov exponent for each clump into pdb file
print_lyap()    ... Prints out Lyapunov exponent for each clump into txt file
plot_diverge()  ... Plot difference b/w particle trajectories and centroids
plot_traj()     ... Poincare' plots of particle trajectories for each clump
plot_spect()    ... Fourier spectrum of particle trajectory
calc_irrev()    ... Calculate Irreversibility factors (coefficients)
save_irrev()    ... Save Irreversibility factors (coefficients) into text file

"""

from rami_scripts import *
from Fitting import *
import __main__, os, sys

postrack_version = "$Id: postrack.py,v 1.4 2007/12/20 00:36:41 dave Exp $"
def postrackdoc():
  import postrack
  print postrack.__doc__

# --- Restore data files if postprocessing

runid = arraytostr(top.runid)
if "traj."+runid+"0.pdb" in os.listdir('.'): restore("traj."+runid+"0.pdb")
if "data."+runid+"0.pdb" in os.listdir('.'): restore("data."+runid+"0.pdb")
try:
    top.ns = eval("hepsx"+runid,__main__.__dict__).shape[0]
    n_steps  = eval("xp_"+runid, __main__.__dict__).shape[1]
    num_part = eval("xp_"+runid, __main__.__dict__).shape[2]
    dzz = eval("zscale"+runid,__main__.__dict__)[21]-eval("zscale"+runid,__main__.__dict__)[20]
except:
    top.ns = 6
    n_steps = 8640
    num_part = 20
    dzz = 0.02

# --- Define Constants and Parameters

colors = ["fg", "red", "blue", "green", "magenta", "cyan"]
ncolors = len(colors)
plots = ['x', 'y', 'xp', 'yp']

try:
    perv = abs(1.515e4*top.ibeam/top.ekin**1.5)
    k0x = sqrt(2.*perv/(top.a0*(top.a0+top.b0))+(top.emitx**2/top.a0**4))
    k0y = sqrt(2.*perv/(top.b0*(top.a0+top.b0))+(top.emity**2/top.b0**4))
    tune_depx = top.emitx/(k0x*(top.a0**2))
    tune_depy = top.emity/(k0y*(top.b0**2))
    lambda_x = 2*pi/(k0x*tune_depx)
    lambda_y = 2*pi/(k0y*tune_depy)
except:
    print "Numbers not right, tunes may be wrong and hence trajects may be wrong"
    k0x = k0y = sqrt(15.775)

# --- Allocate arrays for Lyapunov Exponents

lyexp_x = zeros([top.ns-1, 9],'d')
lyexp_y = zeros([top.ns-1, 9],'d')
init_zx  = ones([top.ns-1, 3],'l')  # Final Cutoffs
init_zy  = ones([top.ns-1, 3],'l')  # Initial Cutoffs
final_zx = ones([top.ns-1, 3],'l')  # Final Cutoffs
final_zy = ones([top.ns-1, 3],'l')  # Final Cutoffs

# --- Allocate Arrays for particle differences ("divergence")

eps_x = zeros([top.ns-1, n_steps, num_part-1],'d')
eps_y = zeros([top.ns-1, n_steps, num_part-1],'d')
eps_xp = zeros([top.ns-1, n_steps, num_part-1],'d')
eps_yp = zeros([top.ns-1, n_steps, num_part-1],'d')
eps_t = zeros([top.ns-1, n_steps, num_part-1],'d')

# --- Allocate Arrays for irreversibility factors

secmoms = ["xrms", "yrms", "epsx", "epsy", "epsg", "epsh"]

for mom in secmoms:
    __main__.__dict__["ir"+mom] = zeros([top.ns],'d')


#============================================================================

def plot_moms(plots = ("env", "emit", "remit", "cent", "vz", "np"), jspec=None, **kw):
    """ plot_moms(plots = ("env", "emit", "remit", "cent", "vz", "np"), jspec=None, **kw)
        Plot Moments for each species
        Plots over all species, except 0, by default, unless a range of species
        is provided via the parameter 'jspec'
    """
    if jspec is None: jspec = range(1, top.ns)
    if 'titlet' in kw.keys(): title = kw['titlet']
    else: title = ''
    for plot in plots:
        kw['titlet'] = title+' '+plot+" from "
        for spec in jspec:
            kw['nwin'] = spec
            kw['color'] = colors[spec % ncolors]
            apply( eval("plot_"+plot), (), kw)
        fma()

def plot_delay(lfull=0, jspec=None):
    """ plot_delay(lfull=0, jspec=None)

        Generate set of delay coordinate plots
        if 'lfull' true, plot delays for each clump, otherwise main beam only
        Plots over all species, except 0, by default, unless a range of species
        is provided via the parameter 'jspec'
    """
    strobex = nint(lambda_x/(wxy.ds*top.nhist))
    strobey = nint(lambda_y/(wxy.ds*top.nhist))
    absc_indx = slice(0,n_steps+1-strobex,strobex)
    oord_indx = slice(strobex,n_steps+1,strobex)
    absc_indy = slice(0,n_steps+1-strobey,strobey)
    oord_indy = slice(strobey,n_steps+1,strobey)

    if lfull: sp_last = top.ns
    else: sp_last = 1
    if jspec is None: jspec = range(0, sp_last)
    for spec in jspec:
        plg(eval("hepsx"+runid,__main__.__dict__)[spec, oord_indx],
            eval("hepsx"+runid,__main__.__dict__)[spec, absc_indx], marker='x')
        plg(eval("hepsy"+runid,__main__.__dict__)[spec, oord_indy],
            eval("hepsy"+runid,__main__.__dict__)[spec, absc_indy], marker='y')
        ptitles("Emittance Delay Plot"+`spec`+' '+colors[spec % ncolors]); fma()
        plg(eval("hxrms"+runid,__main__.__dict__)[spec, oord_indx],
            eval("hxrms"+runid,__main__.__dict__)[spec, absc_indx], marker='x')
        plg(eval("hyrms"+runid,__main__.__dict__)[spec, oord_indy],
            eval("hyrms"+runid,__main__.__dict__)[spec, absc_indy], marker='y')
        ptitles("Envelope Delay Plot for species "+`spec`+' '+colors[spec % ncolors]); fma()


# ----------------------

def calc_irrev():
    """ calc_irrev()

        Calculate irreversibility factors (coefficients) for each clump
    """
    for mom in secmoms:
        __main__.__dict__["ir"+mom][:] = eval('h'+mom+runid, __main__.__dict__
                    )[:, top.it]/eval('h'+mom+runid, __main__.__dict__)[:, 0]
        print mom, __main__.__dict__["ir"+mom]


def save_irrev(simlen=0.0):
    """ save_irrev(simlen=0.0)

        Save irreversibility factors (coefficients) in Text File
        Simlen is the length of simulation before reversal.
    """
    out = open("rev."+runid[0:3]+".txt", "a")
    out.write("Simulation Length = "+`simlen`+" * lambda_p"'\n\n')
    for mom in secmoms:
        out.write(mom+" = %s\n" % eval("ir"+mom, __main__.__dict__))
    out.write('\n\n')
    out.close()


# =========================

def print_lyap(outf=None, sep=' ', fmt="%12.8f"):
    """ print_lyap(outf=None, sep=' ', fmt="%12.8f")
        Prints out Lyapunov Exponents in tabular form into a text file:
        with filename "outf".  If no file provided prints to sys.stdout (default).
        The arrays are printed as a table with separators 'sep' and the
        numbers are formatted according to format string 'fmt'.
    """
    if outf is not None: out = open(outf, 'w')
    else: out = sys.stdout

    nmeth = init_zx.shape[1]
    vect = ones(2*array(init_zx.shape),'l')
    vect[0:top.ns-1, :nmeth] = init_zx
    vect[0:top.ns-1, nmeth:] = init_zy
    vect[top.ns-1:, :nmeth] = final_zx
    vect[top.ns-1:, nmeth:] = final_zy

    out.write(runid+'\n')
    out.write("\nStart X Indices for each Spec.         Y indices\n")
    tablewrite(out, vect[:top.ns-1], sep, "%8d")
    out.write("\nStart X Positions per Spec (m)         Y Positions\n")
    tablewrite(out, vect[:top.ns-1]*dzz, sep, "%8.3f")
    out.write("\nFinal X Indices for each Spec.         Y indices\n")
    tablewrite(out, vect[top.ns-1:], sep, "%8d")
    out.write("\nFinal X Positions per Spec (m)         Y Positions\n")
    tablewrite(out, vect[top.ns-1:]*dzz, sep, "%8.3f")

    out.write("\nLyapunov Exponents (Rows = species; columns = methods)\n")
    out.write("\n  X moms\n")
    tablewrite(out, lyexp_x, sep, fmt)
    out.write("\n")
    out.write("  Y moms\n")
    tablewrite(out, lyexp_y, sep, fmt)
    out.write("\n")

    if outf is not None: out.close()


# ----------------------

def findindices(vect, lowlim, hilim, sense="+"):
    """ (lowind, hiind) = findindices(vect, lowlim, hilim, sense="+")
    Searches for indices of vector 'vect' corresponding to the lower
and upper limits.  The relation between hilim and lowlim determine the
search criterion (e.g., '<' / '>').  If hilim ==lowlim, indices (0, len(vect))
are returned.  'sense' = '+' / '-' determines whether search is done
forward or backward.
    """
    # Determine condition rel for "while" statement
    if lowlim == hilim:  return (0, len(vect))
    if sense=="+":  # Forward Search
        lowind = 0
        if lowlim < hilim:  rel = '<'       # Growth
        else:               rel = '>'       # Decay
        exec("while vect[lowind]" + rel + " lowlim: lowind = lowind " + sense + " 1")
        hiind = lowind
        exec("while vect[hiind]" + rel + " hilim: hiind = hiind " + sense + " 1")
    else:           # Reverse Search
        hiind = len(vect)-1
        if lowlim < hilim:  rel = '>'
        else:               rel = '<'
        exec("while vect[hiind]" + rel + " hilim: hiind = hiind " + sense + " 1")
        lowind = hiind
        exec("while vect[lowind]" + rel + " lowlim: lowind = lowind " + sense + " 1")
    return (lowind, hiind)


def calc_lyap(meth="75%"):
    """ calc_lyap(meth="75%")

        Calculate Lyapunov exponents using various methods.

        The following criteria are used for determining the points that
        go into the calculation:
            - The beginning point is the time when the log of the moment reaches
            1/9 of the final cutoff.
            - The end point can be determined using different methods, depending on
            the variable 'meth':
            "NN%" picks the NN% level of the final value (e.g., 90%)

        The cutoff z_intercepts are determined either by searching forward 'f',
        reverse 'r', or both then averaging the two 'b'.
        The data is then fitted using a least squares fit to determine slope 'a',
        intercept 'b', and fit parameter 'L'.

        The output lyexp_x/y is given for each method:
            0: 'b' a: slope ~ Lyapunov Exponent (average)
            1: 'f' a: slope ~ Lyapunov Exponent (forward)
            2: 'r' a: slope ~ Lyapunov Exponent (reverse)
            3: 'b' L: Ave square of deviation, "goodness" of fit (average)
            4: 'f' L: Ave square of deviation, "goodness" of fit (forward)
            5: 'r' L: Ave square of deviation, "goodness" of fit (reverse)
            6: 'b' b: y-intercept (average)
            7: 'f' b: y-intercept (forward)
            8: 'r' b: y-intercept (reverse)
    """
    # --- Determine Cutoffs

    init_vsx = log(eval("hepsx"+runid,__main__.__dict__)[1:top.ns,0])
    init_vsy = log(eval("hepsy"+runid,__main__.__dict__)[1:top.ns,0])
    final_vsx = log(eval("hepsx"+runid,__main__.__dict__)[1:top.ns,-1])
    final_vsy = log(eval("hepsy"+runid,__main__.__dict__)[1:top.ns,-1])
    #
    if meth == "dev":
        pass
    elif meth[2] == '%':
        frac = float(meth[0:2])/100.0
        final_cfx = init_vsx + frac*(final_vsx-init_vsx)
        final_cfy = init_vsy + frac*(final_vsy-init_vsy)
    #
    init_cfx = init_vsx + (final_cfx-init_vsx)/9.0
    init_cfy = init_vsy + (final_cfy-init_vsy)/9.0

    for sp in range(0, top.ns-1):

    # Forward Search
        (init_zx[sp,1], final_zx[sp,1]) = findindices(
                log(eval("hepsx"+runid,__main__.__dict__)[sp+1]),
                init_cfx[sp], final_cfx[sp], "+")
        (init_zy[sp,1], final_zy[sp,1]) = findindices(
                log(eval("hepsy"+runid,__main__.__dict__)[sp+1]),
                init_cfy[sp], final_cfy[sp], "+")

    # Reverse Search
        (init_zx[sp,2], final_zx[sp,2]) = findindices(
                log(eval("hepsx"+runid,__main__.__dict__)[sp+1]),
                init_cfx[sp], final_cfx[sp], "-")
        (init_zy[sp,2], final_zy[sp,2]) = findindices(
                log(eval("hepsy"+runid,__main__.__dict__)[sp+1]),
                init_cfy[sp], final_cfy[sp], "-")

    # Avergaing Both Results:
        init_zx[sp,0]  = nint((init_zx[sp,1]+init_zx[sp,2])/2)
        init_zy[sp,0]  = nint((init_zy[sp,1]+init_zy[sp,2])/2)
        final_zx[sp,0] = nint((final_zx[sp,1]+final_zx[sp,2])/2)
        final_zy[sp,0] = nint((final_zy[sp,1]+final_zy[sp,2])/2)

    # --- Calculate Lyapunov Exponents from Least Squares Fit
        for i in range(0,3):
            (a, b, L) = lsqfit(eval("zscale"+runid,__main__.__dict__)[init_zx[sp,i]:final_zx[sp,i]],
                        log(eval("hepsx"+runid,__main__.__dict__)[sp+1,init_zx[sp,i]:final_zx[sp,i]]))
            lyexp_x[sp, i] = a
            lyexp_x[sp, i+6] = b
            lyexp_x[sp, i+3] = L

            (a, b, L) = lsqfit(eval("zscale"+runid,__main__.__dict__)[init_zy[sp,i]:final_zy[sp,i]],
                        log(eval("hepsy"+runid,__main__.__dict__)[sp+1,init_zy[sp,i]:final_zy[sp,i]]))
            lyexp_y[sp, i] = a
            lyexp_y[sp, i+6] = b
            lyexp_y[sp, i+3] = L

    print_lyap()


# ----------------------

def save_lyap(crun="0"):
    """ save_lyap(crun="0")
        Saves Lyapunov exponents into file """
    outfile = PW.PW("lyap."+runid+crun+".pdb")
    for vname in ("lyexp_x", "lyexp_y"):
        outfile.write( vname+runid, eval(vname, __main__.__dict__) )
    outfile.close()


def plot_lyap(jspec=None):
    "plot_lyap() retained as dummy for backward compatibility"
    pass


# =========================

def plot_diverge(lfull=0, jspec=None):
    """ plot_diverge(lfull=0, jspec=None)
        Plot difference b/w particle trajectories and centroids.
        If 'lfull' true, plots trajectories in each dimension separately;
        otherwise, only plots total separation.
        Plots over all species, except 0, by default, unless a range of species
        is provided via the parameter 'jspec'
    """
    if jspec is None: jspec = range(1, top.ns)
    absc = eval("zscale"+runid, __main__.__dict__)
    for sp in jspec:
        for part in range(1, num_part):
            eps_x[sp-1, :, part-1] = eval("x_"+runid, __main__.__dict__)[sp,:,part]-eval("x_"+runid, __main__.__dict__)[sp,:,0]
            eps_y[sp-1, :, part-1] = eval("y_"+runid, __main__.__dict__)[sp,:,part]-eval("y_"+runid, __main__.__dict__)[sp,:,0]
            eps_xp[sp-1, :, part-1] = eval("xp_"+runid, __main__.__dict__)[sp,:,part]-eval("xp_"+runid, __main__.__dict__)[sp,:,0]
            eps_yp[sp-1, :, part-1] = eval("yp_"+runid, __main__.__dict__)[sp,:,part]-eval("yp_"+runid, __main__.__dict__)[sp,:,0]

        if lfull:
            title = runid+": species "+`sp`+' '+colors[sp % ncolors]
            for plot in plots:
                for part in range(0, num_part-1):
                    plg(eval("eps_"+plot)[sp-1,:,part], absc)
                ptitles(title, "S (m)", plot+" separation"); fma()

    eps_t = sqrt(eps_x**2 + eps_y**2 + (eps_xp/k0x)**2 + (eps_yp/k0y)**2)

    for sp in jspec:
        title = runid+": species "+`sp`+' '+colors[sp % ncolors]
        for part in range(0, num_part-1):
            plg(eps_t[sp-1,:,part], absc)
        ptitles(title, "S (m)", "Total Separation"); fma()


def plot_traj(jspec=None, beg=0, end=n_steps):
    """  plot_traj(jspec=None, beg=0, end=n_steps)
    Poincare' plots of particle trajectories.
    Plots over all species, by default, unless a range of species
    is provided via the parameter 'jspec'
    """
    if jspec is None: jspec = range(0, top.ns)
    for sp in jspec:
        title = runid+": species "+`sp`+' '+colors[sp % ncolors]
        #
        for part in range(0, num_part):
            plg(eval("y_"+runid, __main__.__dict__)[sp,beg:end,part],
                eval("x_"+runid, __main__.__dict__)[sp,beg:end,part], type='dot')
        ptitles(title, "X", "Y"); fma()
        #
        for part in range(0, num_part):
            plg(eval("xp_"+runid, __main__.__dict__)[sp,beg:end,part],
                 eval("x_"+runid, __main__.__dict__)[sp,beg:end,part], type='dot')
        ptitles(title, "X", "X'"); fma()
        #
        for part in range(0, num_part):
            plg(eval("yp_"+runid, __main__.__dict__)[sp,beg:end,part],
                 eval("y_"+runid, __main__.__dict__)[sp,beg:end,part], type='dot')
        ptitles(title, "Y", "Y'"); fma()
        #
        for part in range(0, num_part):
            plg(eval("yp_"+runid, __main__.__dict__)[sp,beg:end,part],
                eval("xp_"+runid, __main__.__dict__)[sp,beg:end,part], type='dot')
        ptitles(title, "X'", "Y'"); fma()


def plot_spect(jspec=None, plist=[0], beg=0, end=n_steps, endp=None):
    """plot_spect(jspec=None, plist=[0], beg=0, end=n_steps, endp=end/2)
    Fourier Spectra of particle trajectories.
    Plots over all species, by default, unless a range of species
    is provided via the parameter 'jspec'
    'plist' is a list of particles to calculate spectrum over (Default is
    test particle #0).
    """
    import FFT
    if endp is None: endp = nint(end/2)
    if jspec is None: jspec = range(0, top.ns)
    for sp in jspec:
      for part in plist:
        title = runid+": species "+`sp`+' '+colors[sp % ncolors]+'; particle '+`part`
        #absc =
        for plot in plots:
            __main__.__dict__['spect_'+plot] = FFT.fft(
                eval(plot+'_'+runid, __main__.__dict__)[sp,beg:end,part])
            plg(abs(eval("spect_"+plot, __main__.__dict__))[:endp],  marker=plot)
            #plg(eval("spect_"+plot, __main__.__dict__), absc, marker=plot)
            ptitles(plot+' '+title, "frequency", "spectral power"); fma()


"""
Module orbitrack.py
by Rami A. Kishek
Created: May 29, 2001

Last Modified: Jan 30, 2002

Need to define a few variables before importing this function:
    'nrun' = # of steps to run
    'num_part' = number of particles to track
    top.runid

Contains functions for obtaining and saving particle trajectories and moments, and for
generating particle plots per species:

find_trajects() ... Extract trajectories of 'num_part' particles from each species
save_trajects() ... Save trajectories calculated using find_trajects to file
ppsp()          ... Plot set of particle plots with species colored differently
spmom()         ... Calculate Moments species-by-species
save_moms()     ... Save Moments calculated using spmom to file
"""

import __main__
from rami_scripts import *
import getzmom

orbitrack_version = "$Id: orbitrack.py,v 1.2 2002/08/14 21:07:53 ramiak Exp $"
def orbitrackdoc():
  import orbitrack
  print orbitrack.__doc__

# --- Define Constants

runid = arraytostr(top.runid)

tvars = ["x_", "y_", "xp_", "yp_"]
colors = ["fg", "red", "blue", "green", "magenta", "cyan"]

moms = def_vars1 + add_vars
nrun = __main__.__dict__["nrun"]
try: num_part = __main__.__dict__["num_part"]
except: num_part = 0

# --- Allocate arrays for moments and trajectories

for mom in moms:
    __main__.__dict__[mom+runid] = zeros((top.ns,(nrun/top.nhist)+1), 'd')
    __main__.__dict__["zscale"+runid] = zeros(((nrun/top.nhist)+1,), 'd')

for var in tvars:
    __main__.__dict__[var+runid] = zeros((top.ns,(nrun/top.nhist)+1,num_part), 'd')


#============================================================================

def find_trajects():
    """ Extract particle trajectories """
    for var in tvars:
        for spec in range(0,top.ns):
            __main__.__dict__[var+runid][spec][top.it/top.nhist] = eval(
                    "get"+var[:-1]+"(js="+`spec`+",iw=-1)[0:num_part]")


def save_trajects(crun="0"):
    """ save_trajects(crun="0")

        Function that dumps particle trajectories to file
        Trajectories calulated using find_trajects()
    """
    outfile = PW.PW("traj."+runid+crun+".pdb")
    for vname in tvars:
        outfile.write( vname+runid, eval(vname+runid, __main__.__dict__) )
    outfile.close()


# ----------------------

def ppsp(iwmb=-2, iwsp=-1, plots=('ppxy','ppxxp','ppyyp','ppxpyp'), scont=1, ncont=7):
    """ppsp(iwmb=-2, iwsp=-1, plots=('ppxy','ppxxp','ppyyp','ppxpyp'), scont=1, ncont=7)

    Plots a select set of particle plots given by the tuple "plots"
    with each species colored differently.  The particles in the main
    beam are weighted according to the window specified by "iwmb", swhile
    the minor species are weghted according to "iwsp".
        If "scont" then use shaded contour for main beam with 'ncont'
        contours, else use specimen.
    """
    for plot in plots:
        if scont:
            palette("gray.gp")
            apply(eval(plot), (), {'js': 0, 'iw': 0, 'color': 'fg',
                                   'contours':ncont, 'filled':1, 'particles':0})
        else:
            apply(eval(plot), (), {'js': 0, 'iw': iwmb, 'color': 'fg'})

        for spec in range(1,top.ns):
            apply(eval(plot), (), {'js': spec, 'iw': iwsp, 'color': colors[spec]})
        fma()


# ----------------------

def spmom(jslist=None):
    """ spmom(jslist=None)
    Calculate Moments species-by-species
    If jslist specified, combine moments over species in jslist
    """
    ifzmmnt = 1
    __main__.__dict__["zscale"+runid][top.it/top.nhist] = top.zbeam
    if jslist is None:
      for spec in range(top.ns-1, -1, -1):    # Reverse indicing so main beam ends up last calc
        getzmom.zmmnt(js=spec)
        for mom in moms:
            __main__.__dict__[mom+runid][spec, top.it/top.nhist] = eval("top."+mom[1:], __main__.__dict__)[0]
    else:
      getzmom.zmmnt(jslist=jslist)
      for mom in moms:
          __main__.__dict__[mom+runid][jslist[0], top.it/top.nhist] = eval("top."+mom[1:], __main__.__dict__)[0]
    ifzmmnt = 0


def save_moms(crun="0"):
    """ save_moms(crun="0")

        Mock save_data which saves the data file in the same format, but
        uses the moments calculated over a species, with the first index indicating
        the species rather than window number.
    """
    outfile = PW.PW("data."+runid+crun+".pdb")
    for mom in moms:
        outfile.write( mom+runid, eval(mom+runid, __main__.__dict__) )
    outfile.write( "zscale"+runid, eval("zscale"+runid, __main__.__dict__) )
    for var in def_vars2[1:]:
        __main__.__dict__[var+runid] = eval(seek_name(var))
        outfile.write( var+runid, eval(var+runid, __main__.__dict__) )
    outfile.close()


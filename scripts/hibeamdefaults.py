# --- Returns a value from the requested distribution function.
from Numeric import *
import RandomArray
hibeamdefaults_version = "$Id: hibeamdefaults.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"
def errordist(type):
  if type == 'GAUSSIAN':
    s = RandomArray.random(2)
    phi = 2.*pi*RandomArray.random(2)
    sq = sqrt(-2.*log(s))
    return tuple(sq*cos(phi))
    return
  elif type == 'ABSOLUTE':
    return (1.,1.)
  else:
    return tuple(2.*RandomArray.random(2)-1.)

#..HINITnamelistdefaults:
hinit={'current0':4.e-3,
       'current_fac':1.0,
       'emit0':0.,
       'emitn0':0.,
       'lambda':-1.,
       'emev0':0.160,
       'fchrome':1.0,
       'charge':1.0,
       'amu':133.0,
       'enorm0':2.0e-8,
       'npart':4096,
       'l_varydz':0,
       'dvzth':0.,
       'a':0.005,
       'b':0.005,
       'ap':0.,
       'bp':0.,
       'xoffset':0.,
       'yoffset':0.,
       'xpoffset':0.,
       'ypoffset':0.,
       'iseed_offset':0,
       'xdsize':0.05,
       'ydsize':0.05,
       'nx':256,
       'ny':256,
       'beam_name':['Beam#1','Beam#2','Beam#3','Beam#4'],
       'rquad':0.016,
       'rwall':0.020,
       'rapert':0.014,
       'iperiod':0,
       'irot':0,
       'idist':2,
       'icomb':0,
       'icapt':0,
       'ibcf':0,
       'l_use_capmatrix':1,
       'l_neut_charge':1,
       'l_env':0,
       'l_offset_one':0,
       'l_mad_lat_file':0,
       'single_cage':'NOTSET'}


#..ZTIMEnamelistdefaults:
ztime={'zmax':0.5,
       'zstart':0.0,
       'zstep0':0.005,
       'imerge':-1,
       'rzf_merge':2.0,
       'iswitch':-1,
       'irezone':-1,
       'rzf':1.0,
       'zparplot':-1.0,
       'zpardump':-1.0,
       'zphaseplot':-1.0,
       'zfldplot':-1.0,
       'dzfldplot':-1.,
       'zplot':-1.0,
       'dzphaseplot':-1.,
       'dzhist':-1.,
       'zhistnext':0.,
       'l_plot_cap_nodes':1,
       'dump_type':'BINARY',
       'l_print_hist':0}
       #'phasep_info(1:4)%zrange(1)':-1.,
       #'phasep_info(1:4)%zrange(2)':-1.,

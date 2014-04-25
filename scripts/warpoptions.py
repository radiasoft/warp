"""Command line options for Warp.
Users can add their own options by importing warpoptions (before importing
warp) add calling parser.add_argument. Exceptions generated by argument parser
in response to unknown command line arguments may be suppressed by setting
ignoreUnknownArgs (before importing warp)'
"""
import argparse

# --- if set this flag suppresses execetpions for
# --- unkown command line options found during command
# --- line parsing.
ignoreUnknownArgs = False

# --- if set the warp imports quietly
quietImport = False

# --- if set warp should skip argument parsing
lskipoptions = 0

# --- global parser object
parser = None
options = None
args = None

def warpoptionsdoc():
    import warpoptions
    print warpoptions.__doc__
    print warpoptions.warpoptionsstr()

def warpoptionsstr():
    return parser.format_help()

def init_parser():
    global parser
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS if quietImport else None,
        add_help=False if quietImport else True)

    parser.add_argument('-p', '--decomp', type=int, nargs=3, dest='decomp',
                      help='Set the 3-D decomposition, nxprocs nyprocs nzprocs')
    parser.add_argument('-l', '--localflags', dest='localflags',
                      help='localflags file to execfile')
    parser.add_argument('--pnumb', dest='pnumb', type=str, default=None,
                      help='Run number, used in the plot and other output file names.')

def parse_args():
    global options, args

    if lskipoptions:
        return

    if ignoreUnknownArgs:
        options,args = parser.parse_known_args()
    else:
        options = parser.parse_args()
        args = [] # --- currently no positional args are defined

# --- initialize default parser
init_parser()

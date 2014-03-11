"""Command line options for Warp.
Users can add their own options by importing warpoptions (before importing
warp) add calling parser.add_option.
"""
import optparse
warpoptions_version = "$Id: warpoptions.py,v 1.5 2010/11/15 23:01:06 dave Exp $"


def warpoptionsdoc():
    import warpoptions
    print warpoptions.__doc__
    print warpoptions.parser.parser.print_help()

parser = optparse.OptionParser()

parser.add_option('-p', '--decomp', type='int', nargs=3, dest='decomp',
                  help='Set the 3-D decomposition, nxprocs nyprocs nzprocs')
parser.add_option('-l', '--localflags', dest='localflags',
                  help='localflags file to execfile')
parser.add_option('--pnumb', dest='pnumb', type='string', default=None,
                  help='Run number, used in the plot and other output file names.')

lskipoptions = 0


def parse_args():
    global options, args
    if lskipoptions:
        return
    options, args = parser.parse_args()

"""Command line options for Warp.
Users can add their own options by importing warpoptions (before importing
warp) add calling parser.add_option.
"""
import optparse
warpoptions_version = "$Id: warpoptions.py,v 1.3 2009/06/12 19:06:21 dave Exp $"

def warpoptionsdoc():
    import warpoptions
    print warpoptions.__doc__
    print warpoptions.parser.parser.print_help()


parser = optparse.OptionParser()

parser.add_option('-p','--decomp',type='int',nargs=3,dest='decomp',
                  help='Set the 3-D decomposition, nxprocs nyprocs nzprocs')
parser.add_option('-l','--localflags',dest='localflags',
                  help='localflags file to execfile')

def parse_args():
    global options, args
    options, args = parser.parse_args()

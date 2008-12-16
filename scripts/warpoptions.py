"""Command line options for Warp
"""
import optparse
warpoptions_version = "$Id: warpoptions.py,v 1.2 2008/12/16 21:54:40 dave Exp $"

def warpoptionsdoc():
  import warpoptions
  print warpoptions.__doc__
  print warpoptions.parser.parser.print_help()


parser = optparse.OptionParser()

parser.add_option('-p','--decomp',type='int',nargs=3,dest='decomp',
                  help='Set the 3-D decomposition, nxprocs nyprocs nzprocs')
parser.add_option('-l','--localflags',dest='localflags',
                  help='localflags file to execfile')

options, args = parser.parse_args()

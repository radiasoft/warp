"""Command line options for Warp
"""
import optparse
warpoptions_version = "$Id: warpoptions.py,v 1.1 2008/11/19 22:22:45 dave Exp $"

def warpoptionsdoc():
  import warpoptions
  print warpoptions.__doc__
  print warpoptions.parser.parser.print_help()


parser = optparse.OptionParser()

parser.add_option('-p','--decomp',type='int',nargs=3,dest='decomp',
                  help='Set the 3-D decomposition, nxprocs nyprocs nzprocs')

options, args = parser.parse_args()

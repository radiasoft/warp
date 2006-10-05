#!/usr/bin/env python
import os
import sys
import time
t = time.localtime(time.time())
date = "%04d%02d%02d%02d%02d%02d"%tuple(t[:6])

if len(sys.argv) > 1:
  pubhomes = sys.argv[1:]
else:
  pubhomes = ['$HOME/pub']

def executecommand(s):
  print s
  e = os.system(s)

for pubhome in pubhomes:
  executecommand('cp warpC.so '+pubhome+'/warpC'+date+'.so')
  executecommand('cd '+pubhome+';chmod go+rx warpC'+date+'.so')
  executecommand('cd '+pubhome+';ln -sf warpC'+date+'.so warpC.so')
  executecommand('cd '+pubhome+'/source/pywarp90;cvs update -d')
  executecommand('cd '+pubhome+'/scripts;cvs update -d')
  executecommand('chmod -R go+rX '+pubhome)

#for pubhome in pubhomes:
#  executecommand('cp pywarp90 '+pubhome+'/pywarp90'+date)
#  executecommand('cp warpC.so '+pubhome+'/warpC'+date+'.so')
#  executecommand('cd '+pubhome+';chmod go+rx pywarp90'+date+' warpC'+date+'.so')
#  executecommand('cd '+pubhome+';ln -sf pywarp90'+date+' pywarp90')
#  executecommand('cd '+pubhome+';ln -sf warpC'+date+'.so warpC.so')
#  executecommand('cd '+pubhome+'/source/pywarp90;cvs update -d')
#  executecommand('cd '+pubhome+'/scripts;cvs update -d')


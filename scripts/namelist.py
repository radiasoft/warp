from Numeric import *
import re
namelist_version = "$Id: namelist.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

class NameList:

  def __init__(self,filename):
    self.filename = filename
    ff = open(filename,'r')
    self.file = ff.readlines()
    ff.close()
    #self.joincontinuedlines()
    self.getheader()
    self.getnamelists()

  def joincontinuedlines(self):
    n = len(self.file)
    i = 0
    while (i < n):
      if re.sub(r'&','',self.file[i]):
        self.file[i] = self.file[i] + self.file[i+1]
        del self.file[i+1]
        n = n - 1
      else:
        i = i + 1

  def getheader(self):
    self.header = ''
    while not re.search(r'\$',self.file[0]):
      self.header = self.header + self.file[0]
      del self.file[0]
    self.header = self.header + re.sub(r'\$','',self.file[0])

  def getnamelists(self):
    self.namelists = {}
    while 1:
      try:
        while not re.search(r'&',self.file[0]):
          del self.file[0]
      except IndexError:
        break
      m = re.search(r'&(\w*)\s*',self.file[0])
      name = m.group(1)
      del self.file[0]
      nl = {}
      while not re.search(r'/END',self.file[0]):
        m = re.search(r'\s*(\w*)\s*=\s*(.*)',self.file[0])
        nl[m.group(1)] = m.group(2)
        del self.file[0]
      del self.file[0]
      self.namelists[name] = nl



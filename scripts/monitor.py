"""
Sets up a remote run monitor. To use, first create the monitor in the
running job by inserting the following command somewhere in the input deck
before the long running commands:

createmonitor(passwd='mypassword')

Then, from a remotely running python, import this module and run the
following to create the connection

connect(machine="machinename",passwd='mypassword')

where "machinename" is the name of the machine where the original run
is. Note that the password must be the same as in the createmonitor
command. (The password provides some security since it is since using md5
encryption.  Though you should still probably not use a sensitive password!)
A "+++" prompt will be shown and commands typed will be sent to the remote job.

When the "+++" prompt is displayed, the remote job will stop execution
and wait for commands. To continue the job, type the control-d key.
"""
from warp import *
import gist
import socket
import time
import re
import md5
monitor_version = "$Id: monitor.py,v 1.5 2003/12/02 02:39:57 dave Exp $"

def socketsend(sock,s):
  """
Set up send command which first sends the number of bytes in the
message. The random number 39487 is subtracted from the number of bytes
sent as a form of minimal security, making it more difficult to hack
into the running job.
  """
  sock.send('%10d'%(len(s)-39487))
  sock.send(s)

def socketrecv(sock):
  """
Set up receive command which first sends the number of bytes in the
message. The random number 39487 is added from the number of bytes sent
as a form of minimal security, making it more difficult to hack into the
running job.
  """
  nbytes = sock.recv(10)
  nbytes = eval(nbytes) + 39487
  return sock.recv(nbytes)

# --- Create class for server side of run.
class Monitor:
  """
Creates a monitor for a running job.
 - port=50007: port number to use
 - passwd='0': password to provide some security
  """
  def __init__(self,port=50007,passwd='fj39jfgks'):
    global _ismonitored
    _ismonitored = true
    self.md5passwd = md5.new(passwd)
    self.hexdigestpasswd = self.md5passwd.hexdigest()
    self.port = port
    self.sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    self.sock.bind(("",port))
    self.sock.getsockname()
    self.sock.listen(1)
    self.sock.setblocking(0)
    self.client = []
    self.gettingcommands = 0
    installafterfs(self.checksocket)
    installafterstep(self.getcommands)
  def checksocket(self):
    """
Check if anyone has opened the port. If so, get the socket and set up
the code to run a special step so that all data is current when the
commands are received. Also, save maxcalls and ncall in case the
user does some additional step commands. If this is already getting
commands, then don't do anything.
    """
    if self.gettingcommands: return
    try:
      self.client = self.sock.accept()[0]
    except:
      if self.client: self.client.close()
      self.client = []
      return
    self.sock.setblocking(1)
    #self.sock.setblocking(0)
    # --- Not sure why this first send is needed, but it doesn't seem to 
    # --- work if it is not done, especially if the passwd is not ok.
    socketsend(self.client,'0')
    socketsend(self.client,self.hexdigestpasswd)
    passwdok = socketrecv(self.client)
    if passwdok == 'ok':
      self.lspecialsave = top.lspecial
      top.lspecial = 1
      self.maxcallssave = top.maxcalls
      self.ncallsave = top.ncall
      return
    else:
      print "Access not authorized"
     #socketsend(self.client,"Access not authorized")
      self.sock.setblocking(0)
      self.client.close()
      self.client = []
      return
  def getcommands(self):
    """
Get commands from a remote job. Checks if a client socket has been created,
and if so, get commands from it.
    """
    if self.gettingcommands: return
    import __main__
    if self.client:
      print "Ready for input."
      socketsend(self.client,"Ready for input.")
      self.gettingcommands = 1
      done = false
      while not done:
        try:
          comm = socketrecv(self.client)
        except:
          break
        if comm == "EOF":
          done = true
          print "Continuing"
          socketsend(self.client,"Continuing")
        elif comm == "exec":
          try:
            socketsend(self.client,"go")
            comm = socketrecv(self.client)
            print 'exec('+comm+')'
          except:
            break
          try:
            exec(comm,__main__.__dict__)
            socketsend(self.client,"OK")
          except:
            socketsend(self.client,"Error")
        else:
          try:
            print 'eval('+comm+')'
            rrr = eval(comm,__main__.__dict__)
            print rrr
            socketsend(self.client,repr(rrr))
          except:
            socketsend(self.client,"Error")
        # --- Handle any gist events
        try:
          v = gist.__version__
          pyg_pending()
          pyg_idler()
        except:
          ygdispatch()
      # --- When done, release the socket and return some variables
      # --- back to normal.
      self.client.close()
      self.sock.setblocking(0)
      self.client = []
      top.lspecial = self.lspecialsave
      top.maxcalls = self.maxcallssave
      top.ncall = self.ncallsave
      self.gettingcommands = 0

def createmonitor(port=50007,passwd='fj39jfgks'):
  """
Creates a monitor for a running job.
 - port=50007: port number to use
 - passwd='0': password to provide some security
  """
  global _monitor
  _monitor = Monitor(port=port,passwd=passwd)


###########################################################################
# --- Create functions for client side

# --- Send command to eval
def sendeval(s):
  global _sock
  socketsend(_sock,s)
  print socketrecv(_sock)

# --- Send command to exec
def sendexec(s):
  global _sock
  socketsend(_sock,'exec')
  r = socketrecv(_sock)
  socketsend(_sock,s)
  print socketrecv(_sock)

# --- Send continue command
def sendcont():
  sendeval('EOF')

# --- Reads command from stdin and sends them off. Assumes that anything
# --- containing an '=' should be exec'ed, otherwise eval'ed. That
# --- isn't always the case but it seems to work anyway.
def sendcommands():
  while 1:
    try:
      s = raw_input('+++ ')
    except EOFError:
      try:
        sendeval('EOF')
      except socket.error:
        pass
      break
    try:
      if not s:
        pass
      elif s == 'EOF':
        sendeval('EOF')
        break
      elif re.search('=',s):
        sendexec(s)
      elif s[:6] == 'print ':
        sendeval(s[6:])
      else:
        sendeval(s)
    except socket.error:
      break

# --- Create functions for client side
def connect(machine="localhost",port=50007,passwd='fj39jfgks',auto=1): 
  """
Make a connection to a running job with a monitor.
  - machine="localhost": machine on which job is running
  - port=50007: port the running job is monitoring
  - password='0': password the running job expects
  - auto=1: when 1, go directly into the "+++" prompt
  """
  global _sock
  try:
    if _ismonitored:
      print  "A monitored job cannot connect to another job"
      return "A monitored job cannot connect to another job"
  except:
    pass
  try:
    if _sock.fileno() != -1:
      print "Already connected!"
      return
  except:
    pass
  _sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
  _sock.connect((machine,port))
  md5passwd = md5.new(passwd)
  hexdigestpasswd = md5passwd.hexdigest()
  r = socketrecv(_sock)
  r = socketrecv(_sock)
  if r == hexdigestpasswd:
    socketsend(_sock,'ok')
  else:
    socketsend(_sock,'nope')
    print "Access not authorized"
    _sock.close()
    return
  print socketrecv(_sock)
  if auto:
    sendcommands()
    _sock.close()


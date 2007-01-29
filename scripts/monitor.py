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
import socket
import time
import re
import md5
monitor_version = "$Id: monitor.py,v 1.10 2007/01/29 19:05:26 dave Exp $"

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
    self.initializesocket()
    self.client = []
    self.gettingcommands = 0
    if self.sock is not None:
      installafterfs(self.checksocket)
      installafterstep(self.getcommands)
    else:
      raise "Error: quiting since socket could not be created"

  def initializesocket(self):
    try:
      self.sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
      self.sock.bind(("",self.port))
      self.sock.getsockname()
      self.sock.listen(1)
      self.sock.setblocking(0)
    except:
      self.sock = None
      print "Warning: socket could not be created"

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
    passwdok = 0
    try:
      # --- Catch all exceptions so running job isn't interrupted if a problem
      # --- occurs.
      self.sock.setblocking(1)
      #self.sock.setblocking(0)
      # --- Not sure why this first send is needed, but it doesn't seem to 
      # --- work if it is not done, especially if the passwd is not ok.
      socketsend(self.client,'0')
      socketsend(self.client,self.hexdigestpasswd)
      passwdok = socketrecv(self.client)
    except:
      pass
    if passwdok == 'ok':
      self.lspecialsave = top.lspecial
      top.lspecial = 1
      self.maxcallssave = top.maxcalls
      self.ncallsave = top.ncall
      return
    else:
      print "Error"
      try:
        socketsend(self.client,"Error")
      except:
        pass
      try:
        self.sock.setblocking(0)
        self.client.close()
      except:
        pass
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
      done = false
      try:
        socketsend(self.client,"Ready for input.")
        self.gettingcommands = 1
      except:
        print "Problem with socket - continuing job"
        done = true

      # --- Continue receiving commands until user finishes.
      while not done:

        # --- Get the command, check for errors in the socket.
        try:
          comm = socketrecv(self.client)
        except socket.error:
          print "Problem with socket - continuing job"
          break

        if comm == "EOF":
          # --- The user has finished.
          done = true
          result = "Continuing"
        else:
          try:
            # --- First, try to evaluate the command
            rrr = eval(comm,__main__.__dict__)
            print 'eval('+comm+')'
            result = repr(rrr)
          except SyntaxError:
            # --- If there was a syntax error, try to exec it.
            # --- The error might mean that this is an assignment.
            # --- Still, catch any errors in case there really are syntax
            # --- or other errors in the command.
            print 'exec('+comm+')'
            try:
              exec(comm,__main__.__dict__)
              result = "OK"
            except:
              result = "Error"
          except:
            result = "Error"

        # --- Print the result and send it back to the client.
        print result
        try:
          socketsend(self.client,result)
        except socket.error:
          print "Problem with socket - continuing job"
          break

        # --- Handle any gist events
        try:
          pyg_pending()
          pyg_idler()
        except KeyError:
          ygdispatch()

      # --- When done, release the socket and return some variables
      # --- back to normal.
      self.client.close()
      try:
        self.sock.setblocking(0)
      except:
        pass
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

# --- Send continue command
def sendcont():
  socketsend(_sock,'EOF')
  print socketrecv(_sock)

# --- Reads command from stdin and sends them off. Assumes that anything
# --- containing an '=' should be exec'ed, otherwise eval'ed. That
# --- isn't always the case but it seems to work anyway.
def sendcommands():
  while 1:
    try:
      s = raw_input('+++ ')
    except EOFError:
      try:
        socketsend(_sock,'EOF')
        print socketrecv(_sock)
      except socket.error:
        pass
      break
    try:
      if not s:
        pass
      elif s == 'EOF':
        socketsend(_sock,'EOF')
        print socketrecv(_sock)
        break
      elif s[:6] == 'print ':
        # --- Print is special since normally the quantity would only be
        # --- printed to the output of the remotely running job.
        socketsend(_sock,s[6:])
        print socketrecv(_sock)
      else:
        socketsend(_sock,s)
        print socketrecv(_sock)
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
  try:
    _sock.connect((machine,port))
  except socket.error:
    _sock.close()
    raise
  md5passwd = md5.new(passwd)
  hexdigestpasswd = md5passwd.hexdigest()
  r = socketrecv(_sock)
  r = socketrecv(_sock)
  if r == hexdigestpasswd:
    socketsend(_sock,'ok')
  else:
    socketsend(_sock,'nope')
    print "Error"
    _sock.close()
    return
  print socketrecv(_sock)
  if auto:
    sendcommands()
    _sock.close()


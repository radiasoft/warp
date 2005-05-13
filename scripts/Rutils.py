""" Module Rutils.py

by:      Rami A. Kishek
Created: April 4, 2003

Last Modified: 4/4/03

This module contains utilities for simplifying management of WARP simulations:

mailme ... sends e-mail message at that point of the run.
"""
Rutils_version = "$Id: Rutils.py,v 1.1 2005/05/13 06:09:00 ramiak Exp $"
import smtplib, os

def Rutilsdoc():
  import Rutils
  print Rutils.__doc__

def mailme(runid, user=None, addr=None, serv=None, host=None, text=None):
    """ mailme(runid, user=(Your username), addr='localhost', serv='localhost',
               host='localhost', text="Your run 'runid' on 'host' is over!")
        Sends e-mail message to user@addr from server 'serv' with 'text'
        containing 'runid' and 'host' running the simulation.
    """
    if user is None:    user = os.environ['USER']
    if serv is None:    serv = 'localhost'
    if addr is None:    addr = serv
    if host is None:    host = os.environ['HOST']
    if text is None:
        text = "Your run '%s' on '%s' is over!" % (runid, host)

    sender = '%s@%s' % (user, serv)
    recip  = '%s@%s' % (user, addr)
    server = smtplib.SMTP(serv)
    server.sendmail(sender, recip, text)
    server.quit()

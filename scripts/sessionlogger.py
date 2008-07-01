"""SessionLogger writes a log of the input and output of an interactice session

Notes:
  When writing out the prompt, it would be nice if there was a way of knowing
  whether to use the ps1 or ps2 prompt. As far as I can tell, that info is not
  available.

"""
__all__ = ['SessionLogger']

import os
import sys
import time
import readline

class ForkedFile(object):
    """
This class forks the output to both an internal buffer and to the write method
of the original file.
    """
    def __init__(self,file):
        self.file = file
        self.resetbuffer()
    def write(self,string):
        self.buffer += string
        self.file.write(string)
    def resetbuffer(self):
        self.buffer = ''
    def flush(self):
        self.file.flush()

class SessionLogger(object):
    def __init__(self,logfilename,comments=None):
        # --- Open the file
        self.logfilename = logfilename
        self.logfile = open(logfilename,'w')

        # --- Creates forks of both stdout and stderr and replace the original
        # --- ones with the forked ones.
        self.bufferedstdout = ForkedFile(sys.stdout)
        self.bufferedstderr = ForkedFile(sys.stderr)
        sys.stdout = self.bufferedstdout
        sys.stderr = self.bufferedstderr

        # --- Write any initial stuff to the log, including a time stamp
        self.logfile.write('#### '+time.asctime(time.localtime())+' ####\n')
        self.logfile.write('Python '+sys.version+'\n')
        if comments is not None: self.logfile.write(comments+'\n')

        # --- The session log writer will be called just before readline
        # --- writes the prompt.
        readline.set_startup_hook(self._writesessionlog)

    def _writeinput(self):
        n = readline.get_current_history_length()
        input = readline.get_history_item(n)
        self.logfile.write(sys.ps1+input+'\n')

    def _writeoutput(self):
        self.logfile.write(self.bufferedstdout.buffer)
        self.logfile.write(self.bufferedstderr.buffer)
        self.bufferedstdout.resetbuffer()
        self.bufferedstderr.resetbuffer()

    def _writesessionlog(self):
        self._writeinput()
        self._writeoutput()

    def close(self):
        # --- Finish writing the buffers
        self._writesessionlog()

        # --- Close the file
        self.logfile.close()

        # --- Remove readline startup hook and reset stdout and err
        readline.set_startup_hook()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__


from warp import *

class newstdout_accumulate:
    def __init__(self):
        self.stdout = []
    def write(self,s):
        self.stdout.append(s)
    def clear(self):
        self.stdout = []
    def flush(self):
        pass         

class newstdout_withCR:
    def __init__(self,winout):
        self.winout = winout
        self.lastCR = false 
        import curses.ascii
        self.CR = curses.ascii.ctrl('m')
    def write(self,s):
        if(self.lastCR):
            l = self.winout.GetLineLength(self.winout.GetNumberOfLines()-1)
            self.winout.Remove(self.winout.GetLastPosition()-l,self.winout.GetLastPosition())
            self.lastCR = false
        if(s==self.CR):
            self.lastCR = true
        else:
            self.winout.WriteText(s)
    def flush(self):
        pass         

class newstdout:
    def __init__(self,winout):
        self.winout = winout
    def write(self,s):
        self.winout.SetInsertionPointEnd()
        self.winout.write(s)
    def flush(self):
        pass         

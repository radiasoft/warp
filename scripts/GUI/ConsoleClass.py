#Boa:FramePanel:ConsoleClass

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from warp import *

[wxID_CONSOLECLASS, wxID_CONSOLECLASSCONSOLE, 
] = map(lambda _init_ctrls: wxNewId(), range(2))

class ConsoleClass(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_CONSOLECLASS, name='', parent=prnt,
              pos=wxPoint(355, 226), size=wxSize(604, 339),
              style=wxMAXIMIZE_BOX | wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(596, 315))
        self.SetAutoLayout(True)

        self.Console = wxTextCtrl(id=wxID_CONSOLECLASSCONSOLE, name='Console',
              parent=self, pos=wxPoint(0, 0), size=wxSize(596, 315),
              style=wxTE_MULTILINE, value='')
        self.Console.SetFont(wxFont(12, wxMODERN, wxNORMAL, wxNORMAL, false,''))
        self.Console.SetConstraints(LayoutAnchors(self.Console, True, True,
              True, True))
        EVT_CHAR(self.Console, self.OnconsoleChar)
        EVT_RIGHT_DOWN(self.Console, self.OnConsoleRightDown)
        EVT_LEFT_UP(self.Console, self.OnConsoleLeftUp)

    def __init__(self, parent, inter, in_notebook=0):
        self._init_ctrls(parent)
        if in_notebook:parent.InsertPage(index=0, imageId=-1, page=self, select=False, text='Console')
        self.listcommands=[]
        self.numcommand=-1
        self.CursorMin = 0
        self.inter = inter
        self.NoEntry = 1
        self.inter.push('import sys')
        self.inter.push("sys.path.append('')")
        self.inter.push('from warp import *')
        self.Console.WriteText('>>> ')
        self.ins_point = self.Console.GetInsertionPoint()

    def send_command_line(self):
        self.NoEntry = 0
        nline = self.Console.GetNumberOfLines()       
        line = self.Console.GetLineText(nline-1)
        self.Console.WriteText('\n')
        nl = len(line)
        r=self.sendcommand(line[4:nl])
        if(r==0):
          self.Console.WriteText('>>> ')
        else:
          self.Console.WriteText('... ')
        self.CursorMin = self.Console.GetLastPosition()
        self.ins_point = self.CursorMin
       
    def sendcommand(self,c,addlist=1):
        r = self.inter.push(c)
        if(addlist and c<>""):
          self.listcommands = self.listcommands + [c]
          self.numcommand = len(self.listcommands)
        return r

    def OnconsoleChar(self, event):
        ascii_code = event.GetKeyCode()
#        print ascii_code
        if(ascii_code == 317): # arrow up
            if(self.numcommand==len(self.listcommands)):
                nline = self.Console.GetNumberOfLines()       
                line = self.Console.GetLineText(nline-1)
                nl = len(line)
                self.LineBuffer = line[4:nl]
            self.numcommand = max(self.numcommand-1,0)
            self.Console.Remove(self.CursorMin,self.Console.GetLastPosition())
            self.Console.WriteText(self.listcommands[self.numcommand])
            self.Reset_ins_point()
        elif(ascii_code == 319): # arrow down
            self.numcommand = min(self.numcommand+1,len(self.listcommands))
            self.Console.Remove(self.CursorMin,self.Console.GetLastPosition())
            if(self.numcommand<len(self.listcommands)):
              self.Console.WriteText(self.listcommands[self.numcommand])
            else:
              self.Console.WriteText(self.LineBuffer)
            self.Reset_ins_point()
        elif(ascii_code == 316): # arrow left
            self.Set_ins_point()
            self.Console.SetInsertionPoint(max(self.ins_point-1,self.CursorMin))
            self.ins_point -= 1
        elif(ascii_code == 318): # arrow right
            self.Set_ins_point()
            self.Console.SetInsertionPoint(min(self.ins_point+1,self.Console.GetLastPosition()))
            self.ins_point += 1
        elif(ascii_code == 8): # backspace
            self.Set_ins_point()
            rmv_point = max(self.ins_point-1,self.CursorMin)
            self.Console.Remove(rmv_point,rmv_point+1)
            self.ins_point -= 1
        elif(ascii_code == 127): # delete
            self.Set_ins_point()
            rmv_point = min(self.ins_point,self.Console.GetLastPosition())
            self.Console.Remove(rmv_point,rmv_point+1)
        elif(ascii_code == 13): # return
            self.Console.SetInsertionPoint(self.Console.GetLastPosition())
            self.send_command_line()
        elif(ascii_code>=32 and ascii_code<255): # character       
            self.Set_ins_point()
            self.Console.SetInsertionPoint(self.ins_point)
            self.Console.WriteText(chr(ascii_code))
            self.ins_point+=1
        else:
            pass

    def Set_ins_point(self):
        newins_point = self.Console.GetInsertionPoint()
        if(self.ins_point<>newins_point and newins_point>=self.CursorMin):
           self.ins_point = newins_point
        
    def Reset_ins_point(self):
        self.ins_point = self.Console.GetInsertionPoint()
        
    def OnConsoleLeftUp(self, event):
        newins_point = self.Console.GetInsertionPoint()
        if(self.ins_point<>newins_point and newins_point>=self.CursorMin):
           self.ins_point = newins_point
        event.Skip()

    def OnConsoleRightDown(self, event):
        self.GetParent().GetParent().GetParent().OutToMessageWindow()
        print doc(self.Console.GetStringSelection(),printit=0)
        self.GetParent().GetParent().GetParent().OutToConsole()
        event.Skip()



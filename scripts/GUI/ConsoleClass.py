#Boa:FramePanel:ConsoleClass

from wxPython.wx import *

[wxID_CONSOLECLASS, wxID_CONSOLECLASSCONSOLE, 
] = map(lambda _init_ctrls: wxNewId(), range(2))

class ConsoleClass(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_CONSOLECLASS, name='', parent=prnt,
              pos=wxPoint(0, 0), size=wxSize(596, 315),
              style=wxMAXIMIZE_BOX | wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(596, 315))

        self.Console = wxTextCtrl(id=wxID_CONSOLECLASSCONSOLE, name='Console',
              parent=self, pos=wxPoint(0, 0), size=wxSize(596, 315),
              style=wxTE_MULTILINE, value='')
        self.Console.SetFont(wxFont(12, wxMODERN, wxNORMAL, wxNORMAL, false,''))
        EVT_CHAR(self.Console, self.OnconsoleChar)

    def __init__(self, parent, inter):
        self._init_ctrls(parent)
        parent.AddPage(bSelect=true, imageId=-1, pPage=self, strText='Console')
        self.listcommands=[]
        self.numcommand=-1
        self.CursorMin = 0
        self.inter = inter
        self.NoEntry = 1
        self.inter.push('import sys')
        self.inter.push("sys.path.append('')")
        self.inter.push('from warp import *')
        self.Console.WriteText('>>> ')

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
        elif(ascii_code == 319): # arrow down
            self.numcommand = min(self.numcommand+1,len(self.listcommands))
            self.Console.Remove(self.CursorMin,self.Console.GetLastPosition())
            if(self.numcommand<len(self.listcommands)):
              self.Console.WriteText(self.listcommands[self.numcommand])
            else:
              self.Console.WriteText(self.LineBuffer)
        elif(ascii_code == 316): # arrow left
            ins_point = self.Console.GetInsertionPoint()
            self.Console.SetInsertionPoint(max(ins_point-1,self.CursorMin))
        elif(ascii_code == 318): # arrow right
            ins_point = self.Console.GetInsertionPoint()
            self.Console.SetInsertionPoint(min(ins_point+1,self.Console.GetLastPosition()))
        elif(ascii_code == 8): # backspace
            ins_point = self.Console.GetInsertionPoint()
            rmv_point = max(ins_point-1,self.CursorMin)
            self.Console.Remove(rmv_point,rmv_point+1)
        elif(ascii_code == 127): # delete
            ins_point = self.Console.GetInsertionPoint()
            rmv_point = min(ins_point,self.Console.GetLastPosition())
            self.Console.Remove(rmv_point,rmv_point+1)
        elif(ascii_code == 13): # return
            self.Console.SetInsertionPoint(self.Console.GetLastPosition())
            self.send_command_line()
        elif(ascii_code>=32 and ascii_code<255): # character       
            self.Console.WriteText(chr(ascii_code))
        else:
            pass

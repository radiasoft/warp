#Boa:Frame:WarpRun

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
import WarpGUIInfo
import ParticlePlotsGUI
import EnvelopeGUI
import LatticeGUI
import DocGUI
import MatchingGUI
import ConsoleClass
import PzplotsGUI
import txtEditorDialog
import newstdout
import gist
import sys
import code
import __main__

from warp import *
from errorcheck import *

def create(parent):
    return WarpRun(parent)

[wxID_WARPRUN, wxID_WARPRUNCONT, wxID_WARPRUNDOC, wxID_WARPRUNENV, 
 wxID_WARPRUNFMA, wxID_WARPRUNHCP, wxID_WARPRUNLAT, wxID_WARPRUNMESSAGEWINDOW, 
 wxID_WARPRUNNEXT, wxID_WARPRUNNOTEBOOK1, wxID_WARPRUNPANEL1, 
 wxID_WARPRUNSTART, wxID_WARPRUNSTATUSBAR1, wxID_WARPRUNSTEP, 
 wxID_WARPRUNTXTEDITOR, wxID_WARPRUNWINON, 
] = map(lambda _init_ctrls: wxNewId(), range(16))

[wxID_WARPRUNTOOLBAR2TOOLS0, wxID_WARPRUNTOOLBAR2TOOLS1, wxID_WARPRUNTOOLBAR2TOOLS2,
 wxID_WARPRUNTOOLBAR2TOOLS3] = map(lambda _init_coll_toolBar2_Tools: wxNewId(), range(4))

[wxID_WARPRUNTOOLBAR1TOOLS0, wxID_WARPRUNTOOLBAR1TOOLS1, 
 wxID_WARPRUNTOOLBAR1TOOLS2, wxID_WARPRUNTOOLBAR1TOOLS3, 
 wxID_WARPRUNTOOLBAR1TOOLS4, wxID_WARPRUNTOOLBAR1TOOLS5, 
 wxID_WARPRUNTOOLBAR1TOOLS6, 
] = map(lambda _init_coll_toolBar1_Tools: wxNewId(), range(7))

[wxID_WARPRUNMNUERRORCHECKCHECKALL, wxID_WARPRUNMNUERRORCHECKENVELOPE, 
 wxID_WARPRUNMNUERRORCHECKIBPUSH, wxID_WARPRUNMNUERRORCHECKPARTICLELOAD, 
 wxID_WARPRUNMNUERRORCHECKSYMMETRY, 
] = map(lambda _init_coll_mnuErrorCheck_Items: wxNewId(), range(5))

[wxID_WARPRUNMNUPACKAGE3D, wxID_WARPRUNMNUPACKAGEENV, 
 wxID_WARPRUNMNUPACKAGEXY, 
] = map(lambda _init_coll_mnuPackage_Items: wxNewId(), range(3))

[wxID_WARPRUNMNUFILEEXEC, wxID_WARPRUNMNUFILEEXIT, wxID_WARPRUNMNUFILEOPEN, 
 wxID_WARPRUNMNUFILEOPENEXEC, wxID_WARPRUNMNUFILESAVE, 
 wxID_WARPRUNMNUFILESAVEAS, 
] = map(lambda _init_coll_mnuFile_Items: wxNewId(), range(6))

[wxID_WARPRUNMNUDUMPDUMP, wxID_WARPRUNMNUDUMPDUMPAS, wxID_WARPRUNMNUDUMPRESTORE, 
 wxID_WARPRUNMNUDUMPRESTART, 
] = map(lambda _init_coll_mnuFile_Items: wxNewId(), range(4))

[wxID_WARPRUNMNUDUMPDUMP, wxID_WARPRUNMNUDUMPDUMPAS, 
 wxID_WARPRUNMNUDUMPRESTART, wxID_WARPRUNMNUDUMPRESTORE, 
] = map(lambda _init_coll_mnuDump_Items: wxNewId(), range(4))

[wxID_WARPRUNMNUHELPABOUT] = map(lambda _init_coll_mnuHelp_Items: wxNewId(), range(1))

class WarpRun(wxFrame):
    def _init_coll_mnuHelp_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='Display info', id=wxID_WARPRUNMNUHELPABOUT,
              item='About', kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPABOUT, self.OnMnuhelpAboutMenu)

    def _init_coll_mnuFile_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='', id=wxID_WARPRUNMNUFILEOPEN, item='Open',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='Opens and Executes file',
              id=wxID_WARPRUNMNUFILEOPENEXEC, item='Open/Execfile',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUFILESAVE, item='Save',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUFILESAVEAS,
              item='Save As', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUFILEEXEC,
              item='ExecFile', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUFILEEXIT, item='Exit',
              kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUFILEOPEN, self.OnMnuOpenMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILESAVE, self.OnMnufileSaveMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILESAVEAS, self.OnMnufileSaveAsMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEEXIT, self.OnMnufileExitMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEEXEC, self.OnMnufileexecfileMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEOPENEXEC, self.OnMnufileOpenExecMenu)

    def _init_coll_menuBar1_Menus(self, parent):
        # generated method, don't edit

        parent.Append(menu=self.mnuFile, title='File')
        parent.Append(menu=self.mnuDump, title='Dump')
        parent.Append(menu=self.mnuHelp, title='Help')
        parent.Append(menu=self.mnuErrorCheck, title='ErrorCheck')
        parent.Append(menu=self.mnuPackage, title='Package')
        parent.Append(menu=self.mnuPalette, title='Palette')

    def _init_coll_mnuPackage_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='Select 3-D code', id=wxID_WARPRUNMNUPACKAGE3D,
              item='3-D', kind=wxITEM_CHECK)
        parent.Append(helpString='Select slice code',
              id=wxID_WARPRUNMNUPACKAGEXY, item='X-Y', kind=wxITEM_CHECK)
        parent.Append(helpString='Select envelope code',
              id=wxID_WARPRUNMNUPACKAGEENV, item='Envelope', kind=wxITEM_CHECK)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGE3D, self.OnMnupackage3dMenu)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGEXY, self.OnMnupackageXYMenu)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGEENV, self.OnMnupackageEnvMenu)

    def _init_coll_mnuDump_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='', id=wxID_WARPRUNMNUDUMPRESTORE,
              item='Restore', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUDUMPRESTART,
              item='Restart', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUDUMPDUMP, item='Dump',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUDUMPDUMPAS,
              item='Dump As', kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUDUMPRESTORE, self.OnMnudumpRestore)
        EVT_MENU(self, wxID_WARPRUNMNUDUMPRESTART, self.OnMnudumpRestart)
        EVT_MENU(self, wxID_WARPRUNMNUDUMPDUMP, self.OnMnudumpDump)
        EVT_MENU(self, wxID_WARPRUNMNUDUMPDUMPAS, self.OnMnudumpDumpAs)

    def _init_coll_mnuErrorCheck_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='', id=wxID_WARPRUNMNUERRORCHECKSYMMETRY,
              item='Symmetry', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUERRORCHECKPARTICLELOAD,
              item='ParticleLoad', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUERRORCHECKENVELOPE,
              item='Envelope', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUERRORCHECKIBPUSH,
              item='Ibpush', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WARPRUNMNUERRORCHECKCHECKALL,
              item='CheckAll', kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUERRORCHECKSYMMETRY,
              self.OnMnuerrorchecksymmetryMenu)
        EVT_MENU(self, wxID_WARPRUNMNUERRORCHECKPARTICLELOAD,
              self.OnMnuerrorcheckparticleloadMenu)
        EVT_MENU(self, wxID_WARPRUNMNUERRORCHECKENVELOPE,
              self.OnMnuerrorcheckenvelopeMenu)
        EVT_MENU(self, wxID_WARPRUNMNUERRORCHECKIBPUSH,
              self.OnMnuerrorcheckibpushMenu)
        EVT_MENU(self, wxID_WARPRUNMNUERRORCHECKCHECKALL,
              self.OnMnuerrorcheckallMenu)

    def _init_coll_notebook1_Pages(self, parent):
        # generated method, don't edit

        parent.AddPage(imageId=-1, page=self.txtEditor, select=True,
              text='Editor')

    def _init_coll_statusBar1_Fields(self, parent):
        # generated method, don't edit
        parent.SetFieldsCount(1)

        parent.SetStatusText(i=0, text='Status')

        parent.SetStatusWidths([-1])

    def _init_utils(self):
        # generated method, don't edit
        self.mnuFile = wxMenu(title='File')
        self._init_coll_mnuFile_Items(self.mnuFile)

        self.mnuDump = wxMenu(title='Dump')
        self._init_coll_mnuDump_Items(self.mnuDump)

        self.mnuHelp = wxMenu(title='Help')
        self._init_coll_mnuHelp_Items(self.mnuHelp)

        self.menuBar1 = wxMenuBar()

        self.mnuErrorCheck = wxMenu(title='ErrorCheck')
        self._init_coll_mnuErrorCheck_Items(self.mnuErrorCheck)

        self.mnuPackage = wxMenu(title='Package')
        self._init_coll_mnuPackage_Items(self.mnuPackage)

        self.mnuPalette = wxMenu(title='Palette')

        self._init_coll_menuBar1_Menus(self.menuBar1)

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxFrame.__init__(self, id=wxID_WARPRUN, name='WarpRun', parent=prnt,
              pos=wxPoint(520, 164), size=wxSize(600, 604),
              style=wxDEFAULT_FRAME_STYLE, title='WARP')
        self._init_utils()
        self.SetClientSize(wxSize(600, 580))
        self.SetMenuBar(self.menuBar1)
        EVT_IDLE(self, self.OnWxframe1Idle)

        self.statusBar1 = wxStatusBar(id=wxID_WARPRUNSTATUSBAR1,
              name='statusBar1', parent=self, style=0)
        self.statusBar1.SetSize(wxSize(550, 25))
        self.statusBar1.SetPosition(wxPoint(0, 596))
        self._init_coll_statusBar1_Fields(self.statusBar1)
        self.SetStatusBar(self.statusBar1)

        self.panel1 = wxPanel(id=wxID_WARPRUNPANEL1, name='panel1', parent=self,
              pos=wxPoint(0, 0), size=wxSize(600, 24), style=wxTAB_TRAVERSAL)

        self.winon = wxButton(id=wxID_WARPRUNWINON, label='win', name='winon',
              parent=self.panel1, pos=wxPoint(0, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.winon, wxID_WARPRUNWINON, self.OnWinonButton)

        self.fma = wxButton(id=wxID_WARPRUNFMA, label='fma', name='fma',
              parent=self.panel1, pos=wxPoint(40, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.fma, wxID_WARPRUNFMA, self.OnFmaButton)

        self.hcp = wxButton(id=wxID_WARPRUNHCP, label='hcp', name='hcp',
              parent=self.panel1, pos=wxPoint(80, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.hcp, wxID_WARPRUNHCP, self.OnHcpButton)

        self.env = wxButton(id=wxID_WARPRUNENV, label='env', name='env',
              parent=self.panel1, pos=wxPoint(120, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.env, wxID_WARPRUNENV, self.OnEnvButton)

        self.lat = wxButton(id=wxID_WARPRUNLAT, label='lat', name='lat',
              parent=self.panel1, pos=wxPoint(160, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.lat, wxID_WARPRUNLAT, self.OnLatButton)

        self.doc = wxButton(id=wxID_WARPRUNDOC, label='doc', name='doc',
              parent=self.panel1, pos=wxPoint(200, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.doc, wxID_WARPRUNDOC, self.OnDocButton)

        self.notebook1 = wxNotebook(id=wxID_WARPRUNNOTEBOOK1, name='notebook1',
              parent=self, pos=wxPoint(0, 24), size=wxSize(600, 350), style=0)
        EVT_NOTEBOOK_PAGE_CHANGED(self.notebook1, wxID_WARPRUNNOTEBOOK1,
              self.OnNotebook1NotebookPageChanged)

        self.txtEditor = wxTextCtrl(id=wxID_WARPRUNTXTEDITOR, name='txtEditor',
              parent=self.notebook1, pos=wxPoint(0, 0), size=wxSize(596, 316),
              style=wxTE_MULTILINE, value='')
        self.txtEditor.SetToolTipString('Text Editor')

        self.MessageWindow = wxTextCtrl(id=wxID_WARPRUNMESSAGEWINDOW,
              name='MessageWindow', parent=self, pos=wxPoint(2, 376),
              size=wxSize(596, 208), style=wxTE_READONLY | wxTE_MULTILINE,
              value='')
        self.MessageWindow.SetFont(wxFont(12, wxMODERN, wxNORMAL, wxNORMAL,
              false, ''))
        self.MessageWindow.SetBackgroundColour(wxColour(192, 192, 192))

        self.Step = wxButton(id=wxID_WARPRUNSTEP, label='Step', name='Step',
              parent=self.panel1, pos=wxPoint(280, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.Step, wxID_WARPRUNSTEP, self.OnStepButton)

        self.Next = wxButton(id=wxID_WARPRUNNEXT, label='Next', name='Next',
              parent=self.panel1, pos=wxPoint(320, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.Next, wxID_WARPRUNNEXT, self.OnNextButton)

        self.Start = wxButton(id=wxID_WARPRUNSTART, label='Start', name='Start',
              parent=self.panel1, pos=wxPoint(240, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.Start, wxID_WARPRUNSTART, self.OnStartButton)

        self.Cont = wxButton(id=wxID_WARPRUNCONT, label='Cont', name='Cont',
              parent=self.panel1, pos=wxPoint(360, 0), size=wxSize(40, 22),
              style=0)
        EVT_BUTTON(self.Cont, wxID_WARPRUNCONT, self.OnContButton)

        self._init_coll_notebook1_Pages(self.notebook1)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.FileName = None
        #EVT_UPDATE_UI(,self.mnuPackageUpdate)
        self.mnuPackageUpdate()
        self.isgistwindowon = 0
        self.oldline = '#'
        self.linenum = 0
        self.EdPos = 0
        self.startrun = 1
        self.inter = code.InteractiveConsole(__main__.__dict__)
        self.ConsolePanel = ConsoleClass.ConsoleClass(parent=self.notebook1, inter=self.inter)
        self.Console = self.ConsolePanel.Console
        self.prefix = ''
        self.PplotsPanel = ParticlePlotsGUI.ParticlePlotsGUI(self.notebook1)
        self.PzplotsGUI = PzplotsGUI.Pzplots(self.notebook1)
        self.MatchingGUI = MatchingGUI.MatchingGUI(self.notebook1)
        self.notebook1.SetSelection(0) # open notebook on Editor
        self.FileExecDialog = txtEditorDialog.txtEditorDialog(self)      
        self.FileExec = self.FileExecDialog.txtEditor  
#        self.FileExec.Show(1)
        Palettes = ["earth","rainbow","gray","yarg","heat","ncar","cool","rainbowaf","stern","christmas"]
        for i in range(0,len(Palettes)):
            self.AddPalette(Palettes[i])
   
    def HandleGistEvents(self):
      try:
        v = gist.__version__
        pyg_pending()
        pyg_idler()
      except:
        ygdispatch()

    def AddPalette(self,name):
        exec("[wxID_WARPRUNMNUPALLETTE"+name+"] = map(lambda _init_coll_mnuPalette_Items: wxNewId(), range(1))")
        exec("self.mnuPalette.Append(helpString='', id=wxID_WARPRUNMNUPALLETTE"+name+", item='"+name+"',kind=wxITEM_NORMAL)")
        exec("def OnMnuPalette"+name+"(event):palette('"+name+".gp')")
        exec("EVT_MENU(self, wxID_WARPRUNMNUPALLETTE"+name+", OnMnuPalette"+name+")")

    def OnMnuhelpAboutMenu(self, event):
        dlg = WarpGUIInfo.WarpGUIInfo(self)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def OnMnufileOpenExecMenu(self, event):
        self.OnMnuOpenMenu(event)
        self.OnMnufileexecfileMenu(event)

    def OnMnuOpenMenu(self, event):
        dlg = wxFileDialog(self, "Choose a file", ".", "", 
              "PYTHON files (*.py)|*.py|ALL files (*.*)|*.*", wxOPEN)
        try:
            if dlg.ShowModal() == wxID_OK:
                filename = dlg.GetPath()
                self.txtEditor.LoadFile(filename)
                self.FileName = filename
                self.FileExec.SetValue(self.txtEditor.GetValue())
                self.SetInputsLimits()
                try:
                    self.ImportGui(filename)
                except:
                    pass
        finally:
            dlg.Destroy()

    def SetInputsLimits(self):
        self.StartInputs = string.find(self.FileExec.GetValue(),'#Inputs')+len('#Inputs\n')
        self.EndInputs   = string.find(self.FileExec.GetValue(),'#EndInputs')
        if(self.StartInputs>0 and self.EndInputs>0):
            self.Inputs = self.FileExec.GetValue()[self.StartInputs:self.EndInputs]
        
    def ChangeInputs(self,newinputs):
        self.FileExec.Remove(self.StartInputs,self.EndInputs)
        self.FileExec.SetInsertionPoint(self.StartInputs)
        self.FileExec.WriteText(newinputs)
        self.FileExec.SetInsertionPoint(0)
        self.inputs = newinputs
        self.SetInputsLimits()

    def OnMnufileSaveMenu(self, event):
        if self.FileName is None:
            return OnMnufileSaveAsMenu(event)
        else:
            self.txtEditor.SaveFile(self.FileName)
            self.FileExec.SetValue(self.txtEditor.GetValue())
            self.SetInputsLimits()
 
    def OnMnufileSaveAsMenu(self, event):
        dlg = wxFileDialog(self, "Choose a file", ".", "", 
              "PYTHON files (*.py)|*.py|ALL files (*.*)|*.*", wxSAVE|wxOVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wxID_OK:
                filename = dlg.GetPath()
                self.txtEditor.SaveFile(filename)
                self.FileName = filename
                self.FileExec.SetValue(self.txtEditor.GetValue())
                self.SetInputsLimits()
        finally:
            dlg.Destroy()

    def OnMnufileExitMenu(self, event):
        self.Close()

    def OnMnufileexecfileMenu1(self, event):
        if self.FileName is None:
            OnMnufileSaveAsMenu(event)
        self.statusBar1.SetStatusText(i=0,text="Executing file %s"%self.FileName)
        execfile(self.FileName)
        self.statusBar1.SetStatusText(i=0,text="Finished executing file %s"%self.FileName)

    def OnMnufileexecfileMenu(self, event):
        self.notebook1.SetSelection(1) # open notebook on Editor
        sys.stdout = newstdout.newstdout(self.Console)
        sys.stderr = newstdout.newstdout(self.Console)
        if self.FileName is None:
            return
        self.statusBar1.SetStatusText(i=0,text="Executing file %s"%self.FileName)
        if(self.startrun and self.ConsolePanel.NoEntry):
            self.Console.Clear()
            startrun = 0
        self.OnContButton(event)

    def OnStartButton(self, event):
        self.FileExec.SetValue(self.txtEditor.GetValue())
        self.SetInputsLimits()
        self.OnMnufileexecfileMenu(event)

    def OnContButton(self, event):
        self.notebook1.SetSelection(1) # open notebook on Editor
        sys.stdout = newstdout.newstdout(self.Console)
        sys.stderr = newstdout.newstdout(self.Console)
        self.SetStatusText('Running')
        dorun = true
        self.notebook1.SetSelection(1) # open notebook on Console
        while(dorun and self.linenum<=self.FileExec.GetNumberOfLines()):
            dorun = self.AnalyzeNextLine(action='next')
            self.prefix='>>> '
        self.ReturnToPrompt(self.line)
        self.prefix=''

    def OnNextButton(self, event):
        if self.linenum<=self.FileExec.GetNumberOfLines() and self.FileExec.GetNumberOfLines()>0:
            self.SetStatusText('Running')
            self.notebook1.SetSelection(1) # open notebook on Console
            dorun = self.AnalyzeNextLine(action='next')
            self.ReturnToPrompt(self.line)

    def OnStepButton(self, event):
#        self.Console.SetInsertionPoint(self.Console.GetLastPosition())
        if self.linenum<=self.FileExec.GetNumberOfLines() and self.FileExec.GetNumberOfLines()>0:
            self.SetStatusText('Running')
            self.notebook1.SetSelection(1) # open notebook on Console
            dorun = self.AnalyzeNextLine(action='step')
            self.ReturnToPrompt(self.line)

    def ReadNextLine(self):
        self.line=self.FileExec.GetLineText(self.linenum)
        self.EdPos = self.EdPos + len(self.line) + 1
#    try to highlight syntax in Editor window: result not satisfying
#            if(self.linenum>0):
#                self.txtEditor.SetDefaultStyle(wxTextAttr(wxBLACK,wxWHITE))
#                self.txtEditor.SetInsertionPoint(self.EdPos-len(line)-1-len(self.oldline)-1)
#                self.txtEditor.Remove(self.EdPos-len(line)-1-len(self.oldline)-1,self.EdPos-len(line)-1)
#                self.txtEditor.WriteText(self.oldline+'\n')
#            self.txtEditor.SetDefaultStyle(wxTextAttr(wxBLACK,wxCYAN))
#            self.txtEditor.SetInsertionPoint(self.EdPos-len(line)-1)
#            self.txtEditor.Remove(self.EdPos-len(line)-1,self.EdPos)
#            self.txtEditor.WriteText(line+'\n')
        self.linenum = self.linenum+1
 
    def AnalyzeNextLine(self,action='cont'):
        dorun = true
        doraise = false
        docomment = false
        redo = true
        endsection = false
        while redo:
          self.ReadNextLine()
          if(string.lstrip(self.line) <> ''):
            firstword = string.split(self.line)[0]
            if(len(firstword)>=5):
                if(firstword[0:5]=='raise'): 
                    doraise=true
                else:
                    doraise=false
            if(not doraise):
                if(len(firstword)>=1):
                    if(firstword[0]=='#'): 
                        docomment=true
                if self.prefix is '... ':
                    if(self.line[0]<>' '):
                        if(len(self.line)>=4):
                            if(self.line[:4]<>'else' and self.line[:4]<>'elif'):
                                self.inter.push('\n')
                                endsection = true
                                redo = false
                                self.EdPos = self.EdPos - len(self.line) - 1
                                self.linenum = self.linenum-1
                                self.prefix=''
                        else:
                            self.inter.push('\n')
                            endsection = true
                            redo = false
                            self.EdPos = self.EdPos - len(self.line) - 1
                            self.linenum = self.linenum-1
                            self.prefix=''
                    self.oldline = self.line
                if docomment:
                    docomment = false
                    if endsection:
                        endsection = false
                    else: 
                        redo = true
                elif endsection:
                    endsection = false
                else:
                    redo = false
                    if(string.strip(self.line)=='winon()'):
                        self.OnWinonButton()
                    else:
                        self.Console.WriteText(self.prefix+self.line+'\n')
                        r=self.ConsolePanel.sendcommand(self.line,addlist=0)
                        if(r==1):
                          if action is 'next': redo=true
                          self.prefix='... '
                        else:
                          self.prefix=''
                          redo=false
            else:   
                self.prefix=''
                dorun = false
                redo = false
          else:
            if self.linenum<=self.FileExec.GetNumberOfLines():
                r=self.inter.push('\n')
                self.oldline = '#'
                redo = true
            else:
                redo = false
        return dorun

    def ReturnToPrompt(self,line):
        self.Console.WriteText('>>> ')
        self.ConsolePanel.CursorMin = self.Console.GetLastPosition()
        self.FileExec.ShowPosition(self.EdPos-len(line)-1) 
        self.statusBar1.SetStatusText(i=0,text="Ready")

    def OnMnudumpDump(self,event):
        dump()
    
    def OnMnudumpDumpAs(self,event):
        import DumpGUI
        self.DumpGUI = DumpGUI.wxDialog1(self)
        self.DumpGUI.Show(1)
    
    def OnMnudumpRestore(self,event):
        import RestoreGUI
        self.RestoreGUI = RestoreGUI.wxDialog1(self)
        self.RestoreGUI.Show(1)

    def GetFileName(self):
        dlg = wxFileDialog(self, "Choose a file", ".", "", "*.*", wxOPEN)
        try:
            if dlg.ShowModal() == wxID_OK:
                filename = dlg.GetPath()
                return filename
        finally:
            dlg.Destroy()
    
    def OnMnudumpRestart(self,event):
        restart(self.GetFileName())
        
    def OnMnuerrorchecksymmetryMenu(self, event):
        checksymmetry()

    def OnMnuerrorcheckibpushMenu(self, event):
        checkibpush()

    def OnMnuerrorcheckparticleloadMenu(self, event):
        checkparticleload()

    def OnMnuerrorcheckenvelopeMenu(self, event):
        checkenv()

    def OnMnuerrorcheckallMenu(self, event):
        errorcheck()

    def mnuPackageUpdate(self):
        currpkg = package()[0]
        if currpkg == 'w3d':
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGE3D,true)
        else:
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGE3D,false)
        if currpkg == 'wxy':
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGEXY,true)
        else:
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGEXY,false)
        if currpkg == 'env':
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGEENV,true)
        else:
          self.mnuPackage.Check(wxID_WARPRUNMNUPACKAGEENV,false)

    def OnMnupackage3dMenu(self, event):
        package('w3d')
        self.mnuPackageUpdate()

    def OnMnupackageXYMenu(self, event):
        package('wxy')
        self.mnuPackageUpdate()

    def OnMnupackageEnvMenu(self, event):
        package('env')
        self.mnuPackageUpdate()

    def OnWxframe1Idle(self, event):
        if self.isgistwindowon:
            self.HandleGistEvents()
            event.RequestMore(1)
    
    def OnWinonButton(self, event):
        if not self.isgistwindowon:
            winon()
            if sys.platform <> 'win32':
                self.HandleGistEvents()
                self.isgistwindowon = 1
        
    def OnFmaButton(self, event):
        fma()

    def OnHcpButton(self, event):
        hcp()

    def OnEnvButton(self, event):
        self.OnMnupackageEnvMenu(None)
        try:
            self.envelopeDialogOn = not self.envelopeDialogOn
        except AttributeError:
            self.envelopeDialogOn = 0
        if not self.envelopeDialogOn:
            self.envelopeDialogOn = 1
            self.envelopeDialog = EnvelopeGUI.EnvelopeGUI(self)
            try:
                self.envelopeDialog.Show(1)
            except:
                pass
        else:
            self.envelopeDialog.Destroy()
            self.envelopeDialogOn = 0

    def OnLatButton(self, event):
        self.LatticeDialog = LatticeGUI.LatticeGUI(self)
        try:
            self.LatticeDialog.Show(1)
        except:
            pass

    def OnDocButton(self, event):
        try:
            self.DocGUIOn = not self.DocGUIOn
        except AttributeError:
            self.DocGUIOn = 0
        if not self.DocGUIOn:
            self.DocGUIOn = 1
            self.DocGUI = DocGUI.DocGUI(self)
            try:
                self.DocGUI.Show(1)
            except:
                pass
        else:
            self.DocGUI.Destroy()
            self.DocGUIOn = 0

    def ImportGui(self,filename):
        pattern = string.split(os.path.basename(filename),'.')[0]
        exec('import '+pattern+'gui')
        exec('self.'+pattern+'Panel = '+pattern+'gui.'+pattern+'(self.notebook1)')
        exec("self.notebook1.AddPage(imageId=-1, page=self."+pattern+"Panel, text='"+pattern+",select=True')")
        self.notebook1.SetSelection(0) # open notebook on Editor
        
    def OutToConsole(self):
        sys.stdout = newstdout.newstdout(self.Console)
        sys.stderr = newstdout.newstdout(self.Console)

    def OutToMessageWindow(self):
        sys.stdout = newstdout.newstdout(self.MessageWindow)
        sys.stderr = newstdout.newstdout(self.MessageWindow)

    def OnNotebook1NotebookPageChanged(self, event):
        if event.GetSelection() is 1:
            self.OutToConsole()
        else:
            self.OutToMessageWindow()
        event.Skip()
            
        

#Boa:Frame:WarpRun

from wxPython.wx import *
from wxPython.stc import *
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
import wxDialog_proto
import pygistDialog
import WarpPanel
import gist
import sys
import code
import __main__
from warp import *
from errorcheck import *
sys.path=sys.path+[os.path.dirname(warp.__file__)+'/GUI/pype']
import pype

# for debugging purpose, output is not redirected in GUI if true
l_standard_out = 0

def create(parent):
    return WarpRun(parent)

[wxID_WARPRUN, wxID_WARPRUNBOOKMARK, wxID_WARPRUNCONT, wxID_WARPRUNDOC, 
 wxID_WARPRUNENV, wxID_WARPRUNFMA, wxID_WARPRUNHCP, wxID_WARPRUNLAT, 
 wxID_WARPRUNMESSAGEWINDOW, wxID_WARPRUNNEXT, wxID_WARPRUNNEXTBOOKMARK, 
 wxID_WARPRUNNOTEBOOK1, wxID_WARPRUNPANEL1, wxID_WARPRUNPREVBOOKMARK, 
 wxID_WARPRUNREDRAW, wxID_WARPRUNSEPARATE, wxID_WARPRUNSPLITTERWINDOW1, 
 wxID_WARPRUNSTART, wxID_WARPRUNSTATUSBAR1, wxID_WARPRUNSTEP, 
 wxID_WARPRUNTXTEDITOR, wxID_WARPRUNWINON, 
] = map(lambda _init_ctrls: wxNewId(), range(22))

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

    def _init_coll_menuBar1_Menus(self, parent):
        # generated method, don't edit

        parent.Append(menu=self.mnuFile, title='File')
        parent.Append(menu=self.mnuDump, title='Dump')
        parent.Append(menu=self.mnuHelp, title='Help')
        parent.Append(menu=self.mnuErrorCheck, title='ErrorCheck')
        parent.Append(menu=self.mnuPackage, title='Package')
        parent.Append(menu=self.mnuPalette, title='Palette')

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
              pos=wxPoint(522, 101), size=wxSize(628, 670),
              style=wxCLIP_CHILDREN | wxDEFAULT_FRAME_STYLE, title='WARP')
        self._init_utils()
        self.SetClientSize(wxSize(620, 646))
        self.SetMenuBar(self.menuBar1)
        self.SetAutoLayout(True)

        self.statusBar1 = wxStatusBar(id=wxID_WARPRUNSTATUSBAR1,
              name='statusBar1', parent=self, style=0)
        self.statusBar1.SetSize(wxSize(620, 19))
        self.statusBar1.SetPosition(wxPoint(0, 0))
        self._init_coll_statusBar1_Fields(self.statusBar1)
        self.SetStatusBar(self.statusBar1)

        self.panel1 = wxPanel(id=wxID_WARPRUNPANEL1, name='panel1', parent=self,
              pos=wxPoint(0, 0), size=wxSize(616, 24), style=wxTAB_TRAVERSAL)

        self.winon = wxButton(id=wxID_WARPRUNWINON, label='win', name='winon',
              parent=self.panel1, pos=wxPoint(0, 0), size=wxSize(40, 22),
              style=0)
        self.winon.SetBackgroundColour(wxColour(0, 0, 160))
        self.winon.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.winon.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.winon, wxID_WARPRUNWINON, self.OnWinonButton)

        self.fma = wxButton(id=wxID_WARPRUNFMA, label='fma', name='fma',
              parent=self.panel1, pos=wxPoint(40, 0), size=wxSize(40, 22),
              style=0)
        self.fma.SetBackgroundColour(wxColour(0, 0, 160))
        self.fma.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.fma.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.fma, wxID_WARPRUNFMA, self.OnFmaButton)

        self.hcp = wxButton(id=wxID_WARPRUNHCP, label='hcp', name='hcp',
              parent=self.panel1, pos=wxPoint(80, 0), size=wxSize(40, 22),
              style=0)
        self.hcp.SetBackgroundColour(wxColour(0, 0, 160))
        self.hcp.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.hcp.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.hcp, wxID_WARPRUNHCP, self.OnHcpButton)

        self.env = wxButton(id=wxID_WARPRUNENV, label='env', name='env',
              parent=self.panel1, pos=wxPoint(168, 0), size=wxSize(40, 22),
              style=0)
        self.env.SetBackgroundColour(wxColour(0, 128, 0))
        self.env.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.env.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.env, wxID_WARPRUNENV, self.OnEnvButton)

        self.lat = wxButton(id=wxID_WARPRUNLAT, label='lat', name='lat',
              parent=self.panel1, pos=wxPoint(208, 0), size=wxSize(40, 22),
              style=0)
        self.lat.SetBackgroundColour(wxColour(0, 128, 0))
        self.lat.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.lat.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.lat, wxID_WARPRUNLAT, self.OnLatButton)

        self.doc = wxButton(id=wxID_WARPRUNDOC, label='doc', name='doc',
              parent=self.panel1, pos=wxPoint(504, 0), size=wxSize(40, 22),
              style=0)
        self.doc.SetForegroundColour(wxColour(0, 0, 0))
        self.doc.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))
        self.doc.SetBackgroundColour(wxColour(255, 255, 128))
        EVT_BUTTON(self.doc, wxID_WARPRUNDOC, self.OnDocButton)

        self.Step = wxButton(id=wxID_WARPRUNSTEP, label='Step', name='Step',
              parent=self.panel1, pos=wxPoint(296, 0), size=wxSize(40, 22),
              style=0)
        self.Step.SetBackgroundColour(wxColour(128, 0, 64))
        self.Step.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.Step.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.Step, wxID_WARPRUNSTEP, self.OnStepButton)

        self.Next = wxButton(id=wxID_WARPRUNNEXT, label='Next', name='Next',
              parent=self.panel1, pos=wxPoint(336, 0), size=wxSize(40, 22),
              style=0)
        self.Next.SetBackgroundColour(wxColour(128, 0, 64))
        self.Next.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.Next.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.Next, wxID_WARPRUNNEXT, self.OnNextButton)

        self.Start = wxButton(id=wxID_WARPRUNSTART, label='Start', name='Start',
              parent=self.panel1, pos=wxPoint(256, 0), size=wxSize(40, 22),
              style=0)
        self.Start.SetBackgroundColour(wxColour(128, 0, 64))
        self.Start.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.Start.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.Start, wxID_WARPRUNSTART, self.OnStartButton)

        self.Cont = wxButton(id=wxID_WARPRUNCONT, label='Cont', name='Cont',
              parent=self.panel1, pos=wxPoint(376, 0), size=wxSize(40, 22),
              style=0)
        self.Cont.SetBackgroundColour(wxColour(128, 0, 64))
        self.Cont.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.Cont.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.Cont, wxID_WARPRUNCONT, self.OnContButton)

        self.separate = wxButton(id=wxID_WARPRUNSEPARATE, label='separate',
              name='separate', parent=self.panel1, pos=wxPoint(560, 0),
              size=wxSize(56, 22), style=0)
        self.separate.SetBackgroundColour(wxColour(128, 128, 128))
        self.separate.SetForegroundColour(wxColour(255, 255, 255))
        self.separate.SetConstraints(LayoutAnchors(self.separate, True, True,
              True, False))
        EVT_BUTTON(self.separate, wxID_WARPRUNSEPARATE, self.OnSeparateButton)

        self.redraw = wxButton(id=wxID_WARPRUNREDRAW, label='rdw',
              name='redraw', parent=self.panel1, pos=wxPoint(120, 0),
              size=wxSize(40, 22), style=0)
        self.redraw.SetBackgroundColour(wxColour(0, 0, 160))
        self.redraw.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxNORMAL, False,
              'MS Sans Serif'))
        self.redraw.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.redraw, wxID_WARPRUNREDRAW, self.OnRedrawButton)

        self.splitterWindow1 = wxSplitterWindow(id=wxID_WARPRUNSPLITTERWINDOW1,
              name='splitterWindow1', parent=self, point=wxPoint(0, 24),
              size=wxSize(616, 576), style=wxSP_3D)
        self.splitterWindow1.SetConstraints(LayoutAnchors(self.splitterWindow1,
              True, True, True, True))
        self.splitterWindow1.SetAutoLayout(True)
        self.splitterWindow1.SetMinimumPaneSize(30)

        self.MessageWindow = wxTextCtrl(id=wxID_WARPRUNMESSAGEWINDOW,
              name='MessageWindow', parent=self.splitterWindow1, pos=wxPoint(2,
              357), size=wxSize(612, 236),
              style=wxHSCROLL | wxVSCROLL | wxTE_READONLY | wxTE_MULTILINE,
              value='')
        self.MessageWindow.SetFont(wxFont(12, wxMODERN, wxNORMAL, wxNORMAL,
              false, ''))
        self.MessageWindow.SetBackgroundColour(wxColour(192, 192, 192))

        self.notebook1 = wxNotebook(id=wxID_WARPRUNNOTEBOOK1, name='notebook1',
              parent=self.splitterWindow1, pos=wxPoint(2, 2), size=wxSize(612,
              348), style=0)
        EVT_NOTEBOOK_PAGE_CHANGED(self.notebook1, wxID_WARPRUNNOTEBOOK1,
              self.OnNotebook1NotebookPageChanged)
        EVT_SIZE(self.notebook1, self.OnNotebook1Size)
        self.splitterWindow1.SplitHorizontally(self.notebook1,
              self.MessageWindow, 350)

        self.txtEditor = wxTextCtrl(id=wxID_WARPRUNTXTEDITOR, name='txtEditor',
              parent=self.notebook1, pos=wxPoint(0, 0), size=wxSize(604, 322),
              style=wxTE_MULTILINE, value='')
        self.txtEditor.SetToolTipString('Text Editor')

        self.BookMark = wxButton(id=wxID_WARPRUNBOOKMARK, label='->',
              name='BookMark', parent=self.panel1, pos=wxPoint(424, 0),
              size=wxSize(22, 22), style=0)
        self.BookMark.SetBackgroundColour(wxColour(155, 202, 230))
        self.BookMark.SetFont(wxFont(8, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.BookMark, wxID_WARPRUNBOOKMARK, self.OnBookmarkButton)

        self.PrevBookMark = wxButton(id=wxID_WARPRUNPREVBOOKMARK, label='<<',
              name='PrevBookMark', parent=self.panel1, pos=wxPoint(446, 0),
              size=wxSize(22, 22), style=0)
        self.PrevBookMark.SetBackgroundColour(wxColour(155, 202, 230))
        self.PrevBookMark.SetFont(wxFont(8, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.PrevBookMark, wxID_WARPRUNPREVBOOKMARK,
              self.OnPrevbookmarkButton)

        self.NextBookMark = wxButton(id=wxID_WARPRUNNEXTBOOKMARK, label='>>',
              name='NextBookMark', parent=self.panel1, pos=wxPoint(468, 0),
              size=wxSize(22, 22), style=0)
        self.NextBookMark.SetBackgroundColour(wxColour(155, 202, 230))
        self.NextBookMark.SetFont(wxFont(8, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.NextBookMark, wxID_WARPRUNNEXTBOOKMARK,
              self.OnNextbookmarkButton)

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
        # substitute default editor by pype
        self.notebook1.DeletePage(0)
        self.launch_pype()
        # start console
        self.ConsolePanel = ConsoleClass.ConsoleClass(parent=self.notebook1, inter=self.inter)
        self.Console = self.ConsolePanel.Console
        self.prefix = ''
        self.PplotsPanel = ParticlePlotsGUI.ParticlePlotsGUI(self.notebook1)
        self.panels = {}
        self.panels['Pzplots']  = self.show_GUI(PzplotsGUI,  'notebook','Pzplots')
        self.panels['Matching'] = self.show_GUI(MatchingGUI, 'notebook','Matching')
        self.panels['Gist']     = self.show_GUI(pygistDialog,'notebook','Gist')
        self.notebook1.SetSelection(0) # open notebook on Editor
        self.FileExecDialog = txtEditorDialog.txtEditorDialog(self)      
        self.FileExec = self.FileExecDialog.txtEditor  
        self.FileExec.Show(1)
        Palettes = ["earth","rainbow","gray","yarg","heat","ncar","cool","rainbowaf","stern","christmas"]
        for i in range(0,len(Palettes)):
            self.AddPalette(Palettes[i])
        self.gist_timer = wxPyTimer(self.HandleGistEvents)
        self.gist_timer.Start(100)
            
    def launch_pype(self):
        def GetKeyPress(evt):
            keycode = evt.GetKeyCode()
            keyname = pype.keyMap.get(keycode, None)
            modifiers = ""
            for mod, ch in [(evt.ControlDown(), 'Ctrl+'),
                            (evt.AltDown(),     'Alt+'),
                            (evt.ShiftDown(),   'Shift+')]:
                if mod:
                    modifiers += ch
            if keyname is None:
                if 27 < keycode < 256:
                    keyname = chr(keycode)
                else:
                    keyname = "(%s)unknown" % keycode
            return modifiers + keyname
        pype.GetKeyPress = GetKeyPress
        def menuAdd(root, menu, name, desc, funct, id, kind=wxITEM_NORMAL):
            a = wxMenuItem(menu, id, 'TEMPORARYNAME', desc, kind)
            menu.AppendItem(a)
            EVT_MENU(root.GetParent(), id, funct)

            ns, oacc = pype._spl(name)
            heir = pype.recmenu(pype.menuBar, id)[:-13] + ns
            if heir in pype.MENUPREF:
                name, acc = pype.MENUPREF[heir]
            else:
                if heir in pype.OLD_MENUPREF:
                    name, acc = pype.MENUPREF[heir] = pype.OLD_MENUPREF[heir]
                else:
                    name, acc = Mpype.ENUPREF[heir] = (ns, oacc)
                pype.MENULIST.append((heir, name, oacc, acc, kind in [wxITEM_NORMAL, wxITEM_CHECK]))

            if acc:
                pype.HOTKEY_TO_ID[acc] = id

            pype.menuBar.SetLabel(id, '%s\t%s'%(name, acc))
#        pype.menuAdd=menuAdd 

        if sys.executable[:6].lower() != 'python':
         import encodings.cp037
         import encodings.cp1006
         import encodings.cp1026
         import encodings.cp1140
         import encodings.cp1250
         import encodings.cp1251
         import encodings.cp1252
         import encodings.cp1253
         import encodings.cp1254
         import encodings.cp1255
         import encodings.cp1256
         import encodings.cp1257
         import encodings.cp1258
         import encodings.cp424
         import encodings.cp437
         import encodings.cp500
         import encodings.cp737
         import encodings.cp775
         import encodings.cp850
         import encodings.cp852
         import encodings.cp855
         import encodings.cp856
         import encodings.cp857
         import encodings.cp860
         import encodings.cp861
         import encodings.cp862
         import encodings.cp863
         import encodings.cp864
         import encodings.cp865
         import encodings.cp866
         import encodings.cp869
         import encodings.cp874
         import encodings.cp875
        if pype.VS[-1] == 'u':
         import encodings.ascii
         import encodings.utf_7
         import encodings.utf_8
         import encodings.utf_16
         import encodings.utf_16_be
         import encodings.utf_16_le
        opn=0
        if len(sys.argv)>1 and (sys.argv[1] == '--last'):opn=1
        pype.frame = pype.MainWindow(self, wxNewId(), "PyPE %s"%pype.VERSION, sys.argv[1+opn:])
        def resize_dummy(self,e=None):
            pass
        resize=pype.frame.OnResize
        pype.frame.OnResize=resize_dummy
        if self.FileName is not None:
            pype.frame.OnDrop([self.FileName])
        self.pype = pype
        self.menuBar1.Remove(0)
        pype.frame.SetSize(self.GetSize())
        sys.stderr.flush()
        sys.stdout.flush()
        self.SetPosition((0,0))
        pype.frame.SetPosition(self.GetPosition()+(10,10))
        if sys.platform <> 'win32':
          pype.frame.Show(0)
          panel = WarpPanel.panel(self.notebook1)
          self.notebook1.AddPage(imageId=-1, page=panel, select=True, text='Editor')
          pype.frame.menubar.Reparent(panel)
          pype.frame.control.Move(wxPoint(0,25))
          pype.frame.control.Reparent(panel)
          self.SetStatusBar(pype.frame.sb)
          pype.frame.sb.Reparent(self)
          self.statusBar1=pype.frame.sb
          def OnCpSize(evt,win=pype.frame.control):
             size = evt.GetSize()
             size.SetHeight(size.GetHeight()-25)
             win.SetSize(size)
          EVT_SIZE(panel,OnCpSize)
        pype.frame.OnResize=resize_dummy
       
        def testfollowpanel():
         panel = WarpPanel.panel(self.notebook1)
         self.notebook1.AddPage(imageId=-1, page=panel, select=True, text='Editor')
         def OnCpSize(evt,win=pype.frame):
            size = evt.GetSize()
            win.SetSize(size)
            pype.frame.Raise()
            pype.frame.SetFocus()
         EVT_SIZE(panel,OnCpSize)
         def getabspos(win):
            pos = win.GetPosition()
            try:
                pos+=getabspos(win.GetParent())
            except:
                pass
            return pos
         def OnCpMove(evt,win=pype.frame,panel=panel):
            pos = getabspos(panel)
            win.Move(pos+(2,40))
            pype.frame.Raise()
            pype.frame.SetFocus()
         EVT_MOVE(self,OnCpMove)
         # also need to add this into OnNotebook1NotebookPageChanged
#        if event.GetSelection() == 0:
#            self.pype.frame.Raise()
#            self.pype.frame.SetFocus()

        #bookmark support
        self.BOOKMARKNUMBER = pype.BOOKMARKNUMBER+1
        self.BOOKMARKSYMBOL = wxSTC_MARK_ARROW
        self.BOOKMARKMASK = 2**self.BOOKMARKNUMBER
        pype.frame.Old_newTab = pype.frame.newTab
        def newTab(d, fn, switch=0):
            pype.frame.Old_newTab(d, fn, switch)
            wnum, win = self.pype.frame.getNumWin()
            win.MarkerDefine(self.BOOKMARKNUMBER, self.BOOKMARKSYMBOL, 'red', 'red')
            win.SetFocus()
        pype.frame.newTab=newTab
        self.win = None
        
    def OnToggleBookmark (self, e):
            wnum, win = self.pype.frame.getNumWin(e)
            lineNo = win.GetCurrentLine()
            if win.MarkerGet(lineNo) & self.BOOKMARKMASK:
                win.MarkerDelete(lineNo, self.BOOKMARKNUMBER)
            else:
                win.MarkerAdd(lineNo, self.BOOKMARKNUMBER)
      
    def OnNextBookmark  (self, e):
            wnum, win = self.pype.frame.getNumWin(e)
            lineNo = win.GetCurrentLine()
            newLineNo = win.MarkerNext(lineNo + 1, self.BOOKMARKMASK)
            if newLineNo != -1:
                win.GotoLine(newLineNo)
            else:
                lineNo = win.GetLineCount()
                newLineNo = win.MarkerNext(0, self.BOOKMARKMASK)
                if newLineNo != -1:
                    win.GotoLine(newLineNo)
            win.EnsureVisible(win.GetCurrentLine())
            win.EnsureCaretVisible()
    
    def OnPreviousBookmark (self, e):
            wnum, win = self.pype.frame.getNumWin(e)
            lineNo = win.GetCurrentLine()
            newLineNo = win.MarkerPrevious(lineNo - 1, self.BOOKMARKMASK)
            if newLineNo != -1:
                win.GotoLine(newLineNo)
            else:
                lineNo = win.GetLineCount()
                newLineNo = win.MarkerPrevious(lineNo, self.BOOKMARKMASK)
                if newLineNo != -1:
                    win.GotoLine(newLineNo)
            win.EnsureVisible(win.GetCurrentLine())
            win.EnsureCaretVisible()

    def add_panel(self,panel,name,out='notebook'):
        if(self.panels.has_key(name)):return
        self.panels[name] = self.show_GUI(panel,out,name)
        self.OutToMessageWindow()

    def show_GUI(self,gui,winout,title):
        if(winout=='notebook'):
          panel = WarpPanel.panel(self.notebook1)
          self.notebook1.AddPage(imageId=-1, page=panel, select=True, text=title)
          panel.panel = gui.panel(panel)
          return {'panel':panel,'gui':gui,'winout':winout,'title':title}
        else:
          dialog = wxDialog_proto.wxDialog1(self,gui.panel,title)
          dialog.Show(1)
          return {'panel':dialog.panel,'gui':gui,'winout':winout,'title':title}
        
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
        if self.FileName is None:
            return
        self.Run()
        
    def Run(self):
        self.notebook1.SetSelection(0) # open notebook on Editor
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.Console)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.Console)
        self.statusBar1.SetStatusText(i=0,text="Executing file %s"%self.FileName)
        if(self.startrun and self.ConsolePanel.NoEntry):
            self.Console.Clear()
            startrun = 0
        self.OnContButton()

    def GetText(self):
        try:
          text = self.txtEditor.GetValue()
        except:
          if self.win is None:
            wnum, win = self.pype.frame.getNumWin()
            self.wnum = wnum
            self.win = win
            self.LineNo = 0
            self.lastline = self.win.LineFromPosition(self.win.GetLength())
          if self.LineNo>self.lastline:
              print '<End of file>'
              return ''
          newLineNo = self.win.MarkerNext(self.LineNo + 1, self.BOOKMARKMASK)
          if newLineNo==-1:
             newLineNo = self.lastline
          startpos = self.win.PositionFromLine(self.LineNo)
          endpos   = self.win.PositionFromLine(newLineNo+1)
          self.LineNo = newLineNo+1
          text = self.win.GetTextRange(startpos,endpos)
          self.win.GotoLine(newLineNo)
          self.win.EnsureVisible(self.win.GetCurrentLine())
          self.win.EnsureCaretVisible()
        return text

    def OnStartButton(self, event):
        self.FileExec.SetValue(self.GetText())
        self.SetInputsLimits()
        self.Run()

    def OnContButton(self, event=None):
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.Console)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.Console)
        self.SetStatusText('Running')
        dorun = true
        self.notebook1.SetSelection(0) # open notebook on Console
        if(self.linenum>self.FileExec.GetNumberOfLines()):
          self.FileExec.SetValue(self.GetText())
          self.linenum=0
        while(dorun and self.linenum<=self.FileExec.GetNumberOfLines()):
            dorun = self.AnalyzeNextLine(action='next')
            self.prefix='>>> '
        self.ReturnToPrompt(self.line)
        self.prefix=''

    def OnNextButton(self, event):
        if(self.linenum>self.FileExec.GetNumberOfLines()):
          self.FileExec.SetValue(self.GetText())
          self.linenum=0
        if self.linenum<=self.FileExec.GetNumberOfLines() and self.FileExec.GetNumberOfLines()>0:
            self.SetStatusText('Running')
            self.notebook1.SetSelection(0) # open notebook on Console
            dorun = self.AnalyzeNextLine(action='next')
            self.ReturnToPrompt(self.line)

    def OnStepButton(self, event):
#        self.Console.SetInsertionPoint(self.Console.GetLastPosition())
        if(self.linenum>self.FileExec.GetNumberOfLines()):
          self.FileExec.SetValue(self.GetText())
          self.linenum=0
        if self.linenum<=self.FileExec.GetNumberOfLines() and self.FileExec.GetNumberOfLines()>0:
            self.SetStatusText('Running')
            self.notebook1.SetSelection(0) # open notebook on Console
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
#                if self.prefix is '... ':
                elif self.prefix is '... ':
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

    def OnWinonButton(self, event=None):
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
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.Console)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.Console)

    def OutToMessageWindow(self):
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.MessageWindow)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.MessageWindow)

    def OnNotebook1NotebookPageChanged(self, event):
        if event.GetSelection() == 0:
            self.OutToConsole()
        else:
            self.OutToMessageWindow()
        if sys.platform == 'win32':event.Skip()

    def OnGistButton(self, event):
        import pygistDialog
        self.pygistDialog = pygistDialog.wxDialog1(self)
        self.pygistDialog.Show(1)
        if sys.platform <> 'win32':event.Skip()

    def OnSeparateButton(self, event):
        current = self.notebook1.GetPage(self.notebook1.GetSelection())
        for i in self.panels.keys():
            if self.panels[i]['panel'] == current:
                dialog = wxDialog_proto.wxDialog1(self,None,self.panels[i]['title'])
                self.panels[i]['panel'].panel.Reparent(dialog)
                dialog.panel=self.panels[i]['panel'].panel
                dialog.panel.Move(wxPoint(0,25))
                size = self.panels[i]['panel'].panel.GetSize()
                dialog.SetSize((size.GetWidth(),size.GetHeight()+25))
                dialog.Show(1)
                dialog.nbselection = self.notebook1.GetSelection()
                self.notebook1.GetPage(self.notebook1.GetSelection()).Hide()
                # the next 3 lines are needed on Windows
                dialog.Update()
                dialog.panel.Refresh()
                self.notebook1.Update()
                size = dialog.panel.GetSize()
                dialog.SetSize((size[0],size[1]+25))
        event.Skip()

    def OnRedrawButton(self, event):
        redraw()
        event.Skip()

    def OnNotebook1Size(self, event):
        try:
          self.pype.frame.OnResize(event)
        except:
          pass
        event.Skip()

    def OnBookmarkButton(self, event):
        self.OnToggleBookmark(event)
        event.Skip()

    def OnPrevbookmarkButton(self, event):
        self.OnPreviousBookmark(event)
        event.Skip()

    def OnNextbookmarkButton(self, event):
        self.OnNextBookmark(event)
        event.Skip()
            
        

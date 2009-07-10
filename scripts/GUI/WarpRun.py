#Boa:Frame:WarpRun

#from wx import py
from wx import *
from wx.stc import *
from wx.lib.anchors import LayoutAnchors
#from wxPython.wx import *
#from wxPython.stc import *
#from wxPython.lib.anchors import LayoutAnchors
import WarpGUIInfo
import ParticlePlotsGUI
import EnvelopeGUI
import LatticeGUI
import DocGUI
import MatchingGUI
import ConsoleClass
import PzplotsGUI
import ManualDialog
import txtEditorDialog
import newstdout
import wxDialog_proto
import pygistDialog
try:
  import Egun_like_gui
  l_egun=1
except:
  l_egun=0
try:
  import pyOpenDXDialog
  l_opendx=1
except:
  l_opendx=0
import WarpPanel
import gist
import sys
import os
import code
import __main__
from warp import *
from errorcheck import *
warp_path = os.path.dirname(warp.__file__)
if warp_path<>'':warp_path+='/'
sys.path=sys.path+[warp_path+'GUI/pype']
import pype

# for debugging purpose, output is not redirected in GUI if true
l_standard_out = 1
l_PyCrust = 0
l_pype = 0

def create(parent):
    return WarpRun(parent)

[wxID_WARPRUN, wxID_WARPRUNBOOKMARK, wxID_WARPRUNCONT, wxID_WARPRUNDOC, 
 wxID_WARPRUNENV, wxID_WARPRUNFMA, wxID_WARPRUNHCP, wxID_WARPRUNLAT, 
 wxID_WARPRUNMESSAGEWINDOW, wxID_WARPRUNNEXT, wxID_WARPRUNNEXTBOOKMARK, 
 wxID_WARPRUNNOTEBOOK1, wxID_WARPRUNPANEL1, wxID_WARPRUNPREVBOOKMARK, 
 wxID_WARPRUNREDRAW, wxID_WARPRUNSEPARATE, wxID_WARPRUNSPLITTERWINDOW1, 
 wxID_WARPRUNSTART, wxID_WARPRUNSTATUSBAR1, wxID_WARPRUNSTEP, 
 wxID_WARPRUNTXTEDITOR, wxID_WARPRUNWINON, 
] = map(lambda _init_ctrls: wx.NewId(), range(22))

[wxID_WARPRUNTOOLBAR2TOOLS0, wxID_WARPRUNTOOLBAR2TOOLS1, wxID_WARPRUNTOOLBAR2TOOLS2,
 wxID_WARPRUNTOOLBAR2TOOLS3] = map(lambda _init_coll_toolBar2_Tools: wx.NewId(), range(4))

[wxID_WARPRUNTOOLBAR1TOOLS0, wxID_WARPRUNTOOLBAR1TOOLS1, 
 wxID_WARPRUNTOOLBAR1TOOLS2, wxID_WARPRUNTOOLBAR1TOOLS3, 
 wxID_WARPRUNTOOLBAR1TOOLS4, wxID_WARPRUNTOOLBAR1TOOLS5, 
 wxID_WARPRUNTOOLBAR1TOOLS6, 
] = map(lambda _init_coll_toolBar1_Tools: wx.NewId(), range(7))

[wxID_WARPRUNMNUERRORCHECKCHECKALL, wxID_WARPRUNMNUERRORCHECKENVELOPE, 
 wxID_WARPRUNMNUERRORCHECKIBPUSH, wxID_WARPRUNMNUERRORCHECKPARTICLELOAD, 
 wxID_WARPRUNMNUERRORCHECKSYMMETRY, 
] = map(lambda _init_coll_mnuErrorCheck_Items: wx.NewId(), range(5))

[wxID_WARPRUNMNUPACKAGE3D, wxID_WARPRUNMNUPACKAGEENV, 
 wxID_WARPRUNMNUPACKAGEXY, 
] = map(lambda _init_coll_mnuPackage_Items: wx.NewId(), range(3))

[wxID_WARPRUNMNUFILEEXEC, wxID_WARPRUNMNUFILEEXIT, wxID_WARPRUNMNUFILEOPEN, 
 wxID_WARPRUNMNUFILEOPENEXEC, wxID_WARPRUNMNUFILESAVE, 
 wxID_WARPRUNMNUFILESAVEAS, 
] = map(lambda _init_coll_mnuFile_Items: wx.NewId(), range(6))

[wxID_WARPRUNMNUDUMPDUMP, wxID_WARPRUNMNUDUMPDUMPAS, 
 wxID_WARPRUNMNUDUMPRESTART, wxID_WARPRUNMNUDUMPRESTORE, 
] = map(lambda _init_coll_mnuDump_Items: wx.NewId(), range(4))

[wxID_WARPRUNMNUHELPMANUAL,wxID_WARPRUNMNUHELPSCRIPTS,wxID_WARPRUNMNUHELPSOURCE,
wxID_WARPRUNMNUHELPTUTORIAL,wxID_WARPRUNMNUHELPABOUT,
] = map(lambda _init_coll_mnuHelp_Items: wx.NewId(), range(5))

class WarpRun(wx.Frame):
    def _init_coll_menuBar1_Menus(self, parent):
        # generated method, don't edit

        parent.Append(menu=self.mnuFile, title='File')
        parent.Append(menu=self.mnuDump, title='Dump')
        parent.Append(menu=self.mnuErrorCheck, title='ErrorCheck')
        parent.Append(menu=self.mnuPackage, title='Package')
        parent.Append(menu=self.mnuPalette, title='Palette')
        parent.Append(menu=self.mnuHelp, title='Help')

    def _init_coll_mnuHelp_Items(self, parent):
        # generated method, don't edit

        parent.Append(help='browse manual', id=wxID_WARPRUNMNUHELPMANUAL,
              text='Manual', kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPMANUAL, self.OnMnuhelpManualMenu)
        parent.Append(help='browse scripts', id=wxID_WARPRUNMNUHELPSCRIPTS,
              text='Scripts', kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPSCRIPTS, self.OnMnuhelpScriptsMenu)
        parent.Append(help='browse source', id=wxID_WARPRUNMNUHELPSOURCE,
              text='Source', kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPSOURCE, self.OnMnuhelpSourceMenu)
        parent.Append(help='display tutorials', id=wxID_WARPRUNMNUHELPTUTORIAL,
              text='Tutorial', kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPTUTORIAL, self.OnMnuhelpTutorialMenu)
        parent.Append(help='Display info', id=wxID_WARPRUNMNUHELPABOUT,
              text='About', kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUHELPABOUT, self.OnMnuhelpAboutMenu)

    def _init_coll_mnuErrorCheck_Items(self, parent):
        # generated method, don't edit

        parent.Append(help='', id=wxID_WARPRUNMNUERRORCHECKSYMMETRY,
              text='Symmetry', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUERRORCHECKPARTICLELOAD,
              text='ParticleLoad', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUERRORCHECKENVELOPE,
              text='Envelope', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUERRORCHECKIBPUSH,
              text='Ibpush', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUERRORCHECKCHECKALL,
              text='CheckAll', kind=wx.ITEM_NORMAL)
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

        parent.Append(help='', id=wxID_WARPRUNMNUFILEOPEN, text='Open',
              kind=wx.ITEM_NORMAL)
        parent.Append(help='Opens and Executes file',
              id=wxID_WARPRUNMNUFILEOPENEXEC, text='Open/Execfile',
              kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUFILESAVE, text='Save',
              kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUFILESAVEAS,
              text='Save As', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUFILEEXEC,
              text='ExecFile', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUFILEEXIT, text='Exit',
              kind=wx.ITEM_NORMAL)
        EVT_MENU(self, wxID_WARPRUNMNUFILEOPEN, self.OnMnuOpenMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILESAVE, self.OnMnufileSaveMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILESAVEAS, self.OnMnufileSaveAsMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEEXIT, self.OnMnufileExitMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEEXEC, self.OnMnufileexecfileMenu)
        EVT_MENU(self, wxID_WARPRUNMNUFILEOPENEXEC, self.OnMnufileOpenExecMenu)

    def _init_coll_mnuPackage_Items(self, parent):
        # generated method, don't edit

        parent.Append(help='Select 3-D code', id=wxID_WARPRUNMNUPACKAGE3D,
              text='3-D', kind=wx.ITEM_CHECK)
        parent.Append(help='Select slice code',
              id=wxID_WARPRUNMNUPACKAGEXY, text='X-Y', kind=wx.ITEM_CHECK)
        parent.Append(help='Select envelope code',
              id=wxID_WARPRUNMNUPACKAGEENV, text='Envelope', kind=wx.ITEM_CHECK)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGE3D, self.OnMnupackage3dMenu)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGEXY, self.OnMnupackageXYMenu)
        EVT_MENU(self, wxID_WARPRUNMNUPACKAGEENV, self.OnMnupackageEnvMenu)

    def _init_coll_mnuDump_Items(self, parent):
        # generated method, don't edit

        parent.Append(help='', id=wxID_WARPRUNMNUDUMPRESTORE,
              text='Restore', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUDUMPRESTART,
              text='Restart', kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUDUMPDUMP, text='Dump',
              kind=wx.ITEM_NORMAL)
        parent.Append(help='', id=wxID_WARPRUNMNUDUMPDUMPAS,
              text='Dump As', kind=wx.ITEM_NORMAL)
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

        parent.SetStatusText(number=0, text='Status')

        parent.SetStatusWidths([-1])

    def _init_utils(self):
        # generated method, don't edit
        self.mnuFile = wx.Menu(title='File')
        self._init_coll_mnuFile_Items(self.mnuFile)

        self.mnuDump = wx.Menu(title='Dump')
        self._init_coll_mnuDump_Items(self.mnuDump)

        self.mnuHelp = wx.Menu(title='Help')
        self._init_coll_mnuHelp_Items(self.mnuHelp)

        self.menuBar1 = wx.MenuBar()

        self.mnuErrorCheck = wx.Menu(title='ErrorCheck')
        self._init_coll_mnuErrorCheck_Items(self.mnuErrorCheck)

        self.mnuPackage = wx.Menu(title='Package')
        self._init_coll_mnuPackage_Items(self.mnuPackage)

        self.mnuPalette = wx.Menu(title='Palette')

        self._init_coll_menuBar1_Menus(self.menuBar1)

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Frame.__init__(self, id=wxID_WARPRUN, name='WarpRun', parent=prnt,
              pos=wx.Point(522, 101), size=wx.Size(628, 670),
              style=wx.CLIP_CHILDREN | wx.DEFAULT_FRAME_STYLE, title='WARP')
        self._init_utils()
        self.SetClientSize(wx.Size(620, 646))
        self.SetMenuBar(self.menuBar1)
        self.SetAutoLayout(True)
        EVT_CLOSE(self,self.OnFrameClose)

        self.statusBar1 = wx.StatusBar(id=wxID_WARPRUNSTATUSBAR1,
              name='statusBar1', parent=self, style=0)
        self.statusBar1.SetSize(wx.Size(620, 19))
        self.statusBar1.SetPosition(wx.Point(0, 0))
        self._init_coll_statusBar1_Fields(self.statusBar1)
        self.SetStatusBar(self.statusBar1)

        self.panel1 = wx.Panel(id=wxID_WARPRUNPANEL1, name='panel1', parent=self,
              pos=wx.Point(0, 0), size=wx.Size(616, 24), style=wx.TAB_TRAVERSAL)

        self.winon = wx.Button(id=wxID_WARPRUNWINON, label='win', name='winon',
              parent=self.panel1, pos=wx.Point(0, 0), size=wx.Size(40, 22),
              style=0)
        self.winon.SetBackgroundColour(wx.Colour(0, 0, 160))
        self.winon.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.winon.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.winon, wxID_WARPRUNWINON, self.OnWinonButton)

        self.fma = wx.Button(id=wxID_WARPRUNFMA, label='fma', name='fma',
              parent=self.panel1, pos=wx.Point(40, 0), size=wx.Size(40, 22),
              style=0)
        self.fma.SetBackgroundColour(wx.Colour(0, 0, 160))
        self.fma.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.fma.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.fma, wxID_WARPRUNFMA, self.OnFmaButton)

        self.hcp = wx.Button(id=wxID_WARPRUNHCP, label='hcp', name='hcp',
              parent=self.panel1, pos=wx.Point(80, 0), size=wx.Size(40, 22),
              style=0)
        self.hcp.SetBackgroundColour(wx.Colour(0, 0, 160))
        self.hcp.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.hcp.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.hcp, wxID_WARPRUNHCP, self.OnHcpButton)

        self.env = wx.Button(id=wxID_WARPRUNENV, label='env', name='env',
              parent=self.panel1, pos=wx.Point(168, 0), size=wx.Size(40, 22),
              style=0)
        self.env.SetBackgroundColour(wx.Colour(0, 128, 0))
        self.env.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.env.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.env, wxID_WARPRUNENV, self.OnEnvButton)

        self.lat = wx.Button(id=wxID_WARPRUNLAT, label='lat', name='lat',
              parent=self.panel1, pos=wx.Point(208, 0), size=wx.Size(40, 22),
              style=0)
        self.lat.SetBackgroundColour(wx.Colour(0, 128, 0))
        self.lat.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.lat.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.lat, wxID_WARPRUNLAT, self.OnLatButton)

        self.doc = wx.Button(id=wxID_WARPRUNDOC, label='doc', name='doc',
              parent=self.panel1, pos=wx.Point(504, 0), size=wx.Size(40, 22),
              style=0)
        self.doc.SetForegroundColour(wx.Colour(0, 0, 0))
        self.doc.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))
        self.doc.SetBackgroundColour(wx.Colour(255, 255, 128))
        EVT_BUTTON(self.doc, wxID_WARPRUNDOC, self.OnDocButton)

        self.Step = wx.Button(id=wxID_WARPRUNSTEP, label='Step', name='Step',
              parent=self.panel1, pos=wx.Point(296, 0), size=wx.Size(40, 22),
              style=0)
        self.Step.SetBackgroundColour(wx.Colour(128, 0, 64))
        self.Step.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.Step.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.Step, wxID_WARPRUNSTEP, self.OnStepButton)

        self.Next = wx.Button(id=wxID_WARPRUNNEXT, label='Next', name='Next',
              parent=self.panel1, pos=wx.Point(336, 0), size=wx.Size(40, 22),
              style=0)
        self.Next.SetBackgroundColour(wx.Colour(128, 0, 64))
        self.Next.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.Next.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.Next, wxID_WARPRUNNEXT, self.OnNextButton)

        self.Start = wx.Button(id=wxID_WARPRUNSTART, label='Start', name='Start',
              parent=self.panel1, pos=wx.Point(256, 0), size=wx.Size(40, 22),
              style=0)
        self.Start.SetBackgroundColour(wx.Colour(128, 0, 64))
        self.Start.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.Start.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.Start, wxID_WARPRUNSTART, self.OnStartButton)

        self.Cont = wx.Button(id=wxID_WARPRUNCONT, label='Cont', name='Cont',
              parent=self.panel1, pos=wx.Point(376, 0), size=wx.Size(40, 22),
              style=0)
        self.Cont.SetBackgroundColour(wx.Colour(128, 0, 64))
        self.Cont.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.Cont.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.Cont, wxID_WARPRUNCONT, self.OnContButton)

        self.separate = wx.Button(id=wxID_WARPRUNSEPARATE, label='separate',
              name='separate', parent=self.panel1, pos=wx.Point(560, 0),
              size=wx.Size(56, 22), style=0)
        self.separate.SetBackgroundColour(wx.Colour(128, 128, 128))
        self.separate.SetForegroundColour(wx.Colour(255, 255, 255))
        self.separate.SetConstraints(LayoutAnchors(self.separate, True, True,
              True, False))
        EVT_BUTTON(self.separate, wxID_WARPRUNSEPARATE, self.OnSeparateButton)

        self.redraw = wx.Button(id=wxID_WARPRUNREDRAW, label='rdw',
              name='redraw', parent=self.panel1, pos=wx.Point(120, 0),
              size=wx.Size(40, 22), style=0)
        self.redraw.SetBackgroundColour(wx.Colour(0, 0, 160))
        self.redraw.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False,
              'MS Sans Serif'))
        self.redraw.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.redraw, wxID_WARPRUNREDRAW, self.OnRedrawButton)

        self.splitterWindow1 = wx.SplitterWindow(id=wxID_WARPRUNSPLITTERWINDOW1,
              name='splitterWindow1', parent=self, point=wx.Point(0, 24),
              size=wx.Size(616, 576), style=wx.SP_3D)
        self.splitterWindow1.SetConstraints(LayoutAnchors(self.splitterWindow1,
              True, True, True, True))
        self.splitterWindow1.SetAutoLayout(True)
        self.splitterWindow1.SetMinimumPaneSize(30)

#        self.MessageWindow = wx.TextCtrl(id=wxID_WARPRUNMESSAGEWINDOW,
#              name='MessageWindow', parent=self.splitterWindow1, pos=wx.Point(2,
#              357), size=wx.Size(580, 310),
#              style=wx.VSCROLL | wx.TE_READONLY | wx.TE_MULTILINE,
#              value='')
#        self.MessageWindow.SetFont(wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL,
#              false, ''))
#        self.MessageWindow.SetBackgroundColour('CADET BLUE')

        self.notebook1 = wx.Notebook(id=wxID_WARPRUNNOTEBOOK1, name='notebook1',
              parent=self.splitterWindow1, pos=wx.Point(2, 2), size=wx.Size(612,
              348), style=0)
        EVT_NOTEBOOK_PAGE_CHANGED(self.notebook1, wxID_WARPRUNNOTEBOOK1,
              self.OnNotebook1NotebookPageChanged)
        EVT_SIZE(self.notebook1, self.OnNotebook1Size)
#        self.splitterWindow1.SplitHorizontally(self.notebook1,
#              self.MessageWindow, 350)

        self.txtEditor = wx.TextCtrl(id=wxID_WARPRUNTXTEDITOR, name='txtEditor',
              parent=self.notebook1, pos=wx.Point(0, 0), size=wx.Size(604, 322),
              style=wx.TE_MULTILINE, value='')
        self.txtEditor.SetToolTipString('Text Editor')

        self.BookMark = wx.Button(id=wxID_WARPRUNBOOKMARK, label='->',
              name='BookMark', parent=self.panel1, pos=wx.Point(424, 0),
              size=wx.Size(22, 22), style=0)
        self.BookMark.SetBackgroundColour(wx.Colour(155, 202, 230))
        self.BookMark.SetFont(wx.Font(8, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.BookMark, wxID_WARPRUNBOOKMARK, self.OnBookmarkButton)

        self.PrevBookMark = wx.Button(id=wxID_WARPRUNPREVBOOKMARK, label='<<',
              name='PrevBookMark', parent=self.panel1, pos=wx.Point(446, 0),
              size=wx.Size(22, 22), style=0)
        self.PrevBookMark.SetBackgroundColour(wx.Colour(155, 202, 230))
        self.PrevBookMark.SetFont(wx.Font(8, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.PrevBookMark, wxID_WARPRUNPREVBOOKMARK,
              self.OnPrevbookmarkButton)

        self.NextBookMark = wx.Button(id=wxID_WARPRUNNEXTBOOKMARK, label='>>',
              name='NextBookMark', parent=self.panel1, pos=wx.Point(468, 0),
              size=wx.Size(22, 22), style=0)
        self.NextBookMark.SetBackgroundColour(wx.Colour(155, 202, 230))
        self.NextBookMark.SetFont(wx.Font(8, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))
        EVT_BUTTON(self.NextBookMark, wxID_WARPRUNNEXTBOOKMARK,
              self.OnNextbookmarkButton)

        self._init_coll_notebook1_Pages(self.notebook1)

    def __init__(self, parent):
        self._init_ctrls(parent)

    def init(self):
        self.FileName = None
        #EVT_UPDATE_UI(,self.mnuPackageUpdate)
        self.mnuPackageUpdate()
        self.isgistwindowon = 0
        self.oldline = '#'
        self.linenum = 0
        self.EdPos = 0
        self.startrun = 1
        self.panels = {}
        # substitute default editor by pype
        if l_pype:
          self.notebook1.DeletePage(0)
          self.launch_pype()
        self.prefix = ''
        # start console
        if l_PyCrust:
          def shortcuts(): 
            print """
* Key bindings:
Home              Go to the beginning of the command or line.
Shift+Home        Select to the beginning of the command or line.
Shift+End         Select to the end of the line.
End               Go to the end of the line.
Ctrl+C            Copy selected text, removing prompts.
Ctrl+Shift+C      Copy selected text, retaining prompts.
Ctrl+X            Cut selected text.
Ctrl+V            Paste from clipboard.
Ctrl+Shift+V      Paste and run multiple commands from clipboard.
Ctrl+Up Arrow     Retrieve Previous History item.
Alt+P             Retrieve Previous History item.
Ctrl+Down Arrow   Retrieve Next History item.
Alt+N             Retrieve Next History item.
Shift+Up Arrow    Insert Previous History item.
Shift+Down Arrow  Insert Next History item.
F8                Command-completion of History item.
                  (Type a few characters of a previous command and press F8.)
Ctrl+Enter        Insert new line into multiline command.
Ctrl+]            Increase font size.
Ctrl+[            Decrease font size.
Ctrl+=            Default font size.
"""

          __main__.shortcuts = shortcuts
          self.Crust = py.crust.Crust(self,-1,intro='For help on:\n - WARP      - type "warphelp()",\n - shortcuts - type "shorcuts()".\n\n')
          self.shell = self.Crust.shell
          self.crustnotebook = self.Crust.notebook
          self.Crust.notebook.Reparent(self)

          self.Crust.filling.Destroy() # too CPU intensive, might be replaced with WARP data exploration at sopme point
          self.Crust.dispatcherlisting.Destroy() # of no interest
          self.Crust.display.Destroy() # not sure what this window does anyway
          self.crustnotebook.SetPageText(3,'History')
          self.shell.Reparent(self.splitterWindow1)
          self.splitterWindow1.SplitHorizontally(self.notebook1,self.shell, 350)
#          self.splitterWindow1.ReplaceWindow(self.MessageWindow,self.shell)
#          self.shell.Reparent(self.splitterWindow1)
          self.Crust.Destroy()
#          self.MessageWindow.Destroy()
          if sys.platform <> 'cygwin':
              self.panels['Session']  = self.show_GUI(self.crustnotebook,  'notebook','Session','frame',True)
          self.Console = self.shell
          self.inter = self.shell.interp
          def SetInsertionPointEnd():  
            self.shell.SetCurrentPos(self.shell.GetTextLength())
          self.shell.SetInsertionPointEnd=SetInsertionPointEnd

          __main__.autocomp=self.AutoComp
          __main__.calltip=self.CallTip
          self.CallTip() # turns off calltip
        else:
          self.inter = code.InteractiveConsole(__main__.__dict__)
          self.ConsolePanel = ConsoleClass.ConsoleClass(parent=self.splitterWindow1,inter=self.inter)
          if l_pype:self.splitterWindow1.ReplaceWindow(self.MessageWindow,self.ConsolePanel)
          self.MessageWindow=self.ConsolePanel.Console
          self.Console = self.ConsolePanel.Console
# old:         self.ConsolePanel = ConsoleClass.ConsoleClass(parent=self.notebook1, inter=self.inter, in_notebook=1)
        self.PplotsPanel = ParticlePlotsGUI.ParticlePlotsGUI(self)
        self.panels['Pplots']  = self.show_GUI(self.PplotsPanel,  'notebook','Pplots','frame',True)
        self.panels['Pzplots']  = self.show_GUI(PzplotsGUI,  'notebook','Pzplots')
        self.panels['Matching'] = self.show_GUI(MatchingGUI, 'notebook','Matching')
        self.panels['Gist']     = self.show_GUI(pygistDialog,'notebook','Gist')
        if l_opendx:
          self.panels['OpenDX']     = self.show_GUI(pyOpenDXDialog,'notebook','OpenDX')
        if l_egun:
          self.panels['Egun']     = self.show_GUI(Egun_like_gui,'notebook','Egun')
#        self.notebook1.SetSelection(0) # open notebook on Editor
        self.notebook1.SetSelection(1) # open notebook on Pplots for now
        self.FileExecDialog = txtEditorDialog.txtEditorDialog(self)      
        self.FileExec = self.FileExecDialog.txtEditor  
        self.FileExec.Show(1)
        Palettes = ["earth","rainbow","gray","yarg","heat","ncar","cool","rainbowaf","stern","christmas"]
        for i in range(0,len(Palettes)):
            self.AddPalette(Palettes[i])
        self.gist_timer = wx.PyTimer(self.HandleGistEvents)
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
        def menuAdd(root, menu, name, desc, funct, id, kind=wx.ITEM_NORMAL):
            a = wx.MenuItem(menu, id, 'TEMPORARYNAME', desc, kind)
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
                pype.MENULIST.append((heir, name, oacc, acc, kind in [wx.ITEM_NORMAL, wx.ITEM_CHECK]))

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
        pype.frame = pype.MainWindow(self, wx.NewId(), "PyPE %s"%pype.VERSION, sys.argv[1+opn:])
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
        if sys.platform <> 'cygwin':
          pype.frame.Show(0)
 #         panel = WarpPanel.panel(self.notebook1)
 #         self.notebook1.AddPage(imageId=-1, page=panel, select=True, text='Editor')
          self.panels['Editor']  = self.show_GUI(None,'notebook','Editor','frame',True)
          panel = self.panels['Editor']['panel'].panel
          pype.frame.menubar.Reparent(panel)
          pype.frame.control.Move(wx.Point(0,25))
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
        self.BOOKMARKSYMBOL = wx.STC_MARK_ARROW
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

    def show_GUI(self,gui,winout='notebook',title='',type='dialog',newpanel=0):
        if(winout=='notebook'):
          panel = WarpPanel.panel(self.notebook1)
          self.notebook1.AddPage(imageId=-1, page=panel, select=True, text=title)
          if not newpanel: # then create
            panel.panel = gui.panel(panel)
          else:
            panel.panel = WarpPanel.panel(panel)
            if gui is not None:
              gui.Reparent(panel.panel)
              ref=panel.panel
              cs = wx.LayoutConstraints()
              cs.top.SameAs(ref,wx.Top)
              cs.bottom.SameAs(ref,wx.Bottom)
              cs.left.SameAs(ref,wx.Left)
              cs.right.SameAs(ref,wx.Right)
              gui.SetConstraints(cs)
              gui.SetAutoLayout(True)
            cs = wx.LayoutConstraints()
            ref=panel
            cs.top.SameAs(ref,wx.Top)
            cs.bottom.SameAs(ref,wx.Bottom)
            cs.left.SameAs(ref,wx.Left)
            cs.right.SameAs(ref,wx.Right)
            panel.panel.SetConstraints(cs)
            panel.panel.SetAutoLayout(True)
        else:
          if frame:
            frame = wxDialog_proto.wxFrame1(self,gui.panel,title)
            frame.Show(1)
            panel = frame.panel
          else:
            dialog = wxDialog_proto.wxDialog1(self,gui.panel,title)
            dialog.Show(1)
            panel = dialog.panel
        return {'panel':panel,'gui':gui,'winout':winout,'title':title,'type':type}
         
    def HandleGistEvents(self):
      refresh()

    def AddPalette(self,name):
        exec("[wxID_WARPRUNMNUPALLETTE"+name+"] = map(lambda _init_coll_mnuPalette_Items: wx.NewId(), range(1))")
        exec("self.mnuPalette.Append(help='', id=wxID_WARPRUNMNUPALLETTE"+name+", text='"+name+"',kind=wx.ITEM_NORMAL)")
        exec("def OnMnuPalette"+name+"(event):palette('"+name+".gp')")
        exec("EVT_MENU(self, wxID_WARPRUNMNUPALLETTE"+name+", OnMnuPalette"+name+")")

    def OnMnuhelpAboutMenu(self, event):
        dlg = WarpGUIInfo.WarpGUIInfo(self)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def OnMnuhelpManualMenu(self, event):
        try: 
           self.ManualHtml
        except:
           self.panels['Manual']  = self.show_GUI(ManualDialog,  'notebook','Manual','frame') 
           self.ManualHtml = self.panels['Manual']['panel'].panel.html
        self.ManualHtml.GoHome('manual/manual')

    def OnMnuhelpScriptsMenu(self, event):
        try: 
           self.ScriptsHtml
        except:
           self.panels['Scripts']  = self.show_GUI(ManualDialog,  'notebook','Scripts','frame') 
           self.ScriptsHtml = self.panels['Scripts']['panel'].panel.html
        self.ScriptsHtml.GoHome('scripts/doc/index')

    def OnMnuhelpSourceMenu(self, event):
        try: 
           self.SourceHtml
        except:
           self.panels['Source']  = self.show_GUI(ManualDialog,  'notebook','Source','frame') 
           self.SourceHtml = self.panels['Source']['panel'].panel.html
        self.SourceHtml.GoHome('source')

    def OnMnuhelpTutorialMenu(self, event):
        try: 
           self.TutorialHtml
        except:
           self.panels['Tutorial']  = self.show_GUI(ManualDialog,  'notebook','Tutorial','frame') 
           self.TutorialHtml = self.panels['Tutorial']['panel'].panel.html
        self.TutorialHtml.GoHome('tutorials')

    def OnMnufileOpenExecMenu(self, event):
        self.OnMnuOpenMenu(event)
        self.OnMnufileexecfileMenu(event)

    def OnMnuOpenMenu(self, event):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", 
              "PYTHON files (*.py)|*.py|ALL files (*.*)|*.*", wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
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
        dlg = wx.FileDialog(self, "Choose a file", ".", "", 
              "PYTHON files (*.py)|*.py|ALL files (*.*)|*.*", wx.SAVE|wx.OVERWRITE_PROMPT)
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
        self.statusBar1.SetStatusText(number=0,text="Executing file %s"%self.FileName)
        execfile(self.FileName)
        self.statusBar1.SetStatusText(number=0,text="Finished executing file %s"%self.FileName)

    def OnMnufileexecfileMenu(self, event):
        if self.FileName is None:
            return
        self.Run()
        
    def Run(self):
        self.notebook1.SetSelection(0) # open notebook on Editor
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.Console)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.Console)
        self.statusBar1.SetStatusText(number=0,text="Executing file %s"%self.FileName)
	if not l_PyCrust:
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
#                self.txtEditor.SetDefaultStyle(wx.TextAttr(wx.BLACK,wx.WHITE))
#                self.txtEditor.SetInsertionPoint(self.EdPos-len(line)-1-len(self.oldline)-1)
#                self.txtEditor.Remove(self.EdPos-len(line)-1-len(self.oldline)-1,self.EdPos-len(line)-1)
#                self.txtEditor.WriteText(self.oldline+'\n')
#            self.txtEditor.SetDefaultStyle(wx.TextAttr(wx.BLACK,wx.CYAN))
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
                if not docomment and self.prefix is '... ':
                    if(self.line[0]<>' '):
                        if(len(self.line)>=4):
                            if(self.line[:4]<>'else' and self.line[:4]<>'elif'):
                                self.inter.push(os.linesep)
                                endsection = true
                                redo = false
                                self.EdPos = self.EdPos - len(self.line) - 1
                                self.linenum = self.linenum-1
                                self.prefix=''
                        else:
                            self.inter.push(os.linesep)
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
			if l_PyCrust:
			  self.shell.SetCurrentPos(self.shell.GetTextLength())
			  self.shell.write(self.line+os.linesep)
			  more=self.inter.push(self.line)
			else:
                          self.Console.WriteText(self.prefix+self.line+os.linesep)
                          more=self.ConsolePanel.sendcommand(self.line,addlist=0)
                        if(more):
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
                self.oldline = '#'
                redo = true
            else:
                redo = false
        return dorun

    def ReturnToPrompt(self,line):
	if l_PyCrust:
	  self.shell.prompt()
        else:    
	  self.Console.WriteText('>>> ')
          self.ConsolePanel.CursorMin = self.Console.GetLastPosition()
#	self.shell.AddText('>>> ')
#	self.shell.DocumentEnd()
	self.FileExec.ShowPosition(self.EdPos-len(line)-1) 
	self.statusBar1.SetStatusText(number=0,text="Ready")

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
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*.*", wx.OPEN)
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
            if sys.platform <> 'cygwin':
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
        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.Console)
        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.Console)
#        if(not l_standard_out): sys.stdout = newstdout.newstdout(self.MessageWindow)
#        if(not l_standard_out): sys.stderr = newstdout.newstdout(self.MessageWindow)

    def OnNotebook1NotebookPageChanged(self, event):
        if event.GetSelection() == 0:
	  try:
            self.OutToConsole()
	  except:
	    pass
        else:
            self.OutToMessageWindow()
        if sys.platform == 'cygwin':event.Skip()

    def OnGistButton(self, event):
        import pygistDialog
        self.pygistDialog = pygistDialog.wxDialog1(self)
        self.pygistDialog.Show(1)
        if sys.platform <> 'cygwin':event.Skip()

    def OnSeparateButton(self, event):
        current = self.notebook1.GetPage(self.notebook1.GetSelection())
        for i in self.panels.keys():
            if self.panels[i]['panel'] == current:
                if self.panels[i]['type'] == 'frame':
                  dialog = wxDialog_proto.wxFrame1(self,None,self.panels[i]['title'])
                else:
                  dialog = wxDialog_proto.wxDialog1(self,None,self.panels[i]['title'])
                self.panels[i]['panel'].panel.Reparent(dialog)
                dialog.panel=self.panels[i]['panel'].panel
                dialog.panel.Move(wx.Point(0,25))
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
                if self.panels[i]['type'] == 'frame':
                  cs = wx.LayoutConstraints()
                  ref=dialog
                  cs.top.SameAs(ref,wx.Top,25)
                  cs.bottom.SameAs(ref,wx.Bottom)
                  cs.left.SameAs(ref,wx.Left)
                  cs.right.SameAs(ref,wx.Right)
                  dialog.panel.SetConstraints(cs)
                  dialog.panel.SetAutoLayout(True)
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

    def CallTip(self):
        self.shell.autoCallTip=1-self.shell.autoCallTip
        if self.shell.autoCallTip:
          print 'Auto call-tip turned on.'
        else:
          print 'Auto call-tip turned off.'

    def AutoComp(self):
        self.shell.autoComplete=1-self.shell.autoComplete
        if self.shell.autoComplete:
          print 'Auto completion turned on.'
        else:
          print 'Auto completion turned off.'

    def OnFrameClose(self,event):
      self.Hide()
      __main__.wgui.initialized=False
      __main__.wgui.closed=True
      __main__.wgui.ExitMainLoop()
      event.Skip()


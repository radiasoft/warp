#Boa:Frame:wxFrame1

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
import WarpGUIInfo
import ParticlePlotsGUI
import EnvelopeGUI
import LatticeGUI
import DocGUI

from warp import *
from errorcheck import *

def create(parent):
    return wxFrame1(parent)

[wxID_WXFRAME1, wxID_WXFRAME1STATUSBAR1, wxID_WXFRAME1TOOLBAR1, 
 wxID_WXFRAME1TXTEDITOR, 
] = map(lambda _init_ctrls: wxNewId(), range(4))

[wxID_WXFRAME1TOOLBAR2TOOLS0, wxID_WXFRAME1TOOLBAR2TOOLS1, wxID_WXFRAME1TOOLBAR2TOOLS2, wxID_WXFRAME1TOOLBAR2TOOLS3] = map(lambda _init_coll_toolBar2_Tools: wxNewId(), range(4))

[wxID_WXFRAME1TOOLBAR1TOOLS0, wxID_WXFRAME1TOOLBAR1TOOLS1, 
 wxID_WXFRAME1TOOLBAR1TOOLS2, wxID_WXFRAME1TOOLBAR1TOOLS3, 
 wxID_WXFRAME1TOOLBAR1TOOLS4, wxID_WXFRAME1TOOLBAR1TOOLS5, 
 wxID_WXFRAME1TOOLBAR1TOOLS6, 
] = map(lambda _init_coll_toolBar1_Tools: wxNewId(), range(7))

[wxID_WXFRAME1MNUERRORCHECKCHECKALL, wxID_WXFRAME1MNUERRORCHECKENVELOPE, 
 wxID_WXFRAME1MNUERRORCHECKIBPUSH, wxID_WXFRAME1MNUERRORCHECKPARTICLELOAD, 
 wxID_WXFRAME1MNUERRORCHECKSYMMETRY, 
] = map(lambda _init_coll_mnuErrorCheck_Items: wxNewId(), range(5))

[wxID_WXFRAME1MNUPACKAGE3D, wxID_WXFRAME1MNUPACKAGEENV, 
 wxID_WXFRAME1MNUPACKAGEXY, 
] = map(lambda _init_coll_mnuPackage_Items: wxNewId(), range(3))

[wxID_WXFRAME1MNUFILEEXEC, wxID_WXFRAME1MNUFILEEXIT, wxID_WXFRAME1MNUFILEOPEN, 
 wxID_WXFRAME1MNUFILEOPENEXEC, wxID_WXFRAME1MNUFILESAVE, 
 wxID_WXFRAME1MNUFILESAVEAS, 
] = map(lambda _init_coll_mnuFile_Items: wxNewId(), range(6))

[wxID_WXFRAME1MENU1ITEMS0] = map(lambda _init_coll_menu1_Items: wxNewId(), range(1))

[wxID_WXFRAME1MNUHELPABOUT] = map(lambda _init_coll_mnuHelp_Items: wxNewId(), range(1))

class wxFrame1(wxFrame):
    def _init_coll_mnuPackage_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='Select 3-D code',
              id=wxID_WXFRAME1MNUPACKAGE3D, item='3-D', kind=wxITEM_CHECK)
        parent.Append(helpString='Select slice code',
              id=wxID_WXFRAME1MNUPACKAGEXY, item='X-Y', kind=wxITEM_CHECK)
        parent.Append(helpString='Select envelope code',
              id=wxID_WXFRAME1MNUPACKAGEENV, item='Envelope',
              kind=wxITEM_CHECK)
        EVT_MENU(self, wxID_WXFRAME1MNUPACKAGE3D, self.OnMnupackage3dMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUPACKAGEXY, self.OnMnupackageXYMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUPACKAGEENV, self.OnMnupackageEnvMenu)

    def _init_coll_mnuHelp_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='Display info', id=wxID_WXFRAME1MNUHELPABOUT,
              item='About', kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WXFRAME1MNUHELPABOUT, self.OnMnuhelpAboutMenu)

    def _init_coll_mnuFile_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='', id=wxID_WXFRAME1MNUFILEOPEN, item='Open',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='Opens and Executes file',
              id=wxID_WXFRAME1MNUFILEOPENEXEC, item='Open/Execfile',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUFILESAVE, item='Save',
              kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUFILESAVEAS,
              item='Save As', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUFILEEXEC,
              item='ExecFile', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUFILEEXIT, item='Exit',
              kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WXFRAME1MNUFILEOPEN, self.OnMnuOpenMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUFILESAVE, self.OnMnufileSaveMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUFILESAVEAS, self.OnMnufileSaveAsMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUFILEEXIT, self.OnMnufileExitMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUFILEEXEC, self.OnMnufileexecfileMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUFILEOPENEXEC, self.OnMnufileOpenExecMenu)

    def _init_coll_menuBar1_Menus(self, parent):
        # generated method, don't edit

        parent.Append(menu=self.mnuFile, title='File')
        parent.Append(menu=self.mnuHelp, title='Help')
        parent.Append(menu=self.mnuErrorCheck, title='ErrorCheck')
        parent.Append(menu=self.mnuPackage, title='Package')

    def _init_coll_mnuErrorCheck_Items(self, parent):
        # generated method, don't edit

        parent.Append(helpString='', id=wxID_WXFRAME1MNUERRORCHECKSYMMETRY,
              item='Symmetry', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUERRORCHECKPARTICLELOAD,
              item='ParticleLoad', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUERRORCHECKENVELOPE,
              item='Envelope', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUERRORCHECKIBPUSH,
              item='Ibpush', kind=wxITEM_NORMAL)
        parent.Append(helpString='', id=wxID_WXFRAME1MNUERRORCHECKCHECKALL,
              item='CheckAll', kind=wxITEM_NORMAL)
        EVT_MENU(self, wxID_WXFRAME1MNUERRORCHECKSYMMETRY,
              self.OnMnuerrorchecksymmetryMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUERRORCHECKPARTICLELOAD,
              self.OnMnuerrorcheckparticleloadMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUERRORCHECKENVELOPE,
              self.OnMnuerrorcheckenvelopeMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUERRORCHECKIBPUSH,
              self.OnMnuerrorcheckibpushMenu)
        EVT_MENU(self, wxID_WXFRAME1MNUERRORCHECKCHECKALL,
              self.OnMnuerrorcheckallMenu)

    def _init_coll_statusBar1_Fields(self, parent):
        # generated method, don't edit
        parent.SetFieldsCount(1)

        parent.SetStatusText(i=0, text='Status')

        parent.SetStatusWidths([-1])

    def _init_coll_toolBar1_Tools(self, parent):
        # generated method, don't edit

        parent.AddTool(bitmap=wxBitmap('images/winon.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS0, isToggle=false, longHelpString='',
              pushedBitmap=wxNullBitmap, shortHelpString='Open Gist Window')
        parent.AddTool(bitmap=wxBitmap('images/fma.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS1, isToggle=false, longHelpString='',
              pushedBitmap=wxNullBitmap, shortHelpString='Frame Advance')
        parent.AddTool(bitmap=wxBitmap('images/hcp.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS2, isToggle=false, longHelpString='',
              pushedBitmap=wxNullBitmap, shortHelpString='Hard Copy')
        parent.AddTool(bitmap=wxBitmap('images/pplots.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS3, isToggle=false, longHelpString='',
              pushedBitmap=wxNullBitmap, shortHelpString='Particle Plots')
        parent.AddTool(bitmap=wxBitmap('images/envelope.bmp',
              wxBITMAP_TYPE_BMP), id=wxID_WXFRAME1TOOLBAR1TOOLS4,
              isToggle=false, longHelpString='', pushedBitmap=wxNullBitmap,
              shortHelpString='Envelope Solver')
        parent.AddTool(bitmap=wxBitmap('images/lattice.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS5, isToggle=false,
              longHelpString='Opens lattice editting dialog',
              pushedBitmap=wxNullBitmap, shortHelpString='Lattice')
        parent.AddTool(bitmap=wxBitmap('images/doc.bmp', wxBITMAP_TYPE_BMP),
              id=wxID_WXFRAME1TOOLBAR1TOOLS6, isToggle=false,
              longHelpString='Opens doc window', pushedBitmap=wxNullBitmap,
              shortHelpString='Doc window')
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS0,
              self.OnToolbar1OpenWindowTool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS1, self.OnToolbar1fmaTool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS2, self.OnToolbar1hcpTool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS3, self.OnToolbar1pplotsTool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS4, self.OnToolbar1envelopeTool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS5, self.OnToolbar1tools5Tool)
        EVT_TOOL(self, wxID_WXFRAME1TOOLBAR1TOOLS6, self.OnToolbar1tools6Tool)

        parent.Realize()

    def _init_utils(self):
        # generated method, don't edit
        self.mnuFile = wxMenu(title='File')
        self._init_coll_mnuFile_Items(self.mnuFile)

        self.mnuHelp = wxMenu(title='Help')
        self._init_coll_mnuHelp_Items(self.mnuHelp)

        self.menuBar1 = wxMenuBar()

        self.mnuErrorCheck = wxMenu(title='ErrorCheck')
        self._init_coll_mnuErrorCheck_Items(self.mnuErrorCheck)

        self.mnuPackage = wxMenu(title='Package')
        self._init_coll_mnuPackage_Items(self.mnuPackage)

        self._init_coll_menuBar1_Menus(self.menuBar1)

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxFrame.__init__(self, id=wxID_WXFRAME1, name='', parent=prnt,
              pos=wxPoint(186, 358), size=wxSize(949, 535),
              style=wxDEFAULT_FRAME_STYLE, title='WARP')
        self._init_utils()
        self.SetClientSize(wxSize(949, 508))
        self.SetMenuBar(self.menuBar1)
        EVT_IDLE(self, self.OnWxframe1Idle)

        self.statusBar1 = wxStatusBar(id=wxID_WXFRAME1STATUSBAR1,
              name='statusBar1', parent=self, style=0)
        self.statusBar1.SetSize(wxSize(877, 25))
        self.statusBar1.SetPosition(wxPoint(0, 596))
        self._init_coll_statusBar1_Fields(self.statusBar1)
        self.SetStatusBar(self.statusBar1)

        self.toolBar1 = wxToolBar(id=wxID_WXFRAME1TOOLBAR1, name='toolBar1',
              parent=self, pos=wxPoint(0, 27), size=wxSize(266, 38),
              style=wxTB_3DBUTTONS | wxSIMPLE_BORDER | wxTB_HORIZONTAL | wxNO_BORDER)
        self.toolBar1.SetConstraints(LayoutAnchors(self.toolBar1, true, true,
              false, false))
        self.toolBar1.SetBackgroundColour(wxColour(100, 100, 221))
        self._init_coll_toolBar1_Tools(self.toolBar1)
        self.SetToolBar(self.toolBar1)

        self.txtEditor = wxTextCtrl(id=wxID_WXFRAME1TXTEDITOR, name='txtEditor',
              parent=self, pos=wxPoint(0, 0), size=wxSize(949, 445),
              style=wxTE_MULTILINE, value='')
        self.txtEditor.SetToolTipString('Text Editor')

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.FileName = None
        #EVT_UPDATE_UI(,self.mnuPackageUpdate)
        self.mnuPackageUpdate()
        self.isgistwindowon = 0

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
        dlg = wxFileDialog(self, "Choose a file", ".", "", "*.*", wxOPEN)
        try:
            if dlg.ShowModal() == wxID_OK:
                filename = dlg.GetPath()
                self.txtEditor.LoadFile(filename)
                self.FileName = filename
        finally:
            dlg.Destroy()

    def OnMnufileSaveMenu(self, event):
        if self.FileName is None:
            return OnMnufileSaveAsMenu(event)
        else:
            self.txtEditor.SaveFile(self.FileName)

    def OnMnufileSaveAsMenu(self, event):
        dlg = wxFileDialog(self, "Choose a file", ".", "", "*.*", wxOPEN)
        try:
            if dlg.ShowModal() == wxID_OK:
                filename = dlg.GetPath()
                self.txtEditor.SaveFile(filename)
                self.FileName = filename
        finally:
            dlg.Destroy()
        event.Skip()

    def OnMnufileExitMenu(self, event):
        self.Close()

    def OnMnufileexecfileMenu(self, event):
        if self.FileName is None:
            OnMnufileSaveAsMenu(event)
        self.statusBar1.SetStatusText(i=0,text="Executing file %s"%self.FileName)
        execfile(self.FileName)
        self.statusBar1.SetStatusText(i=0,text="Finished executing file %s"%self.FileName)

    def OnToolbar1OpenWindowTool(self, event):
        winon()
#       dlg = wxMessageDialog(self, 'Press enter in the starting terminal window',
#         'Press enter', wxOK | wxICON_INFORMATION)
#       try:
#           dlg.ShowModal()
#       finally:
#           dlg.Destroy()
        ygdispatch()
        self.isgistwindowon = 1

    def OnToolbar1hcpTool(self, event):
        hcp()

    def OnToolbar1fmaTool(self, event):
        fma()

    def OnToolbar1pplotsTool(self, event):
        try:
            self.pplotsDialogOn = not self.pplotsDialogOn
        except AttributeError:
            self.pplotsDialogOn = 0
        if not self.pplotsDialogOn:
            self.pplotsDialogOn = 1
            self.pplotsDialog = ParticlePlotsGUI.ParticlePlotsGUI(self)
            try:
                self.pplotsDialog.Show(1)
            except:
                pass
        else:
            self.pplotsDialog.Destroy()
            self.pplotsDialogOn = 0

    def OnToolbar1envelopeTool(self, event):
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

    def OnToolbar1tools6Tool(self, event):
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
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGE3D,true)
        else:
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGE3D,false)
        if currpkg == 'wxy':
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGEXY,true)
        else:
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGEXY,false)
        if currpkg == 'env':
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGEENV,true)
        else:
          self.mnuPackage.Check(wxID_WXFRAME1MNUPACKAGEENV,false)

    def OnMnupackage3dMenu(self, event):
        package('w3d')
        self.mnuPackageUpdate()

    def OnMnupackageXYMenu(self, event):
        package('wxy')
        self.mnuPackageUpdate()

    def OnMnupackageEnvMenu(self, event):
        package('env')
        self.mnuPackageUpdate()

    def OnToolbar1tools5Tool(self, event):
        self.LatticeDialog = LatticeGUI.LatticeGUI(self)
        try:
            self.LatticeDialog.Show(1)
        except:
            pass

    def OnWxframe1Idle(self, event):
        if self.isgistwindowon:
          ygdispatch()
          event.RequestMore(1)



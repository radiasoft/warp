#Boa:FramePanel:panel

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from wxPython.grid import *
from warp import *
from pyOpenDX import *
import __main__

[wxID_PANEL, wxID_PANELFILENAME, wxID_PANELFORMAT, wxID_PANELMODE, 
 wxID_PANELSAVE, wxID_PANELSTATICBOX1, wxID_PANELSTATICTEXT1, 
 wxID_PANELSTATICTEXT2, wxID_PANELSTATICTEXT3, 
] = map(lambda _init_ctrls: wxNewId(), range(9))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='OpenDX', parent=prnt,
              pos=wxPoint(370, 238), size=wxSize(344, 232),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(336, 208))
        EVT_ENTER_WINDOW(self, self.OnPanelEnterWindow)

        self.staticText1 = wxStaticText(id=wxID_PANELSTATICTEXT1, label='Mode',
              name='staticText1', parent=self, pos=wxPoint(16, 20),
              size=wxSize(44, 16), style=0)

        self.Mode = wxChoice(choices=['rotate', 'translate', 'zoom'],
              id=wxID_PANELMODE, name='Mode', parent=self, pos=wxPoint(76, 16),
              size=wxSize(90, 21), style=0, validator=wxDefaultValidator)
        self.Mode.SetLabel('')
        self.Mode.Show(True)
        EVT_CHOICE(self.Mode, wxID_PANELMODE, self.OnModeChoice)

        self.staticText2 = wxStaticText(id=wxID_PANELSTATICTEXT2,
              label='Filename', name='staticText2', parent=self, pos=wxPoint(16,
              48), size=wxSize(50, 13), style=0)

        self.Format = wxChoice(choices=['eps', 'eps grey', 'miff', 
              'ps', 'ps grey', 'rgb', 'r+g+b', 'tiff', 'yuv'],
              id=wxID_PANELFORMAT, name='Format', parent=self, pos=wxPoint(76,
              68), size=wxSize(90, 21), style=0, validator=wxDefaultValidator)

        self.Filename = wxTextCtrl(id=wxID_PANELFILENAME, name='Filename',
              parent=self, pos=wxPoint(76, 42), size=wxSize(90, 21), style=0,
              value='Image')

        self.staticText3 = wxStaticText(id=wxID_PANELSTATICTEXT3,
              label='Format', name='staticText3', parent=self, pos=wxPoint(16,
              72), size=wxSize(40, 13), style=0)

        self.Save = wxButton(id=wxID_PANELSAVE, label='Save', name='Save',
              parent=self, pos=wxPoint(76, 96), size=wxSize(90, 21), style=0)
        EVT_BUTTON(self.Save, wxID_PANELSAVE, self.OnSaveButton)

        self.staticBox1 = wxStaticBox(id=wxID_PANELSTATICBOX1, label='Display',
              name='staticBox1', parent=self, pos=wxPoint(8, 0),
              size=wxSize(168, 128), style=0)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wxPoint(0,0))

    def OnModeChoice(self, event):
        __main__.interactor = self.Mode.GetSelection()
        event.Skip()

    def OnPanelEnterWindow(self, event):
        self.Mode.SetSelection(__main__.interactor)
        event.Skip()

    def OnSaveButton(self, event):
        self.filename = str(self.Filename.GetValue())
        self.format = str(self.Format.GetStringSelection())
        if self.format=='jpeg':self.format='ImageMagick supported format gamma=2.2 compression=JPEG quality=90'
        DXWriteImage(self.filename,__main__.wgui.dxobject,__main__.dcamera,None,self.format)
        event.Skip()


#Boa:FramePanel:panel

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from wxPython.grid import *
from warp import *
from pyOpenDX import *
import __main__

[wxID_PANEL, wxID_PANELFILENAME, wxID_PANELFORMAT, wxID_PANELMODE, 
 wxID_PANELRENDERING, wxID_PANELSAVE, wxID_PANELSTATICBOX1, 
 wxID_PANELSTATICBOX2, wxID_PANELSTATICTEXT1, wxID_PANELSTATICTEXT2, 
 wxID_PANELSTATICTEXT3, wxID_PANELSTATICTEXT4, 
] = map(lambda _init_ctrls: wxNewId(), range(12))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='OpenDX', parent=prnt,
              pos=wxPoint(393, 242), size=wxSize(365, 127),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(357, 103))
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
              label='Filename', name='staticText2', parent=self,
              pos=wxPoint(192, 20), size=wxSize(50, 13), style=0)

        self.Format = wxChoice(choices=['eps', 'eps grey', 'miff', 'ps',
              'ps grey', 'rgb', 'r+g+b', 'tiff', 'yuv'], id=wxID_PANELFORMAT,
              name='Format', parent=self, pos=wxPoint(244, 40), size=wxSize(90,
              21), style=0, validator=wxDefaultValidator)

        self.Filename = wxTextCtrl(id=wxID_PANELFILENAME, name='Filename',
              parent=self, pos=wxPoint(244, 16), size=wxSize(90, 21), style=0,
              value='Image')

        self.staticText3 = wxStaticText(id=wxID_PANELSTATICTEXT3,
              label='Format', name='staticText3', parent=self, pos=wxPoint(192,
              44), size=wxSize(40, 13), style=0)

        self.Save = wxButton(id=wxID_PANELSAVE, label='Save', name='Save',
              parent=self, pos=wxPoint(244, 64), size=wxSize(90, 21), style=0)
        EVT_BUTTON(self.Save, wxID_PANELSAVE, self.OnSaveButton)

        self.staticBox1 = wxStaticBox(id=wxID_PANELSTATICBOX1, label='Display',
              name='staticBox1', parent=self, pos=wxPoint(8, 0),
              size=wxSize(168, 96), style=0)

        self.staticText4 = wxStaticText(id=wxID_PANELSTATICTEXT4,
              label='Rendering', name='staticText4', parent=self,
              pos=wxPoint(16, 44), size=wxSize(55, 13), style=0)

        self.Rendering = wxChoice(choices=['software', 'hardware'],
              id=wxID_PANELRENDERING, name='Rendering', parent=self,
              pos=wxPoint(76, 40), size=wxSize(90, 21), style=0,
              validator=wxDefaultValidator)
        EVT_CHOICE(self.Rendering, wxID_PANELRENDERING, self.OnRenderingChoice)

        self.staticBox2 = wxStaticBox(id=wxID_PANELSTATICBOX2, label='Output',
              name='staticBox2', parent=self, pos=wxPoint(184, 0),
              size=wxSize(164, 96), style=0)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wxPoint(0,0))

    def OnModeChoice(self, event):
        __main__.interactor = self.Mode.GetSelection()
        event.Skip()

    def OnPanelEnterWindow(self, event):
        self.Mode.SetSelection(__main__.interactor)
        self.Rendering.SetSelection(__main__.l_hardware_acceleration)
        event.Skip()

    def OnSaveButton(self, event):
        self.filename = str(self.Filename.GetValue())
        self.format = str(self.Format.GetStringSelection())
        if self.format=='jpeg':self.format='ImageMagick supported format gamma=2.2 compression=JPEG quality=90'
        DXWriteImage(self.filename,__main__.wgui.dxobject,__main__.dcamera,None,self.format)
        event.Skip()

    def OnRenderingChoice(self, event):
        DXRendering(str(self.Rendering.GetStringSelection()))
        event.Skip()


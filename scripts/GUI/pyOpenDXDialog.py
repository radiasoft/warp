#Boa:FramePanel:panel

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from wxPython.grid import *
from warp import *
from pyOpenDX import *
import __main__

[wxID_PANEL, wxID_PANELMODE, wxID_PANELSTATICTEXT1
] = map(lambda _init_ctrls: wxNewId(), range(3))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='OpenDX', parent=prnt,
              pos=wxPoint(498, 297), size=wxSize(344, 232),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(336, 208))
        EVT_ENTER_WINDOW(self, self.OnPanelEnterWindow)

        self.staticText1 = wxStaticText(id=wxID_PANELSTATICTEXT1, label='Mode',
              name='staticText1', parent=self, pos=wxPoint(16, 20),
              size=wxSize(44, 16), style=0)

        self.Mode = wxChoice(choices=['rotate', 'translate','zoom'], id=wxID_PANELMODE,
              name='Mode', parent=self, pos=wxPoint(68, 16), size=wxSize(90,
              21), style=0, validator=wxDefaultValidator)
        self.Mode.SetLabel('')
        self.Mode.Show(True)
        EVT_CHOICE(self.Mode, wxID_PANELMODE, self.OnModeChoice)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wxPoint(0,0))

    def OnModeChoice(self, event):
        __main__.interactor = self.Mode.GetSelection()
        event.Skip()

    def OnPanelEnterWindow(self, event):
        self.Mode.SetSelection(__main__.interactor)
        event.Skip()


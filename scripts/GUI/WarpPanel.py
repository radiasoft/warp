#Boa:FramePanel:panel

from wxPython.wx import *
from warp import *

[wxID_PANEL] = map(lambda _init_ctrls: wxNewId(), range(1))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='', parent=prnt,
              pos=wxPoint(498, 297), size=wxSize(604, 339),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(596, 315))

    def __init__(self, parent):
        self._init_ctrls(parent)

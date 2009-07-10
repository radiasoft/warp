#Boa:FramePanel:panel

from wx import *
from warp import *

[wxID_PANEL] = map(lambda _init_ctrls: wx.NewId(), range(1))

class panel(wx.Panel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='', parent=prnt,
              pos=wx.Point(498, 297), size=wx.Size(604, 339),
              style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(596, 315))
        self.SetAutoLayout(True)

    def __init__(self, parent):
        self._init_ctrls(parent)

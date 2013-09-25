#Boa:MiniFrame:DocGUI

from wx import *
from warp import *
import sys

def create(parent):
    return DocGUI(parent)

[wxID_DOCGUI, wxID_DOCGUIDOCLABEL, wxID_DOCGUIDOCNAME, wxID_DOCGUIDOCTEXT, wxID_DOCGUIGETNAME] = map(lambda _init_ctrls: wx.NewId(), range(5))

class DocGUI(wx.MiniFrame):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wx.MiniFrame.__init__(self, id = wxID_DOCGUI, name = 'DocGUI', parent = prnt, pos = wx.Point(395, 276), size = wx.Size(625, 338), style = wx.DEFAULT_FRAME_STYLE, title = 'Doc')
        self._init_utils()
        self.SetClientSize(wx.Size(625, 338))
        self.SetToolTipString('Doc')

        self.GetName = wx.TextCtrl(id = wxID_DOCGUIGETNAME, name = 'GetName', parent = self, pos = wx.Point(8, 8), size = wx.Size(184, 22), style = wx.TE_PROCESS_ENTER | wx.TE_PROCESS_TAB, value = '')
        self.GetName.SetToolTipString('Enter name')
        EVT_TEXT_ENTER(self.GetName, wxID_DOCGUIGETNAME, self.OnGetnameTextEnter)

        self.DocText = wx.TextCtrl(id = wxID_DOCGUIDOCTEXT, name = 'DocText', parent = self, pos = wx.Point(8, 64), size = wx.Size(608, 264), style = wx.TE_READONLY | wx.TE_MULTILINE, value = '')
        self.DocText.SetToolTipString('')

        self.DocLabel = wx.StaticText(id = wxID_DOCGUIDOCLABEL, label = 'Enter name to get documentation', name = 'DocLabel', parent = self, pos = wx.Point(200, 10), size = wx.Size(183, 16), style = 0)

        self.DocName = wx.StaticText(id = wxID_DOCGUIDOCNAME, label = 'doc name', name = 'DocName', parent = self, pos = wx.Point(8, 38), size = wx.Size(184, 16), style = 0)

    def __init__(self, parent):
        self._init_ctrls(parent)

    def OnGetnameTextEnter(self, event):
        text = self.GetName.GetValue()
        self.DocName.SetLabel(text)
        self.GetName.SetValue('')
        d = doc(text,printit=0)
        self.DocText.SetValue(d)

#Boa:Dialog:WarpGUIInfo

#from wxPython.wx import *
from wx import *

def create(parent):
    return WarpGUIInfo(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1BUTTON1, wxID_WXDIALOG1STATICTEXT1] = map(lambda _init_ctrls: wx.NewId(), range(3))

class WarpGUIInfo(wx.Dialog):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wx.Dialog.__init__(self, id = wxID_WXDIALOG1, name = '', parent = prnt, pos = wx.Point(329, 258), size = wx.Size(233, 116), style = wxDEFAULT_DIALOG_STYLE, title = 'About Notebook')
        self._init_utils()
        self.SetClientSize(wx.Size(233, 116))

        self.staticText1 = wx.StaticText(id = wxID_WXDIALOG1STATICTEXT1, label = 'Notebook text editor', name = 'staticText1', parent = self, pos = wx.Point(16, 16), size = wx.Size(220, 26), style = wx.ALIGN_CENTRE)
        self.staticText1.SetFont(wx.wFont(20, wx.SWISS, wx.NORMAL, wx.NORMAL, false, ''))

        self.button1 = wx.Button(id = wxID_WXDIALOG1BUTTON1, label = 'Close', name = 'button1', parent = self, pos = wx.Point(24, 56), size = wx.Size(80, 22), style = 0)
        EVT_BUTTON(self.button1, wxID_WXDIALOG1BUTTON1, self.OnButton1Button)

    def __init__(self, parent):
        self._init_ctrls(parent)

    def OnButton1Button(self, event):
        self.Close()

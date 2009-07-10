#Boa:Dialog:txtEditorDialog

from wx import *

def create(parent):
    return txtEditorDialog(parent)

[wxID_TXTEDITORDIALOG, wxID_TXTEDITORDIALOGTXTEDITOR, 
] = map(lambda _init_ctrls: wx.NewId(), range(2))

class txtEditorDialog(wx.Dialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Dialog.__init__(self, id=wxID_TXTEDITORDIALOG, name='txtEditorDialog',
              parent=prnt, pos=wx.Point(318, 222), size=wx.Size(529, 336),
              style=wx.DEFAULT_DIALOG_STYLE, title='wx.Dialog1')
        self._init_utils()
        self.SetClientSize(wx.Size(529, 336))

        self.txtEditor = wx.TextCtrl(id=wxID_TXTEDITORDIALOGTXTEDITOR,
              name='txtEditor', parent=self, pos=wx.Point(0, 0), size=wx.Size(529,
              336), style=wx.TE_MULTILINE, value='')

    def __init__(self, parent):
        self._init_ctrls(parent)

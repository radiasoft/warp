#Boa:Dialog:txtEditorDialog

from wxPython.wx import *

def create(parent):
    return txtEditorDialog(parent)

[wxID_TXTEDITORDIALOG, wxID_TXTEDITORDIALOGTXTEDITOR, 
] = map(lambda _init_ctrls: wxNewId(), range(2))

class txtEditorDialog(wxDialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxDialog.__init__(self, id=wxID_TXTEDITORDIALOG, name='txtEditorDialog',
              parent=prnt, pos=wxPoint(318, 222), size=wxSize(529, 336),
              style=wxDEFAULT_DIALOG_STYLE, title='wxDialog1')
        self._init_utils()
        self.SetClientSize(wxSize(529, 336))

        self.txtEditor = wxTextCtrl(id=wxID_TXTEDITORDIALOGTXTEDITOR,
              name='txtEditor', parent=self, pos=wxPoint(0, 0), size=wxSize(529,
              336), style=wxTE_MULTILINE, value='')

    def __init__(self, parent):
        self._init_ctrls(parent)

#Boa:Dialog:WarpGUIInfo

from wxPython.wx import *

def create(parent):
    return WarpGUIInfo(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1BUTTON1, wxID_WXDIALOG1STATICTEXT1] = map(lambda _init_ctrls: wxNewId(), range(3))

class WarpGUIInfo(wxDialog):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wxDialog.__init__(self, id = wxID_WXDIALOG1, name = '', parent = prnt, pos = wxPoint(329, 258), size = wxSize(233, 116), style = wxDEFAULT_DIALOG_STYLE, title = 'About Notebook')
        self._init_utils()
        self.SetClientSize(wxSize(233, 116))

        self.staticText1 = wxStaticText(id = wxID_WXDIALOG1STATICTEXT1, label = 'Notebook text editor', name = 'staticText1', parent = self, pos = wxPoint(16, 16), size = wxSize(220, 26), style = wxALIGN_CENTRE)
        self.staticText1.SetFont(wxFont(20, wxSWISS, wxNORMAL, wxNORMAL, false, ''))

        self.button1 = wxButton(id = wxID_WXDIALOG1BUTTON1, label = 'Close', name = 'button1', parent = self, pos = wxPoint(24, 56), size = wxSize(80, 22), style = 0)
        EVT_BUTTON(self.button1, wxID_WXDIALOG1BUTTON1, self.OnButton1Button)

    def __init__(self, parent):
        self._init_ctrls(parent)

    def OnButton1Button(self, event):
        self.Close()

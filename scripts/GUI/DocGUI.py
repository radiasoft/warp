#Boa:MiniFrame:DocGUI

from wxPython.wx import *
from warp import *
import sys

def create(parent):
    return DocGUI(parent)

[wxID_DOCGUI, wxID_DOCGUIDOCLABEL, wxID_DOCGUIDOCNAME, wxID_DOCGUIDOCTEXT, wxID_DOCGUIGETNAME] = map(lambda _init_ctrls: wxNewId(), range(5))

class DocGUI(wxMiniFrame):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wxMiniFrame.__init__(self, id = wxID_DOCGUI, name = 'DocGUI', parent = prnt, pos = wxPoint(395, 276), size = wxSize(625, 338), style = wxDEFAULT_FRAME_STYLE, title = 'Doc')
        self._init_utils()
        self.SetClientSize(wxSize(625, 338))
        self.SetToolTipString('Doc')

        self.GetName = wxTextCtrl(id = wxID_DOCGUIGETNAME, name = 'GetName', parent = self, pos = wxPoint(8, 8), size = wxSize(184, 22), style = wxTE_PROCESS_ENTER | wxTE_PROCESS_TAB, value = '')
        self.GetName.SetToolTipString('Enter name')
        EVT_TEXT_ENTER(self.GetName, wxID_DOCGUIGETNAME, self.OnGetnameTextEnter)

        self.DocText = wxTextCtrl(id = wxID_DOCGUIDOCTEXT, name = 'DocText', parent = self, pos = wxPoint(8, 64), size = wxSize(608, 264), style = wxTE_READONLY | wxTE_MULTILINE, value = '')
        self.DocText.SetToolTipString('')

        self.DocLabel = wxStaticText(id = wxID_DOCGUIDOCLABEL, label = 'Enter name to get documentation', name = 'DocLabel', parent = self, pos = wxPoint(200, 10), size = wxSize(183, 16), style = 0)

        self.DocName = wxStaticText(id = wxID_DOCGUIDOCNAME, label = 'doc name', name = 'DocName', parent = self, pos = wxPoint(8, 38), size = wxSize(184, 16), style = 0)

    def __init__(self, parent):
        self._init_ctrls(parent)

    def OnGetnameTextEnter(self, event):
        text = self.GetName.GetValue()
        self.DocName.SetLabel(text)
        self.GetName.SetValue('')
        d = doc(text,printit=0)
        self.DocText.SetValue(d)
        

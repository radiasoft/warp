#Boa:Dialog:wxDialog1

from wxPython.wx import *
import WarpPanel

def create(parent):
    return wxDialog1(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1TONOTEBOOK, wxID_WXDIALOG1WINDOW1, 
] = map(lambda _init_ctrls: wxNewId(), range(3))

class wxDialog1(wxDialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxDialog.__init__(self, id=wxID_WXDIALOG1, name='', parent=prnt,
              pos=wxPoint(419, 165), size=wxSize(638, 434),
              style=wxDEFAULT_DIALOG_STYLE, title='wxDialog1')
        self._init_utils()
        self.SetClientSize(wxSize(630, 410))

        self.window1 = wxWindow(id=wxID_WXDIALOG1WINDOW1, name='window1',
              parent=self, pos=wxPoint(0, 24), size=wxSize(640, 416), style=0)

        self.tonotebook = wxButton(id=wxID_WXDIALOG1TONOTEBOOK,
              label='to notebook', name='tonotebook', parent=self,
              pos=wxPoint(0, 0), size=wxSize(75, 23), style=0)
        self.tonotebook.SetBackgroundColour(wxColour(128, 128, 128))
        self.tonotebook.SetForegroundColour(wxColour(255, 255, 255))
        EVT_BUTTON(self.tonotebook, wxID_WXDIALOG1TONOTEBOOK,
              self.OnTonotebookButton)

    def __init__(self, parent, child, title):
        self._init_ctrls(parent)
        self.SetTitle(title)
        self.parent = parent
        self.title = title
        if child is None: return
        self.panel = child(self.window1)
        self.child = child

    def OnTonotebookButton(self, event):
        self.panel.Reparent(self.parent.notebook1.GetPage(self.nbselection))
        self.parent.notebook1.GetPage(self.nbselection).Show(1)
        self.parent.notebook1.SetSelection(self.nbselection)
        self.panel.Move(wxPoint(0,0))
        self.Destroy()
        event.Skip()

#Boa:Dialog:wxDialog1

from wxPython.wx import *
from warp import *

def create(parent):
    return wxDialog1(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1BROWSER, wxID_WXDIALOG1FILENAME, 
 wxID_WXDIALOG1RESTORE, wxID_WXDIALOG1SKIP, wxID_WXDIALOG1STATICTEXT1, 
 wxID_WXDIALOG1STATICTEXT2, wxID_WXDIALOG1TEXT3, wxID_WXDIALOG1VARSUFFIX, 
 wxID_WXDIALOG1VERBOSE, 
] = map(lambda _init_ctrls: wxNewId(), range(10))

class wxDialog1(wxDialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxDialog.__init__(self, id=wxID_WXDIALOG1, name='', parent=prnt,
              pos=wxPoint(394, 138), size=wxSize(320, 103),
              style=wxDEFAULT_DIALOG_STYLE, title='Restore')
        self._init_utils()
        self.SetClientSize(wxSize(320, 103))
        self.SetFont(wxFont(12, wxSWISS, wxNORMAL, wxNORMAL, false, ''))
        self.SetForegroundColour(wxColour(0, 0, 0))

        self.staticText1 = wxStaticText(id=wxID_WXDIALOG1STATICTEXT1,
              label='Filename:', name='staticText1', parent=self, pos=wxPoint(4,
              10), size=wxSize(54, 16), style=0)

        self.Filename = wxTextCtrl(id=wxID_WXDIALOG1FILENAME, name='Filename',
              parent=self, pos=wxPoint(60, 8), size=wxSize(160, 22),
              style=wxTE_PROCESS_ENTER, value='')

        self.Browser = wxButton(id=wxID_WXDIALOG1BROWSER, label='Browse',
              name='Browser', parent=self, pos=wxPoint(230, 8), size=wxSize(80,
              22), style=0)
        EVT_BUTTON(self.Browser, wxID_WXDIALOG1BROWSER, self.OnBrowserButton)

        self.staticText2 = wxStaticText(id=wxID_WXDIALOG1STATICTEXT2,
              label='Skip:', name='staticText2', parent=self, pos=wxPoint(4,
              40), size=wxSize(27, 16), style=0)

        self.Skip = wxTextCtrl(id=wxID_WXDIALOG1SKIP, name='Skip', parent=self,
              pos=wxPoint(60, 38), size=wxSize(100, 22),
              style=wxTE_PROCESS_ENTER, value='')

        self.Verbose = wxToggleButton(id=wxID_WXDIALOG1VERBOSE, label='Verbose',
              name='Verbose', parent=self, pos=wxPoint(170, 38),
              size=wxSize(140, 22), style=0)
        self.Verbose.SetValue(false)
        EVT_TOGGLEBUTTON(self.Verbose, wxID_WXDIALOG1VERBOSE,
              self.OnVerboseTogglebutton)

        self.Text3 = wxStaticText(id=wxID_WXDIALOG1TEXT3, label='Suffix:',
              name='Text3', parent=self, pos=wxPoint(4, 70), size=wxSize(33,
              16), style=0)

        self.Varsuffix = wxTextCtrl(id=wxID_WXDIALOG1VARSUFFIX,
              name='Varsuffix', parent=self, pos=wxPoint(60, 68),
              size=wxSize(100, 22), style=wxTE_PROCESS_ENTER, value='')

        self.Restore = wxButton(id=wxID_WXDIALOG1RESTORE, label='RESTORE',
              name='Restore', parent=self, pos=wxPoint(170, 68),
              size=wxSize(140, 22), style=0)
        self.Restore.SetFont(wxFont(12, wxSWISS, wxNORMAL, wxBOLD, false, ''))
        self.Restore.SetForegroundColour(wxColour(230, 0, 0))
        EVT_BUTTON(self.Restore, wxID_WXDIALOG1RESTORE, self.OnRestoreButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.verbose = 0
        self.suffix = ''
        self.skip = []

    def OnBrowserButton(self, event):
        dlg = wxFileDialog(self, "Choose a file", ".", "", 
              "DUMP files (*.dump)|*.dump|PDB files (*.pdb)|*.pdb|ALL files (*.*)|*.*", wxOPEN)
        try:
            if dlg.ShowModal() == wxID_OK:
                self.filename = dlg.GetPath()
                self.Filename.SetValue(self.filename)
        finally:
            dlg.Destroy()

    def OnVerboseTogglebutton(self, event):
        if(self.Verbose.GetValue()):
            self.verbose = 1
        else:
            self.verbose = 0

    def OnRestoreButton(self, event):
        self.varsuffix = self.Varsuffix.GetValue()
        self.filename = self.Filename.GetValue()
        self.skip = [self.Skip.GetValue()]
        if self.varsuffix is '': self.varsuffix=None
        restore(filename=self.filename,verbose=self.verbose,skip=self.skip,varsuffix=self.varsuffix)
        self.Destroy()

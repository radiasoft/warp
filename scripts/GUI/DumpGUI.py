#Boa:Dialog:wxDialog1

from wxPython.wx import *
from warp import *

def create(parent):
    return wxDialog1(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1ATTRIBUTES, wxID_WXDIALOG1BROWSER, 
 wxID_WXDIALOG1DUMP, wxID_WXDIALOG1FILENAME, wxID_WXDIALOG1PYVARS, 
 wxID_WXDIALOG1STATICTEXT1, wxID_WXDIALOG1STATICTEXT2, wxID_WXDIALOG1TEXT3, 
 wxID_WXDIALOG1VARSUFFIX, 
] = map(lambda _init_ctrls: wxNewId(), range(10))

class wxDialog1(wxDialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxDialog.__init__(self, id=wxID_WXDIALOG1, name='', parent=prnt,
              pos=wxPoint(394, 138), size=wxSize(320, 103),
              style=wxDEFAULT_DIALOG_STYLE, title='Dump')
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
              label='Attributes:', name='staticText2', parent=self,
              pos=wxPoint(4, 40), size=wxSize(55, 16), style=0)

        self.Attributes = wxTextCtrl(id=wxID_WXDIALOG1ATTRIBUTES,
              name='Attributes', parent=self, pos=wxPoint(60, 38),
              size=wxSize(100, 22), style=wxTE_PROCESS_ENTER, value='dump')

        self.Pyvars = wxToggleButton(id=wxID_WXDIALOG1PYVARS,
              label='Save python variables', name='Pyvars', parent=self,
              pos=wxPoint(170, 38), size=wxSize(140, 22), style=0)
        self.Pyvars.SetValue(true)
        EVT_TOGGLEBUTTON(self.Pyvars, wxID_WXDIALOG1PYVARS,
              self.OnPyvarsTogglebutton)

        self.Text3 = wxStaticText(id=wxID_WXDIALOG1TEXT3, label='Suffix:',
              name='Text3', parent=self, pos=wxPoint(4, 70), size=wxSize(33,
              16), style=0)

        self.Varsuffix = wxTextCtrl(id=wxID_WXDIALOG1VARSUFFIX,
              name='Varsuffix', parent=self, pos=wxPoint(60, 68),
              size=wxSize(100, 22), style=wxTE_PROCESS_ENTER, value='')

        self.Dump = wxButton(id=wxID_WXDIALOG1DUMP, label='DUMP', name='Dump',
              parent=self, pos=wxPoint(170, 68), size=wxSize(140, 22), style=0)
        self.Dump.SetFont(wxFont(12, wxSWISS, wxNORMAL, wxBOLD, false, ''))
        self.Dump.SetForegroundColour(wxColour(230, 0, 0))
        EVT_BUTTON(self.Dump, wxID_WXDIALOG1DUMP, self.OnDumpButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.pyvars = 1
        self.suffix = ''
        self.filename = arraytostr(top.runid)+('%06d'%top.it)+self.suffix+'.dump'
        self.Filename.SetValue(self.filename)
        self.attr = 'dump'
        self.Attributes.SetValue('dump')

    def OnBrowserButton(self, event):
        dlg = wxFileDialog(self, "Choose a file", ".", "", 
              "DUMP files (*.dump)|*.dump|PDB files (*.pdb)|*.pdb|ALL files (*.*)|*.*", wxSAVE|wxOVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wxID_OK:
                self.filename = dlg.GetPath()
                self.Filename.SetValue(self.filename)
        finally:
            dlg.Destroy()

    def OnPyvarsTogglebutton(self, event):
        if(self.Pyvars.GetValue()):
            self.pyvars = 1
        else:
            self.pyvars = 0

    def OnDumpButton(self, event):
        self.attr = self.Attributes.GetValue()
        self.filename = self.Filename.GetValue()
        self.varsuffix = self.Varsuffix.GetValue()
        if self.varsuffix is '': self.varsuffix=None
        dump(filename=self.filename,attr=self.attr,pyvars=self.pyvars,varsuffix=self.varsuffix)
        self.Destroy()

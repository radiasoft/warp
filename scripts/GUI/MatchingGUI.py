#Boa:FramePanel:MatchingGUI

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
import matchenv
import newstdout
import warp
import __main__

[wxID_MATCHINGGUI, wxID_MATCHINGGUIENDMATCHOUTPUT, wxID_MATCHINGGUIFINALA, 
 wxID_MATCHINGGUIFINALALABEL, wxID_MATCHINGGUIFINALAP, 
 wxID_MATCHINGGUIFINALAPLABEL, wxID_MATCHINGGUIFINALB, 
 wxID_MATCHINGGUIFINALBLABEL, wxID_MATCHINGGUIFINALBP, 
 wxID_MATCHINGGUIFINALBPLABEL, wxID_MATCHINGGUIMATCHEND, 
 wxID_MATCHINGGUIMATCHINGTYPES, wxID_MATCHINGGUIPANEL1, 
 wxID_MATCHINGGUIPANEL2, wxID_MATCHINGGUIPERIODICLABEL, 
 wxID_MATCHINGGUIPLOTENDMATCH, wxID_MATCHINGGUISETQUAD0, 
 wxID_MATCHINGGUISETQUAD1, wxID_MATCHINGGUISETQUAD2, wxID_MATCHINGGUISETQUAD3, 
 wxID_MATCHINGGUIVARYQUADS, 
] = map(lambda _init_ctrls: wxNewId(), range(21))

class MatchingGUI(wxPanel):
    def _init_coll_MatchingTypes_Pages(self, parent):
        # generated method, don't edit

        parent.AddPage(imageId=-1, page=self.panel1, select=False,
              text='Periodic')
        parent.AddPage(imageId=-1, page=self.panel2, select=True,
              text='End Conditions')

    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_MATCHINGGUI, name='MatchingGUI',
              parent=prnt, pos=wxPoint(380, 304), size=wxSize(548, 315),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(548, 315))
        EVT_PAINT(self, self.OnMatchingguiPaint)

        self.MatchingTypes = wxNotebook(id=wxID_MATCHINGGUIMATCHINGTYPES,
              name='MatchingTypes', parent=self, pos=wxPoint(0, 0),
              size=wxSize(548, 315), style=0)
        self.MatchingTypes.SetToolTipString('Envelope matching')

        self.panel1 = wxPanel(id=wxID_MATCHINGGUIPANEL1, name='panel1',
              parent=self.MatchingTypes, pos=wxPoint(0, 0), size=wxSize(544,
              281), style=wxTAB_TRAVERSAL)

        self.PeriodicLabel = wxStaticText(id=wxID_MATCHINGGUIPERIODICLABEL,
              label='Matching to a periodic lattice', name='PeriodicLabel',
              parent=self.panel1, pos=wxPoint(4, 4), size=wxSize(185, 18),
              style=0)
        self.PeriodicLabel.SetToolTipString('')
        self.PeriodicLabel.SetFont(wxFont(14, wxSWISS, wxNORMAL, wxNORMAL,
              False, ''))

        self.panel2 = wxPanel(id=wxID_MATCHINGGUIPANEL2, name='panel2',
              parent=self.MatchingTypes, pos=wxPoint(0, 0), size=wxSize(544,
              281), style=wxTAB_TRAVERSAL)

        self.FinalaLabel = wxStaticText(id=wxID_MATCHINGGUIFINALALABEL,
              label='Final a', name='FinalaLabel', parent=self.panel2,
              pos=wxPoint(16, 32), size=wxSize(39, 16), style=0)
        self.FinalaLabel.SetToolTipString('')

        self.FinalbLabel = wxStaticText(id=wxID_MATCHINGGUIFINALBLABEL,
              label='Final b', name='FinalbLabel', parent=self.panel2,
              pos=wxPoint(16, 62), size=wxSize(39, 16), style=0)
        self.FinalbLabel.SetToolTipString('')

        self.Finala = wxTextCtrl(id=wxID_MATCHINGGUIFINALA, name='Finala',
              parent=self.panel2, pos=wxPoint(70, 24), size=wxSize(80, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='afinal')
        self.Finala.SetToolTipString('Final value of a to match to')
        EVT_TEXT_ENTER(self.Finala, wxID_MATCHINGGUIFINALA,
              self.OnFinalaTextEnter)

        self.Finalb = wxTextCtrl(id=wxID_MATCHINGGUIFINALB, name='Finalb',
              parent=self.panel2, pos=wxPoint(70, 54), size=wxSize(80, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='bfinal')
        self.Finalb.SetToolTipString('Final value of b to match to')
        EVT_TEXT_ENTER(self.Finalb, wxID_MATCHINGGUIFINALB,
              self.OnFinalbTextEnter)

        self.FinalapLabel = wxStaticText(id=wxID_MATCHINGGUIFINALAPLABEL,
              label="Final a'", name='FinalapLabel', parent=self.panel2,
              pos=wxPoint(16, 92), size=wxSize(42, 16), style=0)
        self.FinalapLabel.SetToolTipString('')

        self.FinalbpLabel = wxStaticText(id=wxID_MATCHINGGUIFINALBPLABEL,
              label="Final b'", name='FinalbpLabel', parent=self.panel2,
              pos=wxPoint(16, 122), size=wxSize(42, 16), style=0)
        self.FinalbpLabel.SetToolTipString('')

        self.Finalap = wxTextCtrl(id=wxID_MATCHINGGUIFINALAP, name='Finalap',
              parent=self.panel2, pos=wxPoint(70, 84), size=wxSize(80, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='apfinal')
        self.Finalap.SetToolTipString("Final value of a' to match to")
        EVT_TEXT_ENTER(self.Finalap, wxID_MATCHINGGUIFINALAP,
              self.OnFinalapTextEnter)

        self.Finalbp = wxTextCtrl(id=wxID_MATCHINGGUIFINALBP, name='Finalbp',
              parent=self.panel2, pos=wxPoint(70, 114), size=wxSize(80, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='bpfinal')
        self.Finalbp.SetToolTipString("Final value of b' to match to")
        EVT_TEXT_ENTER(self.Finalbp, wxID_MATCHINGGUIFINALBP,
              self.OnFinalbpTextEnter)

        self.VaryQuads = wxStaticText(id=wxID_MATCHINGGUIVARYQUADS,
              label='Quads to vary', name='VaryQuads', parent=self.panel2,
              pos=wxPoint(12, 150), size=wxSize(80, 16), style=0)

        self.SetQuad0 = wxTextCtrl(id=wxID_MATCHINGGUISETQUAD0, name='SetQuad0',
              parent=self.panel2, pos=wxPoint(8, 170), size=wxSize(30, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.SetQuad0, wxID_MATCHINGGUISETQUAD0,
              self.OnSetquad0TextEnter)

        self.SetQuad1 = wxTextCtrl(id=wxID_MATCHINGGUISETQUAD1, name='SetQuad1',
              parent=self.panel2, pos=wxPoint(42, 170), size=wxSize(30, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.SetQuad1, wxID_MATCHINGGUISETQUAD1,
              self.OnSetquad1TextEnter)

        self.SetQuad2 = wxTextCtrl(id=wxID_MATCHINGGUISETQUAD2, name='SetQuad2',
              parent=self.panel2, pos=wxPoint(76, 170), size=wxSize(30, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='')
        self.SetQuad2.SetToolTipString('')
        EVT_TEXT_ENTER(self.SetQuad2, wxID_MATCHINGGUISETQUAD2,
              self.OnSetquad2TextEnter)

        self.SetQuad3 = wxTextCtrl(id=wxID_MATCHINGGUISETQUAD3, name='SetQuad3',
              parent=self.panel2, pos=wxPoint(110, 170), size=wxSize(30, 22),
              style=wxTAB_TRAVERSAL | wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER,
              value='')
        self.SetQuad3.SetToolTipString('')
        EVT_TEXT_ENTER(self.SetQuad3, wxID_MATCHINGGUISETQUAD3,
              self.OnSetquad3TextEnter)

        self.MatchEnd = wxButton(id=wxID_MATCHINGGUIMATCHEND, label='Match',
              name='MatchEnd', parent=self.panel2, pos=wxPoint(12, 210),
              size=wxSize(80, 22), style=wxTAB_TRAVERSAL)
        self.MatchEnd.SetToolTipString('')
        EVT_BUTTON(self.MatchEnd, wxID_MATCHINGGUIMATCHEND,
              self.OnMatchendButton)

        self.PlotEndMatch = wxCheckBox(id=wxID_MATCHINGGUIPLOTENDMATCH,
              label='Plot envelope after match', name='PlotEndMatch',
              parent=self.panel2, pos=wxPoint(12, 240), size=wxSize(168, 24),
              style=wxTAB_TRAVERSAL)
        self.PlotEndMatch.SetValue(True)
        EVT_CHECKBOX(self.PlotEndMatch, wxID_MATCHINGGUIPLOTENDMATCH,
              self.OnPlotendmatchCheckbox)

        self.EndMatchOutput = wxTextCtrl(id=wxID_MATCHINGGUIENDMATCHOUTPUT,
              name='EndMatchOutput', parent=self.panel2, pos=wxPoint(180, 24),
              size=wxSize(360, 250), style=wxTE_MULTILINE, value='')
        self.EndMatchOutput.SetToolTipString('')

        self._init_coll_MatchingTypes_Pages(self.MatchingTypes)

    def __init__(self, parent, id=0, pos=0, size=0, style=0, name=0):
        self._init_ctrls(parent)
        parent.AddPage(imageId=-1, page=self, select=True, text='Matching')
        self.plotafterendmatch = self.PlotEndMatch.GetValue()
        self.endmatchquads = [None,None,None,None]

    def OnFinalaTextEnter(self, event,defval=''):
        try:
            self.afinal = eval(self.Finala.GetValue(),__main__.__dict__)
            self.Finala.SetValue(str(self.afinal))
        except:
            self.Finala.SetValue(defval)

    def OnFinalbTextEnter(self, event,defval=''):
        try:
            self.bfinal = eval(self.Finalb.GetValue(),__main__.__dict__)
            self.Finalb.SetValue(str(self.bfinal))
        except:
            self.Finalb.SetValue(defval)

    def OnFinalapTextEnter(self, event,defval=''):
        try:
            self.apfinal = eval(self.Finalap.GetValue(),__main__.__dict__)
            self.Finalap.SetValue(str(self.apfinal))
        except:
            self.Finalap.SetValue(defval)

    def OnFinalbpTextEnter(self, event,defval=''):
        try:
            self.bpfinal = eval(self.Finalbp.GetValue(),__main__.__dict__)
            self.Finalbp.SetValue(str(self.bpfinal))
        except:
            self.Finalbp.SetValue(defval)

    def OnMatchendButton(self, event):
        savestdout = sys.stdout
        sys.stdout = newstdout.newstdout(self.EndMatchOutput)
        matchenv.matchenv(self.endmatchquads,self.afinal,self.bfinal,self.apfinal,self.bpfinal,usequad=1)
        sys.stdout = savestdout
        if self.plotafterendmatch:
            fma()
            penv()

    def OnPlotendmatchCheckbox(self, event):
        self.plotafterendmatch = self.PlotEndMatch.GetValue()

    def OnSetquad0TextEnter(self, event):
        self.endmatchquads[0] = eval(self.SetQuad0.GetValue())
        if self.endmatchquads[1] is None:
            self.endmatchquads[1] = self.endmatchquads[0] + 1
            self.endmatchquads[2] = self.endmatchquads[1] + 1
            self.endmatchquads[3] = self.endmatchquads[2] + 1
            self.SetQuad1.SetValue(str(self.endmatchquads[1]))
            self.SetQuad2.SetValue(str(self.endmatchquads[2]))
            self.SetQuad3.SetValue(str(self.endmatchquads[3]))

    def OnSetquad1TextEnter(self, event):
        self.endmatchquads[1] = eval(self.SetQuad1.GetValue())

    def OnSetquad2TextEnter(self, event):
        self.endmatchquads[2] = eval(self.SetQuad2.GetValue())

    def OnSetquad3TextEnter(self, event):
        self.endmatchquads[3] = eval(self.SetQuad3.GetValue())

    def OnMatchingguiPaint(self, event):
        warp.package('env')
        self.OnFinalaTextEnter(None,'afinal')
        self.OnFinalbTextEnter(None,'bfinal')
        self.OnFinalapTextEnter(None,'apfinal')
        self.OnFinalbpTextEnter(None,'bpfinal')

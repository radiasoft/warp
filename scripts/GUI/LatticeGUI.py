#Boa:Dialog:LatticeGUI

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from warp import *
import sortlattice

def create(parent):
    return LatticeGUI(parent)

[wxID_LATTICEGUI, wxID_LATTICEGUIDOSTEP, wxID_LATTICEGUIDRFT, wxID_LATTICEGUIDRFTAPLABEL, wxID_LATTICEGUIDRFTAPUNITS, wxID_LATTICEGUIDRFTLABEL, wxID_LATTICEGUIDRFTZELABEL, wxID_LATTICEGUIDRFTZEUNITS, wxID_LATTICEGUIDRFTZSLABEL, wxID_LATTICEGUIDRFTZSUNITS, wxID_LATTICEGUIELEMENTNUM, wxID_LATTICEGUIELEMENTSPIN, wxID_LATTICEGUIELEMNUMLABEL1, wxID_LATTICEGUIELEMNUMLABEL2, wxID_LATTICEGUIGETDRFTAP, wxID_LATTICEGUIGETDRFTZE, wxID_LATTICEGUIGETDRFTZS, wxID_LATTICEGUIGETQUADDB, wxID_LATTICEGUIGETQUADDE, wxID_LATTICEGUIGETQUADZE, wxID_LATTICEGUIGETQUADZS, wxID_LATTICEGUIQUAD, wxID_LATTICEGUIQUADDBLABEL, wxID_LATTICEGUIQUADDBUNITS, wxID_LATTICEGUIQUADDELABEL, wxID_LATTICEGUIQUADDEUNITS, wxID_LATTICEGUIQUADLABEL, wxID_LATTICEGUIQUADZELABEL, wxID_LATTICEGUIQUADZEUNITS, wxID_LATTICEGUIQUADZSLABEL, wxID_LATTICEGUIQUADZSUNITS, wxID_LATTICEGUISETELEMENTNUM] = map(lambda _init_ctrls: wxNewId(), range(32))

class LatticeGUI(wxDialog):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wxDialog.__init__(self, id = wxID_LATTICEGUI, name = 'LatticeGUI', parent = prnt, pos = wxPoint(381, 289), size = wxSize(389, 277), style = wxDEFAULT_DIALOG_STYLE, title = 'Lattice Editor')
        self._init_utils()
        self.SetClientSize(wxSize(389, 277))
        self.SetToolTipString('Lattice editor')

        self.ElementNum = wxPanel(id = wxID_LATTICEGUIELEMENTNUM, name = 'ElementNum', parent = self, pos = wxPoint(8, 8), size = wxSize(136, 264), style = wxRAISED_BORDER | wxTAB_TRAVERSAL)

        self.SetElementNum = wxTextCtrl(id = wxID_LATTICEGUISETELEMENTNUM, name = 'SetElementNum', parent = self.ElementNum, pos = wxPoint(8, 48), size = wxSize(40, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '0')
        self.SetElementNum.SetToolTipString('Sets element number')
        EVT_TEXT_ENTER(self.SetElementNum, wxID_LATTICEGUISETELEMENTNUM, self.OnSetelementnumTextEnter)

        self.ElemNumLabel1 = wxStaticText(id = wxID_LATTICEGUIELEMNUMLABEL1, label = 'Element', name = 'ElemNumLabel1', parent = self.ElementNum, pos = wxPoint(8, 8), size = wxSize(44, 16), style = 0)

        self.ElemNumLabel2 = wxStaticText(id = wxID_LATTICEGUIELEMNUMLABEL2, label = 'Number', name = 'ElemNumLabel2', parent = self.ElementNum, pos = wxPoint(8, 24), size = wxSize(43, 16), style = 0)

        self.ElementSpin = wxSpinButton(id = wxID_LATTICEGUIELEMENTSPIN, name = 'ElementSpin', parent = self.ElementNum, pos = wxPoint(56, 40), size = wxSize(15, 34), style = wxDOUBLE_BORDER | wxSP_HORIZONTAL | wxSIMPLE_BORDER)
        self.ElementSpin.SetToolTipString('Change element number')
        EVT_SPIN_DOWN(self.ElementSpin, wxID_LATTICEGUIELEMENTSPIN, self.OnElementspinSpinDown)
        EVT_SPIN_UP(self.ElementSpin, wxID_LATTICEGUIELEMENTSPIN, self.OnElementspinSpinUp)

        self.Quad = wxPanel(id = wxID_LATTICEGUIQUAD, name = 'Quad', parent = self, pos = wxPoint(144, 8), size = wxSize(240, 264), style = wxSUNKEN_BORDER | wxTAB_TRAVERSAL)
        self.Quad.Show(false)

        self.quadLabel = wxStaticText(id = wxID_LATTICEGUIQUADLABEL, label = 'Hard edged quadrupole', name = 'quadLabel', parent = self.Quad, pos = wxPoint(42, 0), size = wxSize(151, 18), style = wxSIMPLE_BORDER | wxALIGN_CENTRE)
        self.quadLabel.Center(wxHORIZONTAL)
        self.quadLabel.SetFont(wxFont(14, wxSWISS, wxNORMAL, wxNORMAL, false, ''))

        self.quadzsLabel = wxStaticText(id = wxID_LATTICEGUIQUADZSLABEL, label = 'Z start', name = 'quadzsLabel', parent = self.Quad, pos = wxPoint(8, 27), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.quadzeLabel = wxStaticText(id = wxID_LATTICEGUIQUADZELABEL, label = 'Z end', name = 'quadzeLabel', parent = self.Quad, pos = wxPoint(8, 51), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.quaddeLabel = wxStaticText(id = wxID_LATTICEGUIQUADDELABEL, label = 'E gradient', name = 'quaddeLabel', parent = self.Quad, pos = wxPoint(8, 75), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.quaddbLabel = wxStaticText(id = wxID_LATTICEGUIQUADDBLABEL, label = 'B gradient', name = 'quaddbLabel', parent = self.Quad, pos = wxPoint(8, 99), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.quadzsUnits = wxStaticText(id = wxID_LATTICEGUIQUADZSUNITS, label = 'meters', name = 'quadzsUnits', parent = self.Quad, pos = wxPoint(152, 27), size = wxSize(36, 16), style = 0)

        self.quadzeUnits = wxStaticText(id = wxID_LATTICEGUIQUADZEUNITS, label = 'meters', name = 'quadzeUnits', parent = self.Quad, pos = wxPoint(152, 51), size = wxSize(36, 16), style = 0)

        self.quaddeUnits = wxStaticText(id = wxID_LATTICEGUIQUADDEUNITS, label = 'V/m^2', name = 'quaddeUnits', parent = self.Quad, pos = wxPoint(152, 75), size = wxSize(35, 16), style = 0)

        self.quaddbUnits = wxStaticText(id = wxID_LATTICEGUIQUADDBUNITS, label = 'B/m^2', name = 'quaddbUnits', parent = self.Quad, pos = wxPoint(152, 99), size = wxSize(34, 16), style = 0)

        self.Getquadzs = wxTextCtrl(id = wxID_LATTICEGUIGETQUADZS, name = 'Getquadzs', parent = self.Quad, pos = wxPoint(72, 24), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getquadzs, wxID_LATTICEGUIGETQUADZS, self.OnGetelemzsTextEnter)

        self.Getquadze = wxTextCtrl(id = wxID_LATTICEGUIGETQUADZE, name = 'Getquadze', parent = self.Quad, pos = wxPoint(72, 48), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getquadze, wxID_LATTICEGUIGETQUADZE, self.OnGetelemzeTextEnter)

        self.Getquadde = wxTextCtrl(id = wxID_LATTICEGUIGETQUADDE, name = 'Getquadde', parent = self.Quad, pos = wxPoint(72, 72), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getquadde, wxID_LATTICEGUIGETQUADDE, self.OnGetquaddeTextEnter)

        self.Getquaddb = wxTextCtrl(id = wxID_LATTICEGUIGETQUADDB, name = 'Getquaddb', parent = self.Quad, pos = wxPoint(72, 96), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getquaddb, wxID_LATTICEGUIGETQUADDB, self.OnGetquaddbTextEnter)

        self.Drft = wxPanel(id = wxID_LATTICEGUIDRFT, name = 'Drft', parent = self, pos = wxPoint(144, 8), size = wxSize(240, 264), style = wxSUNKEN_BORDER | wxTAB_TRAVERSAL)
        self.Drft.SetToolTipString('Drift elements')

        self.drftLabel = wxStaticText(id = wxID_LATTICEGUIDRFTLABEL, label = 'Drift', name = 'drftLabel', parent = self.Drft, pos = wxPoint(105, 0), size = wxSize(25, 18), style = wxALIGN_CENTRE)
        self.drftLabel.Center(wxHORIZONTAL)
        self.drftLabel.SetFont(wxFont(14, wxSWISS, wxNORMAL, wxNORMAL, false, ''))

        self.drftzsLabel = wxStaticText(id = wxID_LATTICEGUIDRFTZSLABEL, label = 'Z start', name = 'drftzsLabel', parent = self.Drft, pos = wxPoint(8, 27), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.drftzeLabel = wxStaticText(id = wxID_LATTICEGUIDRFTZELABEL, label = 'Z end', name = 'drftzeLabel', parent = self.Drft, pos = wxPoint(8, 51), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.drftapLabel = wxStaticText(id = wxID_LATTICEGUIDRFTAPLABEL, label = 'Aperture', name = 'drftapLabel', parent = self.Drft, pos = wxPoint(8, 75), size = wxSize(60, 16), style = wxALIGN_RIGHT)

        self.drftzsUnits = wxStaticText(id = wxID_LATTICEGUIDRFTZSUNITS, label = 'meters', name = 'drftzsUnits', parent = self.Drft, pos = wxPoint(152, 27), size = wxSize(36, 16), style = 0)

        self.drftzeUnits = wxStaticText(id = wxID_LATTICEGUIDRFTZEUNITS, label = 'meters', name = 'drftzeUnits', parent = self.Drft, pos = wxPoint(152, 51), size = wxSize(36, 16), style = 0)

        self.drftapUnits = wxStaticText(id = wxID_LATTICEGUIDRFTAPUNITS, label = 'meters', name = 'drftapUnits', parent = self.Drft, pos = wxPoint(152, 75), size = wxSize(36, 16), style = 0)

        self.Getdrftzs = wxTextCtrl(id = wxID_LATTICEGUIGETDRFTZS, name = 'Getdrftzs', parent = self.Drft, pos = wxPoint(72, 24), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getdrftzs, wxID_LATTICEGUIGETDRFTZS, self.OnGetelemzsTextEnter)

        self.Getdrftze = wxTextCtrl(id = wxID_LATTICEGUIGETDRFTZE, name = 'Getdrftze', parent = self.Drft, pos = wxPoint(72, 48), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getdrftze, wxID_LATTICEGUIGETDRFTZE, self.OnGetelemzeTextEnter)

        self.Getdrftap = wxTextCtrl(id = wxID_LATTICEGUIGETDRFTAP, name = 'Getdrftap', parent = self.Drft, pos = wxPoint(72, 72), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.Getdrftap, wxID_LATTICEGUIGETDRFTAP, self.OnGetelemapTextEnter)

        self.DoStep = wxCheckBox(id = wxID_LATTICEGUIDOSTEP, label = 'Step on change', name = 'DoStep', parent = self.ElementNum, pos = wxPoint(8, 96), size = wxSize(120, 24), style = 0)
        self.DoStep.SetValue(false)
        self.DoStep.SetHelpText('When checked, execute step command on a change.')
        self.DoStep.SetToolTipString('Turns on code calculation on change')
        EVT_CHECKBOX(self.DoStep, wxID_LATTICEGUIDOSTEP, self.OnDostepCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.elemnum = 0
        self.sortedelems = sortlattice.sortlattice()
        self.numelements = len(self.sortedelems)
        self.shownelem = None
        self.UpdatePanel()
        self.steponchange = self.DoStep.GetValue()

    def Step(self):
        if self.steponchange:
            step()

    def OnSetelementnumTextEnter(self, event):
        try:
            self.elemnum = eval(self.SetElementNum.GetValue())
            self.UpdatePanel()
        except:
            self.SetElementNum.SetValue("%d"%self.elemnum)

    def OnElementspinSpinDown(self, event):
        self.elemnum = (self.elemnum - 1 + self.numelements) % self.numelements
        self.SetElementNum.SetValue("%d"%self.elemnum)
        self.UpdatePanel()

    def OnElementspinSpinUp(self, event):
        self.elemnum = (self.elemnum + 1 + self.numelements) % self.numelements
        self.SetElementNum.SetValue("%d"%self.elemnum)
        self.UpdatePanel()
        
    def OnDostepCheckbox(self, event):
        self.steponchange = self.DoStep.GetValue()

    def GetElemValue(self,val,valtype,valformat):
        el = self.sortedelems[self.elemnum]
        ctl = self.__dict__["Get"+el.type+val]
        try:
            value = valtype(ctl.GetValue())
            el.__setattr__(val,value)
            self.Step()
        except ValueError:
            value = el.__getattr__(val)
            ctl.SetValue(valformat%value)
            
    def GetElemFloatValue(self,val):
        self.GetElemValue(val,float,"%f")

    def GetElemIntValue(self,val):
        self.GetElemValue(val,int,"%d")
        
    def OnGetelemzsTextEnter(self, event):
        self.GetElemFloatValue('zs')

    def OnGetelemzeTextEnter(self, event):
        self.GetElemFloatValue('ze')

    def OnGetelemapTextEnter(self, event):
        self.GetElemFloatValue('ap')

    def OnGetquaddeTextEnter(self, event):
        self.GetElemValue('de',float,"%e")

    def OnGetquaddbTextEnter(self, event):
        self.GetElemFloatValue('db')

    def UpdatePanel(self):
        if self.shownelem is not None: self.shownelem.Show(0)
        el = self.sortedelems[self.elemnum]
        if el.type == 'quad':
            self.Getquadzs.SetValue("%f"%el.zs)
            self.Getquadze.SetValue("%f"%el.ze)
            self.Getquadde.SetValue("%e"%el.de)
            self.Getquaddb.SetValue("%f"%el.db)
            self.Quad.Show(1)
            self.shownelem = self.Quad
        elif el.type == 'drft':
            self.Getdrftzs.SetValue("%f"%el.zs)
            self.Getdrftze.SetValue("%f"%el.ze)
            self.Getdrftap.SetValue("%f"%el.ap)
            self.Drft.Show(1)
            self.shownelem = self.Drft



#Boa:Frame:EnvelopeGUI

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from warp import *

def create(parent):
    return EnvelopeGUI(parent)

[wxID_ENVELOPEGUI, wxID_ENVELOPEGUIAIONLABEL, wxID_ENVELOPEGUIANGLELABEL, wxID_ENVELOPEGUIANGLEUNITS, wxID_ENVELOPEGUIBUTTON1, wxID_ENVELOPEGUIBUTTON2, wxID_ENVELOPEGUIBUTTON3, wxID_ENVELOPEGUIBUTTON4, wxID_ENVELOPEGUIBUTTON5, wxID_ENVELOPEGUICURRENTLABEL, wxID_ENVELOPEGUICURRENTUNITS, wxID_ENVELOPEGUIEMITNLABEL, wxID_ENVELOPEGUIEMITNUNITS, wxID_ENVELOPEGUIENERGYLABEL, wxID_ENVELOPEGUIENERGYUNITS, wxID_ENVELOPEGUIFMA, wxID_ENVELOPEGUIHCP, wxID_ENVELOPEGUIINITIALIZE, wxID_ENVELOPEGUIINPUTPANEL, wxID_ENVELOPEGUIPLOTSPANEL, wxID_ENVELOPEGUIPRINTENVOUTPUT, wxID_ENVELOPEGUIREDRAW, wxID_ENVELOPEGUIRUNONCHANGE, wxID_ENVELOPEGUISETAION, wxID_ENVELOPEGUISETCURRENT, wxID_ENVELOPEGUISETENERGY, wxID_ENVELOPEGUISETFMA, wxID_ENVELOPEGUISETXEMITN, wxID_ENVELOPEGUISETXPSTART, wxID_ENVELOPEGUISETXSTART, wxID_ENVELOPEGUISETYEMITN, wxID_ENVELOPEGUISETYPSTART, wxID_ENVELOPEGUISETYSTART, wxID_ENVELOPEGUISETZION, wxID_ENVELOPEGUISIZELABEL, wxID_ENVELOPEGUISTATICTEXT1, wxID_ENVELOPEGUIXLABEL, wxID_ENVELOPEGUIXPANDYP, wxID_ENVELOPEGUIXSTARTUNITS, wxID_ENVELOPEGUIYLABEL, wxID_ENVELOPEGUIYP, wxID_ENVELOPEGUIZIONLABEL] = map(lambda _init_ctrls: wxNewId(), range(42))

class EnvelopeGUI(wxFrame):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wxFrame.__init__(self, id = wxID_ENVELOPEGUI, name = 'EnvelopeGUI', parent = prnt, pos = wxPoint(435, 313), size = wxSize(597, 340), style = wxDEFAULT_FRAME_STYLE, title = 'wxFrame1')
        self._init_utils()
        self.SetClientSize(wxSize(597, 340))
        EVT_CLOSE(self, self.OnEnvelopeguiClose)

        self.PlotsPanel = wxPanel(id = wxID_ENVELOPEGUIPLOTSPANEL, name = 'PlotsPanel', parent = self, pos = wxPoint(64, 0), size = wxSize(104, 216), style = wxRAISED_BORDER | wxTAB_TRAVERSAL)

        self.button1 = wxButton(id = wxID_ENVELOPEGUIBUTTON1, label = 'Run', name = 'button1', parent = self, pos = wxPoint(8, 8), size = wxSize(48, 56), style = 0)
        EVT_BUTTON(self.button1, wxID_ENVELOPEGUIBUTTON1, self.OnStepButton)

        self.Initialize = wxButton(id = wxID_ENVELOPEGUIINITIALIZE, label = 'Init', name = 'Initialize', parent = self, pos = wxPoint(8, 64), size = wxSize(48, 56), style = 0)
        EVT_BUTTON(self.Initialize, wxID_ENVELOPEGUIINITIALIZE, self.OnInitializeButton)

        self.button2 = wxButton(id = wxID_ENVELOPEGUIBUTTON2, label = 'x and y', name = 'button2', parent = self.PlotsPanel, pos = wxPoint(8, 56), size = wxSize(80, 22), style = 0)
        EVT_BUTTON(self.button2, wxID_ENVELOPEGUIBUTTON2, self.OnPlotXYButton)

        self.staticText1 = wxStaticText(id = wxID_ENVELOPEGUISTATICTEXT1, label = 'Plots', name = 'staticText1', parent = self.PlotsPanel, pos = wxPoint(22, 0), size = wxSize(56, 16), style = wxALIGN_CENTRE)
        self.staticText1.Center(wxHORIZONTAL)

        self.button3 = wxButton(id = wxID_ENVELOPEGUIBUTTON3, label = 'x', name = 'button3', parent = self.PlotsPanel, pos = wxPoint(8, 32), size = wxSize(40, 22), style = 0)
        EVT_BUTTON(self.button3, wxID_ENVELOPEGUIBUTTON3, self.OnPlotXButton)

        self.button4 = wxButton(id = wxID_ENVELOPEGUIBUTTON4, label = 'y', name = 'button4', parent = self.PlotsPanel, pos = wxPoint(48, 32), size = wxSize(40, 22), style = 0)
        EVT_BUTTON(self.button4, wxID_ENVELOPEGUIBUTTON4, self.OnPlotYButton)

        self.button5 = wxButton(id = wxID_ENVELOPEGUIBUTTON5, label = "x'", name = 'button5', parent = self.PlotsPanel, pos = wxPoint(8, 88), size = wxSize(40, 22), style = 0)
        EVT_BUTTON(self.button5, wxID_ENVELOPEGUIBUTTON5, self.OnPlotXPButton)

        self.yp = wxButton(id = wxID_ENVELOPEGUIYP, label = "y'", name = 'yp', parent = self.PlotsPanel, pos = wxPoint(48, 88), size = wxSize(40, 22), style = 0)
        EVT_BUTTON(self.yp, wxID_ENVELOPEGUIYP, self.OnPlotYPButton)

        self.xpandyp = wxButton(id = wxID_ENVELOPEGUIXPANDYP, label = "x' and y'", name = 'xpandyp', parent = self.PlotsPanel, pos = wxPoint(8, 112), size = wxSize(80, 22), style = 0)
        EVT_BUTTON(self.xpandyp, wxID_ENVELOPEGUIXPANDYP, self.OnPlotXPYPButton)

        self.fma = wxButton(id = wxID_ENVELOPEGUIFMA, label = 'fma', name = 'fma', parent = self, pos = wxPoint(8, 160), size = wxSize(48, 32), style = 0)
        EVT_BUTTON(self.fma, wxID_ENVELOPEGUIFMA, self.OnFmaButton)

        self.hcp = wxButton(id = wxID_ENVELOPEGUIHCP, label = 'hcp', name = 'hcp', parent = self, pos = wxPoint(8, 192), size = wxSize(48, 32), style = 0)
        EVT_BUTTON(self.hcp, wxID_ENVELOPEGUIHCP, self.OnHcpButton)

        self.RunOnChange = wxCheckBox(id = wxID_ENVELOPEGUIRUNONCHANGE, label = 'Run on change', name = 'RunOnChange', parent = self, pos = wxPoint(216, 256), size = wxSize(112, 24), style = 0)
        self.RunOnChange.SetValue(true)
        self.RunOnChange.SetHelpText('Sets whether to run envelope calculation on any change of parameters.')
        self.RunOnChange.SetToolTipString('')
        EVT_CHECKBOX(self.RunOnChange, wxID_ENVELOPEGUIRUNONCHANGE, self.OnRunonchangeCheckbox)

        self.SetFma = wxCheckBox(id = wxID_ENVELOPEGUISETFMA, label = 'Frame advance before each plot', name = 'SetFma', parent = self, pos = wxPoint(216, 272), size = wxSize(208, 24), style = 0)
        self.SetFma.SetValue(true)
        EVT_CHECKBOX(self.SetFma, wxID_ENVELOPEGUISETFMA, self.OnSetfmaCheckbox)

        self.InputPanel = wxPanel(id = wxID_ENVELOPEGUIINPUTPANEL, name = 'InputPanel', parent = self, pos = wxPoint(168, 0), size = wxSize(312, 216), style = wxRAISED_BORDER | wxTAB_TRAVERSAL)

        self.SetXstart = wxTextCtrl(id = wxID_ENVELOPEGUISETXSTART, name = 'SetXstart', parent = self.InputPanel, pos = wxPoint(56, 80), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetXstart.SetToolTipString('Set X')
        EVT_TEXT_ENTER(self.SetXstart, wxID_ENVELOPEGUISETXSTART, self.OnSetXstartTextEnter)

        self.CurrentLabel = wxStaticText(id = wxID_ENVELOPEGUICURRENTLABEL, label = 'Current', name = 'CurrentLabel', parent = self.InputPanel, pos = wxPoint(8, 11), size = wxSize(41, 16), style = wxALIGN_RIGHT)

        self.SetCurrent = wxTextCtrl(id = wxID_ENVELOPEGUISETCURRENT, name = 'SetCurrent', parent = self.InputPanel, pos = wxPoint(56, 8), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.SetCurrent, wxID_ENVELOPEGUISETCURRENT, self.OnSetcurrentTextEnter)

        self.CurrentUnits = wxStaticText(id = wxID_ENVELOPEGUICURRENTUNITS, label = 'Amps', name = 'CurrentUnits', parent = self.InputPanel, pos = wxPoint(136, 11), size = wxSize(31, 16), style = 0)

        self.XstartUnits = wxStaticText(id = wxID_ENVELOPEGUIXSTARTUNITS, label = 'mm', name = 'XstartUnits', parent = self.InputPanel, pos = wxPoint(216, 83), size = wxSize(18, 16), style = 0)

        self.SizeLabel = wxStaticText(id = wxID_ENVELOPEGUISIZELABEL, label = 'Size', name = 'SizeLabel', parent = self.InputPanel, pos = wxPoint(24, 83), size = wxSize(24, 16), style = wxALIGN_RIGHT)

        self.SetYstart = wxTextCtrl(id = wxID_ENVELOPEGUISETYSTART, name = 'SetYstart', parent = self.InputPanel, pos = wxPoint(136, 80), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetYstart.SetToolTipString('Set Y')
        EVT_TEXT_ENTER(self.SetYstart, wxID_ENVELOPEGUISETYSTART, self.OnSetystartTextEnter)

        self.EmitnLabel = wxStaticText(id = wxID_ENVELOPEGUIEMITNLABEL, label = 'Emitn', name = 'EmitnLabel', parent = self.InputPanel, pos = wxPoint(16, 131), size = wxSize(30, 16), style = 0)

        self.SetXEmitn = wxTextCtrl(id = wxID_ENVELOPEGUISETXEMITN, name = 'SetXEmitn', parent = self.InputPanel, pos = wxPoint(56, 128), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetXEmitn.SetToolTipString('Set X normalized emittance')
        EVT_TEXT_ENTER(self.SetXEmitn, wxID_ENVELOPEGUISETXEMITN, self.OnSetxemitnTextEnter)

        self.SetYEmitn = wxTextCtrl(id = wxID_ENVELOPEGUISETYEMITN, name = 'SetYEmitn', parent = self.InputPanel, pos = wxPoint(136, 128), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetYEmitn.SetToolTipString('Set Y normalized emittance')
        EVT_TEXT_ENTER(self.SetYEmitn, wxID_ENVELOPEGUISETYEMITN, self.OnSetyemitnTextEnter)

        self.XLabel = wxStaticText(id = wxID_ENVELOPEGUIXLABEL, label = 'X', name = 'XLabel', parent = self.InputPanel, pos = wxPoint(64, 64), size = wxSize(16, 16), style = 0)

        self.YLabel = wxStaticText(id = wxID_ENVELOPEGUIYLABEL, label = 'Y', name = 'YLabel', parent = self.InputPanel, pos = wxPoint(144, 64), size = wxSize(24, 16), style = 0)

        self.EmitnUnits = wxStaticText(id = wxID_ENVELOPEGUIEMITNUNITS, label = 'pi-mm-mrad', name = 'EmitnUnits', parent = self.InputPanel, pos = wxPoint(216, 131), size = wxSize(65, 16), style = 0)

        self.AngleLabel = wxStaticText(id = wxID_ENVELOPEGUIANGLELABEL, label = 'Angle', name = 'AngleLabel', parent = self.InputPanel, pos = wxPoint(16, 107), size = wxSize(33, 16), style = 0)

        self.SetXpstart = wxTextCtrl(id = wxID_ENVELOPEGUISETXPSTART, name = 'SetXpstart', parent = self.InputPanel, pos = wxPoint(56, 104), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetXpstart.SetToolTipString('Set X angle')
        EVT_TEXT_ENTER(self.SetXpstart, wxID_ENVELOPEGUISETXPSTART, self.OnSetxpstartTextEnter)

        self.SetYpstart = wxTextCtrl(id = wxID_ENVELOPEGUISETYPSTART, name = 'SetYpstart', parent = self.InputPanel, pos = wxPoint(136, 104), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetYpstart.SetToolTipString('Set Y angle')
        EVT_TEXT_ENTER(self.SetYpstart, wxID_ENVELOPEGUISETYPSTART, self.OnSetypstartTextEnter)

        self.AngleUnits = wxStaticText(id = wxID_ENVELOPEGUIANGLEUNITS, label = 'mrad', name = 'AngleUnits', parent = self.InputPanel, pos = wxPoint(216, 107), size = wxSize(27, 16), style = 0)

        self.Redraw = wxButton(id = wxID_ENVELOPEGUIREDRAW, label = 'Redraw', name = 'Redraw', parent = self, pos = wxPoint(8, 128), size = wxSize(48, 32), style = 0)
        EVT_BUTTON(self.Redraw, wxID_ENVELOPEGUIREDRAW, self.OnRedrawButton)

        self.EnergyLabel = wxStaticText(id = wxID_ENVELOPEGUIENERGYLABEL, label = 'Energy', name = 'EnergyLabel', parent = self.InputPanel, pos = wxPoint(8, 35), size = wxSize(40, 16), style = 0)

        self.EnergyUnits = wxStaticText(id = wxID_ENVELOPEGUIENERGYUNITS, label = 'MV', name = 'EnergyUnits', parent = self.InputPanel, pos = wxPoint(136, 35), size = wxSize(20, 16), style = 0)

        self.SetEnergy = wxTextCtrl(id = wxID_ENVELOPEGUISETENERGY, name = 'SetEnergy', parent = self.InputPanel, pos = wxPoint(56, 32), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetEnergy.SetToolTipString('Kinetic energy')
        EVT_TEXT_ENTER(self.SetEnergy, wxID_ENVELOPEGUISETENERGY, self.OnSetenergyTextEnter)

        self.AionLabel = wxStaticText(id = wxID_ENVELOPEGUIAIONLABEL, label = 'Aion', name = 'AionLabel', parent = self.InputPanel, pos = wxPoint(184, 11), size = wxSize(26, 16), style = 0)

        self.ZionLabel = wxStaticText(id = wxID_ENVELOPEGUIZIONLABEL, label = 'Zion', name = 'ZionLabel', parent = self.InputPanel, pos = wxPoint(184, 35), size = wxSize(26, 16), style = 0)

        self.SetAion = wxTextCtrl(id = wxID_ENVELOPEGUISETAION, name = 'SetAion', parent = self.InputPanel, pos = wxPoint(216, 8), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        EVT_TEXT_ENTER(self.SetAion, wxID_ENVELOPEGUISETAION, self.OnSetaionTextEnter)

        self.SetZion = wxTextCtrl(id = wxID_ENVELOPEGUISETZION, name = 'SetZion', parent = self.InputPanel, pos = wxPoint(216, 32), size = wxSize(80, 22), style = wxTE_PROCESS_TAB | wxTE_PROCESS_ENTER, value = '')
        self.SetZion.SetToolTipString('Set charge')
        EVT_TEXT_ENTER(self.SetZion, wxID_ENVELOPEGUISETZION, self.OnSetzionTextEnter)

        self.PrintEnvOutput = wxCheckBox(id = wxID_ENVELOPEGUIPRINTENVOUTPUT, label = 'Print envelope output', name = 'PrintEnvOutput', parent = self, pos = wxPoint(216, 288), size = wxSize(200, 24), style = 0)
        self.PrintEnvOutput.SetValue(false)
        EVT_CHECKBOX(self.PrintEnvOutput, wxID_ENVELOPEGUIPRINTENVOUTPUT, self.OnPrintenvoutputCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        package('env')
        self.plotcolor = 'fg'
        self.dorunonchange = self.RunOnChange.GetValue()
        self.dofmabeforeplot = self.SetFma.GetValue()
        self.OnPrintenvoutputCheckbox(None)
        self.SetCurrent.SetValue("%.4f"%top.ibeam)
        self.SetEnergy.SetValue("%.4f"%(top.ekin*1.e-6))
        self.SetAion.SetValue("%.4f"%(top.aion))
        self.SetZion.SetValue("%.1f"%(top.zion))
        self.SetXstart.SetValue("%.4f"%(top.a0*1.e3))
        self.SetYstart.SetValue("%.4f"%(top.b0*1.e3))
        self.SetXpstart.SetValue("%.4f"%(top.ap0*1.e3))
        self.SetYpstart.SetValue("%.4f"%(top.bp0*1.e3))
        self.SetXEmitn.SetValue("%.4f"%(top.emitnx*1.e6))
        if top.emitny != top.emitnx: self.SetYEmitn.SetValue("%.4f"%(top.emitny*1.e6))
        self.xsave = []
        self.ysave = []
        installafterstep(self.MakeEnvPlot)

    def OnEnvelopeguiClose(self, event):
        uninstallafterstep(self.MakeEnvPlot)

    def OnStepButton(self, event):
        self.DoStep(override=1)

    def OnInitializeButton(self, event):
        envgen()
                
    def OnFmaButton(self, event):
        fma()

    def OnHcpButton(self, event):
        hcp()

    def OnRedrawButton(self, event):
        redraw()

    def OnPlotXYButton(self, event):
        self.MakeEnvPlot(["aenv","benv"],["zenv","zenv"])

    def OnPlotXButton(self, event):
        self.MakeEnvPlot(["aenv"],["zenv"])

    def OnPlotYButton(self, event):
        self.MakeEnvPlot(["benv"],["zenv"])

    def OnPlotXPButton(self, event):
        self.MakeEnvPlot(["apenv"],["zenv"])

    def OnPlotYPButton(self, event):
        self.MakeEnvPlot(["bpenv"],["zenv"])

    def OnPlotXPYPButton(self, event):
        self.MakeEnvPlot(["apenv","bpenv"],["zenv","zenv"])

    def OnRunonchangeCheckbox(self, event):
        self.dorunonchange = self.RunOnChange.GetValue()

    def OnPrintenvoutputCheckbox(self, event):
        env.lenvout = self.PrintEnvOutput.GetValue()

    def OnSetfmaCheckbox(self, event):
        self.dofmabeforeplot = self.SetFma.GetValue()

    def OnSetcurrentTextEnter(self, event):
        top.ibeam = eval(self.SetCurrent.GetValue())
        self.DoStep()

    def OnSetXstartTextEnter(self, event):
        top.a0 = 1.e-3*eval(self.SetXstart.GetValue())
        self.DoStep()

    def OnSetystartTextEnter(self, event):
        top.b0 = 1.e-3*eval(self.SetYstart.GetValue())
        self.DoStep()
  
    def OnSetxemitnTextEnter(self, event):
        top.emitnx = 1.e-6*eval(self.SetXEmitn.GetValue())
        if self.SetYEmitn.GetValue() == '': top.emitny = top.emitnx
        self.DoStep()

    def OnSetyemitnTextEnter(self, event):
        top.emitny = 1.e-6*eval(self.SetYEmitn.GetValue())
        self.DoStep()

    def OnSetxpstartTextEnter(self, event):
        top.ap0 = 1.e-3*eval(self.SetXpstart.GetValue())
        self.DoStep()

    def OnSetypstartTextEnter(self, event):
        top.bp0 = 1.e-3*eval(self.SetYpstart.GetValue())
        self.DoStep()

    def OnSetenergyTextEnter(self, event):
        top.ekin = 1.e6*eval(self.SetEnergy.GetValue())
        top.vbeam = 0.
        top.vbeam_s = 0.
        derivqty()
        self.DoStep()

    def OnSetaionTextEnter(self, event):
        top.aion = eval(self.SetAion.GetValue())
        top.vbeam = 0.
        top.vbeam_s = 0.
        derivqty()
        self.DoStep()

    def OnSetzionTextEnter(self, event):
        top.zion = eval(self.SetZion.GetValue())
        derivqty()
        self.DoStep()
      
    def DoStep(self,override=0):
        if self.dorunonchange or override:
            envexe()
            self.MakeEnvPlot()
        
    def MakeEnvPlot(self,y=None,x=None):
        if self.dofmabeforeplot: fma()
        if x is not None: self.xsave = x
        if y is not None: self.ysave = y
        try:
            for x,y in map(None,self.xsave,self.ysave):
                x = env.getpyobject(x)
                y = env.getpyobject(y)
                if x is not None and y is not None:
                  plg(y,x,color=self.plotcolor)
                  redraw()
        except AttributeError:
            pass










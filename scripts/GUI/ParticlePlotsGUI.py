#Boa:Dialog:ParticlePlotsGUI

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from warp import *

def create(parent):
    return ParticlePlotsGUI(parent)

[wxID_WXDIALOG2, wxID_WXDIALOG2CELLARRAY, wxID_WXDIALOG2CONTOURS, wxID_WXDIALOG2DENSITY, wxID_WXDIALOG2FMA, wxID_WXDIALOG2FMABEFOREPLOT, wxID_WXDIALOG2HCP, wxID_WXDIALOG2IZSLIDER, wxID_WXDIALOG2PALETTE, wxID_WXDIALOG2PARTICLES, wxID_WXDIALOG2PLOTCHOICE, wxID_WXDIALOG2PLOTREFRESH, wxID_WXDIALOG2PPTRACE, wxID_WXDIALOG2PPXPYP, wxID_WXDIALOG2PPXXP, wxID_WXDIALOG2PPXY, wxID_WXDIALOG2PPYYP, wxID_WXDIALOG2PPZX, wxID_WXDIALOG2PPZXP, wxID_WXDIALOG2PPZY, wxID_WXDIALOG2PPZYP, wxID_WXDIALOG2STATICLINE1, wxID_WXDIALOG2STATICLINE2, wxID_WXDIALOG2STATICLINE3, wxID_WXDIALOG2STATICLINE4, wxID_WXDIALOG2STATICLINE5, wxID_WXDIALOG2SURFACE] = map(lambda _init_ctrls: wxNewId(), range(27))

class ParticlePlotsGUI(wxDialog):
    def _init_utils(self):
        pass

    def _init_ctrls(self, prnt):
        wxDialog.__init__(self, id = wxID_WXDIALOG2, name = '', parent = prnt, pos = wxPoint(468, 342), size = wxSize(543, 205), style = wxDIALOG_MODELESS | wxDEFAULT_DIALOG_STYLE, title = 'Particle Plots')
        self._init_utils()
        self.SetClientSize(wxSize(543, 205))
        self.SetAutoLayout(true)

        self.ppxy = wxButton(id = wxID_WXDIALOG2PPXY, label = 'ppxy', name = 'ppxy', parent = self, pos = wxPoint(72, 8), size = wxSize(56, 22), style = 0)
        self.ppxy.SetHelpText('Make x-y plot')
        self.ppxy.SetToolTipString('Plots x-y')
        EVT_BUTTON(self.ppxy, wxID_WXDIALOG2PPXY, self.OnPpxyButton)

        self.ppxxp = wxButton(id = wxID_WXDIALOG2PPXXP, label = 'ppxxp', name = 'ppxxp', parent = self, pos = wxPoint(136, 8), size = wxSize(56, 22), style = 0)
        self.ppxxp.SetToolTipString("Plots x-x'")
        self.ppxxp.SetHelpText("Plots x-x'")
        EVT_BUTTON(self.ppxxp, wxID_WXDIALOG2PPXXP, self.OnPpxxpButton)

        self.ppyyp = wxButton(id = wxID_WXDIALOG2PPYYP, label = 'ppyyp', name = 'ppyyp', parent = self, pos = wxPoint(72, 40), size = wxSize(56, 22), style = 0)
        self.ppyyp.SetToolTipString("Plots y-y'")
        self.ppyyp.SetHelpText("Plots y-y'")
        EVT_BUTTON(self.ppyyp, wxID_WXDIALOG2PPYYP, self.OnPpyypButton)

        self.ppxpyp = wxButton(id = wxID_WXDIALOG2PPXPYP, label = 'ppxpyp', name = 'ppxpyp', parent = self, pos = wxPoint(136, 40), size = wxSize(56, 22), style = 0)
        self.ppxpyp.SetToolTipString("Plots x'-y'")
        self.ppxpyp.SetHelpText("Plots x'-y'")
        EVT_BUTTON(self.ppxpyp, wxID_WXDIALOG2PPXPYP, self.OnPpxpypButton)

        self.pptrace = wxButton(id = wxID_WXDIALOG2PPTRACE, label = 'pptrace', name = 'pptrace', parent = self, pos = wxPoint(72, 72), size = wxSize(120, 22), style = 0)
        self.pptrace.SetToolTipString('Plots trace-space')
        self.pptrace.SetHelpText('Plots trace-space')
        EVT_BUTTON(self.pptrace, wxID_WXDIALOG2PPTRACE, self.OnPptraceButton)

        self.izslider = wxSlider(id = wxID_WXDIALOG2IZSLIDER, maxValue = 100, minValue = 0, name = 'izslider', parent = self, point = wxPoint(72, 112), size = wxSize(112, 35), style = wxSL_LABELS | wxSL_AUTOTICKS | wxSL_HORIZONTAL, validator = wxDefaultValidator, value = 0)
        self.izslider.SetLabel('iz')
        self.izslider.SetToolTipString('Grid location')
        EVT_SLIDER(self.izslider, wxID_WXDIALOG2IZSLIDER, self.OnIzsliderSlider)

        self.fma = wxButton(id = wxID_WXDIALOG2FMA, label = 'fma', name = 'fma', parent = self, pos = wxPoint(8, 8), size = wxSize(40, 38), style = 0)
        self.fma.SetToolTipString('Frame Advance')
        EVT_BUTTON(self.fma, wxID_WXDIALOG2FMA, self.OnFmaButton)

        self.hcp = wxButton(id = wxID_WXDIALOG2HCP, label = 'hcp', name = 'hcp', parent = self, pos = wxPoint(8, 56), size = wxSize(40, 38), style = 0)
        self.hcp.SetToolTipString('Hard Copy')
        EVT_BUTTON(self.hcp, wxID_WXDIALOG2HCP, self.OnHcpButton)

        self.plotchoice = wxChoice(choices = ['iw','iz','iy','ix'], id = wxID_WXDIALOG2PLOTCHOICE, name = 'plotchoice', parent = self, pos = wxPoint(8, 112), size = wxSize(60, 25), style = 0, validator = wxDefaultValidator)
        self.plotchoice.SetToolTipString('Choose plot location')
        EVT_CHOICE(self.plotchoice, wxID_WXDIALOG2PLOTCHOICE, self.OnPlotchoiceChoice)

        self.Contours = wxRadioButton(id = wxID_WXDIALOG2CONTOURS, label = 'Contours', name = 'Contours', parent = self, pos = wxPoint(392, 24), size = wxSize(94, 24), style = 0)
        self.Contours.SetValue(true)
        self.Contours.SetToolTipString('Selects contour plots')
        EVT_RADIOBUTTON(self.Contours, wxID_WXDIALOG2CONTOURS, self.OnContoursRadiobutton)

        self.Density = wxRadioButton(id = wxID_WXDIALOG2DENSITY, label = 'Density', name = 'Density', parent = self, pos = wxPoint(392, 40), size = wxSize(94, 24), style = 0)
        self.Density.SetValue(false)
        self.Density.SetToolTipString('Selects density color plots')
        EVT_RADIOBUTTON(self.Density, wxID_WXDIALOG2DENSITY, self.OnDensityRadiobutton)

        self.Cellarray = wxRadioButton(id = wxID_WXDIALOG2CELLARRAY, label = 'Cellarray', name = 'Cellarray', parent = self, pos = wxPoint(392, 56), size = wxSize(94, 24), style = 0)
        self.Cellarray.SetValue(false)
        self.Cellarray.SetToolTipString('Selects density cellarray plots')
        EVT_RADIOBUTTON(self.Cellarray, wxID_WXDIALOG2CELLARRAY, self.OnCellarrayRadiobutton)

        self.Surface = wxRadioButton(id = wxID_WXDIALOG2SURFACE, label = 'Surface', name = 'Surface', parent = self, pos = wxPoint(392, 72), size = wxSize(94, 24), style = 0)
        self.Surface.SetValue(false)
        self.Surface.SetToolTipString('Selects surface plots')
        EVT_RADIOBUTTON(self.Surface, wxID_WXDIALOG2SURFACE, self.OnSurfaceRadiobutton)

        self.Palette = wxChoice(choices = ["earth","rainbow","gray","yarg","heat","ncar","cool","rainbowaf","stern","christmas"], id = wxID_WXDIALOG2PALETTE, name = 'Palette', parent = self, pos = wxPoint(392, 96), size = wxSize(112, 25), style = 0, validator = wxDefaultValidator)
        self.Palette.SetToolTipString('Selects palette')
        EVT_CHOICE(self.Palette, wxID_WXDIALOG2PALETTE, self.OnPaletteChoice)

        self.Particles = wxRadioButton(id = wxID_WXDIALOG2PARTICLES, label = 'Particles', name = 'Particles', parent = self, pos = wxPoint(392, 8), size = wxSize(72, 24), style = 0)
        self.Particles.SetValue(true)
        self.Particles.SetToolTipString('Selects particles plots')
        EVT_RADIOBUTTON(self.Particles, wxID_WXDIALOG2PARTICLES, self.OnParticlesRadiobutton)

        self.PlotRefresh = wxCheckBox(id = wxID_WXDIALOG2PLOTREFRESH, label = 'Refresh on Change', name = 'PlotRefresh', parent = self, pos = wxPoint(8, 152), size = wxSize(128, 24), style = 0)
        self.PlotRefresh.SetValue(false)
        self.PlotRefresh.SetToolTipString('Plot will refresh on change')
        EVT_CHECKBOX(self.PlotRefresh, wxID_WXDIALOG2PLOTREFRESH, self.OnPlotrefreshCheckbox)

        self.staticLine1 = wxStaticLine(id = wxID_WXDIALOG2STATICLINE1, name = 'staticLine1', parent = self, pos = wxPoint(192, 0), size = wxSize(16, 208), style = wxLI_VERTICAL)

        self.staticLine2 = wxStaticLine(id = wxID_WXDIALOG2STATICLINE2, name = 'staticLine2', parent = self, pos = wxPoint(376, 0), size = wxSize(20, 208), style = wxLI_VERTICAL)

        self.staticLine3 = wxStaticLine(id = wxID_WXDIALOG2STATICLINE3, name = 'staticLine3', parent = self, pos = wxPoint(0, 144), size = wxSize(200, 16), style = wxLI_HORIZONTAL)

        self.staticLine4 = wxStaticLine(id = wxID_WXDIALOG2STATICLINE4, name = 'staticLine4', parent = self, pos = wxPoint(48, 0), size = wxSize(20, 104), style = wxLI_VERTICAL)

        self.staticLine5 = wxStaticLine(id = wxID_WXDIALOG2STATICLINE5, name = 'staticLine5', parent = self, pos = wxPoint(0, 96), size = wxSize(200, 16), style = wxLI_HORIZONTAL)

        self.ppzx = wxButton(id = wxID_WXDIALOG2PPZX, label = 'ppzx', name = 'ppzx', parent = self, pos = wxPoint(216, 8), size = wxSize(56, 22), style = 0)
        self.ppzx.SetToolTipString('Plots z-x')
        EVT_BUTTON(self.ppzx, wxID_WXDIALOG2PPZX, self.OnPpzxButton)

        self.ppzxp = wxButton(id = wxID_WXDIALOG2PPZXP, label = 'ppzxp', name = 'ppzxp', parent = self, pos = wxPoint(216, 40), size = wxSize(56, 22), style = 0)
        self.ppzxp.SetToolTipString("Plots z-x'")
        EVT_BUTTON(self.ppzxp, wxID_WXDIALOG2PPZXP, self.OnPpzxpButton)

        self.ppzyp = wxButton(id = wxID_WXDIALOG2PPZYP, label = 'ppzyp', name = 'ppzyp', parent = self, pos = wxPoint(280, 40), size = wxSize(56, 22), style = 0)
        self.ppzyp.SetToolTipString("Plots z-y'")
        EVT_BUTTON(self.ppzyp, wxID_WXDIALOG2PPZYP, self.OnPpzypButton)

        self.ppzy = wxButton(id = wxID_WXDIALOG2PPZY, label = 'ppzy', name = 'ppzy', parent = self, pos = wxPoint(280, 8), size = wxSize(56, 22), style = 0)
        self.ppzy.SetToolTipString('Plots z-y')
        EVT_BUTTON(self.ppzy, wxID_WXDIALOG2PPZY, self.OnPpzyButton)

        self.fmabeforeplot = wxCheckBox(id = wxID_WXDIALOG2FMABEFOREPLOT, label = 'Frame Advance before plot', name = 'fmabeforeplot', parent = self, pos = wxPoint(8, 176), size = wxSize(184, 24), style = 0)
        self.fmabeforeplot.SetValue(true)
        self.fmabeforeplot.SetToolTipString('When checked, do frame advance before each plot')
        EVT_CHECKBOX(self.fmabeforeplot, wxID_WXDIALOG2FMABEFOREPLOT, self.OnFmabeforeplotCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.plotchoiceslidervalue = 0
        self.plotchoicekw = {'iw':self.plotchoiceslidervalue}
        self.plottypekw = {}
        self.doplotrefreshonchange = 0
        self.dofmabeforeplot = 1

    def OnPpxyButton(self, event):
        self.currentplot = ppxy
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpxxpButton(self, event):
        self.currentplot = ppxxp
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpyypButton(self, event):
        self.currentplot = ppyyp
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpxpypButton(self, event):
        self.currentplot = ppxpyp
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPptraceButton(self, event):
        self.currentplot = pptrace
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpzxButton(self, event):
        self.currentplot = ppzx
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpzxpButton(self, event):
        self.currentplot = ppzxp
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpzypButton(self, event):
        self.currentplot = ppzyp
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnPpzyButton(self, event):
        self.currentplot = ppzy
        self.MakeParticlePlot(1,self.dofmabeforeplot)

    def OnIzsliderSlider(self, event):
        self.plotchoiceslidervalue = self.izslider.GetValue()
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnPlotchoiceChoice(self, event):
        plotchoice = self.plotchoice.GetStringSelection()
        if plotchoice == 'iw':
          self.izslider.SetRange(0,top.nzwind)
          self.plotchoicekw = {'iw':self.plotchoiceslidervalue}
        elif plotchoice == 'iz':
          self.izslider.SetRange(0,w3d.nz)
          self.plotchoicekw = {'iz':self.plotchoiceslidervalue}
        elif plotchoice == 'ix':
          self.izslider.SetRange(0,w3d.nx)
          self.plotchoicekw = {'ix':self.plotchoiceslidervalue}
        elif plotchoice == 'iy':
          self.izslider.SetRange(0,w3d.ny)
          self.plotchoicekw = {'iy':self.plotchoiceslidervalue}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnParticlesRadiobutton(self, event):
        self.plottypekw = {}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnContoursRadiobutton(self, event):
        self.plottypekw = {'contours':10}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnDensityRadiobutton(self, event):
        self.plottypekw = {'color':'density'}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnCellarrayRadiobutton(self, event):
        self.plottypekw = {'cellarray':1}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnSurfaceRadiobutton(self, event):
        self.plottypekw = {'surface':1}
        self.MakeParticlePlot(self.doplotrefreshonchange,1)

    def OnPaletteChoice(self, event):
        newpalette = self.Palette.GetStringSelection()
        try:
          palette(newpalette+".gp")
        except gist.error:
          pass

    def OnFmaButton(self, event):
        fma()

    def OnHcpButton(self, event):
        hcp()

    def OnPlotrefreshCheckbox(self, event):
        self.doplotrefreshonchange = self.PlotRefresh.GetValue()

    def OnFmabeforeplotCheckbox(self, event):
        self.dofmabeforeplot = self.fmabeforeplot.GetValue()

    def OnAllowZoom(self, event):
        ygdispatch()

    def MakeParticlePlot(self,refresh,dofma):
        if not refresh: return
        if dofma: fma()
        kw = self.plotchoicekw.copy()
        kw.update(self.plottypekw)
        apply(self.currentplot,[],kw)
        redraw()


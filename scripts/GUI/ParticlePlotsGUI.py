#Boa:FramePanel:ParticlePlotsGUI

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from warp import *
l_opendx=1
try:
    from pyOpenDX import *
except:
    l_opendx=0

[wxID_PARTICLEPLOTSGUI, wxID_PARTICLEPLOTSGUICELLARRAY, 
 wxID_PARTICLEPLOTSGUICONTOURS, wxID_PARTICLEPLOTSGUIDENSITY, 
 wxID_PARTICLEPLOTSGUIFMABEFOREPLOT, wxID_PARTICLEPLOTSGUIIZSLIDER, 
 wxID_PARTICLEPLOTSGUIPALETTE, wxID_PARTICLEPLOTSGUIPARTICLES, 
 wxID_PARTICLEPLOTSGUIPLOTCHOICE, wxID_PARTICLEPLOTSGUIPLOTREFRESH, 
 wxID_PARTICLEPLOTSGUIPPTRACE, wxID_PARTICLEPLOTSGUIPPXPYP, 
 wxID_PARTICLEPLOTSGUIPPXXP, wxID_PARTICLEPLOTSGUIPPXY, 
 wxID_PARTICLEPLOTSGUIPPXYZ, wxID_PARTICLEPLOTSGUIPPXYZ_V, 
 wxID_PARTICLEPLOTSGUIPPYYP, wxID_PARTICLEPLOTSGUIPPZX, 
 wxID_PARTICLEPLOTSGUIPPZXP, wxID_PARTICLEPLOTSGUIPPZY, 
 wxID_PARTICLEPLOTSGUIPPZYP, wxID_PARTICLEPLOTSGUISTATICLINE1, 
 wxID_PARTICLEPLOTSGUISTATICLINE2, wxID_PARTICLEPLOTSGUISTATICLINE3, 
 wxID_PARTICLEPLOTSGUISTATICLINE5, wxID_PARTICLEPLOTSGUISTATICLINE7, 
 wxID_PARTICLEPLOTSGUISURFACE, 
] = map(lambda _init_ctrls: wxNewId(), range(27))

class ParticlePlotsGUI(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PARTICLEPLOTSGUI, name='', parent=prnt,
              pos=wxPoint(428, 291), size=wxSize(479, 227),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(471, 203))

        self.ppxy = wxButton(id=wxID_PARTICLEPLOTSGUIPPXY, label='ppxy',
              name='ppxy', parent=self, pos=wxPoint(8, 8), size=wxSize(56, 22),
              style=0)
        self.ppxy.SetHelpText('Make x-y plot')
        self.ppxy.SetToolTipString('Plots x-y')
        EVT_BUTTON(self.ppxy, wxID_PARTICLEPLOTSGUIPPXY, self.OnPpxyButton)

        self.ppxxp = wxButton(id=wxID_PARTICLEPLOTSGUIPPXXP, label='ppxxp',
              name='ppxxp', parent=self, pos=wxPoint(72, 8), size=wxSize(56,
              22), style=0)
        self.ppxxp.SetToolTipString("Plots x-x'")
        self.ppxxp.SetHelpText("Plots x-x'")
        EVT_BUTTON(self.ppxxp, wxID_PARTICLEPLOTSGUIPPXXP, self.OnPpxxpButton)

        self.ppyyp = wxButton(id=wxID_PARTICLEPLOTSGUIPPYYP, label='ppyyp',
              name='ppyyp', parent=self, pos=wxPoint(8, 40), size=wxSize(56,
              22), style=0)
        self.ppyyp.SetToolTipString("Plots y-y'")
        self.ppyyp.SetHelpText("Plots y-y'")
        EVT_BUTTON(self.ppyyp, wxID_PARTICLEPLOTSGUIPPYYP, self.OnPpyypButton)

        self.ppxpyp = wxButton(id=wxID_PARTICLEPLOTSGUIPPXPYP, label='ppxpyp',
              name='ppxpyp', parent=self, pos=wxPoint(72, 40), size=wxSize(56,
              22), style=0)
        self.ppxpyp.SetToolTipString("Plots x'-y'")
        self.ppxpyp.SetHelpText("Plots x'-y'")
        EVT_BUTTON(self.ppxpyp, wxID_PARTICLEPLOTSGUIPPXPYP,
              self.OnPpxpypButton)

        self.pptrace = wxButton(id=wxID_PARTICLEPLOTSGUIPPTRACE,
              label='pptrace', name='pptrace', parent=self, pos=wxPoint(8, 72),
              size=wxSize(120, 22), style=0)
        self.pptrace.SetToolTipString('Plots trace-space')
        self.pptrace.SetHelpText('Plots trace-space')
        EVT_BUTTON(self.pptrace, wxID_PARTICLEPLOTSGUIPPTRACE,
              self.OnPptraceButton)

        self.izslider = wxSlider(id=wxID_PARTICLEPLOTSGUIIZSLIDER, maxValue=100,
              minValue=0, name='izslider', parent=self, point=wxPoint(115, 104),
              size=wxSize(172, 35),
              style=wxSL_LABELS | wxSL_AUTOTICKS | wxSL_HORIZONTAL,
              validator=wxDefaultValidator, value=0)
        self.izslider.SetLabel('iz')
        self.izslider.SetToolTipString('Grid location')
        EVT_LEFT_UP(self.izslider, self.OnIzsliderLeftUp)

        self.plotchoice = wxChoice(choices=['iw', 'iz', 'iy', 'ix'],
              id=wxID_PARTICLEPLOTSGUIPLOTCHOICE, name='plotchoice',
              parent=self, pos=wxPoint(8, 112), size=wxSize(60, 21), style=0,
              validator=wxDefaultValidator)
        self.plotchoice.SetToolTipString('Choose plot location')
        EVT_CHOICE(self.plotchoice, wxID_PARTICLEPLOTSGUIPLOTCHOICE,
              self.OnPlotchoiceChoice)

        self.Contours = wxRadioButton(id=wxID_PARTICLEPLOTSGUICONTOURS,
              label='Contours', name='Contours', parent=self, pos=wxPoint(354,
              24), size=wxSize(94, 24), style=0)
        self.Contours.SetValue(true)
        self.Contours.SetToolTipString('Selects contour plots')
        EVT_RADIOBUTTON(self.Contours, wxID_PARTICLEPLOTSGUICONTOURS,
              self.OnContoursRadiobutton)

        self.Density = wxRadioButton(id=wxID_PARTICLEPLOTSGUIDENSITY,
              label='Density', name='Density', parent=self, pos=wxPoint(354,
              40), size=wxSize(94, 24), style=0)
        self.Density.SetValue(false)
        self.Density.SetToolTipString('Selects density color plots')
        EVT_RADIOBUTTON(self.Density, wxID_PARTICLEPLOTSGUIDENSITY,
              self.OnDensityRadiobutton)

        self.Cellarray = wxRadioButton(id=wxID_PARTICLEPLOTSGUICELLARRAY,
              label='Cellarray', name='Cellarray', parent=self, pos=wxPoint(354,
              56), size=wxSize(94, 24), style=0)
        self.Cellarray.SetValue(false)
        self.Cellarray.SetToolTipString('Selects density cellarray plots')
        EVT_RADIOBUTTON(self.Cellarray, wxID_PARTICLEPLOTSGUICELLARRAY,
              self.OnCellarrayRadiobutton)

        self.Surface = wxRadioButton(id=wxID_PARTICLEPLOTSGUISURFACE,
              label='Surface', name='Surface', parent=self, pos=wxPoint(354,
              72), size=wxSize(94, 24), style=0)
        self.Surface.SetValue(false)
        self.Surface.SetToolTipString('Selects surface plots')
        EVT_RADIOBUTTON(self.Surface, wxID_PARTICLEPLOTSGUISURFACE,
              self.OnSurfaceRadiobutton)

        self.Palette = wxChoice(choices=["earth", "rainbow", "gray", "yarg",
              "heat", "ncar", "cool", "rainbowaf", "stern", "christmas"],
              id=wxID_PARTICLEPLOTSGUIPALETTE, name='Palette', parent=self,
              pos=wxPoint(354, 96), size=wxSize(112, 21), style=0,
              validator=wxDefaultValidator)
        self.Palette.SetToolTipString('Selects palette')
        EVT_CHOICE(self.Palette, wxID_PARTICLEPLOTSGUIPALETTE,
              self.OnPaletteChoice)

        self.Particles = wxRadioButton(id=wxID_PARTICLEPLOTSGUIPARTICLES,
              label='Particles', name='Particles', parent=self, pos=wxPoint(354,
              8), size=wxSize(72, 24), style=0)
        self.Particles.SetValue(true)
        self.Particles.SetToolTipString('Selects particles plots')
        EVT_RADIOBUTTON(self.Particles, wxID_PARTICLEPLOTSGUIPARTICLES,
              self.OnParticlesRadiobutton)

        self.PlotRefresh = wxCheckBox(id=wxID_PARTICLEPLOTSGUIPLOTREFRESH,
              label='Refresh on Change', name='PlotRefresh', parent=self,
              pos=wxPoint(8, 144), size=wxSize(128, 24), style=0)
        self.PlotRefresh.SetValue(false)
        self.PlotRefresh.SetToolTipString('Plot will refresh on change')
        EVT_CHECKBOX(self.PlotRefresh, wxID_PARTICLEPLOTSGUIPLOTREFRESH,
              self.OnPlotrefreshCheckbox)

        self.staticLine1 = wxStaticLine(id=wxID_PARTICLEPLOTSGUISTATICLINE1,
              name='staticLine1', parent=self, pos=wxPoint(136, 0),
              size=wxSize(2, 104), style=wxLI_VERTICAL)

        self.staticLine2 = wxStaticLine(id=wxID_PARTICLEPLOTSGUISTATICLINE2,
              name='staticLine2', parent=self, pos=wxPoint(344, -8),
              size=wxSize(2, 208), style=wxLI_VERTICAL)

        self.staticLine3 = wxStaticLine(id=wxID_PARTICLEPLOTSGUISTATICLINE3,
              name='staticLine3', parent=self, pos=wxPoint(0, 140),
              size=wxSize(346, 2), style=wxLI_HORIZONTAL)

        self.staticLine5 = wxStaticLine(id=wxID_PARTICLEPLOTSGUISTATICLINE5,
              name='staticLine5', parent=self, pos=wxPoint(0, 104),
              size=wxSize(344, 2), style=wxLI_HORIZONTAL)

        self.ppzx = wxButton(id=wxID_PARTICLEPLOTSGUIPPZX, label='ppzx',
              name='ppzx', parent=self, pos=wxPoint(144, 8), size=wxSize(56,
              22), style=0)
        self.ppzx.SetToolTipString('Plots z-x')
        EVT_BUTTON(self.ppzx, wxID_PARTICLEPLOTSGUIPPZX, self.OnPpzxButton)

        self.ppzxp = wxButton(id=wxID_PARTICLEPLOTSGUIPPZXP, label='ppzxp',
              name='ppzxp', parent=self, pos=wxPoint(144, 40), size=wxSize(56,
              22), style=0)
        self.ppzxp.SetToolTipString("Plots z-x'")
        EVT_BUTTON(self.ppzxp, wxID_PARTICLEPLOTSGUIPPZXP, self.OnPpzxpButton)

        self.ppzyp = wxButton(id=wxID_PARTICLEPLOTSGUIPPZYP, label='ppzyp',
              name='ppzyp', parent=self, pos=wxPoint(208, 40), size=wxSize(56,
              22), style=0)
        self.ppzyp.SetToolTipString("Plots z-y'")
        EVT_BUTTON(self.ppzyp, wxID_PARTICLEPLOTSGUIPPZYP, self.OnPpzypButton)

        self.ppzy = wxButton(id=wxID_PARTICLEPLOTSGUIPPZY, label='ppzy',
              name='ppzy', parent=self, pos=wxPoint(208, 8), size=wxSize(56,
              22), style=0)
        self.ppzy.SetToolTipString('Plots z-y')
        EVT_BUTTON(self.ppzy, wxID_PARTICLEPLOTSGUIPPZY, self.OnPpzyButton)

        self.fmabeforeplot = wxCheckBox(id=wxID_PARTICLEPLOTSGUIFMABEFOREPLOT,
              label='Frame Advance before plot', name='fmabeforeplot',
              parent=self, pos=wxPoint(8, 168), size=wxSize(184, 24), style=0)
        self.fmabeforeplot.SetValue(true)
        self.fmabeforeplot.SetToolTipString('When checked, do frame advance before each plot')
        EVT_CHECKBOX(self.fmabeforeplot, wxID_PARTICLEPLOTSGUIFMABEFOREPLOT,
              self.OnFmabeforeplotCheckbox)

        self.staticLine7 = wxStaticLine(id=wxID_PARTICLEPLOTSGUISTATICLINE7,
              name='staticLine7', parent=self, pos=wxPoint(272, 0),
              size=wxSize(2, 104), style=0)

        self.ppxyz = wxButton(id=wxID_PARTICLEPLOTSGUIPPXYZ, label='ppxyz',
              name='ppxyz', parent=self, pos=wxPoint(280, 8), size=wxSize(56,
              22), style=0)
        EVT_BUTTON(self.ppxyz, wxID_PARTICLEPLOTSGUIPPXYZ, self.OnPpxyzButton)

        self.ppxyz_v = wxButton(id=wxID_PARTICLEPLOTSGUIPPXYZ_V,
              label='ppxyz-v', name='ppxyz_v', parent=self, pos=wxPoint(280,
              40), size=wxSize(56, 22), style=0)
        EVT_BUTTON(self.ppxyz_v, wxID_PARTICLEPLOTSGUIPPXYZ_V,
              self.OnPpxyz_vButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
#        parent.AddPage(imageId=-1, page=self, select=True, text='PPlots')
        self.plotchoiceslidervalue = 0
        self.plottypekw = {}
        self.doplotrefreshonchange = 0
        self.dofmabeforeplot = 1
        self.iw = 0
        self.ix = 0
        self.iy = 0
        self.iz = 0
        self.plotchoicekw = 'iw'
        self.izslider.SetRange(0,top.nzwind)
        self.izslider.SetValue(self.iw)
        self.plotchoicekw = 'iw'

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

#    def OnIzsliderSlider(self, event):
    def OnIzsliderLeftUp(self, event):
        self.plotchoiceslidervalue = self.izslider.GetValue()
        self.MakeParticlePlot(self.doplotrefreshonchange,self.dofmabeforeplot)

    def OnPlotchoiceChoice(self, event):
        plotchoice = self.plotchoice.GetStringSelection()
        if plotchoice == 'iw':
          self.izslider.SetRange(0,top.nzwind)
          self.izslider.SetValue(self.iw)
          self.plotchoicekw = 'iw'
        elif plotchoice == 'iz':
          self.izslider.SetRange(0,w3d.nz)
          self.izslider.SetValue(self.iz)
          self.plotchoicekw = 'iz'
        elif plotchoice == 'ix':
          self.izslider.SetRange(0,w3d.nx)
          self.izslider.SetValue(self.ix)
          self.plotchoicekw = 'ix'
        elif plotchoice == 'iy':
          self.izslider.SetRange(0,w3d.ny)
          self.izslider.SetValue(self.iy)
          self.plotchoicekw = 'iy'
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
        if dofma and not self.plottypekw.has_key('surface'): fma()
        kw = {self.plotchoicekw:self.plotchoiceslidervalue}
        kw.update(self.plottypekw)
        apply(self.currentplot,[],kw)
        redraw()

    def OnPpxyzButton(self, event):
        if l_opendx:ppxyz()
        event.Skip()

    def OnPpxyz_vButton(self, event):
        if l_opendx:ppxyzvxvyvz()
        event.Skip()


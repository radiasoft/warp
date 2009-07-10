#Boa:Dialog:LatticeGUI

from wx import *
from wx.lib.anchors import LayoutAnchors
from warp import *
import sortlattice
import lattice
import __main__

def create(parent):
    return LatticeGUI(parent)

[wxID_LATTICEGUI, wxID_LATTICEGUIDOSTEP, wxID_LATTICEGUIDRFT, 
 wxID_LATTICEGUIDRFTAPLABEL, wxID_LATTICEGUIDRFTAPUNITS, 
 wxID_LATTICEGUIDRFTLABEL, wxID_LATTICEGUIDRFTZELABEL, 
 wxID_LATTICEGUIDRFTZEUNITS, wxID_LATTICEGUIDRFTZSLABEL, 
 wxID_LATTICEGUIDRFTZSUNITS, wxID_LATTICEGUIELEMENTNUM, 
 wxID_LATTICEGUIELEMENTSPIN, wxID_LATTICEGUIELEMNUMLABEL1, 
 wxID_LATTICEGUIELEMNUMLABEL2, wxID_LATTICEGUIGETDRFTAP, 
 wxID_LATTICEGUIGETDRFTZE, wxID_LATTICEGUIGETDRFTZS, wxID_LATTICEGUIGETQUADDB, 
 wxID_LATTICEGUIGETQUADDE, wxID_LATTICEGUIGETQUADZE, wxID_LATTICEGUIGETQUADZS, 
 wxID_LATTICEGUIMAKEUNIQUE, wxID_LATTICEGUIQUAD, wxID_LATTICEGUIQUADDBLABEL, 
 wxID_LATTICEGUIQUADDBUNITS, wxID_LATTICEGUIQUADDELABEL, 
 wxID_LATTICEGUIQUADDEUNITS, wxID_LATTICEGUIQUADLABEL, 
 wxID_LATTICEGUIQUADZELABEL, wxID_LATTICEGUIQUADZEUNITS, 
 wxID_LATTICEGUIQUADZSLABEL, wxID_LATTICEGUIQUADZSUNITS, 
 wxID_LATTICEGUISETELEMENTNUM, wxID_LATTICEGUISETMADLATTICE, 
 wxID_LATTICEGUISETMADLATTICELABEL, 
] = map(lambda _init_ctrls: wx.NewId(), range(35))

class LatticeGUI(wx.Dialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Dialog.__init__(self, id=wxID_LATTICEGUI, name='LatticeGUI',
              parent=prnt, pos=wx.Point(493, 320), size=wx.Size(389, 277),
              style=wx.DEFAULT_DIALOG_STYLE, title='Lattice Editor')
        self._init_utils()
        self.SetClientSize(wx.Size(389, 277))
        self.SetToolTipString('Lattice editor')

        self.ElementNum = wx.Panel(id=wxID_LATTICEGUIELEMENTNUM,
              name='ElementNum', parent=self, pos=wx.Point(8, 8),
              size=wx.Size(136, 264), style=wx.RAISED_BORDER | wx.TAB_TRAVERSAL)

        self.SetElementNum = wx.TextCtrl(id=wxID_LATTICEGUISETELEMENTNUM,
              name='SetElementNum', parent=self.ElementNum, pos=wx.Point(8, 48),
              size=wx.Size(40, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='0')
        self.SetElementNum.SetToolTipString('Sets element number')
        EVT_TEXT_ENTER(self.SetElementNum, wxID_LATTICEGUISETELEMENTNUM,
              self.OnSetelementnumTextEnter)

        self.ElemNumLabel1 = wx.StaticText(id=wxID_LATTICEGUIELEMNUMLABEL1,
              label='Element', name='ElemNumLabel1', parent=self.ElementNum,
              pos=wx.Point(8, 8), size=wx.Size(44, 16), style=0)

        self.ElemNumLabel2 = wx.StaticText(id=wxID_LATTICEGUIELEMNUMLABEL2,
              label='Number', name='ElemNumLabel2', parent=self.ElementNum,
              pos=wx.Point(8, 24), size=wx.Size(43, 16), style=0)

        self.ElementSpin = wx.SpinButton(id=wxID_LATTICEGUIELEMENTSPIN,
              name='ElementSpin', parent=self.ElementNum, pos=wx.Point(56, 40),
              size=wx.Size(15, 34),
              style=wx.DOUBLE_BORDER | wx.SP_HORIZONTAL | wx.SIMPLE_BORDER)
        self.ElementSpin.SetToolTipString('Change element number')
        EVT_SPIN_DOWN(self.ElementSpin, wxID_LATTICEGUIELEMENTSPIN,
              self.OnElementspinSpinDown)
        EVT_SPIN_UP(self.ElementSpin, wxID_LATTICEGUIELEMENTSPIN,
              self.OnElementspinSpinUp)

        self.Quad = wx.Panel(id=wxID_LATTICEGUIQUAD, name='Quad', parent=self,
              pos=wx.Point(144, 8), size=wx.Size(240, 264),
              style=wx.SUNKEN_BORDER | wx.TAB_TRAVERSAL)
        self.Quad.Show(false)

        self.QuadLabel = wx.StaticText(id=wxID_LATTICEGUIQUADLABEL,
              label='Hard edged quadrupole', name='QuadLabel', parent=self.Quad,
              pos=wx.Point(42, 0), size=wx.Size(151, 18),
              style=wx.SIMPLE_BORDER | wx.ALIGN_CENTRE)
        self.QuadLabel.Center(wx.HORIZONTAL)
        self.QuadLabel.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.NORMAL, false,
              ''))

        self.QuadzsLabel = wx.StaticText(id=wxID_LATTICEGUIQUADZSLABEL,
              label='Z start', name='QuadzsLabel', parent=self.Quad,
              pos=wx.Point(8, 27), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.QuadzeLabel = wx.StaticText(id=wxID_LATTICEGUIQUADZELABEL,
              label='Z end', name='QuadzeLabel', parent=self.Quad,
              pos=wx.Point(8, 51), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.QuaddeLabel = wx.StaticText(id=wxID_LATTICEGUIQUADDELABEL,
              label='E gradient', name='QuaddeLabel', parent=self.Quad,
              pos=wx.Point(8, 75), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.QuaddbLabel = wx.StaticText(id=wxID_LATTICEGUIQUADDBLABEL,
              label='B gradient', name='QuaddbLabel', parent=self.Quad,
              pos=wx.Point(8, 99), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.QuadzsUnits = wx.StaticText(id=wxID_LATTICEGUIQUADZSUNITS,
              label='meters', name='QuadzsUnits', parent=self.Quad,
              pos=wx.Point(152, 27), size=wx.Size(36, 16), style=0)

        self.QuadzeUnits = wx.StaticText(id=wxID_LATTICEGUIQUADZEUNITS,
              label='meters', name='QuadzeUnits', parent=self.Quad,
              pos=wx.Point(152, 51), size=wx.Size(36, 16), style=0)

        self.QuaddeUnits = wx.StaticText(id=wxID_LATTICEGUIQUADDEUNITS,
              label='V/m^2', name='QuaddeUnits', parent=self.Quad,
              pos=wx.Point(152, 75), size=wx.Size(35, 16), style=0)

        self.QuaddbUnits = wx.StaticText(id=wxID_LATTICEGUIQUADDBUNITS,
              label='B/m^2', name='QuaddbUnits', parent=self.Quad,
              pos=wx.Point(152, 99), size=wx.Size(34, 16), style=0)

        self.GetQuadzs = wx.TextCtrl(id=wxID_LATTICEGUIGETQUADZS,
              name='GetQuadzs', parent=self.Quad, pos=wx.Point(72, 24),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetQuadzs, wxID_LATTICEGUIGETQUADZS,
              self.OnGetelemzsTextEnter)

        self.GetQuadze = wx.TextCtrl(id=wxID_LATTICEGUIGETQUADZE,
              name='GetQuadze', parent=self.Quad, pos=wx.Point(72, 48),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetQuadze, wxID_LATTICEGUIGETQUADZE,
              self.OnGetelemzeTextEnter)

        self.GetQuadde = wx.TextCtrl(id=wxID_LATTICEGUIGETQUADDE,
              name='GetQuadde', parent=self.Quad, pos=wx.Point(72, 72),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetQuadde, wxID_LATTICEGUIGETQUADDE,
              self.OnGetQuaddeTextEnter)

        self.GetQuaddb = wx.TextCtrl(id=wxID_LATTICEGUIGETQUADDB,
              name='GetQuaddb', parent=self.Quad, pos=wx.Point(72, 96),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetQuaddb, wxID_LATTICEGUIGETQUADDB,
              self.OnGetQuaddbTextEnter)

        self.Drft = wx.Panel(id=wxID_LATTICEGUIDRFT, name='Drft', parent=self,
              pos=wx.Point(144, 8), size=wx.Size(240, 264),
              style=wx.SUNKEN_BORDER | wx.TAB_TRAVERSAL)
        self.Drft.SetToolTipString('Drift elements')
        self.Drft.Show(False)

        self.DrftLabel = wx.StaticText(id=wxID_LATTICEGUIDRFTLABEL,
              label='Drift', name='DrftLabel', parent=self.Drft,
              pos=wx.Point(105, 0), size=wx.Size(25, 18), style=wx.ALIGN_CENTRE)
        self.DrftLabel.Center(wx.HORIZONTAL)
        self.DrftLabel.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.NORMAL, false,
              ''))

        self.DrftzsLabel = wx.StaticText(id=wxID_LATTICEGUIDRFTZSLABEL,
              label='Z start', name='DrftzsLabel', parent=self.Drft,
              pos=wx.Point(8, 27), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.DrftzeLabel = wx.StaticText(id=wxID_LATTICEGUIDRFTZELABEL,
              label='Z end', name='DrftzeLabel', parent=self.Drft,
              pos=wx.Point(8, 51), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.DrftapLabel = wx.StaticText(id=wxID_LATTICEGUIDRFTAPLABEL,
              label='Aperture', name='DrftapLabel', parent=self.Drft,
              pos=wx.Point(8, 75), size=wx.Size(60, 16), style=wx.ALIGN_RIGHT)

        self.DrftzsUnits = wx.StaticText(id=wxID_LATTICEGUIDRFTZSUNITS,
              label='meters', name='DrftzsUnits', parent=self.Drft,
              pos=wx.Point(152, 27), size=wx.Size(36, 16), style=0)

        self.DrftzeUnits = wx.StaticText(id=wxID_LATTICEGUIDRFTZEUNITS,
              label='meters', name='DrftzeUnits', parent=self.Drft,
              pos=wx.Point(152, 51), size=wx.Size(36, 16), style=0)

        self.DrftapUnits = wx.StaticText(id=wxID_LATTICEGUIDRFTAPUNITS,
              label='meters', name='DrftapUnits', parent=self.Drft,
              pos=wx.Point(152, 75), size=wx.Size(36, 16), style=0)

        self.GetDrftzs = wx.TextCtrl(id=wxID_LATTICEGUIGETDRFTZS,
              name='GetDrftzs', parent=self.Drft, pos=wx.Point(72, 24),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetDrftzs, wxID_LATTICEGUIGETDRFTZS,
              self.OnGetelemzsTextEnter)

        self.GetDrftze = wx.TextCtrl(id=wxID_LATTICEGUIGETDRFTZE,
              name='GetDrftze', parent=self.Drft, pos=wx.Point(72, 48),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetDrftze, wxID_LATTICEGUIGETDRFTZE,
              self.OnGetelemzeTextEnter)

        self.GetDrftap = wx.TextCtrl(id=wxID_LATTICEGUIGETDRFTAP,
              name='GetDrftap', parent=self.Drft, pos=wx.Point(72, 72),
              size=wx.Size(80, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        EVT_TEXT_ENTER(self.GetDrftap, wxID_LATTICEGUIGETDRFTAP,
              self.OnGetelemapTextEnter)

        self.DoStep = wx.CheckBox(id=wxID_LATTICEGUIDOSTEP,
              label='Step on change', name='DoStep', parent=self.ElementNum,
              pos=wx.Point(8, 96), size=wx.Size(120, 24), style=wx.TAB_TRAVERSAL)
        self.DoStep.SetValue(false)
        self.DoStep.SetHelpText('When checked, execute step command on a change.')
        self.DoStep.SetToolTipString('Turns on code calculation on change')
        EVT_CHECKBOX(self.DoStep, wxID_LATTICEGUIDOSTEP, self.OnDostepCheckbox)

        self.SetMADLattice = wx.TextCtrl(id=wxID_LATTICEGUISETMADLATTICE,
              name='SetMADLattice', parent=self.ElementNum, pos=wx.Point(4, 168),
              size=wx.Size(104, 22), style=wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        self.SetMADLattice.SetToolTipString('Specify MAD lattice to use')
        EVT_TEXT_ENTER(self.SetMADLattice, wxID_LATTICEGUISETMADLATTICE,
              self.OnSetmadlatticeTextEnter)

        self.SetMADLatticeLabel = wx.StaticText(id=wxID_LATTICEGUISETMADLATTICELABEL,
              label='Set MAD lattice', name='SetMADLatticeLabel',
              parent=self.ElementNum, pos=wx.Point(8, 150), size=wx.Size(88, 16),
              style=0)

        self.MakeUnique = wx.CheckBox(id=wxID_LATTICEGUIMAKEUNIQUE,
              label='Unique elements', name='MakeUnique',
              parent=self.ElementNum, pos=wx.Point(4, 192), size=wx.Size(116, 24),
              style=wx.TAB_TRAVERSAL)
        self.MakeUnique.SetValue(False)
        self.MakeUnique.SetToolTipString('Make all elements of the lattice unique')
        EVT_CHECKBOX(self.MakeUnique, wxID_LATTICEGUIMAKEUNIQUE,
              self.OnMakeuniqueCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.elemnum = 0
        self.sortedelems = sortlattice.sortlattice()
        self.madlattice = None
        self.MADelementsunique = self.MakeUnique.GetValue()
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

    def GetCurrentElement(self):
        if self.madlattice:
            return self.madlattice[self.elemnum]
        else:
            return self.sortedelems[self.elemnum]

    def GetElemValue(self,val,valtype,valformat):
        el = self.GetCurrentElement()
        ctl = self.__dict__["Get"+el.type+val]
        try:
            value = valtype(ctl.GetValue())
            el.setattr(val,value)
            if self.madlattice:
                self.madlattice.derivedquantities()
                self.madlattice.setextent(0.)
                lattice.madtowarp(self.madlattice)
            self.Step()
        except ValueError:
            value = el.getattr(val)
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

    def OnGetQuaddeTextEnter(self, event):
        self.GetElemValue('de',float,"%e")

    def OnGetQuaddbTextEnter(self, event):
        self.GetElemFloatValue('db')

    def UpdatePanel(self):
        if self.shownelem is not None: self.shownelem.Show(0)
        el = self.GetCurrentElement()
        if el.type == 'Quad':
            self.GetQuadzs.SetValue("%f"%el.zs)
            self.GetQuadze.SetValue("%f"%el.ze)
            self.GetQuadde.SetValue("%e"%el.de)
            self.GetQuaddb.SetValue("%f"%el.db)
            self.Quad.Show(1)
            self.shownelem = self.Quad
        elif el.type == 'Drft':
            self.GetDrftzs.SetValue("%f"%el.zs)
            self.GetDrftze.SetValue("%f"%el.ze)
            self.GetDrftap.SetValue("%f"%el.ap)
            self.Drft.Show(1)
            self.shownelem = self.Drft

    def OnSetmadlatticeTextEnter(self, event):
        ll = self.SetMADLattice.GetValue()
        if ll:
            try:
                self.madlattice = __main__.__dict__[ll]
                if self.MADelementsunique:
                    self.madlattice.deepexpand()
                else:
                    self.madlattice.expand(lredo=1)
                self.madlattice.setextent(0.)
                self.numelements = len(self.madlattice)
                self.UpdatePanel()
            except:
                self.madlattice = None
                self.numelements = len(self.sortedelems)
                self.SetMADLattice.SetValue('')
        else:
            self.madlattice = None
            self.SetMADLattice.SetValue('')
            self.numelements = len(self.sortedelems)

    def OnMakeuniqueCheckbox(self, event):
        self.MADelementsunique = self.MakeUnique.GetValue()
            



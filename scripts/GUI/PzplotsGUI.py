#Boa:FramePanel:panel

from wx import *
from warp import *
from StringIO import *

[wxID_PANEL, wxID_PANELCOLOR, wxID_PANELLINETYPE, wxID_PANELMARKER, 
 wxID_PANELMARKERSIZE, wxID_PANELMARKS, wxID_PANELSIZE, wxID_PANELSTATICTEXT1, 
 wxID_PANELSTATICTEXT2, wxID_PANELSTATICTEXT3, 
] = map(lambda _init_ctrls: wx.NewId(), range(10))

class panel(wx.Panel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='Zplots', parent=prnt,
              pos=wx.Point(0, 0), size=wx.Size(604, 349), style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(596, 315))

        self.Color = wx.Choice(choices=['black', 'white', 'red', 'green', 'blue',
              'cyan', 'magenta', 'yellow'], id=wxID_PANELCOLOR, name='Color',
              parent=self, pos=wx.Point(512, 0), size=wx.Size(80, 21), style=0,
              validator=wx.DefaultValidator)
        EVT_CHOICE(self.Color, wxID_PANELCOLOR, self.OnColorChoice)

        self.LineType = wx.Choice(choices=['solid', 'dash', 'dot', 'dashdot',
              'dashdotdot', 'none'], id=wxID_PANELLINETYPE, name='LineType',
              parent=self, pos=wx.Point(512, 24), size=wx.Size(80, 21), style=0,
              validator=wx.DefaultValidator)
        EVT_CHOICE(self.LineType, wxID_PANELLINETYPE, self.OnLinetypeChoice)

        self.Marker = wx.TextCtrl(id=wxID_PANELMARKER, name='Marker',
              parent=self, pos=wx.Point(560, 72), size=wx.Size(22, 22),
              style=wx.TE_PROCESS_ENTER, value='A')
        EVT_TEXT_ENTER(self.Marker, wxID_PANELMARKER, self.OnMarkerTextEnter)

        self.staticText1 = wx.StaticText(id=wxID_PANELSTATICTEXT1,
              label='Marker', name='staticText1', parent=self, pos=wx.Point(516,
              74), size=wx.Size(39, 16), style=0)

        self.staticText2 = wx.StaticText(id=wxID_PANELSTATICTEXT2, label='Size',
              name='staticText2', parent=self, pos=wx.Point(516, 50),
              size=wx.Size(24, 16), style=0)

#        self.Size = wx.SpinCtrl(id=wxID_PANELSIZE, initial=1, max=10, min=1,
#              name='Size', parent=self, pos=wx.Point(550, 48), 
#              size=wx.Size(40,22), style=wx.SP_ARROW_KEYS)
#        EVT_SPINCTRL(self.Size, wxID_PANELSIZE, self.OnSizeSpinctrl)

        self.MarkerSize = wx.SpinCtrl(id=wxID_PANELMARKERSIZE, initial=1, max=10,
              min=1, name='MarkerSize', parent=self, pos=wx.Point(550, 96),
              size=wx.Size(40, 22), style=wx.SP_ARROW_KEYS)
        EVT_SPINCTRL(self.MarkerSize, wxID_PANELMARKERSIZE,
              self.OnMarkersizeSpinctrl)

        self.staticText3 = wx.StaticText(id=wxID_PANELSTATICTEXT3, label='Size',
              name='staticText3', parent=self, pos=wx.Point(516, 98),
              size=wx.Size(24, 16), style=0)

        self.Marks = wx.CheckBox(id=wxID_PANELMARKS, label='', name='Marks',
              parent=self, pos=wx.Point(496, 71), size=wx.Size(20, 20), style=0)
        self.Marks.SetValue(false)
        EVT_CHECKBOX(self.Marks, wxID_PANELMARKS, self.OnMarksCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wx.Point(0,0))
        import pzplots
        import string
        listpzplots = StringIO(pzplots.__doc__)
        self.pltcolor = 'fg'
        self.pltlinetype = 'solid'
        self.pltwidth = 1
        self.pltmarks = false
        self.pltmarker = 'A'
        self.pltmsize = 1
        self.varsuffix = None
        doread = true
        i = 0
        il = 0
        ibloc = 0
        iymin = 0
        nlong = 7
        ixsize = 60
        iysize = 20
        while(doread):
            line = listpzplots.readline()
            ipos = string.find(line,':')
            if(ipos>0):
                name = string.strip(line[:ipos])
                help = string.strip(line[ipos+1:])
                if(help==''):
                    if(il<>0):
                        iymin = iymin + il*iysize + 4
                        il = 0
                    i = 0
                    il = il+1
                    self.AddPzplotsCategory(name, il, 100, 14, iymin)
                else:
                    if(i%nlong==0):
                        il=il+1
                        i=0
                    i = i+1
                    self.AddPzplotsButton(i-1, name[2:], help, il, nlong, ixsize, iysize, iymin) 
            elif(line==''):
                doread = false
                              
    def AddPzplotsCategory(self, n, il, ixsize, iysize, iymin):
        exec("[wxID_PANEL"+n+",] = map(lambda _init_ctrls: wx.NewId(), range(1))")
        ix = 0
        iy = iymin + (il-1)*iysize 
        exec("self."+n+" = wx.StaticText(id=wxID_PANEL"+n+", label='"+n+"',name='"+n+"', parent=self, pos=wx.Point(%g, %g), size=wx.Size(%g,%g), style=0)"%(ix,iy,ixsize,iysize))

    def AddPzplotsButton(self, i, n, h, il, nlong, ixsize, iysize, iymin):
        exec("def On"+n+"Button(self,event):pz"+n+"("+
              "color=self.pltcolor,"+
              "linetype=self.pltlinetype,"+
              "width=self.pltwidth,"+
              "marks=self.pltmarks,"+
              "marker=self.pltmarker,"+
              "msize=self.pltmsize)")
        import new
        exec("self.On"+n+"Button=new.instancemethod(On"+n+"Button,self,panel)")
        exec("[wxID_PANEL"+n+",] = map(lambda _init_ctrls: wx.NewId(), range(1))")
        ix = i*ixsize
        iy = iymin + (il-1)*iysize
        exec("self."+n+" = wx.Button(id=wxID_PANEL"+n+", label='"+n+"',name='"+n+"', parent=self, pos=wx.Point(%g, %g), size=wx.Size(%g,%g), style=0)"%(ix,iy,ixsize,iysize))
        exec('self.'+n+'.SetToolTipString("'+h+'")')
        exec("EVT_BUTTON(self."+n+", wxID_PANEL"+n+", self.On"+n+"Button)")

    def addfunction(self,f1,f2):
        self.f1=f2

    def OnColorChoice(self, event):
        self.pltcolor = str(self.Color.GetStringSelection())

    def OnLinetypeChoice(self, event):
        self.pltlinetype = str(self.LineType.GetStringSelection())

    def OnMarkerTextEnter(self, event):
        self.pltmarker = self.Marker.GetValue()

    def OnSizeSpinctrl(self, event):
        self.pltsize = self.Size.GetValue()

    def OnMarkersizeSpinctrl(self, event):
        self.pltmsize = self.MarkerSize.GetValue()

    def OnMarksCheckbox(self, event):
        self.pltmarks = self.Marks.GetValue()


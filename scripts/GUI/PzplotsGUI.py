#Boa:FramePanel:Pzplots

from wxPython.wx import *
from warp import *
from StringIO import *

[wxID_PZPLOTS, wxID_PZPLOTSCOLOR, wxID_PZPLOTSLINETYPE, wxID_PZPLOTSMARKER, 
 wxID_PZPLOTSMARKERSIZE, wxID_PZPLOTSMARKS, wxID_PZPLOTSSIZE, 
 wxID_PZPLOTSSTATICTEXT1, wxID_PZPLOTSSTATICTEXT2, wxID_PZPLOTSSTATICTEXT3, 
] = map(lambda _init_ctrls: wxNewId(), range(10))

class Pzplots(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PZPLOTS, name='Zplots', parent=prnt,
              pos=wxPoint(591, 167), size=wxSize(596, 315),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(596, 315))

        self.Color = wxChoice(choices=['black', 'white', 'red', 'green', 'blue',
              'cyan', 'magenta', 'yellow'], id=wxID_PZPLOTSCOLOR, name='Color',
              parent=self, pos=wxPoint(512, 0), size=wxSize(80, 24), style=0,
              validator=wxDefaultValidator)
        EVT_CHOICE(self.Color, wxID_PZPLOTSCOLOR, self.OnColorChoice)

        self.LineType = wxChoice(choices=['solid', 'dash', 'dot', 'dashdot',
              'dashdotdot', 'none'], id=wxID_PZPLOTSLINETYPE, name='LineType',
              parent=self, pos=wxPoint(512, 24), size=wxSize(80, 24), style=0,
              validator=wxDefaultValidator)
        EVT_CHOICE(self.LineType, wxID_PZPLOTSLINETYPE, self.OnLinetypeChoice)

        self.Marker = wxTextCtrl(id=wxID_PZPLOTSMARKER, name='Marker',
              parent=self, pos=wxPoint(560, 72), size=wxSize(22, 22),
              style=wxTE_PROCESS_ENTER, value='A')
        EVT_TEXT_ENTER(self.Marker, wxID_PZPLOTSMARKER, self.OnMarkerTextEnter)

        self.staticText1 = wxStaticText(id=wxID_PZPLOTSSTATICTEXT1,
              label='Marker', name='staticText1', parent=self, pos=wxPoint(516,
              74), size=wxSize(39, 16), style=0)

        self.staticText2 = wxStaticText(id=wxID_PZPLOTSSTATICTEXT2,
              label='Size', name='staticText2', parent=self, pos=wxPoint(516,
              50), size=wxSize(24, 16), style=0)

        self.Size = wxSpinCtrl(id=wxID_PZPLOTSSIZE, initial=1, max=10, min=1,
              name='Size', parent=self, pos=wxPoint(550, 48), size=wxSize(40,
              22), style=wxSP_ARROW_KEYS)
        EVT_SPINCTRL(self.Size, wxID_PZPLOTSSIZE, self.OnSizeSpinctrl)

        self.MarkerSize = wxSpinCtrl(id=wxID_PZPLOTSMARKERSIZE, initial=1,
              max=10, min=1, name='MarkerSize', parent=self, pos=wxPoint(550,
              96), size=wxSize(40, 22), style=wxSP_ARROW_KEYS)
        EVT_SPINCTRL(self.MarkerSize, wxID_PZPLOTSMARKERSIZE,
              self.OnMarkersizeSpinctrl)

        self.staticText3 = wxStaticText(id=wxID_PZPLOTSSTATICTEXT3,
              label='Size', name='staticText3', parent=self, pos=wxPoint(516,
              98), size=wxSize(24, 16), style=0)

        self.Marks = wxCheckBox(id=wxID_PZPLOTSMARKS, label='', name='Marks',
              parent=self, pos=wxPoint(496, 71), size=wxSize(20, 20), style=0)
        self.Marks.SetValue(false)
        EVT_CHECKBOX(self.Marks, wxID_PZPLOTSMARKS, self.OnMarksCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        parent.AddPage(bSelect=true, imageId=-1, pPage=self, strText='Pzplots')
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
        exec("[wxID_PZPLOTS"+n+",] = map(lambda _init_ctrls: wxNewId(), range(1))")
        ix = 0
        iy = iymin + (il-1)*iysize 
        exec("self."+n+" = wxStaticText(id=wxID_PZPLOTS"+n+", label='"+n+"',name='"+n+"', parent=self, pos=wxPoint(%g, %g), size=wxSize(%g,%g), style=0)"%(ix,iy,ixsize,iysize))

    def AddPzplotsButton(self, i, n, h, il, nlong, ixsize, iysize, iymin):
        exec("def On"+n+"Button(self,event):pz"+n+"("+
              "color=self.pltcolor,"+
              "linetype=self.pltlinetype,"+
              "width=self.pltwidth,"+
              "marks=self.pltmarks,"+
              "marker=self.pltmarker,"+
              "msize=self.pltmsize)")
        import new
        exec("self.On"+n+"Button=new.instancemethod(On"+n+"Button,self,Pzplots)")
        exec("[wxID_PZPLOTS"+n+",] = map(lambda _init_ctrls: wxNewId(), range(1))")
        ix = i*ixsize
        iy = iymin + (il-1)*iysize
        exec("self."+n+" = wxButton(id=wxID_PZPLOTS"+n+", label='"+n+"',name='"+n+"', parent=self, pos=wxPoint(%g, %g), size=wxSize(%g,%g), style=0)"%(ix,iy,ixsize,iysize))
        exec('self.'+n+'.SetToolTipString("'+h+'")')
        exec("EVT_BUTTON(self."+n+", wxID_PZPLOTS"+n+", self.On"+n+"Button)")

    def addfunction(self,f1,f2):
        self.f1=f2

    def OnColorChoice(self, event):
        self.pltcolor = self.Color.GetStringSelection()

    def OnLinetypeChoice(self, event):
        self.pltlinetype = self.LineType.GetStringSelection()

    def OnMarkerTextEnter(self, event):
        self.pltmarker = self.Marker.GetValue()

    def OnSizeSpinctrl(self, event):
        self.pltsize = self.Size.GetValue()

    def OnMarkersizeSpinctrl(self, event):
        self.pltmsize = self.MarkerSize.GetValue()

    def OnMarksCheckbox(self, event):
        self.pltmarks = self.Marks.GetValue()


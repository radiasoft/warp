#Boa:FramePanel:panel

from wxPython.wx import *
from wxPython.lib.anchors import LayoutAnchors
from wxPython.grid import *
from warp import *

[wxID_PANEL, wxID_PANELBOLDLABEL, wxID_PANELCOLOR, wxID_PANELELEMENT_SPIN, 
 wxID_PANELFONTLABEL, wxID_PANELHEIGHTLABEL, wxID_PANELHIDE, 
 wxID_PANELITALICLABEL, wxID_PANELLABELAXIS, wxID_PANELMARKER_LETTER, 
 wxID_PANELMARKS, wxID_PANELMPHASE, wxID_PANELMSIZE, wxID_PANELMSPACE, 
 wxID_PANELSTATICBOX1, wxID_PANELSTATICBOX2, wxID_PANELSTATICTEXT1, 
 wxID_PANELSTATICTEXT10, wxID_PANELSTATICTEXT11, wxID_PANELSTATICTEXT2, 
 wxID_PANELSTATICTEXT3, wxID_PANELSTATICTEXT3, wxID_PANELSTATICTEXT4, 
 wxID_PANELSTATICTEXT5, wxID_PANELSTATICTEXT6, wxID_PANELSTATICTEXT7, 
 wxID_PANELSTATICTEXT9, wxID_PANELTYPE, wxID_PANELWIDTH_SLIDER, 
] = map(lambda _init_ctrls: wxNewId(), range(29))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='Gist', parent=prnt,
              pos=wxPoint(498, 297), size=wxSize(344, 232),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(336, 208))
        EVT_ENTER_WINDOW(self, self.OnPanelEnterWindow)

        self.staticText1 = wxStaticText(id=wxID_PANELSTATICTEXT1, label='Color',
              name='staticText1', parent=self, pos=wxPoint(16, 68),
              size=wxSize(44, 16), style=0)

        self.Color = wxChoice(choices=['black', 'red', 'green', 'blue', 'cyan',
              'magenta', 'yellow', 'white', 'bg', 'fg'], id=wxID_PANELCOLOR,
              name='Color', parent=self, pos=wxPoint(68, 64), size=wxSize(90,
              21), style=0, validator=wxDefaultValidator)
        self.Color.SetLabel('')
        self.Color.Show(True)
        EVT_CHOICE(self.Color, wxID_PANELCOLOR, self.OnColorChoice)

        self.staticText2 = wxStaticText(id=wxID_PANELSTATICTEXT2,
              label='Number', name='staticText2', parent=self, pos=wxPoint(16,
              20), size=wxSize(48, 12), style=0)

        self.element_spin = wxSpinCtrl(id=wxID_PANELELEMENT_SPIN, initial=0,
              max=100, min=1, name='element_spin', parent=self, pos=wxPoint(68,
              16), size=wxSize(40, 21), style=wxSP_ARROW_KEYS)
        EVT_SPINCTRL(self.element_spin, wxID_PANELELEMENT_SPIN,
              self.OnElementSpinctrl)

        self.StaticText3 = wxStaticText(id=wxID_PANELSTATICTEXT3, label='Width',
              name='StaticText3', parent=self, pos=wxPoint(16, 92),
              size=wxSize(44, 18), style=0)

        self.width_slider = wxSlider(id=wxID_PANELWIDTH_SLIDER, maxValue=10,
              minValue=1, name='width_slider', parent=self, point=wxPoint(68,
              92), size=wxSize(90, 21), style=wxSL_HORIZONTAL,
              validator=wxDefaultValidator, value=0)
        EVT_SCROLL(self.width_slider, self.OnWidthSliderScroll)

        self.Hide = wxCheckBox(id=wxID_PANELHIDE, label='Hide', name='Hide',
              parent=self, pos=wxPoint(112, 18), size=wxSize(50, 20), style=0)
        self.Hide.SetValue(False)
        EVT_CHECKBOX(self.Hide, wxID_PANELHIDE, self.OnHideCheckbox)

        self.staticText3 = wxStaticText(id=wxID_PANELSTATICTEXT3, label='Type',
              name='staticText3', parent=self, pos=wxPoint(16, 44),
              size=wxSize(44, 13), style=0)

        self.Type = wxChoice(choices=['solid', 'dash', 'dot', 'dashdot',
              'dashdotdot', 'none'], id=wxID_PANELTYPE, name='Type',
              parent=self, pos=wxPoint(68, 40), size=wxSize(90, 21), style=0,
              validator=wxDefaultValidator)
        EVT_CHOICE(self.Type, wxID_PANELTYPE, self.OnTypeChoice)

        self.staticText4 = wxStaticText(id=wxID_PANELSTATICTEXT4,
              label='Marker', name='staticText4', parent=self, pos=wxPoint(16,
              120), size=wxSize(44, 13), style=0)

        self.Marks = wxCheckBox(id=wxID_PANELMARKS, label='Marks', name='Marks',
              parent=self, pos=wxPoint(104, 116), size=wxSize(56, 20), style=0)
        self.Marks.SetValue(False)
        EVT_CHECKBOX(self.Marks, wxID_PANELMARKS, self.OnMarksCheckbox)

        self.marker_letter = wxTextCtrl(id=wxID_PANELMARKER_LETTER,
              name='marker_letter', parent=self, pos=wxPoint(68, 116),
              size=wxSize(24, 21), style=0, value='textCtrl1')
        EVT_TEXT(self.marker_letter, wxID_PANELMARKER_LETTER,
              self.OnMarker_letterText)

        self.staticText5 = wxStaticText(id=wxID_PANELSTATICTEXT5, label='Size',
              name='staticText5', parent=self, pos=wxPoint(16, 136),
              size=wxSize(44, 16), style=0)

        self.staticText6 = wxStaticText(id=wxID_PANELSTATICTEXT6, label='Phase',
              name='staticText6', parent=self, pos=wxPoint(16, 152),
              size=wxSize(44, 16), style=0)

        self.staticText7 = wxStaticText(id=wxID_PANELSTATICTEXT7, label='Space',
              name='staticText7', parent=self, pos=wxPoint(16, 168),
              size=wxSize(44, 16), style=0)

        self.msize = wxSlider(id=wxID_PANELMSIZE, maxValue=10, minValue=1,
              name='msize', parent=self, point=wxPoint(68, 136), size=wxSize(90,
              21), style=wxSL_HORIZONTAL, validator=wxDefaultValidator,
              value=0)
        EVT_SCROLL(self.msize, self.OnMsizeScroll)

        self.mspace = wxSlider(id=wxID_PANELMSPACE, maxValue=100, minValue=1,
              name='mspace', parent=self, point=wxPoint(68, 168),
              size=wxSize(90, 21), style=wxSL_HORIZONTAL,
              validator=wxDefaultValidator, value=0)
        EVT_SCROLL(self.mspace, self.OnMspaceScroll)

        self.staticText9 = wxStaticText(id=wxID_PANELSTATICTEXT9, label='Font',
              name='staticText9', parent=self, pos=wxPoint(184, 24),
              size=wxSize(31, 20), style=0)

        self.staticText10 = wxStaticText(id=wxID_PANELSTATICTEXT10,
              label='Height', name='staticText10', parent=self, pos=wxPoint(184,
              72), size=wxSize(41, 16), style=0)

        self.BoldLabel = wxCheckBox(id=wxID_PANELBOLDLABEL, label='Bold',
              name='BoldLabel', parent=self, pos=wxPoint(224, 48),
              size=wxSize(48, 16), style=0)
        self.BoldLabel.SetValue(False)
        EVT_CHECKBOX(self.BoldLabel, wxID_PANELBOLDLABEL,
              self.OnBoldlabelCheckbox)

        self.ItalicLabel = wxCheckBox(id=wxID_PANELITALICLABEL, label='Italic',
              name='ItalicLabel', parent=self, pos=wxPoint(272, 48),
              size=wxSize(48, 16), style=0)
        self.ItalicLabel.SetValue(False)
        EVT_CHECKBOX(self.ItalicLabel, wxID_PANELITALICLABEL,
              self.OnItaliclabelCheckbox)

        self.LabelAxis = wxChoice(choices=['x', 'y', 'all'],
              id=wxID_PANELLABELAXIS, name='LabelAxis', parent=self,
              pos=wxPoint(224, 96), size=wxSize(96, 21), style=0,
              validator=wxDefaultValidator)

        self.FontLabel = wxChoice(choices=['Courier', 'Times', 'Helvetica',
              'Symbol', 'New Century'], id=wxID_PANELFONTLABEL,
              name='FontLabel', parent=self, pos=wxPoint(224, 22),
              size=wxSize(96, 21), style=0, validator=wxDefaultValidator)
        EVT_CHOICE(self.FontLabel, wxID_PANELFONTLABEL, self.OnFontlabelChoice)

        self.HeightLabel = wxSlider(id=wxID_PANELHEIGHTLABEL, maxValue=100,
              minValue=1, name='HeightLabel', parent=self, point=wxPoint(224,
              70), size=wxSize(96, 20), style=wxSL_HORIZONTAL,
              validator=wxDefaultValidator, value=0)
        EVT_SCROLL(self.HeightLabel, self.OnHeightlabelScroll)

        self.staticText11 = wxStaticText(id=wxID_PANELSTATICTEXT11,
              label='Axis', name='staticText11', parent=self, pos=wxPoint(184,
              100), size=wxSize(30, 13), style=0)

        self.staticBox1 = wxStaticBox(id=wxID_PANELSTATICBOX1, label='Labels',
              name='staticBox1', parent=self, pos=wxPoint(176, 0),
              size=wxSize(152, 200), style=0)
        self.staticBox1.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))

        self.staticBox2 = wxStaticBox(id=wxID_PANELSTATICBOX2, label='Element',
              name='staticBox2', parent=self, pos=wxPoint(8, 0),
              size=wxSize(160, 200), style=0)
        self.staticBox2.SetFont(wxFont(10, wxSWISS, wxNORMAL, wxBOLD, False,
              'MS Sans Serif'))

        self.mphase = wxSlider(id=wxID_PANELMPHASE, maxValue=100, minValue=0,
              name='mphase', parent=self, point=wxPoint(68, 152),
              size=wxSize(90, 21), style=wxSL_HORIZONTAL,
              validator=wxDefaultValidator, value=0)
        EVT_SCROLL(self.mphase, self.OnMphaseScroll)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wxPoint(0,0))
        self.element=1
        self.element_spin.SetValue(self.element)
        self.updatenow=1
        self.setcolor=0
        self.setwidth=0
        self.settype =0
        self.setmarker=0
        self.setmarks =0
        self.setmsize =0
        self.setmphase =0
        self.setmspace =0
        self.setfontlabel   = 0
        self.setheightlabel = 0
        self.setboldlabel   = 0
        self.setitaliclabel = 0
        self.getlist()
        self.initoptions()
      
    def getlist(self):
        self.list = aplq()
        self.nelements = len(self.list)    
            
    def initoptions(self):
        self.element_spin.SetRange(1,max(1,self.nelements))
        if(self.nelements==0): return
        plist = self.list[self.element-1]
        try:
          icolor = int(plist['color'])
          if icolor>=246:
            self.Color.SetStringSelection(['bg','fg','black','white','red','green','blue','cyan','magenta','yellow'][255-icolor])
        except:
          pass
        try:
          self.width_slider.SetValue(int(plist['width']))
        except:
          pass
        try:
          self.Type.SetStringSelection(plist['type'])  
        except:
          pass
        try:
          self.Marks.SetValue(int(plist['marks']))
        except:
          pass
        try:
	  marker = plist['marker']
	  if marker=='\\':marker='.'
          self.marker_letter.SetValue(marker)
        except:
          pass
        try:
          self.msize.SetValue(int(plist['msize']))
        except:
          pass
        try:
          self.mspace.SetValue(int(plist['mspace']*100))
        except:
          pass
        try:
          self.mphase.SetValue(int(plist['mphase']*100))
        except:
          pass
        self.Hide.SetValue(int(plist['hide']))
        isys = plsys(plsys())
        font = get_style()['systems'][isys-1]['ticks']['horizontal']['textStyle']['font']
        bold = font%4%2
        italic = (font%4-bold)/2
        font = (font-2*italic-bold)/4
        font = ['Courier','Times','Helvetica','Symbol','New Century'][font]        
        self.FontLabel.SetStringSelection(font)
        self.BoldLabel.SetValue(bold)
        self.ItalicLabel.SetValue(italic)
        height = nint(get_style()['systems'][isys-1]['ticks']['horizontal']['textStyle']['height']/0.0003)
        self.HeightLabel.SetValue(height)

    def OnUpdateButton(self, event):
        if(self.nelements<1):return
        if(self.setcolor): 
	    pledit(self.element,color=str(self.Color.GetStringSelection()))
            self.setcolor=0
        if(self.setwidth): 
            pledit(self.element,width=self.width_slider.GetValue())
            self.setwidth=0  
        if(self.settype):  
            pledit(self.element,type =str(self.Type.GetStringSelection()))
            self.settype=0  
        if(self.setmarker): 
            pledit(self.element,marker=str(self.marker_letter.GetValue()))
            self.setmarker=0
        if(self.setmarks): 
            pledit(self.element,marks=self.Marks.GetValue())
            self.setmarks=0  
        if(self.setmsize): 
            pledit(self.element,msize=self.msize.GetValue())
            self.setmsize=0  
        if(self.setmphase): 
            pledit(self.element,mphase=self.mphase.GetValue()*0.01)
            self.setmphase=0  
        if(self.setmspace): 
            pledit(self.element,mspace=self.mspace.GetValue()*0.01)
            self.setmspace=0  
        if(self.setfontlabel or self.setboldlabel or self.setitaliclabel):
            set_label(font=str(self.FontLabel.GetStringSelection()), \
                      bold=self.BoldLabel.GetValue() , \
                      italic=self.ItalicLabel.GetValue(),\
                      axis=str(self.LabelAxis.GetStringSelection())) 
            self.setfontlabel=0  
            self.setboldlabel=0  
            self.setitaliclabel=0  
        if(self.setheightlabel):
            set_label(height=self.HeightLabel.GetValue()*0.0003, \
                      axis=str(self.LabelAxis.GetStringSelection())) 
            self.setheightlabel=0  
        event.Skip()

    def OnImmediateCheckbox(self, event):
        self.updatenow = self.immediate.GetValue()
        event.Skip()

    def OnColorChoice(self, event):
        self.setcolor=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnElementSpinctrl(self, event):
        self.element=self.element_spin.GetValue()
        self.initoptions()
        event.Skip()

    def OnWidthSliderScroll(self, event):
        self.setwidth=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnPanelEnterWindow(self, event):
        self.getlist()
        self.initoptions()
        event.Skip()

    def OnHideCheckbox(self, event):
        if(self.nelements<1):return
        pledit(self.element,hide=self.Hide.GetValue())
        event.Skip()

    def OnTypeChoice(self, event):
        self.settype=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnMarksCheckbox(self, event):
        self.setmarks=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnMarker_letterText(self, event):
        self.setmarker=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnMsizeScroll(self, event):
        self.setmsize=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnMphaseScroll(self, event):
        self.setmphase=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnMspaceScroll(self, event):
        self.setmspace=1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnBoldlabelCheckbox(self, event):
        self.setboldlabel = 1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnItaliclabelCheckbox(self, event):
        self.setitaliclabel = 1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnFontlabelChoice(self, event):
        self.setfontlabel = 1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()

    def OnHeightlabelScroll(self, event):
        self.setheightlabel = 1
        if self.updatenow: self.OnUpdateButton(event)
        event.Skip()
       
            

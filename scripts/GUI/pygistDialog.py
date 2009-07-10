#Boa:FramePanel:panel

from wx import *
from wx.lib.anchors import LayoutAnchors
from wx.grid import *
from warp import *

[wxID_PANEL, wxID_PANELBOLDLABEL, wxID_PANELCOLOR, wxID_PANELELEMENT_SPIN, 
 wxID_PANELFILENAME, wxID_PANELFONTLABEL, wxID_PANELFORMAT, 
 wxID_PANELHEIGHTLABEL, wxID_PANELHIDE, wxID_PANELITALICLABEL, 
 wxID_PANELLABELAXIS, wxID_PANELMARKER_LETTER, wxID_PANELMARKS, 
 wxID_PANELMPHASE, wxID_PANELMSIZE, wxID_PANELMSPACE, wxID_PANELSAVE, 
 wxID_PANELSTATICBOX1, wxID_PANELSTATICBOX2, wxID_PANELSTATICBOX3, 
 wxID_PANELSTATICTEXT1, wxID_PANELSTATICTEXT10, wxID_PANELSTATICTEXT11, 
 wxID_PANELSTATICTEXT12, wxID_PANELSTATICTEXT2, wxID_PANELSTATICTEXT3, 
 wxID_PANELSTATICTEXT3, wxID_PANELSTATICTEXT4, wxID_PANELSTATICTEXT5, 
 wxID_PANELSTATICTEXT6, wxID_PANELSTATICTEXT7, wxID_PANELSTATICTEXT8, 
 wxID_PANELSTATICTEXT9, wxID_PANELTYPE, wxID_PANELWIDTH_SLIDER, 
] = map(lambda _init_ctrls: wx.NewId(), range(35))

class panel(wx.Panel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='Gist', parent=prnt,
              pos=wx.Point(415, 280), size=wx.Size(342, 232),
              style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(334, 208))
        EVT_ENTER_WINDOW(self, self.OnPanelEnterWindow)

        self.staticText1 = wx.StaticText(id=wxID_PANELSTATICTEXT1, label='Color',
              name='staticText1', parent=self, pos=wx.Point(16, 68),
              size=wx.Size(44, 16), style=0)

        self.Color = wx.Choice(choices=['black', 'red', 'green', 'blue', 'cyan',
              'magenta', 'yellow', 'white', 'bg', 'fg'], id=wxID_PANELCOLOR,
              name='Color', parent=self, pos=wx.Point(68, 64), size=wx.Size(90,
              21), style=0, validator=wx.DefaultValidator)
        self.Color.SetLabel('')
        self.Color.Show(True)
        EVT_CHOICE(self.Color, wxID_PANELCOLOR, self.OnColorChoice)

        self.staticText2 = wx.StaticText(id=wxID_PANELSTATICTEXT2,
              label='Number', name='staticText2', parent=self, pos=wx.Point(16,
              20), size=wx.Size(48, 12), style=0)

        self.element_spin = wx.SpinCtrl(id=wxID_PANELELEMENT_SPIN, initial=0,
              max=100, min=1, name='element_spin', parent=self, pos=wx.Point(68,
              16), size=wx.Size(40, 21), style=wx.SP_ARROW_KEYS)
        EVT_SPINCTRL(self.element_spin, wxID_PANELELEMENT_SPIN,
              self.OnElementSpinctrl)

        self.StaticText3 = wx.StaticText(id=wxID_PANELSTATICTEXT3, label='Width',
              name='StaticText3', parent=self, pos=wx.Point(16, 92),
              size=wx.Size(44, 18), style=0)

        self.width_slider = wx.Slider(id=wxID_PANELWIDTH_SLIDER, maxValue=10,
              minValue=1, name='width_slider', parent=self, point=wx.Point(68,
              92), size=wx.Size(90, 21), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=0)
        EVT_SCROLL(self.width_slider, self.OnWidthSliderScroll)

        self.Hide = wx.CheckBox(id=wxID_PANELHIDE, label='Hide', name='Hide',
              parent=self, pos=wx.Point(112, 18), size=wx.Size(50, 20), style=0)
        self.Hide.SetValue(False)
        EVT_CHECKBOX(self.Hide, wxID_PANELHIDE, self.OnHideCheckbox)

        self.staticText3 = wx.StaticText(id=wxID_PANELSTATICTEXT3, label='Type',
              name='staticText3', parent=self, pos=wx.Point(16, 44),
              size=wx.Size(44, 13), style=0)

        self.Type = wx.Choice(choices=['solid', 'dash', 'dot', 'dashdot',
              'dashdotdot', 'none'], id=wxID_PANELTYPE, name='Type',
              parent=self, pos=wx.Point(68, 40), size=wx.Size(90, 21), style=0,
              validator=wx.DefaultValidator)
        EVT_CHOICE(self.Type, wxID_PANELTYPE, self.OnTypeChoice)

        self.staticText4 = wx.StaticText(id=wxID_PANELSTATICTEXT4,
              label='Marker', name='staticText4', parent=self, pos=wx.Point(16,
              120), size=wx.Size(44, 13), style=0)

        self.Marks = wx.CheckBox(id=wxID_PANELMARKS, label='Marks', name='Marks',
              parent=self, pos=wx.Point(104, 116), size=wx.Size(56, 20), style=0)
        self.Marks.SetValue(False)
        EVT_CHECKBOX(self.Marks, wxID_PANELMARKS, self.OnMarksCheckbox)

        self.marker_letter = wx.TextCtrl(id=wxID_PANELMARKER_LETTER,
              name='marker_letter', parent=self, pos=wx.Point(68, 116),
              size=wx.Size(24, 21), style=0, value='textCtrl1')
        EVT_TEXT(self.marker_letter, wxID_PANELMARKER_LETTER,
              self.OnMarker_letterText)

        self.staticText5 = wx.StaticText(id=wxID_PANELSTATICTEXT5, label='Size',
              name='staticText5', parent=self, pos=wx.Point(16, 136),
              size=wx.Size(44, 16), style=0)

        self.staticText6 = wx.StaticText(id=wxID_PANELSTATICTEXT6, label='Phase',
              name='staticText6', parent=self, pos=wx.Point(16, 152),
              size=wx.Size(44, 16), style=0)

        self.staticText7 = wx.StaticText(id=wxID_PANELSTATICTEXT7, label='Space',
              name='staticText7', parent=self, pos=wx.Point(16, 168),
              size=wx.Size(44, 16), style=0)

        self.msize = wx.Slider(id=wxID_PANELMSIZE, maxValue=10, minValue=1,
              name='msize', parent=self, point=wx.Point(68, 136), size=wx.Size(90,
              21), style=wx.SL_HORIZONTAL, validator=wx.DefaultValidator,
              value=0)
        EVT_SCROLL(self.msize, self.OnMsizeScroll)

        self.mspace = wx.Slider(id=wxID_PANELMSPACE, maxValue=100, minValue=1,
              name='mspace', parent=self, point=wx.Point(68, 168),
              size=wx.Size(90, 21), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=0)
        EVT_SCROLL(self.mspace, self.OnMspaceScroll)

        self.staticText9 = wx.StaticText(id=wxID_PANELSTATICTEXT9, label='Font',
              name='staticText9', parent=self, pos=wx.Point(184, 20),
              size=wx.Size(31, 20), style=0)

        self.staticText10 = wx.StaticText(id=wxID_PANELSTATICTEXT10,
              label='Height', name='staticText10', parent=self, pos=wx.Point(184,
              62), size=wx.Size(41, 16), style=0)

        self.BoldLabel = wx.CheckBox(id=wxID_PANELBOLDLABEL, label='Bold',
              name='BoldLabel', parent=self, pos=wx.Point(224, 40),
              size=wx.Size(48, 16), style=0)
        self.BoldLabel.SetValue(False)
        self.BoldLabel.SetFont(wx.Font(8, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))
        EVT_CHECKBOX(self.BoldLabel, wxID_PANELBOLDLABEL,
              self.OnBoldlabelCheckbox)

        self.ItalicLabel = wx.CheckBox(id=wxID_PANELITALICLABEL, label='Italic',
              name='ItalicLabel', parent=self, pos=wx.Point(272, 40),
              size=wx.Size(48, 16), style=0)
        self.ItalicLabel.SetValue(False)
        self.ItalicLabel.SetFont(wx.Font(8, wx.SWISS, wx.ITALIC, wx.NORMAL, False,
              'MS Sans Serif'))
        EVT_CHECKBOX(self.ItalicLabel, wxID_PANELITALICLABEL,
              self.OnItaliclabelCheckbox)

        self.LabelAxis = wx.Choice(choices=['x', 'y', 'all'],
              id=wxID_PANELLABELAXIS, name='LabelAxis', parent=self,
              pos=wx.Point(224, 84), size=wx.Size(96, 21), style=0,
              validator=wx.DefaultValidator)

        self.FontLabel = wx.Choice(choices=['Courier', 'Times', 'Helvetica',
              'Symbol', 'New Century'], id=wxID_PANELFONTLABEL,
              name='FontLabel', parent=self, pos=wx.Point(224, 16),
              size=wx.Size(96, 21), style=0, validator=wx.DefaultValidator)
        EVT_CHOICE(self.FontLabel, wxID_PANELFONTLABEL, self.OnFontlabelChoice)

        self.HeightLabel = wx.Slider(id=wxID_PANELHEIGHTLABEL, maxValue=100,
              minValue=1, name='HeightLabel', parent=self, point=wx.Point(224,
              60), size=wx.Size(96, 20), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=0)
        EVT_SCROLL(self.HeightLabel, self.OnHeightlabelScroll)

        self.staticText11 = wx.StaticText(id=wxID_PANELSTATICTEXT11,
              label='Axis', name='staticText11', parent=self, pos=wx.Point(184,
              88), size=wx.Size(30, 13), style=0)

        self.staticBox1 = wx.StaticBox(id=wxID_PANELSTATICBOX1, label='Labels',
              name='staticBox1', parent=self, pos=wx.Point(176, 0),
              size=wx.Size(152, 112), style=0)
        self.staticBox1.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))

        self.staticBox2 = wx.StaticBox(id=wxID_PANELSTATICBOX2, label='Element',
              name='staticBox2', parent=self, pos=wx.Point(8, 0),
              size=wx.Size(160, 200), style=0)
        self.staticBox2.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))

        self.mphase = wx.Slider(id=wxID_PANELMPHASE, maxValue=100, minValue=0,
              name='mphase', parent=self, point=wx.Point(68, 152),
              size=wx.Size(90, 21), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=0)
        EVT_SCROLL(self.mphase, self.OnMphaseScroll)

        self.staticBox3 = wx.StaticBox(id=wxID_PANELSTATICBOX3, label='Output',
              name='staticBox3', parent=self, pos=wx.Point(176, 112),
              size=wx.Size(152, 88), style=0)
        self.staticBox3.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD, False,
              'MS Sans Serif'))

        self.staticText8 = wx.StaticText(id=wxID_PANELSTATICTEXT8, label='File',
              name='staticText8', parent=self, pos=wx.Point(184, 130),
              size=wx.Size(32, 13), style=0)

        self.staticText12 = wx.StaticText(id=wxID_PANELSTATICTEXT12,
              label='Format', name='staticText12', parent=self, pos=wx.Point(184,
              152), size=wx.Size(40, 13), style=0)

        self.Filename = wx.TextCtrl(id=wxID_PANELFILENAME, name='Filename',
              parent=self, pos=wx.Point(232, 128), size=wx.Size(88, 21), style=0,
              value='GistPlot')

        self.Format = wx.Choice(choices=['ps', 'cgm'], id=wxID_PANELFORMAT,
              name='Format', parent=self, pos=wx.Point(232, 152), size=wx.Size(88,
              21), style=0, validator=wx.DefaultValidator)

        self.Save = wx.Button(id=wxID_PANELSAVE, label='Save', name='Save',
              parent=self, pos=wx.Point(232, 176), size=wx.Size(88, 20), style=0)
        EVT_BUTTON(self.Save, wxID_PANELSAVE, self.OnSaveButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wx.Point(0,0))
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

    def OnSaveButton(self, event):
        self.filename = str(self.Filename.GetValue())
        self.format = str(self.Format.GetStringSelection())
        if self.format=='ps':
            suffix = '.ps'
        else:
            suffix = '.cgm'
        hcp_file(self.filename+suffix)
        hcp()
        hcp_finish()
        event.Skip()
       
            

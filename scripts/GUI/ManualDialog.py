#Boa:FramePanel:panel

from wx import *
from wx.html import *
from wx.lib.anchors import LayoutAnchors
from wx.grid import *
import os, sys, string
import warp

[wxID_PANEL, wxID_PANELBACK, wxID_PANELFORWARD, wxID_PANELHOME, 
 wxID_PANELHTMLWINDOW1, 
] = map(lambda _init_ctrls: wx.NewId(), range(5))

class panel(wx.Panel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='Manual', parent=prnt,
              pos=wx.Point(393, 242), size=wx.Size(513, 363),
              style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(505, 339))
        self.SetAutoLayout(True)

        self.htmlWindow1 = wx.HtmlWindow(id=wxID_PANELHTMLWINDOW1,
              name='htmlWindow1', parent=self, pos=wx.Point(0, 24),
              size=wx.Size(504, 312))

        self.Back = wx.Button(id=wxID_PANELBACK, label='Back', name='Back',
              parent=self, pos=wx.Point(0, 0), size=wx.Size(75, 23), style=0)
        EVT_BUTTON(self.Back, wxID_PANELBACK, self.OnBackButton)

        self.Forward = wx.Button(id=wxID_PANELFORWARD, label='Forward',
              name='Forward', parent=self, pos=wx.Point(80, 0), size=wx.Size(75,
              23), style=0)
        EVT_BUTTON(self.Forward, wxID_PANELFORWARD, self.OnForwardButton)

        self.Home = wx.Button(id=wxID_PANELHOME, label='Home', name='Home',
              parent=self, pos=wx.Point(160, 0), size=wx.Size(75, 23), style=0)
        EVT_BUTTON(self.Home, wxID_PANELHOME, self.OnHomeButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.html=self.htmlWindow1
	parent.GetParent().html=self.html
	self.html.GoHome=self.GoHome
        self.Move(wx.Point(0,0))
        cs = wx.LayoutConstraints()
        ref=self.GetParent()
        cs.top.SameAs(ref,wx.Top)
        cs.bottom.SameAs(ref,wx.Bottom)
        cs.left.SameAs(ref,wx.Left)
        cs.right.SameAs(ref,wx.Right)
        self.SetConstraints(cs)
        self.SetAutoLayout(True)
        cs = wx.LayoutConstraints()
        ref=self
        cs.top.SameAs(ref,wx.Top,25)
        cs.bottom.SameAs(ref,wx.Bottom)
        cs.left.SameAs(ref,wx.Left)
        cs.right.SameAs(ref,wx.Right)
        self.html.SetConstraints(cs)
        self.html.SetAutoLayout(True)
    
    def GoHome(self,which=None):
	if which is not None:self.which = which
        warp_path = os.path.dirname(warp.__file__)
        if sys.platform=='cygwin':
          cpos = string.find(warp_path,'cygdrive')
          if cpos>=0:
            warp_path = string.upper(warp_path[cpos+9])+':'+warp_path[cpos+10:]
          else:
            warp_path = warp_path[1]+':'+warp_path[3:]
        if warp_path <> '':warp_path+='/'
        self.html.LoadPage(warp_path+'doc/html/'+self.which+'.html')

    def OnBackButton(self, event):
        if not self.html.HistoryBack():
            wx.MessageBox("No more items in history!")
        event.Skip()

    def OnForwardButton(self, event):
        if not self.html.HistoryForward():
            wx.MessageBox("No more items in history!")
        event.Skip()

    def OnHomeButton(self, event):
        self.GoHome()
        event.Skip()
 

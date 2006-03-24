#Boa:FramePanel:panel

from wxPython.wx import *
from wxPython.html import *
from wxPython.lib.anchors import LayoutAnchors
from wxPython.grid import *
import os, sys, string
import warp

[wxID_PANEL, wxID_PANELBACK, wxID_PANELFORWARD, wxID_PANELHOME, 
 wxID_PANELHTMLWINDOW1, 
] = map(lambda _init_ctrls: wxNewId(), range(5))

class panel(wxPanel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wxPanel.__init__(self, id=wxID_PANEL, name='Manual', parent=prnt,
              pos=wxPoint(393, 242), size=wxSize(513, 363),
              style=wxTAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wxSize(505, 339))
        self.SetAutoLayout(True)

        self.htmlWindow1 = wxHtmlWindow(id=wxID_PANELHTMLWINDOW1,
              name='htmlWindow1', parent=self, pos=wxPoint(0, 24),
              size=wxSize(504, 312))

        self.Back = wxButton(id=wxID_PANELBACK, label='Back', name='Back',
              parent=self, pos=wxPoint(0, 0), size=wxSize(75, 23), style=0)
        EVT_BUTTON(self.Back, wxID_PANELBACK, self.OnBackButton)

        self.Forward = wxButton(id=wxID_PANELFORWARD, label='Forward',
              name='Forward', parent=self, pos=wxPoint(80, 0), size=wxSize(75,
              23), style=0)
        EVT_BUTTON(self.Forward, wxID_PANELFORWARD, self.OnForwardButton)

        self.Home = wxButton(id=wxID_PANELHOME, label='Home', name='Home',
              parent=self, pos=wxPoint(160, 0), size=wxSize(75, 23), style=0)
        EVT_BUTTON(self.Home, wxID_PANELHOME, self.OnHomeButton)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.html=self.htmlWindow1
	parent.GetParent().html=self.html
	self.html.GoHome=self.GoHome
        self.Move(wxPoint(0,0))
        cs = wxLayoutConstraints()
        ref=self.GetParent()
        cs.top.SameAs(ref,wxTop)
        cs.bottom.SameAs(ref,wxBottom)
        cs.left.SameAs(ref,wxLeft)
        cs.right.SameAs(ref,wxRight)
        self.SetConstraints(cs)
        self.SetAutoLayout(True)
        cs = wxLayoutConstraints()
        ref=self
        cs.top.SameAs(ref,wxTop,25)
        cs.bottom.SameAs(ref,wxBottom)
        cs.left.SameAs(ref,wxLeft)
        cs.right.SameAs(ref,wxRight)
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
            wxMessageBox("No more items in history!")
        event.Skip()

    def OnForwardButton(self, event):
        if not self.html.HistoryForward():
            wxMessageBox("No more items in history!")
        event.Skip()

    def OnHomeButton(self, event):
        self.GoHome()
        event.Skip()
 

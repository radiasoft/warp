#Boa:Dialog:wxDialog1

from wx import *
import WarpPanel

def create(parent):
    return wx.Dialog1(parent)

[wxID_WXDIALOG1, wxID_WXDIALOG1TONOTEBOOK, wxID_WXDIALOG1WINDOW1,
] = map(lambda _init_ctrls: wx.NewId(), range(3))

class wxDialog1(wx.Dialog):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Dialog.__init__(self, id=wxID_WXDIALOG1, name='', parent=prnt,
              pos=wx.Point(419, 165), size=wx.Size(638, 434),
              style=wx.DEFAULT_DIALOG_STYLE, title='wx.Dialog1')
        self._init_utils()
        self.SetClientSize(wx.Size(630, 410))
        EVT_CLOSE(self, self.OnWxdialog1Close)

        self.window1 = wx.Window(id=wxID_WXDIALOG1WINDOW1, name='window1',
              parent=self, pos=wx.Point(0, 24), size=wx.Size(640, 416), style=0)

        self.tonotebook = wx.Button(id=wxID_WXDIALOG1TONOTEBOOK,
              label='to notebook', name='tonotebook', parent=self,
              pos=wx.Point(0, 0), size=wx.Size(75, 23), style=0)
        self.tonotebook.SetBackgroundColour(wx.Colour(128, 128, 128))
        self.tonotebook.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.tonotebook, wxID_WXDIALOG1TONOTEBOOK,
              self.OnTonotebookButton)

    def __init__(self, parent, child, title):
        self._init_ctrls(parent)
        self.SetTitle(title)
        self.parent = parent
        self.title = title
        if child is None: return
        self.panel = child(self.window1)
        self.child = child

    def OnTonotebookButton(self, event):
        self.panel.Reparent(self.parent.notebook1.GetPage(self.nbselection))
        self.parent.notebook1.GetPage(self.nbselection).Show(1)
        self.parent.notebook1.SetSelection(self.nbselection)
        self.panel.Move(wx.Point(0,0))
#        self.panel.SetSize(self.parent.notebook1.GetSize())
        self.panel.Refresh()
        self.Destroy()
        event.Skip()

    def OnWxdialog1Close(self, event):
        self.OnTonotebookButton(event)
        event.Skip()


class wxFrame1(wx.Frame):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Frame.__init__(self, id=wxID_WXDIALOG1, name='', parent=prnt,
              pos=wx.Point(419, 165), size=wx.Size(638, 434),
              style=wx.DEFAULT_FRAME_STYLE, title='wx.Dialog1')
        self._init_utils()
        self.SetClientSize(wx.Size(630, 410))
        self.SetAutoLayout(True)
        EVT_CLOSE(self, self.OnWxframe1Close)

        self.window1 = wx.Window(id=wxID_WXDIALOG1WINDOW1, name='window1',
              parent=self, pos=wx.Point(0, 24), size=wx.Size(640, 416), style=0)

        self.tonotebook = wx.Button(id=wxID_WXDIALOG1TONOTEBOOK,
              label='to notebook', name='tonotebook', parent=self,
              pos=wx.Point(0, 0), size=wx.Size(75, 23), style=0)
        self.tonotebook.SetBackgroundColour(wx.Colour(128, 128, 128))
        self.tonotebook.SetForegroundColour(wx.Colour(255, 255, 255))
        EVT_BUTTON(self.tonotebook, wxID_WXDIALOG1TONOTEBOOK,
              self.OnTonotebookButton)

    def __init__(self, parent, child, title):
        self._init_ctrls(parent)
        self.SetTitle(title)
        self.parent = parent
        self.title = title
        if child is None: return
        self.panel = child(self.window1)
        self.child = child

    def OnTonotebookButton(self, event):
        self.panel.Reparent(self.parent.notebook1.GetPage(self.nbselection))
        self.parent.notebook1.GetPage(self.nbselection).Show(1)
        self.parent.notebook1.SetSelection(self.nbselection)
        self.panel.Move(wx.Point(0,0))
        cs = wx.LayoutConstraints()
        ref=self.parent.notebook1.GetPage(self.nbselection)
        cs.top.SameAs(ref,wx.Top)
        cs.bottom.SameAs(ref,wx.Bottom)
        cs.left.SameAs(ref,wx.Left)
        cs.right.SameAs(ref,wx.Right)
        self.panel.SetConstraints(cs)
        self.panel.SetAutoLayout(True)
        self.panel.SetSize(self.parent.notebook1.GetSize())
        self.parent.Refresh()
        self.panel.Refresh()
        self.parent.notebook1.Update()
        self.Destroy()
        event.Skip()

    def OnWxframe1Close(self, event):
        self.OnTonotebookButton(event)
        event.Skip()

#!/usr/bin/env python

#-------------- User changable settings are in configuration.py --------------

#------------------------------ System Imports -------------------------------
from __future__ import generators

import sys
#handy parsing and external command execution without external scripts
if __name__ == '__main__':
    if (len(sys.argv) > 1):
        if (sys.argv[1] == '--exec'):
            del sys.argv[1]
            import os
            cwd = os.getcwd()
            if sys.argv[1] != '*':
                os.chdir(sys.argv[1])
            args = []
            ins = 0
            for arg in sys.argv[2:]:
                if ' ' in arg:
                    ins = 1
                    if sys.platform=='win32':
                        args.append('"%s"'%arg)
                    else:
                        args.append(arg.replace(' ', '\\ '))
                else:
                    args.append(arg)
            if sys.platform=='win32':
                if ins:
                    args[0:0] = ['""']
                os.system("start %s"%' '.join(args))
                #args[0:0] = ['start']
                #os.spawnv(os.P_NOWAIT, args[0], args)
            else:
                os.spawnvp(os.P_NOWAIT, args[0], args)
            sys.exit(0)

import os
import keyword, traceback, cStringIO, imp, fnmatch, re
import time, pprint
from wxPython.wx import *
from wxPython.stc import *
from wxPython.lib.rcsizer import RowColSizer
from wxPython.lib.dialogs import wxScrolledMessageDialog
from wxPython.lib.mixins.listctrl import wxListCtrlAutoWidthMixin

#--------------------------- configuration import ----------------------------
from configuration import *
#--------- The two most useful links for constructing this editor... ---------
# http://www.pyframe.com/stc/index.html
# http://personalpages.tds.net/~edream/out2.htm

#---------------------------- Event Declarations -----------------------------
if 1:
    #under an if so that I can collapse the declarations

    VERSION = "1.9.1"
    VREQ = '2.4.2.4'

    try:
        True
    except:
        True = bool(1)
        False = not True

    import string
    STRINGPRINTABLE = string.printable[:]
    STRINGPRINTABLE = dict(zip(map(ord, STRINGPRINTABLE), len(STRINGPRINTABLE)*[None]))
    ST = str(long(time.time()*100))
    OUTF = "%s/.%s.tmp"%(homedir, ST)
    INF = OUTF+".out"

    class cancelled(Exception):
        pass

    #keypresses
    if 1:
        keys = ["BACK", "TAB", "RETURN", "ESCAPE", "SPACE", "DELETE", "START", 
        "LBUTTON", "RBUTTON", "CANCEL", "MBUTTON", "CLEAR", "PAUSE", 
        "CAPITAL", "PRIOR", "NEXT", "END", "HOME", "LEFT", "UP", "RIGHT", 
        "DOWN", "SELECT", "PRINT", "EXECUTE", "SNAPSHOT", "INSERT", "HELP", 
        "NUMPAD0", "NUMPAD1", "NUMPAD2", "NUMPAD3", "NUMPAD4", "NUMPAD5", 
        "NUMPAD6", "NUMPAD7", "NUMPAD8", "NUMPAD9", "MULTIPLY", "ADD", 
        "SEPARATOR", "SUBTRACT", "DECIMAL", "DIVIDE", "F1", "F2", "F3", "F4", 
        "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14", 
        "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23", "F24", 
        "NUMLOCK", "SCROLL", "PAGEUP", "PAGEDOWN", "NUMPAD_SPACE", 
        "NUMPAD_TAB", "NUMPAD_ENTER", "NUMPAD_F1", "NUMPAD_F2", "NUMPAD_F3", 
        "NUMPAD_F4", "NUMPAD_HOME", "NUMPAD_LEFT", "NUMPAD_UP", 
        "NUMPAD_RIGHT", "NUMPAD_DOWN", "NUMPAD_PRIOR", "NUMPAD_PAGEUP", 
        "NUMPAD_NEXT", "NUMPAD_PAGEDOWN", "NUMPAD_END", "NUMPAD_BEGIN", 
        "NUMPAD_INSERT", "NUMPAD_DELETE", "NUMPAD_EQUAL", "NUMPAD_MULTIPLY", 
        "NUMPAD_ADD", "NUMPAD_SEPARATOR", "NUMPAD_SUBTRACT", "NUMPAD_DECIMAL", 
        "NUMPAD_DIVIDE"] 
        
        keyMap = {}
        #RkeyMap = {}
        for i in keys:
            key = eval("WXK_"+i)
            keyMap[key] = i
            #RkeyMap[i] = key
        for i in ["SHIFT", "ALT", "CONTROL", "MENU"]:
            key = eval("WXK_"+i)
            keyMap[key] = ''
        del key
        
        def GetKeyPress(evt):
            keycode = evt.GetKeyCode()
            keyname = keyMap.get(keycode, None)
            modifiers = ""
            for mod, ch in [(evt.ControlDown(), 'Ctrl+'),
                            (evt.AltDown(),     'Alt+'),
                            (evt.ShiftDown(),   'Shift+'),
                            (evt.MetaDown(),    'Meta+')]:
                if mod:
                    modifiers += ch
            if keyname is None:
                if 27 < keycode < 256:
                    keyname = chr(keycode)
                else:
                    keyname = "(%s)unknown" % keycode
            return modifiers + keyname

    MENULIST = []
    MENUPREF = {}
    OLD_MENUPREF= {}
    HOTKEY_TO_ID = {}

    def _spl(n):
        return (n.split('\t', 1) + [''])[:2]

    def recmenu(menu, id):
        #this was a pain in the ass.
        if isinstance(menu, wxMenuItem) or isinstance(menu, wxMenuItemPtr):
            sm = menu.GetSubMenu()
            if sm:
                a = recmenu(sm, id)
                if a:
                    return '%s->%s'%(menu.GetLabel(), a)
                return ''
            elif menu.GetId() == id:
                return menu.GetLabel()
            else:
                return ''
        elif isinstance(menu, wxMenu) or isinstance(menu, wxMenuPtr):
            ITEMS = menu.GetMenuItems()
            for i in ITEMS:
                a = recmenu(i, id)
                if a:
                    return a
            return ''
        else:
            for i in xrange(menu.GetMenuCount()):
                r = menu.GetMenu(i)
                if r.FindItemById(id):
                    return "%s->%s"%(menu.GetLabelTop(i), recmenu(r, id))
        raise Exception("Item not found.")

    def menuAdd(root, menu, name, desc, funct, id, kind=wxITEM_NORMAL):
        
        a = wxMenuItem(menu, id, 'TEMPORARYNAME', desc, kind)
        menu.AppendItem(a)
        EVT_MENU(root, id, funct)

        ns, oacc = _spl(name)
        heir = recmenu(menuBar, id)[:-13] + ns
        if heir in MENUPREF:
            name, acc = MENUPREF[heir]
        else:
            if heir in OLD_MENUPREF:
                name, acc = MENUPREF[heir] = OLD_MENUPREF[heir]
            else:
                name, acc = MENUPREF[heir] = (ns, oacc)
            MENULIST.append((heir, name, oacc, acc, kind in [wxITEM_NORMAL, wxITEM_CHECK]))

        if acc:
            HOTKEY_TO_ID[acc] = id

        menuBar.SetLabel(id, '%s\t%s'%(name, acc))

    def menuAddM(parent, menu, name, help=''):
        if isinstance(parent, wxMenu) or isinstance(parent, wxMenuPtr):
            id = wxNewId()
            parent.AppendMenu(id, "TEMPORARYNAME", menu, help)
            heir = recmenu(menuBar, id) + name
            name, toss = MENUPREF.setdefault(heir, (name, ''))

            menuBar.SetLabel(id, name)
        else:
            heir = name
            name, toss = MENUPREF.setdefault(name, (name, ''))
            
            parent.Append(menu, name)
        
        MENULIST.append((heir, name, '', '', 0))

    def getIcon():
        data = getData()
        stream = cStringIO.StringIO(data)
        image = wxImageFromStream(stream)
        bitmap = wxBitmapFromImage(image)
        icon = wxEmptyIcon()
        icon.CopyFromBitmap(bitmap)
        return icon

    NEWDOCUMENT = 0L
    
    #required ids
    if 1:
        #style ids
        PY_S = wxNewId()
        HT_S = wxNewId()
        CC_S = wxNewId()
        XM_S = wxNewId()
        TX_S = wxNewId()
        lexers = dict(zip([PY_S, HT_S, CC_S, XM_S, TX_S], ['python', 'html', 'cpp', 'xml', 'text']))
        lexers2 = dict(zip([wxSTC_LEX_PYTHON, wxSTC_LEX_HTML, wxSTC_LEX_CPP, wxSTC_LEX_XML, wxSTC_LEX_NULL], [PY_S, HT_S, CC_S, XM_S, TX_S]))

        PY_DS = wxNewId()
        HT_DS = wxNewId()
        CC_DS = wxNewId()
        XM_DS = wxNewId()
        TX_DS = wxNewId()
        lexers.update(dict(zip([PY_DS, HT_DS, CC_DS, XM_DS, TX_DS], ['python', 'html', 'cpp', 'xml', 'text'])))
        lexers3 = dict(zip(['python', 'html', 'cpp', 'xml', 'text'], [PY_DS, HT_DS, CC_DS, XM_DS, TX_DS]))
    
        #checkbox ids
        SNIPT = wxNewId()
        AUTO = wxNewId()
        NUM = wxNewId()
        MARGIN = wxNewId()
        USETABS = wxNewId()
        INDENTGUIDE = wxNewId()
        WRAPL = wxNewId()
        SORTBY = wxNewId()
        TODOT = wxNewId()
        SAVE_CURSOR = wxNewId()

        ZI = wxNewId()
        
        #margin ids
        LL_BACK = wxNewId()# wxSTC_EDGE_BACKGROUND
        LL_LINE = wxNewId()# wxSTC_EDGE_LINE
        LL_NONE = wxNewId()# wxSTC_EDGE_NONE
        LL_MAPPING = {LL_BACK:wxSTC_EDGE_BACKGROUND,
                      LL_LINE:wxSTC_EDGE_LINE,
                      LL_NONE:wxSTC_EDGE_NONE}
        LL_RMAPPING = {}
        for i,j in LL_MAPPING.iteritems():
            LL_RMAPPING[j] = i

        #line ending ids
        LE_CRLF = wxNewId()
        LE_LF   = wxNewId()
        LE_CR   = wxNewId()
        LE_MAPPING = {LE_CRLF:wxSTC_EOL_CRLF,
                        LE_LF:wxSTC_EOL_LF,
                        LE_CR:wxSTC_EOL_CR}
        LE_RMAPPING = {}
        for i,j in LE_MAPPING.iteritems():
            LE_RMAPPING[j] = i


    
        #pathmark ids
        BM1 = wxNewId()
        BM2 = wxNewId()
        BM3 = wxNewId()
        BM4 = wxNewId()
        BM5 = wxNewId()
        BM6 = wxNewId()
        BM7 = wxNewId()
        BM8 = wxNewId()
        BM9 = wxNewId()
    
        pm = [BM1,BM2,BM3,BM4,BM5,BM6,BM7,BM8,BM9]
        kp = range(49,58)
        #The range(49,58) represents the key codes for numbers 1...9 inclusve.
        #                                key codes            49...57
        pathmarks = dict(zip(kp, 9*[0]))
        bmId2Keypress = dict(zip(pm, kp))
        bmPm2Id = dict(zip(kp, pm))
        del pm;del kp

    #bookmark support
    BOOKMARKNUMBER = 1
    BOOKMARKSYMBOL = wxSTC_MARK_CIRCLE
    BOOKMARKMASK = 2

    #unicode BOM stuff
    if 1:
        BOM = [('+/v8-', 'utf-7'),
               ('\xef\xbb\xbf', 'utf-8'),
               ('\xfe\xff', 'utf-16-be'),
               ('\xff\xfe\xff\xfe', 'utf-16'),
               ('\xff\xfe', 'utf-16-le'),
               ('', 'ascii')]
        ADDBOM = {}
        ENCODINGS = {}
        for i,j in BOM:
            ADDBOM[j] = i
            ENCODINGS[j] = wxNewId()
        del i;del j;

    #font stuff
    if 1:
        cn = 'Courier New'
        if wxPlatform == '__WXMSW__':
            faces = {'times': cn, 'mono' : cn, 'helv' : cn, 'other': cn,
                     'size' : 10, 'size2': 9}
        else:
            faces = {'times': 'Courier', 'mono' : 'Courier',
                     'helv' : 'Courier', 'other': 'Courier', 'size' : 10,
                     'size2': 10 }
        del cn

#---------------------- Frame that contains everything -----------------------
class MainWindow(wxFrame):
    def __init__(self,parent,id,title,fnames):
        wxFrame.__init__(self,parent,id, title, size = ( 800, 600 ),
                       style=wxDEFAULT_FRAME_STYLE|wxNO_FULL_REPAINT_ON_RESIZE)
  #      wxFrame.__init__(self,parent,id, title, size = ( 800, 600 ),style=wxNO_FULL_REPAINT_ON_RESIZE)

#        self.SetIcon(getIcon())
        self.FINDSHOWN = 0
        path = os.path.join(homedir, 'menus.txt')
        try:    OLD_MENUPREF.update(self.readAndEvalFile(path))
        except: pass

        #recent menu relocated to load configuration early on.
        recentmenu = wxMenu()
#--------------------------- cmt-001 - 08/06/2003 ----------------------------
#----------------- Adds opened file history to the File menu -----------------
        self.fileHistory = wxFileHistory()
        self.fileHistory.UseMenu(recentmenu)
        self.configPath = homedir
        self.loadHistory()
        self.restart = 0
        self.restartt = 0
        EVT_MENU_RANGE(self, wxID_FILE1, wxID_FILE9, self.OnFileHistory)
        self.lastused = lastused(64+len(lastopen), LASTUSED)
#------------------------- end cmt-001 - 08/06/2003 --------------------------
        self.toparse = []
        self.parsing = 0
        
        #EVT_IDLE(self, self.SetPos)
        #a = wxNewId()
        #self.T = wxTimer(self, a)
        #EVT_TIMER(self, a, self.SetPos)
        #self.T.Start(100)

        self.sb = wxStatusBar(self, -1)
        if VS[-1] == 'u':
            #to display encoding in unicode supporting platforms
            self.sb.SetFieldsCount(3)
            self.sb.SetStatusWidths([-1, 95, 60])
        else:
            self.sb.SetFieldsCount(2)
            self.sb.SetStatusWidths([-1, 95])
        self.SetStatusBar(self.sb)

        #support for the optional snippet bar.
        if self.config['usesnippets']:
            self.split = wxSplitterWindow(self, -1)

            self.control = MyNB(self, -1, self.split)
            self.snippet = CodeSnippet(self, -1, self.split, self.config['display2code'], self.config['displayorder'])

            self.split.SetMinimumPaneSize(1)
            self.split.SplitVertically(self.snippet, self.control, 1)
        else:
            self.control = MyNB(self, -1, self)

        EVT_CLOSE(self, self.OnExit)

        # Setting up the menu.

#------------------------- Insert menus into Menubar -------------------------
        global menuBar
        menuBar = wxMenuBar()

        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        self.menubar = menuBar

#--------------------------------- File Menu ---------------------------------

        filemenu= wxMenu()
        menuAddM(menuBar, filemenu, "&File")
        menuAdd(self, filemenu, "&New\tCtrl+N",         "New file", self.OnNew, wxID_NEW)
        menuAdd(self, filemenu, "&Open\tCtrl+O",        "Open a file", self.OnOpen, wxID_OPEN)
        menuAdd(self, filemenu, "Open &Module\tAlt+M",  "Open a module for editing using the same path search as import would", self.OnOpenModule, wxNewId())
        menuAdd(self, filemenu, "Open &Last\t",         "Open all the documents that were opening before last program exit", self.OnOpenPrevDocs, wxNewId())
        menuAddM(filemenu, recentmenu, "Open Recent")
        filemenu.AppendSeparator()
        menuAdd(self, filemenu, "&Save\tCtrl+S",        "Save a file", self.OnSave, wxID_SAVE)
        menuAdd(self, filemenu, "Save &As",             "Save a file as...", self.OnSaveAs, wxID_SAVEAS)
        menuAdd(self, filemenu, "Sa&ve All",            "Save all open files...", self.OnSaveAll, wxNewId())
        filemenu.AppendSeparator()
        menuAdd(self, filemenu, "Add Module Search Path", "Add a path to search during subsequent 'Open Module' executions", self.AddSearchPath, wxNewId())
        menuAdd(self, filemenu, "&Reload",              "Reload the current document from disk", self.OnReload, wxID_REVERT)
        menuAdd(self, filemenu, "&Close\tCtrl+W",        "Close the file in this tab", self.OnClose, wxNewId())
        filemenu.AppendSeparator()
        menuAdd(self, filemenu, "E&xit\tAlt+F4",        "Terminate the program", self.OnExit, wxNewId())

#--------------------------------- Edit Menu ---------------------------------

        editmenu= wxMenu()
        menuAddM(menuBar, editmenu, "&Edit")
        menuAdd(self, editmenu, "Undo\tCtrl+Z",         "Undo last modifications", self.OnUndo, wxNewId())
        menuAdd(self, editmenu, "Redo\tCtrl+Y",         "Redo last modifications", self.OnRedo, wxNewId())
        editmenu.AppendSeparator()
        menuAdd(self, editmenu, "Select All\tCtrl+A",   "Select all text", self.OnSelectAll, wxNewId())
        menuAdd(self, editmenu, "Cut\tCtrl+X",          "Cut selected text", self.OnCut, wxNewId())
        menuAdd(self, editmenu, "Copy\tCtrl+C",         "Copy selected text", self.OnCopy, wxNewId())
        menuAdd(self, editmenu, "Paste\tCtrl+V",        "Paste selected text", self.OnPaste, wxNewId())
        editmenu.AppendSeparator()
        menuAdd(self, editmenu, "Indent Region\tCtrl+]", "Indent region %i spaces%indent", self.OnIndent, wxNewId())
        menuAdd(self, editmenu, "Dedent Region\tCtrl+[", "Dedent region %i spaces%indent", self.OnDedent, wxNewId())
        editmenu.AppendSeparator()
        menuAdd(self, editmenu, "Find...\tCtrl+F",      "Find text in a file", self.OnShowFind, wxNewId())
        menuAdd(self, editmenu, "Find again\tF3",       "Find again same text in a file", self.OnShowFindAgain, wxNewId())
        menuAdd(self, editmenu, "Find In Files",        "Find a string or regular expression in a group of files", self.OnFindInFiles, wxNewId())
        menuAdd(self, editmenu, "Replace\tCtrl+R",      "Replace text in a file", self.OnShowFindReplace, wxNewId())
        menuAdd(self, editmenu, "Wrap Selected Text\tAlt+W", "Wrap selected text to a specified width", self.OnWrap, wxNewId())

        editmenu.AppendSeparator()
        if self.config['usesnippets']:
            menuAdd(self, editmenu, "Insert Snippet\tCtrl+return", "Insert the currently selected snippet into the document", self.snippet.OnListBoxDClick, wxNewId())
        editmenu.AppendSeparator()
        menuAdd(self, editmenu, "Insert Comment\tCtrl+I", "Insert a centered comment", self.OnInsertComment, wxNewId())
        menuAdd(self, editmenu, "Comment Selection\tAlt+8", "Comment selected lines", self.OnCommentSelection, wxNewId())
        menuAdd(self, editmenu, "Uncomment Selection\tAlt+9", "Uncomment selected lines", self.OnUncommentSelection, wxNewId())

        EVT_COMMAND_FIND(self, -1, self.OnFind)
        EVT_COMMAND_FIND_NEXT(self, -1, self.OnFind)
        EVT_COMMAND_FIND_REPLACE(self, -1, self.OnFind)
        EVT_COMMAND_FIND_REPLACE_ALL(self, -1, self.OnFind)
        EVT_COMMAND_FIND_CLOSE(self, -1, self.OnFindClose)

#--------------------------------- View Menu ---------------------------------
        
        viewmenu= wxMenu()
        menuAddM(menuBar, viewmenu,"&View")
        menuAdd(self, viewmenu, "Zoom In\tCtrl++", "Make the text in the editing component bigger", self.OnZoom, ZI)
        menuAdd(self, viewmenu, "Zoom Out\tCtrl+-", "Make the text in the editing component smaller", self.OnZoom, wxNewId())
        viewmenu.AppendSeparator()
        menuAdd(self, viewmenu, "Go to line number\tAlt+G", "Advance to the given line in the currently open document", self.OnGoto, wxNewId())
        viewmenu.AppendSeparator()
        menuAdd(self, viewmenu, "Toggle Bookmark\tCtrl+M", "Create/remove bookmark for this line", self.OnToggleBookmark, wxNewId())
        menuAdd(self, viewmenu, "Next Bookmark\tF2", "Hop to the next bookmark in this file", self.OnNextBookmark, wxNewId())
        menuAdd(self, viewmenu, "Previous Bookmark\tShift+F2", "Hop to the previous bookmark in this file", self.OnPreviousBookmark, wxNewId())

#------------------------------- Document menu -------------------------------
        
        setmenu= wxMenu()
        menuAddM(menuBar, setmenu, "&Document")
        menuAdd(self, setmenu, "Use Snippets (req restart)", "Enable or disable the use of snippets, requires restart for change to take effect", self.OnSnipToggle, SNIPT, wxITEM_CHECK)
        menuAdd(self, setmenu, "Use Todo (req restart)", "Enable or disable the use of todo, requires restart for change to take effect", self.OnTodoToggle, TODOT, wxITEM_CHECK)
        setmenu.AppendSeparator()

        #-------------------------------- Style subenu -----------------------
        stylemenu= wxMenu()
        menuAddM(setmenu, stylemenu, "Syntax Highlighting", "Change the syntax highlighting for the currently open document")
        try:
            typ = wxITEM_RADIO
            menuAdd(self, stylemenu, "Python",      "Highlight for Python syntax", self.OnStyleChange, PY_S, typ)
            self.HAS_RADIO = 1
        except:
            #to handle the cases when a platform lacks radio items in menus
            typ = wxITEM_NORMAL
            menuAdd(self, stylemenu, "Python",      "Highlight for Python syntax", self.OnStyleChange, PY_S, typ)
            self.HAS_RADIO = 0

        menuAdd(self, stylemenu, "HTML",        "Highlight for HTML syntax", self.OnStyleChange, HT_S, typ)
        menuAdd(self, stylemenu, "XML",         "Highlight for XML syntax", self.OnStyleChange, XM_S, typ)
        menuAdd(self, stylemenu, "C/C++",       "Highlight for C/C++ syntax", self.OnStyleChange, CC_S, typ)
        menuAdd(self, stylemenu, "Text",        "No Syntax Highlighting", self.OnStyleChange, TX_S, typ)

        #---------------------------- Default Style submenu ------------------
        stylemenu2= wxMenu()
        menuAddM(setmenu, stylemenu2, "Default Highlighting", "Set the default syntax highlighting for new or unknown documents")
        menuAdd(self, stylemenu2, 'Python', '', self.OnDefaultStyleChange, PY_DS, typ)
        menuAdd(self, stylemenu2, 'HTML', '', self.OnDefaultStyleChange, HT_DS, typ)
        menuAdd(self, stylemenu2, 'XML', '', self.OnDefaultStyleChange, XM_DS, typ)
        menuAdd(self, stylemenu2, 'C/C++', '', self.OnDefaultStyleChange, CC_DS, typ)
        menuAdd(self, stylemenu2, 'Text', '', self.OnDefaultStyleChange, TX_DS, typ)

        #------------------------------ Encodings submenu --------------------
        if VS[-1] == 'u':
            encmenu= wxMenu()
            menuAddM(setmenu, encmenu, "Encodings", "Change text encoding")
            menuAdd(self, encmenu, 'ascii', "Change encoding for the current file to ascii (will use utf-8 if unicode characters found)", self.OnEncChange, ENCODINGS['ascii'], typ)
            for bom, enc in BOM[:-1]:
                menuAdd(self, encmenu, enc, "Change encoding for the current file to %s"%enc, self.OnEncChange, ENCODINGS[enc], typ)

        #----------------------------- Line ending menu ----------------------
        endingmenu = wxMenu()
        menuAddM(setmenu, endingmenu, "Line Ending", "Change the line endings on the current document")
        menuAdd(self, endingmenu, "CRLF (windows)", "", self.OnLineEndChange, LE_CRLF, typ)
        menuAdd(self, endingmenu, "LF (*nix)", "", self.OnLineEndChange, LE_LF, typ)
        menuAdd(self, endingmenu, "CR (mac)", "", self.OnLineEndChange, LE_CR, typ)
        #
        
        setmenu.AppendSeparator()
        menuAdd(self, setmenu, "Show Autocomplete", "Show the autocomplete dropdown while typing", self.OnAutoCompleteToggle, AUTO, wxITEM_CHECK)
        menuAdd(self, setmenu, "Show line numbers", "Show or hide the line numbers on the current document", self.OnNumberToggle, NUM, wxITEM_CHECK)
        menuAdd(self, setmenu, "Show margin", "Show or hide the bookmark signifier margin on the current document", self.OnMarginToggle, MARGIN, wxITEM_CHECK)
        menuAdd(self, setmenu, "Show Indentation Guide", "Show or hide gray indentation guides in indentation", self.OnIndentGuideToggle, INDENTGUIDE, wxITEM_CHECK)
        menuAdd(self, setmenu, "Save Position", "Remember or forget the last position of the cursor when the current document is closed", self.OnSavePositionToggle, SAVE_CURSOR, wxITEM_CHECK)
        setmenu.AppendSeparator()
        menuAdd(self, setmenu, "Show/hide tree\tCtrl+Shift+G", "Show/hide the hierarchical source tree for the currently open document", self.OnTree, wxNewId())
        menuAdd(self, setmenu, "Hide all trees", "Hide the browsable source tree for all open documents", self.OnTreeHide, wxNewId())
        menuAdd(self, setmenu, "Refresh\tF5", "Refresh the browsable source tree, autocomplete listing, and the tooltips (always accurate, but sometimes slow)", self.OnRefresh, wxNewId())
        setmenu.AppendSeparator()
        menuAdd(self, setmenu, "Sort Tree by Name", "If checked, will sort the items in the browsable source tree by name, otherwise by line number", self.OnTreeSortToggle, SORTBY, wxITEM_CHECK)
        menuAdd(self, setmenu, "Expand all", "Expand all folded code through the entire document", self.OnExpandAll, wxNewId())
        menuAdd(self, setmenu, "Fold all", "Fold all expanded code through the entire document", self.OnFoldAll, wxNewId())
        menuAdd(self, setmenu, "Use Tabs", "New indentation will include tabs", self.OnSetTabToggle, USETABS, wxITEM_CHECK)
        menuAdd(self, setmenu, "Wrap Long Lines", "Visually continue long lines to the next line", self.OnWrapL, WRAPL, wxITEM_CHECK)
        setmenu.AppendSeparator()
        menuAdd(self, setmenu, "Indent Width", "Set the number of spaces per indentation level", self.OnSetIndent, wxNewId())
        menuAdd(self, setmenu, "Set Tab Width", "Set the visual width of tabs in the current open document", self.OnSetTabWidth, wxNewId())
        menuAdd(self, setmenu, "Set Long Line Column", "Set the column number for the long line indicator", self.OnSetLongLinePosition, wxNewId())

        #------------------------------ Long line submenu --------------------
        longlinemenu = wxMenu()
        menuAddM(setmenu, longlinemenu, "Set Long Line Indicator", "Change the mode that signifies long lines")
        menuAdd(self, longlinemenu, "Background", "Long lines will have a different background color beyond the column limit", self.OnSetLongLineMode, LL_BACK, typ)
        menuAdd(self, longlinemenu, "Line", "Long lines will have a vertical line at the column limit", self.OnSetLongLineMode, LL_LINE, typ)
        menuAdd(self, longlinemenu, "None", "Show no long line indicator", self.OnSetLongLineMode, LL_NONE, typ)
        #
        
        setmenu.AppendSeparator()
        menuAdd(self, setmenu, "Save settings", "Save the settings for the current document as the default for all documents (excluding syntax and encodings)", self.OnSaveSettings, wxNewId())        

#--------------------------------- Tab Menu ----------------------------------
        
        tabmenu= wxMenu()
        menuAddM(menuBar, tabmenu, "&Tabs")
        menuAdd(self, tabmenu, "Previous Tab\tAlt+,",          "View the tab to the left of the one you are currently", self.OnLeft, wxNewId())
        menuAdd(self, tabmenu, "Next Tab\tAlt+.",              "View the tab to the right of the one you are currently", self.OnRight, wxNewId())
        tabmenu.AppendSeparator()
        menuAdd(self, tabmenu, "Move tab left\tCtrl+Alt+,",     "Swap the current tab with the one on the left", self.MoveLeft, wxNewId())
        menuAdd(self, tabmenu, "Move tab right\tCtrl+Alt+.",    "Swap the current tab with the one on the right", self.MoveRight, wxNewId())

#------------------------------- Snippet Menu --------------------------------
        if self.config['usesnippets']:
            snippetmenu= wxMenu()
            self.OnSnipToggle(None)
            menuAddM(menuBar, snippetmenu, "Sn&ippets")
            menuAdd(self, snippetmenu, "Previous snippet\tCtrl+Shift+,",        "Select the previous code snippet", self.snippet.OnSnippetP, wxNewId())
            menuAdd(self, snippetmenu, "Next snippet\tCtrl+Shift+.",            "Select the next code snippet", self.snippet.OnSnippetN, wxNewId())
            menuAdd(self, snippetmenu, "Insert Snippet\tCtrl+return",           "Insert the currently selected snippet into the document", self.snippet.OnListBoxDClick, wxNewId())
            menuAdd(self, snippetmenu, "Add copied text to snippets\tCtrl+Shift+V", "Insert the currently copied text (in the clipboard) into the snippet bar", self.snippet.onpaste, wxNewId())
            menuAdd(self, snippetmenu, "Show/hide snippet bar\tCtrl+Shift+B",   "Show/hide the global code snippet bar on the left", self.OnSnippet, wxNewId())

#------------------------------- Pathmark Menu -------------------------------

        pathmarkmenu = wxMenu()
        self.pm = pathmarkmenu
        menuAddM(menuBar, pathmarkmenu, "&Pathmarks")
        menuAdd(self, pathmarkmenu, "View pathmarks",           "Edit your Pathmarks", self.ViewPathmarks, wxNewId())
        menuAdd(self, pathmarkmenu, "Edit pathmarks\tCtrl+B",   "Add a path to your bookmarks", self.AddPathmark, wxNewId())
        menuAdd(self, pathmarkmenu, "Remove pathmark",          "Remove a bookmarked path", self.RemovePathmark, wxNewId())
        pathmarkmenu.AppendSeparator()

        pmk = [BM1,BM2,BM3,BM4,BM5,BM6,BM7,BM8,BM9]
        for i in xrange(49, 58):
            if pathmarks.get(i, 0) != 0:
                menuAdd(self, pathmarkmenu,
                        "Ctrl+%i\t%s"%(i-48, pathmarks[i]),
                        "Change the current working directory to %s"%pathmarks[i],
                        self.OnPathmark, pmk[i-49])
            else:
                EVT_MENU(self, pmk[i-49], self.OnPathmark)

#-------------------------------- Shell Menu ---------------------------------

        self.shell = RunShell(self, menuBar, "&Shell")

#--------------------------------- Help Menu ---------------------------------
        
        helpmenu= wxMenu()
        menuAddM(menuBar, helpmenu, "&Help")
        menuAdd(self, helpmenu, "About...", "About this piece of software", self.OnAbout, wxID_ABOUT)
        menuAdd(self, helpmenu, "Change Menus (req restart)", "Change the name of menu items and their hotkeys", self.OnChangeMenu, wxNewId())
        helpmenu.AppendSeparator()
        menuAdd(self, helpmenu, "PyPE Help\tF1", "View the help", self.OnHelp, wxID_HELP)

#------------------------ A ...few... state variables ------------------------

        self.Show(true)
        self.dirname = '.'
        self.closing = 0
        self.openfiles = {}
        self.dpm = 0
        self.menubar.Check(SNIPT, usesnippets)
        self.menubar.Check(TODOT, usetodo)
        self.menubar.Check(AUTO, showautocomp)
        self.menubar.Check(WRAPL, wrapmode != wxSTC_WRAP_NONE)
        self.menubar.Check(USETABS, use_tabs)
        self.menubar.Check(INDENTGUIDE, indent_guide)
        self.menubar.Check(SORTBY, sortmode)
        self.menubar.Check(lexers3[DEFAULTLEXER], 1)
        self.menubar.Check(SAVE_CURSOR, save_cursor)

#------------------------ Drag and drop file support -------------------------
        self.SetDropTarget(FileDropTarget(self))

#------------------ Open files passed as arguments to PyPE -------------------
        self.OnDrop(fnames, 0)
        EVT_SIZE(self, self.OnResize)
        EVT_ACTIVATE(self, self.OnActivation)
        EVT_KEY_DOWN(self, self.OnKeyPressed)

    def dialog(self, message, title, styl=wxOK):
        d= wxMessageDialog(self,message,title,styl)
        retr = d.ShowModal()
        d.Destroy()
        return retr

    def exceptDialog(self, title="Error"):
        k = cStringIO.StringIO()
        traceback.print_exc(file=k)
        k.seek(0)
        dlg = wxScrolledMessageDialog(self, k.read(), title)
        dlg.ShowModal()
        dlg.Destroy()

    def OnDrop(self, fnames, error=1):
        cwd = os.getcwd()
        for i in fnames:
            i = os.path.normcase(os.path.normpath(os.path.realpath(os.path.join(cwd, i))))
            if self.isAbsOpen(i):
                if len(fnames)==1:
                    self.selectAbsolute(i)
            else:
                self.makeAbsOpen(i)
                d,f = os.path.split(i)
                try:
                    a = self.newTab(d,f, len(fnames)==1)
                    if VS[-1] == 'u': a = "%s as %s"%(i, a)
                    else:             a = i
                    self.SetStatusText("Opened %s"%a)
                except:
                    if error:
                        self.exceptDialog("File open failed")

    def SetStatusText(self, text, number=0):
        if (number == 0) and text:
            text = "[%s] %s"%(time.asctime(), text)
        self.sb.SetStatusText(text, number)

    def OnActivation(self, e):
        try:
            self.control.__iter__
        except wxPyDeadObjectError:
            return e.Skip()
        try:
            self.iter
            return e.Skip()
        except:
            self.iter = None
        for document in self.control.__iter__():
            if document.mod != None:
                fn = self.getAbsolute(document.filename, document.dirname)
                try:
                    mod = os.stat(fn)[8]
                except OSError:
                    document.MakeDirty(None)
                    self.dialog("%s\n"\
                                "has been deleted from the disk by an\n"
                                "external program.\n"\
                                "Unless the file is saved again, data\n"
                                "loss may occur."%fn,
                                "WARNING!")
                    document.mod = None
                    continue
                if mod != document.mod:
                    document.MakeDirty(None)
                    a = self.dialog("%s\n"\
                                    "has been modified by an external\n"\
                                    "program.\n"\
                                    "Would you like to reload the file\n"\
                                    "from disk?"%fn,
                                    "WARNING!", wxYES_NO)
                    if a == wxID_NO:
                        document.mod = None
                    else:
                        self.OnReload(None, document)
        del self.iter
        e.Skip()
                    
#--------------------------- cmt-001 - 08/06/2003 ----------------------------
#--------------------- Saving and loading of preferences ---------------------
    def loadHistory(self):
        if not os.path.exists(self.configPath):
            os.mkdir(self.configPath)
        path = os.path.join(self.configPath, 'history.txt')
        try:    self.config = self.readAndEvalFile(path)
        except: self.config = {}
        if 'history' in self.config:
            history = self.config['history']
            history.reverse()
            for h in history:
                self.fileHistory.AddFileToHistory(h)
        if 'lastpath' in self.config:
            self.lastpath = self.config['lastpath']

        for (nam, dflt) in [('modulepaths', []),
                            ('usesnippets', 0),
                            ('usetodo', 0),
                            ('paths', {}),
                            ('display2code', {}),
                            ('displayorder', []),
                            ('shellcommands', []),
                            ('lastopen', []),
                            ('LASTUSED', []),
                            ('DEFAULTLEXER', 'python'),
                            ('splittersize', 40),
                            ('FIF_STATE', ([], [], [], 0, 0, 0))]:
            if not (nam in self.config):
                self.config[nam] = dflt
            globals()[nam] = self.config[nam]

        #insert global defaults here
        dct =   {'use_tabs':0,
           'spaces_per_tab':8,
                   'indent':4,
                 'collapse':1,
            'marker_margin':1,
              'line_margin':1,
                 'col_line':78,
                 'col_mode':wxSTC_EDGE_LINE,
             'indent_guide':0,
             'showautocomp':0,
                 'wrapmode':wxSTC_WRAP_NONE,
                 'sortmode':1,
            'SAVEDPOSITION':splittersize,
              'save_cursor':0,
              'cursor_posn':0}
        
        globals().update(dct)
        globals().update(self.config.setdefault('DOCUMENT_DEFAULTS', dct))

        pathmarks.update(paths)

    def saveHistory(self):
        history = []
        for i in range(self.fileHistory.GetNoHistoryFiles()):
            history.append(self.fileHistory.GetHistoryFile(i))
        self.config['history'] = history
        a = []
        for i in self.shell.order:
            a.append(self.shell.menu[i])
        self.config['shellcommands'] = a
        self.config['paths'] = pathmarks
        if self.config['usesnippets'] and (not self.restart):
            self.config['display2code'] = self.snippet.display2code
            self.config['displayorder'] = self.snippet.displayorder
        self.config['lastpath'] = self.config.get('lp', os.getcwd())
        try:
            path = os.sep.join([self.configPath, 'history.txt'])
            f = open(path, "w")
            f.write(pprint.pformat(self.config))
            f.close()
            path = os.sep.join([self.configPath, 'menus.txt'])
            f = open(path, "w")
            f.write(pprint.pformat(MENUPREF))
            f.close()
        except:
            self.exceptDialog("Could not save preferences to %s"%path)

    def readAndEvalFile(self, filename):
        f = open(filename)
        txt = f.read().replace('\r\n','\n')
        f.close()
        return eval(txt)
#------------------------- end cmt-001 - 08/06/2003 --------------------------

    def redrawvisible(self, win=None):
        if (not win) and self.control.GetPageCount() > 0:
            num, win = self.getNumWin()
        if win:
            a = win.parent.GetSashPosition()
            win.parent.SetSashPosition(a-1)
            win.parent.SetSashPosition(a)
            if win.todo:
                wtp = win.todo.parent
                a = wtp.GetSashPosition()
                wtp.SetSashPosition(a-1)
                wtp.SetSashPosition(a)

#---------------------------- File Menu Commands -----------------------------
    def selectAbsolute(self, path):
        if self.isAbsOpen(path):
            dn, fn = self.splitAbsolute(path)
            for i in xrange(self.control.GetPageCount()):
                win = self.control.GetPage(i).GetWindow1()
                if (win.filename == fn) and (win.dirname == dn):
                    self.control.SetSelection(i)
                    return

    def isOpen(self, fn, dn):
        return dn and fn and (self.getAbsolute(fn, dn) in self.openfiles)
    def isAbsOpen(self, path):
        return path in self.openfiles

    def makeOpen(self, fn, dn):
        if fn and dn:
            a = self.getAbsolute(fn, dn)
            self.openfiles[a] = self.splitAbsolute(a)
    def makeAbsOpen(self, path):
        if path:
            self.openfiles[path] = self.splitAbsolute(path)
    def closeOpen(self, fn, dn):
        if fn and dn:
            del self.openfiles[self.getAbsolute(fn, dn)]

    def getAbsolute(self, fn, dn):
        return os.path.normcase(os.path.normpath(os.path.realpath(os.path.join(dn, fn))))
    def splitAbsolute(self, path):
        return os.path.split(os.path.normcase(path))

    def OnNew(self,e):
        self.newTab('', ' ', 1)
        self.control.GetPage(self.control.GetSelection()).GetWindow1().opened = 1
        

    def OnSave(self,e):
        wnum, win = self.getNumWin(e)
        if win.dirname:
            try:
                ofn = os.path.join(win.dirname, win.filename)
                fil = open(ofn, 'wb')
                txt = win.GetText()
                fil.write(txt)
                fil.close()
                if VS[-1] == 'u': a = "%s as %s"%(ofn, win.enc)
                else:             a = ofn
                win.mod = os.stat(ofn)[8]
                self.SetStatusText("Correctly saved %s"%a)
                win.MakeClean()
            except:
                self.exceptDialog("Save Failed")
                raise cancelled
        else:
            self.OnSaveAs(e)

    def OnSaveAs(self,e):
        wnum, win = self.getNumWin(e)
        
        dlg = wxFileDialog(self, "Save file as...", os.getcwd(), "", "All files (*.*)|*.*", wxSAVE|wxOVERWRITE_PROMPT)
        rslt = dlg.ShowModal()
        if rslt == wxID_OK:
            fn=dlg.GetFilename()
            dn=dlg.GetDirectory()
            if self.isOpen(fn, dn):
                self.dialog("Another file with that name and path is already open.\nSave aborted to prevent data corruption.", "Save Aborted!")
                raise cancelled
            if self.isOpen(win.filename, win.dirname):
                self.closeOpen(win.filename, win.dirname)
            dn, fn = self.splitAbsolute(self.getAbsolute(fn, dn))
            win.filename = fn
            win.dirname = dn
            self.makeOpen(fn, dn)
            self.fileHistory.AddFileToHistory(self.getAbsolute(fn, dn))
            self.OnSave(e)
            win.MakeDirty()
            win.MakeClean()
        else:
            raise cancelled

    def OnSaveAll(self, e):
        sel = self.control.GetSelection()
        cnt = self.control.GetPageCount()
        for i in xrange(cnt):
            self.control.SetSelection(i)
            try:
                self.OnSave(e)
            except cancelled:
                pass
        if cnt:
            self.control.SetSelection(sel)
    
    def OnOpen(self,e):
        wd = self.config.get('lastpath', os.getcwd())
        dlg = wxFileDialog(self, "Choose a/some file(s)...", wd, "", wildcard, wxOPEN|wxMULTIPLE|wxHIDE_READONLY)
        if dlg.ShowModal() == wxID_OK:
            self.OnDrop(dlg.GetPaths())
            self.config['lp'] = dlg.GetDirectory()
        dlg.Destroy()

    def FindModule(self, mod):
        fndmod = mod.split('.')
        lf = len(fndmod)        
        pth = sys.path[:]
        pth[1:1] = [os.getcwd()]
        for j in range(lf):
            i = fndmod[j]
            fdsc = imp.find_module(i, pth)
            if not (fdsc[0] is None):
                fdsc[0].close()
                mod = fdsc[1]
                if fdsc[2][2] != imp.PY_SOURCE:
                    return self.dialog("%s is not python source"%mod, "not correct format for editing")
                return mod
            elif fdsc[1]:
                pth[1:1] = [fdsc[1]]
            else:
                raise cancelled
        #If we are outside the loop, this means that the current 'module'
        #that we are on is a folder-module.  We can easily load the __init__.py
        #from this folder.
        return os.sep.join([pth[1], '__init__.py'])

    def OnOpenModule(self,e):
        dlg = wxTextEntryDialog(self, 'Enter the module name you would like to open', 'Open module...')
        if dlg.ShowModal() == wxID_OK:
            mod = dlg.GetValue()
        else:
            mod = ''
        dlg.Destroy()
        if mod:
            sp = sys.path[:]
            sys.path.extend(self.config['modulepaths'])
            try:
                self.OnDrop([self.FindModule(mod)])
            except:
                self.dialog("module %s not found"%mod, "not found")
            sys.path = sp

    def OnOpenPrevDocs(self, e):
        if "lastopen" in self.config:
            self.OnDrop(self.config['lastopen'], 1)                
    def AddSearchPath(self, e):
        dlg = wxDirDialog(self, "Choose a path", "", style=wxDD_DEFAULT_STYLE|wxDD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wxID_OK:
            path = os.path.normcase(os.path.normpath(dlg.GetPath()))
            if not (path in self.config['modulepaths']) and not (path in sys.path):
                self.config['modulepaths'].append(path)

    def newTab(self, d, fn, switch=0):
        if 'lastpath' in self.config:
            del self.config['lastpath']
        #ctrlwidth, ctrlh = self.control.GetSizeTuple()

        split = wxSplitterWindow(self.control, wxNewId())
        split.parent = self
        EVT_SPLITTER_SASH_POS_CHANGED(split, split.GetId(), self.SetPos)
        nwin = PythonSTC(self.control, wxNewId(), split)
        nwin.split = split
        nwin.filename = fn
        nwin.dirname = d
        nwin.changeStyle(stylefile, self.style(fn))
        if usetodo:
            spl = wxSplitterWindow(split, wxNewId())
            nwin.tree = hierCodeTreePanel(self, spl)
            nwin.todo = VirtualTodo(spl, self, nwin)
            
            spl.SetMinimumPaneSize(3)
            spl.SplitHorizontally(nwin.tree, nwin.todo)
            
        else:
            nwin.tree = hierCodeTreePanel(self, split)

            nwin.todo = None
            spl = nwin.tree
        
        split.SetMinimumPaneSize(30)
        split.SplitVertically(nwin, spl)
        state = self.config.get('DOCUMENT_DEFAULTS', {})
        if d:
            FN = self.getAbsolute(nwin.filename, nwin.dirname)
            f=open(FN,'rb')
            txt = f.read()
            f.close()
            nwin.mod = os.stat(FN)[8]
            nwin.format = detectLineEndings(txt)
            nwin.SetText(txt)
            
            #if FN in self.config['LASTOPEN']:
            #    state = self.config['LASTOPEN'][FN]
            if FN in self.lastused:
                state = self.lastused[FN]
                del self.lastused[FN]

        else:
            FN = ''
            nwin.mod = None
            nwin.format = eol
            nwin.SetText('')
        if not ((d == '') and (fn == ' ')):
            self.fileHistory.AddFileToHistory(FN)
        nwin.SetSaveState(state)
        f = nwin.filename
        if f == ' ':
            globals()['NEWDOCUMENT'] += 1
            nwin.NEWDOCUMENT = NEWDOCUMENT
            f = '<untitled %i>'%NEWDOCUMENT
            
        self.control.AddPage(split, f, switch)
        self.OnRefresh(None, nwin)
        return nwin.enc

    def OnReload(self, e, win=None):
        if win == None:
            num, win = self.getNumWin(e)
        if not e is None:
            dlg = wxMessageDialog(self, "%s was modified after last save.\nReloading from disk will destroy all changes.\n\nContinue anyway?"%win.filename, 'File was modified, data loss may occur!', wxYES_NO|wxCANCEL)
            a = dlg.ShowModal()
            if a == wxID_CANCEL:
                raise cancelled
            elif a == wxID_NO:
                return
        try:
            FN = self.getAbsolute(win.filename, win.dirname)
            fil = open(FN, 'rb')
            txt = fil.read()
            fil.close()
            win.BeginUndoAction()
            win.SetText(txt, 0)
            win.EndUndoAction()
            win.mod = os.stat(FN)[8]
            win.SetSavePoint()
            win.MakeClean()
            self.OnRefresh(None, win)
        except:
            self.dialog("Error encountered while trying to reload from disk.", "Reload failed")

    def sharedsave(self, win):
        nam = win.filename
        if not win.dirname:
            nam = "<untitled %i>"%win.NEWDOCUMENT
        a = self.dialog("%s was modified after last save.\nSave changes before closing?"%nam,\
                        "Save changes?", wxYES_NO|wxCANCEL)
        if a == wxID_CANCEL:
            raise cancelled
        elif a == wxID_NO:
            pass
        else:
            self.OnSave(None)

    def OnClose(self, e):
        wnum, win = self.getNumWin(e)
        if win.dirty:
            self.sharedsave(win)
        if self.isOpen(win.filename, win.dirname):
            fn = self.getAbsolute(win.filename, win.dirname)
            self.lastused[fn] = win.GetSaveState()
            self.closeOpen(win.filename, win.dirname)
            self.SetStatusText("Closed %s"%self.getAbsolute(win.filename, win.dirname))
        else:
            self.SetStatusText("Closed unnamed file without saving")
        self.control.DeletePage(wnum)
                    
    def OnExit(self,e):
        if self.closing:
            return e.Skip()
        self.closing = 1
        sel = self.control.GetSelection()
        cnt = self.control.GetPageCount()
        try:
            for i in xrange(cnt):
                win = self.control.GetPage(i).GetWindow1()
                if win.dirty:
                    self.control.SetSelection(i)
                    self.sharedsave(win)
        except cancelled:
            self.closing = 0
            try:    return e.Veto()
            except: return e.Skip()

        LASTOPEN = []
        sav = []
        for win in self.control:
            if win.dirname:
                fn = self.getAbsolute(win.filename, win.dirname)
                sav.append(fn)
                LASTOPEN.append((fn, win.GetSaveState()))
        
        if 'LASTOPEN' in self.config:
            del self.config['LASTOPEN']
        #saving document state
        self.config["lastopen"] = sav
        self.config["LASTUSED"] = self.lastused.items() + LASTOPEN

        self.saveHistory()
        if sel > -1:
            self.control.SetSelection(sel)
        try:    os.remove(OUTF)
        except: pass
        try:    os.remove(INF)
        except: pass
        return self.Close(true)

#--------------------------- cmt-001 - 08/06/2003 ----------------------------
#------------- Open the file selected from the file history menu -------------
    def OnFileHistory(self, e):
        fileNum = e.GetId() - wxID_FILE1
        path = self.fileHistory.GetHistoryFile(fileNum)
        self.OnDrop([path])
#------------------------- end cmt-001 - 08/06/2003 --------------------------

#---------------------------- Edit Menu Commands -----------------------------
    def OneCmd(self, funct_name,evt):
        wnum, win = self.getNumWin(evt)
        getattr(win, funct_name)()

    def OnUndo(self, e):
        self.OneCmd('Undo',e)
    def OnRedo(self, e):
        self.OneCmd('Redo',e)
    def OnSelectAll(self, e):
        self.OneCmd('SelectAll',e)
    def OnCut(self, e):
        self.OneCmd('Cut',e)
    def OnCopy(self, e):
        self.OneCmd('Copy',e)
    def OnPaste(self, e):
        self.OneCmd('Paste',e)
    def OnWrap(self, e):
        wnum, win = self.getNumWin(e)
        
        dlg = wxTextEntryDialog(self, '', 'Wrap to how many columns?', str(col_line))
        resp = dlg.ShowModal()
        valu = dlg.GetValue()
        dlg.Destroy()
        try:
            if resp != wxID_OK:
                raise
            valu = int(valu)
        except:
            return
        win.MakeDirty()
        x,y = win.GetSelection()
        if x==y:
            return
        lnstart = win.LineFromPosition(x)
        lnend = win.LineFromPosition(y-1)
        lines = []
        for ln in xrange(lnstart, lnend+1):
            lines.append(win.GetLine(ln))
        x = win.GetLineEndPosition(lnstart)-len(lines[0])
        y = win.GetLineEndPosition(lnend)
        #win.SetSelection(x, y)
        win.ReplaceSelection(wrap_lines(''.join(lines), valu, win.format))

    def Dent(self, win, incr):
        x,y = win.GetSelection()
        if x==y:
            lnstart = win.GetCurrentLine()
            lnend = lnstart
            if incr < 0:
                a = win.GetLineIndentation(lnstart)%(abs(incr))
                if a:
                    incr = -a
            pos = win.GetCurrentPos()
            col = win.GetColumn(pos)
            linestart = pos-col
            a = max(linestart+col+incr, linestart)
        else:
            lnstart = win.LineFromPosition(x)
            lnend = win.LineFromPosition(y-1)
        win.BeginUndoAction()
        for ln in xrange(lnstart, lnend+1):
            count = win.GetLineIndentation(ln)
            m = (count+incr)
            m += cmp(0, incr)*(m%incr)
            m = max(m, 0)
            win.SetLineIndentation(ln, m)
        if x==y:
            pos = pos + (m-count) - min(0, col + (m-count))
            win.SetSelection(pos, pos)
        win.EndUndoAction()
    def OnIndent(self, e):
        wnum, win = self.getNumWin(e)
        self.Dent(win, win.GetIndent())
    def OnDedent(self, e):
        wnum, win = self.getNumWin(e)
        self.Dent(win, -win.GetIndent())
    def OnInsertComment(self, e):
        wnum, win = self.getNumWin(e)
        dlg = wxTextEntryDialog(self, '', 'Enter a comment.', '')
        resp = dlg.ShowModal()
        valu = dlg.GetValue()
        dlg.Destroy()
        if resp == wxID_OK:
            k = len(valu)
            a = col_line-3-k
            b = a*'-'
            st = '%s%s %s %s%s'%('#', b[:a/2], valu, b[a/2:], win.format)
            lin = win.GetCurrentLine()
            if lin>0:
                win.InsertText(win.GetLineEndPosition(lin-1)+len(win.format), st)
            else:
                win.InsertText(0, st)
            win.MakeDirty()
        else:
            e.Skip()

    def OnCommentSelection(self, e):
        wnum, win = self.getNumWin(e)
        sel = win.GetSelection()
        start = win.LineFromPosition(sel[0])
        end = win.LineFromPosition(sel[1])
        if end > start and win.GetColumn(sel[1]) == 0:
            end = end - 1
        win.MakeDirty()
        win.BeginUndoAction()
        for lineNumber in range(start, end + 1):
            firstChar = win.GetLineIndentPosition(lineNumber)
            lastChar = win.GetLineEndPosition(lineNumber)
            ranga = win.GetTextRange(firstChar,lastChar)
            if len(ranga.strip()) != 0:
                win.InsertText(firstChar, '## ')
        win.SetCurrentPos(win.PositionFromLine(start))
        win.SetAnchor(win.GetLineEndPosition(end))
        win.EndUndoAction()

    def OnUncommentSelection(self, e):
        wnum, win = self.getNumWin(e)
        sel = win.GetSelection()
        start = win.LineFromPosition(sel[0])
        end = win.LineFromPosition(sel[1])
        if end > start and win.GetColumn(sel[1]) == 0:
            end = end - 1
        win.MakeDirty()
        win.BeginUndoAction()
        for lineNumber in range(start, end + 1):
            firstChar = win.GetLineIndentPosition(lineNumber)
            texta = win.GetTextRange(firstChar,firstChar+3)
            lengtha = 0
            if texta[0:3] == '## ':
                lengtha = 3
            elif texta[0:2] == '##':
                lengtha = 2
            elif texta[0:1] == '#':
                lengtha = 1
            win.SetSelection(firstChar,firstChar+lengtha)
            win.ReplaceSelection("")

        win.SetCurrentPos(win.PositionFromLine(start))
        win.SetAnchor(win.GetLineEndPosition(end))
        win.EndUndoAction()

#--------------------- Find and replace dialogs and code ---------------------
    def OnFindInFiles(self, e):
        findinfiles(self, self).Show()

    def getNumWin(self, e=None):
        num = self.control.GetSelection()
        if num >= 0:
            return num, self.control.GetPage(num).GetWindow1()
        return num,0
        if e:
            e.Skip()
        raise cancelled

    def makedata(self, evt):
        data = wxFindReplaceData()
        num, win = self.getNumWin(evt)
        gcp = win.GetCurrentPos()
        st = win.WordStartPosition(gcp, 1)
        end = win.WordEndPosition(gcp, 1)
        if st != end:
            data.SetFindString(win.GetTextRange(st, end))
        return data

    def OnShowFind(self, evt):
        if self.FINDSHOWN: return
        else: self.FINDSHOWN = 1
        wcount = self.control.GetPageCount()
        if not wcount:
            return evt.Skip()
        for wnum in xrange(wcount):
            win = self.control.GetPage(wnum).GetWindow1()
            win.gcp = win.GetCurrentPos()
            #print win.gcp, "found gcp"
            win.last = 0
        data = self.makedata(evt)
        dlg = wxFindReplaceDialog(self, data, "Find")
        dlg.data = data
        dlg.Show(True)

    def OnShowFindReplace(self, evt):
        if self.FINDSHOWN: return
        else: self.FINDSHOWN = 1
        wcount = self.control.GetPageCount()
        if not wcount:
            return evt.Skip()
        for wnum in xrange(wcount):
            win = self.control.GetPage(wnum).GetWindow1()
            win.gcp = win.GetCurrentPos()
            win.last = 0
        data = self.makedata(evt)
        dlg = wxFindReplaceDialog(self, data, "Find & Replace", wxFR_REPLACEDIALOG)
        dlg.data = data
        return dlg.Show(True)

    def OnShowFindAgain(self, evt):
        try: self.findTxt
        except: return
        a = wxFindDialogEvent(wxEVT_COMMAND_FIND_NEXT)
        a.SetFlags(self.flags)
        a.SetFindString(self.findTxt)
        self.OnFind(a)

    def OnFind(self, evt):
        wnum, win = self.getNumWin(evt)
        et = evt.GetEventType()
        m = {
            wxEVT_COMMAND_FIND : "FIND",
            wxEVT_COMMAND_FIND_NEXT : "FIND_NEXT",
            wxEVT_COMMAND_FIND_REPLACE : "REPLACE",
            wxEVT_COMMAND_FIND_REPLACE_ALL : "REPLACE_ALL"}

        if not et in m:
            return evt.Skip()

        findTxt = evt.GetFindString()

        #the next couple lines deal with python strings in find
        if findTxt and (findTxt[-1] in ['"', "'"]):
            try:    findTxt = eval(findTxt)
            except: pass
        self.findTxt = findTxt

        #the next couple lines deal with python strings in replace
        if et == wxEVT_COMMAND_FIND_REPLACE or et == wxEVT_COMMAND_FIND_REPLACE_ALL:
            replaceTxt = evt.GetReplaceString()
            if replaceTxt and (replaceTxt[0] in ['"', "'"]):
                try:    replaceTxt = eval(replaceTxt)
                except: pass
        else:
            replaceTxt = ""

        #next line makes ininite loops not happen
        incr = len(replaceTxt)
        #next line moves cursor for replacements
        diff = incr-len(findTxt)

        flags = evt.GetFlags()
        if VS[-1] == 'u':
            #searching up causes a crash when searching for unicode strings
            flags = flags|wxFR_DOWN

        self.flags = flags

        if et == wxEVT_COMMAND_FIND_REPLACE_ALL:
            totl = 0
            win.BeginUndoAction()
            while 1:
                win.last = win.FindText(win.last, win.GetTextLength(), findTxt, flags)
                if win.last > -1:
                    totl += 1
                    win.SetSelection(win.last, win.last+len(findTxt))
                    win.ReplaceSelection(replaceTxt)
                    win.last += incr
                    if win.last < win.gcp:
                        win.gcp += diff
                    win.MakeDirty()
                else:
                    win.SetSelection(win.gcp, win.gcp)
                    win.ScrollToColumn(0)
                    win.EnsureCaretVisible()
                    win.EndUndoAction()
                    return self.dialog("%i replacements made."%totl, "Finished replacing")
        elif et == wxEVT_COMMAND_FIND_REPLACE:
            win.last = win.FindText(win.last, win.GetTextLength(), findTxt, flags)
            if win.last > -1:
                win.SetSelection(win.last, win.last+len(findTxt))
                win.ReplaceSelection(replaceTxt)
                win.last += incr
                if win.last < win.gcp:
                    win.gcp += diff
                win.MakeDirty()
                a = win.FindText(win.last, win.GetTextLength(), findTxt, flags)
                win.SetSelection(a, a+len(findTxt))
        elif et == wxEVT_COMMAND_FIND:
            if wxFR_DOWN&flags:
                win.last = win.FindText(0, win.GetTextLength(), findTxt, flags)
            else:
                win.last = win.FindText(win.GetTextLength(), 0, findTxt, flags)
            if win.last > -1:
                win.SetSelection(win.last, win.last+len(findTxt))
        else:# et == wxEVT_COMMAND_FIND_NEXT:
            if wxFR_DOWN&flags:
                if (win.last >= win.GetTextLength()) or (win.last == -1):
                    win.last = 0
                else:
                    win.last += len(findTxt)
                win.last = win.FindText(win.last, win.GetTextLength(), findTxt, flags)
            else:
                if (win.last >= win.GetTextLength()) or (win.last == -1):
                    win.last = win.GetTextLength()
                else:
                    win.last -= len(findTxt)
                win.last = win.FindText(win.last, 0, findTxt, flags)
            if win.last > -1:
                win.SetSelection(win.last, win.last+len(findTxt))
        win.EnsureVisible(win.GetCurrentLine())
        win.EnsureCaretVisible()
        if win.last == -1:
            self.dialog("Reached the end of the document.", "")

        return None


    def OnFindClose(self, evt):
        #self.log.write("wxFindReplaceDialog closing...\n")
        evt.GetDialog().Destroy()
        self.FINDSHOWN = 0
        for win in self.control:
            try:
                if win.last == -1:
                    #print win.gcp, "setting gcp"
                    win.SetSelection(win.gcp, win.gcp)
                    win.ScrollToColumn(0)
            except:
                continue

#---------------------------- View Menu Commands -----------------------------

    def Maximize(self, b):
        wxFrame.Maximize(b)
        wxPostEvent(wxSizeEvent((0,0), self.GetId()))

    def OnZoom(self, e):
        wnum, win = self.getNumWin(e)
        if e.GetId() == ZI:incr = 1
        else:              incr = -1
        win.SetZoom(win.GetZoom()+incr)

    def OnGoto(self, e):
        wnum, win = self.getNumWin(e)
        dlg = wxTextEntryDialog(self, '', 'Which line would you like to advance to?', '')
        resp = dlg.ShowModal()
        valu = dlg.GetValue()
        dlg.Destroy()
        if resp == wxID_OK:
            try:
                valu = int(valu)-1
            except:
                return
            if valu < win.GetLineCount():
                linepos = win.GetLineEndPosition(valu)
                win.EnsureVisible(valu)
                win.SetSelection(linepos-len(win.GetLine(valu))+len(win.format), linepos)
                win.ScrollToColumn(0)

    def OnTree(self, e):
        num, win = self.getNumWin(e)
        split = self.control.GetPage(num)
        width = self.control.GetClientSize()[0]-10
        if (width-split.GetSashPosition())<10:
            split.SetSashPosition(max(width-200, int(width/2)))
        else:
            split.SetSashPosition(width-split.GetMinimumPaneSize())
        win.SAVEDPOSITION = split.GetSashPosition()-self.control.GetClientSize()[0]

    def OnTreeHide(self, e):
        width = self.control.GetClientSize()[0]-10
        for stc in self.control:
            split = stc.parent
            split.SetSashPosition(width-split.GetMinimumPaneSize())
            stc.SAVEDPOSITION = split.GetSashPosition()-width

    def OnResize(self, e):
        num, stc = self.getNumWin(e)
        if num<0: return
        width = e.GetSize()[0]
        split = stc.parent
        split.SetSashPosition(max(width-stc.SAVEDPOSITION, int(width/2)))
        if e:
            e.Skip()

    def SetPos(self, e):
        try:    n, w = self.getNumWin(e)
        except: return
        sp = w.parent.GetSashPosition()
        width = self.control.GetClientSize()[0]
        if sp < 0:
            w.SAVEDPOSITION = abs(sp)
        else:
            w.SAVEDPOSITION = width-sp
        e.Skip()

    def startnext(self):
        self.parsing = 1
        done = 0
        while self.toparse:
            win = self.toparse.pop(0)
            if win in self.control:
                start = time.time()

                win.ConvertEOLs(fmt_mode[win.format])
                out = wxStyledTextCtrl.GetText(win).replace('\t', win.GetTabWidth()*' ')

                tpl = fast_parser(out, win.format, 3, wxYield)
                if not win in self.control:
                    continue
                h1, win.kw, win.tooltips, todo = tpl

                win.kw.sort()
                win.kw = ' '.join(win.kw)
                win.tree.new_hierarchy(h1)
                t = ' and tooltips'
                if win.todo:
                    win.todo.NewItemList(todo)
                    t = ', tooltips and todo'
                self.SetStatusText(("Browsable source tree, autocomplete%s"+
                                    " updated for %s in %.1f seconds.")%(t,
                                    win.filename, time.time()-start))
                win.refresh = 0
        self.parsing = 0

    def OnRefresh(self, e, win=None):
        if win is None:
            num, win = self.getNumWin(e)
        if win.refresh:
            return
        win.refresh = 1
        self.toparse.append(win)
        if not self.parsing:
            self.startnext()

    def OnToggleBookmark (self, e):
        wnum, win = self.getNumWin(e)
        lineNo = win.GetCurrentLine()
        if win.MarkerGet(lineNo) & BOOKMARKMASK:
            win.MarkerDelete(lineNo, BOOKMARKNUMBER)
        else:
            win.MarkerAdd(lineNo, BOOKMARKNUMBER)
      
    def OnNextBookmark  (self, e):
        wnum, win = self.getNumWin(e)
        lineNo = win.GetCurrentLine()
        newLineNo = win.MarkerNext(lineNo + 1, BOOKMARKMASK)
        if newLineNo != -1:
            win.GotoLine(newLineNo)
        else:
            lineNo = win.GetLineCount()
            newLineNo = win.MarkerNext(0, BOOKMARKMASK)
            if newLineNo != -1:
                win.GotoLine(newLineNo)
        win.EnsureVisible(win.GetCurrentLine())
        win.EnsureCaretVisible()
    
    def OnPreviousBookmark (self, e):
        wnum, win = self.getNumWin(e)
        lineNo = win.GetCurrentLine()
        newLineNo = win.MarkerPrevious(lineNo - 1, BOOKMARKMASK)
        if newLineNo != -1:
            win.GotoLine(newLineNo)
        else:
            lineNo = win.GetLineCount()
            newLineNo = win.MarkerPrevious(lineNo, BOOKMARKMASK)
            if newLineNo != -1:
                win.GotoLine(newLineNo)
        win.EnsureVisible(win.GetCurrentLine())
        win.EnsureCaretVisible()

#----------------------------- Set menu commands -----------------------------
    def OnSnipToggle(self, event):
        if event is None:
            self.config['usesnippets'] = 1
            self.restart = 0
        else:
            self.config['usesnippets'] = (self.config['usesnippets'] + 1) % 2
            self.restart = (self.restart + 1) % 2
            self.dialog("You have just %sabled code snippets.\n"\
                        "Restart is%s required for change to take effect."%\
                        (["dis", "en"][self.config['usesnippets']],
                         [" not", ''][self.restart]), "Preference change")

    def OnTodoToggle(self, event):
        if event is None:
            self.config['usetodo'] = 1
            self.restartt = 0
        else:
            self.config['usetodo'] = (self.config['usetodo'] + 1) % 2
            self.restartt = (self.restartt + 1) % 2
            self.dialog("You have just %sabled todos.\n"\
                        "Restart is%s required for change to take effect."%\
                        (["dis", "en"][self.config['usetodo']],
                         [" not", ''][self.restartt]), "Preference change")

    def OnChangeMenu(self, e):
        MenuItemDialog(self, self).ShowModal()

    def OnStyleChange(self,e):
        wnum, win = self.getNumWin(e)
        win.changeStyle(stylefile, lexers[e.GetId()])

    def OnDefaultStyleChange(self, e):
        dl = lexers[e.GetId()]
        self.config['DEFAULTLEXER'] = dl
        globals()['DEFAULTLEXER'] = dl

    def style(self, fn):
        ext = fn.split('.')[-1].lower()
        return extns.get(ext, DEFAULTLEXER)
    
    def OnEncChange(self, e):
        num, win = self.getNumWin(e)
        mid = e.GetId()
        newenc = self.menubar.GetLabel(mid)
        oldenc = win.enc
        if oldenc != newenc:
            win.enc = newenc
            self.SetStatusText("encoding changed to %s for %s"%(newenc, win.filename))
        self.SetStatusText(win.enc, 2)

    def OnLineEndChange(self, e):
        n, win = self.getNumWin(e)
        endid = e.GetId()
        newend = LE_MAPPING[endid]
        oldend = win.GetEOLMode()
        if oldend != newend:
            win.format = fmt_Rmode[newend]
            win.ConvertEOLs(newend)
            win.SetEOLMode(newend)

    def OnAutoCompleteToggle(self, event):
        # Images are specified with a appended "?type"
        #for i in range(len(kw)):
        #    if kw[i] in keyword.kwlist:
        #        kw[i] = kw[i]# + "?1"
        n, win = self.getNumWin(event)
        win.showautocomp = (win.showautocomp+1)%2

    def OnNumberToggle(self, e):
        n, win = self.getNumWin(e)
        win.SetMarginWidth(0, (win.GetMarginWidth(0)+40)%80)
        
    def OnMarginToggle(self, e):
        n, win = self.getNumWin(e)
        win.SetMarginWidth(1, (win.GetMarginWidth(1)+16)%32)

    def OnIndentGuideToggle(self, e):
        n, win = self.getNumWin(e)
        win.SetIndentationGuides((win.GetIndentationGuides()+1)%2)
    def OnSavePositionToggle(self, e):
        n, win = self.getNumWin(e)
        win.save_cursor = (1+win.save_cursor)%2

    def OnExpandAll(self, e):
        n, win = self.getNumWin(e)
        lc = win.GetLineCount()
        win.ShowLines(0, lc-1)
        for line in xrange(lc):
            if win.GetFoldLevel(line) & wxSTC_FOLDLEVELHEADERFLAG:
                win.SetFoldExpanded(line, 1)

    def OnFoldAll(self, e):
        n, win = self.getNumWin(e)
        #these next two lines are to allow
        #win.GetFoldLevel() to be accurate
        win.HideLines(0, win.GetLineCount()-1)
        try: wxYield()
        except: pass
        
        #toss all the old folds
        self.OnExpandAll(e)
        
        lc = win.GetLineCount()
        lines = []
        for line in xrange(lc):
            if win.GetFoldLevel(line) & wxSTC_FOLDLEVELHEADERFLAG:
                lines.append(line)
        lines.reverse()
        for line in lines:
            a = win.GetLastChild(line, -1)
            win.HideLines(line+1,a)
            win.SetFoldExpanded(line, 0)

    def OnSetTabToggle(self, e):
        n, win = self.getNumWin(e)
        a = win.GetUseTabs()
        win.SetUseTabs((a+1)%2)
        win.SetProperty("tab.timmy.whinge.level", str(a))

    def OnWrapL(self, e):
        wnum, win = self.getNumWin(e)
        self.WrapToggle(win)

    def WrapToggle(self, win):
        #tossing the current home/end key configuration, as it will change
        win.CmdKeyClear(wxSTC_KEY_HOME,0)
        win.CmdKeyClear(wxSTC_KEY_HOME,wxSTC_SCMOD_SHIFT)
        win.CmdKeyClear(wxSTC_KEY_HOME,wxSTC_SCMOD_ALT)
        win.CmdKeyClear(wxSTC_KEY_HOME,wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT)
        win.CmdKeyClear(wxSTC_KEY_END,0)
        win.CmdKeyClear(wxSTC_KEY_END,wxSTC_SCMOD_SHIFT)
        win.CmdKeyClear(wxSTC_KEY_END,wxSTC_SCMOD_ALT)
        win.CmdKeyClear(wxSTC_KEY_END,wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT)

        if win.GetWrapMode() == wxSTC_WRAP_NONE:
            win.SetWrapMode(wxSTC_WRAP_WORD)
            
            #making home and end work like we expect them to when lines are wrapped
            win.CmdKeyAssign(wxSTC_KEY_HOME, 0, wxSTC_CMD_HOMEDISPLAY)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_SHIFT, wxSTC_CMD_HOMEDISPLAYEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_ALT, wxSTC_CMD_VCHOME)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT, wxSTC_CMD_VCHOMEEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_END, 0, wxSTC_CMD_LINEENDDISPLAY)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_SHIFT, wxSTC_CMD_LINEENDDISPLAYEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_ALT, wxSTC_CMD_LINEEND)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT, wxSTC_CMD_LINEENDEXTEND)
        else:
            win.SetWrapMode(wxSTC_WRAP_NONE)

            #making home and end work like we expect them to when lines are not wrapped
            win.CmdKeyAssign(wxSTC_KEY_HOME, 0, wxSTC_CMD_VCHOME)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_SHIFT, wxSTC_CMD_VCHOMEEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_ALT, wxSTC_CMD_HOMEDISPLAY)
            win.CmdKeyAssign(wxSTC_KEY_HOME, wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT, wxSTC_CMD_HOMEDISPLAYEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_END, 0, wxSTC_CMD_LINEEND)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_SHIFT, wxSTC_CMD_LINEENDEXTEND)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_ALT, wxSTC_CMD_LINEENDDISPLAY)
            win.CmdKeyAssign(wxSTC_KEY_END, wxSTC_SCMOD_ALT|wxSTC_SCMOD_SHIFT, wxSTC_CMD_LINEENDDISPLAYEXTEND)

    def OnTreeSortToggle(self, e):
        n, win = self.getNumWin(e)
        win.tree.tree.SORTTREE = (win.tree.tree.SORTTREE+1)%2
        self.OnRefresh(e, win)

    def OnSetIndent(self, e):
        n, win = self.getNumWin(e)
        dlg = wxTextEntryDialog(self, "Enter an integer > 0.",
                                "How many spaces per indent level?")
        orig = win.GetIndent()
        dlg.SetValue(str(orig))
        rslt = dlg.ShowModal()
        if rslt == wxID_OK:
            win.SetIndent(validate(dlg, orig))
        dlg.Destroy()

    def OnSetTabWidth(self, e):
        n, win = self.getNumWin(e)
        dlg = wxTextEntryDialog(self,\
                                "How many spaces should a tab character,\n"\
                                "'\\t' represent?  Enter an integer > 0.",\
                                "How many spaces per tab?")
        orig = win.GetTabWidth()
        dlg.SetValue(str(orig))
        rslt = dlg.ShowModal()
        if rslt == wxID_OK:
            win.SetTabWidth(validate(dlg, orig))
        dlg.Destroy()

    def OnSetLongLinePosition(self, e):
        n, win = self.getNumWin(e)
        dlg = wxTextEntryDialog(self,
                                "At what column would you like a long line\n"\
                                "signifier to be displayed?  Enter an integer > 0.",
                                "Long Line Indicator")
        orig = win.GetEdgeColumn()
        dlg.SetValue(str(orig))
        rslt = dlg.ShowModal()
        if rslt == wxID_OK:
            win.SetEdgeColumn(validate(dlg, orig))
        dlg.Destroy()

    def OnSetLongLineMode(self, e):
        n, win = self.getNumWin(e)
        eid = e.GetId()
        win.SetEdgeMode(LL_MAPPING[eid])

    def OnSaveSettings(self, e):
        n, win = self.getNumWin()
        dct = win.GetSaveState()
        del dct['BM']
        del dct['FOLD']

        globals().update(dct)
        self.config['DOCUMENT_DEFAULTS'] = dct

#----------------------------- Tab menu commands -----------------------------
    def OnLeft(self, e):
        self.control.AdvanceSelection(False)
    def OnRight(self, e):
        self.control.AdvanceSelection(True)
    def MoveLeft(self, e):
        wnum = self.control.GetSelection()
        pagecount = self.control.GetPageCount()
        if (wnum >= 0) and (pagecount>1):
            l, r = (wnum-1)%pagecount, wnum
            self.control.swapPages(l, r, 1)
            self.control.SetSelection(l)
        else:
            e.Skip()
    def MoveRight(self, e):
        wnum = self.control.GetSelection()
        pagecount = self.control.GetPageCount()
        if (wnum >= 0) and (pagecount>1):
            l, r = (wnum+1)%pagecount, wnum
            self.control.swapPages(l, r, 0)
            self.control.SetSelection(l)
        else:
            e.Skip()

#----------------------------- Snippet commands ------------------------------
    def OnSnippet(self, e, resize=0):
        split = self.split
        num = self.control.GetSelection()
        if num != -1:
            resize = 1
            width = self.control.GetClientSize()[0]
            orig = self.control.GetPage(num).GetSashPosition()
        if split.GetSashPosition() < 10:
            split.SetSashPosition(100)
        else:
            split.SetSashPosition(split.GetMinimumPaneSize())
        if resize:
            delta = self.control.GetClientSize()[0]-width
            self.control.GetPage(num).SetSashPosition(orig+delta)

#------------------------- Bookmarked Path Commands --------------------------
    def OnPathmark(self, e, st=type('')):
        if 'lastpath' in self.config:
            del self.config['lastpath']
        try:
            key = e.GetKeyCode()
            pth = pathmarks.get(key)
        except:
            pth = pathmarks[bmId2Keypress[e.GetId()]]
        if (type(pth) is st) and os.path.isdir(pth):
            os.chdir(pth)
            self.SetStatusText("Changed path to %s"%pth)

    def ViewPathmarks(self, e, titl="Pathmarks", styl=wxOK, sel=0, st=type('')):
        kys = pathmarks.keys()
        kys.sort()
        out = []
        for i in kys:
            if type(pathmarks[i]) is st:
                out.append("%i            %s    "%(i-48, pathmarks[i]))
            else:
                out.append("%i            <unused>    "%(i-48))

        dlg = wxSingleChoiceDialog(self, "Pathmark       Path",titl, out, styl|wxDD_DEFAULT_STYLE)
        dlg.SetSelection(sel)
        if dlg.ShowModal() == wxID_OK:
            retr = dlg.GetStringSelection()
        else:
            retr = None
        dlg.Destroy()
        return retr
        
    def AddPathmark(self, e):
        pth = self.ViewPathmarks(e, "Select a pathmark to update", wxOK|wxCANCEL)
        while pth != None:
            posn = int(pth[0])-1
            opth = pathmarks[posn+1+48]
            if opth == 0: opth = os.getcwd()
            dlg = wxDirDialog(self, "Choose a path", opth, style=wxDD_DEFAULT_STYLE|wxDD_NEW_DIR_BUTTON)
            if dlg.ShowModal() == wxID_OK:
                pmn = posn+1+48
                path = dlg.GetPath()
                self.dpm = 1
                self.remPos(pmn)
                self.addPos(pmn, path)
            else:
                self.SetStatusText("Pathmark modification aborted")
            pth = self.ViewPathmarks(e, "Select a pathmark to update", wxOK|wxCANCEL, posn)
        self.SetStatusText("Finished updating pathmarks")

    def itemPos(self, pmn):
        cnt = 0
        for i in xrange(49, pmn):
            if pathmarks[i] != 0:
                cnt += 1
        return cnt

    def remPos(self, pmn):
        if pathmarks[pmn] != 0:
            posn = self.itemPos(pmn)
            a = pathmarks[pmn]
            self.pm.Delete(self.pm.GetMenuItems()[posn+4].GetId())
            pathmarks[pmn] = 0
            self.SetStatusText("Deleted %s from bookmark %i"%(a, pmn-48))

    def addPos(self, pmn, path):
        posn = self.itemPos(pmn)
        self.pm.Insert(posn+4, bmPm2Id[pmn], "Ctrl+%i\t%s"%(pmn-48, path), "Change the current working directory to %s"%path)
        pathmarks[pmn] = path
        self.SetStatusText("Set bookmark %i to %s"%(pmn-48, path))

    def RemovePathmark(self, e):
        pth = self.ViewPathmarks(e, "Delete pathmark (cancel to finish)", wxOK|wxCANCEL)
        while pth != None:
            self.dpm = 1
            posn = int(pth[0])-1
            pmn = posn+1+48
            if pathmarks[pmn] != 0:
                self.remPos(pmn)
            pth = self.ViewPathmarks(e, "Delete pathmark (cancel to finish)", wxOK|wxCANCEL, posn)

#---------------------------- Help Menu Commands -----------------------------
    def OnAbout(self, e):
        txt = """
        You're wondering what this editor is all about, right?  Easy, this edior was
        written to scratch an itch.  I (Josiah Carlson), was looking for an editor
        that had the features I wanted, I couldn't find one.  So I wrote one.  And
        here it is.
        
        PyPE %s (Python Programmers Editor)
        http://come.to/josiah
        PyPE is copyright 2003 Josiah Carlson.
        
        This software is licensed under the GPL (GNU General Public License) as it
        appears here: http://www.gnu.org/copyleft/gpl.html  It is also included with
        this software as gpl.txt.
        
        If you do not also receive a copy of gpl.txt with your version of this
        software, please inform the me of the violation at the web page near the top
        of this document."""%VERSION
        self.dialog(txt.replace('        ', ''), "About...")
    def OnHelp(self, e):
        a = open(os.path.join(runpath, 'readme.txt'), 'rb')
        txt = a.read()
        a.close()
        dlg = wxScrolledMessageDialog(self, txt, "Help!")
        dlg.ShowModal()
        dlg.Destroy()
#-------------------------- Hot Key support madness --------------------------
    def OnKeyPressed(self, event):
        showpress=0
        
        keypressed = GetKeyPress(event)
        
        if showpress: print "keypressed", keypressed
        key = event.KeyCode()
        wnum = self.control.GetSelection()
        pagecount = self.control.GetPageCount()
        if wnum > -1:
            win = self.control.GetPage(wnum).GetWindow1()
        
        if keypressed in HOTKEY_TO_ID:
            menuid = HOTKEY_TO_ID[keypressed]
            wxPostEvent(self, wxMenuEvent(wxEVT_COMMAND_MENU_SELECTED, menuid))

        else:
            if pagecount:
                if (key==13):
                    #when 'enter' is pressed, indentation needs to happen.
                    #
                    #will indent the current line to be equivalent to the line above
                    #unless a ':' is at the end of the previous, then will indent
                    #configuration.indent more.
                    if win.AutoCompActive():
                        return win.AutoCompComplete()
                    if win.CallTipActive():
                        win.CallTipCancel()
                    #get information about the current cursor position
                    linenum = win.GetCurrentLine()
                    pos = win.GetCurrentPos()
                    col = win.GetColumn(pos)
                    linestart = win.PositionFromLine(linenum)
                    line = win.GetLine(linenum)[:pos-linestart]
                    
                    #get info about the current line's indentation
                    ind = win.GetLineIndentation(linenum)
                    
                    colon = ord(':')
                    if col <= ind:
                        if win.GetUseTabs():
                            win.ReplaceSelection(win.format+(col*' ').replace(win.GetTabWidth()*' ', '\t'))
                        else:
                            win.ReplaceSelection(win.format+(col*' '))
                    elif pos:
                        xtra = 0
                        if (line.find(':')>-1):
                            for i in xrange(linestart, min(pos, win.GetTextLength())):
                                styl = win.GetStyleAt(i)
                                #print styl, win.GetCharAt(i)
                                if not xtra:
                                    if (styl==10) and (win.GetCharAt(i) == colon):
                                        xtra = 1
                                elif (styl == 1):
                                    #it is a comment, ignore the character
                                    pass
                                elif (styl == 0) and (win.GetCharAt(i) in [ord(i) for i in ' \t\r\n']):
                                    #is not a comment, but is the space before a comment
                                    #or the end of a line, ignore the character
                                    pass
                                else:
                                    #this is not a comment or some innocuous other character
                                    #this is a docstring or otherwise, no additional indent
                                    xtra = 0
                                    #commenting the break fixes stuff like this
                                    # for i in blah[:]:
                                    #break
                            if xtra:
                                #This deals with ending single and multi-line definitions properly.
                                while linenum >= 0:
                                    found = []
                                    for i in ['def', 'class', 'if', 'else', 'elif', 'while',
                                              'for', 'try', 'except', 'finally']:
                                        a = line.find(i)
                                        if (a > -1):
                                            found.append(a)
                                    #print 'fnd', found
                                    if found: found = min(found)
                                    else:     found = -1
                                    if (found > -1) and\
                                       (win.GetStyleAt(win.GetLineEndPosition(linenum)-len(line)+found)==5) and\
                                       (win.GetLineIndentation(linenum) == found):
                                           ind = win.GetLineIndentation(linenum)
                                           break
                                    linenum -= 1
                                    line = win.GetLine(linenum)
                        #if we were to do indentation for ()[]{}, it would be here
                        if not xtra:
                            #yep, right here.
                            fnd = 0
                            for i in "(){}[]":
                                if (line.find(i) > -1):
                                    fnd = 1
                                    break
                            if fnd:
                                seq = []
                                #print "finding stuff"
                                for i in "(){}[]":
                                    a = line.find(i)
                                    start = 0
                                    while a > -1:
                                        start += a+1
                                        if win.GetStyleAt(start+linestart-1)==10:
                                            seq.append((start, i))
                                        a = line[start:].find(i)
                                seq.sort()
                                cl = {')':'(', ']': '[', '}': '{',
                                      '(':'',  '[': '',  '{': ''}
                                stk = []
                                #print "making tree"
                                for po, ch in seq:
                                    #print ch,
                                    if not cl[ch]:
                                        #standard opening
                                        stk.append((po, ch))
                                    elif stk:
                                        if cl[ch] == stk[-1][1]:
                                            #proper closing of something
                                            stk.pop()
                                        else:
                                            #probably a syntax error
                                            #does it matter what is done?
                                            stk = []
                                            break
                                    else:
            #Probably closing something on another line, should probably find
            #the indent level of the opening, but that looks to be a hard
            #problem, and multiple lines would need to be checked for the
            #location of the closing item.
            #This would require an incremental line-by-line stepping back.
            #Almost like what was done for ['def', 'class', 'if', 'while', 'for']
            #and ':'.  However, there is no guarantee that there is a proper
            #start, and we could end up parsing the entire file...a few times
            #over.  Incremental backwards things end up doing worst case
            #O(n^2) line parsings.
            #Fuck that.
            #People can use the improved single-line dedent with crtl+[.
                                        stk = []
                                        break
                                if stk:
                                    #print "stack remaining", stk
                                    ind = stk[-1][0]
                        if not xtra:
                            ls = line.lstrip()
                            if (ls[:6] == 'return') or (ls[:4] == 'pass'):
                                xtra = -1
                        a = ind*' '
                        if win.GetUseTabs():
                            a = a.replace(win.GetTabWidth()*' ', '\t')
                        win.ReplaceSelection(win.format+a)
                        if xtra == 1:
                            self.OnIndent(event)
                        elif xtra == -1:
                            self.OnDedent(event)
                    else:
                        win.ReplaceSelection(win.format)
                else:
                    event.Skip()
                    if not (win.GetStyleAt(win.GetCurrentPos())):
                        #3, 13, 6, 7
                        if win.CallTipActive():
                            good = STRINGPRINTABLE
                            if (key in good) or (key in (WXK_SHIFT, WXK_CONTROL, WXK_ALT)):
                                pass
                                #if key in (48, 57) and event.ShiftDown():
                                #    win.CallTipCancel()
                                #else it is something in the arguments that is OK.
                            else:
                                win.CallTipCancel()
                        if (not win.CallTipActive()) and event.ShiftDown() and (key == ord('9')):
                            win.CallTipSetBackground(wxColour(255, 255, 232))
                            cur, colpos, word = self.getLeftFunct(win)
                            tip = '\n'.join(win.tooltips.get(word, []))
                            if tip:
                                win.CallTipShow(win.GetCurrentPos(),tip)
                        elif win.showautocomp and bool(win.kw) and (not win.AutoCompActive()) and (ord('A') <= key) and (key <= ord('Z')):
                            #if win.CallTipActive():
                            #    win.CallTipCancel()
                            cur, colpos, word = self.getLeftFunct(win)
                            indx = win.kw.find(word)
                            if (not word) or ((indx > -1) and ((win.kw[indx-1] == ' ') or (indx==0))):
                                win.AutoCompSetIgnoreCase(False)
                                win.AutoCompShow(colpos-cur, win.kw)
            else:
                return event.Skip()

    def getLeftFunct(self, win):
        t = ' .,;:([)]}\'"\\<>%^&+-=*/|`'
        bad = dict(zip(t, [0]*len(t)))
        line = win.GetLine(win.GetCurrentLine())
        colpos = win.GetColumn(win.GetCurrentPos())
        cur = colpos-1
        while (cur >= 0) and not (line[cur:cur+1] in bad):
            cur -= 1
        cur += 1
        return cur, colpos, line[cur:colpos]

#------------- Ahh, Styled Text Control, you make this possible. -------------
class PythonSTC(wxStyledTextCtrl):
    def __init__(self, notebook, ID, parent):
        wxStyledTextCtrl.__init__(self, parent, ID, style = wxNO_FULL_REPAINT_ON_RESIZE)

        self.MarkerDefine(BOOKMARKNUMBER, BOOKMARKSYMBOL, 'blue', 'blue')

        self.hierarchy = []
        self.kw = []
        self.tooltips = {}

        self.parent = parent
        self.notebook = notebook #should also be equal to self.parent.parent
        self.root = self.notebook.root
        self.dirty = 0
        self.refresh = 0
        
        #Text is included in the original, but who drags text?  Below for dnd file support.
        if dnd_file: self.SetDropTarget(FileDropTarget(self.root))

        #for command comlpetion
        self.SetKeyWords(0, " ".join(keyword.kwlist))
        
        self.SetMargins(0,0)
        self.SetViewWhiteSpace(False)
        self.SetBackSpaceUnIndents(1)
        #self.SetBufferedDraw(False)
        #self.SetViewEOL(True)

        self.SetEdgeMode(col_mode)
        self.SetEdgeColumn(col_line)
        

        self.SetMarginType(0, wxSTC_MARGIN_NUMBER)
        if line_margin:
            self.SetMarginWidth(0, 40)
        else:
            self.SetMarginWidth(0, 0)
        #self.StyleSetSpec(wxSTC_STYLE_LINENUMBER, "size:%(size)d,face:%(mono)s" % faces)

        # Setup a margin to hold fold markers
        #I agree, what is this value?
        #self.SetFoldFlags(16)  ###  WHAT IS THIS VALUE?  WHAT ARE THE OTHER FLAGS?  DOES IT MATTER?
        if marker_margin:
            self.SetMarginWidth(1, 16)
        else:
            self.SetMarginWidth(1, 0)
        
        self.SetProperty("fold", "1")
        self.SetMarginType(2, wxSTC_MARGIN_SYMBOL)
        self.SetMarginMask(2, wxSTC_MASK_FOLDERS)
        self.SetMarginSensitive(2, True)
        self.SetMarginWidth(2, 12)

        if collapse_style: # simple folder marks, like the old version
            self.MarkerDefine(wxSTC_MARKNUM_FOLDER, wxSTC_MARK_BOXPLUS, "navy", "white")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDEROPEN, wxSTC_MARK_BOXMINUS, "navy", "white")
            # Set these to an invisible mark
            self.MarkerDefine(wxSTC_MARKNUM_FOLDEROPENMID, wxSTC_MARK_BACKGROUND, "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERMIDTAIL, wxSTC_MARK_BACKGROUND, "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERSUB, wxSTC_MARK_BACKGROUND, "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERTAIL, wxSTC_MARK_BACKGROUND, "white", "black")

        else: # more involved "outlining" folder marks
            self.MarkerDefine(wxSTC_MARKNUM_FOLDEREND,     wxSTC_MARK_BOXPLUSCONNECTED,  "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDEROPENMID, wxSTC_MARK_BOXMINUSCONNECTED, "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERMIDTAIL, wxSTC_MARK_TCORNER,  "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERTAIL,    wxSTC_MARK_LCORNER,  "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDERSUB,     wxSTC_MARK_VLINE,    "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDER,        wxSTC_MARK_BOXPLUS,  "white", "black")
            self.MarkerDefine(wxSTC_MARKNUM_FOLDEROPEN,    wxSTC_MARK_BOXMINUS, "white", "black")
        

        # Make some styles,  The lexer defines what each style is used for, we
        # just have to define what each style looks like.  This set is adapted from
        # Scintilla sample property files.
        
        #And good wxPython users who have seen some demos know the above was copied
        #and pasted, along with alot of sample code, right out of the demo.  The
        #demo r0XX0rs.
        
        # Global default styles for all languages
        self.StyleSetSpec(wxSTC_STYLE_DEFAULT,     "fore:#000000,face:%(mono)s,back:#FFFFFF,size:%(size)d" % faces)
        self.StyleSetSpec(wxSTC_STYLE_LINENUMBER,  "back:#C0C0C0,face:Lucida Console,size:%(size2)d" % faces)
        self.StyleSetSpec(wxSTC_STYLE_CONTROLCHAR, "face:%(other)s" % faces)
        self.StyleSetSpec(wxSTC_STYLE_BRACELIGHT,  "fore:#003000,face:%(mono)s,back:#80E0E0"% faces)
        self.StyleSetSpec(wxSTC_STYLE_BRACEBAD,    "fore:#E0FFE0,face:%(mono)s,back:#FF0000"% faces)

        #various settings
        if use_tabs:
            self.SetProperty("tab.timmy.whinge.level", "0")
        else:
            self.SetProperty("tab.timmy.whinge.level", "1")
        self.SetSelBackground(1, '#B0B0FF')
        self.SetIndent(indent)
        self.SetUseTabs(use_tabs)
        self.SetTabWidth(spaces_per_tab)

        #again, some state variables
        self.filename = ''
        self.dirname = ''
        self.opened = 0
        self.AutoCompStops(' .,;:()[]{}\'"\\<>%^&+-=*/|`')

        EVT_STC_UPDATEUI(self,    ID, self.OnUpdateUI)
        EVT_STC_MARGINCLICK(self, ID, self.OnMarginClick)
        EVT_KEY_DOWN(self, self.root.OnKeyPressed)
        #EVT_CHAR(self, self.root.OnKeyPressed)
        EVT_KEY_UP(self, self.key_up)
        
        EVT_STC_CHARADDED(self, ID, self.key_up)
        EVT_STC_CHANGE(self, ID, self.cha)
        EVT_STC_POSCHANGED(self, ID, self.pos)
        EVT_STC_SAVEPOINTREACHED(self, ID, self.MakeClean)
        EVT_STC_SAVEPOINTLEFT(self, ID, self.MakeDirty)
        self.SetModEventMask(wxSTC_MOD_INSERTTEXT|wxSTC_MOD_DELETETEXT|wxSTC_PERFORMED_USER|wxSTC_PERFORMED_UNDO|wxSTC_PERFORMED_REDO)

        if REM_SWAP:
            self.CmdKeyClear(ord('T'), wxSTC_SCMOD_CTRL)
        if wrapmode != wxSTC_WRAP_NONE:
            self.root.WrapToggle(self)
        
        self.CmdKeyClear(ord('Z'), wxSTC_SCMOD_CTRL)
        self.CmdKeyClear(ord('Y'), wxSTC_SCMOD_CTRL)
        self.CmdKeyClear(ord('X'), wxSTC_SCMOD_CTRL)
        self.CmdKeyClear(ord('C'), wxSTC_SCMOD_CTRL)
        self.CmdKeyClear(ord('V'), wxSTC_SCMOD_CTRL)
        self.CmdKeyClear(ord('A'), wxSTC_SCMOD_CTRL)
        
        
    def pos(self, e):
        #print "position changed event"
        #for some reason I never get called
        self.pos_ch(e)

    def cha(self, e):
        #print "changed event"
        self.key_up(e)
    
    def save(self, e):
        #print "save reached"
        self.key_up(e)
    
    def savel(self, e):
        #print "save left"
        self.key_up(e)
        
    def key_up(self, e):
        if self.GetModify(): self.MakeDirty()
        else:                self.MakeClean()
        self.pos_ch(e)
        #e.Skip()

    def pos_ch(self, e):
        lin = self.GetCurrentLine()+1
        col = self.GetColumn(self.GetCurrentPos())
        self.root.SetStatusText("L%i C%i"%(lin, col), 1)
        if e:
            e.Skip()

#------------------------- persistant document state -------------------------
    def GetSaveState(self):
        BM = []
        FOLD = []
        ret =  {'BM':BM,
                'FOLD':FOLD,
                'use_tabs':self.GetUseTabs(),
                'spaces_per_tab':self.GetTabWidth(),
                'indent':self.GetIndent(),
                'marker_margin':bool(self.GetMarginWidth(1)),
                'line_margin':bool(self.GetMarginWidth(0)),
                'col_line':self.GetEdgeColumn(),
                'col_mode':self.GetEdgeMode(),
                'indent_guide':self.GetIndentationGuides(),
                'showautocomp':self.showautocomp,
                'wrapmode':self.GetWrapMode(),
                'sortmode':self.tree.tree.SORTTREE,
                'SAVEDPOSITION':self.SAVEDPOSITION,
                'save_cursor':self.save_cursor,
                'cursor_posn':self.GetCurrentPos()
               }
        for line in xrange(self.GetLineCount()):
            if self.MarkerGet(line) & BOOKMARKMASK:
                BM.append(line)
            
            if (self.GetFoldLevel(line) & wxSTC_FOLDLEVELHEADERFLAG) and\
            (not self.GetFoldExpanded(line)):
                FOLD.append(line)
        FOLD.reverse()
        return ret

    def SetSaveState(self, saved):
        self.SetUseTabs(saved.get('use_tabs', use_tabs))
        self.SetTabWidth(saved.get('space_per_tab', spaces_per_tab))
        self.SetIndent(saved.get('indent', indent))
        self.SetMarginWidth(0, 40*saved.get('line_margin', line_margin))
        self.SetMarginWidth(1, 16*saved.get('marker_margin', marker_margin))
        self.SetEdgeColumn(saved.get('col_line', col_line))
        self.SetEdgeMode(saved.get('col_mode', col_mode))
        self.SetIndentationGuides(saved.get('indent_guide', indent_guide))
        self.showautocomp = saved.get('showautocomp', showautocomp)
        if self.GetWrapMode() != saved.get('wrapmode', wrapmode):
            self.root.WrapToggle(self)
        self.tree.tree.SORTTREE = saved.get('sortmode', sortmode)
        self.SAVEDPOSITION = saved.get('SAVEDPOSITION', SAVEDPOSITION)
        self.save_cursor = saved.get('save_cursor', save_cursor)
        
        for bml in saved.get('BM', []):
            self.MarkerAdd(bml, BOOKMARKNUMBER)
        
        for exl in saved.get('FOLD', []):
            a = self.GetLastChild(exl, -1)
            self.HideLines(exl+1,a)
            self.SetFoldExpanded(exl, 0)
        
        try: wxYield()
        except: pass
        
        if self.save_cursor:
            a = saved.get('cursor_posn', 0)
            self.SetSelection(a,a)
            self.EnsureCaretVisible()
            self.ScrollToColumn(0)

#-------------------- fix for SetText for the 'dirty bit' --------------------
    def SetText(self, txt, emptyundo=1):
        self.SetEOLMode(fmt_mode[self.format])
        self.enc = 'ascii'
        if VS[-1] == 'u':
            for bom, enc in BOM:
                if txt[:len(bom)] == bom:
                    self.enc = enc
                    txt = txt[len(bom):]
                    #print "chose", enc
                    break
            if self.enc != 'ascii':
                try:    txt = txt.decode(self.enc)
                except:
                    #print "failed text decoding"
                    self.root.dialog("There has been a unicode decoding error."
                                     "The cause of this error is unknown to PyPE."
                                     "To prevent loss or corruption of data, it"
                                     "is suggested that you close the document,"
                                     "do not save.  Then try to open the document"
                                     "with the application you originally created"
                                     "it with.  If PyPE was the original creator,"
                                     "and only editor of the document, please"
                                     "contact the author and submit a bug report.",
                                     "Unicode decoding error.")
                    self.enc = ascii
        wxStyledTextCtrl.SetText(self, txt)
        self.ConvertEOLs(fmt_mode[self.format])
        self.opened = 1
        if emptyundo:
            self.EmptyUndoBuffer()
            self.SetSavePoint()

    def GetText(self):
        self.ConvertEOLs(fmt_mode[self.format])
        if VS[-1] == 'u':
            if self.enc == 'ascii':
                try:
                    return wxStyledTextCtrl.GetText(self).encode(self.enc)
                except:
                    #Previously non-unicode ascii file has had unicode characters
                    #inserted.  Must encode into some sort of unicode format.
                    self.enc = 'utf-8'
                    self.root.SetStatusText(self.enc, 2)
            return ADDBOM[self.enc] + wxStyledTextCtrl.GetText(self).encode(self.enc)
        return wxStyledTextCtrl.GetText(self)

#----- Takes care of the little '*' modified next to the open file name ------
    def MakeDirty(self, e=None):
        if (not self.dirty) and self.opened:
            self.dirty = 1
            f = self.filename
            if f == ' ':
                f = '<untitled %i>'%self.NEWDOCUMENT
            c = 0
            for i in self.notebook:
                if i == self:
                    break
                c += 1
            self.notebook.SetPageText(c, '* '+f)
            self.root.redrawvisible(self)
        if e:
            e.Skip()

    def MakeClean(self, e=None):
        if self.dirty:
            self.dirty = 0
            f = self.filename
            if f == ' ':
                f = '<untitled %i>'%self.NEWDOCUMENT
            c = 0
            for i in self.notebook:
                if i == self:
                    break
                c += 1
            self.notebook.SetPageText(c, f)
            self.SetSavePoint()
            self.root.redrawvisible(self)
        if e:
            e.Skip()


    def nada(self, e):
        pass

    def do(self, funct, dirty=1):
        if dirty: self.MakeDirty(None)
        funct(self)
        self.root.redrawvisible(self)

    def Cut(self):      self.do(wxStyledTextCtrl.Cut)
    def Paste(self):    self.do(wxStyledTextCtrl.Paste)
    def Undo(self):     self.do(wxStyledTextCtrl.Undo)
    def Redo(self):     self.do(wxStyledTextCtrl.Redo)

#--------- Ahh, the style change code...isn't it great?  Not really. ---------
    def changeStyle(self, stylefile, language):
        try:
            #from StyleSupport import initSTC
            from STCStyleEditor import initSTC
            initSTC(self, stylefile, language)

        except:
            #self.root.exceptDialog("Style Change failed, assuming plain text")
            self.root.SetStatusText("Style Change failed, assuming plain text")
            
#----------------- Defaults, in case the other code was bad. -----------------
            #for some default font styles

            self.SetLexer(wxSTC_LEX_NULL)
    
            ### Python styles
            ##
            ### White space
            ##self.StyleSetSpec(wxSTC_P_DEFAULT, "fore:#000000,face:%(helv)s,size:%(size)d" % faces)
            ### Comment
            ##self.StyleSetSpec(wxSTC_P_COMMENTLINE, "fore:#007F00,face:%(mono)s,back:#E0FFE0,size:%(size)d" % faces)
            ### Number
            ##self.StyleSetSpec(wxSTC_P_NUMBER, "fore:#007F7F,face:%(times)s,size:%(size)d" % faces)
            ### String
            ##self.StyleSetSpec(wxSTC_P_STRING, "fore:#7F007F,face:%(times)s,size:%(size)d" % faces)
            ### Single quoted string
            ##self.StyleSetSpec(wxSTC_P_CHARACTER, "fore:#7F007F,face:%(times)s,size:%(size)d" % faces)
            ### Keyword
            ##self.StyleSetSpec(wxSTC_P_WORD, "fore:#F0B000,face:%(mono)s,size:%(size)d,bold" % faces)
            ### Triple quotes
            ##self.StyleSetSpec(wxSTC_P_TRIPLE, "fore:#603000,face:%(times)s,back:#FFFFE0,size:%(size)d" % faces)
            ### Triple double quotes
            ##self.StyleSetSpec(wxSTC_P_TRIPLEDOUBLE, "fore:#603000,face:%(times)s,back:#FFFFE0,size:%(size)d" % faces)
            ### Class name definition
            ##self.StyleSetSpec(wxSTC_P_CLASSNAME, "fore:#0000FF,face:%(times)s,size:%(size)d,bold" % faces)
            ### Function or method name definition
            ##self.StyleSetSpec(wxSTC_P_DEFNAME, "fore:#0000FF,face:%(times)s,size:%(size)d,bold" % faces)
            ### Operators
            ##self.StyleSetSpec(wxSTC_P_OPERATOR, "face:%(times)s,size:%(size)d" % faces)
            ### Identifiers
            ##self.StyleSetSpec(wxSTC_P_IDENTIFIER, "fore:#000000,face:%(helv)s,size:%(size)d" % faces)
            ### Comment-blocks
            ##self.StyleSetSpec(wxSTC_P_COMMENTBLOCK, "fore:#7F7F7F,face:%(times)s,size:%(size)d" % faces)
            ### End of line where string is not closed
            ##self.StyleSetSpec(wxSTC_P_STRINGEOL, "fore:#000000,face:%(mono)s,back:#E0C0E0,eol,size:%(size)d" % faces)
            ##
            ##self.SetCaretForeground("BLACK")
            ##
            ### prototype for registering some images for use in the AutoComplete box.
            ###self.RegisterImage(1, images.getSmilesBitmap())

#------------ copied and pasted from the wxStyledTextCtrl_2 demo -------------
    def OnUpdateUI(self, evt):
        # check for matching braces
        braceAtCaret = -1
        braceOpposite = -1
        charBefore = None
        caretPos = self.GetCurrentPos()
        if caretPos > 0:
            charBefore = self.GetCharAt(caretPos - 1)
            styleBefore = self.GetStyleAt(caretPos - 1)

        # check before
        if charBefore and chr(charBefore) in "[]{}()" and styleBefore == wxSTC_P_OPERATOR:
            braceAtCaret = caretPos - 1

        # check after
        if braceAtCaret < 0:
            charAfter = self.GetCharAt(caretPos)
            styleAfter = self.GetStyleAt(caretPos)
            if charAfter and chr(charAfter) in "[]{}()" and styleAfter == wxSTC_P_OPERATOR:
                braceAtCaret = caretPos

        if braceAtCaret >= 0:
            braceOpposite = self.BraceMatch(braceAtCaret)

        if braceAtCaret != -1  and braceOpposite == -1:
            self.BraceBadLight(braceAtCaret)
        else:
            self.BraceHighlight(braceAtCaret, braceOpposite)
            #pt = self.PointFromPosition(braceOpposite)
            #self.Refresh(True, wxRect(pt.x, pt.y, 5,5))
            #print pt
            #self.Refresh(False)

    def OnMarginClick(self, evt):
        # fold and unfold as needed
        if evt.GetMargin() == 2:
            lineClicked = self.LineFromPosition(evt.GetPosition())
            if self.GetFoldLevel(lineClicked) & wxSTC_FOLDLEVELHEADERFLAG:
                if evt.GetShift():
                    self.SetFoldExpanded(lineClicked, True)
                    self.Expand(lineClicked, True, True, 1)
                elif evt.GetControl():
                    if self.GetFoldExpanded(lineClicked):
                        self.SetFoldExpanded(lineClicked, False)
                        self.Expand(lineClicked, False, True, 0)
                    else:
                        self.SetFoldExpanded(lineClicked, True)
                        self.Expand(lineClicked, True, True, 100)
                else:
                    self.ToggleFold(lineClicked)

    def FoldAll(self):
        lineCount = self.GetLineCount()
        expanding = True

        # find out if we are folding or unfolding
        for lineNum in range(lineCount):
            if self.GetFoldLevel(lineNum) & wxSTC_FOLDLEVELHEADERFLAG:
                expanding = not self.GetFoldExpanded(lineNum)
                break;

        lineNum = 0
        while lineNum < lineCount:
            level = self.GetFoldLevel(lineNum)
            if level & wxSTC_FOLDLEVELHEADERFLAG and \
               (level & wxSTC_FOLDLEVELNUMBERMASK) == wxSTC_FOLDLEVELBASE:

                if expanding:
                    self.SetFoldExpanded(lineNum, True)
                    lineNum = self.Expand(lineNum, True)
                    lineNum = lineNum - 1
                else:
                    lastChild = self.GetLastChild(lineNum, -1)
                    self.SetFoldExpanded(lineNum, False)
                    if lastChild > lineNum:
                        self.HideLines(lineNum+1, lastChild)

            lineNum = lineNum + 1

    def Expand(self, line, doExpand, force=False, visLevels=0, level=-1):
        lastChild = self.GetLastChild(line, level)
        line = line + 1
        while line <= lastChild:
            if force:
                if visLevels > 0:
                    self.ShowLines(line, line)
                else:
                    self.HideLines(line, line)
            else:
                if doExpand:
                    self.ShowLines(line, line)

            if level == -1:
                level = self.GetFoldLevel(line)

            if level & wxSTC_FOLDLEVELHEADERFLAG:
                if force:
                    if visLevels > 1:
                        self.SetFoldExpanded(line, True)
                    else:
                        self.SetFoldExpanded(line, False)
                    line = self.Expand(line, doExpand, force, visLevels-1)

                else:
                    if doExpand and self.GetFoldExpanded(line):
                        line = self.Expand(line, True, force, visLevels-1)
                    else:
                        line = self.Expand(line, False, force, visLevels-1)
            else:
                line = line + 1;

        return line
#-------------- end of copied code from wxStyledTextCtrl_2 demo --------------
#(really I copied alot more, but that part I didn't modify at all, I 
#wanted to understand it, but it just worked, so there was no need)

#------------------------- Ahh, the tabbed notebook --------------------------
class MyNB(wxNotebook):
    def __init__(self, root, id, parent):
        #the left-tab, while being nice, turns the text sideways, ick.
        wxNotebook.__init__(self, parent, id, style=
                            wxNB_TOP
                            #wxNB_BOTTOM
                            #wxNB_LEFT
                            #wxNB_RIGHT
                            )
        self.root = root

        #for some reason, the notebook needs the next line...the text control
        #doesn't.
        EVT_KEY_DOWN(self, self.root.OnKeyPressed)

        EVT_NOTEBOOK_PAGE_CHANGED(self, self.GetId(), self.OnPageChanged)
        
        self.SetDropTarget(FileDropTarget(self.root))

    def __iter__(self):
        count = self.GetPageCount()
        cur = 0
        while cur < count:
            r = self.GetPage(cur).GetWindow1()
            yield r
            cur += 1

    def OnPageChanged(self, event):
        new = event.GetSelection()
        if new > -1:
            win = self.GetPage(new).GetWindow1()
            #fix for dealing with current paths.  They are wonderful.
            if win.dirname:
                os.chdir(win.dirname)
            if VS[-1] == 'u':
                self.root.SetStatusText(win.enc, 2)
                if self.root.HAS_RADIO:
                    self.root.menubar.Check(ENCODINGS[win.enc], 1)
            if self.root.HAS_RADIO:
                self.root.menubar.Check(lexers2[win.GetLexer()], 1)
                self.root.menubar.Check(LL_RMAPPING[win.GetEdgeMode()], 1)
            for m,id in ((0, NUM), (1, MARGIN)):
                self.root.menubar.Check(id, bool(win.GetMarginWidth(m)))
            self.root.menubar.Check(INDENTGUIDE, win.GetIndentationGuides())
            self.root.menubar.Check(USETABS, win.GetUseTabs())
            self.root.menubar.Check(AUTO, win.showautocomp)
            self.root.menubar.Check(WRAPL, win.GetWrapMode() != wxSTC_WRAP_NONE)
            self.root.menubar.Check(SORTBY, win.tree.tree.SORTTREE)
            self.root.menubar.Check(LE_RMAPPING[win.GetEOLMode()], 1)
            self.root.menubar.Check(SAVE_CURSOR, win.save_cursor)
            width = self.GetClientSize()[0]
            split = win.parent
            split.SetSashPosition(max(width-win.SAVEDPOSITION, int(width/2)))
            #if win.GetWrapMode() == wxSTC_WRAP_NONE:
            #    self.parent.SetStatusText("", 1)
            #else:
            #    self.parent.SetStatusText("WRAP",1)
        if event:
            event.Skip()

#----------------- This deals with the tab swapping support. -----------------
    def swapPages(self, p1, p2, moveleft):
        mn = min(p1,p2)
        mx = max(p1,p2)
        if not (moveleft or (p1-p2==1)):
            mn, mx = mx, mn
        page = self.GetPage(mn)
        text = self.GetPageText(mn)[:]
        self.RemovePage(mn)
        self.InsertPage(mx, page, text, 1)

class FileDropTarget(wxFileDropTarget):
    def __init__(self, root):
        wxFileDropTarget.__init__(self)
        self.root = root
    def OnDropFiles(self, x, y, filenames):
        self.root.OnDrop(filenames)

class TextDropTarget(wxTextDropTarget):
    def __init__(self, parent):
        wxTextDropTarget.__init__(self)
        self.parent = parent

    def OnDropText(self, x, y, text):
        self.parent.OnDropText(text)
        #This is so that you can keep your document unchanged while adding
        #code snippets.
        return False

class CodeSnippet(wxPanel):
    def __init__(self, root, id, parent, display2code, displayorder):
        wxPanel.__init__(self, parent, id, style=wxWANTS_CHARS)
        self.root = root

        self.display2code = display2code
        self.displayorder = displayorder
        self.dirty = 0
        
        LB_ID = wxNewId()

        self.lb = wxListBox(self, LB_ID, choices = self.displayorder, style=wxLB_SINGLE|wxLB_NEEDED_SB)
        EVT_LISTBOX_DCLICK(self, LB_ID, self.OnListBoxDClick)
        EVT_KEY_DOWN(self.lb, self.OnKeyPressed)
        EVT_KEY_DOWN(self, self.OnKeyPressed)

        self.SetDropTarget(TextDropTarget(self))
        self.lb.SetDropTarget(TextDropTarget(self))

        sizer = RowColSizer()
        sizer.Add(self.lb, flag=wxEXPAND, row=0, col=0)
        sizer.AddGrowableRow(0)
        sizer.AddGrowableCol(0)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def lb_refresh(self):
        self.lb.Clear()
        self.lb.Set(self.displayorder)

    def onpaste(self, e):
            wxTheClipboard.Open()
            do = wxTextDataObject()
            tmp = wxTheClipboard.GetData(do)
            wxTheClipboard.Close()
            if tmp:
                self.OnDropText(do.GetText())
    def OnKeyPressed(self,e):
        key = e.KeyCode()
        if key in [WXK_DELETE, WXK_BACK] and self.displayorder:
            self.OnListBoxDelete(e)
        elif key == 86 and e.ControlDown() and e.ShiftDown():
            self.onpaste(e)
        elif key == 13:
            self.OnListBoxDClick(e)
        else:
            e.Skip()

    def OnListBoxDClick(self, e):
        num, win = self.root.getNumWin(e)
        sel = self.lb.GetSelection()
        if sel < 0:
            return e.Skip()
        code = self.display2code.get(self.displayorder[sel], '')
        if code != '':
            win.ReplaceSelection(code)
            win.MakeDirty()
            win.SetFocus()

    def OnListBoxDelete(self, e):
        disp = self.lb.GetSelection()
        if disp >= 0:
            self.dirty = 1
            a = self.displayorder[disp]
            del self.display2code[a]
            self.lb.Delete(disp)
            self.displayorder.pop(disp)
            cnt = self.lb.GetCount()
            #if this is the last snippet in the list, select the previous
            #snippet (if one exists)
            if disp == cnt:
                if cnt > 0:
                    self.lb.SetSelection(disp-1)
            #otherwise select the snippet that is just after this one
            else:
                self.lb.SetSelection(disp)
        else:   e.Skip()
    
    def OnDropText(self, text):
        rpt = 0
        str1 = 'What would you like this code snippet to be called?'
        err = {0:'',
               1:'Snippet name already used, please choose another.\n\n',
               2:'Snippet is composed entirely of blank characters,\nplease choose another.\n\n'
               }
        while 1:
            dlg = wxTextEntryDialog(self, err[rpt]+str1, 'Snippet Name')
            if dlg.ShowModal() == wxID_OK:
                nam = dlg.GetValue().strip()
            else:
                nam = None
            dlg.Destroy()
            if nam != None:
                if not nam: rpt = 2;self.root.SetStatusText("Snippet name will not display")
                elif not nam in self.display2code:
                    self.display2code[nam] = text[:]
                    self.displayorder.append(nam)
                    self.lb.Append(nam)
                    self.dirty = 1
                    return self.root.SetStatusText("Snippet %s added"%nam)
                else: rpt = 1;self.root.SetStatusText("Snippet name already used")
            else: self.root.SetStatusText("Snippet addition aborted");return

    def seladj(self, incr):
        cnt = self.lb.GetCount()
        if cnt:
            self.lb.SetSelection((self.lb.GetSelection()+incr)%cnt)

    def OnSnippetP(self, e):
        self.seladj(-1)
    def OnSnippetN(self, e):
        self.seladj(1)

class RunShell(wxMenu):
    def __init__(self, root, parent, label):
        wxMenu.__init__(self)
        menuAddM(parent, self, label)
        self.root = root
        
        menuAdd(self.root, self, "Run current file",    "Run the current open file", self.runScript, wxNewId())
        menuAdd(self.root, self, "View shell commands", "View the list of shell command working paths and the command that would be executed", self.viewShellCommands, wxNewId())
        menuAdd(self.root, self, "Edit shell command",  "Edit a pre-existing shell command", self.editShellCommand, wxNewId())
        menuAdd(self.root, self, "Add shell command",   "Add a shell command to this menu", self.addShellCommand, wxNewId())
        menuAdd(self.root, self, "Remove shell command", "Remove a shell command from this menu", self.removeShellCommand, wxNewId())
        self.AppendSeparator()
        
        #parse settings file, adding the menu
        #format will be:
        #[("name in menu", "working path", "command args"), ...]
        self.menu = {}
        self.order = []
        for i in shellcommands:
            self.appendShellCommand(*i)
        self.dirty = 0
        self.cnt = 6

    def runScript(self, evt):
        num, win = self.root.getNumWin(evt)
        use = self.makeuse(win.dirname, win.filename)
        #print use
        #os.spawnv(os.P_NOWAIT, spawnargs[0], spawnargs + ['--exec', win.dirname, win.filename])
        x = (runme, use['path'], use['full'])
        a = "%s --exec %s %s"%x
        #print a.encode('utf-8')
        #self.root.dialog('path: %s\nfile: %s'%(win.dirname, win.filename), 'running this command')
        try:
            os.system(a)
        except UnicodeError:
            try:
                for i in x:
                    i.encode('ascii')
            except:
                self.root.dialog("The path to: %s\ncontains unicode characters that Python's\nprocess spawning cannot handle."%i, "Error")
        except:
            self.root.exceptDialog()

    def viewShellCommands(self, event, titl="Shell Commands", styl=wxOK, sel=0):
        count = 1
        out = []
        for i in self.order:
            out.append("%i          %s    "%(count, repr(self.menu[i])))
            count += 1
        dlg = wxSingleChoiceDialog(self.root, "order      (name, working dir, command, arguments)", titl, out, styl|wxDD_DEFAULT_STYLE)
        if self.order:
            dlg.SetSelection(min(len(self.order)-1, sel))
        if dlg.ShowModal() == wxID_OK:
            a = dlg.GetStringSelection()
            if a:
                retr = int(dlg.GetStringSelection()[:].split(' ', 1)[0])-1
            else:
                retr = None
        else:
            retr = None
        dlg.Destroy()
        return retr

    def editShellCommand(self, event):
        title = "Choose a shell command to edit"
        style = wxOK|wxCANCEL
        tmp = self.viewShellCommands(event, title, style)
        while tmp != None:
            a = self.order[tmp]
            a1, a2, a3 = self.menu[a]
            if self.addShellCommand(event, a1, a2, a3, tmp) != None:
                del self.order[tmp+1]
                del self.menu[a]
                self.Delete(a)
            self.root.SetStatusText("Edited shell command named %s"%a1)
            tmp = self.viewShellCommands(event, title, style, min(tmp, len(self.order)))

    def addShellCommand(self, event, a1='', a2='', a3='', posn=-1):
        commandname = a1
        recentcwd = a2 
        recentexec = a3
        scaa = "Shell command addition aborted"
        ex = {'path':runpath, 'file':os.path.split(sys.argv[0])[1]}
        ex['full'] = os.path.join(ex['path'], ex['file'])
        hlp = """What is given is completely unchanged, except for those that include the
        below string expansions:
        %(path)s  ->  The path of the file that is currently open
        %(file)s    ->  The name of the file that is currently open
        %(full)s    ->  The equivalent of os.path.join(path, file)
        All other expansions will cause an exception dialog, and your shell
        command will not run.
        
        Spaces in the above expansions, are handled properly for your system.
        Spaces in the paths you enter manually are not repaired.
        If you have paths or arguments that require spaces, other than the above
        expansions, you must repair them yourself.
        For example:
        If I wanted to ssh to another machine to run a command, like ps mayhaps...
        In *nix arguments: hostname ps -u user
        In windows arguments: hostname "ps -u user"
        
        If we were passing files as arguments while connecting to a *nix host...
        In *nix arguments: hostname cat file\\ name\\ with\\ spaces
        In windows arguments: hostname "cat file\\ name\\ with\\ spaces"
        
        Play around, you'll figure it out.
        """.replace('        ', '')

        while 1:
            dlg = wxTextEntryDialog(self.root, "What would you like your command to be called?", "Command name")
            dlg.SetValue(commandname)
            rslt = dlg.ShowModal()
            commandname = dlg.GetValue()
            dlg.Destroy()
            if rslt != wxID_OK:
                return self.root.SetStatusText(scaa)
            
            dlg = wxTextEntryDialog(self.root, hlp, "Choose a working directory (if any)")
            dlg.SetValue(recentcwd)
            rslt = dlg.ShowModal()
            recentcwd = dlg.GetValue()
            dlg.Destroy()
            if rslt != wxID_OK:
                return self.root.SetStatusText(scaa)
            
            dlg = wxTextEntryDialog(self.root, hlp, "Enter the command as you would in a console (including arguments)")
            dlg.SetValue(recentexec)
            rslt = dlg.ShowModal()
            recentexec = dlg.GetValue()
            dlg.Destroy()
            if rslt != wxID_OK:
                return self.root.SetStatusText(scaa)
            
            cmd = recentexec
            
            out = "The functional equivalent to the following will be\nevalutated when the menu item is started, if the\nsource file for this editor was currently being viewed:\n\n"
            if recentcwd:
                out += 'cd %s\n'%recentcwd%ex
            out += cmd%ex
            if self.root.dialog(out, "Is this correct?", styl=wxYES_NO) == wxID_YES:
                self.appendShellCommand(commandname, recentcwd, recentexec, posn)
                self.root.SetStatusText("Added shell command")
                return 1
            else:
                self.root.SetStatusText("Editing shell command")

    def removeShellCommand(self, event):
        title = "Choose a shell command to remove"
        style = wxOK|wxCANCEL
        tmp = self.viewShellCommands(event, title, style)
        while tmp != None:
            a = self.order.pop(tmp)
            b = self.menu[a]
            del self.menu[a]
            self.Delete(a)
            self.root.SetStatusText("Removed shell command named %s"%b[0])
            tmp = self.viewShellCommands(event, title, style, min(tmp, len(self.order)))

    def OnMenu(self, evt):
        try:
            (name, path, command) = self.menu[evt.GetId()]
            try:
                num, win = self.root.getNumWin(evt)
                dn = win.dirname
                fn = win.filename
            except:
                dn = '.'
                fn = ''

            use = self.makeuse(dn, fn)
            use2 = self.makeuse(path, command)
            path = use2['path']
            #command = use2['full']
            #os.spawnv(os.P_NOWAIT, spawnargs[0], spawnargs + ['--exec', path or '*', command])
            if sys.platform=='win32':
                os.system("start %s --exec %s %s"%(runme, path or '*', command)%use)
            else:
                os.system("%s --exec %s %s"%(runme, path or '*', command)%use)
        except:
            self.root.exceptDialog()

    def makeuse(self, path, command):
        use = {}
        use['path'] = fixpath(path)
        use['file'] = fixpath(command)
        use['full'] = fixpath(os.path.join(path, command))
        return use

    def appendShellCommand(self, name, path, command, posn=-1):
        self.dirty = 1
        a = wxNewId()
        if posn == -1:
            self.Append(a, name)
            self.order.append(a)
        else:
            self.Insert(posn+self.cnt, a, name)
            self.order.insert(posn, a)
        EVT_MENU(self.root, a, self.OnMenu)
        self.menu[a] = (name, path, command)

class hierCodeTreePanel(wxPanel):
    class TreeCtrl(wxTreeCtrl):
        def __init__(self, parent, tid):
            wxTreeCtrl.__init__(self, parent, tid, style=wxTR_DEFAULT_STYLE|wxTR_HAS_BUTTONS|wxTR_HIDE_ROOT)
            
            isz = (16,16)
            il = wxImageList(isz[0], isz[1])
            self.images = [wxArtProvider_GetBitmap(wxART_FOLDER, wxART_OTHER, isz),
                           wxArtProvider_GetBitmap(wxART_FILE_OPEN, wxART_OTHER, isz),
                           wxArtProvider_GetBitmap(wxART_HELP_PAGE, wxART_OTHER, isz)]
            for i in self.images:
                il.Add(i)
            self.SetImageList(il)
            self.il = il

        def OnCompareItems(self, item1, item2):
            if self.SORTTREE:
                return cmp(self.GetItemData(item1).GetData(),\
                           self.GetItemData(item2).GetData())
            else:
                return cmp(self.GetItemData(item1).GetData()[1],\
                           self.GetItemData(item2).GetData()[1])

        def getchlist(self, parent):
            #when using ItemHasChildren(parent)
            #an assert exception saying "can't retrieve virtual root item"
            #kept coming up.
            #for some reason GetChildrenCount(parent, 0) doesn't.
            numchildren = self.GetChildrenCount(parent, 0) 
            lst = {}
            if numchildren>0:
                ch, cookie = self.GetFirstChild(parent, wxNewId())
                lst[self.GetItemText(ch)] = [ch]
                for i in xrange(numchildren-1):
                    ch, cookie = self.GetNextChild(parent, cookie)
                    txt = self.GetItemText(ch)
                    if txt in lst:
                        lst[txt].append(ch)
                    else:
                        lst[txt] = [ch]
            return lst

    def __init__(self, root, parent):
        # Use the WANTS_CHARS style so the panel doesn't eat the Return key.
        wxPanel.__init__(self, parent, -1, style=wxWANTS_CHARS)
        EVT_SIZE(self, self.OnSize)

        self.root = root

        tID = wxNewId()

        self.tree = self.TreeCtrl(self, tID)
        self.tree.troot = self.tree.AddRoot("Unseen Root")

        #self.tree.Expand(self.root)
        EVT_LEFT_DCLICK(self, self.OnLeftDClick)
        EVT_TREE_ITEM_ACTIVATED(self, tID, self.OnActivate)

    def new_hierarchy(self, hier):
        #self.tree.DeleteAllItems()
        root = [self.tree.troot]
        stk = [(self.tree.getchlist(self.tree.troot), hier)]
        #D = {'c':wxColour(0, 0, 200),
        #     'd':wxColour(200, 0, 0)}
        while stk:
            chlist, cur = stk.pop()
            while cur:
                name, line_no, leading, children = cur.pop()
                if name in chlist:
                    item_no = chlist[name].pop()
                    if not chlist[name]:
                        del chlist[name]
                    self.tree.SetPyData(item_no, line_no)
                    icons = 0
                else:
                    item_no = self.tree.AppendItem(root[-1], name)
                    self.tree.SetPyData(item_no, line_no)
                    #self.tree.SetItemTextColour(item_no, D[name[0]])
                    icons = 1
                if children:
                    if icons:
                        self.tree.SetItemImage(item_no, 0, wxTreeItemIcon_Normal)
                        self.tree.SetItemImage(item_no, 1, wxTreeItemIcon_Expanded)
                    stk.append((chlist, cur))
                    root.append(item_no)
                    chlist = self.tree.getchlist(item_no)
                    cur = children
                elif icons:
                    self.tree.SetItemImage(item_no, 2, wxTreeItemIcon_Normal)
                    self.tree.SetItemImage(item_no, 2, wxTreeItemIcon_Selected)
            for j in chlist.itervalues():
                for i in j:
                    self.tree.DeleteChildren(i)
                    self.tree.Delete(i)
            self.tree.SortChildren(root[-1])
            root.pop()

    def OnSize(self, event):
        w,h = self.GetClientSizeTuple()
        self.tree.SetDimensions(0, 0, w, h)

    def OnLeftDClick(self, event):
        #pity this doesn't do what it should.
        num, win = self.root.getNumWin(event)
        win.SetFocus()

    def OnActivate(self, event):
        num, win = self.root.getNumWin(event)
        dat = self.tree.GetItemData(event.GetItem()).GetData() 
        if dat == None:
            return event.Skip()
        ln = dat[1]
        #print ln
        #print dir(win)
        linepos = win.GetLineEndPosition(ln)
        win.EnsureVisible(ln)
        win.SetSelection(linepos-len(win.GetLine(ln))+len(win.format), linepos)
        win.ScrollToColumn(0)
        win.SetFocus()

class VirtualTodo(wxPanel):
    class virtualTodo(wxListCtrl):
        def __init__(self, parent, root, stc):
            wxListCtrl.__init__(self, parent, -1,
                                style=wxLC_REPORT|wxLC_VIRTUAL|wxLC_HRULES|wxLC_VRULES)
            self.parent = parent
            self.root = root
            self.stc = stc
    
            self.InsertColumn(0, "Line")
            self.InsertColumn(1, "!")
            self.InsertColumn(2, "Todo")
            self.SetColumnWidth(0, 40)
            self.SetColumnWidth(1, 30)
            self.SetColumnWidth(2, 400)
    
            self.items = []
            self.SetItemCount(0)
    
            EVT_LIST_ITEM_ACTIVATED(self, self.GetId(), self.OnItemActivated)
    
        def OnItemActivated(self, event):
            sel = self.items[event.m_itemIndex][0]
            win = self.stc
            if sel < win.GetLineCount():
                linepos = win.GetLineEndPosition(sel)
                win.EnsureVisible(sel)
                win.SetSelection(linepos-len(win.GetLine(sel))+len(win.format), linepos)
                win.ScrollToColumn(0)
    
        def getColumnText(self, index, col):
            return str(self.items[index][col])
    
        def OnGetItemText(self, item, col):
            return str(self.items[item][col])
    
        def OnGetItemImage(self, item):
            return -1
    
        def OnGetItemAttr(self, item):
            return None
    
    def __init__(self, parent, root, stc):
        wxPanel.__init__(self, parent, -1, style=wxWANTS_CHARS)
        self.vtd = self.virtualTodo(self, root, stc)
        
        self.parent = parent
        self.root = root
        
        EVT_SIZE(self, self.OnSize)

    def NewItemList(self, items):
        self.vtd.items = []
        self.vtd.SetItemCount(0)
        self.vtd.items = items
        self.vtd.SetItemCount(len(items))

    def OnSize(self, event):
        w,h = self.GetClientSizeTuple()
        self.vtd.SetDimensions(0, 0, w, h)




class KeySink(wxWindow):
    def __init__(self, parent, defa, curr):
        wxWindow.__init__(self, parent, -1, style=wxWANTS_CHARS)
        self.SetBackgroundColour(wxBLUE)
        self.haveFocus = False
        self.defa = defa
        self.curr = curr
        

        EVT_PAINT(self, self.OnPaint)
        EVT_SET_FOCUS(self, self.OnSetFocus)
        EVT_KILL_FOCUS(self, self.OnKillFocus)
        EVT_MOUSE_EVENTS(self, self.OnMouse)
        EVT_KEY_DOWN(self, self.OnKey)

    def OnPaint(self, evt):
        dc = wxPaintDC(self)
        rect = self.GetClientRect()
        dc.SetTextForeground(wxWHITE)
        dc.DrawLabel("Click here and enter your key combination.\n"
                     "Default: %s  Current: %s"%(self.defa, self.curr),
                     rect, wxALIGN_CENTER | wxALIGN_TOP)
        if self.haveFocus:
            dc.SetTextForeground(wxGREEN)
            dc.DrawLabel("Have Focus", rect, wxALIGN_RIGHT | wxALIGN_BOTTOM)
        else:
            dc.SetTextForeground(wxRED)
            dc.DrawLabel("Need Focus!", rect, wxALIGN_RIGHT | wxALIGN_BOTTOM)

    def OnSetFocus(self, evt):
        self.haveFocus = True
        self.Refresh()

    def OnKillFocus(self, evt):
        self.haveFocus = False
        self.Refresh()

    def OnMouse(self, evt):
        if evt.ButtonDown():
            self.SetFocus()

    def OnKey(self, evt):
        self.GetParent().text.SetValue(GetKeyPress(evt))

class GetKeyDialog(wxDialog):
    def __init__(self, parent, defa, curr):
        wxDialog.__init__(self, parent, -1, "Enter key combination",
                          style=wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER,
                          size=(320, 240))
        self.parent = parent
        
        self.keysink = KeySink(self, defa, curr)
        
        self.text = wxTextCtrl(self, -1, curr, style=wxTE_READONLY|wxTE_LEFT)
        
        buttons = wxBoxSizer(wxHORIZONTAL)
        ok = wxButton(self, wxOK, "OK")
        cancel = wxButton(self, wxCANCEL, "Cancel")
        buttons.Add(ok, 0)
        buttons.Add(cancel, 0)
        
        sizer = wxBoxSizer(wxVERTICAL)
        sizer.Add(self.keysink, 1, wxGROW)
        sizer.Add(self.text, 0, wxGROW)
        sizer.Add(buttons, 0, wxALIGN_RIGHT)

        self.SetSizer(sizer)
        EVT_BUTTON(self, wxOK, self.OnOK)
        EVT_BUTTON(self, wxCANCEL, self.OnCancel)
        EVT_CLOSE(self, self.OnClose)

    def OnClose(self, evt):
        self.OnCancel(evt)

    def OnOK(self, evt):
        self.parent.accelerator = self.text.GetValue()
        self.Destroy()
    
    def OnCancel(self, evt):
        self.parent.accelerator = ""
        self.Destroy()

class MenuItemDialog(wxDialog):
    class VirtualMenu(wxPanel):
        class virtualMenu(wxListCtrl, wxListCtrlAutoWidthMixin):
            def __init__(self, parent):
                wxListCtrl.__init__(self, parent, -1,
                                    style=wxLC_REPORT|wxLC_VIRTUAL|wxLC_HRULES|wxLC_VRULES|wxLC_SINGLE_SEL)
                wxListCtrlAutoWidthMixin.__init__(self)
                
                self.parent = parent
                self.sel = None
    
                self.InsertColumn(0, "Default Menu Heirarchy")
                self.InsertColumn(1, "Current Name")
                self.InsertColumn(2, "Default Hotkey")
                self.InsertColumn(3, "Current Hotkey")
                self.SetColumnWidth(0, 250)
                self.SetColumnWidth(1, 150)
                self.SetColumnWidth(2, 100)
        
                self.items = []
                self.SetItemCount(0)
                
                EVT_LIST_ITEM_SELECTED(self, self.GetId(), self.OnItemSelected)
                EVT_LIST_ITEM_ACTIVATED(self, self.GetId(), self.OnItemActivated)
        
            def OnItemActivated(self, event):
                inum = event.GetIndex()
                item = self.items[inum]
                dlg = wxTextEntryDialog(self,
                                        "Enter the new name of the menu item\n"\
                                        "Default: %s  Current: %s"%(item[0].split('->')[-1], item[1]),\
                                        "What should the item be called?")
                dlg.SetValue(item[1])
                rslt = dlg.ShowModal()
                if rslt == wxID_OK:
                    name = dlg.GetValue()
                else:
                    name = None
                dlg.Destroy()
                if not name:
                    return
                if item[0].find('->') == -1 or not item[4]:
                    self.items[inum] = (item[0], name, '', '', 0)
                    return self.RefreshItem(inum)
                dlg = GetKeyDialog(self, item[2], item[3])
                dlg.ShowModal()
                if not self.accelerator:
                    return
                self.items[inum] = (item[0], name, item[2], self.accelerator, 1)
                self.RefreshItem(inum)
        
            def getColumnText(self, index, col):
                return self.items[index][col]
        
            def OnGetItemText(self, item, col):
                return self.items[item][col]
        
            def OnGetItemImage(self, item):
                return -1
        
            def OnGetItemAttr(self, item):
                return None
            def OnItemSelected(self, evt):
                self.sel = evt.GetIndex()
        
        def __init__(self, parent):
            wxPanel.__init__(self, parent, -1, style=wxWANTS_CHARS)
            self.vm = self.virtualMenu(self)
            self.refresh(MENULIST)
            
            self.parent = parent
            
            EVT_SIZE(self, self.OnSize)

        def refresh(self, items):
            self.vm.items = []
            self.vm.SetItemCount(0)
            self.vm.items = items
            self.vm.SetItemCount(len(items))
            if self.vm.sel is not None:
                self.vm.EnsureVisible(self.vm.sel)

        def OnSize(self, event):
            w,h = self.GetClientSizeTuple()
            self.vm.SetDimensions(0, 0, w, h)
    
    def __init__(self, parent, root):
        wxDialog.__init__(self, parent, -1, "Menu item names and hotkey bindings.",
                          size=(640,480), style=wxRESIZE_BORDER|wxDEFAULT_DIALOG_STYLE)
        
        self.root = root
        self.vmu = self.VirtualMenu(self)
        sizer = wxBoxSizer(wxVERTICAL)
        
        sizer.Add(wxStaticText(self, -1,
                               "Double click on an entry to change the name and/or hotkey.  Most any hotkey should work fine.  Give it a try.\n"
                               "Name changes with accelerators should work fine, as long as there are no accelerator key collisions."))
        sizer.Add(self.vmu, 1, wxGROW)
        
        s2 = wxBoxSizer(wxHORIZONTAL)

        ok = wxButton(self, wxOK, "OK")
        clear = wxButton(self, wxNewId(), "Clear Hotkey")
        revert = wxButton(self, wxNewId(), "Revert Name and Hotkey")
        cancel = wxButton(self, wxCANCEL, "Cancel")
        
        s2.Add(ok)
        s2.Add(clear)
        s2.Add(revert)
        s2.Add(cancel)
        sizer.Add(s2, 0, wxALIGN_RIGHT)
        
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        
        EVT_BUTTON(self, wxOK, self.OnOK)
        EVT_BUTTON(self, clear.GetId(), self.OnClear)
        EVT_BUTTON(self, revert.GetId(), self.OnRevert)
        EVT_BUTTON(self, wxCANCEL, self.OnCancel)
    
    def OnOK(self, evt):
        global MENUPREF
        global MENULIST
        nmu = {}
        changed = 0
        for heir, cn, da, ca, hk in self.vmu.vm.items:
            if MENUPREF[heir] != (cn, ca):
                changed = 1
            nmu[heir] = (cn, ca)
        if changed:
            MENUPREF.clear()
            MENUPREF.update(nmu)
            MENULIST[:] = self.vmu.vm.items
            self.root.dialog("You must restart in order for your\n"
                             "menu item changes to go into effect.",
                             "Restart Required")
        self.Destroy()

    def OnClear(self, evt):
        IL = self.vmu.vm
        items = IL.items
        indx = IL.sel
        items[indx] = items[indx][:3] + ('', items[indx][4])
        IL.RefreshItem(indx)

    def OnRevert(self, evt):
        IL = self.vmu.vm
        items = IL.items
        indx = IL.sel
        item = items[indx]
        items[indx] = (item[0], item[0].split('->')[-1], item[2], item[2], item[4])
        IL.RefreshItem(indx)

    def OnCancel(self, evt):
        self.Destroy()

class lastused:
    #Implementation of a length-limited O(1) LRU cache
    class Node:
        def __init__(self, prev, me):
            self.prev = prev
            self.me = me
            self.next = None
    def __init__(self, count, lst=[]):
        self.count = max(count, 2)
        self.d = {}
        self.first = None
        self.last = None
        for key, value in lst:
            self[key] = value
    def __contains__(self, obj):
        return obj in self.d
    def __getitem__(self, obj):
        return self.d[obj].me[1]
    def __setitem__(self, obj, val):
        if obj in self.d:
            del self[obj]
        nobj = self.Node(self.last, (obj, val))
        if self.first is None:
            self.first = nobj
        if self.last:
            self.last.next = nobj
        self.last = nobj
        self.d[obj] = nobj
        if len(self.d) > self.count:
            if self.first == self.last:
                self.first = None
                self.last = None
                return
            a = self.first
            a.next.prev = None
            self.first = a.next
            a.next = None
            del self.d[a.me[0]]
            del a
    def __delitem__(self, obj):
        nobj = self.d[obj]
        if nobj.prev:
            nobj.prev.next = nobj.next
        else:
            self.first = nobj.next
        if nobj.next:
            nobj.next.prev = nobj.prev
        else:
            self.last = nobj.prev
        del self.d[obj]
    def __iter__(self):
        cur = self.first
        while cur != None:
            cur2 = cur.next
            yield cur.me[1]
            cur = cur2
        raise StopIteration
    def iteritems(self):
        cur = self.first
        while cur != None:
            cur2 = cur.next
            yield cur.me
            cur = cur2
        raise StopIteration
    def iterkeys(self):
        return iter(self.d)
    def itervalues(self):
        for i,j in self.iteritems():
            yield j
    def keys(self):
        return [i for i,j in self.iteritems()]
    def values(self):
        return [j for i,j in self.iteritems()]
    def items(self):
        return [i for i in self.iteritems()]
        
class findinfiles(wxDialog):
    def __init__(self, parent, choices):
        wxDialog.__init__(self, parent, -1, "Find in files...", size=(640,480),
                          style=wxRESIZE_BORDER|wxDEFAULT_DIALOG_STYLE)
        
        self.parent = parent
        self.root = parent
        
        self.running = 0
        self.stopping = 0
        self.starting = 0
        
        def static(text, style=wxALIGN_LEFT, self=self):
            return wxStaticText(self, -1, text, style=style)

        def checkb(text, val=0, style=wxNO_BORDER):
            a = wxCheckBox(self, -1, text, style=style)
            a.SetValue(val)
            return a
        
        def combobox(ch, d=''):
            a = wxComboBox(self, -1, choices=ch, style=wxCB_DROPDOWN)
            if ch: a.SetSelection(0)
            elif d: a.SetValue(d)
            return a
        
        sizer = RowColSizer()
        
        def sAdd(obj, pos, size=(1,1), flag=wxALIGN_CENTER_VERTICAL|wxEXPAND|wxALL, sizer=sizer):
            sizer.Add(obj, pos=pos, size=size, flag=flag, border=3)
        
        sizer.AddGrowableCol(2)
        sizer.AddGrowableRow(5)
        
        sAdd(static("Multiple directories or extensions "\
                    "should be separated by semicolons   ;"),
                    (0, 1), (1,3))
        
        fc = self.root.config['FIF_STATE']
        
        sAdd(static("Search for:"), (1,1))
        self.search = combobox(fc[0])
        sAdd(self.search, (1,2), (1,2))

        sAdd(static("Directories:"), (2,1))
        self.sdirs = combobox(fc[1])
        sAdd(self.sdirs, (2,2))
        b = wxButton(self, wxNewId(), "Add Path")
        sAdd(b, (2,3))
        EVT_BUTTON(self, b.GetId(), self.OnDirButtonClick)
        
        sAdd(static("Extensions:"), (3,1))
        self.extns = combobox(fc[2], "*.*")
        sAdd(self.extns, (3,2), (1,2))
        

        self.cs = checkb("Case sensitive", fc[3])
        self.ss = checkb("Search Subdirectories", fc[4])
        self.re = checkb("Regular Expression", fc[5])
        hsize = wxBoxSizer(wxHORIZONTAL)
        hsize.Add(self.cs, 5, wxEXPAND)
        hsize.Add(self.ss, 6, wxEXPAND)
        hsize.Add(self.re, 6, wxEXPAND)
        sAdd(hsize, (4,1), (1,3))

        self.results = wxListBox(self, wxNewId(), style=wxLB_SINGLE|wxLB_NEEDED_SB)
        EVT_LISTBOX_DCLICK(self, self.results.GetId(), self.OpenFound)
        sAdd(self.results, (5,1), (1,3))
        
        sAdd(static('Status:'), (6, 1))
        self.status = wxTextCtrl(self, -1)
        self.status.SetValue("Ready.")
        sAdd(self.status, (6, 2))
        
        self.b = wxButton(self, wxNewId(), "Start Search")
        sAdd(self.b, (6,3))
        EVT_BUTTON(self, self.b.GetId(), self.OnFindButtonClick)
        
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        tid = wxNewId()
        self.timer = wxTimer(self, tid)
        EVT_TIMER(self, tid, self.OnFindButtonClick)

    def OpenFound(self, e):
        selected = self.results.GetSelection()
        if selected < 0:
            return
        cur = selected
        while (cur > 0) and (self.results.GetString(cur)[0] == ' '):
            cur -= 1
        fn = self.results.GetString(cur)
        line = 1
        a = self.results.GetString(selected) 
        if a[0] == ' ':
            line = int(a.split(':', 1)[0])
        
        self.root.OnDrop([fn], 1)
        a = self.root.control.GetSelection()
        if a > -1:
            stc = self.root.control.GetPage(a).GetWindow1()
            stc.GotoLine(line)

    def OnDirButtonClick(self, e):
        dlg = wxDirDialog(self, "Choose a directory", style=wxDD_DEFAULT_STYLE)
        a = dlg.ShowModal()
        if a == wxID_OK:
            a = self.sdirs.GetValue()
            if a:
                self.sdirs.SetValue(a+';'+dlg.GetPath())
            else:
                self.sdirs.SetValue(dlg.GetPath())
        dlg.Destroy()

    def OnFindButtonClick(self, e):
        if self.stopping:
            #already stopping
            return
        
        elif self.running:
            #running, we want to try to stop
            self.stopping = 1
            self.status.SetValue("Stopping...please wait.")
            self.b.SetLabel("Start Search")
            return

        if e.GetId() == self.b.GetId():
            if self.starting:
                #previously was waiting to start due to an
                #external yield, abort waiting and cancel
                #search
                self.starting = 0
                self.status.SetValue("Search cancelled.")
                self.b.SetLabel("Start Search")
                return
        elif not self.starting:
            #I was a waiting timer event, but I don't
            #need to start anymore
            return

        #try to start
        self.starting = 1
        self.b.SetLabel("Stop Search")

        try:
            wxYield()
        except:
            #would have tried to start while in another function's
            #wxYield() call.  Will wait 100ms and try again
            self.status.SetValue("Waiting for another process to stop, please wait.")
            self.timer.Start(100, wxTIMER_ONE_SHOT)
            return

        #am currently the topmost call, so will continue.
        self.starting = 0
        self.running = 1
        wxYield()

        def getlist(c):
            cc = c.GetCount()
            e = [c.GetString(i) for i in xrange(cc)]
            a = c.GetValue()
            if a in e:
                e.remove(a)
            e = [a] + e
            e = e[:10]
            if len(e) > cc:
                c.Append(e[-1])
            for i in xrange(len(e)):
                c.SetString(i, e[i])
            c.SetSelection(0)
            return e

        FIFS = (getlist(self.search),
              getlist(self.sdirs),
              getlist(self.extns),
              self.cs.IsChecked(),
              self.ss.IsChecked(),
              self.re.IsChecked())
        self.root.config['FIF_STATE'] = FIFS

        search = self.search.GetValue()
        paths = self.sdirs.GetValue().split(';')
        extns = self.extns.GetValue().split(';') or ['*.*']
        case = self.cs.IsChecked()
        subd = self.ss.IsChecked()
        
        sfunct = self.searchST
        if self.re.IsChecked():
            sfunct = self.searchRE
            if case:
                search = re.compile(search)
            else:
                search = re.compile(search, re.IGNORECASE)

        results = self.results
        results.Clear()

        def file_iterator(path, subdirs, extns):
            try:    lst = os.listdir(path)
            except: return
            d = []
            for file in lst:
                a = os.path.join(path, file)
                if os.path.isfile(a):
                    for extn in extns:
                        if fnmatch.fnmatch(file, extn):
                            yield a
                            break
                elif subdirs and os.path.isdir(a):
                    d.append(a)
            if not subdirs:
                return
            for p in d:
                for f in file_iterator(p, subdirs, extns):
                    yield f

        filecount = 0
        filefcount = 0
        foundcount = 0
        ss = "Found %i instances in %i files out of %i files checked."
        for path in paths:
            wxYield()
            if not self.running:
                break
            for filename in file_iterator(path, subd, extns):
                wxYield()
                filecount += 1
                if not self.stopping:
                    r = sfunct(filename, search, case)
                    if r:
                        try:
                            for a in r:
                                results.Append(a)
                        except:
                            #for platforms with limited sized
                            #wxListBox controls
                            pass
                        filefcount += 1
                        foundcount += len(r)-1
                    self.status.SetValue((ss%(foundcount, filefcount, filecount)) + '...searching...')
                else:
                    break
        if self.stopping:
            self.stopping = 0
            ex = '...cancelled.'
            #was stopped by a button press
        else:
            self.running = 0
            ex = '...done.'
        self.b.SetLabel("Start Search")
        self.status.SetValue((ss%(foundcount, filefcount, filecount)) + ex)
        self.stopping = 0
        self.running = 0
        self.starting = 0

    def searchST(self, filename, pattern, casesensitive):
        try:
            lines = open(filename, 'r').readlines()
        except:
            lines = []
        found = []
        if not casesensitive:
            pattern = pattern.lower()
        spt = spaces_per_tab*' '
        for i in range(len(lines)):
            try:
                if not casesensitive:
                    line = lines[i].lower()
                else:
                    line = lines[i]
                if line.find(pattern) > -1:
                    if not found:
                        found.append(filename)
                    found.append('  '+`i+1` + ': '+lines[i].rstrip().replace('\t', spt))
            except: pass
            wxYield()
        return found

    def searchRE(self, filename, pattern, toss):
        try:
            lines = open(filename, 'r').readlines()
        except:
            lines = []
        found = []
        spt = spaces_per_tab*' '
        for i in range(len(lines)):
            try:
                line = lines[i]
                if pattern.search(line) is not None:
                    if not found:
                        found.append(filename)
                    found.append('  '+`i+1` + ': '+line.rstrip().replace('\t', spt))
            except: pass
            wxYield()
        return found

#--------------------------- And the main...*sigh* ---------------------------
import wx
VS = wx.VERSION_STRING
del wx

def main():
    app = wxPySimpleApp()
    if VS < VREQ:
        app.frame = wxFrame(None, -1, "Improper version of wxPython")
        dlg = wxMessageDialog(app.frame, "PyPE is only known to run on wxPython %s or later.\nYou are running version %s of wxPython.\nYou should upgrade your version of wxPython.\nNewer versions of wxPython are available at\nwww.wxpython.org"%(VREQ,VS), "Improper version of wxPython", wxOK|wxICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        return
    #imports for encodings support in the binary windows distribution
    if sys.executable[:6].lower() != 'python':
        import encodings.cp037
        import encodings.cp1006
        import encodings.cp1026
        import encodings.cp1140
        import encodings.cp1250
        import encodings.cp1251
        import encodings.cp1252
        import encodings.cp1253
        import encodings.cp1254
        import encodings.cp1255
        import encodings.cp1256
        import encodings.cp1257
        import encodings.cp1258
        import encodings.cp424
        import encodings.cp437
        import encodings.cp500
        import encodings.cp737
        import encodings.cp775
        import encodings.cp850
        import encodings.cp852
        import encodings.cp855
        import encodings.cp856
        import encodings.cp857
        import encodings.cp860
        import encodings.cp861
        import encodings.cp862
        import encodings.cp863
        import encodings.cp864
        import encodings.cp865
        import encodings.cp866
        import encodings.cp869
        import encodings.cp874
        import encodings.cp875
    if VS[-1] == 'u':
        import encodings.ascii
        import encodings.utf_7
        import encodings.utf_8
        import encodings.utf_16
        import encodings.utf_16_be
        import encodings.utf_16_le
    opn=0
    if len(sys.argv)>1 and (sys.argv[1] == '--last'):
        opn=1
    app.frame = MainWindow(None, -1, "PyPE %s"%VERSION, sys.argv[1+opn:])
    app.SetTopWindow(app.frame)
    app.frame.Show(1)
    if opn:
        app.frame.OnOpenPrevDocs(None)
    app.MainLoop()

if __name__ == '__main__':
    main()

#----------------------------- configuration.py ------------------------------
#----------------------- Settings for configuring PyPE -----------------------

collapse_style = 0

#if the following is set to 1, then drag and drop TEXT support will not work
#in the editor control, but will enable drag and drop FILE support in the
#editor control.  Enabled by default.  Why?  Because I drag files, I cut and
#paste text.
dnd_file = 1

#CTRL-T swaps two lines.  Setting the below to 1 disables this hotkey.
REM_SWAP = 1

import os
import sys
#from wxPython.wx import *
from wxPython.stc import wxSTC_EOL_CRLF, wxSTC_EOL_LF, wxSTC_EOL_CR
from parsers import *

fmt_mode = {"\r\n":wxSTC_EOL_CRLF,
              "\n":wxSTC_EOL_LF,
              "\r":wxSTC_EOL_CR}
fmt_Rmode = {}
for i,j in fmt_mode.items():
    fmt_Rmode[j] = i

eol = os.linesep
eolmode = fmt_mode[eol]

def fp(path):
    if sys.platform=='win32':
        if ' ' in path:
            return '"%s"'%path
        return path
    return path.replace(' ', '\\ ')

def fixpath(path):
    return os.path.normpath(fp(path))

se = fixpath(sys.executable)
spawnargs = [se]

if sys.executable[-8:].lower() == 'pype.exe':
    runme = se
    runpath = se
else:
    runpath = os.path.dirname(os.path.normpath(os.path.abspath(__file__)))
    b = fixpath(os.path.join(runpath, sys.argv[0]))
    runme = "%s %s"%(se,b)
    spawnargs.append(b)
    runpath = __file__
runpath = os.path.dirname(os.path.normpath(os.path.abspath(runpath)))

stylefile = os.path.join(runpath, 'stc-styles.rc.cfg')

def command_parser(command):
    args = []
    beg = 0
    cur = 0
    if os.name == 'nt':
        while cur < len(command):
            if command[cur] == ' ':
                if command[beg] == '"':
                    if cur > (beg+1):
                        if command[cur-1] == '"':
                            args.append(command[beg:cur])
                            cur += 1
                            beg = cur
                        else:
                            cur += 1
                    else:
                        cur += 1
                elif cur > beg:
                    args.append(command[beg:cur])
                    cur += 1
                    beg = cur
                else:
                    cur += 1
                    beg = cur
            else:
                cur += 1
    else:
        while cur < len(command):
            if command[cur:cur+2] == '\\ ':
                cur = cur + 2
            elif command[cur] == ' ':
                args.append(command[beg:cur])
                cur += 1
                beg = cur
            else:
                cur += 1
    if beg != cur:
        args.append(command[beg:cur])
    return args

#for open/save dialogs
wildcard = "All python files (*.py *.pyw)|*.py;*.pyw;*.PY;*.PYW|"\
           "C/C++ files (*.c* *.h)|*.c*;*.h;*.C*;*.H|"\
           "HTML/XML files (*.htm* *.shtm* *.xml)|*.htm*;*.shtm*;*.xml;*.HTM*;*.SHTM*;*.XML|"\
           "All files (*.*)|*.*"

#for style mappings from extensions
extns = {'py' : 'python',
        'pyw' : 'python',
          'c' : 'cpp',
         'cc' : 'cpp',
        'cpp' : 'cpp',
        'c++' : 'cpp',
          'h' : 'cpp',
        'htm' : 'html',
       'html' : 'html',
       'shtm' : 'html',
      'shtml' : 'html',
        'xml' : 'xml',
        'txt' : 'text'}

default_homedir = os.path.dirname(os.path.abspath(__file__))

try:
    #all user-based OSes
    thd = os.path.expanduser("~")
    if thd == "~": raise
    homedir = os.path.join(thd, ".pype")
except:
    try:
        #*nix fallback
        homedir = os.path.join(os.environ['HOME'], ".pype")
    except:
        try:
            #windows NT,2k,XP,etc. fallback
            homedir = os.path.join(os.environ['USERPROFILE'], ".pype")
        except:
            #What os are people using?
            homedir = os.path.join(default_homedir, ".pype")
try:
    # create the config directory if it
    # doesn't already exist
    def expandfull(var, rem=3):
        if not rem:
            return os.path.expandvars(var)
        a = os.path.expandvars(var)
        b = []
        d = [b.extend(i.split('\\')) for i in a.split('/')]
        c = []
        for i in b:
            if '%' in i:
                c.append(expandfull(i, rem-1))
            else:
                c.append(i)
        return '\\'.join(c)
    if eol == "\r\n" and '%' in homedir:
        homedir = expandfull(homedir)
    if not os.path.exists(homedir):
        os.mkdir(homedir)
except:
    #print "unable to create config directory", homedir
    homedir = default_homedir

for fil in os.listdir(homedir):
    if fil.find('.tmp.') > -1:
        try: os.remove(os.path.join(homedir, fil))
        except: pass

def get_paragraphs(text, l_sep):
    in_lines = text.split(l_sep)
    lines = []
    cur = []
    for line in in_lines:
        cur.append(line)
        if line:
            if line[-1] != ' ':
                cur.append(' ')
        else:
            if cur:
                lines.append(cur)
                cur = []
            lines.append([])
    if cur:
        lines.append(cur)
    return [''.join(i) for i in lines]

def wrap_paragraph(text, width):
    words = text.split(' ')
    lines = []
    cur = []
    l = 0
    for word in words:
        lw = len(word)
        if not lw:
            cur.append(word)
        elif (l + len(cur) + len(word)) <= width:
            cur.append(word)
            l += lw
        else:
            if cur[-1]:
                cur.append('')
            lines.append(cur)
            cur = [word]
            l = lw
    if cur:
        lines.append(cur)
    return '\n'.join([' '.join(i) for i in lines])

def wrap_lines(text, width, lend):
    paragraphs = get_paragraphs(text, lend)
    retr = lend.join([wrap_paragraph(i, width) for i in paragraphs])
    return retr

def validate(dlg, orig):
    try:
        a = int(dlg.GetValue())
    except:
        a = 0
    if a < 1:
        return orig
    return a

def getData():
    import zlib
    a = zlib.decompress(
'x\xda\x01\xf4\x04\x0b\xfb\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00 \
\x00\x00\x00 \x08\x06\x00\x00\x00szz\xf4\x00\x00\x00\x04sBIT\x08\x08\x08\x08\
|\x08d\x88\x00\x00\x04\xabIDATx\x9c\xb5\x97\xcdk\x1bG\x18\xc6\x7f\x93f\x9d \
\xc5\xd2\xacI\x94\x14Y\x91\xe3B\x1a\xb0\x14(-!\x10\xd5\xd4\xa4\x1fid\x8a/=\
\xd4\xb9\xf4Tj\x9f\n=\xd4\xfe\x0b\xd2K\xae=\xe4n\xe5TP\x8b\x95@\xd3\x90b\xcb\
\xd4\x18\xd2\x84:\x854\x07G\xb6%bKP\xad\x94nh\xbd\x84\xe9a\xbd\xab]}Y\xf9\
\xe8\x03\x02\xed\xec\xec<\xcf\xfb\xcc;\xef\xcc\x08)%\x9d\x106\x0cU\x03\xc2@M\
J\xd1\xb1\xe3K`\x7f7\xf2\x82R\xee\xf3\x90\x10\xea\xff\x10\xd1V@39@A)\x86\x84\
\xf05>\xaf \xcb\n\xa9Xt\x9c\xc8\x91\x14\xba~\x1a\xd3\\o\x15\x106\x0c\x97D\
\x17\xf6\xf8\xd5]15\xcf\x7f\xe8\xdd\x15\xcb\n\xa9\xe4\xc8,\xc3C\x93\x9e\xb6:\
Ec\xb5!\xc0%\xceBm\x02\x1f\xb1#\xa6\xda\xc1\x95n"\x06dZ\x9d=\xf3\x1d\x00\xe5\
J\x9e\xb5G\x19\xb6+\x8bhZ]\xc0\xee\x14\x84\rC\x91\x85\xf5c\xf6G\xf2f{BGH\xb3\
\xb8n\x08\x06\x8f\x03\xb0\xbc2\xcd_FN\x00\xf4i1\xa5\xcbw\xd5ve\xb1\xe1\x80C\
\x0e`\xf47D4\xa3\x17\xe2\x03}\t\xa5\xcb\xd3l\x953\xbe\x01N\x9d\xbc\xacb\xd1q\
4-\x04\xc0\xad_\xd2\xecw\xa2o\x86\xd1\x0f,7\x9e\xe5\x93VWt!\x08g!<a\xa8\x9a\
\x94\xc2\xb2Bj\xf4\xdc5t\x99\xd8%\xf8\xdd\xed\xecL\x83\xf9t\x93\xcd\xd2<\x86\
\xb1J\xb5\xba\xday\x19\xb6\x13\xd4\xe2\xcaMX\xef\x87x\x16\x0e|\x91P\x1f\x8c]\
C\xd3BXV\x1dM\x0b\x11\x0c\xc4\xdd\xae\xe5J\x9e\xd5?.\xf3\xef\xce}w\x80@P\xb0\
\xaf&\xa5`\xa2w\x11\xc6\xb2\xe7\xd7o\xb7\xbf\xf6\xd5%\xce\xbf\x97c\xc7\xaaq\
\xe3\xa7\x14\x0f\x1e\xda\xd1nW\x16\xddo\xbd\xe4\x82\x98\x1a\x90i5 \xd3\xcau \
\xbe\xe5\xcf\x83^\xf1\xf6\xd4\xb7\xbcyf\x8ab)\xc7\xf2\xca\x14\x00\xa7NNS,\
\xe5\x9cLW\x00\xc1@\x9c\x13C\x97|9P\xae\xe4\xed)\xa8I)\xc2\x13\x86\x8a\x03\
\xde\xd5\xe0\xc3Y|9\xe1@\x97I\x00\x16\xaeL\xa2_H2\x9a\xba\x06\xc0\x9d{3h\x9a\
\xe7\xf369\xb0]\xf6\xac\x02g-\xbbB\x9a\x91\x85\xf5.N\x1c\xfd$\xc5\xe8\x98M\
\xfe\xf3\xed\x8b\xee:w\xd0.\x07\x00D\xb7\xcd\xc8\x0bo\x85$\x0b\x8f?\r\x030\
\x9a\xca\x109\x92\x02\xa0j\xdcga\xe93\x97<l\x18\xeaYt\x1c\xcb2Z\x88\x1d\xec\
\xeb\x89\x1d\xdb!\xe7W\xff<\xc9\xc7\x1f-\xf9\xde\x17K9n\xdd\xbe\xc8as\xa3!8\
\x0b\xc5\xef\xe79x5\xef\x0f\xe0E\x04x\x91\x1c\x99eg\xa7\xc6\xebV\xcdm[\xb82\
\xc9\xceb\xe3\xd9Y]\xf1\xb3\xbb\rY\xda\x8ax!\x01\x83\xd14\xc5\xd2<\x0f7\xda\
\xbc\xf4\x109\x8e1\xd1y\x85\xf5$\xc0\xb2B\xae\xf2C\xc1\x94\x02(\x96\xe6Q\x91\
\xd6\xbe\xeb\xc7h\x1bm|\xab}{\x8b\x00AL\xbdqbF\xbd\xf3\xd6\x9c\x8aE\xa7\xd4\
\xa1`J]\xfcp\x89\x01\x99V\x00G#),\xab\xce\xe3\xadU\x8f\xc0z\xab\x08\x0f\xba\
\xb9\xb0\x1f\xec-3\x12I\xa1\xcb\xd3n\x1d\x07\x88\x1cI\xb1V\xc8\xa0i!\x92\x89\
Y\xee\xcd\xcf\xa9\xc1\xe88\x9b\xa5y\xb7O\xd5Xe!?\t\xdf4\x8d\xec\xd9#Z\\\xf0\
\xb4\xef\x03\xbbH\x0c\x0fM\xa2\xcb\x04\xc5R\x8e\x85\xa5K\x98O7\x01\xbb\x84\
\x02\x04\x031B#\xb3\x04\x031\xca\xe5\xbc\xcb\xb3\xf4\xeb\x0c\xf5\'\x8d\xe4{^\
\x17\xdcB\xb4V\xc8\xf0\xdb\xdd\x19\x02A!\x9e\x9a\xf6\x96W\xae\xe49ln\xb0V\
\xc80<4Ird\x06\xb0\xe7\xffU\xc1\xcd\x01\xd3\xdc \x10\x14\xe2@_B\x9d\x1f\xbbN\
0\x10\xb3\xa3\xcf6\\\x00\xbb\xd8\xb4\x8b\xf8E\xe1:0\x18\x1dg0:\xae\x9c\x1cX^\
\x99\xe6\xe0\xd5\xbcm\xd9\xad\r\x8a_\xe7\x18\x8c\xa6){v8\x07\x07\xefu&\x08\
\x1bF\xd7#\x9b\xeb\x80.\x13\xe82\xc1Z!\xc3\x0f\xb9$\xcf\ns\xbe\xf9\xaa\x1av\
\xd6\xaf=\x9aC\x94{\x8bn\x1d\xdc\xa5\xd7\xe9\xe0\xe3\xcb\x81\x07\x0fg\x05\
\xc0a\xd3\xdfY\xbf\x90\xe4\xd4\xd84\xe5J\x1e\xc5\xa68y\x1c\xf5\xe7?\xf6\xbbn\
\xd1\xc3n\xd2yv\xd1\xf8V\x07\x01\xa6\xe9/kN\xf4\xfa\x85$\xef\x8f]\x07\xe0\
\xce\xdd\x19_\x9fn\xe4n\t\xf6\x9c3\xe2[\xc0\x84\xff>\xb1\xe7\x91\xac\xaf/\
\xcc\xdf\xe6\x06\xcb+_\x12\xb2\xa7A\x91\xdd\x83|7J\xdf9\x83F\x9b\xb7oG\x07\
\x9c\x13\xd2\xf6\x8fynp\x0e\x80\xea^ji\x8dr\xaf\x8bK\xc7\xf3@\xf3]\xa1\x17R\
\x07\xcfse\xebz \xf1\xde\x96\x9a\x85\xbc\x0ci\xcf\x02\xbcB\x9aK\xcf\xab\xba\
\xb2\xff\x07\n$<\xcc\x9c\x11r\xaf\x00\x00\x00\x00IEND\xaeB`\x82_YY\xda' )
    del zlib
    return a

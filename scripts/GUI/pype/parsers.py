#!/usr/bin/python
import time
import os


def detectLineEndings(text):
    crlf_ = text.count('\r\n')
    lf_ = text.count('\n')
    cr_ = text.count('\r')
    mx = max(lf_, cr_)
    if not mx:
        return eol
    elif crlf_ >= mx/2:
        return '\r\n'
    elif lf_ is mx:
        return '\n'
    else:# cr_ is mx:
        return '\r'

def leading(line):
    return len(line)-len(line.lstrip())

def fast_parser(source, line_ending, flat, wxYield):
    lines = source.split(line_ending)
    docstring = {} #new_kwl()
    todo = []
    
    out = []
    stk = []
    line_no = -1
##    SEQ = ('def ','class ')
    
    FIL = lambda A:A[1][2]

    def fun(i, line, ls, line_no, stk):
        try: wxYield()
        except: pass
        na = ls.find('(')
        ds = ls.find(':')
        if na == -1:
            na = ds
        if na != -1:
            if ds == -1:
                ds = na
            fn = ls[len(i):ds].strip()
            if fn:
                lead = len(line)-len(ls)
                while stk and (stk[-1][2] >= lead):
                    prev = stk.pop()
                    if stk: stk[-1][-1].append(prev)
                    else:   out.append(prev)
                nam = i+fn
                nl = nam.lower()
                f = ls[len(i):na].strip()
                stk.append((nam, (f.lower(), line_no, f), lead, []))
                docstring.setdefault(f, []).append(" ".join([fn, '.'.join(map(FIL, stk))]))
    
    for line in lines:
        line_no += 1
        ls = line.lstrip()
        
        #this method is actually the fastest for the
        #single-pass method, but only by ~1%
        #the other versions are easier to maintain
##        if ls[:4] == 'def ':
##            i = 'def '
##            na = ls.find('(')
##            if na == -1:
##                na = ls.find(':')
##            if na != -1:
##                fn = ls[len(i):na].strip()
##                if fn:
##                    lead = len(line)-len(ls)#Leading(None, line, i)
##                    while stk and (stk[-1][2] >= lead):
##                        prev = stk.pop()
##                        if stk: stk[-1][-1].append(prev)
##                        else:   out.append(prev)
##                    nam = i+fn
##                    nl = nam.lower()
##                    docstring[fn] = []
##                    stk.append((nam, (nl, line_no), lead, []))
##        elif ls[:6] == 'class ':
##            i = 'class '
##            na = ls.find('(')
##            if na == -1:
##                na = ls.find(':')
##            if na != -1:
##                fn = ls[len(i):na].strip()
##                if fn:
##                    lead = len(line)-len(ls)#Leading(None, line, i)
##                    while stk and (stk[-1][2] >= lead):
##                        prev = stk.pop()
##                        if stk: stk[-1][-1].append(prev)
##                        else:   out.append(prev)
##                    nam = i+fn
##                    nl = nam.lower()
##                    docstring[fn] = []
##                    stk.append((nam, (nl, line_no), lead, []))

        if ls[:4] == 'def ':
            fun('def ', line, ls, line_no, stk)
        elif ls[:6] == 'class ':
            fun('class ', line, ls, line_no, stk)
        elif ls[:1] == '#':
            a = ls.lower().find('todo:')
            if a+1:
                todo.append((line_no, ls.count('!'), ls[a+5:].strip()))
        #elif ls[:3] == '#>>':
        #    fun('#>>', line, ls, line_no, stk)

##        if ls.startswith('def '):
##            fun('def ', line, ls, line_no, stk)
##        elif ls.startswith('class '):
##            fun('class ', line, ls, line_no, stk)

##        for i in SEQ:
##            if ls[:len(i)] == i:
##                fun(i, line, ls, line_no, stk)
##                break
        #else non-function or non-class definition lines
    while len(stk)>1:
        a = stk.pop()
        stk[-1][-1].append(a)
    out.extend(stk)
    if flat == 0:
        return out, docstring.keys()
    elif flat==1:
        return docstring
    elif flat==2:
        return out, docstring.keys(), docstring
    else:
        return out, docstring.keys(), docstring, todo

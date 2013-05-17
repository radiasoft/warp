"""Run all of the available tests.
This will run all files that have a name with the format *_test.py.
It will print out each file name and whether all of the tests passed or not.
"""
import os
import re

allfiles = os.listdir('.')

# --- Get only the python script files that match the format *_test.py
testfiles = []
for f in allfiles:
    froot,fext = os.path.splitext(f)
    if fext == '.py' and re.match('.+_test.py',f):
        testfiles.append(f)

# --- Run each of the files independently.
for f in testfiles:
    fin,fout,ferr = os.popen3('python %s'%f)
    serr = ferr.readlines()
    if serr[-1] == 'OK\n':
        print('%s OK\n'%f)
    else:
        print('\033[1;31m%s\033[0m'%(f))
        for e in serr[:-1]:
            print(e)
        # --- Use ASCII codes to print this is red
        print('\033[1;31m%s %s\033[0m'%(f,serr[-1]))


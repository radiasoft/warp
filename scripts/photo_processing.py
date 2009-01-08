"""
# Module photo_processing.py
#
# by: Agust Valfells and Lilya Kharevych
# created: Sep. 22, 2000
#
#	Last Modified: Nov. 12, 2001
#
# Additional photo manipulation functions to create tif montages
# and gif animations from previously saved tif photos in WARP
# The module and the functions in it can be invoked anytime after the
# end of a simulation for postprocessing the data.
#
# ==================
# PHOTO PROCESSING:
#   make_mont   ... Assembles tif files at end of run to make montage
# *  make_mov    ... Assembles tif files at end of run to make AVI movie
# ==================
"""
from warp import *
import Image
import os
import string
import numpy
import mphoto
import gifmaker

photo_processing_version = "$Id: photo_processing.py,v 1.4 2009/01/08 19:30:23 dave Exp $"
def photo_processingdoc():
  import photo_processing
  print photo_processing.__doc__

def make_montage(runid=None, label = "z", Rows = None, Columns = None):
    """ make_montage(runid=None, label = "z", Rows = None, Columns = None)
    This function looks in the current directory, assembles all the tif files
prefixed by the runid provided (top.runid by default), and arranges them
in a montage of given dimensions.  If "Rows" and/or "Columns" are not
given, they are calculated to produce the best fit, depending on the number
of photos available.
    The result is saved as a tif file in the same directory, and with the
same name but with "montage" appended.
    """
    if (runid is None):
        runid = arraytostr(top.runid)

    L = []
    T0 = []
    T1 = []
    for file in os.listdir(os.curdir):
        fname = string.split(file, '.')
        if len(fname) == 3:
            if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0] == label):
                L.append(file)
                im = Image.open(file)
                T0.append(im.size[0])
                T1.append(im.size[1])
                del(im)

    print "L =", L
    N_photo = len(L)

    if Rows is None:
        if Columns is None:
            Rows = int(numpy.sqrt(N_photo))
            Columns = (N_photo + (N_photo % Rows)) / Rows
        else:
            if Columns <= N_photo:    Rows =  (N_photo + (N_photo % Columns)) / Columns
            else:
                Columns = 1
                Rows = N_photo
    else:
        Columns = (N_photo + (N_photo % Rows)) / Rows

    D0 = max(T0)
    D1 = max(T1)

    print "D0 =", D0
    print "D1 =", D1
    print "Rows =", Rows
    print "Columns =", Columns

    mphoto.save_tif(numpy.zeros([Columns*D1,Rows*D0],'l'),"background.tif")
    B = Image.open("background.tif")

    print "Background size =", B.size

    L.sort()
    for H in L:
    	im = Image.open(H)
    	n = L.index(H)
    	R = (n / Columns)
    	C = (n - R*Columns)
    	M0 = (D0 - im.size[0]) / 2
    	M1 = (D1 - im.size[1]) / 2
    	B.paste(im,(C*D1 + M1, R*D0 + M0, C*D1 + M1 + im.size[1], R*D0 + M0 + im.size[0]))

    filename = runid + "_" + label + "_montage.tif"
    B.save(filename)
    os.remove("background.tif")


def make_movie(runid=None, num_movie=None, tz=1,ty=0,tx=0):
	""" make_montage(runid=None, num_movie=None, tz=1,ty=0,tx=0):
    This function looks in the current directory, assembles all the tif files
prefixed by the runid provided (top.runid by default), and arranges them
in a gif animation.  <Needs More Work (to change fps, etc)>
	"""
	if runid is None:
		runid = arraytostr(top.runid)

	if tz:
	# ---- Reading z pictures ----
    	# -- Looking for biggest number of picture file name --
		max_num = '00';
		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'z'):
						if (fname[1][1:]>max_num):
							max_num = fname[1][1:]
			except IndexError: continue

		print "maxnum: ", max_num

		# -- Making movie file name --
		m_num='00'
		if num_movie==None:
			for mfile in os.listdir(os.curdir):
				mfname = string.split(mfile, '.')

				try:
					if len(mfname)== 4:
						if (mfname[0]==runid) and (mfname[3]=='gif') and (mfname[2]=='m'):
							if (mfname[1][1:] > m_num):
								m_num = mfname[1][1:]
				except IndexError: continue
        	m_num = `(int(m_num)+1)`

#	else:
#		m_num = `num_movie`

		if len(m_num) < 2:    m_num = '0'+m_num

		filename = runid+'.'+'z'+m_num+'.'+'m'+".gif"

    	# -- Opening needed files in the directory and saving them in frames --
		max_num = `(int(max_num)+1)`
		sequence = []

		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'z'):
						if (fname[1][1:]<max_num) and (fname[1][1:]>'0'):
							try:
								im = Image.open(file)
								sequence.append(im)
							except IOError:
								print "!!Wrong input output"

			except IndexError: continue

		#write GIF animation
		fp = open(filename, "wb")
		gifmaker.makedelta(fp, sequence)
		fp.close()



	if ty:
	# ---- Reading y pictures ----
    	# -- Looking for biggest number of picture file name --
		max_num = '00';
		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'y'):
						if (fname[1][1:]>max_num):
							max_num = fname[1][1:]
			except IndexError: continue

		print "maxnum: ", max_num

		# -- Making movie file name --
		m_num='00'
		if num_movie==None:
			for mfile in os.listdir(os.curdir):
				mfname = string.split(mfile, '.')

				try:
					if len(mfname)== 4:
						if (mfname[0]==runid) and (mfname[3]=='gif') and (mfname[2]=='m'):
							if (mfname[1][1:] > m_num):
								m_num = mfname[1][1:]
				except IndexError: continue
        	m_num = `(int(m_num)+1)`

#		else:
#			m_num = `num_movie`

		if len(m_num) < 2:    m_num = '0'+m_num

		filename = runid+'.'+'y'+m_num+'.'+'m'+".gif"

    	# -- Opening needed files in the directory and saving them in frames --
		max_num = `(int(max_num)+1)`
		sequence = []

		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:	
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'y'):
						if (fname[1][1:]<max_num) and (fname[1][1:]>'0'):
							try:
								im = Image.open(file)
								sequence.append(im)
							except IOError:
								print "!!Wrong input output"
							
			except IndexError: continue
		
		#write GIF animation
		fp = open(filename, "wb")
		gifmaker.makedelta(fp, sequence)
		fp.close()
	
	
	if tx:
	# ---- Reading x pictures ----
 	# -- Looking for biggest number of picture file name --
		max_num = '00';
		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'x'):
						if (fname[1][1:]>max_num):
							max_num = fname[1][1:]
			except IndexError: continue
		
		print "maxnum: ", max_num
	
		# -- Making movie file name --
		m_num='00'
		if num_movie==None:
			for mfile in os.listdir(os.curdir):
				mfname = string.split(mfile, '.')
			
				try:
					if len(mfname)== 4:
						if (mfname[0]==runid) and (mfname[3]=='gif') and (mfname[2]=='m'):
							if (mfname[1][1:] > m_num):
								m_num = mfname[1][1:]
				except IndexError: continue
        	m_num = `(int(m_num)+1)`

	#	else:	
	#		m_num = `num_movie`
		
		if len(m_num) < 2:    m_num = '0'+m_num

		filename = runid+'.'+'x'+m_num+'.'+'m'+".gif"

    	# -- Opening needed files in the directory and saving them in frames --
		max_num = `(int(max_num)+1)`
		sequence = []

		for file in os.listdir(os.curdir):
			fname = string.split(file,'.')
			try:
				if len(fname)== 3:
					if (fname[0]==runid) and (fname[2]=='tif') and (fname[1][0]== 'x'):
						if (fname[1][1:]<max_num) and (fname[1][1:]>'0'):
							try:
								im = Image.open(file)
								sequence.append(im)
							except IOError:
								print "!!Wrong input output"
							
			except IndexError: continue
		
		#write GIF animation
		fp = open(filename, "wb")
		gifmaker.makedelta(fp, sequence)
		fp.close()

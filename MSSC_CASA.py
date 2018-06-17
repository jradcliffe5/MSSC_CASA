#### MSSC_CASA v0.1 ###
### If you use MSSC in your work please cite either Radcliffe+16 ##
### Or the new version which will be in Moldon, Radcliffe et al. +18 ###

### Imports ###
import numpy as np
import sys
import os
import shutil
import stat
import time
from math import sqrt
from taskinit import *


### Inputs ###
### TO do: make a input file ###

vis = ['1252+5634.ms']


## Imaging options
manual = False

###

def uvdiv(vis):
    ## Somehow try to parallelise this?
    ## Not sure how to do that apart from using TaQL or the casa mpi?
    t = tbtool()
    casalog.post('Dividing visibilities by MODEL', 'INFO')

    mycb = cbtool()
    if (type(vis) == str) & os.path.exists(vis):
        # add CORRECTED_DATA column
        mycb.open(filename=vis, compress=False, addcorr=True,
                  addmodel=False)
    else:
        raise Exception, \
            'Visibility data set not found - please verify the name'
    mycb.close()


    t.open(vis, nomodify=False)
    for colname in ['DATA']:
        if (colname in t.colnames()) and (t.iscelldefined(colname,0)):
            for j in xrange(0,t.nrows()):
                a = t.getcell(colname, j)
                model = t.getcell('MODEL_DATA',j)
                a = a/model
                t.putcell('CORRECTED_DATA', j, a)
    t.close()

def separate_sources(msfile):
    ### Generate a copy of these data
    for i in range(len(msfile)):
        os.system('rm -r %s_temp_MSSC%d.ms' % (msfile[i],i))
        os.system('cp -r %s %s_temp_MSSC%d.ms' % (msfile[i],msfile[i],i))
        t = tbtool()
        t.open('%s_temp_MSSC%d.ms/FIELD' % (msfile[i],i), nomodify=False)
        x = t.getcol('NAME')
        if len(x) > 1:
            print('Number of fields needs to be 1 per measurement set')
            exit()
        else:
            x = 'MSSC_%d' % i
        t.putcell('NAME', 0, x)
        t.close()

def initial_image(msfile):
    for i in range(len(msfile)):
        vis = '%s_temp_MSSC%d.ms' % (msfile[i],i)
        tclean(vis=vis, )
### MSSC steps ###
### 1. Generate initial image ###
#separate_sources(vis)
#initial_image(vis)
uvdiv('1252+5634.ms_temp.ms')

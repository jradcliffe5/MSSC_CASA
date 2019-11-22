#### MSSC_CASA v0.1 ###
### If you use MSSC in your work please cite either Radcliffe+16 ##
### Or the new version which will be in Moldon, Radcliffe et al. +18 ###

### Imports ###
import numpy as np
import numpy.ma as ma
import sys
import os
import shutil
import stat
import time
from math import sqrt
from taskinit import *


### Inputs ###
### TO do: make a input file ###

msfile = ['evn_SRC0.ms','evn_SRC1.ms']
chanaverage = 1 # channels
timeaverage= 0 # seconds

## Imaging options
manual = False

###

def uvdiv(vis):
    ## Somehow try to parallelise this?
    ## Not sure how to do that apart from using TaQL or the casa mpi?
    t = tbtool()
    casalog.post('Dividing visibilities by MODEL', 'INFO')

    t.open(vis, nomodify=False)
    ram_restrict = 100000
    ranger = list(range(0,t.nrows(),ram_restrict))
    
    for colname in ['DATA']:
        if (colname in t.colnames()) and (t.iscelldefined(colname,0)):
            for j in ranger:
                if j == ranger[-1]:
                    ram_restrict = t.nrows()%ram_restrict
                a = t.getcol('DATA',startrow=j, nrow=ram_restrict, rowincr=1)
                model = t.getcol('MODEL_DATA',startrow=j, nrow=ram_restrict, rowincr=1)
                model_mask = model == 0
                model = ma.array(data=model,mask=model_mask)
                #model = ma.array(data=model,mask=np.isin(model,0))
                ### note that the use of np.isin cannot be used in CASA as it still uses numpy version v.1.11 whereas this came in v1.14
                ### Had to use old horrible version (see line 51)
                a = a/model
                t.putcol('CORRECTED_DATA',a.filled(0),startrow=j, nrow=ram_restrict, rowincr=1)
                #t.putcell('CORRECTED_DATA', j, a.filled(0))
    
    '''
    for colname in ['DATA']:
        if (colname in t.colnames()) and (t.iscelldefined(colname,0)):
            for j in xrange(0,t.nrows()):
                a = t.getcell(colname, j)
                model = t.getcell('MODEL_DATA',j)
                model_mask = model == 0
                model = ma.array(data=model,mask=model_mask)
                #model = ma.array(data=model,mask=np.isin(model,0))
                ### note that the use of np.isin cannot be used in CASA as it still uses numpy version v.1.11 whereas this came in v1.14
                ### Had to use old horrible version (see line 51)
                a = a/model
                t.putcell('CORRECTED_DATA', j, a.filled(0))
    '''
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

def initial_image(msfile,datacolumn='data',position=''):
    ## First one is to get the rms values, we use the fit-half algorithm
    if datacolumn=='corrected':
        appendix = '_uvdiv'
    else:
        appendix = ''
    os.system('rm -r %s_dirty%s.*'%(msfile,appendix))
    tclean(vis=msfile,
           deconvolver='clark',
           datacolumn=datacolumn,
           imagename='%s_dirty%s'%(msfile,appendix),
           imsize=[2548,2548],
           cell='0.001arcsec',
           niter=0)
    threshold = 3.5*imstat(imagename='%s_dirty%s.image'%(msfile,appendix),algorithm='fit-half')['rms'][0]
    os.system('rm -r %s_IM%s.*'%(msfile,appendix))
    tclean(vis=msfile,
           deconvolver='clark',
           imagename='%s_IM%s'%(msfile,appendix),
           imsize=[2548,2548],
           datacolumn=datacolumn,
           cell='0.001arcsec',
           gain=0.05,
           niter=10000,
           phasecenter=position,
           threshold=threshold,
           savemodel='modelcolumn')

def adjust_phase_centre(ms,position):
    tb.open('%s/FIELD'%ms,nomodify=False)
    for i in ['DELAY_DIR','PHASE_DIR','REFERENCE_DIR']:
        pdir = tb.getcol(i)
        pdir[0] = [[position[0]]]
        pdir[1] = [[position[1]]]
        tb.putcol(i,pdir)
    tb.close()

def add_columns(ms):
    mycb = cbtool()
    if (type(ms) == str) & os.path.exists(ms):
        # add CORRECTED_DATA column
        mycb.open(filename=ms, compress=False, addcorr=True,
                  addmodel=True)
    else:
        raise Exception, \
            'Visibility data set not found - please verify the name'
    mycb.close()
### MSSC steps ###
### 1. Generate initial image ###
#separate_sources(vis)
#initial_image(vis)
for i, ms in enumerate(msfile):
    add_columns(ms)
    initial_image(msfile=ms,datacolumn='data')
    uvdiv(ms)
    initial_image(msfile=ms,datacolumn='corrected')
    os.system('rm -r MSSC_%s.ms'%i)
    split(vis=ms,outputvis='MSSC_%s.ms'%i, width=chanaverage,timebin='%ss'%timeaverage)
    adjust_phase_centre('MSSC_%s.ms'%i,[0,1.04])

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
           phasecenter=position,
           cell='0.001arcsec',
           niter=0,
           pblimit=1e-10,
           parallel=True)
    threshold = 1.0*imstat(imagename='%s_dirty%s.image'%(msfile,appendix),algorithm='fit-half')['rms'][0]
    os.system('rm -r %s_IM%s.*'%(msfile,appendix))
    tclean(vis=msfile,
           deconvolver='clark',
           imagename='%s_IM%s'%(msfile,appendix),
           imsize=[2548,2548],
           datacolumn=datacolumn,
           cell='0.001arcsec',
           pblimit=1e-10,
           gain=0.05,
           niter=10000,
           phasecenter=position,
           threshold=threshold,
           savemodel='modelcolumn',
           parallel=True)

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

ms='VLBA_SRC005_sp.ms'
initial_image(msfile=ms,datacolumn='data',position=phasecenter[i])
uvdiv(ms)
initial_image(msfile=ms,datacolumn='corrected')
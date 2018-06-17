import numpy as np
import sys
import os
import shutil
import stat
import time
from math import sqrt
import numpy.ma as ma
from taskinit import *

vis = '1252+5634.ms'

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
            model = ma.array(data=model,mask=np.isin(model,0))
            a = a/model
            t.putcell('CORRECTED_DATA', j, a.filled(0))
t.close()

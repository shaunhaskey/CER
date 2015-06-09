import numpy as np
import scipy.io.netcdf as netcdf
import os
import OMFITtree
import matplotlib.pyplot as pt
import numpy
import scipy
from scipy.optimize import curve_fit
import cer_funcs as CER

#  params=MPFITFUN('MTANH_MPFIT_EVALUATE',x,y,err,guess,$
#                  FUNCTARGS={NCORE:ncore,$
#                             NEDGE:nedge})
#  yfit=MTANH_MPFIT_EVALUATE(x,params,NCORE=ncore,NEDGE=nedge)

#for ident, width_mult, top_val in zip(['ne', 'te', 'ti'], [1., 1., 2.],[4.,3.,3.]):
offset = 0.1
xsym = 0.95
core = [0.02, 0.001]
core = [0.02]#, 0.001]
edge = [-0.02]
width = 0.1
npts = 121
max_val = 1.21
top_val = 4.
npts=121;max_val=1.21;x = np.linspace(0,max_val, npts)

y = CER.mtanh(x, top_val, offset, xsym, width, len(core), len(edge), *(core + edge))

fig, ax = pt.subplots()
ax.plot(x,y)
#guess = [top_val,offset,xsym,width]
guess = [1.5, 0.2, 0.9, 0.05,]
def mtanh_wrapper(x,*p):
    core = [0.02, 0.001]
    core = [0.02]#, 0.001]
    edge = [-0.02]
    p = [i for i in p]
    p = p + [len(core), len(edge)] + core + edge
    print p
    y = CER.mtanh(x,*p)
    return y

coeff, var_matrix = curve_fit(mtanh_wrapper, x, y, p0=guess,)
y2 = mtanh_wrapper(x, *coeff)
ax.plot(x,y2)
fig.canvas.draw();fig.show()

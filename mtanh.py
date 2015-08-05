import numpy as np
import matplotlib.pyplot as pt
offset = 0.1
top = 3.0
xsym = 0.95
hwid = 0.05
core = [0.02, 0.001]
core = [0.02, 0.001]
edge = [-0.02]

from scipy.io import netcdf
f = netcdf.netcdf_file('test.nc','w')


def mtanh(top, offset,xsym,hwid, core = None, npts=121, max_val=1.21, edge = None):
    '''Copied from /scp:venus:/u/grierson/idlpros/utilities/mpfit_mtanh.pro

    SRH : 02June2015
    '''
    x = np.linspace(0,max_val, npts)
    #x = np.arange(121)/100.
    a = (top-offset)/2.0
    b = (3*offset+top)/2.0
    z = (xsym-x)/hwid
    if core==None: core = []
    if edge==None: edge = []
    if len(core)>0:
        poly_core = 1.0
        for i in range(len(core)):
            poly_core += core[i] * z**(i+1)
        num_core = poly_core*np.exp(z)

    #;; Edge with a negative slope
    if len(edge)>0:
        poly_edge = 1.0
        for i in range(len(edge)):
            poly_edge += edge[i] * z**(i+1)
        num_edge = poly_edge*np.exp(-z)

        num = num_core - num_edge
        den = np.exp(z) + np.exp(-z)
        y = a * num/den + b
    return x, y
fig, ax = pt.subplots()
for i, hwid in enumerate(np.linspace(0.05,0.15,5)):
    x, y = mtanh(top, offset,xsym, hwid, core = core, npts=121, max_val=1.21, edge = edge)
    ax.plot(x,y)
    ax.axhline(offset)
    ax.axhline(top)
    ax.axvline(xsym)
    ax.axvline(xsym+hwid)
    ax.axvline(xsym-hwid)
    length = len(x)
    if i==0:
        f.createDimension('rho',length)
        rho = f.createVariable('rho','float',('rho',))
        rho[:] = +x
        rho.units = ''
    cur_var = 'ne{}'.format(i)
    f.createDimension(cur_var,length)
    prof = f.createVariable(cur_var,'float',(cur_var,))
    prof[:] = +y
f.close()

fig.canvas.draw();fig.show()
1/0

# Y(Z) = A*MTANH(Z,core,edge) + B

# MTANH = ((1+core[0]*Z + core[1]*Z^2 + ...)*np.exp(Z) - (1+edge[0]*Z + edge[1]*Z^2 + ...)*np.exp(-Z)) / (np.exp(Z) + np.exp(-Z))


# Z = (Xsym-X)/Hwid
# Pedestal = A+B
# Offset = B-A
# Width = 2*Hwid
# Knee = Xsym-Hwid
# Foot = Xsym+Hwid

# x = FINDGEN(121)/100.
# offset = 0.1
# top = 3.0
# #For a reasonable set of parameters, the best guess for
# #the parameters A and B is that B-A is the min and A+B is the
# #value at the top of the pedestal, so b-a = min and a+b=top
# #so a+(min +a) = top, or 
# #    a = (top-min)/2
# # similarly b = a+min = (3min+top)/2
# a = (top-offset)/2.0
# b = (3*offset+top)/2.0

# #Symmetry and half-width
# xsym = 1.0
# hwid = 0.05

# #;; Core with a negative slope and curvature
# core = [0.02, 0.001]

# #;; Edge with a negative slope
# edge = [-0.02]

# #Evaluate the function
# y = MTANH_MPFIT_EVALUATE(x,[a,b,xsym,hwid,core,edge],$
#                            NCORE=N_ELEMENTS(core),$
#                            NEDGE=N_ELEMENTS(edge))

# #Now fit the data
# params = MTANH_MPFIT(x,y,0.1*y,[a,b,xsym,hwid,core,edge],$
#                        NCORE=N_ELEMENTS(core),$
#                        NEDGE=N_ELEMENTS(edge))
# print,params
# BAG_PLOT_SETUP,/DIR,WIN=0
# PLOT,x,y

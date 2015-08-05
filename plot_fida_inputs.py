import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

from mpl_toolkits.mplot3d import Axes3D

loc = '/home/shaskey/01010/MAIN_ION330/'
fname = 'def_inputs.cdf'
net = netcdf.netcdf_file(loc + fname,mmap=False)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(net.variables['x_grid'][0,0,:],net.variables['y_grid'][0,:,0],net.variables['z_grid'][:,0,0])
fig.canvas.draw();fig.show()

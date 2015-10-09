import scipy.io.netcdf as net
import numpy as np
import matplotlib.pyplot as pt


shot = 900000
a = net.netcdf_file('{}.nc'.format(shot))
fig, ax = pt.subplots()
#a = netcdf_file
d = a.variables['m01'].data
im = ax.imshow(d,aspect='auto',cmap='spectral', interpolation = 'nearest')
ax.set_title('m01')
im.set_clim([0,16000])
print im.get_clim()
for i in range(0,728,128): ax.axvline(i,color='k')
#im.set_clim([0,4000])
fig.canvas.draw();fig.show()

fig, ax = pt.subplots()
mult = [3,0,-3,6]
labels = ['+3dB','0dB','-3dB','+6dB']
for shot,lab in zip([900000,900001,900002,900004], labels):
    a = net.netcdf_file('{}.nc'.format(shot))
    d = a.variables['m01'].data
    ax.plot(np.mean(d[0:20,:],axis=0), label = lab)

# a = net.netcdf_file('900000.nc'.format(shot))
# d = a.variables['m01'].data*0.707
# ax.plot(d[20,:])

# a = net.netcdf_file('900002.nc'.format(shot))
# d = a.variables['m01'].data/0.707
# ax.plot(d[20,:])

a = net.netcdf_file('163601.nc'.format(shot))
d = a.variables['m01white'].data
ax.plot(d, label='Bad Data 163601')

a = net.netcdf_file('163600.nc'.format(shot))
d = a.variables['m01white'].data
ax.plot(d, label='Good data 163600')

#ax.set_ylim([0,10000])
ax.set_ylim([6000,8000])
ax.set_xlim([0,769])
for i in range(0,728,128): ax.axvline(i,color='k')
ax.legend(loc='best')
fig.suptitle('m01 Bad white light calibration tests')
fig.savefig('m01_bad_whitecal_tests2.png')
fig.canvas.draw();fig.show()

import pidly, os, time
import numpy as np
import matplotlib.pyplot as pt
import gc
import scipy.io.netcdf as net
#idl = pidly.IDL('cerview')
#shot = 162288
#chord = 'm13'
#idl("get_chord,{},'{}'".format(shot, chord))
#idl("openw,1,'/u/haskeysr/hello.txt'")
#idl("printf,1,chord_data.data")
#idl("close,1")

#a = np.loadtxt('/u/haskeysr/hello.txt')


def get_data(chord_list, shot, idl = None, raw = True, net_cdf = True):
    if idl==None:idl = pidly.IDL('cerview')
    for i in range(len(chord_list)):
        os.system('rm /u/haskeysr/hello{}.txt'.format(i))
    for i in range(len(chord_list)):
        while os.path.isfile('/u/haskeysr/hello{}.txt'.format(i)):
            print 'waiting for file to go', i
            time.sleep(1)
    if net_cdf:
        idl("id=ncdf_create('{}.nc',/CLOBBER)".format(shot))
        for i,chord in enumerate(chord_list):
            idl("get_chord,{},'{}'".format(shot, chord))
            print 'm{:02d}'.format(chord)
            idl("xid = ncdf_dimdef(id,'x',768)".format(chord))
            idl("yid = ncdf_dimdef(id,'y',2032)".format(chord))
            idl("vid = ncdf_vardef(id,'{}',[xid,yid],/LONG)".format(chord))
            idl("ncdf_varput, id, 'm{}', chord_data.raw".format(chord))
        idl("ncdf_close, id")
    else:
        for i,chord in enumerate(chord_list):
            idl("get_chord,{},'{}'".format(shot, chord))
            idl("openw,1,'/u/haskeysr/hello{}.txt'".format(i))
            data = 'raw' if raw else 'data'
            idl("printf,1,chord_data.{}".format(data))
            idl("close,1")

    for i in range(len(chord_list)):
        while not os.path.isfile('/u/haskeysr/hello{}.txt'.format(i)):
            print 'waiting for file to exist', i
            time.sleep(1)
    time.sleep(2)
    print 'loading file'
    output_data = []
    for i in range(len(chord_list)):
        a = np.loadtxt('/u/haskeysr/hello{}.txt'.format(i))
        if raw:
            d = a.flatten().reshape((768,2032), order = 'F')
        else:
            d = a.flatten().reshape((768,2000), order = 'F')
        output_data.append(d.copy())
    return idl, output_data

shot = 162288
def write_shot_nc(shot):
    ch_list = [1,2,3,4,5,6,7,8,10,11,12,13,17,24,25,26,27,28,29,31]
    #chord_list = ['m10','m11','m12','m13']
    chord_list = ['m{:02d}'.format(i) for i in ch_list]
    cmd = "id=ncdf_create('{}.nc',/CLOBBER)".format(shot); print cmd; idl(cmd)
    cmd = "xid = ncdf_dimdef(id,'x',768)"; print cmd; idl(cmd)
    cmd = "yid = ncdf_dimdef(id,'y',2032)"; print cmd; idl(cmd)
    for chord in chord_list:
        cmd = "vid = ncdf_vardef(id,'{}',[xid,yid],/FLOAT)".format(chord); print cmd; idl(cmd)
    cmd = "NCDF_CONTROL, id, /ENDEF"; print cmd; idl(cmd)
    for chord in chord_list:
        cmd = "get_chord,{},'{}'".format(shot, chord);print cmd; idl(cmd)
        cmd = "ncdf_varput, id, '{}', chord_data.raw".format(chord); print cmd; idl(cmd)
    cmd = "ncdf_close, id"; print cmd; idl(cmd)


def clr_plot(system, shot, netcdf_file = None):
    if system=='MU1': chord_list = ['m{:02d}'.format(i) for i in [1,2,3,4]]
    if system=='MU2': chord_list = ['m{:02d}'.format(i) for i in [5,6,7,8]]
    if system=='MU3': chord_list = ['m{:02d}'.format(i) for i in [10,11,12,13]]
    if system=='MU4': chord_list = ['m{:02d}'.format(i) for i in [17,24,25,26,27,28,29,31]]
    if len(chord_list)==8:
        fig, ax = pt.subplots(ncols = 4, nrows = 2)
        ax = ax.flatten()
    else:
        fig, ax = pt.subplots(ncols = 4)
        ax = ax.flatten()
    if netcdf_file == None:
        a = net.netcdf_file('{}.nc'.format(shot))
    else:
        a = netcdf_file
    for ax_tmp, chord in zip(ax,chord_list):
        d = a.variables[chord].data
        im = ax_tmp.imshow(d,aspect='auto',cmap='spectral', interpolation = 'nearest')
        ax_tmp.set_title(chord)
        for i in range(0,728,128): ax_tmp.axvline(i,color='k')
        #ims.append(im)
        im.set_clim([0,100])
    fig.suptitle('{}'.format(shot))
    #fig.canvas.draw();fig.show()
    fig.savefig('{}_{}.png'.format(system, shot))
    fig.clf()
    pt.close()
    gc.collect()


for shot in range(162330,162348):
    a = net.netcdf_file('{}.nc'.format(shot))
    clr_plot('MU1', shot, netcdf_file = a)
    clr_plot('MU2', shot, netcdf_file = a)
    clr_plot('MU3', shot, netcdf_file = a)
    clr_plot('MU4', shot, netcdf_file = a)

for shot in range(162295,162319):
    a = net.netcdf_file('{}.nc'.format(shot))
    clr_plot('MU1', shot, netcdf_file = a)
    clr_plot('MU2', shot, netcdf_file = a)
    clr_plot('MU3', shot, netcdf_file = a)
    clr_plot('MU4', shot, netcdf_file = a)

1/0
for shot in range(162330,163348):
    write_shot_nc(shot)

a = net.netcdf_file('162288.nc',mmap=False)
print a.variables

1/0

raw = True
master_dict = {}

ch_list = [1,2,3,4,5,6,7,8,10,11,12,13,17,24,25,26,27,28,29,31]
#chord_list = ['m10','m11','m12','m13']
chord_list = ['m{:02d}'.format(i) for i in ch_list]

idl, output_data = get_data(chord_list, shot, idl = idl, raw = raw)
master_dict[shot]={}
for chord,dat in zip(chord_lis,output_data):
    master_dict[shot][chord] = dat.copy()
pickle.dump


#c = a.flatten().reshape((768,2000))

#d = a.flatten().reshape((768,2000), order = 'F')
import matplotlib.pyplot as pt
fig, ax = pt.subplots(ncols = len(chord_list), sharex = True, sharey = True)
ims = []
for ax_tmp, d in zip(ax,output_data):
    im = ax_tmp.imshow(d.transpose(),aspect='auto',cmap='spectral', interpolation = 'nearest')
    for i in range(0,728,128): ax_tmp.axvline(i,color='k')
    ims.append(im)
    if raw:
        im.set_clim([0,100])
    else:
        im.set_clim([0,50])
fig.canvas.draw();fig.show()

def change_clim(val):
    for i in ims:
        i.set_clim(val)
    fig.canvas.draw()
    fig.tight_layout()

import scipy.signal as sig
sig.correlate(d[0],d[1])
fig, ax = pt.subplots()
col = ['b','k','r','y']
for d,c in zip(output_data,col): ax.plot(np.sum(d,axis=0),'-{}'.format(c))
ax.plot(np.sum(output_data[-1] + output_data[-2],axis=0)/2., '--')
fig.canvas.draw();fig.show()

fig2,ax2 =pt.subplots();
for i in range(20,30):
    ax2.plot(output_data[2][i,:]);
#ax2.plot(output_data[2][180,:]);
fig2.canvas.draw();fig2.show()

fig3,ax3 =pt.subplots(ncols = 2, sharey = True);
im_list  = []
im = ax3[0].imshow(d.transpose()[:,0:66],aspect='auto',cmap='spectral', interpolation = 'nearest')
im_list.append(im)
im.set_clim([0,200])
im = ax3[1].imshow(d.transpose()[:,370:430],aspect='auto',cmap='spectral', interpolation = 'nearest')
im_list.append(im)
im.set_clim([0,200])
fig3.canvas.draw();fig3.show()
#ax2.plot(output_data[2][180,:]);

1/0
idl.close()
import scipy.io.netcdf


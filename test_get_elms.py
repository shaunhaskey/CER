import os
import numpy as np
import data as MDSdata
os.environ['VPN_ACTIVE']='yes'
import data as MDSdata
import matplotlib.pyplot as pt

#input data: shot, fs, thresh,plus_dt, minus_dt
#time_window
#buttons: plot, write
shot = 160409
channel = 'fs04f'
dat = MDSdata.Data(channel,shot)
#return dat.x[0], dat.y
# if plot:
#fs04 = data.Data('fs04',160820)
#'\pedestal::top.elm_times.'+fs+'.elmsmarker'
fig, ax = pt.subplots()
ax.plot(dat.x[0], dat.y,)
dat = []
for i,style,mult in zip(['elmsmarker','elmpeak','elmstart','elmend'],['o','x','d','s'], [0,1,0,0]):
    dat.append(MDSdata.Data('\pedestal::top.elm_times.{}.{}'.format(channel,i),shot))
    ax.plot(dat[-1].x[0], dat[-1].y*mult,style)
simple_x= np.zeros(dat[2].x[0].shape[0]*3,dtype=float)
simple_y= np.zeros(dat[2].x[0].shape[0]*3,dtype=float)
simple_x[0::3] = dat[2].x[0]
simple_x[1::3] = dat[1].x[0]
simple_x[2::3] = dat[3].x[0]
simple_y[0::3] = dat[2].y*0
simple_y[1::3] = dat[1].y
simple_y[2::3] = dat[3].y*0
ax.plot(simple_x, simple_y,'-')

#MDSdata.Data('\electrons::tstime_core',160409)
ax.plot(simple_x, simple_y,'-')
tstime_core = MDSdata.Data('tstime_core',160409)
tstimes = tstime_core.x[0]
remove = np.zeros(len(tstimes),dtype=bool)
thresh =  6.9e14
thresh =  2.64e16
ax.axhline(thresh)
plus_time = 3.
minus_time = 1.
for start,peak,end in zip(dat[2].x[0],dat[1].y,dat[3].x[0]):
    cur_vals = ((tstimes>(start-minus_time)) * (tstimes<(end+plus_time)) * (peak>thresh))
    remove = remove + cur_vals
    print np.sum(cur_vals),np.sum(remove)

for x in tstimes[remove]:
    ax.axvline(x,color='r')
for x in tstimes[np.invert(remove)]:
    ax.axvline(x,color='k')
#    ax.axvline(tstimes[good_bad],'k')
ax.set_xlim([2900,3050])
fig.canvas.draw(); fig.show()

  # elmsmarker=mdsvalue('\pedestal::top.elm_times.'+fs+'.elmsmarker')
  # elmsmarker_t=mdsvalue('dim_of(\pedestal::top.elm_times.'+fs+'.elmsmarker)')

  # elmpeak=mdsvalue('\pedestal::top.elm_times.'+fs+'.elmpeak')
  # elmpeak_t=mdsvalue('dim_of(\pedestal::top.elm_times.'+fs+'.elmpeak)')      

  # elmstart=mdsvalue('\pedestal::top.elm_times.'+fs+'.elmstart')
  # elmstart_t=mdsvalue('dim_of(\pedestal::top.elm_times.'+fs+'.elmstart)')

  # elmend=mdsvalue('\pedestal::top.elm_times.'+fs+'.elmend')
  # elmend_t=mdsvalue('dim_of(\pedestal::top.elm_times.'+fs+'.elmend)')
'/u/smithsp/Komodo8/bin/komodo'

#for tmp_dat in ['R','Z','CDPOLYBOX', 'LFORDER', 'DENSITY', 'DENSITY_E', 'TEMP', 'TEMP_E', 'REDCHISQ']:
thom_dat = {}
for tmp_dat in ['R','Z', 'DENSITY', 'DENSITY_E', 'TEMP', 'TEMP_E']:
    rev = 'REVISION02'
    #rev = MDSdata.Data('\TOP.TS:BLESSED_ID',shot)
    node = r'\electrons::top.TS.REVISIONS.{rev}.{branch}.{dat}'.format(rev=rev, branch='CORE', dat=tmp_dat)
    print node
    thom_dat[tmp_dat] = MDSdata.Data(node,shot)
#    mdssetdefault,'\TOP.TS.REVISIONS.'+rev+'.'+branch,quiet=quiet,status=status
#    infonodes = {time:'TIME', chan:'CHANNEL', cdpulse:'CDPULSE'}
#    nodes = 
#  endelse
Z_array = thom_dat['DENSITY'].y*0
for tmp in range(thom_dat['Z'].y.shape[0]):
    Z_array[tmp,:] = thom_dat['Z'].y[tmp]

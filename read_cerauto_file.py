import time
import cer_funcs as CER
import os
import matplotlib.pyplot as pt
import numpy as np
shot = 160414

print ''
shot_imp = []
#for shot in range(160403,160423):
for shot in range(160818,160840):
    shot_data = CER.read_shot_status(shot)
    for tmp_shot in [shot]:#shot_data.keys():
        print '#### {} ####'.format(tmp_shot)
        for i in np.sort(shot_data.keys()):
            print '{:8s} {:15s} {:6s} {}'.format(i, shot_data[i]['TIMI'], shot_data[i]['WAVE'], shot_data[i]['CHOR'])
        if shot_data['CT1']['WAVE']=='5290.5' and shot_data['Unit1']['WAVE']=='5167.0':shot_imp.append(tmp_shot)

1/0

chord = 't01'
for i in shot_data.keys():
    if shot_data[i]['CHOR'].lower().find(chord)>=0:
        spec = i
print spec, shot_data[spec]['CHOR']
stri = shot_data[spec]['TIMI']
print stri
#stri = '0.0:2000@5.0/1'
start_time, n_slices, t_int, start_pts, end_pts = CER.split_timing_string(stri)
nbi_x, nbi_y = CER.get_nbi_data('30l', plot = False)
#dat.x[0], dat.y
fig, ax = pt.subplots()
ax.plot(nbi_x,nbi_y/np.max(nbi_y))
fig.canvas.draw();fig.show()
ax.plot(start_pts, start_pts*0+1,'o')
ax.plot(end_pts, end_pts*0+1,'x')
ax.set_ylim([0,1.5])
fig.canvas.draw();fig.show()

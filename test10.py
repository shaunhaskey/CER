import numpy as np
lens1 = np.array([-58.0452,238.6632,0.6822])
loc1 = np.array([-133.64747, 172.84416, -1.4165587])
lens2 = np.array([-58.045200, 238.66320, 0.68220000])
loc2 = np.array([-134.95745, 186.51464, -1.0956602])
ids = ['m17','m24','m25','m26','m27','m28','m29','m31','m32']
n = len(ids)
lens = [];loc = []
axes = ['x','y','z']
for i in range(3):
    lens.append(np.linspace(lens1[i],lens2[i],n))
    loc.append(np.linspace(loc1[i],loc2[i],n))
for i in range(3):
    print '{}lens = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in lens[i]]))
for i in range(3):
    print '{}los = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in loc[i]]))

import numpy as np
import matplotlib.pyplot as pt
xlens = np.array([-58.045200,-58.045200,-58.045200,-58.045200])
ylens = np.array([238.66320, 238.66320, 238.66320, 238.66320])
zlens = np.array([0.68220000, 0.68220000, 0.68220000, 0.68220000])
################
xlos = np.array([-133.64747, -134.55281, -134.71830, -134.95745])
ylos = np.array([172.84416, 182.29196, 184.01889, 186.51464])
zlos = np.array([-1.4165587, -1.2403410, -1.1361711, -1.0956602])
fig, ax = pt.subplots()
ax.plot(xlens,ylens,'x')
ax.plot(xlos,ylos,'o')
fig.canvas.draw();fig.show()

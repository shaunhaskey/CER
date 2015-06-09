import numpy as np
import cer_funcs as CER
import matplotlib.pyplot as pt
x = np.linspace(0,1.21,100)
offset = 0.1
top = 3
sym = 0.95
width = 0.05
core = [0.02]#, 0.001]
edge = []#-0.02]
#core = []
#edge = []
ncore = len(core)
nedge = len(edge)
fig, ax = pt.subplots()
def mtanh_plot_wrap():
    y = CER.mtanh(x,top,offset,sym,width,ncore,nedge,False,*(core + edge))
    ax.plot(x,y)
for top in np.linspace(0.5,4,5):
    for width in np.linspace(0.05,0.15,10):
        mtanh_plot_wrap()
fig.canvas.draw();fig.show()

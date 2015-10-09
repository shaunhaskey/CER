import matplotlib.pyplot as pt
fig, ax = pt.subplots()
import pmds
X = pmds.mdsvalue('ptdata2("cer3tang01",163257,0)')
X = pmds.mdsvalue('ptdata2("cer3mion01",163257,0)')

Y = X.reshape(2032,X.shape[0]/2032)
im = ax.imshow(Y, aspect = 'auto',interpolation = 'nearest',cmap='spectral')
im.set_clim([0,3000])
fig.canvas.draw();fig.show()
import Ptdata
Ptdata.ptgetheader("cer3tang01",163257)
print Ptdata._real64[:20]
print Ptdata._int16[:20]

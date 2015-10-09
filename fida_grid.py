import matplotlib.pyplot as pt
import numpy as np
import scipy.io.netcdf as netcdf

directory = '/u/haskeysr/FIDASIM/RESULTS/D3D/158676/01000/MAIN_ION330/'
id_name = 'def'
def fida_plot_grid(directory, id_name, RMIDOUT = 2.253, RMIDIN=1.1043):
    fig, ax = pt.subplots()
    inputs = netcdf.netcdf_file('{}/{}_inputs.cdf'.format(directory, id_name), mmap = False)
    for x in inputs.variables['x_grid'].data[0,0,:]:
        ax.axvline(x)
    for y in inputs.variables['y_grid'].data[0,:,0]:
        ax.axhline(y)
    yrange = [np.min(inputs.variables['y_grid'].data[0,:,0]), np.max(inputs.variables['y_grid'].data[0,:,0])]
    xrange = [np.min(inputs.variables['x_grid'].data[0,0,:]), np.max(inputs.variables['x_grid'].data[0,0,:])]
    circle2=pt.Circle((RMIDOUT*100,RMIDOUT*100),.5,color='b',fill=False)
    circle2=pt.Circle((RMIDIN*100,RMIDIN*100),.5,color='b',fill=False)
    for i, (x1,y1, x2, y2) in enumerate(zip(inputs.variables['xlos'].data, inputs.variables['ylos'].data, inputs.variables['xlens'].data, inputs.variables['ylens'].data)):
        poly_vals = np.polyfit([x1,x2],[y1,y2],1)
        poly_vals = np.polyfit([y1,y2],[x1,x2],1)
        y3  = -20
        x3 = np.polyval(poly_vals,y3)
        if i==60:
            ax.plot([x1,x2],[y1,y2],'k-')
            ax.plot([x2,x3],[y2,y3],'k-')
        else:
            ax.plot([x1,x2],[y1,y2],'b-')
            ax.plot([x2,x3],[y2,y3],'b-')
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    fig.canvas.draw();fig.show()

fida_plot_grid(directory, id_name)
directory = '/u/haskeysr/FIDASIM/RESULTS/D3D/158676/01010/MAIN_ION330/'
id_name = 'def'
fida_plot_grid(directory, id_name)

import matplotlib
matplotlib.use('TkAgg')

#from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar


#import gtk
import numpy as np
import time

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

#import datetime, sys
import scipy.io.netcdf as netcdf
#import matplotlib.pyplot as pt
import os
HOME = os.environ['HOME']
shot = 160409
netcdf_files = []
temp_list = []
fit_list = []
intensity_list = []
# 'location' peak location, each column is a peak

for tang_ch in range(1,6):
    netcdf_files.append(netcdf.netcdf_file(HOME + '/cerfit/{shot}/d{shot}_tang{}.nc'.format(tang_ch,shot=shot)))
    temp_list.append(netcdf_files[-1].variables['temperature'].data[:,0])
    fit_list.append(netcdf_files[-1].variables['fit'].data)
    intensity_list.append(netcdf_files[-1].variables['intensity'].data)

# uncomment to select /GTK/GTKAgg/GTKCairo
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas

class gui():
    def __init__(self,):
        self.root = Tk.Tk()
        #self.win = gtk.Window()
        self.root.wm_title("CER analysis")
        self.win.connect("destroy", lambda x: gtk.main_quit())
        #self.win.set_default_size(400,300)
        #self.win.set_title("Embedding in GTK")
        self.f = Figure(figsize=(5,4), dpi=100)
        self.axis_images = []
        self.axis_plot = []
        for i in range(4):
            self.axis_images.append(self.f.add_subplot(4,2,i+1))
        for i in range(4,4+4):
            self.axis_plot.append(self.f.add_subplot(4,2,i+1))
        #self.a = self.f.add_subplot(111)
        self.im_quant = 'fit'
        self.im_quant = 'residuals'
        for i in range(4):
            self.axis_images[i].imshow(netcdf_files[i].variables[self.im_quant].data, aspect = 'auto',cmap='spectral')
        #t = np.arange(0.0,3.0,0.01)
        #s = np.sin(2*pi*t)
        #self.a.plot(t,s)

        self.canvas = FigureCanvas(self.f, master = root)  # a gtk.DrawingArea
        self.toolbar = NavigationToolbar(self.canvas, self.root)
        self.mainframe = Tk.ttk.Frame(root)
        #self.vbox = gtk.VBox(gtk.FALSE, 0)
        #self.win.add(self.vbox)
        
        #self.vbox.pack_start(self.canvas)
        #self.vbox.pack_start(self.toolbar, False, False)
        #self.vbox.pack_start(self.canvas2)
        #vbox.add(canvas)
        #vbox.add(toolbar)
        #vbox.show()
        self.x_val = 0
        #self.win.connect("key-press-event",self.on_window_key_press_event)

        self.startup = True
        self.clicked()
        #self.win.show_all()
        self.startup = False
    def on_window_key_press_event(self,window, event):
        print event.state, event.keyval
        value = gtk.gdk.keyval_name(event.keyval)
        if value == 'f':
            self.x_val += netcdf_files[0].variables['intensity'].data.shape[0]/2
            self.clicked()

        
    def clicked(self,):
        start_time = time.time()
        if self.startup:
            self.int_lines = []
            self.fit_lines = []
        #if self.f.canvas.manager.toolbar._active is None:
        if True:
            for i in range(4):
                if self.startup:
                    #self.axis_plot[i].cla()
                    self.int_lines.append(self.axis_plot[i].plot(netcdf_files[i-4].variables['intensity'].data[self.x_val,:],linestyle='-',color='b')[0])
                    self.fit_lines.append(self.axis_plot[i].plot(netcdf_files[i-4].variables['fit'].data[self.x_val,:],linestyle='--',color='r')[0])
                else:
                    #self.axis_plot[i].cla()
                    self.int_lines[i].set_ydata(netcdf_files[i-4].variables['intensity'].data[self.x_val,:])
                    self.fit_lines[i].set_ydata(netcdf_files[i-4].variables['fit'].data[self.x_val,:])
            mid_time = time.time()
            #if self.startup:
            if True:
                self.f.canvas.draw()
            else:
                for i in range(4):
                    self.axis_plot[i].draw_artist(self.axis_plot[i].patch)
                    self.axis_plot[i].draw_artist(self.int_lines[i])
                    self.axis_plot[i].draw_artist(self.fit_lines[i])
            
            self.f.canvas.flush_events()
            print -start_time+ mid_time, time.time() - mid_time
        else:
            print 'not doing anything....'
def onclick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata)
    gui1.x_val = int(np.round(event.xdata))
    gui1.y_val = int(np.round(event.ydata))
    if event.inaxes in gui1.axis_images:
        gui1.clicked()
    else:
        print 'not in image axis'
gui1 = gui()
cid = gui1.f.canvas.mpl_connect('button_press_event', onclick)
gtk.main()


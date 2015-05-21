from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import gtk
import numpy as np
import socket, time
import datetime, sys
import scipy.io.netcdf as netcdf
import matplotlib.pyplot as pt
import os
HOME = os.environ['HOME']
shot = 160409
netcdf_files = []
temp_list = []
fit_list = []
intensity_list = []
# 'location' peak location, each column is a peak

#for tang_ch in range(1,9):
#    netcdf_files.append(netcdf.netcdf_file(HOME + '/cerfit/{shot}/d{shot}_tang{}.nc'.format(tang_ch,shot=shot)))
#    temp_list.append(netcdf_files[-1].variables['temperature'].data[:,0])
#    fit_list.append(netcdf_files[-1].variables['fit'].data)
#    intensity_list.append(netcdf_files[-1].variables['intensity'].data)

# uncomment to select /GTK/GTKAgg/GTKCairo
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas

class gui():
    def __init__(self,):
        self.win = gtk.Window()
        self.win.connect("destroy", lambda x: gtk.main_quit())
        self.win.connect("key-press-event",self.on_window_key_press_event)
        self.win.set_default_size(400,300)
        self.win.set_title("Embedding in GTK")

        self.cold_line_select = False

        #for i in range(4,4+4):
        #self.a = self.f.add_subplot(111)

        self.im_quant = 'residuals'

        self.notebook = gtk.Notebook()
        self.win.add(self.notebook,)

        self.plot_dict = {'resid':{'plot':1,'axes':[]},
                          'spec':{'plot':0,'axes':[]},
                          'resid_dot':{'plot':2,'axes':[]},
                          'temp':{'plot':3,'axes':[]},
                          'vel':{'plot':-1,'axes':[]},
                          'amp':{'plot':-1,'axes':[]}}

        self.get_netcdf()
        self.chan_sels, self.chans_vbox, self.plot_checkboxes = self.selection_frame()
        self.f, self.canvas, self.vbox, self.toolbar = self.add_matplotlib()
        self.resids_figure, self.resids_canvas, self.resids_parent, self.resids_toolbar = self.add_matplotlib()
        self.figure_parent = self.vbox
        self.status_label = gtk.Label()
        self.action_status = gtk.Label()

        entry_hbox = gtk.HBox()
        self.entry_text = gtk.Entry()
        self.entry_button = gtk.Button(label='Enter')
        self.entry_button.connect("clicked", self.entry_button_pressed, "cool button")
        entry_hbox.pack_start(self.entry_text)
        entry_hbox.pack_start(self.entry_button, expand = False, fill = False, padding = 0)
        self.entry_text.set_editable(False)
        vbox = gtk.VBox()
        self.action_list_label = gtk.Label()
        vbox.pack_start(self.action_list_label)
        self.action_list = []

        hbox = gtk.HBox(homogeneous = False)
        self.figure_parent.pack_end(entry_hbox, expand = False, fill = False, padding = 0)
        self.figure_parent.pack_end(self.status_label, expand = False, fill = False, padding = 0)
        self.figure_parent.pack_end(self.action_status, expand = False,fill = False, padding = 0)
        hbox.pack_start(self.figure_parent, expand = True, fill = True)
        hbox.pack_end(self.chans_vbox, expand = False,fill = False, padding = 0)

        self.notebook.append_page(hbox)
        self.notebook.append_page(self.resids_parent)
        self.notebook.append_page(vbox)
        #self.notebook.append_page(self.chans_vbox)

        self.x_val = 1000
        self.get_active_channels()

        self.plot_all_residuals()

        #self.notebook.show()


        self.win.show_all()

    def get_active_channels(self,*args):
        #self.plot_channels = ['tang{}'.format(i) for i in range(1,9)]
        #self.plot_channels.append('tang13')
        #print 'hello world'
        self.plot_channels = []
        for i in self.chan_sels:
            print i.get_active(), i.get_label()
            if i.get_active(): self.plot_channels.append(i.get_label())
        self.n_chans = len(self.plot_channels)
        print self.plot_channels, self.n_chans
        count = 0
        for i in self.plot_checkboxes:
            if i.get_active(): 
                self.plot_dict[i.get_label()]['plot'] = count
                print 'active plot:',i
                count+=1
            else:
                self.plot_dict[i.get_label()]['plot'] = -1
                print 'inactive plot:',i
        self.t = self.netcdf_dict[self.plot_channels[0]].variables['time'].data
        for i in self.plot_channels:
            self.t = np.append(self.t, self.netcdf_dict[i].variables['time'].data)
        print self.t, self.t.shape
        self.t = np.sort(np.unique(self.t))
        print self.t
        self.generate_figure_layout()
        self.startup_sequence()
        #def update_plots(self,*args):

    def get_active_plots(self,*args):
        self.plot_channels = []
        for i in self.chan_sels:
            print i.get_active(), i.get_label()
            if i.get_active(): self.plot_channels.append(i.get_label())
        self.n_chans = len(self.plot_channels)
        print self.plot_channels, self.n_chans
        self.generate_figure_layout()
        self.startup_sequence()
        #def update_plots(self,*args):


    def selection_frame(self,):
        vbox = gtk.VBox()
        #window.add(vbox)
        checkbox_list = []
        count = 0
        for row in self.avail_chans:
            hbox = gtk.HBox(homogeneous=True)
            #hbox.pack_start(gtk.Label(row))
            enabled_checkbox = gtk.CheckButton(row)
            checkbox_list.append(enabled_checkbox)
            if count <=5:
                enabled_checkbox.set_active(True)#row[2] is '1')
            else:
                enabled_checkbox.set_active(False)#row[2] is '1')
            hbox.pack_start(enabled_checkbox)
            vbox.pack_start(hbox)
            count += 1
        plot_checkboxes_list = []
        self.plot_lists = ['spec','resid','resid_dot','temp','vel','amp']

        vbox.pack_start(gtk.HSeparator())
        for tmp_key in self.plot_lists:
            hbox = gtk.HBox(homogeneous=True)
            enabled_checkbox = gtk.CheckButton(tmp_key)
            plot_checkboxes_list.append(enabled_checkbox)
            if tmp_key in ['spec','resid','temp']:
                enabled_checkbox.set_active(True)
            else:
                enabled_checkbox.set_active(False)
            hbox.pack_start(enabled_checkbox)
            vbox.pack_start(hbox)
        vbox.pack_start(gtk.HSeparator())
        go_button = gtk.Button(label='Update')
        go_button.connect("clicked", self.get_active_channels, "cool button")
        vbox.pack_start(go_button)
        return checkbox_list, vbox, plot_checkboxes_list

    def plot_all_residuals(self,):
        ncols = 4
        nrows = np.ceil(len(self.netcdf_dict.keys())/float(ncols))
        resids_all_axes = [self.resids_figure.add_subplot(nrows, ncols, 1)]
        print 'length', len(self.netcdf_dict.keys())
        for j in range(1,len(self.netcdf_dict.keys())):
            resids_all_axes.append(self.resids_figure.add_subplot(nrows, ncols, j, sharex = resids_all_axes[0], sharey = resids_all_axes[0]))
        
        for i,tmp_key in enumerate(self.netcdf_dict.keys()):
            #rel_axes[i].imshow(netcdf_files[i].variables[self.im_quant].data, aspect = 'auto',cmap='spectral')
            im = resids_all_axes[i].pcolormesh(self.netcdf_dict[tmp_key].variables['time'].data,
                                               np.arange(self.netcdf_dict[tmp_key].variables[self.im_quant].data.shape[1]),
                                               self.netcdf_dict[tmp_key].variables[self.im_quant].data.transpose(), cmap='RdBu',)
            im.set_clim([-4,4])
        self.resids_figure.canvas.draw()

    def add_matplotlib(self,):
        f = Figure(figsize=(5,4), dpi=100)
        canvas = FigureCanvas(f)  # a gtk.DrawingArea
        vbox = gtk.VBox(gtk.FALSE, 0)
        toolbar = NavigationToolbar(canvas, vbox)
        vbox.pack_start(canvas)
        vbox.pack_start(toolbar, False, False)
        return f, canvas, vbox, toolbar
        
    def startup_sequence(self,):
        self.startup = True
        #self.clicked()

        if self.plot_dict['resid']['plot']>=0:self.plot_residuals()
        if self.plot_dict['temp']['plot']>=0:self.plot_temperatures()
        if self.plot_dict['spec']['plot']>=0:self.plot_spectra()
        if self.plot_dict['resid_dot']['plot']>=0:self.plot_residual_dots()


        self.f.tight_layout(pad=0)
        self.f.canvas.draw()
        self.startup = False

    def get_netcdf(self,):
        #for tang_ch in range(1,9):
        self.netcdf_files = []
        self.netcdf_dict = {}
        dir_loc = HOME + '/cerfit/{shot}'.format(shot = shot)
        dir_list = os.listdir(dir_loc)
        filt_list = []
        chan_list = []; vert_list = []; tang_list = []
        for i in dir_list:
            file_start = 'd{shot}_'.format(shot = shot)
            if i.find(file_start)==0:
                print i
                if i.find('tang')>0:
                    tang_list.append(i.replace(file_start, '',).replace('.nc',''))
                elif i.find('vert')>0:
                    vert_list.append(i.replace(file_start, '',).replace('.nc',''))
                else:
                    print 'not vert or tang'
        #for i in range(100):
        #    if tang_list.find(')
        print tang_list, vert_list, sorted(tang_list), sorted(vert_list)
        chan_list = sorted(tang_list)
        for i in sorted(vert_list): chan_list.append(i) 
        #chan_list.append(i.replace(file_start, '',).replace('.nc',''))#.replace('tang','t').replace('vert','v'))
        self.avail_chans = chan_list
        print 'hello',self.avail_chans
        for i in self.avail_chans:
            f = netcdf.netcdf_file(dir_loc +  '/d{shot}_{}.nc'.format(i,shot=shot))
            self.netcdf_dict[i] = f

    def plot_residuals(self,):
        rel_axes = self.plot_dict['resid']['axes']
        self.plot_dict['resid']['im'] = []
        #for i in range(self.n_chans):
        for i,id in enumerate(self.plot_channels):
            
            #rel_axes[i].imshow(netcdf_files[i].variables[self.im_quant].data, aspect = 'auto',cmap='spectral')
            im = rel_axes[i].pcolormesh(self.netcdf_dict[id].variables['time'].data,
                                        np.arange(self.netcdf_dict[id].variables[self.im_quant].data.shape[1]),
                                        self.netcdf_dict[id].variables[self.im_quant].data.transpose(), cmap='RdBu',)
            self.plot_dict['resid']['im'].append(im)
            im.set_clim([-4,4])

    def plot_temperatures(self,):
        rel_axes = self.plot_dict['temp']['axes']
        if self.startup:
            self.plot_dict['temp']['lines'] = []
            self.plot_dict['temp']['vlines'] = []
            for i,id in enumerate(self.plot_channels):
            #for i in range(len(self.plot_channels)):
                self.plot_dict['temp']['vlines'].append(rel_axes[i].axvline(self.x_val))
                self.plot_dict['temp']['lines'].append(rel_axes[i].plot(self.netcdf_dict[id].variables['time'].data, self.netcdf_dict[id].variables['temperature'].data[:,0],'bx')[0])
        else:
            #for i in range(len(self.plot_channels)):
            for i,(id,xval) in enumerate(zip(self.plot_channels,self.time_values)):
            #for i in range(self.n_chans):
                self.plot_dict['temp']['vlines'][i].set_xdata([xval,xval])
                rel_axes[i].draw_artist(rel_axes[i].patch)
                rel_axes[i].draw_artist(self.plot_dict['temp']['vlines'][i])
                rel_axes[i].draw_artist(self.plot_dict['temp']['lines'][i])
                self.f.canvas.blit(rel_axes[i].bbox)

    # def plot_vlines(self,):
    #     rel_axes = self.plot_dict['temp']['axes']
    #     if self.startup:
    #     else:
    #         for i in range(self.n_chans):
    #             self.plot_dict['temp']['vlines'][i].set_xdata([self.x_val,self.x_val])
    #             rel_axes[i].draw_artist(rel_axes[i].patch)
    #             rel_axes[i].draw_artist(self.plot_dict['temp']['vlines'][i])
    #             self.f.canvas.blit(rel_axes[i].bbox)


    def generate_figure_layout(self,):
        self.f.clf()
        self.ncols = self.n_chans
        self.nrows = np.sum(np.array([self.plot_dict[i]['plot'] for i in self.plot_dict.keys()])>=0)

        for i in self.plot_dict.keys():
            if self.plot_dict[i]['plot']>=0:

                if i=='resid':
                    self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, 0+1 + self.plot_dict[i]['plot'] *self.n_chans)]
                    for j in range(1,self.n_chans):
                        self.plot_dict[i]['axes'].append(self.f.add_subplot(self.nrows, self.ncols, j+1 + self.plot_dict[i]['plot'] *self.n_chans, sharex = self.plot_dict[i]['axes'][0], sharey = self.plot_dict[i]['axes'][0]))
                elif i=='temp':
                    if self.plot_dict['resid']['plot']>=0: 
                        ref_ax = self.plot_dict['resid']['axes'][0]
                        self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, 0+1 + self.plot_dict[i]['plot'] *self.n_chans,sharex = ref_ax)]
                    else:
                        self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, 0+1 + self.plot_dict[i]['plot'] *self.n_chans)]
                        ref_ax = self.plot_dict[i]['axes'][0]
                    for j in range(1,self.n_chans):
                        self.plot_dict[i]['axes'].append(self.f.add_subplot(self.nrows, self.ncols, j+1 + self.plot_dict[i]['plot'] *self.n_chans, sharex = ref_ax))
                else:
                    self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, j+1 + self.plot_dict[i]['plot'] *self.n_chans) for j in range(self.n_chans)]

            #self.axis_plot.append(self.f.add_subplot(self.nrows,self.ncols,i+1+(self.n_chans)))
        

    def on_window_key_press_event(self,window, event):
        print event.state, event.keyval
        print self.x_val
        value = gtk.gdk.keyval_name(event.keyval)
        best = 0
        for i in self.plot_channels:
            if len(self.netcdf_dict[i].variables['time'].data)>best:
                tmp_key = i
                best = len(self.netcdf_dict[i].variables['time'].data)
        print 'using ', tmp_key
        #t = self.netcdf_dict[tmp_key].variables['time'].data
        if value == 'f': 
            loc = np.argmin(np.abs(self.t - self.x_val))
            print loc, self.x_val
            self.x_val = np.min([self.t[loc+1],np.max(self.t)])
            print loc, self.x_val
            self.clicked()
        elif value == 'b':
            loc = np.argmin(np.abs(self.t - self.x_val))
            #loc = np.argmin(np.abs(self.netcdf_dict[tmp_key].variables['time'].data - self.x_val))
            self.x_val = np.min([self.t[loc-1],np.max(self.t)])
            self.clicked()
        elif value=='r':
            self.startup_sequence()
            self.clicked()
        if value == 'c': 
            self.cold_line_select = True
            self.action_status.set_label('RIGHT CLICK TO SELECT A COLD LINE LOCATION')
        elif value =='k':
            print 'kill a certain datapoint'
            possible_list = []
            for i in self.netcdf_dict.keys():
                if self.x_val in self.netcdf_dict[i].variables['time']:
                    possible_list.append(i)
            print possible_list, len(possible_list)
            self.action_list.append('kill {} {}'.format(','.join(possible_list), self.x_val))
            self.action_list_label.set_label('\n'.join(self.action_list))

        elif value =='m':
            print 'modify a certain datapoint'
            possible_list = []
            for i in self.netcdf_dict.keys():
                if self.x_val in self.netcdf_dict[i].variables['time']:
                    possible_list.append(i)
            print possible_list, len(possible_list)
            self.entry_text.set_text('modify {} {}'.format(','.join(possible_list), self.x_val))
            self.entry_text.set_editable(True)
        #self.x_val += self.netcdf_dict[tmp_key].variables['intensity'].data.shape[tmp_key]/2

    def entry_button_pressed(self,*args):
        self.action_list.append(self.entry_text.get_text())
        self.action_list_label.set_label('\n'.join(self.action_list))
        self.entry_text.set_text('')
        self.entry_text.set_editable(False)


    def create_status_string(self,):
        tmp_str = ','.join(['{}:{:.2f}'.format(ch,tmp_time) for tmp_time, ch in zip(self.time_values,self.plot_channels)])
        self.status_label.set_label(tmp_str)
        print 'tmp_str', tmp_str

    def plot_spectra(self,):
        rel_axes = self.plot_dict['spec']['axes']
        if self.startup:
            self.plot_dict['spec']['line1'] = []
            self.plot_dict['spec']['line2'] = []
        #for i in range(self.n_chans):
        self.time_values = []
        for i,id in enumerate(self.plot_channels):
            #rel_axes[i].cla()
            loc = np.argmin(np.abs(self.netcdf_dict[id].variables['time'].data - self.x_val))
            print self.x_val, loc, self.netcdf_dict[id].variables['intensity'].data.shape
            self.time_values.append(self.netcdf_dict[id].variables['time'].data[loc])
            if self.startup:
                fit_dat = self.netcdf_dict[id].variables['fit'].data[loc,:]
                max_val = np.max(fit_dat)*1.2
                min_val = -np.max(fit_dat)*0.1
                self.plot_dict['spec']['line1'].append(rel_axes[i].plot(self.netcdf_dict[id].variables['intensity'].data[loc,:],linestyle='-',color='b')[0])
                self.plot_dict['spec']['line2'].append(rel_axes[i].plot(fit_dat,linestyle='--',color='k')[0])
                rel_axes[i].set_title(id)
                rel_axes[i].set_ylim([min_val,max_val])
            else:
                self.plot_dict['spec']['line1'][i].set_ydata(self.netcdf_dict[id].variables['intensity'].data[loc,:])
                self.plot_dict['spec']['line2'][i].set_ydata(self.netcdf_dict[id].variables['fit'].data[loc,:])
                rel_axes[i].draw_artist(rel_axes[i].patch)
                for tmp_spine in rel_axes[i].spines.values(): rel_axes[i].draw_artist(tmp_spine)
                rel_axes[i].draw_artist(self.plot_dict['spec']['line1'][i])
                rel_axes[i].draw_artist(self.plot_dict['spec']['line2'][i])
                self.f.canvas.blit(rel_axes[i].bbox)
        self.create_status_string()

    def plot_residual_dots(self,):
        rel_axes = self.plot_dict['resid_dot']['axes']
        if self.startup:
            self.plot_dict['resid_dot']['dots'] = []
        #for i in range(self.n_chans):
        for i,id in enumerate(self.plot_channels):
            #rel_axes[i].cla()
            loc = np.argmin(np.abs(self.netcdf_dict[id].variables['time'].data - self.x_val))
            print self.x_val, loc, self.netcdf_dict[id].variables['intensity'].data.shape
            data = self.netcdf_dict[id].variables['intensity'].data[loc,:] - self.netcdf_dict[id].variables['fit'].data[loc,:]
            if self.startup:
                self.plot_dict['resid_dot']['dots'].append(rel_axes[i].plot(data,'o')[0])
                #self.plot_dict['spec']['line2'].append(,linestyle='--',color='k')[0])
                #rel_axes[i].set_title(id)
                rel_axes[i].set_ylim([-10,10])
            else:
                self.plot_dict['resid_dot']['dots'][i].set_ydata(data)
                rel_axes[i].draw_artist(rel_axes[i].patch)
                for tmp_spine in rel_axes[i].spines.values(): rel_axes[i].draw_artist(tmp_spine)
                rel_axes[i].draw_artist(self.plot_dict['resid_dot']['dots'][i])
                self.f.canvas.blit(rel_axes[i].bbox)


    def clicked(self,):
        start_time = time.time()
        #if.self.plot_dict['resid']['plot']>0:self.plot_residuals()
        if self.plot_dict['spec']['plot']>=0:self.plot_spectra()
        if self.plot_dict['temp']['plot']>=0:self.plot_temperatures()
        if self.plot_dict['resid_dot']['plot']>=0:self.plot_residual_dots()
        #self.plot_spectra()
        #self.plot_temperatures()
        self.f.canvas.flush_events()

def onclick(event):
    #print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #    event.button, event.x, event.y, event.xdata, event.ydata)
    print 'button',event.button
    if gui1.cold_line_select:
        if (event.inaxes in gui1.plot_dict['resid']['axes']) and event.button==3:
            line_at = int(np.round(event.ydata))
            ch_name = gui1.plot_channels[gui1.plot_dict['resid']['axes'].index(event.inaxes)]
        elif (event.inaxes in gui1.plot_dict['resid_dot']['axes']) and event.button==3:
            line_at = int(np.round(event.xdata))
            ch_name = gui1.plot_channels[gui1.plot_dict['resid']['axes'].index(event.inaxes)]
        elif (event.inaxes in gui1.plot_dict['spec']['axes']) and event.button==3:
            line_at = int(np.round(event.xdata))
            ch_name = gui1.plot_channels[gui1.plot_dict['spec']['axes'].index(event.inaxes)]
        print '####### line : ', line_at, ch_name
        gui1.action_list.append('cold_line {} {}'.format(ch_name, line_at))
        gui1.cold_line_select = False
        gui1.action_list_label.set_label('\n'.join(gui1.action_list))
        gui1.action_status.set_label('')
    if (event.inaxes in gui1.plot_dict['resid']['axes']) and event.button==3:
        gui1.x_val = gui1.t[np.argmin(np.abs(gui1.t - event.xdata))]
        print '****', gui1.x_val, event.xdata
        #int(np.round(event.xdata))
        gui1.y_val = int(np.round(event.ydata))
        gui1.clicked()
    elif (event.inaxes in gui1.plot_dict['temp']['axes']) and event.button==3:
        gui1.x_val = gui1.t[np.argmin(np.abs(gui1.t - event.xdata))]
        print '****', gui1.x_val, event.xdata
        #gui1.x_val = int(np.round(event.xdata))
        gui1.y_val = int(np.round(event.ydata))
        #gui1.startup = True
        # for i in gui1.plot_dict.keys():
        #     try:
        #         print gui1.plot_dict[i]['axes']
        #         val = True
        #     except:
        #         val = False

        #     if val:
        #         for j in gui1.plot_dict[i]['axes']:j.cla()
        # #gui1.f.clf()
        gui1.clicked()
        #gui1.f.canvas.draw()
        #gui1.startup = False

gui1 = gui()
cid = gui1.f.canvas.mpl_connect('button_press_event', onclick)
gtk.main()

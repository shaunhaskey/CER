'''
Todo list:
-> Run CERFIT from within this widget for minor changes like modifying the time subtraction slice
-> Plot the profiles as a function of time, and allow stepping through this data
-> Modification, removal etc... of action items from the widget
'''
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import datetime, sys,socket, time, gtk, os
import scipy.io.netcdf as netcdf
import scipy.interpolate as interp
import numpy as np
import matplotlib.pyplot as pt
import sys
import cer_funcs as CER
HOME = os.environ['HOME']
shot = int(sys.argv[1])

#shot = 160409

def change_name(ch_name):
    if ch_name.find('tang')>=0:
        ch_name = 't{:02d}'.format(int(ch_name[4:]))
    elif ch_name.find('vert')>=0:
        ch_name = 'v{:02d}'.format(int(ch_name[4:]))
    return ch_name

class gui():
    def __init__(self,shot):
        self.shot = shot
        self.win = gtk.Window()
        self.win.connect("destroy", lambda x: gtk.main_quit())
        self.win.connect("key-press-event",self.on_window_key_press_event)
        self.win.set_default_size(400,300)
        self.win.set_title("CERFIT helper")
        self.cold_line_select = False
        self.im_quant = 'residuals'
        self.notebook = gtk.Notebook()
        self.overall_vbox = gtk.VBox()
        self.overall_vbox.pack_start(self.notebook, expand = True, fill = True, padding = 0)
        #self.win.add(self.notebook,)
        self.win.add(self.overall_vbox,)
        self.plot_dict = {'resid':{'plot':1,'axes':[],'plot_func':self.plot_residuals, 'refresh':False, 'args':[],'xax':''},
                          'spec':{'plot':0,'axes':[],'plot_func':self.plot_spectra, 'refresh':True, 'args':[],'xax':'pixel'},
                          'resid_dot':{'plot':2,'axes':[],'plot_func':self.plot_residual_dots, 'refresh':True, 'args':[],'xax':'pixel'},
                          'temp':{'plot':3,'axes':[],'plot_func':self.plot_time_varying_1D, 'refresh':True, 'args':['temperature', 'temp'],'xax':'time'},
                          'vel':{'plot':-1,'axes':[],'plot_func':self.plot_time_varying_1D, 'refresh':True, 'args':['location','vel'],'xax':'time'},
                          'amp':{'plot':-1,'axes':[],'plot_func':self.plot_time_varying_1D, 'refresh':True, 'args':['amplitude','amp'],'xax':'time'}}

        self.fig_radial, self.canvas_radial, self.vbox_radial, self.toolbar_radial = self.add_matplotlib()
        self.fig_timing, self.canvas_timing, self.vbox_timing, self.toolbar_timing = self.add_matplotlib()

        self.ax_radial = [self.fig_radial.add_subplot(2,2,i) for i in range(1,5)]
        self.ax_timing = [self.fig_timing.add_subplot(1,1,i) for i in range(1,2)]
        #Get the available netCDF data 
        self.get_netcdf()
        #Setup the plot selection frame
        self.chan_sels, self.chans_vbox, self.plot_checkboxes = self.selection_frame()
        #Add the main matplotlib 'canvas'
        self.f, self.canvas, self.vbox, self.toolbar = self.add_matplotlib()
        #Another tab for the notebook that all of the residuals plotted - useful?
        self.resids_figure, self.resids_canvas, self.resids_parent, self.resids_toolbar = self.add_matplotlib()
        #setup an alias
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
        self.action_list_scroll = gtk.ScrolledWindow()
        self.action_list_scroll.set_border_width(10)
        self.action_list_scroll.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)

        #self.get_timing_info()
        self.action_list_button = gtk.Button(label='Write File')
        self.timing_button = gtk.Button(label='Timing')
        self.action_list_button.connect("clicked", self.action_list_button_pressed, "cool button")
        self.timing_button.connect("clicked", self.update_timing, "cool button")

        self.action_list_scroll.add_with_viewport(self.action_list_label)
        #vbox.pack_start(self.action_list_label)
        vbox_timing = gtk.VBox()
        self.timing_combobox_chan = gtk.combo_box_new_text()
        self.timing_combobox_beam = gtk.combo_box_new_text()
        for i in ['30','33','21']:
            self.timing_combobox_beam.append_text(i+'r')
            self.timing_combobox_beam.append_text(i+'l')
        hbox_timing = gtk.HBox()
        hbox_timing.pack_start(self.timing_combobox_chan)
        hbox_timing.pack_start(self.timing_combobox_beam)
        hbox_timing.pack_start(self.timing_button)
        
        for i in self.chan_sels: self.timing_combobox_chan.append_text(i.get_label())

        vbox_timing.pack_start(hbox_timing, expand = False, fill = False, padding = 0)
        #vbox_timing.pack_start(self.timing_button, expand = False, fill = False, padding = 0)
        vbox_timing.pack_start(self.vbox_timing, expand = True, fill = True, padding = 0)
        vbox.pack_start(self.action_list_scroll)
        vbox.pack_start(self.action_list_button, expand = False, fill = False, padding = 0)
        #vbox.pack_start(self.timing_button)


        tmp_txt = self.action_list_label.get_label()
        fname = HOME + '/cerfit/{shot}/extra_cmds2.txt'.format(shot = self.shot)
        if os.path.isfile(fname):
            with file(fname,'r') as filehandle: 
                self.action_list = [i.rstrip('\n') for i in filehandle.readlines()]
        else:
            self.action_list = []
        self.action_list_label.set_label('\n'.join(self.action_list))

        hbox = gtk.HBox(homogeneous = False)
        #self..pack_start(self.notebook, expand = True, fill = True, padding = 0)
        self.overall_vbox.pack_end(entry_hbox, expand = False, fill = False, padding = 0)
        self.figure_parent.pack_end(self.status_label, expand = False, fill = False, padding = 0)
        self.figure_parent.pack_end(self.action_status, expand = False,fill = False, padding = 0)
        hbox.pack_start(self.figure_parent, expand = True, fill = True)
        hbox.pack_end(self.chans_vbox, expand = False,fill = False, padding = 0)

        self.notebook.append_page(hbox,gtk.Label("Data"))
        self.notebook.append_page(self.resids_parent,gtk.Label("Residuals"))
        self.notebook.append_page(vbox, gtk.Label("ActionItems"))
        self.notebook.append_page(self.vbox_radial,gtk.Label('RadialProfiles'))
        self.notebook.append_page(vbox_timing,gtk.Label('Timing'))

        #Set an initial time value
        self.x_val = 1000
        self.get_active_channels()

        self.plot_all_residuals()
        self.win.show_all()

    def get_timing_info(self,beam_chan, chord):
        shot_data = CER.read_shot_status(shot)
        for i in shot_data.keys():
            print i, shot_data[i]['CHOR'].lower()
            if shot_data[i]['CHOR'].lower().find(chord)>=0:
                spec = i
        print spec, shot_data[spec]['CHOR']
        stri = shot_data[spec]['TIMI']
        print stri
        #stri = '0.0:2000@5.0/1'
        self.start_time, self.n_slices, self.t_int, self.start_pts, self.end_pts = CER.split_timing_string(stri)
        self.nbi_x, self.nbi_y = CER.get_nbi_data(self.shot, beam_chan)
        self.fscope_x, self.fscope_y = CER.get_fscope(self.shot,)
        if beam_chan.find('33')>=0: 
            beam_chan_mod = beam_chan.replace('33','330')
        elif beam_chan.find('21')>=0: 
            beam_chan_mod = beam_chan.replace('21','210')
        else:
            beam_chan_mod = beam_chan
        ts_fname =  HOME + '/cerfit/{shot}/'.format(shot = self.shot) + 'timesub_{}_{}t_{}.dat'.format(chord[0],beam_chan_mod,chord)
        print ts_fname
        self.ts_data = CER.read_subtract_file(ts_fname)
        self.ts_data = np.array([[float(i[0]),float(i[1])] for i in self.ts_data])
        print self.ts_data.shape

    def update_timing(self,*args):
        chord = self.timing_combobox_chan.get_active_text()
        beam = self.timing_combobox_beam.get_active_text()
        chord_num = int(chord[4:])
        chord_tmp = chord[:4].replace('tang','t').replace('vert','v')
        chord = '{}{:02d}'.format(chord_tmp,chord_num)
        print 'beam:{},chord{}'.format(beam, chord)
        #chord = 't01'
        self.get_timing_info(beam, chord)
        self.ax_timing[0].cla()
        self.ax_timing[0].plot(self.nbi_x,self.nbi_y/np.max(self.nbi_y))
        self.ax_timing[0].plot(self.start_pts, self.start_pts*0+1,'o')
        self.ax_timing[0].plot(self.fscope_x,self.fscope_y/np.max(self.fscope_y))
        valid = (self.ts_data[:,0]<(self.x_val+100)) * (self.ts_data[:,0]>(self.x_val-100))
        tmp_x = np.zeros((len(self.ts_data[valid,0]),2),dtype=float)
        tmp_y = np.zeros((len(self.ts_data[valid,0]),2),dtype=float)
        print tmp_x.shape
        tmp_x[:,0] = self.ts_data[valid,0]
        tmp_x[:,1] = self.ts_data[valid,0] + self.ts_data[valid,1]*self.t_int
        tmp_y[:,0] = self.ts_data[valid,0]*0+0.8
        tmp_y[:,1] = self.ts_data[valid,0]*0+0.6
        
        for i in range(tmp_x.shape[0]):
            print tmp_x[i,:], tmp_y[i,:]
            self.ax_timing[0].plot([tmp_x[i,0],tmp_x[i,0]+self.t_int], [0.8,0.8],'xb-')
            self.ax_timing[0].plot([tmp_x[i,1],tmp_x[i,1]+self.t_int], [0.6,0.6],'xb-')
            self.ax_timing[0].plot(tmp_x[i,:]+self.t_int/2., tmp_y[i,:], '-b')

        #self.ax_timing[0].plot(tmp_x, tmp_y, 'k-x')
        #self.ax_timing[0].plot(self.ts_data[valid,0], self.ts_data[valid,0]*0+0.8,'x')
        #self.ax_timing[0].plot(self.ts_data[valid,0]+self.ts_data[valid,1]*self.t_int, self.ts_data[valid,0]*0+0.6,'d')
        
        self.ax_timing[0].set_xlim([self.x_val-50, self.x_val+50])
        self.ax_timing[0].set_ylim([0,1.3])
        self.ax_timing[0].axvline(self.x_val)
        self.fig_timing.canvas.draw()

    def action_list_button_pressed(self,*args):
        tmp_txt = self.action_list_label.get_label()
        fname = HOME + '/cerfit/{shot}/extra_cmds2.txt'.format(shot = self.shot)
        with file(fname,'w') as filehandle: filehandle.write(tmp_txt+'\n')

    def get_active_channels(self,*args):
        '''Figures out which channels should be plotted, and which plots should be shown

        SRH: 21May2015
        '''
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
        '''Figures out which channels should be plotted, and which plots should be shown
        SRH: 21May2015
        '''
        self.plot_channels = []
        for i in self.chan_sels:
            print i.get_active(), i.get_label()
            if i.get_active(): self.plot_channels.append(i.get_label())
        self.n_chans = len(self.plot_channels)
        print self.plot_channels, self.n_chans
        self.generate_figure_layout()
        self.startup_sequence()

    def selection_frame(self,):
        vbox = gtk.VBox()
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
        self.calc_closest_times()
        for i in self.plot_dict.keys():
            if self.plot_dict[i]['plot']>=0:
                self.plot_dict[i]['plot_func'](*self.plot_dict[i]['args'])
        self.f.tight_layout(pad=0)
        self.f.canvas.draw()
        self.startup = False

    def get_netcdf(self,):
        self.netcdf_files = []
        self.netcdf_dict = {}
        dir_loc = HOME + '/cerfit/{shot}'.format(shot = self.shot)
        dir_list = os.listdir(dir_loc)
        filt_list = []
        chan_list = []; vert_list = []; tang_list = []
        for i in dir_list:
            file_start = 'd{shot}_'.format(shot = self.shot)
            if i.find(file_start)==0:
                print i
                if i.find('tang')>0:
                    tang_list.append(i.replace(file_start, '',).replace('.nc',''))
                elif i.find('vert')>0:
                    vert_list.append(i.replace(file_start, '',).replace('.nc',''))
                else:
                    print 'not vert or tang'
        chan_list_new = []
        for i in range(100):
            cur_name = 'tang{}'.format(i)
            if cur_name in tang_list:
                chan_list_new.append(cur_name)
        for i in range(100):
            cur_name = 'vert{}'.format(i)
            if cur_name in vert_list:
                chan_list_new.append(cur_name)
        print tang_list, vert_list, sorted(tang_list), sorted(vert_list)
        chan_list = sorted(tang_list)
        for i in sorted(vert_list): chan_list.append(i) 
        self.avail_chans = chan_list
        print '##########', len(chan_list), len(chan_list_new)
        self.avail_chans = chan_list_new
        print 'hello',self.avail_chans
        for i in self.avail_chans:
            f = netcdf.netcdf_file(dir_loc +  '/d{shot}_{}.nc'.format(i,shot=self.shot))
            self.netcdf_dict[i] = f

        #Make an overall interpolated time base
        # self.read_in_cer_files()
        # self.overall_t = self.netcdf_dict[self.avail_chans[0]].variables['time'].data
        # for i in self.avail_chans:
        #     self.overall_t = np.append(self.overall_t, self.netcdf_dict[i].variables['time'].data)
        # print self.overall_t, self.overall_t.shape
        # self.overall_t = np.sort(np.unique(self.overall_t))
        # print '##############'
        # self.overall_temp_data = np.zeros((len(self.avail_chans),self.overall_t.shape[0]),dtype=float)
        # for i, chan in enumerate(self.avail_chans):
        #     tmp = interp.interp1d(self.netcdf_dict[chan].variables['time'].data, self.netcdf_dict[chan].variables['temperature'].data[:,0], bounds_error = False, fill_value = 0.,kind='nearest')
        #     self.overall_temp_data[i,:] = tmp(self.overall_t)
        # print self.overall_t, self.overall_temp_data
        # #fig, ax = pt.subplots()
        # arg_ord = np.argsort(self.chan_R)
        # self.t_mesh, self.R_mesh = np.meshgrid(self.overall_t, self.chan_R[arg_ord])
        # im = self.ax_radial[0].pcolormesh(self.t_mesh, self.R_mesh, self.overall_temp_data[arg_ord,:],cmap = 'spectral')
        # im.set_clim([0,5])
        # self.fig_radial.canvas.draw()

    def read_in_cer_files(self,):
        print 'reading cer files'
        dir_loc = HOME + '/cerfit/{shot}'.format(shot = self.shot)
        self.chan_R = []
        for i in self.avail_chans:
            num = int(i[4:])
            if i.find('tang')>=0:
                letter='t'
            else:
                letter = 'v'
            if num<=7 and letter == 't':
                append_letter = 'a'
            if num>=8 and letter == 't':
                append_letter = 'd'
            if letter == 'v':
                append_letter = ''
            fname = '{dir_loc}/{shot}.8{sys}{num}{ext}'.format(dir_loc = dir_loc, shot = self.shot, sys=letter,num=num,ext=append_letter)
            print fname
            with file(fname,'r') as filehandle: lines = filehandle.readlines()
            find_txt = 'RADIUS ='
            start = lines[1].find(find_txt)
            end = lines[1].find('cm')
            self.chan_R.append(float(lines[1][start+len(find_txt):end]))
        self.chan_R = np.array(self.chan_R)
        print self.chan_R

    def plot_residuals(self,):
        rel_axes = self.plot_dict['resid']['axes']
        self.plot_dict['resid']['im'] = []
        for i,id in enumerate(self.plot_channels):
            im = rel_axes[i].pcolormesh(self.netcdf_dict[id].variables['time'].data,
                                        np.arange(self.netcdf_dict[id].variables[self.im_quant].data.shape[1]),
                                        self.netcdf_dict[id].variables[self.im_quant].data.transpose(), cmap='RdBu',)
            self.plot_dict['resid']['im'].append(im)
            im.set_clim([-4,4])


    def plot_time_varying_1D(self,var_name, plot_key):
        #var_name = 'amplitude'
        #plot_key = 'amp'
        #print 'in plot_time_varying_1D ', var_name, plot_key
        rel_axes = self.plot_dict[plot_key]['axes']
        if self.startup:
            self.plot_dict[plot_key]['lines'] = []
            self.plot_dict[plot_key]['vlines'] = []
            for i,id in enumerate(self.plot_channels):
                self.plot_dict[plot_key]['vlines'].append(rel_axes[i].axvline(self.x_val))
                self.plot_dict[plot_key]['lines'].append(rel_axes[i].plot(self.netcdf_dict[id].variables['time'].data, self.netcdf_dict[id].variables[var_name].data[:,0],'bx')[0])
        else:
            for i,(id,xval) in enumerate(zip(self.plot_channels,self.time_values)):
                self.plot_dict[plot_key]['vlines'][i].set_xdata([xval,xval])
                rel_axes[i].draw_artist(rel_axes[i].patch)
                rel_axes[i].draw_artist(self.plot_dict[plot_key]['vlines'][i])
                rel_axes[i].draw_artist(self.plot_dict[plot_key]['lines'][i])
                self.f.canvas.blit(rel_axes[i].bbox)

    def generate_figure_layout(self,):
        self.f.clf()
        self.ncols = self.n_chans
        self.nrows = np.sum(np.array([self.plot_dict[i]['plot'] for i in self.plot_dict.keys()])>=0)
        self.clickable_axes = []
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
                elif i=='vel':
                    if self.plot_dict['resid']['plot']>=0: 
                        ref_ax = self.plot_dict['resid']['axes'][0]
                        self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, 0+1 + self.plot_dict[i]['plot'] *self.n_chans,sharex = ref_ax)]
                    else:
                        self.plot_dict[i]['axes'] = [self.f.add_subplot(self.nrows, self.ncols, 0+1 + self.plot_dict[i]['plot'] *self.n_chans)]
                        ref_ax = self.plot_dict[i]['axes'][0]
                    for j in range(1,self.n_chans):
                        self.plot_dict[i]['axes'].append(self.f.add_subplot(self.nrows, self.ncols, j+1 + self.plot_dict[i]['plot'] *self.n_chans, sharex = ref_ax))
                elif i=='amp':
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
            if self.plot_dict[i]['xax'] =='time':
                self.clickable_axes.extend(self.plot_dict[i]['axes'])

        print len(self.clickable_axes), self.clickable_axes

    def on_window_key_press_event(self,window, event):
        value = gtk.gdk.keyval_name(event.keyval)
        best = 0
        for i in self.plot_channels:
            if len(self.netcdf_dict[i].variables['time'].data)>best:
                tmp_key = i
                best = len(self.netcdf_dict[i].variables['time'].data)
        #t = self.netcdf_dict[tmp_key].variables['time'].data
        if value == 'f': 
            loc = np.argmin(np.abs(self.t - self.x_val))
            #print loc, self.x_val
            self.x_val = np.min([self.t[loc+1],np.max(self.t)])
            #print loc, self.x_val
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
            possible_list = [change_name(ch_name) for ch_name in possible_list]
            print possible_list, len(possible_list)
            self.entry_text.set_text('kill:{}:{:.3f}'.format(','.join(possible_list), self.x_val))
            self.entry_text.set_editable(True)
            #self.action_list_label.set_label('\n'.join(self.action_list))
        elif value =='m':
            print 'modify a certain datapoint'
            possible_list = []
            for i in self.netcdf_dict.keys():
                if self.x_val in self.netcdf_dict[i].variables['time']:
                    possible_list.append(i)
            possible_list = [change_name(ch_name) for ch_name in possible_list]
            print possible_list, len(possible_list)
            self.entry_text.set_text('modify:{}:{},{}'.format(','.join(possible_list), self.x_val,'VAL'))
            self.entry_text.set_editable(True)
        elif value =='d':
            self.entry_text.set_text('remove:')
            self.entry_text.set_editable(True)
        elif value=='h':
            print 'r: reload'
            print 'c: add cold line, cold_line:chords:location'
            print 'k: kill datapoint, kill:chords:time'
            print 'm: modify datapoint modify:chords:time,new_tss'
            print 'd: Delete a chord completely remove:chords'
            print 'lower_time:time'
            print 'upper_time:time'
            print 'd: Delete a chord completely remove:chords'
        #self.x_val += self.netcdf_dict[tmp_key].variables['intensity'].data.shape[tmp_key]/2

    def entry_button_pressed(self,*args):
        txt = self.entry_text.get_text()
        if len(txt)>=1:
            self.action_list.append(self.entry_text.get_text())
            self.action_list_label.set_label('\n'.join(self.action_list))
        self.entry_text.set_text('')
        self.entry_text.set_editable(False)


    def create_status_string(self,):
        tmp_str = ','.join(['{}:{:.2f}'.format(ch,tmp_time) for tmp_time, ch in zip(self.time_values,self.plot_channels)])
        self.status_label.set_label(tmp_str)
        #print 'tmp_str', tmp_str

    def calc_closest_times(self,):
        self.time_values = []; self.locs = []
        for i,id in enumerate(self.plot_channels):
            loc = np.argmin(np.abs(self.netcdf_dict[id].variables['time'].data - self.x_val))
            self.locs.append(np.argmin(np.abs(self.netcdf_dict[id].variables['time'].data - self.x_val)))
            #print self.x_val, loc, self.netcdf_dict[id].variables['intensity'].data.shape
            self.time_values.append(self.netcdf_dict[id].variables['time'].data[loc])

    def plot_spectra(self,):
        rel_axes = self.plot_dict['spec']['axes']
        if self.startup:
            self.plot_dict['spec']['line1'] = []
            self.plot_dict['spec']['line2'] = []
            self.plot_dict['spec']['vlines'] = []
        for i,id in enumerate(self.plot_channels):
            loc = self.locs[i]
            #tmp_time = self.netcdf_dict[id].variables['time'][loc]
            xdata= self.netcdf_dict[id].variables['location'][loc,0] - 1
            if self.startup:
                fit_dat = self.netcdf_dict[id].variables['fit'].data[loc,:]
                max_val = np.max(fit_dat)*1.2
                min_val = -np.max(fit_dat)*0.1
                #self.plot_dict['spec']['line1'].append(rel_axes[i].plot(self.netcdf_dict[id].variables['intensity'].data[loc,:],linestyle='-',color='b')[0])
                self.plot_dict['spec']['line1'].append(rel_axes[i].plot(self.netcdf_dict[id].variables['intensity'].data[loc,:],'b.')[0])
                self.plot_dict['spec']['line2'].append(rel_axes[i].plot(fit_dat,linestyle='--',color='k')[0])
                self.plot_dict['spec']['vlines'].append(rel_axes[i].axvline(xdata))
                rel_axes[i].set_title(id)
                rel_axes[i].set_ylim([min_val,max_val])
            else:
                self.plot_dict['spec']['line1'][i].set_ydata(self.netcdf_dict[id].variables['intensity'].data[loc,:])
                self.plot_dict['spec']['line2'][i].set_ydata(self.netcdf_dict[id].variables['fit'].data[loc,:])
                self.plot_dict['spec']['vlines'][i].set_xdata([xdata,xdata])
                rel_axes[i].draw_artist(rel_axes[i].patch)
                rel_axes[i].draw_artist(self.plot_dict['spec']['vlines'][i])
                for tmp_spine in rel_axes[i].spines.values(): rel_axes[i].draw_artist(tmp_spine)
                rel_axes[i].draw_artist(self.plot_dict['spec']['line1'][i])
                rel_axes[i].draw_artist(self.plot_dict['spec']['line2'][i])
                self.f.canvas.blit(rel_axes[i].bbox)
        self.create_status_string()

    def plot_residual_dots(self,):
        rel_axes = self.plot_dict['resid_dot']['axes']
        if self.startup:
            self.plot_dict['resid_dot']['dots'] = []
        for i,id in enumerate(self.plot_channels):
            loc = np.argmin(np.abs(self.netcdf_dict[id].variables['time'].data - self.x_val))
            print self.x_val, loc, self.netcdf_dict[id].variables['intensity'].data.shape
            data = self.netcdf_dict[id].variables['intensity'].data[loc,:] - self.netcdf_dict[id].variables['fit'].data[loc,:]
            if self.startup:
                self.plot_dict['resid_dot']['dots'].append(rel_axes[i].plot(data,'.')[0])
                rel_axes[i].set_ylim([-10,10])
            else:
                self.plot_dict['resid_dot']['dots'][i].set_ydata(data)
                rel_axes[i].draw_artist(rel_axes[i].patch)
                for tmp_spine in rel_axes[i].spines.values(): rel_axes[i].draw_artist(tmp_spine)
                rel_axes[i].draw_artist(self.plot_dict['resid_dot']['dots'][i])
                self.f.canvas.blit(rel_axes[i].bbox)

    def clicked(self,):
        self.calc_closest_times()
        cur_xlim = self.plot_dict['temp']['axes'][0].get_xlim()
        for i in self.plot_dict.keys():
            if self.plot_dict[i]['plot']>=0 and self.plot_dict[i]['refresh']:
                self.plot_dict[i]['plot_func'](*self.plot_dict[i]['args'])
                #self.plot_dict[i]['plot_func']()
        self.f.canvas.flush_events()
        if (self.x_val>np.max(cur_xlim)) or (self.x_val<np.min(cur_xlim)):
            xlim_range = cur_xlim[1] - cur_xlim[0]
            print 'changing x lim'
            if self.x_val>np.max(cur_xlim):
                self.plot_dict['temp']['axes'][0].set_xlim([cur_xlim[0]+xlim_range*0.9,cur_xlim[1]+xlim_range*0.9])
            if self.x_val<np.min(cur_xlim):
                self.plot_dict['temp']['axes'][0].set_xlim([cur_xlim[0]-xlim_range*0.9,cur_xlim[1]-xlim_range*0.9])
            self.f.canvas.draw()


def onclick(event):
    print 'button',event.button
    if gui1.cold_line_select:
        if (event.inaxes in gui1.plot_dict['resid']['axes']) and event.button==3:
            line_at = int(np.round(event.ydata))
            ch_name = gui1.plot_channels[gui1.plot_dict['resid']['axes'].index(event.inaxes)]
        elif (event.inaxes in gui1.plot_dict['resid_dot']['axes']) and event.button==3:
            line_at = int(np.round(event.xdata))
            ch_name = gui1.plot_channels[gui1.plot_dict['resid']['axes'].index(event.inaxes)]
        elif (event.inaxes in gui1.clickable_axesplot_dict['spec']['axes']) and event.button==3:
            print 'hello!!'
            line_at = int(np.round(event.xdata))
            ch_name = gui1.plot_channels[gui1.plot_dict['spec']['axes'].index(event.inaxes)]
        #self.clickable_axes
        print '####### line : ', line_at, ch_name
        ch_name = change_name(ch_name)
        gui1.action_list.append('cold_line:{}:{}'.format(ch_name, line_at))
        gui1.cold_line_select = False
        gui1.action_list_label.set_label('\n'.join(gui1.action_list))
        gui1.action_status.set_label('')
    if (event.inaxes in gui1.plot_dict['resid']['axes']) and event.button==3:
        gui1.x_val = gui1.t[np.argmin(np.abs(gui1.t - event.xdata))]
        print '****', gui1.x_val, event.xdata
        gui1.y_val = int(np.round(event.ydata))
        gui1.clicked()
    elif (event.inaxes in gui1.clickable_axes) and event.button==3:
    #elif (event.inaxes in gui1.clickable_axesgui1.plot_dict['temp']['axes']) and event.button==3:
        gui1.x_val = gui1.t[np.argmin(np.abs(gui1.t - event.xdata))]
        print '****', gui1.x_val, event.xdata
        gui1.y_val = int(np.round(event.ydata))
        gui1.clicked()
        #elif (event.inaxes in gui1.clickable_axes) and event.button==3:

gui1 = gui(shot)
cid = gui1.f.canvas.mpl_connect('button_press_event', onclick)
gtk.main()
#sd

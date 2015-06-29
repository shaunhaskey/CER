'''Script to ease the processing on GA computers and PPPL computers

Runs everything the whole way through including idl process at the start, and submitting the job to the PPPL cluster at the end

SRH : 30March2015
/u/stagnerl/FIDA/CGRAD/
'''
import pidly
import os, shutil
import time as time_mod
import numpy as np
import itertools as iter
import fida_funcs as FIDA
import cer_funcs as CER
#import matplotlib.pyplot as pt
import subprocess as sub
from scipy.io import netcdf
import cPickle as pickle

HOME = os.environ['HOME']
#Should I remove those directories before doing anything else?
#'/u/grierson/transp/ACFILE/blank_fi_1.cdf'

#Default pedestal parameters
#core = [0.02, 0.001]
HOST='venus' if os.environ['HOST'].find('venus')>=0 else 'portal'

#Open the file to write the new profiles to (for passing to IDL) -> move this to pIDLy?
model_dir = '/u/haskeysr/gaprofiles/model/'
model_shot = 158676; model_time = 3220
new_shot = 158676; base_time = 1400
single = False

#model_dir = '/u/haskeysr/gaprofiles/155196/haskeysr/'
##beam155196.02400
#model_shot = 155196; model_time = 2400
#new_shot = 155196; base_time = 500
#single = True

# model_dir = '/u/haskeysr/gaprofiles/model/'
# model_shot = 158676; model_time = 3220
# new_shot = 158676; base_time = 8000
# single = True

#Check for existing filenames
model_files = os.listdir(model_dir)
model_fnames = []
model_shot_str = '{}'.format(model_shot)
model_time_str = '{:05d}'.format(model_time)
file_prefix = ['dti','dte','dne','dtrot','dimp','beam']
master_dict = {'sims':{},'settings':{'model_dir':model_dir,'model_shot':model_shot,'new_shot':new_shot}}
for tmp_prefix in file_prefix:
    for i in model_files:
        if (i.find(tmp_prefix)==0) and (i.find(model_shot_str)>=0) and (i.find(model_time_str)>=0):
            model_fnames.append([tmp_prefix, i])
#idl_strs = []
widths = np.linspace(0.01, 0.12, 5)
te_val_top_list = np.linspace(0.5, 2.5, 5)
#0.5 -> 10 are the 'possible' ranges for ne at the top of the pedestal
ne_val_top_list = np.linspace(1, 8, 5)
ti_val_top_list = None

fnames = ['beam','dimp','dne','dte','dti','dtrot','g']
total_simuls = len(widths)*len(te_val_top_list)*len(ne_val_top_list)
for i in range(base_time, base_time+total_simuls):
    cmd = 'rm {home}/gaprofiles/{shot}/{ftype}{shot}.{time_string:05d}'.format(home=HOME,shot=new_shot,ftype='{'+','.join(fnames) + '}',time_string = i)
    #cmd = 'rm {home}/gaprofiles/{shot}/{ftype}{shot}.{range_string}'.format(home=HOME,shot=new_shot,ftype=i,range_string = '{'+'{:05}..{:05}'.format(base_time, base_time + total_simuls)+'}')
    print cmd
    os.system(cmd)
    cmd = 'rm -r {home}/gaprofiles/f90fidasim/{shot}/{time_string:05d}'.format(home=HOME,shot=new_shot, time_string = i)
    print cmd
    os.system(cmd)
    cmd = 'rm -r {home}/FIDASIM/RESULTS/D3D/{shot}/{time_string:05d}'.format(home=HOME,shot=new_shot, time_string = i)
    print cmd
    os.system(cmd)
#new_times = np.arange(len(widths)) + base_time
#ne_val_top_list = np.linspace(0.5,2.5,5)

shot_list, time_list = FIDA.generate_dirs(widths, base_time, master_dict, new_shot, model_fnames, model_dir, model_shot, model_time, single = single, te_val_top_list = te_val_top_list, ti_val_top_list = ti_val_top_list, ne_val_top_list = ne_val_top_list)

idl_string = FIDA.generate_profiles_nc(widths, base_time, master_dict, new_shot, model_fnames, model_time_str, single = single, te_val_top_list = te_val_top_list, ti_val_top_list = ti_val_top_list, ne_val_top_list = ne_val_top_list)

idl = pidly.IDL('idl')

# f1 = '/u/haskeysr/idl_test.pro'
# f2 = '/u/haskeysr/idl_test2.pro'
# with file(f1,'w') as filehandle:filehandle.write(idl_string)
# idl_text2 = '@{}\nexit\n'.format('/u/haskeysr/idl_test.pro')
# with file(f2,'w') as filehandle:
#     filehandle.write(idl_text2)
# FIDA.check_file_exists(f1, max_time = 10, interval = 0.1)
# FIDA.check_file_exists(f2, max_time = 10, interval = 0.1)
# time_mod.sleep(4)
# print '***===***'

idl_string_new = idl_string.replace('exit','')
for i in idl_string_new.split('\n'):
    print i
    idl(i)
print '***===***'
#os.system('idl < /u/haskeysr/idl_test2.pro')

#Finished with the profiles at this point, now move on to the actual FIDASIM part

#Create and idl process
grid_settings = {'xmin':35,'xmax':95,'nr_halo':5000000}
for i in master_dict['sims'].keys():master_dict['sims'][i]['prefida_changes'] = grid_settings
FIDA.generate_run_FIDASIM(master_dict, time_list, shot_list, grid_settings, idl = idl, HOST='venus', setpoint = 15)

1/0

import numpy as np
lens1 = np.array([-58.0452,238.6632,0.6822])
loc1 = np.array([-133.64747, 172.84416, -1.4165587])
lens2 = np.array([-58.045200, 238.66320, 0.68220000])
loc2 = np.array([-134.95745, 186.51464, -1.0956602])
dr_loc = loc1 - loc2
dr_lens = lens1 - lens2
n = 80
frac_vec = np.linspace(-0.25,2.2, n)

loc3 = [dr_loc[i]*frac_vec+loc2[i] for i in range(3)]
lens3 = [dr_lens[i]*frac_vec+lens2[i] for i in range(3)]
use_lens = lens3
use_loc = loc3
ids = ['m17','m24','m25','m26','m27','m28','m29','m31','m32']
ids = ['m17','m24','m25','m26','m27','m28','m29','m31','m32']
ids = ['m{}'.format(i) for i in range(1,n)]
lens = [];loc = []
axes = ['x','y','z']
for i in range(3):
    lens.append(np.linspace(lens1[i],lens2[i],n))
    loc.append(np.linspace(loc1[i],loc2[i],n))
for i in range(3):
    print '{}lens = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in use_lens[i]]))
for i in range(3):
    print '{}los = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in use_loc[i]]))
print 'mchords = [{}]'.format(','.join(["'{}'".format(i) for i in ids]))
fig, ax = pt.subplots()
ax.plot(lens3[0],lens3[1],'x')
ax.plot(loc3[0],loc3[1],'x')
ax.plot(lens[0],lens[1],'o')
ax.plot(loc[0],loc[1],'o')
phi = np.linspace(0,2.*np.pi)
ax.plot(171*np.cos(phi),171*np.sin(phi))
ax.set_xlim([-240,240])
ax.set_ylim([-240,240])
fig.canvas.draw();fig.show()



setup_profiles = r'''.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/haskeysr/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile writeg
.compile /u/haskeysr/gaprofiles_compute_profile_beam.pro

'''


'''
% Compiled module: BST_PARAM_WAVELENGTH.
% LOAD_ALL_CCD_DATA: Invalid chord identifier: m17
% GET_CCD_CHORD_DATA: Error in external call to LOAD_ALL_CCD_DATA
% Tag name WL is undefined for structure <Anonymous>.
% Execution halted at:  BST_CERVIEW_GET_LAMBDA0  512 /u/grierson/idlpros/b-stark/bst_cerview/bst_cerview.pro
%                       BST_CHORD_PARAM_WAVELENGTH  263 /u/grierson/idlpros/b-stark/bst_chord_param/bst_chord_param.pro
%                       BST_CHORD_PARAM   814 /u/grierson/idlpros/b-stark/bst_chord_param/bst_chord_param.pro
%                       GET_MAINION_GEOM   61 /u/haskeysr/FIDASIM/D3D/d3d_chords.pro
%                       D3D_CHORDS        276 /u/haskeysr/FIDASIM/D3D/d3d_chords.pro
%                       D3D_ROUTINES       21 /u/grierson/FIDASIM/D3D/d3d_routines.pro
%                       PREFIDA           939 /u/haskeysr/FIDASIM/prefida.pro
%                       $MAIN$          
% Tag name LENS is undefined for structure NULL.
% Execution halted at: GET_MAINION_GEOM   62
   /u/haskeysr/FIDASIM/D3D/d3d_chords.pro
%                      D3D_CHORDS        276
   /u/haskeysr/FIDASIM/D3D/d3d_chords.pro
%                      D3D_ROUTINES       21
   /u/grierson/FIDASIM/D3D/d3d_routines.pro
%                      PREFIDA           939 /u/haskeysr/FIDASIM/prefida.pro
%                      $MAIN$          
IDL>  
'''

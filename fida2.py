'''Script to ease the processing on GA computers and PPPL computers

Runs everything the whole way through including idl process at the start, and submitting the job to the PPPL cluster at the end

SRH : 30March2015
/u/stagnerl/FIDA/CGRAD/
'''
import os, shutil
import time as time_mod
import os
import numpy as np
import itertools as iter
import fida_funcs as FIDA
import cer_funcs as CER
import numpy as np
import matplotlib.pyplot as pt
import subprocess as sub
from scipy.io import netcdf
#Should I remove those directories before doing anything else?
#'/u/grierson/transp/ACFILE/blank_fi_1.cdf'

offset = 0.1
te_top = 3.0
ti_top = 3.0
ne_top = 4.0
xsym = 0.95
hwid = 0.05
#core = [0.02, 0.001]
core = [0.02]#, 0.001]
edge = []#-0.02]

if os.environ['HOST'].find('venus')>=0: 
    HOST='venus'
else:
    HOST='portal'
#Open the file to write the new profiles to (for passing to IDL) -> move this to pIDLy?
f = netcdf.netcdf_file('/u/haskeysr/test.nc','w')
idl_string = FIDA.idl_header
model_dir = '/u/haskeysr/gaprofiles/model/'
model_shot = 158676
model_time = 3220
new_shot = 158676
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
idl_strs = []
widths = np.linspace(0.01,0.12,5)
base_time = 900
#new_times = np.arange(len(widths)) + base_time
shot_list = []; time_list = []
max_val = 1.21
npts = 121
count = 0

te_val_top_list = np.linspace(0.5,2.5,5)
#ne_val_top_list = np.linspace(0.5,2.5,5)
ne_val_top = 4.; te_val_top = 3.; ti_val_top = 3.
for i, width in enumerate(widths):
    for te_val_top in te_val_top_list:
        new_time = count + base_time
        master_dict['sims'][new_time] = {'dir_dict':{},'shot':new_time, 'profiles':{}}
        print new_time, width
        #cur_dir = '/u/haskeysr/gaprofiles/{}/{:05d}/'.format(new_shot,new_time)
        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)
        master_dict['sims'][new_time]['dir_dict']['profiles'] = cur_dir
        os.system('mkdir -p {}'.format(cur_dir))
        #Copy the model shot but with the modified filenames. Data inside the files will be modified by idl soon....
        for (tmp_prefix, tmp_fname) in model_fnames:
            new_name = tmp_fname.replace(model_shot_str,'{}'.format(new_shot)).replace(model_time_str,'{:05d}'.format(new_time))
            print tmp_prefix, tmp_fname, new_name
            os.system('cp {}/{} {}/{}'.format(model_dir, tmp_fname, cur_dir,new_name ))
        #Do the same for the g-file
        os.system('cp {}/{} {}/{}'.format(model_dir, 'g{}.{:05d}'.format(model_shot,model_time), cur_dir,'g{}.{:05d}'.format(new_shot,new_time) ))

        idl_string+=FIDA.idl_ind.format(dir=cur_dir,shot=new_shot,time=new_time,time_str='{:05d}'.format(new_time),tag='{}'.format(new_time))
        idl_string+=FIDA.idl_bulk
        ti_val_top = te_val_top
        for ident, width_mult, top_val in zip(['ne', 'te', 'ti'], [1., 1., 2.], [ne_val_top, te_val_top, ti_val_top]):
        #Generate the profiles
            print ident, width_mult, top_val
            x = np.linspace(0,max_val, npts)
            y = CER.mtanh(x, top_val, offset, xsym, width*width_mult, len(core), len(edge), False, *(core + edge))
            #y = CER.mtanh(x, top_val, offset, xsym, width*width_mult, core = core, edge = edge)
            cur_var = '{}{}'.format(ident,new_time)
            f.createDimension(cur_var,len(x))
            prof = f.createVariable(cur_var,'float',(cur_var,))
            prof[:] = +y
            master_dict['sims'][new_time]['profiles'][ident] = {'data_y':+y,'data_x':+x,'settings':{'width':width,'xsym':xsym,'offset':offset,'pedestal_top':top_val,'core':core,'npts':npts, 'max_val':max_val,'edge':edge}}

        if count==0:
            f.createDimension('rho',len(x))
            rho = f.createVariable('rho','float',('rho',))
            rho[:] = +x
            rho.units = ''
        shot_list.append('{}'.format(new_shot))
        time_list.append('{:05d}'.format(new_time))
        count += 1
        #idl_strs.append('''gap_model_dir = '{}'\ngap_model_shot = '{}'\ngap_model_time = {}\n'''.format(cur_dir,new_shot,new_time))

f.close()
with file('/u/haskeysr/idl_test.pro','w') as filehandle:filehandle.write(idl_string)
idl_text2 = '@{}\nexit\n'.format('/u/haskeysr/idl_test.pro')
with file('/u/haskeysr/idl_test2.pro','w') as filehandle:
    filehandle.write(idl_text2)
time_mod.sleep(4)
os.system('idl < /u/haskeysr/idl_test2.pro')

#Finished with the profiles at this point, now move on to the actual FIDASIM part

#Settings
diag = 'MAIN_ION330'
beam = '330lt'
comment = 'helloworld'
master_dict['settings']['diag'] = diag
master_dict['settings']['beam'] = beam
master_dict['settings']['comment'] = comment

#Important directories
pppl_base_dir = r'/p/fida/shaskey/RESULTS/D3D/'.replace('//','/')
d3d_base_dir = r'/u/haskeysr/FIDASIM/RESULTS/D3D/'.replace('//','/')

master_dict['settings']['pppl_base_dir'] = pppl_base_dir
master_dict['settings']['d3d_base_dir'] = d3d_base_dir

1/0
def create_dirs(shot_list, time_list, diag):
    pppl_actual_dir_list = []
    d3d_actual_dir_list = []
    for shot,time in zip(shot_list, [int(tmp1) for tmp1 in time_list]):
        pppl_actual_dir = r'{}/{}/{:05d}/{}/'.format(pppl_base_dir,shot, time, diag).replace('//','/') 
        d3d_actual_dir = r'{}/{}/{:05d}/{}/'.format(d3d_base_dir,shot, time, diag).replace('//','/')
        pppl_actual_dir_list.append(pppl_actual_dir)
        d3d_actual_dir_list.append(d3d_actual_dir)

        if os.path.exists(d3d_actual_dir):
            os.system('rm -r {}'.format(d3d_actual_dir))

        tmp = 'mkdir -p {}'.format(d3d_actual_dir)
        master_dict['sims'][time]['dir_dict']['d3d_actual_dir'] = d3d_actual_dir
        master_dict['sims'][time]['dir_dict']['pppl_actual_dir'] = pppl_actual_dir
        print '==> making the directory, ', tmp
        os.system(tmp)
    return pppl_actual_dir_list, d3d_actual_dir_list

def make_run_idl(shot_list, time_list, diag):
    idl_input_file = r'{}/{}/idl_input_file'.format(d3d_base_dir,shot_list[0]).replace('//','/')
    idl_output_log = r'{}/{}/idl_output_log'.format(d3d_base_dir,shot_list[0]).replace('//','/')
    f90_overall = ''
    d3d_prefida_overall = ''
    for shot, time in zip(shot_list, time_list):
        time_str = '{:05d}'.format(int(time))
        f90_overall+=FIDA.f90_txt.format(shot=shot, time = time_str, beam=beam,comment=comment, diag = diag)
        d3d_prefida_overall+=FIDA.d3d_input_prefida.format(shot=shot, time = time_str,)

    print f90_overall
    print d3d_prefida_overall
    
    with file(idl_input_file,'w') as filehandle:
        tmp_txt = FIDA.setup_txt2.format(f90_txt = f90_overall, d3d_input_prefida = d3d_prefida_overall)
        print tmp_txt
        filehandle.write(tmp_txt)

    idl_text2 = "REPEAT BEGIN\nresult = FILE_TEST('{}')\nprint,result\nwait,1\nENDREP UNTIL result EQ 1\n@{}\nexit\n".format(idl_input_file,idl_input_file)
    idl_input_file2 = idl_input_file+'2'
    with file(idl_input_file2,'w') as filehandle:
        filehandle.write(idl_text2)
    #Run IDL to set everything up
    file_is_there = False
    time_mod.sleep(10)
    while not os.path.isfile(idl_input_file):
        time_mod.sleep(0.01)
        print 'hello'

    while not os.path.isfile(idl_input_file2):
        time_mod.sleep(0.01)
        print 'hello2'
    time_mod.sleep(2)
        
    #time_mod.sleep(2)
    os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))

def write_job_file(d3d_actual_dir_list, pppl_actual_dir_list, job_id):
    if HOST=='venus':
        job_template = FIDA.pbs_file_venus
        run_dir_list = d3d_actual_dir_list
    else:
        job_template = FIDA.pbs_file
        run_dir_list = pppl_actual_dir_list
    #Write the jobfile for qsub
    count = 0
    for d3d_actual_dir, run_dir in zip(d3d_actual_dir_list, run_dir_list):
        with file('{}/jobfile.pbs'.format(d3d_actual_dir),'w') as filehandle:
            filehandle.write(job_template.format(PPPL_dir=run_dir,input_file='def_inputs.dat',tmp_work_name='fida_{:d}'.format(count),job_id=job_id))
        wrapper_fname = '{}/wrapper'.format(d3d_actual_dir)
        if HOST=='portal':
            with file(wrapper_fname,'w') as filehandle:
                filehandle.writelines(FIDA.wrapper.format(pppl_dir = run_dir))
            os.system('chmod +x {}'.format(wrapper_fname))
        count+=1

def modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list):
    if HOST=='venus':
        dir_list = d3d_actual_dir_list
    else:
        dir_list = pppl_actual_dir_list
    for d3d_actual_dir, actual_dir in zip(d3d_actual_dir_list, dir_list):
        with file('{}def_inputs.dat'.format(d3d_actual_dir),'r') as filehandle:
            fidasim_input = filehandle.readlines()
        for i,line in enumerate(fidasim_input):
            if line.find("result_dir")>=0:
                line_num = +i
        fidasim_input[line_num] = r"result_dir = '{}'".format(actual_dir)
        with file('{}/def_inputs.dat'.format(d3d_actual_dir),'w') as filehandle:
            filehandle.writelines(fidasim_input)

def copy_files(shot):
    #Copy everything across
    os.system('rsync -avz --progress /u/haskeysr/FIDASIM/RESULTS/D3D/{} shaskey@portal.pppl.gov:{}'.format(shot,pppl_base_dir))

def execute(pppl_actual_dir_list, HOST):
    #Submit the job on the PPPL cluster
    if HOST=='portal':
        dir_list = pppl_actual_dir_list
        cmd_template = r"ssh -t shaskey@portal.pppl.gov '{}/wrapper'"
    else:
        dir_list = d3d_actual_dir_list
        cmd_template = r"qsub {}/jobfile.pbs"
    for actual_dir in dir_list:
        #cmd = r"ssh -t shaskey@portal.pppl.gov '{}/wrapper'".format(pppl_actual_dir).replace('//','/')
        cmd = cmd_template.format(actual_dir).replace('//','/')
        print '==>',cmd
        os.system(cmd)

pppl_actual_dir_list, d3d_actual_dir_list = create_dirs(shot_list, time_list,diag)
make_run_idl(shot_list, time_list, diag)
job_id = 'fidasim'
write_job_file(d3d_actual_dir_list, pppl_actual_dir_list, job_id)
if HOST=='portal':
    modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list)
#What to do about the results that may already exist in the remote directory?
if HOST=='portal': 
    for i in set(shot_list):copy_files(i)

fida_runs_fname = '/u/haskeysr/fida_runs'
setpoint = 13
with file(fida_runs_fname,'w') as filehandle: filehandle.write('{}\n'.format(setpoint))
import cPickle as pickle
pickle.dump(master_dict,file('/u/haskeysr/fida_sim_dict.pickle','w'))
FIDA.batch_launch_fida(master_dict['sims'], fida_runs_fname, setpoint = setpoint, id_string = job_id)


1/0


import numpy as np
lens1 = np.array([-58.0452,238.6632,0.6822])
loc1 = np.array([-133.64747, 172.84416, -1.4165587])
lens2 = np.array([-58.045200, 238.66320, 0.68220000])
loc2 = np.array([-134.95745, 186.51464, -1.0956602])
dr_loc = loc1 - loc2
dr_lens = lens1 - lens2
n = 120
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

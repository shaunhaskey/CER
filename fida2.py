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

#Should I remove those directories before doing anything else?

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

setup_txt2 = '''.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/haskeysr/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim

{f90_txt}
.compile /u/haskeysr/FIDASIM/prefida
.compile /u/haskeysr/FIDASIM/D3D/d3d_beams
.compile /u/haskeysr/FIDASIM/D3D/d3d_chords
.compile /u/haskeysr/FIDASIM/D3D/d3d_profiles
{d3d_input_prefida}
exit
'''

f90_txt = '''f90fidasim_setup,{shot},{time},'{diag}','{beam}',DIR='/u/haskeysr/gaprofiles/{shot}/',COMMENT='{comment}',EINJ=75.0,PINJ=2.0\n'''
#'/u/grierson/transp/ACFILE/blank_fi_1.cdf'

d3d_input_prefida = '''.compile /u/haskeysr/gaprofiles/f90fidasim/{shot}/{time}/MAIN_ION330/def/d3d_input
prefida,'input_template'\n'''

setup_txt = r''';spawn,'rm -r /u/haskeysr/gaprofiles/f90fidasim/{shot}'
;spawn,'rm -r /u/haskeysr/FIDASIM/RESULTS/D3D/{shot}'
;@/e/fida/FIDASIM/D3D/d3d_startup.pro; This can probably replace stuff below
print,'################################'
.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
get_chord,155196,'t01'
addf90fidasim
setenv,'FIDASIM_DIR=/u/haskeysr/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile /u/haskeysr/f90fidasim_setup.pro

f90fidasim_setup,{shot},{time},'{diag}','{beam}',DIR='/u/haskeysr/gaprofiles/{shot}/',COMMENT='{comment}',cdf_file='/u/grierson/transp/ACFILE/blank_fi_1.cdf'
.compile /u/haskeysr/FIDASIM/prefida
.compile /u/haskeysr/FIDASIM/D3D/d3d_beams
.compile /u/haskeysr/gaprofiles/f90fidasim/{shot}/{time}/MAIN_ION330/def/d3d_input
prefida,'input_template'
exit
'''

#This is everything for the qsub job
pbs_file='''#!/bin/tcsh
# PBS script for f90fidasim
# One CPU for 4 hours
#PBS -l nodes=1
#PBS -l walltime=24:00:00
#PBS -l mem=2gbb
#PBS -N f90fidasim
##PBS -e pbs.err
##PBS -o pbs.log
#PBS -r n
#PBS -M shaskey@pppl.gov
#PBS -m aeb

echo -n 'Started job at : ' ; date

limit stacksize unlimited

cd {PPPL_dir}
rm log
unbuffer /u/shaskey/bin/FIDASIM {input_file} >& log
echo -n 'Ended job at : ' ; date
echo " "
exit
'''

#This is everything for the qsub job
pbs_file_venus='''#!/bin/tcsh
#$ -N {job_id}
#$ -q all.q
#$ -o sge_output.dat
#$ -e sge_error.dat
#$ -cwd
#$ -l mem_free=4G,h_vmem=4G
setenv FIDASIM_DIR /u/haskeysr/FIDASIM/
limit stacksize unlimited
#Maybe store these environment variables before changing them, so they can be restored after fidasim runs
setenv LD_LIBRARY_PATH /task/imd/local64/lib
setenv LD_LIBRARY_PATH /task/imd/local64/include:$LD_LIBRARY_PATH
setenv NETCDF_INCLUDE /task/imd/local64/include
setenv NETCDF_LIB /task/imd/local64/lib



#cd {PPPL_dir}
#Remove, make, copy, goto, run, copy back
rm /tmp/{tmp_work_name}/def_*.cdf
mkdir /tmp/{tmp_work_name}
cp -a {PPPL_dir}* /tmp/{tmp_work_name}/
cd /tmp/{tmp_work_name}/
rm log
unbuffer /u/haskeysr/FIDASIM/fidasim {input_file} >& log
echo -n 'Ended job at : ' ; date
echo " "
echo -n 'Copying files back to proper place'
cp -a /tmp/{tmp_work_name}/* {PPPL_dir}
echo -n 'Finished'
exit
'''

#This is a wrapper around qsub to allow remote execution
wrapper='''#!/bin/bash -l
cd {pppl_dir}
qsub {pppl_dir}jobfile.pbs

'''

idl_header='''.compile writeg
file = '/u/haskeysr/test.nc'
fileID = ncdf_open(file, /NOWRITE)

'''

idl_bulk = '''
dne_file=dir+'dne'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dne_file
ne_str.rho_dens=rho
ne_str.dens=dens
dne_newfile=dir+'dne'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
PRINT,dne_newfile
SAVE,ne_str,FILE=dne_newfile,/VERB

dte_file=dir+'dte'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dte_file
te_str.rho_te=rho
te_str.te=te
dte_newfile=dir+'dte'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
PRINT,dte_newfile
SAVE,te_str,FILE=dte_newfile,/VERB

dti_file=dir+'dti'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dti_file
ti_str.rho_ti=rho
ti_str.ti=ti
dti_newfile=dir+'dti'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
PRINT,dti_newfile
SAVE,ti_str,FILE=dti_newfile,/VERB

dtrot_file=dir+'dtrot'+STRTRIM(shot,2)+'.'+ STRING(time,FOR='(I05)')
RESTORE,dtrot_file
tor_rot_str.tor_rot_local *= 0 ;tor_rot_local ;; Omega = V/R
tor_rot_str.v_tor_rot *= 0 ;tor_rot_local*tor_rot_str.r_tor_rot ;; V = Omega*R
dtrot_newfile=dir+'dtrot'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
PRINT,dtrot_newfile
SAVE,tor_rot_str,FILE=dtrot_newfile,/VERB

g=READG(g_file,MODE='FILE')
g.shot=shot
g.time=time
WRITEG_FILE,g,FILE=g_file_new

'''

idl_ind='''dir = '{dir}'
shot = '{shot}'
time = {time}
g_file = '{dir}/g{shot}.{time:05d}'
g_file_new = '{dir}/g{shot}.{time:05d}'

rho_tag=NCDF_VARID(fileID,'rho')
NCDF_VARGET,fileID, rho_tag, rho

ne_tag=NCDF_VARID(fileID,'ne{tag}')
NCDF_VARGET,fileID, ne_tag, dens

te_tag=NCDF_VARID(fileID,'te{tag}')
NCDF_VARGET,fileID, te_tag, te

ti_tag=NCDF_VARID(fileID,'ti{tag}')
NCDF_VARGET,fileID, ti_tag, ti

'''

import cer_funcs as CER
import numpy as np
import matplotlib.pyplot as pt
import subprocess as sub
offset = 0.1

te_top = 3.0
ti_top = 3.0
ne_top = 4.0

xsym = 0.95
hwid = 0.05

#core = [0.02, 0.001]
core = [0.02]#, 0.001]
edge = []#-0.02]

from scipy.io import netcdf


def check_job_number_file(file_name):
    file_name = open(file_name,'r')
    simul_jobs_content = file_name.read()
    file_name.close()
    stab_lines = simul_jobs_content.split('\n')
    simul_jobs = int(stab_lines[0].rstrip('\n'))
    return simul_jobs
def launch_job(job_name):
    os.system('qsub ' + job_name)

def running_jobs(id_string):
    #p = sub.Popen(['qstat', '-u', username],stdout=sub.PIPE,stderr=sub.PIPE)
    p = sub.Popen(['qstat'],stdout=sub.PIPE,stderr=sub.PIPE)
    output, errors = p.communicate()
    number_of_jobs = output.count(id_string)
    #number_of_jobs = output.count(count_string) - 2
    #if number_of_jobs<0:number_of_jobs = 0
    #print 'Running Jobs : ' + str(number_of_jobs)
    return number_of_jobs

def batch_launch_fida(master_dict, job_num_filename, setpoint = 5, id_string = 'fida'):
    start_time = time_mod.time()
    submitted_jobs = 0; startup = 1; startup_runs = 0

    total_jobs = len(master_dict.keys())

    for i in master_dict.keys():
        setpoint = check_job_number_file(job_num_filename)
        print i
        if startup==1:
            launch_job(master_dict[i]['dir_dict']['d3d_actual_dir'] + 'jobfile.pbs')
            startup_runs += 1
            if startup_runs == setpoint:
                startup = 0
                time_mod.sleep(20)
        else:
            while running_jobs(id_string)>setpoint:
                time_mod.sleep(10)
                setpoint = check_job_number_file(job_num_filename)
            launch_job(master_dict[i]['dir_dict']['d3d_actual_dir'] + 'jobfile.pbs')
        #master_dict[i]['chease_run'] = 1
        submitted_jobs += 1
        print 'Submitted %d of %d, Waiting for last %d jobs to finish, time so far : %.2fmins'%(submitted_jobs, total_jobs, running_jobs(id_string), (time_mod.time()-start_time)/60)
        print 'Submitted %d jobs of %d, time %.2fmins'%(submitted_jobs,total_jobs,(time_mod.time()-start_time)/60)
    while running_jobs(id_string)>0:
        print 'Submitted %d of %d, Waiting for last %d jobs to finish, time so far : %.2fmins'%(submitted_jobs, total_jobs, running_jobs(id_string), (time_mod.time()-start_time)/60)
        time_mod.sleep(15)
    return master_dict

if os.environ['HOST'].find('venus')>=0: 
    HOST='venus'
else:
    HOST='portal'

f = netcdf.netcdf_file('/u/haskeysr/test.nc','w')
idl_string = idl_header
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
base_time = 800
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

        idl_string+=idl_ind.format(dir=cur_dir,shot=new_shot,time=new_time,time_str='{:05d}'.format(new_time),tag='{}'.format(new_time))
        idl_string+=idl_bulk
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


'''
profile_ranges = {'NE':[[0.05], np.linspace(1,8,n).tolist(), [8.0]],
                  'TE':[[0.05], np.linspace(0.5,3,n).tolist(), [3.0]],
                  'TI':[[0.05], [1.], [2.0]],
                  'TROT':[[0.], [0.], [0.]]}
profiles = {}
for i in profile_ranges.keys():profiles[i] = list(iter.product(*profile_ranges[i]))
key_ordering = profiles.keys()
unique_profiles = list(iter.product(*[profiles[i] for i in key_ordering]))
print key_ordering
print unique_profiles
if force_ti_te_equal:
    key_ordering.append('TI')
    te_key = key_ordering.index('TE')
    unique_profiles = [list(i) for i in unique_profiles]

    for i in range(len(unique_profiles)): unique_profiles[i].append(unique_profiles[i][te_key])
print unique_profiles
'''


with file('/u/haskeysr/idl_test.pro','w') as filehandle:filehandle.write(idl_string)
f.close()
time_mod.sleep(4)
idl_text2 = '@{}\nexit\n'.format('/u/haskeysr/idl_test.pro')
#idl_input_file2 = idl_input_file+'2'
with file('/u/haskeysr/idl_test2.pro','w') as filehandle:
    filehandle.write(idl_text2)
os.system('idl < /u/haskeysr/idl_test2.pro')

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

def create_dirs(shot_list, time_list, diag):

    #pppl_actual_dir_list = [r'{}/{}/{:05d}/{}/'.format(pppl_base_dir,shot, int(time), diag).replace('//','/') for shot,time in zip(shot_list, time_list)]
    #d3d_actual_dir_list = [r'{}/{}/{:05d}/{}/'.format(d3d_base_dir,shot, int(time), diag).replace('//','/') for shot, time in zip(shot_list, time_list)]
    #for pppl_actual_dir, d3d_actual_dir in zip(pppl_actual_dir_list, d3d_actual_dir_list):
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
        f90_overall+=f90_txt.format(shot=shot, time = time_str, beam=beam,comment=comment, diag = diag)
        d3d_prefida_overall+=d3d_input_prefida.format(shot=shot, time = time_str,)

    print f90_overall
    print d3d_prefida_overall
    
    with file(idl_input_file,'w') as filehandle:
        tmp_txt = setup_txt2.format(f90_txt = f90_overall, d3d_input_prefida = d3d_prefida_overall)
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
        job_template = pbs_file_venus
        run_dir_list = d3d_actual_dir_list
    else:
        job_template = pbs_file_venus
        run_dir_list = pppl_actual_dir_list
    #Write the jobfile for qsub
    count = 0
    for d3d_actual_dir, run_dir in zip(d3d_actual_dir_list, run_dir_list):
        with file('{}/jobfile.pbs'.format(d3d_actual_dir),'w') as filehandle:
            filehandle.write(job_template.format(PPPL_dir=run_dir,input_file='def_inputs.dat',tmp_work_name='fida_{:d}'.format(count),job_id=job_id))
        wrapper_fname = '{}/wrapper'.format(d3d_actual_dir)
        if HOST=='portal':
            with file(wrapper_fname,'w') as filehandle:
                filehandle.writelines(wrapper.format(pppl_dir = run_dir))
            os.system('chmod +x {}'.format(wrapper_fname))

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
with file(fida_runs_fname,'w') as filehandle: filehandle.write('13\n')
import cPickle as pickle

pickle.dump(master_dict,file('/u/haskeysr/fida_sim_dict.pickle','w'))
batch_launch_fida(master_dict['sims'], fida_runs_fname, setpoint = 5, id_string = job_id)


1/0
1/0
execute(pppl_actual_dir_list, HOST)


1/0


ref_shot = 158676; ref_time = 3220; new_shot = 158676
run_information = {}
n = 5
dnep_vals = np.linspace(0.2,0.6,n)
count = 300
force_ti_te_equal = True
profile_ranges = {'NE':[[0.05], np.linspace(1,8,n).tolist(), [8.0]],
                  'TE':[[0.05], np.linspace(0.5,3,n).tolist(), [3.0]],
                  'TI':[[0.05], [1.], [2.0]],
                  'TROT':[[0.], [0.], [0.]]}
dNE = 1.
dTE = 3.5
dTI = 3.5
dTROT = False

#TI, TE in keV
#ne in 10^13 cm**-3
if force_ti_te_equal:del(profile_ranges['TI'])

profiles = {}
for i in profile_ranges.keys():profiles[i] = list(iter.product(*profile_ranges[i]))
key_ordering = profiles.keys()
unique_profiles = list(iter.product(*[profiles[i] for i in key_ordering]))
print key_ordering
print unique_profiles
if force_ti_te_equal:
    key_ordering.append('TI')
    te_key = key_ordering.index('TE')
    unique_profiles = [list(i) for i in unique_profiles]

    for i in range(len(unique_profiles)): unique_profiles[i].append(unique_profiles[i][te_key])
print unique_profiles


for tmp in unique_profiles:
    arg_string = []
    for dat, key_name in zip(tmp,key_ordering):
        dat = list(dat)
        modded = False
        if eval('d{}'.format(key_name))!=False:
            dat[2] = dat[1] + eval('d{}'.format(key_name))
            modded = True
        print key_name, dat, modded
        arg_string.append('D{key}EDGE={},D{key}PED={},D{key}CORE={}'.format(*dat,key=key_name))
        run_information[count] = {'D{}EDGE'.format(key_name):dat[0],'D{}PED'.format(key_name):dat[1],'D{}CORE'.format(key_name):dat[2]}
    profile_args = ','.join(arg_string)
    setup_profiles+= "gaprofiles_modify_profiles,{},{},{},{},{},/write\n".format(new_shot,count,ref_shot,ref_time, profile_args)
    #setup_profiles+="gaprofiles_modify_profiles,{},{},{},{},DNEEDGE=0.05,DNEPED={:.4f},/write\n".format(new_shot,count,ref_shot,ref_time,dnep)
    #run_information[count] = {}
    count+=1
#setup_txt+='exit\n'
print setup_profiles
clean_gaprofiles_dir = True
clean_f90_dir = True
ga_profiles_dir = '/u/haskeysr/gaprofiles/{}'.format(new_shot)

if clean_gaprofiles_dir:shutil.rmtree(ga_profiles_dir,ignore_errors=True)
if not os.path.exists(ga_profiles_dir):os.makedirs(ga_profiles_dir)
    
idl_input_file = 'test.txt'
idl_input_file2 = 'test2.txt'
idl_output_log = 'test.log'
with file(idl_input_file,'w') as filehandle:
    filehandle.write(setup_profiles)
idl_text2 = '@{}\nexit\n'.format(idl_input_file)
with file(idl_input_file2,'w') as filehandle:
    filehandle.write(idl_text2)
#Run IDL to set everything up
#1/0
time_mod.sleep(1.)
os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))

use_new_ones = True
a = os.listdir(ga_profiles_dir)
time_list = []
shot_list = []
for i in a:
    if i[0] == 'g':# and i[1:1+6]=='{}'.format(shot):
        time_tmp = int(i.split('.')[-1])
        if ((use_new_ones) and (time_tmp in run_information.keys())):
            time_list.append(i.split('.')[-1])
            shot_list.append(i[1:7])
        elif use_new_ones==False:
            time_list.append(i.split('.')[-1])
            shot_list.append(i[1:7])

print time_list
print shot_list

#count = 0
#os.system('rm -r 
new_shot_dir = '/u/haskeysr/gaprofiles/f90fidasim/{shot}'.format(shot=new_shot)
if clean_f90_dir:shutil.rmtree(new_shot_dir,ignore_errors=True)

time_str_list, pppl_actual_dir_list, d3d_actual_dir_list = create_dirs(shot_list, time_list,diag)

make_run_idl(shot_list, time_list, diag)

write_job_file(d3d_actual_dir_list, pppl_actual_dir_list)

modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list)

#What to do about the results that may already exist in the remote directory?
for i in set(shot_list):copy_files(i)
execute(pppl_actual_dir_list, HOST)

1/0

for shot, time in zip(shot_list, time_list):
    print shot, time
    #if count!=0:
    run_case(int(shot), int(time), diag)
    count +=1
    if count==6:
        break

prefix = ['dne','dti','dte','dtrot']



'''
spawn,'rm -r /u/haskeysr/gaprofiles/f90fidasim/158673'
.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
get_chord,155196,'t01'
addf90fidasim
setenv,'FIDASIM_DIR=/u/grierson/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile /u/haskeysr/f90fidasim_setup.pro
f90fidasim_setup,158673,2,'MAIN_ION330','330lt',DIR='/u/haskeysr/gaprofiles/158673/',COMMENT='helloworld',cdf_file='/u/grierson/transp/ACFILE/blank_fi_1.cdf'


;f90fidasim_setup,900003,00004,'MAIN_ION330','330lt',DIR='/u/grierson/gaprofiles/900003/transp/00004',COMMENT='helloworld'
print,'################################'
.compile /u/grierson/FIDASIM/prefida
print,'################################'
.compile /u/grierson/FIDASIM/D3D/d3d_beams
print,'################################'
.compile /u/haskeysr/gaprofiles/f90fidasim/900003/00004/MAIN_ION330/def/d3d_input
.compile /u/haskeysr/gaprofiles/f90fidasim/158673/00002/MAIN_ION330/def/d3d_input
print,'################################'
prefida,'input_template'
;Need to loop over the last two things to do multiple times and shots?
print,'################################'
exit
'''

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
print ''
print ids
for i in range(3):
    print '{}lens = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in lens[i]]))
for i in range(3):
    print '{}los = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in loc[i]]))



import numpy as np
lens1 = np.array([-58.0452,238.6632,0.6822])
loc1 = np.array([-133.64747, 172.84416, -1.4165587])
lens2 = np.array([-58.045200, 238.66320, 0.68220000])
loc2 = np.array([-134.95745, 186.51464, -1.0956602])
dr_loc = loc1 - loc2
dr_lens = lens1 - lens2
n = 120
frac_vec = np.linspace(-0.25,2.5, n)

loc3 = [dr_loc[i]*frac_vec+loc2[i] for i in range(3)]
lens3 = [dr_lens[i]*frac_vec+lens2[i] for i in range(3)]
ids = ['m17','m24','m25','m26','m27','m28','m29','m31','m32']
ids = ['m17','m24','m25','m26','m27','m28','m29','m31','m32']
ids = ['m{}'.format(i) for i in range(1,n)]
lens = [];loc = []
axes = ['x','y','z']
for i in range(3):
    lens.append(np.linspace(lens1[i],lens2[i],n))
    loc.append(np.linspace(loc1[i],loc2[i],n))
for i in range(3):
    print '{}lens = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in lens[i]]))
for i in range(3):
    print '{}los = [{}]'.format(axes[i], ','.join(['{:.4f}'.format(j) for j in loc[i]]))
print 'mchords = [{}]'.format(','.join(["'{}'".format(i) for i in ids]))
fig, ax = pt.subplots()
ax.plot(lens3[0],lens3[1],'x')
ax.plot(loc3[0],loc3[1],'x')
ax.plot(lens[0],lens[1],'o')
ax.plot(loc[0],loc[1],'o')
ax.set_xlim([-240,240])
ax.set_ylim([-240,240])
fig.canvas.draw();fig.show()

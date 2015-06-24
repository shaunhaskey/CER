import OMFITtree
import pidly
import os
import time as time_mod
import subprocess as sub
import cPickle as pickle
from scipy.io import netcdf
import numpy as np
import cer_funcs as CER
import scipy
import shutil
import itertools as iter
import cer_funcs as CER

idl_setup = '''.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/haskeysr/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim

'''

idl_setup2 = '''.compile /u/haskeysr/FIDASIM/prefida
.compile /u/haskeysr/FIDASIM/D3D/d3d_beams
.compile /u/haskeysr/FIDASIM/D3D/d3d_chords
.compile /u/haskeysr/FIDASIM/D3D/d3d_profiles

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

d3d_input_prefida = '''.compile /u/haskeysr/gaprofiles/f90fidasim/{shot}/{time}/MAIN_ION330/def/d3d_input
prefida,'input_template'\n'''


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

idl_bulk_single = '''
dne_file=dir+'dne'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dne_file
dne_newfile=dir+'dne'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
SAVE,ne_str,FILE=dne_newfile,/VERB

dte_file=dir+'dte'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dte_file
dte_newfile=dir+'dte'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
SAVE,te_str,FILE=dte_newfile,/VERB

dti_file=dir+'dti'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
RESTORE,dti_file
dti_newfile=dir+'dti'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
SAVE,ti_str,FILE=dti_newfile,/VERB

dtrot_file=dir+'dtrot'+STRTRIM(shot,2)+'.'+ STRING(time,FOR='(I05)')
RESTORE,dtrot_file
dtrot_newfile=dir+'dtrot'+STRTRIM(shot,2)+'.'+STRING(time,FOR='(I05)')
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

idl_ind_single='''dir = '{dir}'
shot = '{shot}'
time = {time}
g_file = '{dir}/g{shot}.{time:05d}'
g_file_new = '{dir}/g{shot}.{time:05d}'

'''

idl_header='''.compile writeg
file = '/u/haskeysr/test.nc'
fileID = ncdf_open(file, /NOWRITE)

'''


################ JOB LAUNCHING DETAILS ###################
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

################ JOB LAUNCHING DETAILS ###################
##
def check_file_exists(fname, max_time = 10, interval = 0.1):
    max_count = max_time / interval
    count = 0
    while (not os.path.isfile(fname)) and (count <max_count):
        time_mod.sleep(interval)
        print 'hello', fname, count, max_count
        count += 1

def create_dirs(shot_list, time_list, diag, pppl_base_dir, d3d_base_dir, master_dict):
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


def make_run_idl_f90setup(shot_list, time_list, diag, d3d_base_dir, beam, comment, idl = None):
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
        tmp_txt = setup_txt2.format(f90_txt = f90_overall, d3d_input_prefida = '\n')
        print tmp_txt
        filehandle.write(tmp_txt)
    if idl!=None:
        for i in idl_setup.split('\n'):
            print i
            idl(i)
        for shot, time in zip(shot_list, time_list):
            time_str = '{:05d}'.format(int(time))
            idl.pro('f90fidasim_setup',shot,time,diag,beam,DIR='/u/haskeysr/gaprofiles/{shot}/'.format(shot=shot),COMMENT=comment,EINJ=75.0,PINJ=2.0)
    else:
        idl_text2 = "REPEAT BEGIN\nresult = FILE_TEST('{}')\nprint,result\nwait,1\nENDREP UNTIL result EQ 1\n@{}\nexit\n".format(idl_input_file,idl_input_file)
        idl_input_file2 = idl_input_file+'2'
        with file(idl_input_file2,'w') as filehandle:
            filehandle.write(idl_text2)
        #Run IDL to set everything up
        check_file_exists(idl_input_file, max_time = 10, interval = 0.1)
        check_file_exists(idl_input_file2, max_time = 10, interval = 0.1)
        time_mod.sleep(4)

        #time_mod.sleep(2)
        os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))



def make_run_idl_prefida(shot_list, time_list, diag, d3d_base_dir, beam, comment, idl = None):
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
    if idl!=None:
        for i in idl_setup2.split('\n'):
            print i
            idl(i)
        for shot, time in zip(shot_list, time_list):
            time_str = '{:05d}'.format(int(time))
            tmp_str = '.compile /u/haskeysr/gaprofiles/f90fidasim/{shot}/{time}/MAIN_ION330/def/d3d_input'.format(shot = shot, time = time_str)
            idl(tmp_str)
            idl.pro('prefida','input_template')
    else:
        with file(idl_input_file,'w') as filehandle:
            tmp_txt = setup_txt2.format(f90_txt = '\n', d3d_input_prefida = d3d_prefida_overall)
            print tmp_txt
            filehandle.write(tmp_txt)
        idl_text2 = "REPEAT BEGIN\nresult = FILE_TEST('{}')\nprint,result\nwait,1\nENDREP UNTIL result EQ 1\n@{}\nexit\n".format(idl_input_file,idl_input_file)
        idl_input_file2 = idl_input_file+'2'
        with file(idl_input_file2,'w') as filehandle:
            filehandle.write(idl_text2)
        #Run IDL to set everything up
        check_file_exists(idl_input_file, max_time = 10, interval = 0.1)
        check_file_exists(idl_input_file2, max_time = 10, interval = 0.1)
        time_mod.sleep(4)
        os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))





def write_job_file(d3d_actual_dir_list, pppl_actual_dir_list, job_id, HOST):
    if HOST=='venus':
        job_template = pbs_file_venus
        run_dir_list = d3d_actual_dir_list
    else:
        job_template = pbs_file
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
        count+=1

def modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list, HOST):
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

def replace_value(input_lines, name, new_val):
    success = False
    for i,line in enumerate(input_lines):
        if (line.find(name)>=0) and (line.replace(name,'').lstrip(' ')[0] == '='):
            line_num = +i
            success = True
            break
    if success:
        new_line = "{} = {}\n".format(name, new_val)
        print '-> {} replaced with {}'.format(input_lines[line_num], new_line) 
        input_lines[line_num] = new_line
    else:
        print " WARNING line not found!!!! ", name 
    return input_lines

def modify_d3d_input(HOME, master_dict):
    shot = master_dict['settings']['new_shot']
    for time, cur_dict in master_dict['sims'].iteritems():
        print time, shot
        fname = '{}/gaprofiles/f90fidasim/{shot}/{time}/MAIN_ION330/def/d3d_input.pro'.format(HOME, shot=shot,time='{:05d}'.format(time))
        print fname
        with file(fname,'r') as filehandle: input_lines = filehandle.readlines()
        for name, new_val in cur_dict['prefida_changes'].iteritems():
            input_lines = replace_value(input_lines, name, new_val)
        with file(fname,'w') as filehandle: filehandle.writelines(input_lines)

def generate_run_FIDASIM(master_dict, time_list, shot_list, grid_settings, idl = None, HOST='venus'):
    if idl==None: idl = pidly.IDL('idl')

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

    HOME = os.environ['HOME']
    pppl_actual_dir_list, d3d_actual_dir_list = create_dirs(shot_list, time_list, diag, pppl_base_dir, d3d_base_dir, master_dict)

    make_run_idl_f90setup(shot_list, time_list, diag, d3d_base_dir, beam, comment, idl = idl)

    #Grid settings
    modify_d3d_input(HOME, master_dict,)
    make_run_idl_prefida(shot_list, time_list, diag, d3d_base_dir, beam, comment, idl = idl)
    idl.close()

    job_id = 'fidasim'
    write_job_file(d3d_actual_dir_list, pppl_actual_dir_list, job_id, HOST)
    if HOST=='portal':
        modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list, HOST)
    #What to do about the results that may already exist in the remote directory?
    if HOST=='portal': 
        for i in set(shot_list):copy_files(i)
    fida_runs_fname = HOME + '/fida_runs'
    setpoint = 13
    with file(fida_runs_fname,'w') as filehandle: filehandle.write('{}\n'.format(setpoint))
    pickle.dump(master_dict,file('/u/haskeysr/fida_sim_dict.pickle','w'))
    1/0
    batch_launch_fida(master_dict['sims'], fida_runs_fname, setpoint = setpoint, id_string = job_id)

def generate_dirs(widths, base_time, master_dict, new_shot, model_fnames, model_dir, model_shot, model_time, single = False, te_val_top_list = None, ti_val_top_list = None, ne_val_top_list = None):
    model_time_str = '{:05d}'.format(model_time)
    model_shot_str = '{}'.format(model_shot)
    new_time = base_time
    shot_list = []; time_list = []
    def copy_files(new_time, cur_dir,):
        for (tmp_prefix, tmp_fname) in model_fnames:
            new_name = tmp_fname.replace(model_shot_str,'{}'.format(new_shot)).replace(model_time_str,'{:05d}'.format(new_time))
            print tmp_prefix, tmp_fname, new_name
            os.system('cp {}/{} {}/{}'.format(model_dir, tmp_fname, cur_dir,new_name))
        #Do the same for the g-file
        os.system('cp {}/{} {}/{}'.format(model_dir, 'g{}.{:05d}'.format(model_shot,model_time), cur_dir,'g{}.{:05d}'.format(new_shot,new_time) ))
       
    if single:
        shot_list = [new_shot]; time_list = [base_time]
        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)
        os.system('mkdir -p {}'.format(cur_dir))
        copy_files(base_time, cur_dir,)
        master_dict['sims'][new_time] = {'dir_dict':{},'shot':new_shot, 'profiles':{}}
        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)
        master_dict['sims'][new_time]['dir_dict']['profiles'] = cur_dir
    else:
        if ne_val_top_list == None:
            ne_val_top_list = [4]
        if ti_val_top_list == None:
            ti_val_top_list = [3]
        if te_val_top_list == None:
            te_val_top_list = [3]

        for i, width in enumerate(widths):
            for te_val_top in te_val_top_list:
                for ti_val_top in ti_val_top_list:
                    for ne_val_top in ne_val_top_list:
                        master_dict['sims'][new_time] = {'dir_dict':{},'shot':new_shot, 'profiles':{}}
                        print new_time, width
                        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)
                        master_dict['sims'][new_time]['dir_dict']['profiles'] = cur_dir
                        os.system('mkdir -p {}'.format(cur_dir))

                        copy_files(base_time, cur_dir,)
                        #Copy model profiles but modify names
                        #idl code to modify the shot details
                        shot_list.append('{}'.format(new_shot))
                        time_list.append('{:05d}'.format(new_time))
                        new_time += 1
    return shot_list, time_list

def generate_profiles_nc(widths, base_time, master_dict, new_shot, model_fnames, model_time_str, single = False, force_ti_eq_te = True, te_val_top_list = None, ti_val_top_list = None, ne_val_top_list = None):
    if single:
        idl_string = '.compile writeg\n'
        new_time = base_time
        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)
        idl_string+=idl_ind_single.format(dir=cur_dir,shot=new_shot,time=new_time,time_str='{:05d}'.format(new_time),tag='{}'.format(new_time))
        idl_string+=idl_bulk_single
    else:
        idl_string = FIDA.idl_header
        new_time = base_time
        f = netcdf.netcdf_file('/u/haskeysr/test.nc','w', mmap=False)
        if ne_val_top_list == None:
            ne_val_top_list = [4]
        if ti_val_top_list == None:
            ti_val_top_list = [3]
        if te_val_top_list == None:
            te_val_top_list = [3]
        max_val = 1.21; npts = 121
        offset = 0.1; xsym = 0.95; hwid = 0.05
        core = [0.02]#, 0.001]
        edge = []#-0.02]
        for i, width in enumerate(widths):
            for j, te_val_top in enumerate(te_val_top_list):
                for k, ti_val_top in enumerate(ti_val_top_list):
f                    for l, ne_val_top in enumerate(ne_val_top_list):
                        #Fix the te and ti temperatures
                        cur_dir = '/u/haskeysr/gaprofiles/{}/'.format(new_shot)

                        #Force te_top = ti_top
                        if force_ti_eq_te: ti_val_top = te_val_top
                        #Create the new profiles and put in he .nc file
                        for ident, width_mult, top_val in zip(['ne', 'te', 'ti'], [1., 1., 2.], [ne_val_top, te_val_top, ti_val_top]):
                            #Generate the profiles
                            print ident, width_mult, top_val
                            x = np.linspace(0,max_val, npts)
                            y = CER.mtanh(x, top_val, offset, xsym, width*width_mult, len(core), len(edge), False, *(core + edge))
                            cur_var = '{}{}'.format(ident,new_time)
                            f.createDimension(cur_var,len(x))
                            prof = f.createVariable(cur_var,'float',(cur_var,))
                            prof[:] = +y
                            master_dict['sims'][new_time]['profiles'][ident] = {'data_y':+y,'data_x':+x,'settings':{'width':width,'xsym':xsym,'offset':offset,'pedestal_top':top_val,'core':core,'npts':npts, 'max_val':max_val,'edge':edge}}
                        #Special case to include rho
                        if new_time==base_time:
                            f.createDimension('rho',len(x))
                            rho = f.createVariable('rho','float',('rho',))
                            rho[:] = +x
                            rho.units = ''
                        idl_string+=idl_ind.format(dir=cur_dir,shot=new_shot,time=new_time,time_str='{:05d}'.format(new_time),tag='{}'.format(new_time))
                        idl_string+=idl_bulk
                        #Create the shot and time list that needs to be executed
                        #idl_strs.append('''gap_model_dir = '{}'\ngap_model_shot = '{}'\ngap_model_time = {}\n'''.format(cur_dir,new_shot,new_time))
                        new_time += 1
        time_mod.sleep(1)
        try_count = 0
        failed = True
        while failed or try_count>10:
            try:
                f.close()
                failed = False
                print 'success'
            except IndexError:
                try_count += 1
                print 'failed'
                time_mod.sleep(0.5)
    return idl_string

def get_los_data(dir='/u/haskeysr/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/', run_id = 'def', plot = False):
    #if 'inputs' not in locals():
    inputs = netcdf.netcdf_file(dir + '{}_inputs.cdf'.format(run_id),'r',mmap = False)
    los_wght = inputs.variables['los_wght'].data
    #spectra = netcdf.netcdf_file(dir + '{}_spectra.cdf'.format(run_id),'r',mmap = False)
    neutrals = netcdf.netcdf_file(dir + '{}_neutrals.cdf'.format(run_id),'r',mmap = False)
    halo_dens = neutrals.variables['halodens'].data
    #weights = netcdf.netcdf_file(dir + '{}_fida_weights.cdf'.format(run_id),'r',mmap = False)
    return halo_dens, los_wght


class interpolate_grid_profiles():
    def __init__(self, shot, time, run_id = 'def'):
        self.shot = shot
        self.time = time
        self.run_id = run_id
        HOME = os.environ['HOME']
        dir =  HOME + '/FIDASIM/RESULTS/D3D/{}/{:05d}/MAIN_ION330/'.format(shot, time) 
        self.eqdsk = OMFITtree.OMFITeqdsk(filename='/u/haskeysr/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/{}/g{}.{:05d}'.format(shot, time, run_id, shot, time))
        self.rgrid, self.zgrid = np.meshgrid(self.eqdsk['AuxQuantities']['R'],self.eqdsk['AuxQuantities']['Z'])

    def flux_values(self, r_values, z_values, use_rho = True,):
        self.r_values = r_values
        self.z_values = z_values
        if use_rho:
            self.interp_key = 'RHOpRZ'
            self.sav_key = 'RHO'
        else:
            self.interp_key = 'PSIRZ_NORM'
            self.sav_key = 'PSI'
        print 'interpolating to find rho at R,Z'
        self.flux_values = scipy.interpolate.griddata((self.rgrid.flatten(),self.zgrid.flatten()),self.eqdsk['AuxQuantities'][self.interp_key].flatten(),(self.r_values,self.z_values))

    def get_profile_RZ(self,):
        self.data_out = {}
        for plot_key in (['ti','te','ne'],ax):
            print plot_key
            fname = '/u/haskeysr/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/{}/d{}{}.{:05d}'.format(self.shot, self.time,self.run_id, plot_key, self.shot, self.time)
            dat_obj = OMFITtree.OMFITidlSav(fname)['{}_str'.format(plot_key)]
            if plot_key=='ne':
                tmp = 'DENS'
            else:
                tmp = plot_key.upper()
            print dat_obj.keys()
            print dat_obj['RHO_{}'.format(tmp)]
            dat_psi = dat_obj['{}_{}'.format(self.sav_key,tmp)]
            dat_dat = dat_obj[tmp]
            dat_interp = np.interp(flux_new, dat_psi, dat_dat,)
            self.data_out[plot_key] = data_interp.copy()
        return self.data_out


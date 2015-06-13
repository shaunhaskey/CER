import os
import time as time_mod
import subprocess as sub
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


def make_run_idl(shot_list, time_list, diag, d3d_base_dir, beam, comment):
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


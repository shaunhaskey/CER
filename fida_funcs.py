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

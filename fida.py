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

# setup_txt2 = r'''.compile /u/grierson/idlpros/add/addanon
# .compile /u/grierson/idlpros/add/addbrian
# addanon
# addbrian
# @/u/grierson/public/idl/fida/bag_fida.idl
# get_chord,155196,'t01'
# addf90fidasim
# setenv,'FIDASIM_DIR=/u/grierson/FIDASIM'
# .compile /u/grierson/idlpros/add/addf90fidasim
# addf90fidasim

# .compile /u/haskeysr/f90fidasim_setup.pro

# {f90_txt}
# .compile /u/grierson/FIDASIM/prefida
# .compile /u/grierson/FIDASIM/D3D/d3d_beams
# {d3d_input_prefida}
# exit
# '''

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

#This is a wrapper around qsub to allow remote execution
wrapper='''#!/bin/bash -l
cd {pppl_dir}
qsub {pppl_dir}jobfile.pbs

'''

#Settings
diag = 'MAIN_ION330'
beam = '330lt'
comment = 'helloworld'

#Important directories
pppl_base_dir = r'/p/fida/shaskey/RESULTS/D3D/'.replace('//','/')
d3d_base_dir = r'/u/haskeysr/FIDASIM/RESULTS/D3D/'.replace('//','/')

def create_dirs(shot_list, time_list, diag):
    time_str_list = ['{:05d}'.format(int(time)) for time in time_list]
    pppl_actual_dir_list = [r'{}/{}/{:05d}/{}/'.format(pppl_base_dir,shot, int(time), diag).replace('//','/') for shot,time in zip(shot_list, time_list)]
    d3d_actual_dir_list = [r'{}/{}/{:05d}/{}/'.format(d3d_base_dir,shot, int(time), diag).replace('//','/') for shot, time in zip(shot_list, time_list)]
    for pppl_actual_dir, d3d_actual_dir in zip(pppl_actual_dir_list, d3d_actual_dir_list):
        #print shot, time, diag, time
        tmp = 'mkdir -p {}'.format(d3d_actual_dir)
        print '==> making the directory, ', tmp
        os.system(tmp)
    return time_str_list, pppl_actual_dir_list, d3d_actual_dir_list

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

    idl_text2 = '@{}\nexit\n'.format(idl_input_file)
    idl_input_file2 = idl_input_file+'2'
    with file(idl_input_file2,'w') as filehandle:
        filehandle.write(idl_text2)
    #Run IDL to set everything up
    file_is_there = False
    while not os.path.isfile(idl_input_file):
        time_mod.sleep(0.01)
        print 'hello'

    while not os.path.isfile(idl_input_file2):
        time_mod.sleep(0.01)
        print 'hello2'
    time_mod.sleep(2)
        
    #time_mod.sleep(2)
    os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))

        
    #Run IDL to set everything up
    #os.system('idl < {} | tee {}'.format(idl_input_file, idl_output_log))

def write_job_file(d3d_actual_dir_list, pppl_actual_dir_list):
    #Write the jobfile for qsub
    for d3d_actual_dir, pppl_actual_dir in zip(d3d_actual_dir_list, pppl_actual_dir_list):
        with file('{}/jobfile.pbs'.format(d3d_actual_dir),'w') as filehandle:
            filehandle.write(pbs_file.format(PPPL_dir=pppl_actual_dir,input_file='def_inputs.dat'))
        wrapper_fname = '{}/wrapper'.format(d3d_actual_dir)
        with file(wrapper_fname,'w') as filehandle:
            filehandle.writelines(wrapper.format(pppl_dir = pppl_actual_dir))
        os.system('chmod +x {}'.format(wrapper_fname))

def modify_dat_file(d3d_actual_dir_list, pppl_actual_dir_list):
    for d3d_actual_dir, pppl_actual_dir in zip(d3d_actual_dir_list, pppl_actual_dir_list):
        with file('{}def_inputs.dat'.format(d3d_actual_dir),'r') as filehandle:
            fidasim_input = filehandle.readlines()
        for i,line in enumerate(fidasim_input):
            if line.find("result_dir")>=0:
                line_num = +i
        fidasim_input[line_num] = r"result_dir = '{}'".format(pppl_actual_dir)
        with file('{}/def_inputs.dat'.format(d3d_actual_dir),'w') as filehandle:
            filehandle.writelines(fidasim_input)

def copy_files(shot):
    #Copy everything across
    os.system('rsync -avz --progress /u/haskeysr/FIDASIM/RESULTS/D3D/{} shaskey@portal.pppl.gov:{}'.format(shot,pppl_base_dir))

def execute(pppl_actual_dir_list):
    #Submit the job on the PPPL cluster
    for pppl_actual_dir in pppl_actual_dir_list:
        cmd = r"ssh -t shaskey@portal.pppl.gov '{}/wrapper'".format(pppl_actual_dir).replace('//','/')
        print '==>',cmd
        os.system(cmd)

# def run_case(shot, time, diag):
#     time_str = '{:05d}'.format(time)
#     print shot, time, diag, time_str

#     #setup the actual directories
#     pppl_actual_dir = r'{}/{}/{:05d}/{}/'.format(pppl_base_dir,shot, time, diag).replace('//','/')
#     d3d_actual_dir = r'{}/{}/{:05d}/{}/'.format(d3d_base_dir,shot, time, diag).replace('//','/')
#     tmp = 'mkdir -p {}'.format(d3d_actual_dir)
#     print '==> making the directory, ', tmp
#     os.system(tmp)
#     idl_input_file = r'{}/{}/{:05d}/{}/idl_input_file'.format(d3d_base_dir,shot, time, diag).replace('//','/')
#     idl_output_log = r'{}/{}/{:05d}/{}/idl_output_log'.format(d3d_base_dir,shot, time, diag).replace('//','/')
#     with file(idl_input_file,'w') as filehandle:
#         tmp_txt = setup_txt.format(shot=shot,time=time_str, diag=diag, beam=beam,comment=comment)
#         print tmp_txt
#         filehandle.write(tmp_txt)
        
#     #Run IDL to set everything up
#     os.system('idl < {} | tee {}'.format(idl_input_file, idl_output_log))

#     #Write the jobfile for qsub
#     with file('{}/jobfile.pbs'.format(d3d_actual_dir),'w') as filehandle:
#         filehandle.write(pbs_file.format(PPPL_dir=pppl_actual_dir,input_file='def_inputs.dat'))

#     #Set the results dir properly
#     time_mod.sleep(1)
#     with file('{}def_inputs.dat'.format(d3d_actual_dir),'r') as filehandle:
#         fidasim_input = filehandle.readlines()
#     for i,line in enumerate(fidasim_input):
#         if line.find("result_dir")>=0:
#             line_num = +i
#     fidasim_input[line_num] = r"result_dir = '{}'".format(pppl_actual_dir)
#     with file('{}/def_inputs.dat'.format(d3d_actual_dir),'w') as filehandle:
#         filehandle.writelines(fidasim_input)

#     #Create a wrapper so that we can remotely execute it
#     wrapper_fname = '{}/wrapper'.format(d3d_actual_dir)
#     with file(wrapper_fname,'w') as filehandle:
#         filehandle.writelines(wrapper.format(pppl_dir = pppl_actual_dir))
#     os.system('chmod +x {}'.format(wrapper_fname))

#     #Copy everything across
#     os.system('rsync -avz --progress /u/haskeysr/FIDASIM/RESULTS/D3D/{} shaskey@portal.pppl.gov:{}'.format(shot,pppl_base_dir))

#     #Submit the job on the PPPL cluster
#     cmd = r"ssh -t shaskey@portal.pppl.gov '{}/wrapper'".format(pppl_actual_dir).replace('//','/')
#     os.system(cmd)

ref_shot = 158676; ref_time = 3220; new_shot = 158676
run_information = {}
n = 5
dnep_vals = np.linspace(0.2,0.6,n)
count = 200
force_ti_te_equal = True
profile_ranges = {'NE':[[0.05], np.linspace(1,8,n).tolist(), [8.0]],
                  'TE':[[0.05], np.linspace(0.5,3,n).tolist(), [3.0]],
                  'TI':[[0.05], [1.], [2.0]],
                  'TROT':[[0.], [0.], [0.]]}
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
execute(pppl_actual_dir_list)

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


import os
import numpy as np
#import matplotlib.pyplot as pt

setup_txt = r'''.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/grierson/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile writeg
.compile /u/haskeysr/gaprofiles_compute_profile_beam.pro

'''

ref_shot = 158676
ref_time = 3220
new_shot = 158676

n = 5
dnep_vals = np.linspace(0.2,0.6,n)
count = 1
for dnep in dnep_vals:
    setup_txt+="gaprofiles_modify_profiles,{},{},{},{},DNEE=0.05,DNEP={:.4f},/write\n".format(new_shot,count,ref_shot,ref_time,dnep)
    count+=1
#setup_txt+='exit\n'
print setup_txt

idl_input_file = 'test.txt'
idl_input_file2 = 'test2.txt'
idl_output_log = 'test.log'
with file(idl_input_file,'w') as filehandle:
    filehandle.write(setup_txt)
idl_text2 = '@{}\nexit\n'.format(idl_input_file)
with file(idl_input_file2,'w') as filehandle:
    filehandle.write(idl_text2)
#Run IDL to set everything up
os.system('idl < {} | tee {}'.format(idl_input_file2, idl_output_log))

'''
;========= Compile stuff
.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/grierson/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile writeg
.compile /u/haskeysr/gaprofiles_compute_profile_beam.pro
;===

;======== Modify profiles
gaprofiles_modify_profiles,158676,1,158676,3220,/write

;======== Run fidasim
f90fidasim_setup,158676,1,'MAIN_ION330','330lt',EINJ=75.0,PINJ=2.0

;======== Compile stuff for prefida
.compile /u/grierson/FIDASIM/prefida
.compile /u/grierson/FIDASIM/D3D/d3d_beams
.compile /u/haskeysr/gaprofiles/f90fidasim/158676/00001/MAIN_ION330/def/d3d_input
prefida,'input_template'

'''

'''
.compile /u/grierson/idlpros/add/addanon
.compile /u/grierson/idlpros/add/addbrian
addanon
addbrian
@/u/grierson/public/idl/fida/bag_fida.idl
addf90fidasim
get_chord,155196,'t01'
setenv,'FIDASIM_DIR=/u/grierson/FIDASIM'
.compile /u/grierson/idlpros/add/addf90fidasim
addf90fidasim
.compile writeg
.compile /u/haskeysr/gaprofiles_compute_profile_beam.pro
dir = '/u/grierson/FIDASIM/RESULTS/D3D/158676/00001/MAIN_ION330/'
f90fidasim_plot_profiles,dir,'def'
'''

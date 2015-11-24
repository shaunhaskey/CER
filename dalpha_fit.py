#!/task/imd/anaconda/bin/python
'''
dalpha_fit.py shot system
system can be either core or edge or both

'''

import pexpect, time, glob, os
import pidly
import sys, os, time
import multiprocessing
shot = int(sys.argv[1])
system = sys.argv[2].lower()
print shot, system
if system!='core' and system!='edge' and system!='both':
    raise(ValueError('Need select core or edge both system'))
idl_executable = 'idl'

def single_job(chord_list, shot, beam_list, ident):
    print chord_list, shot, beam_list
    cmds = ["@/u/haskeysr/code/idl/cc_cerview.idl\n"]
    cmds.append(".compile /u/haskeysr/bst_inst_curvefit__define.pro\n")
    cmds.append("bst_spectral_fit\n")
    cmds.append(".compile /u/haskeysr/dalpha_run_bst_spectral_fit.pro\n")
    for chord, beam in zip(chord_list,beam_list):
        cmds.append("dalpha_setup_bst_spectral_fit,{},'{}','{}'\n".format(shot, chord, beam))
        #print cmd;idl(cmd)
        cmds.append("dalpha_run_bst_spectral_fit,{},'{}','{}',/double,/guess\n".format(shot, chord, beam))
        #print cmd;idl(cmd)
    cmds.append('exit\n')
    fname = 'test_{}'.format(ident)
    with file(fname,'w') as filehandle:
        filehandle.writelines(cmds)
    time.sleep(5)
    cmd = '{} -e @{} >{}_log'.format(idl_executable, fname, fname)
    #cmd = '{} -e @{}'.format(idl_executable, fname, fname)
    os.system(cmd)

#Some level of intelligence - i.e check that
# * Wavecal has been done, and do it if it hasn't
# * Check which beams tssub are available
# * 
# system = 'edge'
# shot = 163209
# shot = 163257
# system = 'core'
# shot = 163240
# shot = 163257

print "SHOT:{}, SYSTEM:{}".format(shot,system)
chord_list_overall =  []
beam_list_overall = []
if system=='edge' or system=='both':
    chord_list = ['m17','m24','m25','m26','m27','m28','m29','m31']
    for chord in chord_list:
        chord_list_overall.append(chord)
        beam_list_overall.append('330lt')
if system=='core' or system=='both':
    chord_list = ['m02','m03','m04','m05','m06','m07','m08']
    for chord in chord_list:
        chord_list_overall.append(chord)
        beam_list_overall.append('30lt')

#chord_list = ['m24']

proc_list = []
items_per_process = 1
HOME = os.environ["HOME"]
def check_wavecal():
    wavecal_dir = "{}/cerfit/wavecal/{}/".format(HOME,shot)
    files = os.listdir(wavecal_dir)
    chords_wavecal = []
    for i in files:
        try:
            chords_wavecal.append(int(i.split('.')[1].lstrip('7').lstrip('m')))
        except:
            print('trouble with file:{}'.format(i))
    print files
    print chords_wavecal
    all_files = True

    for i in chord_list_overall:
        tmp = int(i.lstrip('m'))
        if tmp not in chords_wavecal:
            all_files = False
    if not all_files:
        raise ValueError('wavecal not there....')

check_wavecal()
def check_tssub():
    cerfit_dir = "{}/cerfit/{}/".format(HOME,shot)
    all_files = True
    for chord,beam in zip(chord_list_overall,beam_list_overall):
        tssub_file = '{}/{}/{}/tssub.dat'.format(cerfit_dir, chord, beam)
        exists = os.path.isfile(tssub_file)
        print tssub_file, exists
        if not exists:
            all_files = False
    if not all_files:
        raise ValueError('wavecal not there....')

check_tssub()

for i in range(0, len(chord_list_overall), items_per_process):
    chunk = chord_list_overall[i:i + items_per_process]
    beam = beam_list_overall[i:i + items_per_process]
    P = multiprocessing.Process(target = single_job, args=(chunk,shot, beam, i))
    proc_list.append(P)
    print chunk
    P.start()

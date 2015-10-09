import pexpect, time, glob, os
import pidly
import sys, os, time
import multiprocessing

idl_executable = 'idl'

def single_job(chord_list, shot, beam, ident):
    print chord_list, shot, beam
    cmds = ["@/u/haskeysr/code/idl/cc_cerview.idl\n"]
    cmds.append(".compile /u/haskeysr/bst_inst_curvefit__define.pro\n")
    cmds.append("bst_spectral_fit\n")
    cmds.append(".compile /u/haskeysr/dalpha_run_bst_spectral_fit.pro\n")
    for chord in chord_list:
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
system = 'edge'
shot = 163209
system = 'core'
shot = 163240
if system=='edge':
    chord_list = ['m17','m24','m25','m26','m27','m28','m29','m31']
    beam = '330lt'
else:
    chord_list = ['m02','m03','m04','m05','m06','m07','m08']
    beam = '30lt'
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

    for i in chord_list:
        tmp = int(i.lstrip('m'))
        if tmp not in chords_wavecal:
            all_files = False
    if not all_files:
        raise ValueError('wavecal not there....')

check_wavecal()
def check_tssub():
    cerfit_dir = "{}/cerfit/{}/".format(HOME,shot)
    all_files = True
    for i in chord_list:
        tssub_file = '{}/{}/{}/tssub.dat'.format(cerfit_dir, i, beam)
        exists = os.path.isfile(tssub_file)
        print tssub_file, exists
        if not exists:
            all_files = False
    if not all_files:
        raise ValueError('wavecal not there....')

check_tssub()

for i in range(0, len(chord_list), items_per_process):
    chunk = chord_list[i:i + items_per_process]
    P = multiprocessing.Process(target = single_job, args=(chunk,shot, beam, i))
    proc_list.append(P)
    print chunk
    P.start()

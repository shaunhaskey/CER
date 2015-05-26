'''
-----------------------------------
Chords and beams being fit
-----------------------------------
Beam: 30LT or 30RT Chords T01-T07, T17-T22, T33-T36
 System: CT1 @ x.0 ms
  Chords T01-T04
 System: CT2 @ x.0 ms
  Chords: T05-T07
 System: U1 @x.0 ms
  Chords: T17-T22
 System: U2 @x.0 ms
  Chords: T33-T36

Beam: 210LT or 210RT Chords T25-T32, T37-T40
 System: U4 @ x.0 ms
  Chords: T25-T32
 System: U2 @ x.0 ms
  Chords: T37-T40

Beam: 330LT or 330RT Chords T08-T16, T23-T24 from 345 R-0
Beam: 330LT Chords V1-V6, V17-V23 from 315 R-2
Beam: 330LT Chords V14-V16 from 330 R+1
Beam: 330RT Chords V7-V16, V24 from 330 R+1
 System: CT2 @ x.0 ms
  Chords: T08
 System: ET1 @ x.0 ms
  Chords: T09, T11, T13, T15
 System: ET2 @ x.0 ms
  Chords: T10, T12, T14, T16
 System: U1 @ x.0 ms
  Chords: T23, T24
 System: U3 @ x.0 ms
  Chords: V01, V02, V06, V08, V21-V24
 System: CV1 @ x.0 ms
  Chords: V17, V18, V19, V20
 System: CV2 @ x.0 ms
  Chords: V03, V04, V05, V07
 System: EV1 @ x.0 ms
  Chords: V09, V11, V13, V15 (V13 half,V15 dead)
 System: EV1 @ x.0 ms
  Chords: V10, V12, V14, V16

-----------------------------------
Description of 30LT, 30RT fit method
-----------------------------------
tssub/auto

Some examples are
cerfit_fit_chords.bash [shot] "t01 t02 t03 t04 t05 t06 t07" "tssub_30lt.dat" "beam=30lt"
cerfit_fit_chords.bash [shot] "t01 t02 t03 t04 t05 t06 t07" "tssub_30rt.dat" "beam=30rt"

cerfit_fit_chords.bash [shot] "t17 t18 t19 t20 t21 t22" "tssub_30lt.dat" "beam=30lt"
cerfit_fit_chords.bash [shot] "t17 t18 t19 t20 t21 t22" "tssub_30rt.dat" "beam=30rt"


-----------------------------------
Description of 210RT, 210LT fit method
-----------------------------------

Some examples are
cerfit_fit_chords.bash [shot] "t25 t26 t27 t28 t29 t30 t31 t32" "tssub_210lt.dat" "beam=210lt"
cerfit_fit_chords.bash [shot] "t25 t26 t27 t28 t29 t30 t31 t32" "tssub_210rt.dat" "beam=210rt"

-----------------------------------
Description of 330LT, 330RT fit method
-----------------------------------

Some examples are
From 345R-0
cerfit_fit_chords.bash [shot] "t08 t23 t09 t24 t10 t11 t12 t13 t14 t15 t16" "tssub_330lt.dat" "beam=330lt"
cerfit_fit_chords.bash [shot] "t08 t23 t09 t24 t10 t11 t12 t13 t14 t15 t16" "tssub_330rt.dat" "beam=330rt"

From 315R-2
cerfit_fit_chords.bash [shot] "v01 v02 v03 v04 v05 v06" "tssub_330lt.dat" "beam=330lt"
cerfit_fit_chords.bash [shot] "v17 v18 v19 v20 v21 v22 v23" "tssub_330lt.dat" "beam=330lt"

From 330R+1
cerfit_fit_chords.bash [shot] "v13 v14 v16" "tssub_330lt.dat" "beam=330lt"
These can get good 330LT signal but are outside FWHM
cerfit_fit_chords.bash [shot] "v07 v24 v08 v09 v10 v11 v12" "tssub_330lt.dat" "beam=330lt"
cerfit_fit_chords.bash [shot] "v07 v24 v08 v09 v10 v11 v12 v14 v16" "tssub_330rt.dat" "beam=330rt"
'''

import os, time
import multiprocessing

def make_input_file(chords, tssubs, beams, output_fname, show_plot=True):
    #show_plot_txt = 'X11' if show_plot==True else 'NONE'
    if show_plot==True: 
        show_plot_txt='X11'
    else:
        show_plot_txt = 'NONE'
    input_file = ''
    for chord, tssub, beam in zip(chords, tssubs, beams):
        #chord = 't{:02d}'.format(i)
        with file('orig_inputs/in_{}.dat'.format(chord),'r') as filehandle:
            config = filehandle.read()
        #with file('tssub_{}.dat'.format(beam),'r') as filehandle:
        with file(tssub,'r') as filehandle:
            in_lines = filehandle.read()
        in_lines = in_lines.replace('exit','')
        if chord in extra_cmds:
            extra_txt=extra_cmds[chord]
        else:
            extra_txt=''

        new_line='\n{}\n\nchord={},shot={},beam={},tplot={},{}\n\nwrite in_{}.dat\n\n{}'.format(config, chord, shot, beam, show_plot_txt, extra_txt, chord, in_lines)
        input_file += new_line
    with file(output_fname,'w') as filehandle:
        filehandle.write(input_file)

def beam_to_chord(beam, exclude_systems = None, tangential = False, vertical = False):
    if exclude_systems == None: exclude_systems = []
    chords = ''
    if ((beam=='30lt' or beam=='30rt') and tangential):
        if 'CT1' not in exclude_systems:
            chords += " t01 t02 t03 t04"
        if 'CT2' not in exclude_systems:
            chords += " t05 t06 t07"
        if 'U1' not in exclude_systems:
            chords += " t17 t18 t19 t20 t21 t22"
        if 'U2' not in exclude_systems:
            chords += " t33 t34 t35 t36"
    if ((beam=='330lt' or beam=='330rt') and tangential):
        if 'CT2' not in exclude_systems:
            chords += " t08"
        if 'ET1' not in exclude_systems:
            chords += " t09 t11 t13 t15"
        if 'ET2' not in exclude_systems:
            chords += " t10 t12 t14 t16"
        if 'U1' not in exclude_systems:
            chords += " t23 t24"
    if ((beam=='330lt') and vertical):
        if 'U3' not in exclude_systems:
            chords += " v01 v02 v06 v21 v22 v23"
        if 'CV1' not in exclude_systems:
            chords += " v17 v18 v19 v20"
        if 'CV2' not in exclude_systems:
            chords += " v03 v04 v05"
        if 'EV2' not in exclude_systems:
            chords += " v14 v16"
        #if 'EV1' not in exclude_systems:
        #    chords += " v15"
    if ((beam=='330rt') and vertical):
        if 'CV2' not in exclude_systems:
            chords += " v07"
        if 'U3' not in exclude_systems:
            chords += " v08"
        if 'EV1' not in exclude_systems:
            chords += " v09 v11 v13" # v15" is dead
        if 'EV2' not in exclude_systems:
            chords += " v10 v12 v14 v16"
        if 'U3' not in exclude_systems:
            chords += " v24"
#        chords = ' '.join([" t09 t24 t10 t11 t12 t13 t14 t15 t16",
#                           "v01 v02 v03 v04 v05 v06",
#                           "v17 v18 v19 v20 v21 v22 v23",])
##                           "v13 v14 v16",
##                           "v07 v24 v08 v09 v10 v11 v12"])
#    if beam=='330rt':
#        chords = ' '.join(["t08 t23 t09 t24 t10 t11 t12 t13 t14 t15 t16",
#                  
#         "v07 v24 v08 v09 v10 v11 v12 v14 v16"])
    if ((beam=='210lt' or beam == '210rt') and tangential):
        if 'U4' not in exclude_systems:
            chords += " t25 t16 t27 t28 t29 t30 t31 t32"
        if 'U2' not in exclude_systems:
            chords += " t33 t34 t35 t36"
        #chords = "t25 t26 t27 t28 t29 t30 t31 t32"
    chords = chords.lstrip(' ').rstrip(' ')
    if len(chords)>0:
        chords = chords.split(' ')
    else:
        chords = []
    return chords

def run_cerfit(output_fname, chords, tssubs, beams,):
    #output_fname = 'tmp2.dat'
    make_input_file(chords, tssubs, beams, output_fname)
    print '======='
    print chords
    print tssubs
    print beams
    print output_fname
    print '======='
    cmd = '''cerfit {} >{}.log'''.format(output_fname, output_fname)
    start_time = time.time()
    print cmd;os.system(cmd)
    print time.time() - start_time
    #Go back to the original directory

def make_list_of_inputs(chords, tssubs, beams, number_of_proc):
    n = len(chords)/number_of_proc
    n_float = len(chords)/float(number_of_proc)
    if float(n_float) > n:
        n = n+1
    chords_chunks = [chords[x:x+n] for x in xrange(0, len(chords), n)]
    tssubs_chunks = [tssubs[x:x+n] for x in xrange(0, len(chords), n)]
    beams_chunks = [beams[x:x+n] for x in xrange(0, len(chords), n)]
    input_data_list = []
    for i,(chords_cur, tssubs_cur, beams_cur) in enumerate(zip(chords_chunks,tssubs_chunks,beams_chunks)):
        input_data_list.append(['tmp{}.dat'.format(i), chords_cur, tssubs_cur, beams_cur])
    return input_data_list


def run_cerfit_wrapper(input_data):
    run_cerfit(*input_data)
    #output_fname, chords, tssubs, beams, output_fname = input_data


extra_cmds = {'t06':'lowerc 190\n',
              't07':'location 2=280, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't08':'location 2=260, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't09':'location 2=300, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't10':'location 2=290, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't11':'location 2=300, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't12':'location 2=290, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't13':'location 2=290, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't15':'location 2=290, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',
              't16':'location 2=290, npeak 2,temp 2=0, trange 2=0.2,freeze 2=F\n',}

shot = 160409
cerfit_dir =  '/u/haskeysr/cerfit/{}/'.format(shot)
cur_dir = os.getcwd()

try:
    with file('{}/extra_cmds.txt'.format(cerfit_dir), 'r') as filehandle:
        lines = filehandle.readlines()
except:
    lines = []
print lines
exra_cmds = {}
for i in lines:
    if len(i)>=3:
        i = i.split(':')
        print i
        lab = i[0].strip(' ')
        val = i[1].strip(' ').replace('\n','').strip(' ').strip(',')
        extra_cmds[lab] = val
for i in extra_cmds.keys():print i,extra_cmds[i]
#This copies all of the standard input files from before fy14
cmd = 'mkdir {}/orig_inputs/'.format(cerfit_dir)
print cmd;os.system(cmd)

cmd = 'cp -a /u/grierson/cerfit/setup/in/tssub/in_* {}/orig_inputs/'.format(cerfit_dir)
print cmd;os.system(cmd)

#This overwrites any of those other files with data from fy14
cmd = 'cp -a /u/grierson/cerfit/setup/in/tssub/fy14/in_* {}/orig_inputs/'.format(cerfit_dir)
print cmd;os.system(cmd)


with file('/cerbl/shotstatus/shot_{:04d}.status'.format(shot/100)) as filehandle:
    lines = filehandle.readlines()

start1 = time.time()
shot_data = {}
shot_tmp = -1
for line_id, i in enumerate(lines):
    if i[:4]=='SHOT':
        shot_tmp = int(i.split(' ')[1])
        shot_data[shot_tmp] = {}
        shot_data[shot_tmp]['lines'] = []
    elif shot_tmp==shot:
        shot_data[shot_tmp]['lines'].append(i)

lines = shot_data[shot]['lines']
spec = ''
for line_id, i in enumerate(lines):
    if i[:4]=='SPEC':
        spec = i.replace('SPEC','').replace(' ','').rstrip('\n')
        shot_data[shot][spec] = {}
    elif spec!='':
        if i[0]!='*':
            data = i.split(' ')
            shot_data[shot][spec][data[0]] = ' '.join(data[1:]).strip().rstrip('\n')
print start1 - time.time()
print shot_data[shot].keys()
b = shot_data[shot].keys()
b.remove('lines')
for i in b:
    print i, shot_data[shot][i]['TIMI'], shot_data[shot][i]['TIMI']

os.chdir(cerfit_dir)
on_beams = ['30lt', '330lt']
t_330l_ref = 't08'
v_330l_ref = 'v04'
t_30l_ref = 't01'

exclude_list = ['t14','t15','t16','v14','v16','v17','v18']
#Need an exclude list along with the times that they should be excluded to be as general as possible
#This would mean that some chords require special tssub files
#How to best implement this? Or run everything and chop at the end?

def setup_tssub_links():
    dest = '{}/{}/30lt/tssub.dat'.format(cerfit_dir,t_30l_ref)
    if os.path.isfile(dest):os.system('ln -sf {} {}/tssub_t_30lt.dat'.format(dest,cerfit_dir))
    dest = '{}/{}/330lt/tssub.dat'.format(cerfit_dir,t_330l_ref)
    if os.path.isfile(dest):os.system('ln -sf {} {}/tssub_t_330lt.dat'.format(dest,cerfit_dir))
    dest = '{}/{}/330lt/tssub.dat'.format(cerfit_dir,v_330l_ref)
    if os.path.isfile(dest):os.system('ln -sf {} {}/tssub_v_330lt.dat'.format(dest,cerfit_dir))
    #os.system('ln -sf {}/{}/330lt/tssub.dat {}/tssub_t_330lt.dat'.format(cerfit_dir,t_330l_ref,cerfit_dir))
    #os.system('ln -sf {}/{}/330lt/tssub.dat {}/tssub_v_330lt.dat'.format(cerfit_dir,v_330l_ref,cerfit_dir))


setup_tssub_links()

test_run = False
number_of_proc = 5
if test_run: number_of_proc = 1
tangential = True
vertical = False
if tangential: os.system('rm 160409.8t* d160409_tang*.nc')
if vertical: os.system('rm 160409.8v* d160409_vert*.nc')
exclude_systems = ['U1', 'U2']

if test_run:
    chords = ['t01', 't08', 'v04']
    beams = ['30lt', '330lt', '330lt']
    tssubs = ['tssub_{}_{}.dat'.format(chord[0], beam) for chord,beam in zip(chords,beams)]
    #tssubs = ['tssub_30lt.dat','tssub_330lt.dat']
else:
    chords = [];tssubs = [];beams = []
    for on_beam in on_beams:
        chords_tmp = beam_to_chord(on_beam, exclude_systems = exclude_systems, tangential = tangential, vertical = vertical)
        for i in exclude_list: 
            try:
                chords_tmp.remove(i)
            except ValueError:
                pass
                
        print chords_tmp
        if len(chords_tmp)>0:
            tssubs_tmp = ['tssub_{}_{}.dat'.format(i[0], on_beam) for i in chords_tmp]
            beams_tmp = [on_beam for i in chords_tmp]
            chords.extend(chords_tmp)
            tssubs.extend(tssubs_tmp)
            beams.extend(beams_tmp)

input_data_list = make_list_of_inputs(chords, tssubs, beams, number_of_proc)

if number_of_proc>1:
    p = multiprocessing.Pool(number_of_proc)
    p.map(run_cerfit_wrapper, input_data_list)
else:
    map(run_cerfit_wrapper, input_data_list)

#output_fname = 'tmp2.dat'
#input_data = [output_fname, chords, tssubs, beams]
#run_cerfit_wrapper(input_data)

##Go back to the original directory
os.chdir(cur_dir)

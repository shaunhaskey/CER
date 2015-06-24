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

import os, time, copy
import multiprocessing
import cer_funcs as CER
#read_subtract_file(fname)

def assemble_subtracts(sub_list, fname = False):
    out_txt = ''.join(['time={}\ntssub={}\ngo\n\n'.format(i[0],i[1]) for i in sub_list])
    if fname != False:
        with file(fname,'w') as filehandle:
            filehandle.write(out_txt)
    return out_txt

def remove_times(sub_data, chord):
    new_sub_data = []
    kill_times_int = [int(float(i)) for i in  kill_times[chord]]
    print '##############'
    print kill_times[chord], lower_time, upper_time
    for i in range(len(sub_data)):
        cur_time = int(float(sub_data[i][0]))
        if (cur_time not in kill_times_int) and (cur_time<upper_time) and (cur_time>lower_time):# kill_times[chord]: 
            new_sub_data.append(sub_data[i])
        else:
            print 'removing', chord, sub_data[i]
    print '##############'
    return new_sub_data

def mod_times(sub_data, chord):
    new_sub_data = []
    modify_times_int = [int(float(i)) for i in  modify_times[chord]['times']]
    for i in range(len(sub_data)):
        cur_time = int(float(sub_data[i][0]))
        if (cur_time not in modify_times_int):
            #if sub_data[i][0] not in modify_times[chord]['times']: 
            new_sub_data.append(sub_data[i])
        else:
            index = modify_times_int.index(cur_time)
            new_sub_data.append([sub_data[i][0],modify_times[chord]['subs'][index]])
            print 'modifying chord {}, old: {}, new: {}'.format(chord, sub_data[i],new_sub_data[-1])
    return new_sub_data

def make_input_file(chords, tssubs, beams, output_fname, show_plot=True):
    #show_plot_txt = 'X11' if show_plot==True else 'NONE'
    if show_plot==True: 
        show_plot_txt='X11'
    else:
        show_plot_txt = 'NONE'
    input_file = ''
    for chord, tssub, beam in zip(chords, tssubs, beams):
        with file('orig_inputs/in_{}.dat'.format(chord),'r') as filehandle:
            config = filehandle.read()
        sub_data = CER.read_subtract_file(tssub)
        tmp_name = '{}_{}.dat'.format(tssub.split('.')[0], chord)
        sub_data = remove_times(sub_data, chord)
        sub_data = mod_times(sub_data, chord)
        sub_lines = assemble_subtracts(sub_data, fname = tmp_name)
        extra_txt=extra_cmds[chord] if chord in extra_cmds else ''
        new_line='\n{}\n\nchord={},shot={},beam={},tplot={}\n{}\n\nwrite in_{}.dat\n\n{}'.format(config, chord, shot, beam, show_plot_txt, extra_txt, chord, sub_lines)
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
    print lower_time, upper_time
    print chords
    print tssubs
    print beams
    print output_fname
    print '======='
    cmd = '''cerfit {} >{}.log'''.format(output_fname, output_fname)
    start_time = time.time()

    while not os.path.isfile(output_fname):
        time.sleep(1)
    time.sleep(1.5)
    print cmd;os.system(cmd)
    print time.time() - start_time, cmd
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

shot = 160409
shot = 160414
cerfit_dir =  '/u/haskeysr/cerfit/{}/'.format(shot)
cur_dir = os.getcwd()

def read_in_npeaks(chords):
    npeaks = {}
    for chord in chords:
        fname = 'orig_inputs/in_{}.dat'.format(chord)
        with file(fname,'r') as filehandle: lines = filehandle.readlines()
        for i in lines:
            if i.find('NPEAK')>=0:
                npeaks[chord] = int(i.split('=')[1].rstrip('\n').strip(' '))
    return npeaks

# def parse_line(line):
#     line = line.split(':')
#     print line
#     if i[0]=='cold_line':
#         chrd = i[1].strip(' ')
#         loc = i[2].strip(' ')
#         npeaks[chrd]= npeaks[chrd]+1
#         cmd = 'location {npeak}={loc},npeak {npeak}, temp {npeak}=0, trange {npeak} = 0.2, freeze {npeak}=F'.format(npeak=npeaks[chrd], loc=loc)
#     if i[0]=='command':
#         chrd = i[1].strip(' ')
#         cmd = i[2]
#     if i[0]=='kill':
#         pass
#     if i[0]=='modify':
#         pass
#     #lab = i[0].strip(' ')
#     #val = i[1].strip(' ').replace('\n','').strip(' ').strip(',')
#     if lab not in extra_cmds.keys():
#         extra_cmds[lab] = val
#     else:
#         extra_cmds[lab] += ',{}'.format(cmd)
#         print extra_cmds[lab]

# def read_extra_cmds():
#     try:
#         with file('{}/extra_cmds2.txt'.format(cerfit_dir), 'r') as filehandle:
#             lines = filehandle.readlines()
#     except:
#         lines = []
#     print lines
#     extra_cmds = {}
#     for i in lines:
#         if len(i)>=3:
#             parse_line(i)
#             i = i.split(':')
#             print i
#             lab = i[0].strip(' ')
#             val = i[1].strip(' ').replace('\n','').strip(' ').strip(',')
#             extra_cmds[lab] = val
#     for i in extra_cmds.keys():print i,extra_cmds[i]
#     kill_times = {}
#     modify_times = {}
#     for i in chords: 
#         kill_times[i] = []
#         modify_times[i] = {'times':[],'subs':[]}

#     return extra_cmds, kill_times, modify_times

def read_extra_cmds():
    try:
        with file('{}/extra_cmds2.txt'.format(cerfit_dir), 'r') as filehandle:
            lines_orig = filehandle.readlines()
    except:
        lines_orig = []
    lines = copy.deepcopy(lines_orig)
    print lines
    extra_cmds = {}
    kill_times = {}
    modify_times = {}
    #for i in chords: 
    for prefix in ['t','v']:
        for i in range(50): 
            kill_times['{}{:02d}'.format(prefix,i)] = []
            modify_times['{}{:02d}'.format(prefix,i)] = {'times':[],'subs':[]}

    lower_time = 0
    upper_time = 10000
    all_chans = ','.join(chords).rstrip(',')
    print all_chans
    for i in range(len(lines)):
        if lines[i].find('all')>=0:
            lines[i] = lines[i].replace('all',all_chans)
            print lines[i]
    for i in lines:
        if len(i)>=3:
            i = i.split(':')
            i[0] = i[0].strip(' ')
            print i
            if i[0]=='cold_line':
                chrd = i[1].strip(' ')
                loc = i[2].rstrip('\n').strip(' ')
                if chrd in npeaks.keys():
                    npeaks[chrd]= npeaks[chrd]+1
                else:
                    npeaks[chrd] = 2
                    
                print 'cold_line', i, chrd, loc, npeaks
                cmd = 'location {npeak}={loc}, temp {npeak}=0\ntrange {npeak} = 0.2, freeze {npeak}=F\n'.format(npeak=npeaks[chrd], loc=loc)
                if chrd not in extra_cmds.keys():
                    extra_cmds[chrd] = cmd
                else:
                    extra_cmds[chrd] += '{}'.format(cmd)
                rerun_whole_shot[chrd] = True
            if i[0]=='command':
                chrds = i[1].strip(' ').split(',')
                #chrd = i[1].strip(' ')
                cmd = i[2].rstrip('\n')
                cmd += '\n'
                for chrd in chrds:
                    if chrd not in extra_cmds.keys():
                        extra_cmds[chrd] = cmd
                    else:
                        extra_cmds[chrd] += '{}'.format(cmd)
                    rerun_whole_shot[chrd] = True
            if i[0]=='kill':
                chrds = i[1].strip(' ').split(',')
                kill_time = i[2].rstrip('\n').strip(' ')
                for tmp_chrd in chrds:
                    if tmp_chrd in kill_times.keys():
                        kill_times[tmp_chrd].append(kill_time)
                    else:
                        kill_times[tmp_chrd] = [kill_time]
            if i[0]=='modify':
                chrds = i[1].strip(' ').split(',')
                modify_time = i[2].rstrip('\n').strip(' ')
                new_sub = i[3].rstrip('\n').strip(' ')
                for tmp_chrd in chrds:
                    if tmp_chrd in modify_times.keys():
                        modify_times[tmp_chrd]['times'].append(modify_time)
                        modify_times[tmp_chrd]['subs'].append(new_sub)
                    else:
                        modify_times[tmp_chrd]['times'] = [modify_time]
                        modify_times[tmp_chrd]['subs'] = [new_sub]
            if i[0]=='lower_time':
                lower_time = int(float(i[1].strip(' ').rstrip('\n')))
            if i[0]=='upper_time':
                upper_time = int(float(i[1].strip(' ').rstrip('\n')))
            #lab = i[0].strip(' ')
            #val = i[1].strip(' ').replace('\n','').strip(' ').strip(',')
            #print extra_cmds[chrd]
            #parse_line(i)
            #i = i.split(':')
            #print i
            #lab = i[0].strip(' ')
            #val = i[1].strip(' ').replace('\n','').strip(' ').strip(',')
            #extra_cmds[lab] = val
    for i in extra_cmds.keys():extra_cmds[i] += 'npeak {}\n'.format(npeaks[i])
    for i in extra_cmds.keys():print i,extra_cmds[i]
    return extra_cmds, kill_times, modify_times, lower_time, upper_time, lines_orig


def read_remove_list():
    try:
        with file('{}/extra_cmds2.txt'.format(cerfit_dir), 'r') as filehandle:
            lines = filehandle.readlines()
    except:
        lines = []

    #Strip out any comments in the commands
    lines_new = []
    for i in lines:
        if i.strip(' ').find('#')==0:
            pass
        elif i.find('#')>0:
            lines_new.append(i[:i.find('#')])
        else:
            lines_new.append(i)
    print lines
    print lines_new
    lines = lines_new
    print lines
    remove_list = []
    for i in lines:
        if i.find('remove')>=0:
            i = i.split(':')
            remove_list.extend(i[1].strip(' ').rstrip('\n').split(','))
    return remove_list

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


def read_references():
    with file('references.txt','r') as filehandle: lines = filehandle.readlines()
    tssub_refs = {}
    print lines
    for i in lines:
        print i,i.find(',')
        if i.find(',')>=0:
            print i
            dat = i.split(',')
            print dat
            name = dat[0].strip(' ')
            val = dat[1].rstrip('\n').strip(' ')
            print name, val
            tssub_refs[name] = val
    return tssub_refs
os.chdir(cerfit_dir)
tssub_refs = read_references()
on_beams = ['30lt', '330lt']
#t_330l_ref = 't08'
#v_330l_ref = 'v04' 160409
#v_330l_ref = 't08'
#t_30l_ref = 't01'


exclude_list = read_remove_list()
#exclude_list = ['t14','t15','t16','v14','v16','v17','v18']
#Need an exclude list along with the times that they should be excluded to be as general as possible
#This would mean that some chords require special tssub files
#How to best implement this? Or run everything and chop at the end?

def setup_tssub_links(subtract_name):
    #subtract_name = 'timesub' if timesumb else 'tssub'
    dest = '{}/{}/30lt/{}.dat'.format(cerfit_dir,tssub_refs['t_30l_ref'],subtract_name)
    if os.path.isfile(dest):os.system('ln -sf {} {}/{}_t_30lt.dat'.format(dest,cerfit_dir, subtract_name))
    dest = '{}/{}/330lt/{}.dat'.format(cerfit_dir,tssub_refs['t_330l_ref'], subtract_name)
    if os.path.isfile(dest):os.system('ln -sf {} {}/{}_t_330lt.dat'.format(dest,cerfit_dir, subtract_name))
    dest = '{}/{}/330lt/{}.dat'.format(cerfit_dir,tssub_refs['v_330l_ref'], subtract_name)
    if os.path.isfile(dest):os.system('ln -sf {} {}/{}_v_330lt.dat'.format(dest,cerfit_dir, subtract_name))
    #os.system('ln -sf {}/{}/330lt/tssub.dat {}/tssub_t_330lt.dat'.format(cerfit_dir,t_330l_ref,cerfit_dir))
    #os.system('ln -sf {}/{}/330lt/tssub.dat {}/tssub_v_330lt.dat'.format(cerfit_dir,v_330l_ref,cerfit_dir))

subtract_name = 'timesub'

setup_tssub_links(subtract_name)

test_run = False
number_of_proc = 4
if test_run: number_of_proc = 1
tangential = True
vertical = True
if tangential: os.system('rm {shot}.8t* d{shot}_tang*.nc'.format(shot = shot))
if vertical: os.system('rm {shot}.8v* d{shot}_vert*.nc'.format(shot=shot))
#exclude_systems = ['U1', 'U2']
exclude_systems = []

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
            tssubs_tmp = ['{}_{}_{}.dat'.format(subtract_name,i[0], on_beam) for i in chords_tmp]
            beams_tmp = [on_beam for i in chords_tmp]
            chords.extend(chords_tmp)
            tssubs.extend(tssubs_tmp)
            beams.extend(beams_tmp)

npeaks =  read_in_npeaks(chords)
rerun_whole_shot = {i:False for i in chords}

extra_cmds, kill_times, modify_times, lower_time, upper_time, orig_lines = read_extra_cmds()
print "Lower time: {}, upper time: {}".format(lower_time, upper_time)
#1/0

input_data_list = make_list_of_inputs(chords, tssubs, beams, number_of_proc)

number_of_proc = 5
if number_of_proc>1:
    p = multiprocessing.Pool(number_of_proc)
    p.map(run_cerfit_wrapper, input_data_list)
else:
    map(run_cerfit_wrapper, input_data_list)

with file('{}/last_used_commands.txt'.format(cerfit_dir), 'w') as filehandle:
    lines_orig = filehandle.writelines(orig_lines)

#output_fname = 'tmp2.dat'
#input_data = [output_fname, chords, tssubs, beams]
#run_cerfit_wrapper(input_data)

##Go back to the original directory
os.chdir(cur_dir)

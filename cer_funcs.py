import os
import numpy as np

def split_timing_string(stri):
    stri.strip(' ').strip('\n')
    start_time = float(stri.split(':')[0])
    n_slices = int(stri.split(':')[1].split('@')[0])
    t_int = float(stri.split(':')[1].split('@')[1].split(r'/')[0])
    start_pts = np.arange(n_slices)*t_int + start_time
    end_pts  = start_pts + t_int
    return start_time, n_slices, t_int, start_pts, end_pts

def get_nbi_data(shot, channel,):
    os.environ['VPN_ACTIVE']='yes'
    import data as MDSdata
    dat = MDSdata.Data('nbvolt_{}'.format(channel),shot)
    return dat.x[0], dat.y

def get_fscope(shot, channel = 'fs04f'):
    os.environ['VPN_ACTIVE']='yes'
    import data as MDSdata
    dat = MDSdata.Data(channel,shot)
    return dat.x[0], dat.y
    # if plot:

def read_shot_status(shot):
    with file('/cerbl/shotstatus/shot_{:04d}.status'.format(shot/100)) as filehandle:lines = filehandle.readlines()
    shot_data = {}; shot_tmp = -1
    #Scan through the file looking for all the shots and keep the relevant part of text for each one
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
    print shot_data[shot].keys()
    b = shot_data[shot].keys()
    del shot_data[shot]['lines']#b.remove('lines')
    return shot_data[shot]

def read_subtract_file(fname, tssub=False):
    '''If tssub is true then it expects a tssub file else a time file
    
    '''
    with file(fname,'r') as filehandle: in_lines = filehandle.readlines()
    data = []
    txt = 'ts=' if tssub else 'time='
    for i in range(len(in_lines)):
        cur_line = in_lines[i]
        if cur_line.find(txt)==0:
            cur_line_plus = in_lines[i+1]
            data.append([cur_line.split('=')[-1].rstrip('\n'), cur_line_plus.split('=')[-1].rstrip('\n')])
    return data


def mtanh(*args):
    '''Copied from /scp:venus:/u/grierson/idlpros/utilities/mpfit_mtanh.pro

    SRH : 02June2015
    '''

    x = args[0] 
    top = args[1]
    offset = args[2] 
    xsym = args[3] 
    hwid = args[4] 
    ncore = args[5] 
    nedge = args[6] 
    plot = args[7] 
    core = args[8:8+ncore] 
    edge = args[8+ncore:8+ncore+nedge] 


    print core, edge
    #x = np.arange(121)/100.
    a = (top-offset)/2.0
    b = (3*offset+top)/2.0
    z = (xsym-x)/hwid
    if core==None: core = []
    if edge==None: edge = []
    if len(core)>0:
        poly_core = 1.0
        for i in range(len(core)):
            poly_core += core[i] * z**(i+1)
        num_core = poly_core*np.exp(z)
    else:
        num_core = np.exp(z)
    #;; Edge with a negative slope
    if len(edge)>0:
        poly_edge = 1.0
        for i in range(len(edge)):
            poly_edge += edge[i] * z**(i+1)
        num_edge = poly_edge*np.exp(-z)
    else:
        num_edge = np.exp(-z)
    num = num_core - num_edge
    den = np.exp(z) + np.exp(-z)
    y = a * num/den + b
    print plot
    if plot:
        import matplotlib.pyplot as pt
        print 'hello'
        fig, ax = pt.subplots()
        ax.plot(x,y)
        fig.canvas.draw();fig.show()
    return y



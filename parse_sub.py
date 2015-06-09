import numpy as np
import matplotlib.pyplot as pt

def read_subtract_file(fname):
    with file(fname,'r') as filehandle: in_lines = filehandle.readlines()
    data = []
    for i in range(len(in_lines)):
        cur_line = in_lines[i]
        if cur_line.find('time')==0:
            cur_line_plus = in_lines[i+1]
            data.append([cur_line.split('=')[-1].rstrip('\n'), cur_line_plus.split('=')[-1].rstrip('\n')])
    return data

def assemble_subtracts(sub_list):
    out_txt = ''
    for i in sub_list:
        out_txt+='time={}\ntssub={}\ngo\n\n'.format(i[0],i[1])
    return out_txt

shot = 160409
cerfit_dir =  '/u/haskeysr/cerfit/{}/'.format(shot)
fname = cerfit_dir + 'timesub_t_30lt_t03.dat'
data = read_subtract_file(fname)
out_txt = assemble_subtracts(data)

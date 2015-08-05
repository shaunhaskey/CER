import numpy as np
import pidly, os
import MDSplus as MDS
idl_setup = '''addanon
@/u/kaplan/cerview/cerview.idl
get_chord,149661,'t01'
.compile /u/grierson/idlpros/add/addbrian
.compile /u/grierson/idlpros/add/addnovi
addbrian
addnovi
.compile bst_dispersion
common bst_chord_param, chord_param
'''

def get_mainion_geom(shot=163110):
    os.environ['ions_path'] = 'atlas.gat.com::'
    T = MDS.Tree('ions',shot)
    data = {}
    beam_order = T.getNode('.cer.calibration:beam_order').data()
    for i in range(1,33):
        lens_phi = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_phi'.format(i)).data()
        lens_R = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_R'.format(i)).data()
        lens_Z = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_Z'.format(i)).data()
        plasma_phi = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_phi'.format(i)).data()
        plasma_R = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_R'.format(i)).data()
        plasma_Z = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_Z'.format(i)).data()
        lens = np.array([lens_R * np.sin(np.deg2rad(lens_phi)), lens_R * np.cos(np.deg2rad(lens_phi)), lens_Z])
        data['m{:02d}'.format(i)] = {'lens_xyz': lens, 'lens_rzphi': np.array([lens_R,lens_Z, lens_phi])}
        print 'hello'
        #print i, n1.data(), n2.data(), n3.data()
        for phi, R, Z, beam in zip(plasma_phi, plasma_R, plasma_Z, beam_order):
            tmp = np.array([R * np.sin(np.deg2rad(phi)), R * np.cos(np.deg2rad(phi)), Z])
            data['m{:02d}'.format(i)]['plasma_{}_xyz'.format(beam.lower())] = +tmp
            data['m{:02d}'.format(i)]['plasma_{}_rzphi'.format(beam.lower())] = np.array([R,Z,phi])
        #print ' ', i, n1.data(), n2.data(), n3.data()
    return data
    
a = get_mainion_geom(shot=163110)

 # 'm28': {'lens_rzphi': array([  2.45620370e+00,   6.82199979e-03,   3.46330505e+02], dtype=float32),
 #  'lens_xyz': array([-0.58045208,  2.38663197,  0.006822  ], dtype=float32),
 #  'plasma_150lt_rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_150lt_xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_150rt_rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_150rt_xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_210lt_rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_210lt_xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_210rt_rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_210rt_xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_30lt _rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_30lt _xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_30rt _rzphi': array([   3.        ,    9.99989986,  999.98999023], dtype=float32),
 #  'plasma_30rt _xyz': array([-2.9545145 ,  0.52042705,  9.99989986], dtype=float32),
 #  'plasma_330lt_rzphi': array([  2.27322340e+00,  -1.24828452e-02,   3.23681702e+02], dtype=float32),
 #  'plasma_330lt_xyz': array([-1.34636307,  1.83162534, -0.01248285], dtype=float32),
330lt 28 [ -58.0452  238.6632    0.6822] [-134.63624032  183.16258173   -1.24841603]

 #  'plasma_330rt_rzphi': array([  2.27415943e+00,  -1.07990215e-02,   3.25760895e+02], dtype=float32),
 #  'plasma_330rt_xyz': array([-1.27955043,  1.88004041, -0.01079902], dtype=float32)},
330rt 28 [ -58.0452  238.6632    0.6822] [-127.95487098  188.00414189   -1.08      ]

1/0
idl = pidly.IDL('/usr/local/bin/idl')
for i in idl_setup.split('\n'):
    print i
    idl(i)
chord_geom = {}
chord_geom['330lt']={'chords':[1,17,24,25,26,27,28,29,31],'lens':[],'location':[]}
chord_geom['330rt']={'chords':[1,2,17,24,25,26,27,28,29,31],'lens':[],'location':[]}
chord_geom['30lt']={'chords':[1,2,3,4,5,6,7,8],'lens':[],'location':[]}
chord_geom['30rt']={'chords':[1,2,3,4,5,6,7,8],'lens':[],'location':[]}
chord_geom['210lt']={'chords':[9,10,11,12,13,14,15,16],'lens':[],'location':[]}
chord_geom['210rt']={'chords':[9,10,11,12,13,14,15,16],'lens':[],'location':[]}
chord_geom['210lt']={'chords':[10,11,12,13,],'lens':[],'location':[]}
chord_geom['210rt']={'chords':[10,11,12,13,],'lens':[],'location':[]}

for beam in chord_geom.keys():
    print 'beam', beam
    for i in chord_geom[beam]['chords']:
        idl("bst_chord_param,163257,'m{}','{}'".format(i, beam))
        idl('lens=chord_param.geometry.lens')
        idl('location=chord_param.geometry.location')
        chord_geom[beam]['location'].append(idl.locatio)
        chord_geom[beam]['lens'].append(idl.lens)
        print beam, i, idl.lens, idl.location


idl('x=chord_param.geometry.lens')


#idl = pidly.idl()
import matplotlib.pyplot as pt
fig, ax = pt.subplots()
for beam in chord_geom.keys():
    print 'beam', beam
    for i in range(len(chord_geom[beam]['chords'])):
        ax.plot([chord_geom[beam]['location'][i][0], chord_geom[beam]['lens'][i][0]], [chord_geom[beam]['location'][i][1], chord_geom[beam]['lens'][i][1]],'-o')
fig.canvas.draw();fig.show()


import MDSplus as MDS
import os
os.environ['ions_path'] = 'atlas.gat.com::'
shot = 163110
T = MDS.Tree('ions',shot)
for i in range(1,33):
    n1 = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_phi'.format(i))
    n2 = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_R'.format(i))
    n3 = T.getNode('.cermain.calibration.tangential.channel{:02d}:lens_Z'.format(i))
    print i, n1.data(), n2.data(), n3.data()
    n1 = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_phi'.format(i))
    n2 = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_R'.format(i))
    n3 = T.getNode('.cermain.calibration.tangential.channel{:02d}:plasma_Z'.format(i))
    print ' ', i, n1.data(), n2.data(), n3.data()

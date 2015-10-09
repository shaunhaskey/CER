from OMFITtree import OMFITmdsValue
import numpy as np
shot = 163257
def get_nbi(shot, plot=False, beams=None):
    out_dict = {'names':[],'name_map':[],'volts':[],'power':[], 'shot':shot}
    print ''
    for deg, deg_long in zip(['30','15','21','33'], ['30','150','210','330']):
        for orient in ['L','R']:
            out_dict['names'].append('nbvac{}{}t'.format(deg,orient.lower()))
            out_dict['name_map'].append('{}{}t'.format(deg_long,orient.lower()))
            b = OMFITmdsValue('atlas',treename=None,shot=shot,TDI='NBVAC{}{}T'.format(deg,orient))
            c = OMFITmdsValue('atlas',treename=None,shot=shot,TDI='NBVOLT_{}{}'.format(deg,orient))
            d = OMFITmdsValue('atlas',treename='d3d',shot=shot,TDI='.nb.nb{}{}:nbvac_scalar'.format(deg,orient.lower()))

            out_dict['volts'].append(np.max(d.data()[0]))
            e = OMFITmdsValue('atlas',treename='d3d',shot=shot,TDI='.nb.nb{}{}:pinj_scalar'.format(deg,orient.lower()))
            out_dict['power'].append(e.data()[0])
            print('{}{}: max NBVAC {:.4f}, max NBVOLT {:.4f}'.format(deg,orient, np.max(b.data()),np.max(c.data())))
            print d.data()[0], e.data()[0]
    return out_dict
out_dict = get_nbi(shot)
c = OMFITmdsValue('atlas',treename='d3d',shot=shot,TDI='.nb.nb15l:nbvac_scalar')

# Can get the chord geometry as follows:
# common bst_chord_param, chord_param1
# @/u/haskeysr/code/idl/cc_cerview.idl
# #dalpha_cerfit_wavecal2,163019,'m02'
# #[1,2,3,4,5,6,7,8,17,24,25,26,27,28,29,30,31]
# bst_chord_param,158676,'m24','330lt'
# bst_chord_param,163257,'m24','330lt'
# bst_chord_param,158676,'m20','330lt'
# print,chord_param1.geometry.location
# print,chord_param1.geometry.location

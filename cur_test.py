import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pt
#import OMFITtree
import pidly
import os
import time as time_mod
import subprocess as sub
import cPickle as pickle
from scipy.io import netcdf
import numpy as np
import cer_funcs as CER
import scipy
import shutil
import itertools as iter
import cer_funcs as CER
import os
import pyMARS.generic_funcs as gen
import fida_funcs as FIDA




def get_los_data(dir='/u/haskeysr/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/', run_id = 'def', plot = False):
    #if 'inputs' not in locals():
    inputs = netcdf.netcdf_file(dir + '{}_inputs.cdf'.format(run_id),'r',mmap = False)
    los_wght = inputs.variables['los_wght'].data
    spectra = netcdf.netcdf_file(dir + '{}_spectra.cdf'.format(run_id),'r',mmap = False)
    neutrals = netcdf.netcdf_file(dir + '{}_neutrals.cdf'.format(run_id),'r',mmap = False)
    halo_dens = neutrals.variables['halodens'].data
    #weights = netcdf.netcdf_file(dir + '{}_fida_weights.cdf'.format(run_id),'r',mmap = False)

    origin = inputs.variables['origin'].data
    # x_grid = inputs.variables['x_grid'].data + origin[0]
    # y_grid = inputs.variables['y_grid'].data + origin[1]
    # z_grid = inputs.variables['z_grid'].data + origin[2]

    x_grid = inputs.variables['x_grid'].data
    y_grid = inputs.variables['y_grid'].data
    z_grid = inputs.variables['z_grid'].data

    r_grid = inputs.variables['r_grid'].data
    u_grid = inputs.variables['u_grid'].data
    v_grid = inputs.variables['v_grid'].data


    return halo_dens, los_wght, x_grid, y_grid, z_grid, r_grid, u_grid, v_grid, inputs, origin, neutrals, spectra



#b = FIDA.interpolate_grid_profiles(155196,500)
HOME = os.environ['HOME']
shot = 155196
time = 500
test = FIDA.fidasim_results(HOME+ '/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/',shot, time)
#halo_dens, los_wght, x_grid, y_grid, z_grid, u_grid,v_grid, r_grid,inputs,origin,neutrals, spectra = test.get_los_data()
#test.transform_data()
test.fit_gaussians(decimation = 1)
test.plot_LOS_densities()
test.plot_profiles_responses(plot_items = ['ti', 'TOR_ROT'])

#halo_dens, los_wght, x_grid, y_grid, z_grid, u_grid,v_grid, r_grid,inputs,origin,neutrals, spectra  = get_los_data(dir=HOME + '/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/', run_id = 'def', plot = False)

#halo_dens, los_wght, x_grid, y_grid, z_grid, u_grid,v_grid, r_grid,inputs,origin,neutrals, spectra  = get_los_data(dir=HOME + '/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/', run_id = 'def', plot = False)

# alpha = inputs.variables['alpha'].data

# def transform(x_grid, y_grid, z_grid, origin, alpha, unit_conversion=1./100):
#     phi_tmp = np.arctan2(y_grid, x_grid)
#     r_tmp = np.sqrt(x_grid**2 + y_grid**2)

#     #Values in m
#     x_grid2 = (r_tmp*np.cos(phi_tmp + alpha) + origin[0])*unit_conversion
#     y_grid2 = (r_tmp*np.sin(phi_tmp + alpha) + origin[1])*unit_conversion
#     z_grid2 = (z_grid + origin[2])*unit_conversion
#     return x_grid2, y_grid2,z_grid2


# x_grid2, y_grid2,z_grid2 =  transform(x_grid, y_grid, z_grid, origin, alpha)

# b.get_flux_values(np.sqrt(x_grid2**2 + y_grid2**2),z_grid)
# c, d = b.get_profile_RZ()
# xlens = inputs.variables['xlens'].data
# ylens = inputs.variables['ylens'].data
# zlens = inputs.variables['zlens'].data

# xlos = inputs.variables['xlos'].data
# ylos = inputs.variables['ylos'].data
# zlos = inputs.variables['zlos'].data

# xlens2, ylens2, zlens2 = transform(xlens, ylens, zlens, origin, alpha)
# xlos2, ylos2, zlos2 =  transform(xlos, ylos, zlos, origin, alpha)

# n_chrds = 16

# xlos_chords,ylos_chords,zlos_chords = [np.linspace(-133.64747, -134.95745, n_chrds)/100,
#                                        np.linspace(172.84416, 186.51464 , n_chrds)/100,
#                                        np.linspace(-1.4165587, -1.0956602,  n_chrds)/100]

# xlens_chords,ylens_chords,zlens_chords = [np.linspace(-58.045200, -58.045200, n_chrds)/100,
#                                           np.linspace(238.66320, 238.66320, n_chrds)/100,
#                                           np.linspace(0.68220000, 0.68220000, n_chrds)/100]


# xlos_chords2 = np.interp(xlos_chords,xlos2,xlos)
# ylos_chords2 = np.interp(ylos_chords,ylos2,ylos)
# zlos_chords2 = np.interp(zlos_chords,zlos2,zlos)

# fig, ax = pt.subplots(ncols = 4,nrows = 6, sharex = True)
# ax6 = ax[5,:]
# ax5 = ax[4,:]
# ax4 = ax[3,:]
# ax3 = ax[2,:]
# ax2 = ax[1,:]
# ax1 = ax[0,:]

# clr_cycle = gen.new_color_cycle(0,len(xlens))
# plot_quants = ['rho_grid','te','dene','ti']
# peaks = {i:{'x':[],'y':[]} for i in plot_quants}


# halo = +spectra.variables['halo'].data
# wave = +spectra.variables['lambda'].data * 10 #Angstroms
# ti_list = []
# vel_list = []

#ti_list = test.ti_list
#vel_list = test.vel_list

#1/0
# mean_L_list = test.mean_L_list
# mean_dens_list = test.mean_dens_list
# mean_L_std_list = test.mean_L_std_list
# rho_vals_list = test.rho_vals_list

# rho_chrds = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
#                                       inputs.variables['rho_grid'].data, 
#                                       (zlos_chords2, ylos_chords2, xlos_chords2),
#                                       bounds_error = False, method='linear')

# rho_chrds2 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
#                                        inputs.variables['rho_grid'].data, 
#                                        (zlos, ylos, xlos),
#                                        bounds_error = False, method='linear')
# plot_items = ['te','ti','ne', 'TOR_ROT']
# plot_items = ['ti', 'TOR_ROT']
# fig10, ax10 = pt.subplots(ncols = len(plot_items), sharex = True, figsize=(20,4))
# for ind, name in enumerate(plot_items):
#     rho_vals = getattr(b,name+'_rho')
#     dat_vals = getattr(b,name)
#     for i in range(0, len(xlens)):
#         #ax10.errorbar(i, rho_vals_list[i][1], xerr=[rho_vals_list[i][0], rho_vals_list[i][2]], color=clr_cycle(i))
#         tmp_loc = np.argmin(np.abs(rho_vals - rho_vals_list[i][1]))
#         ch_te = np.interp(rho_vals_list[i][1], rho_vals,dat_vals)
#         tmp_val = [ch_te, ch_te, ch_te]
#         ax10[ind].plot(rho_vals_list[i],tmp_val,  color=clr_cycle(i), linewidth=2)
#         ax10[ind].plot(rho_vals_list[i][1], ch_te, 'o',  color=clr_cycle(i),)
#         i = len(xlens) - i-1
#         ax10[ind].plot(rho_chrds2[i], ch_te, 'xk')
#         ax10[ind].set_title(name)
#         ax10[ind].set_xlabel('rho')
#         #if ind==0:
#         #    ax10[0].plot(rho_vals_list[i],[i,i,i],  color=clr_cycle(i))
#     for rho_cur in rho_chrds:
#         ax10[ind].axvline(rho_cur)
#     ax10[ind].plot(rho_vals, dat_vals, '-', linewidth=0.2)
#     if name=='ti':
#         ax10[ind].plot(rho_chrds2,np.array(ti_list)*1., color='k',marker='x',linestyle='-')
#     if name=='TOR_ROT':
#         ax10[ind].plot(rho_chrds2,np.array(vel_list)*1000., color='k',marker='x',linestyle='-')
# ax10[1].set_xlim([0,1.2])
# fig10.canvas.draw();fig10.show()
# pt.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.05)
# #mlab.mesh(x_grid, y_grid, z_grid, scalars=S, colormap='YlGnBu', )
# #mlab.mesh(x_grid, y_grid, z_grid, colormap='YlGnBu', )
# for tmp in c.keys():
#     mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
#     dat = c[tmp]
#     for i in range(0,x_grid.shape[0],5):
#         print i
#         mlab.mesh(x_grid[i,:,:], y_grid[i,:,:], z_grid[i,:,:], scalars = dat[i,:,:], colormap='jet', vmin=0,vmax = 1.5)
#     for phi in np.linspace(0,2.*np.pi,10):
#         mlab.mesh(b.rgrid*np.sin(phi)*100,b.rgrid*np.cos(phi)*100, b.zgrid*100, scalars = d[tmp], colormap='jet', vmin=0,vmax = 1.5, opacity = 0.5)
# mlab.show()



# def calc_mean_L(L_frm_lens, val):
#     dL = L_frm_lens[1]-L_frm_lens[0]
#     rel_ind = np.invert(np.isnan(val))
#     rel_pts = L_frm_lens[rel_ind]
#     rel_vals = val[rel_ind]
#     tot_area = np.sum(dL*rel_vals)
#     rel_vals_prob = val[rel_ind]/ tot_area
#     mean_L = np.sum(rel_pts*rel_vals_prob*dL)
#     mean_dens = np.interp(mean_L, rel_pts, rel_vals)
#     #mean_dens = rel_vals[np.argmin(np.abs(rel_pts - mean_L))]
#     std_L = np.sqrt(np.sum(rel_vals_prob * (rel_pts - mean_L)**2))
#     print mean_L, mean_dens, tot_area
#     return mean_L, mean_dens, std_L

# def convert_L_to_rho(rho, L, vals):
#     rel_ind = np.invert(np.isnan(rho))
#     rel_L = L[rel_ind]
#     rel_rho = rho[rel_ind]
#     #for i in vals:

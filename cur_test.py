import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pt
import OMFITtree
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


    return halo_dens, los_wght, x_grid, y_grid, z_grid, r_grid, u_grid, v_grid, inputs, origin, neutrals



class interpolate_grid_profiles():
    def __init__(self, shot, time, run_id = 'def'):
        self.shot = shot
        self.time = time
        self.run_id = run_id
        HOME = os.environ['HOME']
        self.HOME = HOME
        self.eqdsk = OMFITtree.OMFITeqdsk(filename=HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/{}/g{}.{:05d}'.format(shot, time, run_id, shot, time))
        self.rgrid, self.zgrid = np.meshgrid(self.eqdsk['AuxQuantities']['R'],self.eqdsk['AuxQuantities']['Z'])

    def get_flux_values(self, r_values, z_values, use_rho = True,):
        self.r_values = r_values
        self.z_values = z_values
        self.inp_shape = r_values.shape
        self.r_values_flat = self.r_values.flatten()
        self.z_values_flat = self.z_values.flatten()

        if use_rho:
            self.interp_key = 'RHOpRZ'
            self.sav_key = 'RHO'
        else:
            self.interp_key = 'PSIRZ_NORM'
            self.sav_key = 'PSI'
        print 'interpolating to find rho at R,Z'
        self.eqdsk_flux = self.eqdsk['AuxQuantities'][self.interp_key]
        self.flux_values = scipy.interpolate.griddata((self.rgrid.flatten(),self.zgrid.flatten()),self.eqdsk_flux.flatten(),(self.r_values_flat,self.z_values_flat))
        self.flux_values = self.flux_values.reshape(self.inp_shape)
        
    def get_profile_RZ(self,):
        self.data_out = {}
        self.data_out_surf = {}
        for plot_key in ['ti','te','ne', 'trot']:
            print plot_key
            fname = self.HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/{}/d{}{}.{:05d}'.format(self.shot, self.time,self.run_id, plot_key, self.shot, self.time)
            if plot_key=='trot':
                plot_key = 'TOR_ROT'
            dat_obj = OMFITtree.OMFITidlSav(fname)['{}_str'.format(plot_key)]
            if plot_key=='ne':
                tmp1 = 'DENS'
                tmp2 = 'DENS'
            elif plot_key=='TOR_ROT':
                tmp1 = 'TOR_ROT'
                tmp2 = 'V_TOR_ROT'
            else:
                tmp1 = plot_key.upper()
                tmp2 = plot_key.upper()
            print dat_obj.keys()
            #print dat_obj['RHO_{}'.format(tmp1)]
            dat_psi = dat_obj['{}_{}'.format(self.sav_key,tmp1)]
            dat_dat = dat_obj[tmp2]
            setattr(self,plot_key,dat_dat)
            setattr(self,plot_key+'_rho',dat_psi)
            dat_interp = np.interp(self.flux_values, dat_psi, dat_dat,)
            dat_interp_surf = np.interp(self.eqdsk_flux, dat_psi, dat_dat,)
            self.data_out[plot_key] = dat_interp.copy()
            self.data_out_surf[plot_key] = dat_interp_surf.copy()
        return self.data_out, self.data_out_surf

b = interpolate_grid_profiles(155196,500)
HOME = os.environ['HOME']
halo_dens, los_wght, x_grid, y_grid, z_grid, u_grid,v_grid, r_grid,inputs,origin,neutrals  = get_los_data(dir=HOME + '/FIDASIM/RESULTS/D3D/155196/00500/MAIN_ION330/', run_id = 'def', plot = False)

alpha = inputs.variables['alpha'].data

def transform(x_grid, y_grid, z_grid, origin, alpha, unit_conversion=1./100):
    phi_tmp = np.arctan2(y_grid, x_grid)
    r_tmp = np.sqrt(x_grid**2 + y_grid**2)

    #Values in m
    x_grid2 = (r_tmp*np.cos(phi_tmp + alpha) + origin[0])*unit_conversion
    y_grid2 = (r_tmp*np.sin(phi_tmp + alpha) + origin[1])*unit_conversion
    z_grid2 = (z_grid + origin[2])*unit_conversion
    return x_grid2, y_grid2,z_grid2


x_grid2, y_grid2,z_grid2 =  transform(x_grid, y_grid, z_grid, origin, alpha)

b.get_flux_values(np.sqrt(x_grid2**2 + y_grid2**2),z_grid)
c, d = b.get_profile_RZ()
xlens = inputs.variables['xlens'].data
ylens = inputs.variables['ylens'].data
zlens = inputs.variables['zlens'].data

xlos = inputs.variables['xlos'].data
ylos = inputs.variables['ylos'].data
zlos = inputs.variables['zlos'].data

xlens2, ylens2, zlens2 = transform(xlens, ylens, zlens, origin, alpha)
xlos2, ylos2, zlos2 =  transform(xlos, ylos, zlos, origin, alpha)

n_chrds = 16

xlos_chords,ylos_chords,zlos_chords = [np.linspace(-133.64747, -134.95745, n_chrds)/100,
                                       np.linspace(172.84416, 186.51464 , n_chrds)/100,
                                       np.linspace(-1.4165587, -1.0956602,  n_chrds)/100]

xlens_chords,ylens_chords,zlens_chords = [np.linspace(-58.045200, -58.045200, n_chrds)/100,
                                          np.linspace(238.66320, 238.66320, n_chrds)/100,
                                          np.linspace(0.68220000, 0.68220000, n_chrds)/100]



xlos_chords2 = np.interp(xlos_chords,xlos2,xlos)
ylos_chords2 = np.interp(ylos_chords,ylos2,ylos)
zlos_chords2 = np.interp(zlos_chords,zlos2,zlos)

fig, ax = pt.subplots(ncols = 4,nrows = 6, sharex = True)
ax6 = ax[5,:]
ax5 = ax[4,:]
ax4 = ax[3,:]
ax3 = ax[2,:]
ax2 = ax[1,:]
ax1 = ax[0,:]

clr_cycle = gen.new_color_cycle(0,len(xlens))
plot_quants = ['rho_grid','te','dene','ti']
peaks = {i:{'x':[],'y':[]} for i in plot_quants}

def calc_mean_L(L_frm_lens, val):
    dL = L_frm_lens[1]-L_frm_lens[0]
    rel_ind = np.invert(np.isnan(val))
    rel_pts = L_frm_lens[rel_ind]
    rel_vals = val[rel_ind]
    tot_area = np.sum(dL*rel_vals)
    rel_vals_prob = val[rel_ind]/ tot_area
    mean_L = np.sum(rel_pts*rel_vals_prob*dL)
    mean_dens = np.interp(mean_L, rel_pts, rel_vals)
    #mean_dens = rel_vals[np.argmin(np.abs(rel_pts - mean_L))]
    std_L = np.sqrt(np.sum(rel_vals_prob * (rel_pts - mean_L)**2))
    print mean_L, mean_dens, tot_area
    return mean_L, mean_dens, std_L

def convert_L_to_rho(rho, L, vals):
    rel_ind = np.invert(np.isnan(rho))
    rel_L = L[rel_ind]
    rel_rho = rho[rel_ind]
    #for i in vals:

    
    pass

mean_L_list = []
mean_dens_list = []
mean_L_std_list = []
rho_vals_list = []
for i in range(0, len(xlens)):
    i = len(xlens) - i-1
    print i
    dz = zlos[i]- zlens[i]
    dy = ylos[i]- ylens[i]
    dx = xlos[i]- xlens[i]
    z_pts = np.linspace(zlens[i],zlens[i] + 2.5*dz,800)
    y_pts = np.linspace(ylens[i],ylens[i] + 2.5*dy,800)
    x_pts = np.linspace(xlens[i],xlens[i] + 2.5*dx,800)

    rho = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                    inputs.variables['rho_grid'].data, 
                                    (z_pts, y_pts, x_pts),
                                    bounds_error = False, method='linear')
    #ax_tmp.plot(rho,q)
    q2 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                  neutrals.variables['halodens'].data[0,:,:,:],
                                  (z_pts, y_pts, x_pts),
                                  bounds_error = False, method='linear')
    q3 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                  neutrals.variables['halodens'].data[1,:,:,:],
                                  (z_pts, y_pts, x_pts),
                                  bounds_error = False, method='linear')
    q4 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                  neutrals.variables['halodens'].data[2,:,:,:],
                                  (z_pts, y_pts, x_pts),
                                  bounds_error = False, method='linear')
    q5 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                  neutrals.variables['halodens'].data[3,:,:,:],
                                  (z_pts, y_pts, x_pts),
                                  bounds_error = False, method='linear')
    #mean_L_list.append(tmp[0])
    #mean_L_list_loc.append(tmp[0])
    max_val = np.max(q5[np.invert(np.isnan(q5))])
    R = np.sqrt((x_pts-x_pts[0])**2 + (y_pts-y_pts[0])**2 + (z_pts-z_pts[0])**2)
    tmp_x = R - R[0]

    tmp = calc_mean_L(tmp_x, q5)
    mean_L_list.append(tmp[0])
    mean_dens_list.append(tmp[1])
    mean_L_std_list.append(tmp[2])

    ax5[0].errorbar(mean_L_list[-1], mean_dens_list[-1], xerr=tmp[2]/2,color=clr_cycle(i))
    ax5[0].plot(mean_L_list[-1], mean_dens_list[-1], 'x',color=clr_cycle(i))

    vals = [mean_L_list[-1] - tmp[2]/2,mean_L_list[-1], mean_L_list[-1] + tmp[2]/2]
    rho_vals_list.append(np.interp(vals, tmp_x, rho))

    for quant,ax_tmp, ax_tmp2, ax_tmp3, ax_tmp4, ax_tmp5, ax_tmp6 in zip(plot_quants,ax1, ax2, ax3, ax4, ax5, ax6):
        q = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                      inputs.variables[quant].data, 
                                      (z_pts, y_pts, x_pts),
                                      bounds_error = False, method='linear')
        tmp_y = q*q5/max_val
        linewidth=0.3
        ax_tmp.plot(R - R[0], q, color = clr_cycle(i), linewidth=linewidth)
        ax_tmp.set_title(quant)
        ax_tmp2.plot(R - R[0], q2, color=clr_cycle(i), linewidth=linewidth)
        ax_tmp3.plot(R - R[0], q3, color=clr_cycle(i), linewidth=linewidth)
        ax_tmp4.plot(R - R[0], q4, color=clr_cycle(i), linewidth=linewidth)
        ax_tmp5.plot(R - R[0], q5, color=clr_cycle(i), linewidth=linewidth)
        ax_tmp6.plot(tmp_x, tmp_y, color=clr_cycle(i), linewidth=linewidth)
        valid = np.invert(np.isnan(q5))
        peaks[quant]['y'].append(np.max(tmp_y[valid]))
        peaks[quant]['x'].append(tmp_x[valid][np.argmax(tmp_y[valid])])
    for i in range(len(plot_quants)):
        ax[-1,i].plot(peaks[plot_quants[i]]['x'],peaks[plot_quants[i]]['y'],'o-')
        ax[-1,i].set_xlabel('Distance frm lens')

#max_loc, min_loc = np.argmin(rho), np.argmax(rho)
        
        #ax_tmp.plot([rho[min_loc],rho[max_loc]],[q[min_loc],q[max_loc]],'o')
pt.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.05)
fig.canvas.draw();fig.show()

fig10, ax10 = pt.subplots(ncols = 4, sharex = True)
rho_chrds = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                      inputs.variables['rho_grid'].data, 
                                      (zlos_chords2, ylos_chords2, xlos_chords2),
                                      bounds_error = False, method='linear')

rho_chrds2 = scipy.interpolate.interpn((z_grid[:,0,0],y_grid[0,:,0], x_grid[0,0,:]),
                                       inputs.variables['rho_grid'].data, 
                                       (zlos, ylos, xlos),
                                       bounds_error = False, method='linear')

for ind, name in enumerate(['te','ti','ne', 'TOR_ROT']):
    rho_vals = getattr(b,name+'_rho')
    dat_vals = getattr(b,name)
    for i in range(0, len(xlens)):
        #ax10.errorbar(i, rho_vals_list[i][1], xerr=[rho_vals_list[i][0], rho_vals_list[i][2]], color=clr_cycle(i))
        tmp_loc = np.argmin(np.abs(rho_vals - rho_vals_list[i][1]))
        ch_te = np.interp(rho_vals_list[i][1], rho_vals,dat_vals)
        tmp_val = [ch_te, ch_te, ch_te]
        ax10[ind].plot(rho_vals_list[i],tmp_val,  color=clr_cycle(i), linewidth=2)
        ax10[ind].plot(rho_vals_list[i][1], ch_te, 'o',  color=clr_cycle(i),)
        i = len(xlens) - i-1
        ax10[ind].plot(rho_chrds2[i], ch_te, 'xk')
        ax10[ind].set_title(name)
        ax10[ind].set_xlabel('rho')
        #if ind==0:
        #    ax10[0].plot(rho_vals_list[i],[i,i,i],  color=clr_cycle(i))
    for rho_cur in rho_chrds:
        ax10[ind].axvline(rho_cur)
    ax10[ind].plot(rho_vals, dat_vals, '-', linewidth=0.2)
ax10[1].set_xlim([0,1.2])
fig10.canvas.draw();fig10.show()

plot_mayavi = False
if plot_mayavi:
    from mayavi import mlab

    mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    #mlab.mesh(x_grid, y_grid, z_grid, scalars=S, colormap='YlGnBu', )
    #mlab.mesh(x_grid, y_grid, z_grid, colormap='YlGnBu', )
    for i in range(0,x_grid2.shape[0],5):
        print i
        #mlab.mesh(x_grid[i,:,:], y_grid[i,:,:], z_grid[i,:,:], scalars = b.flux_values[i,:,:], colormap='jet', vmin=0,vmax = 1.5 , opacity = 0.5)
        clr_dat = inputs.variables['rho_grid'].data[i,:,:]
        #clr_dat = b.flux_values[i,:,:]
        #mlab.mesh(u_grid[i,:,:]/100., v_grid[i,:,:]/100., z_grid[i,:,:]/100., scalars = clr_dat, colormap='jet', vmin=0,vmax = 1.5 , opacity = 0.5)
        mlab.mesh(x_grid2[i,:,:], y_grid2[i,:,:], z_grid2[i,:,:], scalars = clr_dat, colormap='jet', vmin=0,vmax = 1.5 , opacity = 0.5)
    #for phi in np.linspace(0,2.*np.pi,10):
    #for phi in np.linspace(0,2.*np.pi,10):
    #for phi in [0,90,180,270,330]:
    for phi in [330, 330+180]:
        phi = np.deg2rad(-360+90-phi)
        mlab.mesh(b.rgrid*np.cos(phi),b.rgrid*np.sin(phi), b.zgrid, scalars = b.eqdsk_flux, colormap='jet', vmin=0,vmax = 1.5, opacity = 0.5)
    x_orig,y_orig,z_orig = inputs.variables['origin'].data


    b = inputs.variables['rho_grid'].data.copy()





    for i in range(0,inputs.variables['xlens'].shape[0],5):
       mlab.plot3d([xlens2[i],xlos2[i]],[ylens2[i], ylos2[i]],[zlens2[i],zlos2[i]], tube_radius = 0.5/100)

    mlab.show()
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

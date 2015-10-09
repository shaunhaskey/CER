import numpy as np
import scipy.io.netcdf as netcdf
import os, copy
import OMFITtree
import matplotlib.pyplot as pt
import scipy
from scipy.optimize import curve_fit

shot = 158676
#shot = 155196
start_time = 520
end_time = 530
start_time = 600
start_time = 900
end_time = 925
start_time = 1200
end_time = 1225
start_time =1400
end_time = start_time + 125
offset = 0
jump = 5
offset = 0
jump = 1
rel_times = range(start_time+offset, end_time,jump)
#import cPickle as pickle
#master_dict = pickle.load(file('/u/haskeysr/fida_sim_dict.pickle','r'))
#master_dict = pickle.load(file('/u/haskeysr/fida_sim_dict.pickle','r'))

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def get_data(dir, run_id = 'def', plot = False):
    #if 'inputs' not in locals():
    inputs = netcdf.netcdf_file(dir + '{}_inputs.cdf'.format(run_id),'r',mmap = False)
    spectra = netcdf.netcdf_file(dir + '{}_spectra.cdf'.format(run_id),'r',mmap = False)
    #neutrals = netcdf.netcdf_file(dir + '{}_neutrals.cdf'.format(run_id),'r',mmap = False)
    #weights = netcdf.netcdf_file(dir + '{}_fida_weights.cdf'.format(run_id),'r',mmap = False)
    neutrals=weights=None
    fida = +spectra.variables['fida'].data
    wave = +spectra.variables['lambda'].data
    halo = +spectra.variables['halo'].data
    BE = [+spectra.variables['full'].data, +spectra.variables['half'].data, +spectra.variables['third'].data]
    if plot:
        fig, ax = pt.subplots(nrows = len(BE)+1)
        for i in range(fida.shape[0]):
            for ax_tmp, dat_tmp in zip(ax,BE):
                ax_tmp.plot(wave,dat_tmp[i,:])
            ax[len(BE)].plot(wave,halo[i,:])
        fig.canvas.draw();fig.show()
    #    ax[0].plot(fida[i,:])
    i = np.argmin(np.abs([np.mean(inputs.variables['z_grid'].data[i,:,:]) for i in range(inputs.variables['z_grid'].data.shape[0])]))
    z = inputs.variables['z_grid'].data[i,:,:]
    x_grid = inputs.variables['x_grid'].data[i,:,:]
    y_grid = inputs.variables['y_grid'].data[i,:,:]
    if plot:
        dat_grid = neutrals.variables['halodens'].data[0,i,:,:]
        n_halos = neutrals.variables['halodens'].shape[0]
        fig, ax = pt.subplots(nrows = n_halos)
        for j in range(n_halos):
            im = ax[j].pcolormesh(x_grid, y_grid, neutrals.variables['halodens'].data[j,i,:,:])
            pt.colorbar(im,ax=ax[j])
        fig.canvas.draw();fig.show
    return inputs, neutrals, spectra, weights

def calc_ti_vel(data, wave, plot = True, ax_tmp = None):
    hist = data/10. # convert to ph/s-m^2/sR-A
    lambda0 = 6561.0
    # ;; Now convert to CCD counts for fitting
    # ;; We measure counts in a time interval.
    # ;; The total count rate CR = TOTAL(counts[pix])/tinteg
    # ;; The total brightness B = CR/calib
    # ;; Then radiance is counts[pix]/TOTAL(counts[pix])*B/D.
    # ;; Or    R = counts/(tinteg*calib*disp)
    # ;; 
    # ;; We then do this backwards with 
    # ;; counts = R*tinteg*calib*disp

    tinteg = 2.5e-3 # s
    calib = 1.0e-7  # (count rate) / ph/s-cm**2-sR
    calib /= (100.0)**2 # (count rate) / ph/s-m**2-sR
    disp = 0.18 #;; A/pix
    r_to_counts = tinteg*calib*disp
    counts = hist*r_to_counts

    #Initial guess for Gaussian fit
    p0 = [np.max(counts), wave[np.argmax(counts)], 20]
    print p0
    #Fit with a Gaussian
    coeff, var_matrix = curve_fit(gauss, wave, counts, p0=p0)

    # Get the fitted curve
    hist_fit = gauss(wave, *coeff)
    print np.max(hist_fit), np.min(hist_fit), coeff
    if plot:
        if ax_tmp==None:
            fig, ax = pt.subplots()
        else:
            ax = ax_tmp
        ax.plot(wave, counts)#, label='Test data')
        ax.plot(wave, hist_fit)#, label='Fitted data')
        #ax.legend(loc='best')
        if ax_tmp==None:
            fig.canvas.draw();fig.show()
    # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:

    print 'Fitted mean = ', coeff[1]
    print 'Fitted standard deviation = ', coeff[2]

    #Now extract temperature and rotation velocity out of the data
    c = 2.9979e8
    # Temperature
    eV_to_J = 1.6022e-19#6.24150934e18
    md = 3.3452e-27
    sigma = coeff[2]
    ti = (sigma/lambda0)**2 * md * c**2
    ti /= eV_to_J
    ti /= 1.e3
    vlos = c*(coeff[1]-lambda0)/lambda0 *1.e-3

    #1.e-3 * eV_to_J #;; J->keV
    #ti /= 1.e3 ;; keV
    print ti, vlos
    print '########## finished ###############'
    return ti, vlos, coeff

print 'hello world'
nrows = 5
ncols = 5
n_plots = nrows * ncols
fig_probes, ax_probes = pt.subplots(ncols = 5, sharex = True, sharey=True)
fig_flux, ax_flux = pt.subplots(ncols = ncols, nrows = nrows, sharex = True, sharey = True)
fig_flux2, ax_flux2 = pt.subplots(ncols = ncols, nrows = nrows, sharex = True, sharey = True)
ax_flux2 = ax_flux2.flatten()
ax_flux = ax_flux.flatten()
def plot_profiles():
    fig, ax = pt.subplots(nrows = 3, sharex = True)
    for num in rel_times:#range(start_time,end_time):
        print num
        HOME = os.environ['HOME']
        dir =  HOME + '/FIDASIM/RESULTS/D3D/{}/{:05d}/MAIN_ION330/'.format(shot, num) 
        run_id = 'def'
        #inputs, neutrals, spectra, weights = get_data(dir, run_id = 'def', plot = False)
        #halo = +spectra.variables['halo'].data
        a = OMFITtree.OMFITeqdsk(filename=HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/def/g{}.{:05d}'.format(shot, num, shot, num))
        #a = OMFITtree.OMFITeqdsk(filename='/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/g158676.{:05d}'.format(900,900 ))
        b = a['AuxQuantities']
        rgrid,zgrid = np.meshgrid(b['R'],b['Z'])
        r_new = np.linspace(1.4, 2.5, 300)
        z_new = r_new * 0
        use_rho = True

        if use_rho:
            interp_key = 'RHOpRZ'
            sav_key = 'RHO'
        else:
            interp_key = 'PSIRZ_NORM'
            sav_key = 'PSI'
        flux_new = scipy.interpolate.griddata((rgrid.flatten(),zgrid.flatten()),b[interp_key].flatten(),(r_new,z_new))
        print num
        for plot_key,ax_tmp in zip(['ti','te','ne'],ax):
            print plot_key, ax_tmp
            fname = HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/def/d{}{}.{:05d}'.format(shot, num,plot_key, shot, num)
            #fname = '/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/d{}158676.{:05d}'.format(900,plot_key, 900)
            dat_obj = OMFITtree.OMFITidlSav(fname)['{}_str'.format(plot_key)]
            if plot_key=='ne':
                tmp = 'DENS'
            else:
                tmp = plot_key.upper()
            print dat_obj.keys()
            print dat_obj['RHO_{}'.format(tmp)]
            dat_psi = dat_obj['{}_{}'.format(sav_key,tmp)]
            dat_dat = dat_obj[tmp]
            dat_interp = np.interp(flux_new, dat_psi, dat_dat,)
            ax_tmp.plot(r_new, dat_interp,'b-')
        #ax_tmp[1].pcolor(R_array, Z_array, psi_array)
    fig.canvas.draw();fig.show()


import cer_funcs as CER
def mtanh_wrapper(x,*p):
    core = [0.02]#, 0.001]
    edge = []#-0.02]
    p = [i for i in p]
    p = p + [len(core), len(edge), False] + core + edge
    print p
    y = CER.mtanh(x,*p)
    return y

#plot_profiles()

fit_coeffs = []
fit_coeffs_orig = []
results = {}
for num in rel_times:
    HOME = os.environ['HOME']
    #num = 200
    dir =  HOME + '/158676/00001/MAIN_ION330/'
    dir =  HOME + '/FIDASIM/RESULTS/D3D/{}/{:05d}/MAIN_ION330/'.format(shot, num)
    run_id = 'def'
    try:
        inputs, neutrals, spectra, weights = get_data(dir, run_id = 'def', plot = False)
        read_data = True
    except IOError:
        print 'COULDN"T READ INPUT FOR ',num
        read_data = False
    if read_data:
        use_rho = True
        if use_rho:
            interp_key = 'RHOpRZ'
            interp_key = 'RHORZ'
            sav_key = 'RHO'
        else:
            interp_key = 'PSIRZ_NORM'
            sav_key = 'PSI'

        halo = +spectra.variables['halo'].data
        plot_key = 'ti'
        a = OMFITtree.OMFITeqdsk(filename=HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/def/g{}.{:05d}'.format(shot, num,shot, num))
        #a = OMFITtree.OMFITeqdsk(filename='/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/g158676.{:05d}'.format(900, 900))
        prof_name = HOME + '/gaprofiles/f90fidasim/{}/{:05d}/MAIN_ION330/def/d{}{}.{:05d}'.format(shot, num,plot_key, shot, num)
        #prof_name = '/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/d{}158676.{:05d}'.format(900,plot_key, 900)
        prof_dat = OMFITtree.OMFITidlSav(prof_name)['{}_str'.format(plot_key)]
        prof_flux = prof_dat['{}_{}'.format(sav_key, plot_key.upper())]
        prof = prof_dat[plot_key.upper()]
        psi_list = []
        Z_list = []
        R_list = []
        plot = False
        if plot:
            fig, ax = pt.subplots()
            for i in a['fluxSurfaces']['flux'].keys():
                ax.plot(a['fluxSurfaces']['flux'][i]['R'], a['fluxSurfaces']['flux'][i]['Z'],'--')

            fig.canvas.draw();fig.show()
        r_new = np.linspace(1.4, 2.5, 300)
        z_new = r_new * 0
        b = a['AuxQuantities']
        rgrid,zgrid = np.meshgrid(b['R'],b['Z'])
        r_probes = spectra.variables['radius'].data/100.
        flux_new = scipy.interpolate.griddata((rgrid.flatten(),zgrid.flatten()),b[interp_key].flatten(),(r_new,z_new))
        flux_probe = scipy.interpolate.griddata((rgrid.flatten(),zgrid.flatten()),b[interp_key].flatten(),(r_probes,r_probes * 0))
        #flux_new = scipy.interpolate.griddata((rgrid.flatten(),zgrid.flatten()),b['PSIRZ_NORM'].flatten(),(r_new,z_new))
        prof_interp = np.interp(flux_new, prof_flux, prof,)


        if plot:
            fig_tmp, ax_tmp = pt.subplots(nrows = 2)
            ax_tmp[0].plot(r_new,flux_new)
            #ax_tmp[1].pcolor(R_array, Z_array, psi_array)
            fig_tmp.canvas.draw();fig_tmp.show()

        wave = +spectra.variables['lambda'].data
        bin_centres = wave
        wave = wave*10.  #Angstroms
        #halo[0,:]
        plot_spectrum = False
        plot_spectra = False
        ti_list = []; vel_list = []; coeffs_list = []
        if plot_spectra:
            fig_tmp, ax_tmp = pt.subplots(nrows=5,ncols=5, sharex = True)
            ax_tmp = ax_tmp.flatten()
        for i in range(halo.shape[0]):
            if (i%5==0) and plot_spectra:
                ax_tmp_in = ax_tmp[i/5]
                plot_spectrum_tmp = True
            else:
                ax_tmp_in = None
                plot_spectrum_tmp = False
            ti, vel, coeffs = calc_ti_vel(halo[i,:], wave, plot=plot_spectrum_tmp, ax_tmp = ax_tmp_in)
            if (i%5==0) and plot_spectra:
                ax_tmp[i/5].text(ax_tmp[i/5].get_xlim()[0],0,'{:.3f}'.format(r_probes[i]),verticalalignment='bottom',horizontalalignment='left')
            ti_list.append(ti)
            vel_list.append(vel)
            coeffs_list.append(coeffs)


        ax_probes[0].plot(r_probes,ti_list, marker='.',linestyle='-')
        ax_flux[num%n_plots].plot(flux_probe,ti_list, marker='.',linestyle='-')
        ax_flux[num%n_plots].plot(prof_flux, prof,'b-')
        ax_flux2[int((num-start_time)/n_plots)].plot(flux_probe,ti_list, marker='.',linestyle='-')
        ax_flux2[int((num-start_time)/n_plots)].plot(prof_flux, prof,'b-')
        guess = [1.5, 0.2, 0.9, 0.05,]
        coeff, var_matrix = curve_fit(mtanh_wrapper, flux_probe, ti_list, p0=guess,)
        guess_orig = [1.5, 0.2, 0.9, 0.05,]
        coeff_orig, var_matrix_orig = curve_fit(mtanh_wrapper, prof_flux, prof, p0=guess_orig,)
        results[num] = {'ti_list':copy.deepcopy(ti_list), 
                        'vel_list':copy.deepcopy(vel_list),
                        'coeffs_list':copy.deepcopy(coeffs_list),
                        'mtanh_diag':copy.deepcopy(coeff),
                        'mtanh_real':copy.deepcopy(coeff_orig)}
        fit_coeffs_orig.append(coeff_orig)
        fit_coeffs.append(coeff)
        if plot_spectra:
            fig_tmp.suptitle('Top:{:.3f},Width:{:.3f}'.format(fit_coeffs[-1][0], fit_coeffs[-1][3]))
            fig_tmp.tight_layout()
            fig_tmp.canvas.draw();fig_tmp.show()
        y2 = mtanh_wrapper(flux_probe, *coeff)
        y3 = mtanh_wrapper(prof_flux, *coeff_orig)
        ax_flux[num%n_plots].plot(flux_probe, y2,'-.')
        ax_flux[num%n_plots].plot(prof_flux, y3,'r--')
        ax_flux2[int((num-start_time)/n_plots)].plot(flux_probe, y2)
        ax_flux2[int((num-start_time)/n_plots)].plot(prof_flux, y3,'r--')
        ax_probes[1].plot(r_probes,vel_list,'-o')
        ax_probes[0].set_ylim([0,4])
        ax_probes[0].set_xlim([1.45, 2.32])
        #ax_probes[0].plot(int1*100,r_new2,'o')
        ax_probes[0].plot(r_new, prof_interp,'b-')
        #ax[].set_ylim([0,4])

        if plot:
            fig1, ax1 = pt.subplots()
            ax1.plot(prof_flux, prof,'x')
            ax1.plot(flux_new, prof_interp,'o')
            fig1.canvas.draw();fig1.show()
    if num%5==0:
        fig_flux.canvas.draw();fig_flux.show()
        fig_flux2.canvas.draw();fig_flux2.show()

fig_probes.canvas.draw();fig_probes.show()
fig_flux.canvas.draw();fig_flux.show()
fig_flux2.canvas.draw();fig_flux2.show()
top = [i[0] for i in fit_coeffs]
offset = [i[1] for i in fit_coeffs]
sym = [i[2] for i in fit_coeffs]
width = [i[3] for i in fit_coeffs]
width_orig = [i[3] for i in fit_coeffs_orig]
top_orig = [i[0] for i in fit_coeffs_orig]
fig, ax = pt.subplots(ncols = 2)
ax[0].plot(width,'-o')
ax[0].plot(width_orig,'-o')
ax[1].plot(top,'-o')
ax[1].plot(top_orig,'-o')
fig.canvas.draw();fig.show()
ax[1].set_ylim([0,ax[1].get_ylim()[1]])
fig.canvas.draw();fig.show()

fig, ax = pt.subplots(nrows = 2)
max_real = []; max_diag = []
#for i in range(1400,1425):#results.keys():
for i in results.keys():
    mtanh_real = mtanh_wrapper(prof_flux, *results[i]['mtanh_real'])
    mtanh_diag = mtanh_wrapper(prof_flux, *results[i]['mtanh_diag'])
    max_real.append(np.max(np.abs(np.diff(mtanh_real))))
    max_diag.append(np.max(np.abs(np.diff(mtanh_diag))))
    np.diff(mtanh_real)
    ax[0].plot(prof_flux, mtanh_real)
    ax[0].plot(prof_flux, mtanh_diag, '--')
ax[1].plot(max_real)
ax[1].plot(max_diag)
fig.canvas.draw();fig.show()

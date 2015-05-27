import numpy as np
import scipy.io.netcdf as netcdf
import os
import OMFITtree
#import OMFITidlSav
import matplotlib.pyplot as pt
import numpy
import scipy
from scipy.optimize import curve_fit


def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))



def get_data(dir, run_id = 'def', plot = False):
    #if 'inputs' not in locals():
    inputs = netcdf.netcdf_file(dir + '{}_inputs.cdf'.format(run_id),'r',mmap = False)
    neutrals = netcdf.netcdf_file(dir + '{}_neutrals.cdf'.format(run_id),'r',mmap = False)
    spectra = netcdf.netcdf_file(dir + '{}_spectra.cdf'.format(run_id),'r',mmap = False)
    weights = netcdf.netcdf_file(dir + '{}_fida_weights.cdf'.format(run_id),'r',mmap = False)
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
    dat_grid = neutrals.variables['halodens'].data[0,i,:,:]
    n_halos = neutrals.variables['halodens'].shape[0]
    if plot:
        fig, ax = pt.subplots(nrows = n_halos)
        for j in range(n_halos):
            im = ax[j].pcolormesh(x_grid, y_grid, neutrals.variables['halodens'].data[j,i,:,:])
            pt.colorbar(im,ax=ax[j])
        fig.canvas.draw();fig.show
    return inputs, neutrals, spectra, weights


def calc_ti_vel(data, wave, plot = True):
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
        fig, ax = pt.subplots()
        ax.plot(wave, counts, label='Test data')
        ax.plot(wave, hist_fit, label='Fitted data')
        ax.legend(loc='best')
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
    return ti, vlos

fig_probes, ax_probes = pt.subplots(nrows = 2, sharex = True)
for num in range(200,225):
    HOME = os.environ['HOME']
    #num = 200
    dir =  HOME + '/158676/00001/MAIN_ION330/'
    dir =  HOME + '/FIDASIM/RESULTS/D3D/158676/{:05d}/MAIN_ION330/'.format(num)
    run_id = 'def'
    inputs, neutrals, spectra, weights = get_data(dir, run_id = 'def', plot = False)
    halo = +spectra.variables['halo'].data
    a = OMFITtree.OMFITeqdsk(filename='/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/g158676.{:05d}'.format(num, num))
    te_name = '/u/haskeysr/gaprofiles/f90fidasim/158676/{:05d}/MAIN_ION330/def/dte158676.{:05d}'.format(num, num)
    te_dat = OMFITtree.OMFITidlSav(te_name)['te_str']
    #te = scipy.io.readsav(te_name,python_dict=True)
    #te_psi_norm = te_dat['te_str']['TE_DATA'][0]['PSI_NORM'][0]
    te_psi = te_dat['PSI_TE']
    te = te_dat['TE']
    #te_psi_norm = te['te_str']['PSI_TE'][0]
    #edge_loc = np.argmin(np.abs(te['te_str']['RHO_TE'][0] - te['te_str']['MAXRHO'][0]))
    #te_psi_norm = te_psi_norm - np.min(te_psi_norm)
    #te_psi_norm /= te_psi_norm[edge_loc]
    psi_list = []
    Z_list = []
    R_list = []
    R_array = np.array([])
    Z_array = np.array([])
    psi_array = np.array([])
    #new_vals = np.linspace()

    plot = False
    if plot:
        fig, ax = pt.subplots()
        for i in a['fluxSurfaces']['flux'].keys():
            ax.plot(a['fluxSurfaces']['flux'][i]['R'], a['fluxSurfaces']['flux'][i]['Z'],'--')

        fig.canvas.draw();fig.show()
    for i in a['fluxSurfaces']['flux'].keys():
        psi_list.append(np.array(a['fluxSurfaces']['flux'][i]['psi']* a['fluxSurfaces']['flux'][i]['Z']))
        R_array = np.append(np.array(a['fluxSurfaces']['flux'][i]['R']),R_array)
        Z_array = np.append(np.array(a['fluxSurfaces']['flux'][i]['Z']),Z_array)
        psi_array = np.append(a['fluxSurfaces']['flux'][i]['psi']+0* np.array(a['fluxSurfaces']['flux'][i]['Z']),psi_array)
        Z_list.append(np.array(a['fluxSurfaces']['flux'][i]['Z']))
        R_list.append(np.array(a['fluxSurfaces']['flux'][i]['R']))

    #Normalise psi_array
    psi_array = psi_array  - np.min(psi_array)
    psi_array /= np.max(psi_array)
    #Values to interpolate the te profile onto
    r_new = np.linspace(np.min(R_array), np.max(R_array),300)
    z_new = r_new * 0
    psi_new = scipy.interpolate.griddata((R_array,Z_array),psi_array,(r_new,z_new))
    te_interp = np.interp(psi_new, te_psi, te,)
    if plot:
        fig_tmp, ax_tmp = pt.subplots(nrows = 2)
        ax_tmp[0].plot(r_new,psi_new)
        #ax_tmp[1].pcolor(R_array, Z_array, psi_array)
        fig_tmp.canvas.draw();fig_tmp.show()

    wave = +spectra.variables['lambda'].data
    bin_centres = wave
    wave = wave*10.  #Angstroms
    #halo[0,:]
    plot_spectrum = False
    ti_list = []; vel_list = []
    for i in range(halo.shape[0]):
        ti, vel = calc_ti_vel(halo[i,:], wave, plot=plot_spectrum)
        ti_list.append(ti)
        vel_list.append(vel)

    r_probes = spectra.variables['radius'].data/100.
    ax_probes[0].plot(r_probes,ti_list,'-o')
    ax_probes[1].plot(r_probes,vel_list,'-o')
    ax_probes[0].set_ylim([0,4])
    ax_probes[0].set_xlim([1.45, 2.32])
    #ax_probes[0].plot(int1*100,r_new2,'o')
    ax_probes[0].plot(r_new, te_interp,'b-')
    #ax[].set_ylim([0,4])
    fig_probes.canvas.draw();fig_probes.show()

    if plot:
        fig1, ax1 = pt.subplots()
        ax1.plot(te_psi, te,'x')
        ax1.plot(psi_new, te_interp,'o')
        fig1.canvas.draw();fig1.show()
# # Define some test data which is close to Gaussian
# data = numpy.random.normal(size=100000)

# hist, bin_edges = numpy.histogram(data, density=True,bins = 300)
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# # Define model function to be used to fit to the data above:
# # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
# p0 = [1., 0., 1.]

# coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# # Get the fitted curve
# hist_fit = gauss(bin_centres, *coeff)
# plt.figure()
# plt.plot(bin_centres, hist, label='Test data')
# plt.plot(bin_centres, hist_fit, label='Fitted data')

# # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
# print 'Fitted mean = ', coeff[1]
# print 'Fitted standard deviation = ', coeff[2]


# plt.show()


#b['ne_str']['DENS']
#b['ne_str']['RHO_DENS']
#b['ne_str']['PSI_DENS']
#ax.plot(b['ne_str']['RHO_DENS'][0],b['ne_str']['DENS'][0],'-')
#te['ne_str']['PSI_TE']
#te['ne_str']['PSI_TE']

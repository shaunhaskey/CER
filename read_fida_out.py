import numpy as np
import scipy.io.netcdf as netcdf
import os
import matplotlib.pyplot as pt

def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))


HOME = os.environ['HOME']
dir =  HOME + '/158676/00001/MAIN_ION330/'
run_id = 'def'
if 'inputs' not in locals():
    inputs = netcdf.netcdf_file(dir + '{}_inputs.cdf'.format(run_id),'r')
    neutrals = netcdf.netcdf_file(dir + '{}_neutrals.cdf'.format(run_id),'r')
    spectra = netcdf.netcdf_file(dir + '{}_spectra.cdf'.format(run_id),'r')
    weights = netcdf.netcdf_file(dir + '{}_fida_weights.cdf'.format(run_id),'r')
fida = +spectra.variables['fida'].data
wave = +spectra.variables['lambda'].data
halo = +spectra.variables['halo'].data
BE = [+spectra.variables['full'].data, +spectra.variables['half'].data, +spectra.variables['third'].data]
fig, ax = pt.subplots(nrows = len(BE)+1)

#inputs.close()
#neutrals.close()
#spectra.close()

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
fig, ax = pt.subplots(nrows = n_halos)
for j in range(n_halos):
    im = ax[j].pcolormesh(x_grid, y_grid, neutrals.variables['halodens'].data[j,i,:,:])
    pt.colorbar(im,ax=ax[j])
fig.canvas.draw();fig.show




import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define some test data which is close to Gaussian
data = numpy.random.normal(size=100000)

hist, bin_edges = numpy.histogram(data, density=True,bins = 300)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# Define model function to be used to fit to the data above:
# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [1., 0., 1.]

coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
hist_fit = gauss(bin_centres, *coeff)
plt.figure()
plt.plot(bin_centres, hist, label='Test data')
plt.plot(bin_centres, hist_fit, label='Fitted data')

# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
print 'Fitted mean = ', coeff[1]
print 'Fitted standard deviation = ', coeff[2]


plt.show()
bin_centres = wave
hist = halo[0,:] / 10 # convert to ph/s-m^2/sR-A
wave = wave*10.  #Angstroms
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


p0 = [np.max(hist), wave[np.argmax(hist)], 20]
print p0
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
hist_fit = gauss(bin_centres, *coeff)
plt.figure()
plt.plot(bin_centres, hist, label='Test data')
plt.plot(bin_centres, hist_fit, label='Fitted data')

# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
print 'Fitted mean = ', coeff[1]
print 'Fitted standard deviation = ', coeff[2]

c = 2.99792458e8
# Temperature
eV_to_J = 6.24150934e18
md = 3.34358348e-27
sigma = coeff[1]
ti = (sigma/lambda0)**2 * md * c**2
ti *= 1.e-3 * eV_to_J #;; J->keV
#ti /= 1.e3 ;; keV
print ti

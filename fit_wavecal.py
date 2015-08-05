import scipy 
import numpy as np
import matplotlib.pyplot as pt
from sklearn import mixture
import pyfusion.clustering.clustering as clust
a = np.loadtxt('wavecal.txt')
a = a.flatten()

subrange = np.abs(a[400:600])#[98:150]
x_axis = np.arange(len(subrange))
q = scipy.stats.rv_discrete(name='custm',values=(x_axis,subrange/np.sum(subrange)))
fig, ax = pt.subplots()
ax.hist(q.rvs(size=100000),len(x_axis))
fig.canvas.draw();fig.show()
gmm = mixture.GMM(n_components=5, covariance_type = 'diag', n_init=2, n_iter = 150, tol=0.0001)

#x_axis_new = np.linspace(np.min(x_axis),np.max(x_axis),len(x_axis)*2)
#interp = np.interp(x_axis_new, x_axis, subrange)

dat = []
for i in x_axis:
    print subrange[i], i
    if subrange[i]>0:dat.extend([i]*int(subrange[i]))

# dat = []
# for i in x_axis_new:
#     print subrange[i], i
#     if subrange[i]>0:dat.extend([i]*int(subrange[i]))

#dat.extend(dat)
dat = np.array(dat)
dat += 1*(np.random.rand(len(dat))-0.5)


dat = q.rvs(size=1000000)
dat = np.array(dat)[:,np.newaxis]

gmm.fit(dat)
weights = gmm.weights_
means = gmm.means_
cov = gmm.covars_
#cov = gmm._get_covars()
#Y_ = gmm.predict(dat)
from scipy.stats.distributions import norm
overall_tot = np.zeros(len(x_axis),dtype=float)
vals = []
fig, ax = pt.subplots(nrows = 3, sharex = True)
for i in range(len(weights)):#zip(weights, means[:.0], cov[:,0]):
    weight = weights[i]
    mean = means[i,0]
    std  = np.sqrt(cov[i,0])
    print weight, mean, std
    vals.append(weight * norm(loc=mean, scale=std).pdf(x_axis))
    #vals.append(weight * norm(loc=mean, scale=std).pdf(x_axis))
    overall_tot+=vals[-1]
    ax[1].plot(x_axis, vals[-1],'--')
ax[0].plot(subrange/np.sum(subrange))
#ax[0].plot(x_axis_new, interp)
z = ax[0].hist(dat,len(x_axis), normed=True)
ax[1].plot(x_axis, overall_tot)#*np.sum(subrange)/np.sum(overall_tot),'kx')
ax[0].plot(x_axis, overall_tot)#*np.sum(subrange)/np.sum(overall_tot),'kx')
#ax[2].plot(x_axis, subrange - overall_tot*np.sum(subrange)/np.sum(overall_tot),'kx')
ax[2].plot(x_axis, z[0] - subrange/np.sum(subrange))
ax[2].plot(x_axis, overall_tot - subrange/np.sum(subrange),'x--')
fig.canvas.draw();fig.show()


from scipy.stats import kde
import matplotlib.pyplot as plt     

density = kde.gaussian_kde(subrange/np.sum(subrange))  # your data
#xgrid = np.linspace(x.min(), x.max(), 1024)   
fig, ax = pt.subplots()
ax.plot(x_axis, density(x_axis))
#pt.plot(x_axis, subrange)
fig.canvas.draw();fig.show()

1/0

b = clust.EM_GMM_GMM_clustering_class_self(dat, n_clusters = 3, n_iterations = 20, n_cpus=1, start='random', kappa_calc='approx', hard_assignments = 0, kappa_converged = 0.1, mu_converged = 0.01, min_iterations=10, LL_converged = 1.e-4, verbose = 0, seed=None, norm_method = 'sum')


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize


data = np.array([-69,3, -68, 2, -67, 1, -66, 1, -60, 1, -59, 1,
                 -58, 1, -57, 2, -56, 1, -55, 1, -54, 1, -52, 1,
                 -50, 2, -48, 3, -47, 1, -46, 3, -45, 1, -43, 1,
                 0, 1, 1, 2, 2, 12, 3, 18, 4, 18, 5, 13, 6, 9,
                 7, 7, 8, 5, 9, 3, 10, 1, 13, 2, 14, 3, 15, 2,
                 16, 2, 17, 2, 18, 2, 19, 2, 20, 2, 21, 3, 22, 1,
                 24, 1, 25, 1, 26, 1, 28, 2, 31, 1, 38, 1, 40, 2])
data = dat
x, y = data.reshape(-1, 2).T

def tri_norm(x, *args):
    m1, m2, m3, s1, s2, s3, k1, k2, k3 = args
    ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
    ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
    ret += k3*scipy.stats.norm.pdf(x, loc=m3 ,scale=s3)
    return ret


params = [-50, 3, 20, 1, 1, 1, 1, 1, 1]

fitted_params,_ = scipy.optimize.curve_fit(tri_norm,x, y, p0=params)

fig, ax = pt.subplots()
#ax.plot(x_axis, density(x_axis))
#pt.plot(x_axis, subrange)
ax.plot(x, y, 'o')
xx = np.linspace(np.min(x), np.max(x), 1000)
ax.plot(xx, tri_norm(xx, *fitted_params))
#plt.show()
fig.canvas.draw();fig.show()

##############
### MP MUS ###
##############

### Importing libraries
import numpy as np
import datetime
import matplotlib.pyplot as plt

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

# Importing my functions:
from functions import deproject, hyperbolic, cut, test_map, rotate, add_noise, test_map_mie, hg_map
from functions import e_plot, e_best, e_opt, e_evo, e_mask, e_score
from functions import a_plot, a_best, a_opt, a_surf_evo, a_surf_opt, a_gau_opt, a_gau_evo, a_mie_opt, a_hg_opt
from errors import e_best_err, e_evo_err, a_gau_opt_err


### COMPUTATION TIME START

t_start = datetime.datetime.now()


### IMPORTING DATA

# Import an rstokesdc:                                                                                                             
filename        = 'MPMus-J_S20160306S0198_combined_rstokesdc_phot.fits'
i,qphi,uphi,v   = getfitsdata(filename)

# Get useful keywords from fits header:                                                                                            
target  = getfitskeywords(filename, 'OBJECT')
itime   = getfitskeywords(filename, 'ITIME', HEADER='SCI')
# print('target, itime', target,itime)

qphi = np.nan_to_num(qphi)

vu = np.quantile(qphi, 0.99)
vl = np.quantile(qphi, 0.01)
# print("upper =",vu, " lower=",vl)


### HYPERBOLIC FUNCTION
# qphi = hyperbolic(qphi, 10, vu, vl)

### PLOTTING IMAGE AND BEST FIT ELLIPSE
# qphi = deproject(qphi, 31)

qphi = cut(35, 0, 0, 141, 141, 0, qphi)

plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

# t = test_map(55.37, 24.36, 23.32, 98.48, 140.2, 141.2, 9.954, 0, 282)
# t = test_map_mie(55, 24, 23, 97, 141, 142, 7, 5, 35, 0, 282)
# t = hg_map(50, 20, 45, 45, 141, 141, 0.5, 20, 0, 282)
# t = hg_map(54.41, 23.87, 18.19, -517.70, 140.49, 142.39, -0.2524, 149.09, -2.03275376e-06, 282)
# t = add_noise(t, 1)
# plt.imshow(t, cmap='seismic', origin='lower')

plt.colorbar(shrink=0.8)

### ELLIPSE FITTING
# e = e_best(55, 65, 30, 35, 90, 110, 141, 142, 141, 142, qphi)
# e = e_opt(60, 30, 110, 141, 141, qphi)
# e = e_evo(55, 65, 25, 40, 80, 120, 138, 144, 138, 144, qphi)

### ELLIPSE PLOTTING
# e = 50, 30, 0, 141, 141
# e_plot(e, 'k')


### ANNULUS FITTING
# a = a_best(58, 62, 18, 22, 28, 32, 90, 120, 141, 142, 141, 142, t)
# a = a_opt(55, 20, 30, 90, 141, 141, t)

# a = a_surf_evo(55, 60, 1, 30, 30, 35, 80, 90, 141, 143, 141, 143, 5, 15, 0, 2, qphi)
# a = a_surf_opt(50, 25, 30, 90, 141, 141, 10, 1, qphi)

# a = a_gau_opt(60, 20, 30, 90, 141, 141, 10, 1, qphi)
# a = a_gau_evo(55, 60, 10, 40, 30, 35, 80, 90, 141, 143, 141, 143, 5, 15, 0, 2, qphi)

# a = a_mie_opt(60, 20, 30, 90, 141, 141, 10, 10, 45, 0, qphi)

a = a_hg_opt(60, 20, 30, 10, 141, 141, 0.5, 10, 0, qphi)

### ANNULUS PLOTTING
a_plot(a[0], a[1], a[2], a[3], a[4], a[5], 'k', 0.4)
# a_plot(55.37, 24.36, 33.18, 98.03, 142.5, 141.3, 'k', 0.4)
# e_plot(a[0], a[2], a[3], a[4], a[5], 'k')

plt.show()


### ERROR FINDING
# limits = (50, 60), (20, 30), (25, 35), (90, 100), (135, 145), (135, 145), (5, 15), (-5, 5)

# errors = e_best_err(3, limits)
# errors = e_evo_err(3, limits)
# errors = a_gau_opt_err(2, limits)

# print(errors)


### COMPUTATION TIME END
t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)

# t_start = datetime.datetime.now()
# e = e_evo(55, 65, 25, 40, 80, 120, 138, 144, 138, 144, qphi)
# e_plot(e, 'k')
# t_end = datetime.datetime.now()
# print("Computation time:", t_end - t_start)


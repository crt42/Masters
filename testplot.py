### MP MUS

import numpy as np
import datetime
import matplotlib.pyplot as plt

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

# Importing my functions:
from functions import deproject, hyperbolic, cut, test_map, rotate
from functions import e_plot, e_best, e_opt, e_evo, e_mask, e_score
from functions import a_plot, a_best, a_opt, a_surf_evo, a_surf_opt, a_gau_opt, a_gau_evo

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

# t = test_map(57, 20, 30, 100, 141, 141, 10, 0, 282, False)
# plt.imshow(t, cmap='seismic', origin='lower')

plt.colorbar(shrink=0.8)

### ELLIPSE FITTING
# e = e_best(55, 65, 30, 35, 90, 110, 141, 142, 141, 142, qphi)
# e = e_opt(60, 30, 110, 141, 141, qphi)
e = e_evo(55, 65, 25, 40, 80, 120, 138, 144, 138, 144, qphi)

### ELLIPSE PLOTTING
# e = 60, 30, 100, 141, 141
e_plot(e, 'k')

### ANNULUS FITTING
# a = a_best(58, 62, 18, 22, 28, 32, 90, 120, 141, 142, 141, 142, t)
# a = a_opt(55, 20, 30, 90, 141, 141, t)

# a = a_surf_evo(55, 60, 1, 30, 30, 35, 80, 90, 141, 143, 141, 143, 5, 15, 0, 2, qphi)
# a = a_surf_opt(50, 25, 30, 90, 141, 141, 10, 1, qphi)

# a = a_gau_opt(60, 20, 30, 90, 141, 141, 10, 1, qphi)
# a = a_gau_evo(55, 60, 10, 40, 30, 35, 80, 90, 141, 143, 141, 143, 5, 15, 0, 2, qphi)


### ANNULUS PLOTTING
# a_plot(a[0], a[1], a[2], a[3], a[4], a[5], 'k', 0.4)
# a_plot(57, 20, 30, 100, 141, 141, 'k', 0.4)
# e_plot(a[0], a[2], a[3], a[4], a[5], 'k')

plt.show()

### COMPUTATION TIME END
t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)

# t_start = datetime.datetime.now()
# e = e_evo(55, 65, 25, 40, 80, 120, 138, 144, 138, 144, qphi)
# e_plot(e, 'k')
# t_end = datetime.datetime.now()
# print("Computation time:", t_end - t_start)


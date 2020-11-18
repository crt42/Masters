### MP MUS

import numpy as np
import datetime
import matplotlib.pyplot as plt

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

# Importing my functions:
from functions import deproject, hyperbolic, test_map
from functions import e_plot, e_best, e_opt, e_evo
from functions import a_plot, a_best, a_opt, a_surf_evo, a_surf_opt

### COMPUTATION TIME START

t_start = datetime.datetime.now()

### IMPORTING DATA

# Import an rstokesdc:                                                                                                             
filename        = 'MPMus-J_S20160306S0198_combined_rstokesdc_phot.fits'
i,qphi,uphi,v   = getfitsdata(filename)

# Get useful keywords from fits header:                                                                                            
target  = getfitskeywords(filename, 'OBJECT')
itime   = getfitskeywords(filename, 'ITIME', HEADER='SCI')
print('target, itime', target,itime)

qphi = np.nan_to_num(qphi)

vu = np.quantile(qphi, 0.99)
vl = np.quantile(qphi, 0.01)
print("upper =",vu, " lower=",vl)

### HYPERBOLIC FUNCTION
# qphi = hyperbolic(qphi, 10, vu, vl)

### PLOTTING IMAGE AND BEST FIT ELLIPSE
# qphi = deproject(qphi, 31)

plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

t = test_map(60, 20, 30, 110, 141, 141, 10, 1, 282, True)
# plt.imshow(t, cmap='seismic', origin='lower')

plt.colorbar(shrink=0.8)

### ELLIPSE FITTING
# e = e_best(60, 65, 25, 35, 85, 115, 135, 145, 135, 145, t)
# e = e_opt(40, 0, 90, 141, 141, t)
# e = e_evo(40, 100, 0, 60, 45, 135, 131, 151, 131, 151, t)

### ELLIPSE PLOTTING
# e_plot(e[0], e[1], e[2], e[3], e[4], 'k')
# e_plot(61, 30, 90, 141, 141, 'k')

### ANNULUS FITTING6
# a = a_best(58, 62, 18, 22, 28, 32, 90, 120, 141, 142, 141, 142, t)
# a = a_opt(55, 20, 30, 90, 141, 141, t)
# a = a_surf_evo(50, 100, 1, 30, 20, 50, 45, 135, 141, 142, 141, 142, 0, 15, 1, 2, t)
# a = a_surf_opt(40, 25, 30, 90, 141, 141, 10, 1, qphi)

### ANNULUS PLOTTING
a_plot(a[0], a[1], a[2], a[3], a[4], a[5], 'k', 0.4)
# a_plot(55, 20, 30, 70, 141, 141, 'k', 0.4)
# e_plot(a[0], a[2], a[3], a[4], a[5], 'k')

plt.show()

### COMPUTATION TIME END
t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)
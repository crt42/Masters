### MP MUS

import numpy as np
import datetime
import matplotlib.pyplot as plt

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

# Importing my functions:
from functions import deproject, hyperbolic
from functions import e_plot, e_best, e_opt, e_evo
from functions import a_plot, a_best, a_opt

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

vu = np.quantile(qphi, 0.9)
vl = np.quantile(qphi, 0.02)
print("upper =",vu, " lower=",vl)

### HYPERBOLIC FUNCTION

qphi = hyperbolic(qphi, 10, vu, vl)
    
### PLOTTING IMAGE AND BEST FIT ELLIPSE

# qphi = deproject(qphi, 31)

plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)
    
# e = e_best(50, 60, 30, 35, 85, 95, 140, 145, 140, 145, qphi)
# e = e_opt(60, 0, 90, 141, 141, qphi)
# e = e_evo(40, 100, 0, 60, 0, 180, 101, 181, 101, 181, qphi)

# e_plot(e[0], e[1], e[2], e[3], e[4], 'k')
# e_plot(61, 30, 90, 141, 141, 'k')

# a = a_best(58, 62, 1, 20, 30, 31, 90, 91, 141, 142, 141, 142, qphi)
# a = a_opt(55, 20, 30, 90, 141, 141, qphi)

# a_plot(a[0], a[1], a[2], a[3], a[4], a[5], 'k', 0.2)
# a_plot(55, 20, 30, 90, 141, 141, 'k', 0.1)
# e_plot(a[0], a[2], a[3], a[4], a[5], 'k')

plt.show()

### COMPUTATION TIME END
t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)
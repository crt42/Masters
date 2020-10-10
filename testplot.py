import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.convolution import convolve, Gaussian2DKernel

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

# Import an rstokesdc:                                                                                                             
filename        = 'MPMus-J_S20160306S0198_combined_rstokesdc_phot.fits'
i,qphi,uphi,v   = getfitsdata(filename)

# Get useful keywords from fits header:                                                                                            
target  = getfitskeywords(filename, 'OBJECT')
itime   = getfitskeywords(filename, 'ITIME', HEADER='SCI')
print('target, itime', target,itime)

qphi = np.nan_to_num(qphi)

# SMOOTH IMAGE
STD = 2.1 / (2.0*np.sqrt(2.0*np.log10(2.0)))
qphi_smooth = convolve(qphi, Gaussian2DKernel(STD))
qphi = qphi_smooth/qphi

vu = np.quantile(qphi, 0.99)
vl = np.quantile(qphi, 0.01)

print("upper =",vu, " lower=",vl)
plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', norm=LogNorm(vmin=0.01, vmax=100), origin='lower')

plt.show()

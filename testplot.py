import numpy as np
import matplotlib.pyplot as plt

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

vu = np.quantile(qphi, 0.9)
vl = np.quantile(qphi, 0.02)

# HYPERBOLIC FUNCTIONS
beta = 10
qphi = (np.arcsinh((qphi - vl)/beta))/(np.arcsinh((vu - vl)/beta))

print("upper =",vu, " lower=",vl)
plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

# ELLIPSE PLOTTING
x_m = 130
y_m = 130
r = 70
i = np.radians(32)
t_rot = np.radians(95)

t = np.linspace(0, 2*np.pi, 100)
Ell = np.array([r*np.cos(t), r*np.cos(i)*np.sin(t)])  
R_rot = np.array([[np.cos(t_rot), -np.sin(t_rot)],
                  [np.sin(t_rot), np.cos(t_rot)]])
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
    
plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], 'k')

plt.show()
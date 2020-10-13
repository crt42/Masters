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

# ELLIPSE PLOTTING
x_m = 145
y_m = 141
r = 70
i = np.radians(32)
t_rot = np.radians(95)

t = np.linspace(0, 2*np.pi, 400)
Ell = np.array([r*np.cos(t), r*np.cos(i)*np.sin(t)])  
R_rot = np.array([[np.cos(t_rot), -np.sin(t_rot)],
                  [np.sin(t_rot), np.cos(t_rot)]])
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

# MASK    

mask = np.full((len(qphi), len(qphi[0])), False)

for i in range(Ell.shape[1]):
    sq_x = x_m + np.int(Ell_rot[0,i])
    sq_y = y_m + np.int(Ell_rot[1,i])
    mask[sq_y, sq_x] = True

for i in range(len(qphi)):
    for j in range(len(qphi[0])):
        if (mask[i,j]):
            qphi[i,j] = 100;

# DATA PLOTTING

print("upper =",vu, " lower=",vl)
plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

#plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], 'k')

plt.show()
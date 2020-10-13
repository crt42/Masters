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
print("upper =",vu, " lower=",vl)

# HYPERBOLIC FUNCTIONS
beta = 10
qphi = (np.arcsinh((qphi - vl)/beta))/(np.arcsinh((vu - vl)/beta))

# ELLIPSE PLOTTING
x_m = 141
y_m = 141
# r = 50
# inc = np.radians(0)
t_rot = np.radians(95)

r_min = 40
r_max = 100
inc_min = 0
inc_max = 40
score_list = np.zeros((r_max - r_min, inc_max - inc_min))

for r in range(r_min, r_max):
    for inc in range(inc_min, inc_max):
    
        t = np.linspace(0, 2*np.pi, 400)
        Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
        R_rot = np.array([[np.cos(t_rot), -np.sin(t_rot)],
                          [np.sin(t_rot), np.cos(t_rot)]])
        Ell_rot = np.zeros((2,Ell.shape[1]))
        for i in range(Ell.shape[1]):
            Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
        
        # MASK    
        mask = np.full((len(qphi), len(qphi[0])), False)
        score = 0;
        
        # Producing the mask
        for i in range(Ell.shape[1]):
            sq_x = x_m + np.int(Ell_rot[0,i])
            sq_y = y_m + np.int(Ell_rot[1,i])
            mask[sq_y, sq_x] = True
        
        # Calculating the score
        for i in range(len(qphi)):
            for j in range(len(qphi[0])):
                if (mask[i,j]):
                    score += qphi[i, j]
                    
        score_list[r - r_min, inc - inc_min] = score
        
print(np.max(score_list))
max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
print(max_coords)

r = max_coords[0] + r_min
inc = max_coords[1] + inc_min

t = np.linspace(0, 2*np.pi, 400)
Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
R_rot = np.array([[np.cos(t_rot), -np.sin(t_rot)],
                  [np.sin(t_rot), np.cos(t_rot)]])
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

# PLOTTING

# Plotting the image
plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

# Plotting the ellipse
plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], 'k')

plt.show()
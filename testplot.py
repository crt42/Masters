### MP MUS

import numpy as np
import matplotlib.pyplot as plt
import datetime

# Custom functions from our pypeline:
from gpi_analysis.plot      import imshow_fancy, get_vlims, scale_colourbar
from gpi_analysis.inputs    import getfitsdata, getfitskeywords
from gpi_analysis.analysis  import make_radialstokes, make_linpolint

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

beta = 10
qphi = (np.arcsinh((qphi - vl)/beta))/(np.arcsinh((vu - vl)/beta))

### BEST FITTING ELLIPSE SEARCH

# Initialsing the centre coordinates

# Setting boundaries for r, inc, and t_rot
r_min = 55
r_max = 65

inc_min = 25
inc_max = 35

t_rot_min = 85
t_rot_max = 95

x_m_min = 138
x_m_max = 150

y_m_min = 138
y_m_max = 150

# Initialising the score array
score_list = np.zeros((r_max - r_min,
                       inc_max - inc_min,
                       t_rot_max - t_rot_min,
                       x_m_max - x_m_min,
                       y_m_max - y_m_min))

# Iterating over all values of r, inc, t_rot, x_m, and y_m
for r in range(r_min, r_max):
    # Printing time and progress
    print(datetime.datetime.now(), "--->", (r - r_min)/(r_max - r_min) * 100, "%")
    for inc in range(inc_min, inc_max):
        for t_rot in range(t_rot_min, t_rot_max):
            for x_m in range(x_m_min, x_m_max):
                for y_m in range(y_m_min, y_m_max):
    
                    # Initialising a mask and score    
                    mask = np.full((len(qphi), len(qphi[0])), False)
                    score = 0;
                    
                    # Drawning an ellipse with this r, inc, t_rot, x_m, and y_m
                    t = np.linspace(0, 2*np.pi, 400)
                    Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
                    R_rot = np.array([[np.cos(np.radians(t_rot)), -np.sin(np.radians(t_rot))],
                                      [np.sin(np.radians(t_rot)), np.cos(np.radians(t_rot))]])
                    Ell_rot = np.zeros((2,Ell.shape[1]))
                    for i in range(Ell.shape[1]):
                        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
                        
                        # Producing the mask and calculating score
                        sq_x = x_m + np.int(Ell_rot[0,i])
                        sq_y = y_m + np.int(Ell_rot[1,i])
                        if (mask[sq_y, sq_x] == False):
                            score += qphi[sq_y, sq_x]
                            mask[sq_y, sq_x] = True
                    
                    # Adding this ellipse's score to the score array
                    score_list[r - r_min,
                               inc - inc_min,
                               t_rot - t_rot_min,
                               x_m - x_m_min,
                               y_m - y_m_min] = score
  
# Printing the maximum score
print("Max score = ", np.max(score_list))

# Finding the maximum score coordinates in the score array
max_coords = np.unravel_index(score_list.argmax(), score_list.shape)

# Finding the maximum score ellipse's r, inc, t_rot, x_m, and y_m
r = max_coords[0] + r_min
inc = max_coords[1] + inc_min
t_rot = max_coords[2] + t_rot_min
x_m = max_coords[3] + x_m_min
y_m = max_coords[4] + y_m_min

print("Best radius = ", r)
print("Best inclination = ", inc)
print("Best rotation angle = ", t_rot)
print("Best centre coordinates = ", "(", x_m, ",", y_m, ")")

### BEST FITTING ELLIPSE DRAWING

t = np.linspace(0, 2*np.pi, 400)
Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
R_rot = np.array([[np.cos(np.radians(t_rot)), -np.sin(np.radians(t_rot))],
                  [np.sin(np.radians(t_rot)), np.cos(np.radians(t_rot))]])
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

### PLOTTING IMAGE AND ELLIPSE

# Plotting the image
plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)

# Plotting the ellipse
plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], 'k')

plt.show()

### COMPUTATION TIME END

t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)
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
def best_ellipse(r_min, r_max,
                 inc_min, inc_max,
                 t_rot_min, t_rot_max,
                 x_m_min, x_m_max,
                 y_m_min, y_m_max):
    
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
                        mask_no = 0;
                        score = 0;
                        
                        # Drawning an ellipse with this r, inc, t_rot, x_m, and y_m
                        t = np.linspace(0, 2*np.pi, 400)
                        Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
                        R_rot = np.array([[np.cos(np.radians(t_rot)), -np.sin(np.radians(t_rot))],
                                          [np.sin(np.radians(t_rot)), np.cos(np.radians(t_rot))]])
                        Ell_rot = np.zeros((2,Ell.shape[1]))
                        
                        for i in range(Ell.shape[1]):
                            
                            # Rotating the ellipse
                            Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
                            
                            # Producing the mask
                            sq_x = x_m + np.int(Ell_rot[0,i])
                            sq_y = y_m + np.int(Ell_rot[1,i])
                            
                            # Checking to not count a pixel's score multiple times
                            if (mask[sq_y, sq_x] == False):
                                
                                # Calculating the score
                                score += qphi[sq_y, sq_x]
                                mask_no += 1
                                mask[sq_y, sq_x] = True
                        
                        # Adding this ellipse's score to the score array
                        score_list[r - r_min,
                                   inc - inc_min,
                                   t_rot - t_rot_min,
                                   x_m - x_m_min,
                                   y_m - y_m_min] = score/mask_no
      
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, t_rot, x_m, and y_m
    
    ell_params = np.empty(5)   
    ell_params[0] = max_coords[0] + r_min
    ell_params[1] = max_coords[1] + inc_min
    ell_params[2] = max_coords[2] + t_rot_min
    ell_params[3] = max_coords[3] + x_m_min
    ell_params[4] = max_coords[4] + y_m_min
    
    print("Best radius = ", ell_params[0])
    print("Best inclination = ", ell_params[1])
    print("Best rotation angle = ", ell_params[2])
    print("Best centre coordinates = ", "(", ell_params[3], ",", ell_params[4], ")")
    
    return ell_params

### ELLIPSE DRAWING
def plot_ellipse(r, inc, t_rot, x_m, y_m, col):

    t = np.linspace(0, 2*np.pi, 100)
    
    Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
    R_rot = np.array([[np.cos(np.radians(t_rot)), -np.sin(np.radians(t_rot))],
                      [np.sin(np.radians(t_rot)), np.cos(np.radians(t_rot))]])
    Ell_rot = np.zeros((2,Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
    
    # Plotting the ellipse
    plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], col)

### BEST FITTING ANNULUS SEARCH
def best_annulus(r_min, r_max,
                 th_min, th_max,
                 inc_min, inc_max,
                 t_rot_min, t_rot_max,
                 x_m_min, x_m_max,
                 y_m_min, y_m_max):
    
    # Initialising the score array
    score_list = np.zeros((r_max - r_min,
                           th_max - th_min,
                           inc_max - inc_min,
                           t_rot_max - t_rot_min,
                           x_m_max - x_m_min,
                           y_m_max - y_m_min))
    
    # Iterating over all values of r, th, inc, t_rot, x_m, and y_m
    for r in range(r_min, r_max):
        # Printing time and progress
        print(datetime.datetime.now(), "--->", (r - r_min)/(r_max - r_min) * 100, "%")
        for th in range(th_min, th_max):
            for inc in range(inc_min, inc_max):
                for t_rot in range(t_rot_min, t_rot_max):
                    for x_m in range(x_m_min, x_m_max):
                        for y_m in range(y_m_min, y_m_max):
            
                            # Initialising a mask and score    
                            mask = np.full((len(qphi), len(qphi[0])), False)
                            mask_no = 0;
                            score = 0;
                            
                            for i in range(0, 282):
                                for j in range(0, 282):
                                    z = (i*np.cos(np.radians(t_rot)) - j*np.sin(np.radians(t_rot)) - x_m*np.cos(np.radians(t_rot)) + y_m*np.sin(np.radians(t_rot)))**2 + ((j*np.cos(np.radians(t_rot)) + i*np.sin(np.radians(t_rot)) - y_m*np.cos(np.radians(t_rot)) - x_m*np.sin(np.radians(t_rot)))/np.cos(np.radians(inc)))**2
                                    if (z < r - th/2 and z < r + th/2):
                                        if (mask[i,j] == False):
                                            # Calculating the score
                                            score += qphi[i, j]
                                            mask_no += 1
                                            mask[i,j] = True

                            # Adding this ellipse's score to the score array
                            score_list[r - r_min,
                                       th - th_min,
                                       inc - inc_min,
                                       t_rot - t_rot_min,
                                       x_m - x_m_min,
                                       y_m - y_m_min] = score/mask_no
      
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, t_rot, x_m, and y_m
    
    a_params = max_coords[0] + r_min, max_coords[1] + th_min, max_coords[2] + inc_min, max_coords[3] + t_rot_min, max_coords[4] + x_m_min, max_coords[5] + y_m_min
    
    return a_params
    

### ANNULUS DRAWING
def plot_annulus(r, th, inc, t_rot, x_m, y_m, col, a):
        
    # Plotting the annulus
    
    x_l = np.linspace(0, 282, 100)
    y_l = np.linspace(0, 282, 100)
    
    x, y = np.meshgrid(x_l,y_l)

    z = (x*np.cos(np.radians(t_rot)) - y*np.sin(np.radians(t_rot)) - x_m*np.cos(np.radians(t_rot)) + y_m*np.sin(np.radians(t_rot)))**2 + ((y*np.cos(np.radians(t_rot)) + x*np.sin(np.radians(t_rot)) - y_m*np.cos(np.radians(t_rot)) - x_m*np.sin(np.radians(t_rot)))/np.cos(np.radians(inc)))**2
    plt.contourf(x, y, z, levels=[(r - th/2)**2, (r + th/2)**2],
                 colors = col, alpha = a)
    
### PLOTTING IMAGE AND BEST FIT ELLIPSE

plt.figure(figsize=(12,12))
plt.imshow(qphi, cmap='seismic', origin='lower', vmin = vl, vmax = vu)
    
# e = best_ellipse(50, 60, 25, 35, 82, 92, 140, 145, 140, 145)
# plot_ellipse(e[0], e[1], e[2], e[3], e[4], 'k')

plot_ellipse(60, 32, 90, 141, 141, 'k')

a = best_annulus(50, 60, 10, 20, 30, 31, 90, 91, 141, 142, 141, 142)

plot_annulus(a[0], a[1], a[2], a[3], a[4], a[5], 'k', 0.2)

plt.show()

### COMPUTATION TIME END
t_end = datetime.datetime.now()
print("Computation time:", t_end - t_start)
### FUNCTIONS FOR PROTOPLANETARY DISC ANALYSIS

import numpy as np
import datetime
import scipy.optimize
import matplotlib.pyplot as plt

### DEPROJECTION
def deproject(data, inc):
    w = len(data[0])
    h = len(data)
    inc = np.radians(inc)
    
    new_w = np.int(w/np.cos(inc))
    new_data = np.zeros((h, new_w))
    print(new_w)
    for i in range(0, h):
        for j in range(0, new_w):

            new_data[i, j] = data[i, np.int(j*np.cos(inc))]
            
    return new_data

### ELLIPSE FUNCTIONS

### BEST FITTING ELLIPSE SEARCH
def best_ellipse(r_min, r_max, inc_min, inc_max, t_rot_min, t_rot_max, x_m_min, x_m_max, y_m_min, y_m_max, data):
    
    # Initialising the score array
    score_list = np.zeros((r_max - r_min, inc_max - inc_min,
                           t_rot_max - t_rot_min, x_m_max - x_m_min,
                           y_m_max - y_m_min))
    
    # Iterating over all values of r, inc, t_rot, x_m, and y_m
    for r in range(r_min, r_max):
        # Printing time and progress
        print(datetime.datetime.now(), "--->", (r - r_min)/(r_max - r_min) * 100, "%")
        for inc in range(inc_min, inc_max):
            for t_rot in range(t_rot_min, t_rot_max):
                for x_m in range(x_m_min, x_m_max):
                    for y_m in range(y_m_min, y_m_max):
                        # Adding this ellipse's score to the score array
                        this_score = score_ellipse(r, inc, t_rot, x_m, y_m, data)
                        score_list[r - r_min, inc - inc_min,
                                    t_rot - t_rot_min, x_m - x_m_min,
                                    y_m - y_m_min] = this_score
    
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, t_rot, x_m, and y_m
    
    ell_params = max_coords[0] + r_min, max_coords[1] + inc_min, max_coords[2] + t_rot_min, max_coords[3] + x_m_min, max_coords[4] + y_m_min
    
    print("Best radius = ", ell_params[0])
    print("Best inclination = ", ell_params[1])
    print("Best rotation angle = ", ell_params[2])
    print("Best centre coordinates = ", "(", ell_params[3],
                                        ",", ell_params[4], ")")
    
    return ell_params  

### SCORE ELLIPSE
def score_ellipse(r, inc, t_rot, x_m, y_m, data):
    
    # Initialising a mask and score    
    mask = np.full((len(data), len(data[0])), False)
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
        sq_x = np.int(x_m) + np.int(Ell_rot[0,i])
        sq_y = np.int(y_m) + np.int(Ell_rot[1,i])
                            
        # Checking to not count a pixel's score multiple times
        if (mask[sq_y, sq_x] == False):
                                
            # Calculating the score
            score += data[sq_y, sq_x]
            mask_no += 1
            mask[sq_y, sq_x] = True
            
    return score

### RECIPROCAL SCORE ELLIPSE
def recip_score_ellipse(init, data):
    r, inc, t_rot, x_m, y_m = init
    # Finds the score of that ellipse
    score = score_ellipse(r, inc, t_rot, x_m, y_m, data)
    # Returns the reciprocal of the score
    return 1/score

### OPTIMISED ELLIPSE
def opt_ellipse(r, inc, t_rot, x_m, y_m, data):
    init = (r, inc, t_rot, x_m, y_m)
    # Finds the minimum reciprocal scored ellipse
    params = scipy.optimize.minimize(recip_score_ellipse, init,
                                     args = data, method = 'Powell')
    print(params['x'])
    return params['x']

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
    
### ANNULUS FUNCTIONS

### BEST FITTING ANNULUS SEARCH
def best_annulus(r_min, r_max, th_min, th_max, inc_min, inc_max, t_rot_min, t_rot_max, x_m_min, x_m_max, y_m_min, y_m_max, data):
    
    # Initialising the score array
    score_list = np.zeros((r_max - r_min,
                           th_max - th_min,
                           inc_max - inc_min,
                           t_rot_max - t_rot_min,
                           x_m_max - x_m_min,
                           y_m_max - y_m_min))
    
    # Iterating over all values of r, th, inc, t_rot, x_m, and y_m
    for r in range(r_min, r_max):
        for th in range(th_min, th_max):
            # Printing time and progress
            print(datetime.datetime.now(), "--->", ((r - r_min + (th - th_min)/(th_max - th_min)))/(r_max - r_min) * 100, "%")
            for inc in range(inc_min, inc_max):
                for t_rot in range(t_rot_min, t_rot_max):
                    for x_m in range(x_m_min, x_m_max):
                        for y_m in range(y_m_min, y_m_max):
                            
                            this_score = score_annulus(r, th, inc, t_rot, x_m, y_m, data)
                            score_list[r - r_min, inc - inc_min, th - th_min,
                                       t_rot - t_rot_min, x_m - x_m_min,
                                       y_m - y_m_min] = this_score
    
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, t_rot, x_m, and y_m
    
    a_params = max_coords[0] + r_min, max_coords[1] + th_min, max_coords[2] + inc_min, max_coords[3] + t_rot_min, max_coords[4] + x_m_min, max_coords[5] + y_m_min
    
    print("Best radius = ", a_params[0])
    print("Best thickness = ", a_params[1])
    print("Best inclination = ", a_params[2])
    print("Best rotation angle = ", a_params[3])
    print("Best centre coordinates = ", "(", a_params[4],
                                        ",", a_params[5], ")")
    
    return a_params

### SCORE ANNULUS
def score_annulus(r, th, inc, t_rot, x_m, y_m, data):
    # Initialising a mask and score    
    mask = np.full((len(data), len(data[0])), False)
    mask_no = 0
    score = 0
                            
    for i in range(0, 282):
        for j in range(0, 282):
            z = (i*np.cos(np.radians(t_rot)) - j*np.sin(np.radians(t_rot)) - x_m*np.cos(np.radians(t_rot)) + y_m*np.sin(np.radians(t_rot)))**2 + ((j*np.cos(np.radians(t_rot)) + i*np.sin(np.radians(t_rot)) - y_m*np.cos(np.radians(t_rot)) - x_m*np.sin(np.radians(t_rot)))/np.cos(np.radians(inc)))**2
            if (z < r - th/2 and z < r + th/2):
                if (mask[i,j] == False):
                    # Calculating the score
                    score +=  data[i, j]
                    mask_no += 1
                    mask[i,j] = True
    if mask_no == 0:
        mask_no = 1
    return score/mask_no

### RECIPROCAL SCORE ANNULUS
def recip_score_annulus(init, data):
    r, th, inc, t_rot, x_m, y_m = init
    # Finds the score of that ellipse
    score = score_annulus(r, th, inc, t_rot, x_m, y_m, data)
    # Returns the reciprocal of the score
    return 1/score    

### OPTIMISED ANNULUS
def opt_annulus(r, th, inc, t_rot, x_m, y_m, data):
    init = (r, th, inc, t_rot, x_m, y_m)
    # Finds the minimum reciprocal scored ellipse
    params = scipy.optimize.minimize(recip_score_annulus, init,
                                     args = data, method = 'Powell')
    print(params['x'])
    return params['x']

### ANNULUS DRAWING
def plot_annulus(r, th, inc, t_rot, x_m, y_m, col, a):
        
    # Plotting the annulus
    
    x_l = np.linspace(0, 282, 100)
    y_l = np.linspace(0, 282, 100)
    
    x, y = np.meshgrid(x_l,y_l)

    z = (x*np.cos(np.radians(t_rot)) - y*np.sin(np.radians(t_rot)) - x_m*np.cos(np.radians(t_rot)) + y_m*np.sin(np.radians(t_rot)))**2 + ((y*np.cos(np.radians(t_rot)) + x*np.sin(np.radians(t_rot)) - y_m*np.cos(np.radians(t_rot)) - x_m*np.sin(np.radians(t_rot)))/np.cos(np.radians(inc)))**2
    plt.contourf(x, y, z, levels=[(r - th/2)**2, (r + th/2)**2],
                 colors = col, alpha = a)
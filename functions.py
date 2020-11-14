### FUNCTIONS FOR PROTOPLANETARY DISC ANALYSIS

import numpy as np
import datetime
import scipy.optimize
import matplotlib.pyplot as plt



#########################
### GENERAL FUNCTIONS ###
#########################

### HYPERBOLIC SCALING
### Takes a set of data, beta value, and upper and lower limits, and scales
### the data with a hyperbolic function.
def hyperbolic(data, beta, vu, vl):
    return (np.arcsinh((data - vl)/beta))/(np.arcsinh((vu - vl)/beta))  

### DEPROJECTION
### Takes a set of data, and stretches this depending on the inclination,
### resulting in it being deprojected.
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

### TEST MAP
### Creates a map of an artificial ring using input parameters, in order to
### test fitting methods on.
def test_map(r, th, inc, rot, x_m, y_m, surf, back, size):
    test_map = np.full((size, size), back)
    for i in range(size):
        for j in range(size):
            R = np.sqrt(pow(i - x_m, 2) + pow(j - y_m, 2))
            test_map[i,j] = back + surf*np.exp((-pow(R - r, 2))/(0.5*pow(th, 2)))
    return test_map



#########################
### ELLIPSE FUNCTIONS ###
#########################

### BEST FITTING ELLIPSE SEARCH
### Searches through all ellipses in the given ranges for parameters,
### returning the best scoring one.
def e_best(r_min, r_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, data):
    
    # Initialising the score array
    score_list = np.zeros((r_max - r_min, inc_max - inc_min,
                           rot_max - rot_min, x_m_max - x_m_min,
                           y_m_max - y_m_min))
    
    # Iterating over all values of r, inc, rot, x_m, and y_m
    for r in range(r_min, r_max):
        # Printing time and progress
        print(datetime.datetime.now(), "--->", (r - r_min)/(r_max - r_min) * 100, "%")
        for inc in range(inc_min, inc_max):
            for rot in range(rot_min, rot_max):
                for x_m in range(x_m_min, x_m_max):
                    for y_m in range(y_m_min, y_m_max):
                        # Adding this ellipse's score to the score array
                        this_score = e_score(r, inc, rot, x_m, y_m, data)
                        score_list[r - r_min, inc - inc_min,
                                    rot - rot_min, x_m - x_m_min,
                                    y_m - y_m_min] = this_score
    
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, rot, x_m, and y_m
    
    ell_params = max_coords[0] + r_min, max_coords[1] + inc_min, max_coords[2] + rot_min, max_coords[3] + x_m_min, max_coords[4] + y_m_min
    
    print("Best radius = ", ell_params[0])
    print("Best inclination = ", ell_params[1])
    print("Best rotation angle = ", ell_params[2])
    print("Best centre coordinates = ", "(", ell_params[3],
                                        ",", ell_params[4], ")")
    
    return ell_params  

### SCORE ELLIPSE
### Scores a particular ellipse on the data given.
def e_score(r, inc, rot, x_m, y_m, data):
    
    # Initialising a mask and score    
    mask = np.full((len(data), len(data[0])), False)
    mask_no = 0;
    score = 0;
                        
    # Drawning an ellipse with this r, inc, rot, x_m, and y_m
    t = np.linspace(0, 2*np.pi, 400)
    Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
    R_rot = np.array([[np.cos(np.radians(rot)), -np.sin(np.radians(rot))],
                      [np.sin(np.radians(rot)), np.cos(np.radians(rot))]])
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
### Finds the reciprocal of the score of an ellipse.
def e_r_score(init, data):
    r, inc, rot, x_m, y_m = init
    # Finds the score of that ellipse
    score = e_score(r, inc, rot, x_m, y_m, data)
    # Returns the reciprocal of the score
    return 1/score

### OPTIMISED ELLIPSE
### Performs an optimised search for the best fitting ellipse, using the
### parameters as a starting point.
def e_opt(r, inc, rot, x_m, y_m, data):
    init = (r, inc, rot, x_m, y_m)
    # Finds the minimum reciprocal scored ellipse
    params = scipy.optimize.minimize(e_r_score, init,
                                     args = data, method = 'Powell')
    print(params['x'])
    return params['x']

### DIFFERENTIAL EVOLUTION ELLIPSE
### Performs a differential evolution algorithm search for the best
### fitting ellipse in the range of parameters given.
def e_evo(r_min, r_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, data):
    init = [(r_min, r_max), (inc_min, inc_max), (rot_min, rot_max),
            (x_m_min, x_m_max), (y_m_min, y_m_max)]   
    params = scipy.optimize.differential_evolution(e_r_score, init, args = (data,))
    print(params['x'])
    return params['x'] 

### ELLIPSE DRAWING
### Draws an ellipse with these parameters and colour.
def e_plot(r, inc, rot, x_m, y_m, col):

    t = np.linspace(0, 2*np.pi, 100)
    
    Ell = np.array([r*np.cos(t), r*np.cos(np.radians(inc))*np.sin(t)])  
    R_rot = np.array([[np.cos(np.radians(rot)), -np.sin(np.radians(rot))],
                      [np.sin(np.radians(rot)), np.cos(np.radians(rot))]])
    Ell_rot = np.zeros((2,Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
    
    # Plotting the ellipse
    plt.plot(x_m + Ell_rot[0,:], y_m + Ell_rot[1,:], col)



#########################
### ANNULUS FUNCTIONS ###
#########################

### ANNULUS Z FUNCTION
### Returns the radius of the ellipse at that point
def a_z(i, j, inc, rot, x_m, y_m):
    z = np.sqrt((i*np.cos(np.radians(rot)) - j*np.sin(np.radians(rot)) - x_m*np.cos(np.radians(rot)) + y_m*np.sin(np.radians(rot)))**2 + ((j*np.cos(np.radians(rot)) + i*np.sin(np.radians(rot)) - y_m*np.cos(np.radians(rot)) - x_m*np.sin(np.radians(rot)))/np.cos(np.radians(inc)))**2)
    return z

### BEST FITTING ANNULUS SEARCH
### Searches through all annuli in the given ranges for parameters,
### returning the best scoring one.
def a_best(r_min, r_max, th_min, th_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, data):
    
    # Initialising the score array
    score_list = np.zeros((r_max - r_min, th_max - th_min, inc_max - inc_min, rot_max - rot_min, x_m_max - x_m_min, y_m_max - y_m_min))
    
    # Iterating over all values of r, th, inc, rot, x_m, and y_m
    for r in range(r_min, r_max):
        for th in range(th_min, th_max):
            # Printing time and progress
            print(datetime.datetime.now(), "--->", ((r - r_min + (th - th_min)/(th_max - th_min)))/(r_max - r_min) * 100, "%")
            for inc in range(inc_min, inc_max):
                for rot in range(rot_min, rot_max):
                    for x_m in range(x_m_min, x_m_max):
                        for y_m in range(y_m_min, y_m_max):
                            
                            this_score = a_score(r, th, inc, rot, x_m, y_m, data)
                            score_list[r - r_min, th - th_min, inc - inc_min, rot - rot_min, x_m - x_m_min, y_m - y_m_min] = this_score
    
    # Printing the maximum score
    print("Max score = ", np.max(score_list))
    
    # Finding the maximum score coordinates in the score array
    max_coords = np.unravel_index(score_list.argmax(), score_list.shape)
    
    # Finding the maximum score ellipse's r, inc, rot, x_m, and y_m
    
    a_params = max_coords[0] + r_min, max_coords[1] + th_min, max_coords[2] + inc_min, max_coords[3] + rot_min, max_coords[4] + x_m_min, max_coords[5] + y_m_min
    
    print("Best radius = ", a_params[0])
    print("Best thickness = ", a_params[1])
    print("Best inclination = ", a_params[2])
    print("Best rotation angle = ", a_params[3])
    print("Best centre coordinates = ", "(", a_params[4],
                                        ",", a_params[5], ")")
    
    return a_params

### SCORE ANNULUS
### Scores a particular annulus on the data given.
def a_score(r, th, inc, rot, x_m, y_m, data):
    # Initialising a mask and score    
    mask = np.full((len(data), len(data[0])), False)
    mask_no = 0
    score = 0
    if th < 1:
        th = 1
                            
    for i in range(0, 282):
        for j in range(0, 282):
            z = a_z(i, j, inc, rot, x_m, y_m)
            if (z > r - th/2 and z < r + th/2):
                if (mask[i,j] == False):
                    # Calculating the score
                    score +=  data[i, j]
                    mask_no += 1
                    mask[i,j] = True
    if mask_no == 0:
        mask_no = 1
    return score/mask_no

### RECIPROCAL SCORE ANNULUS
### Finds the reciprocal of the score of an annulus.
def a_r_score(init, data):
    r, th, inc, rot, x_m, y_m = init
    # Finds the score of that ellipse
    score = a_score(r, th, inc, rot, x_m, y_m, data)
    # Returns the reciprocal of the score
    if score == 0:
        score = 0.001
    return 1/score    

### OPTIMISED ANNULUS
### Performs an optimised search for the best fitting annulus, using the
### parameters as a starting point.
def a_opt(r, th, inc, rot, x_m, y_m, data):
    init = (r, th, inc, rot, x_m, y_m)
    # Finds the minimum reciprocal scored annulus
    params = scipy.optimize.minimize(a_r_score, init,
                                     args = data, method = 'Powell')
    print(params['x'])
    return params['x']

### SURFACE BRIGHTNESS ANNULUS MAP
### Returns a 2D array plotting an annulus of given parameters.
def a_surf_map(r, th, inc, rot, x_m, y_m, surf, back, data):
    map = np.full((len(data), len(data[0])), back)
    for i in range (0, len(data)):
        for j in range(0, len(data[0])):
            z = a_z(i, j, inc, rot, x_m, y_m)
            if (z > r - th/2 and z < r + th/2):
                map[j, i] = surf
    return map

### SURFACE BRIGHTNESS ANNULUS SCORE
### Returns the score of the data minus a surface brightness annulus map.
def a_surf_score(init, data):
    r, th, inc, rot, x_m, y_m, surf, back = init
    score = 0
    map = a_surf_map(r, th, inc, rot, x_m, y_m, surf, back, data)
    for i in range (0, len(data)):
        for j in range(0, len(data[0])):
            score += (data[i,j] - map[i,j])**2
    return score

###
def a_surf_opt(r, th, inc, rot, x_m, y_m, surf, back, data):
    init = (r, th, inc, rot, x_m, y_m, surf, back)
    # Finds the minimum difference between data and annulus
    params = scipy.optimize.minimize(a_surf_score, init,
                                     args = data, method = 'Powell')
    print(params['x'])
    return params['x']

### DIFFERENTIAL EVOLUTION SURFACE BRIGHTNESS ANNULUS
### Performs a differential evolution algorithm search for the best
### fitting annulus with a given surface brightness, in the range of
### parameters given.
def a_surf_evo(r_min, r_max, th_min, th_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, surf_min, surf_max, back_min, back_max, data):
    init = [(r_min, r_max), (th_min, th_max), (inc_min, inc_max), (rot_min, rot_max), (x_m_min, x_m_max), (y_m_min, y_m_max), (surf_min, surf_max), (back_min, back_max)]   
    
    params = scipy.optimize.differential_evolution(a_surf_score, init, args = (data,), popsize=100, mutation=(1, 1.5))
    print(params['x'])
    return params['x'] 

### ANNULUS DRAWING
### Draws an annulus with these parameters and colour.
def a_plot(r, th, inc, rot, x_m, y_m, col, a):
    
    # Plotting the annulus
    x_l = np.linspace(0, 282, 100)
    y_l = np.linspace(0, 282, 100)
    x, y = np.meshgrid(x_l,y_l)
    
    z = a_z(x, y, inc, rot, x_m, y_m)
    plt.contourf(x, y, z, levels=[r - th/2, r + th/2], colors = col, alpha = a)
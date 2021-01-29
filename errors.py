###########################################
### ERROR FINDING FOR FITTING FUNCTIONS ###
###########################################

### Importing libraries
import numpy as np

### Importing my functions:
from functions import test_map, add_noise, e_best, e_evo, a_gau_opt

### ERROR FINDING FOR E_BEST
### This function finds the errors on the fitting function e_best
def e_best_err(no, limits):
     
    # Limits for the 8 parameters
    (r_min, r_max), (th_min, th_max), (inc_min, inc_max), (rot_min, rot_max), (x_m_min, x_m_max), (y_m_min, y_m_max), (surf_min, surf_max), (back_min, back_max) = limits
    
    # Creating an array to store differences in values
    result = np.empty([no, 5])
    
    for i in range(no):

        # Setting the actual values
        r = np.random.randint(r_min, r_max)
        th = np.random.randint(th_min, th_max)
        inc = np.random.randint(inc_min, inc_max)
        rot = np.random.randint(rot_min, rot_max)
        x_m = np.random.randint(x_m_min, x_m_max)
        y_m = np.random.randint(y_m_min, y_m_max)
        surf = np.random.randint(surf_min, surf_max)
        back = np.random.randint(back_min, back_max)
        noise = np.random.randint(1, 5)
        
        print("r, th, inc, rot, x_m, y_m, surf, back, noise")
        print(r, th, inc, rot, x_m, y_m, surf, back, noise)
        
        # Creating the test map
        t = test_map(r, th, inc, rot, x_m, y_m, surf, back, 282)
        t = add_noise(t, noise)        
        
        # Fitting to this test map
        e = e_best(r_min, r_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, t)
        
        print(e)
        
        result[i, 0] = e[0] - r 
        result[i, 1] = e[1] - inc
        result[i, 2] = e[2] - rot
        result[i, 3] = e[3] - x_m
        result[i, 4] = e[4] - y_m
        
    return result

### ERROR FINDING FOR E_EVO
### This function finds the errors on the fitting function e_evo
def e_evo_err(no, limits):
     
    # Limits for the 8 parameters
    (r_min, r_max), (th_min, th_max), (inc_min, inc_max), (rot_min, rot_max), (x_m_min, x_m_max), (y_m_min, y_m_max), (surf_min, surf_max), (back_min, back_max) = limits
    
    # Creating an array to store differences in values
    result = np.empty([no, 5])
    
    for i in range(no):

        # Setting the actual values for the test map
        r = np.random.randint(r_min, r_max)
        th = np.random.randint(th_min, th_max)
        inc = np.random.randint(inc_min, inc_max)
        rot = np.random.randint(rot_min, rot_max)
        x_m = np.random.randint(x_m_min, x_m_max)
        y_m = np.random.randint(y_m_min, y_m_max)
        surf = np.random.randint(surf_min, surf_max)
        back = np.random.randint(back_min, back_max)
        noise = np.random.randint(1, 5)
        
        print("r, th, inc, rot, x_m, y_m, surf, back, noise")
        print(r, th, inc, rot, x_m, y_m, surf, back, noise)
        
        # Creating the test map
        t = test_map(r, th, inc, rot, x_m, y_m, surf, back, 282)
        t = add_noise(t, noise)        
        
        # Fitting an ellipse to this test map
        e = e_evo(r_min, r_max, inc_min, inc_max, rot_min, rot_max, x_m_min, x_m_max, y_m_min, y_m_max, t)
        
        print(e)
        
        # Adding the results to the array
        result[i, 0] = e[0] - r
        result[i, 1] = e[1] - inc
        result[i, 2] = e[2] - rot
        result[i, 3] = e[3] - x_m
        result[i, 4] = e[4] - y_m
        
    return result

### ERROR FINDING FOR A_GAU_OPT
### This function finds the errors on the fitting function a_gau_opt       
def a_gau_opt_err(no, limits):
     
    # Limits for the 8 parameters
    (r_min, r_max), (th_min, th_max), (inc_min, inc_max), (rot_min, rot_max), (x_m_min, x_m_max), (y_m_min, y_m_max), (surf_min, surf_max), (back_min, back_max) = limits
    
    # Creating an array to store differences in values
    result = np.empty([no, 6])
    
    for i in range(no):

        # Setting the actual values for the test map
        r = np.random.randint(r_min, r_max)
        th = np.random.randint(th_min, th_max)
        inc = np.random.randint(inc_min, inc_max)
        rot = np.random.randint(rot_min, rot_max)
        x_m = np.random.randint(x_m_min, x_m_max)
        y_m = np.random.randint(y_m_min, y_m_max)
        surf = np.random.randint(surf_min, surf_max)
        back = np.random.randint(back_min, back_max)
        noise = np.random.randint(1, 5)
        
        print("r, th, inc, rot, x_m, y_m, surf, back, noise")
        print(r, th, inc, rot, x_m, y_m, surf, back, noise)
        
        # Setting the initial values for fitting function
        r_n = np.random.randint(r_min, r_max)
        th_n = np.random.randint(th_min, th_max)
        inc_n = np.random.randint(inc_min, inc_max)
        rot_n = np.random.randint(rot_min, rot_max)
        x_m_n = np.random.randint(x_m_min, x_m_max)
        y_m_n = np.random.randint(y_m_min, y_m_max)
        surf_n = np.random.randint(surf_min, surf_max)
        back_n = np.random.randint(back_min, back_max)
        
        # Creating the test map
        t = test_map(r, th, inc, rot, x_m, y_m, surf, back, 282)
        t = add_noise(t, noise)        
        
        # Fitting an ellipse to this test map
        a = a_gau_opt(r_n, th_n, inc_n, rot_n, x_m_n, y_m_n, surf_n, back_n, t)
        
        print(a)
        
        # Adding the results to the array
        result[i, 0] = a[0] - r
        result[i, 1] = a[1] - th
        result[i, 2] = a[2] - inc
        result[i, 3] = a[3] - rot
        result[i, 4] = a[4] - x_m
        result[i, 5] = a[5] - y_m
        
    return result
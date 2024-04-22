# initials.py
# version: 1.0.1 (Apr 22, 2024)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com 
#
import OmpApd as apd
import numpy as np

# Note: apd.phi_m(x) is the positive equilibrium value of the order parameter at a temperature of x [K].

def noise():
    
    if np.random.rand() >= 0.5:
        return 1
        
    else:
        return -1
        

initial_PF_2d = np.random.normal(loc = 0, scale = 0.0001, size = (60, 60))


for i in np.arange(0,60,10):
    initial_PF_2d[55:60,i+1:i+9] = noise() * apd.phi_m(673)
   
#initial_PF_2d = np.load("pfsnap2d_100.npy")    
   
#initial_PF_2d = np.loadtxt("pfsnap2d_100.txt")  
    
#####


initial_PF_3d = np.random.normal(loc = 0, scale = 0.0001, size = (60, 60, 60))


for i in np.arange(0,60,10):
    for j in np.arange(0,60,10):
    
        initial_PF_3d[55:60, i+1:i+9, j+1:j+9] = noise() * apd.phi_m(673)


#initial_PF_3d = np.load("pfsnap3d_100.npy")
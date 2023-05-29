# params.py
# version: 1.0.0 (May 29, 2023)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com
#
import OmpApd as apd
import numpy as np
import initials
import importlib

importlib.reload(initials)

CalcParams = {
    
            "omp_v": (6.040 + 6.619)/2*10**(-5), # molar volume of omp in m^3. Holland & Powell (2011) average of jd and di
            
            "nx": 400,  # numbers of the grid. ex) 2D: 400; 3D (normal): 60; 3d (cross-section mode): 100
            "ny": 400,
            "nz": 60, # only used for 3D (normal/cross-section mode) calculation
            "dx": 5.0e-9, # grid size in m. ex) 5.0e-9
            "dy": 5.0e-9,
            "dz": 5.0e-9, # only used for 3D (normal/cross-section mode) calculation
            
            "temp": 823, # temperature in K
            
            "aaa": 0, # gradient coef in m(J/mol)^0.5. If zero, aaa will be calculated based on intE and molar volume.
            "intE": 0.50, # APB energy in J/m2
            
            "ini": 0, # -1: user setting  for the initial field, 0: initial turbulence (Gaussian), 1: uniform distribution
            "turb": 0.0001, # standard deviation (if ini = 0) or the range (if ini = 1) of the initial turbulence
            "inifield": initials.initial_PF_2d, # initial phase field (valid when ini = -1)
            
            "periodi": 1, # 1/0 when i boundary condition is periodic/no-flow
            "periodj": 1, # 1/0 when j boundary condition is periodic/no-flow
            "periodk": 1, # 1/0 when k boundary condition is periodic/no-flow
            
            "dt": 0.01, # step size in normalized time (t* = mobility * rrr * temp * time(yr)). fixed to 0.01
            
            "nsteps": 2000, # numbers of time steps
            "repeat": 1, # numbers of loop. 1 or larger value
        
        	"slicefor3d_cs": 5, # for 3d_cs calculation. Enter the z coordinate.
        
            "mob": 0, # mobility in mol/J/yr. If zero, this will be calculated with the Arrhenius-type function. e.g., 1.0e-9
            "premob": 9.31092e14, # pre-exponential factor of the mobility in mol/J/yr ex. 9.31092e14
            "actEmob": 3.74389e5, # activation energy of ordering in J/mol. ex. 3.74389e5
       
            }
        
        
VisualParams = {
    
            # not required for PF calculation itself, but required for summary calculation
            
            "ApdMeasureWidth": 10, # the number of mesh of sieving. tentatively fixed to 10
            "MaxOffset": 100, # for calculation of autocorrelation. 1/4 of the axis length?
            
            # not required for PF calculation itself, but required for summary calculation and visualization
            
            "slicefor3d": 5, # for 3d to 2d convert. Enter the z coordinate.
    		"binthres": 0.4, # order parameter threshold for binarization. tentatively fixed to 0.4
            "plotNo": np.arange(0,2001,100), #　list of time steps for making plots
    		
    		# required only for DF-image like visualization
    		
    		# Conversion of the absolute value of degree of order (phi) into brightness (B) is defined as: 
    		# a convex downward function where "phiadjust" < phi < 1 (B = "Icoef" * (phi - phiadjust)^2 + "y1"),
    		# and a convex upward function where 0 < phi < "phiadjust" (B = "y1" - "y1" * (phi - "phiadjust")^2 / "phiadjust"^2).
    		# "phiadjust" and "Icoef" values should be directly fixed below, and "y1" will be calculated as "y1" = "Ifactor" * B(phi = 1).
    		# "saturation" designates the saturation level as a ratio against the B(phi = 1) value.
    		
            "saturation": 0.4, # 0-1 value
            
            "phiadjust": 0.1, # tentatively fixed to 0.1. (Do not enter "0" here.)
            "Ifactor": 0.03, # 0-1 value. tentatively fixed to 0.03.
            "Icoef": 5, # tentatively fixed to 5.
            
             
            
            }

ParamsSummary = {"Calc": CalcParams, "Visual": VisualParams}


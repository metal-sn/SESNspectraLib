import numpy as np

# initial parameter values

# v/1000 in km/s, sigma/10 in angstrom, y-amplitude, wave-range in angstrom
P0 = {'Fe' : np.array([11.0, 1.0, 1.0, 1.0])}

# prior for v/1000, sigma/10, y-amplitude
Prior = {'Fe' : np.array([0, 30, 0, 5, 0, 3])}

# region to find initial template fit region
X0 = {'Fe' : np.array([4200, 4800, 5000, 5600])}



from discstress import Disc
import numpy as np


""" 
Disc Properties
"""
boreRadius=0.05 # m
rimRadius=0.4 # m
rimStress= 0.0#Pa
boreStress= 0.0 # Pa (Only arises if "shrunk-fitting")
rho=7990 # kg/m3
nu=0.3
omega=3750*np.pi/30 # rad/s
youngsMod = lambda TdegC: (-0.0727*TdegC+212.04)*1E9 # Pa. Youngs modulus as a function of temperature (could be constant as well).
thermExp = lambda TdegC: (0.0073*TdegC + 9.6201)*1E-6 # [1/°C]. Thermal expansion coeff. as a func. of temp. (could be constant as well).

"""
Thickness and temperature distributions.
"""
# Note that the radii for the two distributions does not need to be the same. 
radii = np.array([0.05, 0.4])
thick = np.array([0.1, 0.1])
thickDist = np.column_stack((radii, thick))

tempDist=(125,500)

"""
Discretization
"""
numrings=100

"""
Instantiate disc, solve and write output to file.
"""

disc=Disc(boreRadius, rimRadius, rho, omega, rimStress, boreStress, youngsMod, thermExp, nu, thickDist, tempDist, numrings)
disc.write("test.csv")


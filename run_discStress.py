"""
This script uses the Disc class from the discstress module to calculate
the stress distribution on a turbomachinery disc. 

Inputs required are:
- Bore radius.
- Rim radius.
- Rim radial stress.
    This is from the blading centrifugal stress.
- Bore radial stress
    This is usually zero for a hollow bore, but could be non-zero 
    if the disk is shrunk-fitted onto a shaft.
- Disk material density.
- Material Poisson's ratio.
- Disk angular velocity. 
- Material Young's Modulus.
    This can be a value or given as a function of temperature.
- Material thermal expansion coefficient.
    This can be a value or given as a function of temperature.
- Disk thickness distribution with radius.
    Given as a numpy 2-column array. First column is radius, 
    second column is thickness. Radius in ascending order (lower first). 
- Disk temperature distribution with radius. 
    This can be given in the same format as the thickness distribution, 
    but it does not be to defined at the same radii as it will be interpolated.
    If no temperature distribution is known but the bore and rim temperatures are known, 
    these can be given as tuple, (boreTemp, rimTemp), and a fifth-law temperature law
    will be assumed. This is a common design assumption. 
- Number of rings to be used for discretization of the disk. 

Ferrer-Vidal, 2026
ferrervl@erau.edu
"""

from discstress import Disc
import numpy as np

"""
$$$$$$$$$$$$$$
    INPUTS
$$$$$$$$$$$$$$
"""

"""
Disc properties
"""
boreRadius=0.05 # m
rimRadius=0.4 # m
rimStress= 66.774*1E6 #Pa
boreStress= 0.0 # Pa (Only arises if "shrink-fitting")
rho=7990 # kg/m3
nu=0.3 # Poisson's ratio
omega=3750*np.pi/30 # rad/s
youngsMod = lambda TdegC: (-0.0727*TdegC+212.04)*1E9 # Pa. Youngs modulus as a function of temperature (could be constant as well).
thermExp = lambda TdegC: (0.0073*TdegC + 9.6201)*1E-6 # [1/°C]. Thermal expansion coeff. as a func. of temp. (could be constant as well).

"""
Thickness and temperature distributions.
"""
# Note that the radii for the two distributions does not need to be the same. 
radii = np.array([0.05, 0.1, 0.2, 0.3, 0.375, 0.376, 0.4])
thick = np.array([0.07, 0.07, 0.0225, 0.0125, 0.01, 0.04, 0.04])

temp = np.array([125, 125, 147, 185, 252, 437, 550])

thickDist = np.column_stack((radii, thick))
tempDist = np.column_stack((radii, temp))

"""
Discretization & filename for output.
"""
numrings=100
outFileName="test.csv"


"""
$$$$$$$$$$$$$$
    INSTANTIATE DISC, SOLVE, AND WRITE DATA
$$$$$$$$$$$$$$
"""


disc=Disc(boreRadius, rimRadius, rho, omega, rimStress, boreStress, youngsMod, thermExp, nu, thickDist, tempDist, numrings)
disc.write(outFileName)


# Copyright (c) 2026 Luis E. Ferrer-Vidal España-Heredia
# SPDX-License-Identifier: BSD-2-Clause
"""
Axisymmetric rotating annular-disc stress analysis.

This module calculates radial and hoop stresses in a rotating disc with
specified bore and rim radii. The disc is discretized into a user-defined
number of annular rings. Each ring is treated as having constant thickness,
temperature, density, Young's modulus, Poisson's ratio, and thermal expansion
coefficient.

The solution is obtained by propagating the stress state from the bore to the
rim. The unknown bore hoop stress is adjusted until the specified radial
stress boundary condition at the rim is satisfied. Ring interfaces are treated
using radial force continuity and circumferential strain compatibility.

The module can also calculate plane-stress von Mises stress at the inner and
outer radii of each ring and export the calculated ring data to a CSV file for
post-processing.

Assumptions
-----------
- Axisymmetric geometry and loading
- Small deformation
- Linear elastic material response
- Plane stress with negligible axial stress
- Constant properties within each ring
- Piecewise-constant ring thickness
- Centrifugal body loading due to rotation
- Linear interpolation of supplied temperature and thickness distributions
- Thermal strain based on the specified thermal expansion model

Input data
----------
Temperature and thickness distributions are supplied as two-column arrays.
The first column contains radius and the second column contains temperature
or thickness, respectively.

Units
-----
All inputs must use a consistent unit system. For example, using SI units
requires radius and thickness in meters, density in kg/m^3, angular velocity
in rad/s, Young's modulus and stress in Pa, and temperature in a consistent
temperature scale.

Limitations
-----------
The current implementation assumes plane stress and does not include axial
stress, plasticity, creep, residual stress, or temperature-dependent
properties within an individual ring. The accuracy of the solution depends on
the number of rings used for discretization.

Attribution
----------
ferrervl@erau.edu
AI (ChatGPT GPT-5.6 Luna) has been used for docstrings and the "write" function. 
"""


from scipy.interpolate import interp1d
from scipy.optimize import minimize
import numpy as np
import csv


class Disc:
    """
    Axisymmetric rotating-disc stress model discretized into annular rings.

    Using the constructor for this class automatically solves the stress in the disc.

    Parameters
    ----------
    r_bore : float
        Inner bore radius of the disc.

    r_rim : float
        Outer rim radius of the disc. Must be greater than `r_bore`.

    rho : float
        Disc material density.

    omega : float
        Disc angular velocity in radians per second.

    st_r_rim : float
        Prescribed radial stress at the rim.

    st_r_bore : float
        Prescribed radial stress at the bore.

    youngs : float or callable
        Young's modulus. If callable, it must accept temperature as its
        argument and return Young's modulus at that temperature.

    thermExp : float or callable
        Coefficient of thermal expansion. If callable, it must accept
        temperature as its argument and return the coefficient at that
        temperature.

    poisson : float or callable
        Poisson's ratio. If callable, it must accept temperature as its
        argument and return Poisson's ratio at that temperature.

    thickDist : array-like
        Two-column radial thickness distribution. Column 0 contains radius
        and column 1 contains thickness.

    tempDist : numpy.ndarray, callable, tuple, or scalar
        Temperature specification. A two-column NumPy array contains radius
        in column 0 and temperature in column 1. A callable must accept
        radius and return temperature. A two-element tuple specifies the
        bore and rim temperatures for the built-in fifth-power temperature
        distribution.

    numRings : int
        Number of annular rings used to discretize the disc.

    Attributes
    ----------
    boreHoopGuess : float
        Initial characteristic estimate of the bore hoop stress used when
        solving for the bore hoop stress.

    omega : float
        Disc angular velocity in radians per second.

    st_r_bore : float
        Prescribed radial stress at the bore.

    st_r_rim : float
        Prescribed radial stress at the rim.

    st_h_bore : float or None
        Calculated hoop stress at the bore after the disc has been solved.

    r_bore : float
        Inner bore radius of the disc.

    r_rim : float
        Outer rim radius of the disc.

    rho : float
        Disc material density.

    radii : numpy.ndarray
        Radial coordinates defining the inner and outer boundaries of the
        discretized rings.

    rings : list of Ring
        List containing the annular Ring objects that make up the disc.

    numRings : int
        Number of annular rings used in the discretization.

    TrefAlpha : float
        Reference temperature used to calculate thermal strain.

    AveHoopStress : float
        Area-weighted mean hoop stress over the disc.

    AveVonMisesStress : float
        Area-weighted mean von Mises stress over the disc.

    st_vm_max : float
        Maximum von Mises stress calculated at the inner and outer boundaries
        of the rings.

    volume : float
        Total disc volume.

    mass : float
        Total disc mass.

    Methods
    -------
    Disc(
        r_bore,
        r_rim,
        rho,
        omega,
        st_r_rim,
        st_r_bore,
        youngs,
        thermExp,
        poisson,
        thickDist,
        tempDist,
        numRings,
    )
        Construct and solve an axisymmetric rotating-disc model.

        The disc is discretized into annular rings using the specified radial
        geometry. The constructor evaluates the temperature and thickness
        distributions, assigns material properties to each ring, solves the
        stress distribution, and calculates the resulting disc properties.

    recalculate(thk)
        Recalculate the disc stresses and properties using a new radial
        thickness distribution.

    write(filename="disc_results.csv")
        Write the ring geometry, material properties, temperatures, and
        calculated stresses to a CSV file.

    Notes
    -----
    All quantities should use a consistent unit system.

    The write function assumes all units to be SI. The temperature is used only
    to interpolate material properties and is expected in Celsius.
    """
    def __init__(self, r_bore, r_rim, rho, omega, st_r_rim, st_r_bore, youngs, thermExp, poisson, thickDist, tempDist, numRings):

        if r_rim <= r_bore:
            raise ValueError("Rim radius cannot be lower than bore radius.")
        
        if numRings<1:
            raise ValueError("The number of rings should be at least 1.")

        centrif_stress_scale=rho * omega**2 * r_rim**2
        self.boreHoopGuess=max(abs(centrif_stress_scale), abs(st_r_rim), abs(st_r_bore), 100E6)

        self.omega=omega    
        self.st_r_bore= st_r_bore # bore radial stress
        self.st_r_rim=st_r_rim # rim radial stress
        self.st_h_bore=None # hoop stress at bore, to be solved. 
        self.r_bore=r_bore
        self.r_rim=r_rim
        self.rho=rho
        self.radii=np.linspace(r_bore, r_rim, numRings+1)
        self.rings=[]
        self.numRings=numRings
        self.TrefAlpha=20 # Reference temperature for thermal expansion.

        self.AveHoopStress = 0
        self.volume = 0
        self.mass = 0

        ThickInterp = interp1d(
            thickDist[:,0], # Radii for thickness dist.
            thickDist[:,1], # Thicknesses. 
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate"
        )

        if isinstance(tempDist, np.ndarray): # If temperature is given as an array, create interpolation function.
            TempInterp = interp1d(
                tempDist[:,0], # Radii for temperature dist.
                tempDist[:,1], # Temperatures. 
                kind="linear",
                bounds_error=False,
                fill_value="extrapolate"
            )
        elif isinstance(tempDist, tuple) and len(tempDist) == 2: # If temperature is given as a tuple (T0, T1), assume 5th power law distribution.
            a=(tempDist[1]-tempDist[0])/ (self.r_rim-self.r_bore)**5
            b=tempDist[0]
            def tempFunc(r):
                return float(b + a * (r-self.r_bore)**5)
            tempDist=tempFunc

        for ri, ro in zip(self.radii, self.radii[1:]):
            rmid=(ri+ro)/2
            h=float(ThickInterp(rmid))
            if isinstance(tempDist, np.ndarray):
                T=float(TempInterp(rmid))
            elif callable(tempDist):
                T=tempDist(rmid)
            else:
                T=tempDist
            if callable(youngs):
                E=youngs(T)
            else:
                E=youngs
            if callable(thermExp):
                alpha=thermExp(T)
            else:
                alpha=thermExp
            if callable(poisson):
                nu=poisson(T)
            else:
                nu=poisson
            ring=Ring(ri, ro, h, rho, self.omega, E, nu, T, alpha)
            ring.TrefAlpha=self.TrefAlpha
            if self.rings:
                ring.prior=self.rings[-1]
            self.rings.append(ring)
        
        self.__solve__()

    def __shootHoop__(self, boreHoop):
        # Returns the radial stress at the rim for the given bore hoop stress.
        self.rings[0].st_h_i=boreHoop
        self.rings[0].st_r_i=self.st_r_bore

        for ring in self.rings: 
            ring() 

        return self.rings[-1].st_r_o

    def __solve__(self, verbose=True):
        # Solve for the bore hoop stress that matches the rim radial stress.
        def rimStressObj(boreHoop):
            rimStressCalc=self.__shootHoop__(boreHoop)
            error=(self.st_r_rim-rimStressCalc)
            return error

        res0=rimStressObj(0)
        res1=rimStressObj(self.boreHoopGuess)

        self.st_h_bore=-res0*self.boreHoopGuess/(res1-res0)
        self.__shootHoop__(self.st_h_bore)  # Update the ring stresses with the solved bore hoop stress.
        self.__discProps__()  # Calculate and store disc properties after solving.
        if verbose:
            self.__print__()

    def __print__(self):
            print("Disc properties calculated.")
            print("----------------------------")
            print(f"Maximum von Mises Stress: {self.st_vm_max/1E6:.2f} MPa")
            print(f"Bore Hoop Stress: {self.st_h_bore/1E6:.2f} MPa")
            print(f"Rim Radial Stress: {self.st_r_rim/1E6:.2f} MPa")
            print(f"Average Hoop Stress: {self.AveHoopStress/1E6:.2f} MPa")
            print(f"Volume: {self.volume:.4f} m^3") 
            print(f"Mass：{self.mass:.2f} kg") 
            print("----------------------------")

    def __discProps__(self):
        discArea = 0
        discVol = 0
        maxVonMises = 0
        discVonMisesA = 0
        discHoopA = 0
        for ring in self.rings:
            Aring = ring.h * (ring.ro - ring.ri)
            discArea += Aring
            ringHoop = (ring.st_h_i+ring.st_h_o)/2
            ringAHoop = Aring * ringHoop
            discHoopA += ringAHoop
            ringVonMisesA = Aring*(ring.st_vm_i+ ring.st_vm_o)/2
            discVonMisesA += ringVonMisesA
            ringVol=np.pi*(ring.ro**2 - ring.ri**2) * ring.h
            discVol += ringVol
            maxVonMises = max(ring.st_vm_i, ring.st_vm_o, maxVonMises)
        self.AveHoopStress = discHoopA / discArea if discArea > 0 else 0
        self.AveVonMisesStress = discVonMisesA / discArea if discArea > 0 else 0
        self.volume = discVol
        self.mass = discVol * self.rho if discVol > 0 else 0
        self.st_vm_max = maxVonMises

    def recalculate(self, thk):
        # Recalculate stresses with new thickness distribution.
        ThkInterp = interp1d(
            thk[:,0], # Radii for thickness dist.
            thk[:,1], # Thicknesses. 
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate"
        )

        for ring in self.rings:
            ring.h=float(ThkInterp((ring.ro+ring.ri)/2))

        self.__solve__(verbose=False)

    def write(self, filename="disc_results.csv"):
    # Write ring properties and stresses to a CSV file.

        if self.st_h_bore is None:
            raise RuntimeError("The disc must be solved before writing results.")

        fieldnames = [
            "Ring",
            "Radius [m]",
            "Thickness [mm]",
            "Density [kg/m^3]",
            "omega [rad/s]",
            "Temperature [°C]",
            "Young's Modulus [MPa]",
            "Poisson's ratio",
            "Thermal expansion [1/°C]",
            "Stress_Hoop [MPa]",
            "Stress_Radial [MPa]",
            "Stress_vonMises [MPa]"
        ]

        rows = []

        for i, ring in enumerate(self.rings):
            points = [
                (ring.ri, ring.st_h_i, ring.st_r_i, ring.st_vm_i),
                (ring.ro, ring.st_h_o, ring.st_r_o, ring.st_vm_o),
            ]

            for  radius, hoop_stress, radial_stress, vonMises_stress in points:
                rows.append({
                    "Ring": i + 1,
                    "Radius [m]": radius,
                    "Thickness [mm]": ring.h*1e3,
                    "Density [kg/m^3]": ring.rho,
                    "omega [rad/s]": ring.omega,
                    "Temperature [°C]": ring.T,
                    "Young's Modulus [MPa]": ring.E/1E6,
                    "Poisson's ratio": ring.nu,
                    "Thermal expansion [1/°C]": ring.alpha,
                    "Stress_Hoop [MPa]": hoop_stress/1E6,
                    "Stress_Radial [MPa]": radial_stress/1E6,  
                    "Stress_vonMises [MPa]": vonMises_stress/1E6,
                })

        with open(filename, "w", newline="", encoding="utf-8") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return filename

class Ring:
    """
    Constant-property annular region within a discretized rotating disc.

    A `Ring` represents the region between inner radius `ri` and outer
    radius `ro`. Its thickness, density, temperature, Young's modulus,
    Poisson's ratio, and thermal expansion coefficient are treated as
    constant throughout the ring.

    The ring receives radial and hoop stresses at its inner radius and
    calculates the corresponding stresses at its outer radius using the
    analytical solution for a rotating annular region.

    For rings other than the innermost ring, the inner boundary stresses
    are obtained from the preceding ring. Radial force continuity and
    circumferential strain compatibility are applied at the interface.

    The current formulation assumes axisymmetric plane stress, small
    deformation, linear elasticity, and no axial stress.

    Parameters
    ----------
    ri : float
        Inner radius of the ring.

    ro : float
        Outer radius of the ring. Must be greater than `ri`.

    thickness : float
        Constant ring thickness.

    rho : float
        Material density within the ring.

    omega : float
        Angular velocity in radians per second.

    E : float
        Young's modulus within the ring.

    nu : float
        Poisson's ratio within the ring.

    Temp : float, optional
        Ring temperature.

    alpha : float, optional
        Coefficient of thermal expansion within the ring.

    Attributes
    ----------
    ri : float
        Inner radius of the ring.

    ro : float
        Outer radius of the ring.

    rho : float
        Material density within the ring.

    omega : float
        Angular velocity in radians per second.

    T : float or None
        Ring temperature.

    alpha : float or None
        Coefficient of thermal expansion within the ring.

    E : float
        Young's modulus within the ring.

    h : float
        Constant ring thickness.

    nu : float
        Poisson's ratio within the ring.

    TrefAlpha : float
        Reference temperature used to calculate thermal strain.

    st_h_i : float or None
        Hoop stress at the inner radius of the ring.

    st_r_i : float or None
        Radial stress at the inner radius of the ring.

    st_h_o : float or None
        Hoop stress at the outer radius of the ring.

    st_r_o : float or None
        Radial stress at the outer radius of the ring.

    st_vm_i : float or None
        Plane-stress von Mises stress at the inner radius of the ring.

    st_vm_o : float or None
        Plane-stress von Mises stress at the outer radius of the ring.

    prior : Ring or None
        Reference to the preceding ring in the radial sequence. This is
        `None` for the innermost ring.
    """
    def __init__(self, ri, ro, thickness, rho, omega, E, nu, Temp=None, alpha=None):
        
        # Ring data
        self.ri=ri
        self.ro=ro
        self.rho=rho
        self.omega=omega
        self.T=Temp
        self.alpha=alpha
        self.E=E
        self.h=thickness
        self.nu=nu
        self.TrefAlpha=20 # Reference temperature for thermal expansion.

        # Ring hoop and radial stresses at inner and outer radii
        self.st_h_i=None
        self.st_r_i=None
        self.st_h_o=None
        self.st_r_o=None

        # Previous ring reference, to pass an update inner stresses.
        self.prior=None 

    def __call__(self):
        if self.prior is not None: # This does not apply for innermost ring. 
            # Inner radial stress for this ring, from equilibrium. 
            self.st_r_i=self.prior.st_r_o*self.prior.h/self.h  
            
            # Inner hoop stress for this ring, from strain compatibility.
            self.st_h_i=(
                self.E*(
                    (1/self.prior.E)*(self.prior.st_h_o-self.prior.nu*self.prior.st_r_o)
                    +self.prior.alpha*(self.prior.T-self.prior.TrefAlpha)
                    -self.alpha*(self.T-self.TrefAlpha)
                    +self.st_r_i*self.nu/self.E
                )
            )
        else:
            if self.st_h_i is None or self.st_r_i is None:
                raise RuntimeError("Bore stresses have not be initialized.")

        # Calculate inner radius von mises stress:
        self.st_vm_i=(self.st_r_i**2-self.st_r_i*self.st_h_i+self.st_h_i**2)**0.5
        
        # Calculate sums and differences at inner radius of ring
        Si=self.st_h_i+self.st_r_i
        Di=self.st_h_i-self.st_r_i

        # Calculate sums and differences at outer radius of ring. 
        So=Si-(1+self.nu)/2*self.rho*self.omega**2*(self.ro**2-self.ri**2)
        Do=Di*(self.ri**2/self.ro**2)-((1-self.nu)*self.rho*self.omega**2)/4*(((self.ri**4)/self.ro**2)-self.ro**2)

        #Calculate outer radius hoop and radial stresses:
        self.st_h_o=(So+Do)/2
        self.st_r_o=(So-Do)/2

        # Calculate outer radius von mises stress:
        self.st_vm_o=(self.st_r_o**2-self.st_r_o*self.st_h_o+self.st_h_o**2)**0.5






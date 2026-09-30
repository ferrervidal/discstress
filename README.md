# Turbomachinery Disc Stress Analysis

This module calculates radial, hoop, and plane-stress von Mises stresses in
an axisymmetric rotating turbomachinery disc. The disc is divided into a
user-defined number of annular rings. Thickness, temperature, density, and
material properties are treated as constant within each ring.

## Files

```text
discstress/
├── discstress.py
├── run_discStress.py
├── run_calcConstDisc.py
├── run_validate.py
├── README.md
└── Post.xlsx
```

`discstress.py` contains the `Disc` and `Ring` classes. `run_discStress.py`
is an example calculation script. `run_calcConstDisc.py` is an another, less-detailed
example for a constant thickness disc with simple temperature distribution. `run_validate.py` 
is a simple validation case to test the code is working well. `Post.xlsx` has sample output
and plots.

## Requirements and setup

The code requires Python, NumPy, and SciPy:

```bash
python -m pip install numpy scipy
```

## Running the example

Place `run_discStress.py` and `discstress.py` in the same folder, then run:

```bash
python run_discStress.py
```

The script constructs and solves a `Disc` object and writes the results to
`test.csv`. The output filename can be changed with:

```python
outFileName = "test.csv"
```

## Basic usage

```python
from discstress import Disc

disc = Disc(
    boreRadius, rimRadius, rho, omega, rimStress, boreStress,
    youngsMod, thermExp, nu, thickDist, tempDist, numrings,
)

disc.write("disc_results.csv")
```

The disc is solved automatically when the `Disc` object is constructed.

## Inputs

All quantities must use a consistent unit system. The example uses SI units:

```python
boreRadius = 0.05          # m
rimRadius = 0.40           # m
rimStress = 66.774e6       # Pa
boreStress = 0.0           # Pa
rho = 7990                 # kg/m^3
nu = 0.3                   # dimensionless
omega = 3750 * np.pi / 30  # rad/s
```

`rimStress` is the prescribed radial stress at the rim, for example from
blade or blade-attachment loading. `boreStress` is the prescribed radial
stress at the bore. It is normally zero for a free hollow bore, but may be
nonzero for a shrink fit.

Angular velocity must be supplied in radians per second. To convert rpm:

```python
omega = rpm * np.pi / 30
```

### Material properties

Young's modulus and thermal expansion coefficient may be constants or
functions of temperature:

```python
youngsMod = 200e9
thermExp = 12e-6
```

For temperature-dependent properties:

```python
def youngsMod(TdegC):
    return (-0.0727 * TdegC + 212.04) * 1e9


def thermExp(TdegC):
    return (0.0073 * TdegC + 9.6201) * 1e-6
```

Temperature is supplied in degrees Celsius in the example and is passed to
the material-property functions in that form.

### Thickness distribution

Thickness is supplied as a two-column NumPy array. Column 0 contains radius
in metres and column 1 contains thickness in metres:

```python
radii = np.array([0.05, 0.10, 0.20, 0.30, 0.375, 0.376, 0.40])
thick = np.array([0.07, 0.07, 0.0225, 0.0125, 0.01, 0.04, 0.04])
thickDist = np.column_stack((radii, thick))
```

The radius values must be in ascending order and should cover the bore-to-rim
range. Thickness is linearly interpolated at the midpoint of each analysis
ring, and the interpolated value is treated as constant within that ring.

### Temperature distribution

Temperature can be supplied in the same two-column format:

```python
temp = np.array([125, 125, 147, 185, 252, 437, 550])
tempDist = np.column_stack((radii, temp))
```

The temperature-distribution radii do not need to be the same as the
thickness-distribution radii.

Alternatively, provide bore and rim temperatures as a tuple:

```python
tempDist = (125, 550)
```

The model then assumes a fifth-power radial temperature distribution between
the bore and rim.

### Discretization

Set the number of annular rings with:

```python
numrings = 100
```

Increase this value and repeat the calculation to check mesh convergence,
especially where thickness or temperature varies rapidly.

## CSV output

The CSV file contains two rows for each ring, corresponding to the inner and
outer ring radii. It includes radius, thickness, material properties,
temperature, hoop stress, radial stress, and plane-stress von Mises stress.

Output units are radius in metres, thickness in millimetres, density in
kg/m^3, angular velocity in rad/s, temperature in degrees Celsius, Young's
modulus in MPa, and stresses in MPa.

A "Post.xlsx" file shows some sample output and graphs.

## Model assumptions and limitations

The current model assumes:

- Axisymmetric geometry and loading
- Small deformation
- Linear elastic material behaviour
- Plane stress with negligible axial stress
- Constant properties within each ring
- Piecewise-constant ring thickness
- Centrifugal loading due to rotation
- Radial force continuity at ring interfaces
- Circumferential strain compatibility at ring interfaces

The model does not currently include plasticity, creep, residual stress,
fatigue, local geometric stress concentrations, or a full three-dimensional
stress state. The elastic von Mises stress should not by itself be treated as
a complete prediction of disc burst or service life.
```

import numpy as np

from discstress import Disc


def analytical_stresses(r, r_bore, r_rim, rho, omega, nu,
                        st_r_bore, st_r_rim):
    """
    Analytical stresses for a constant-property rotating annulus.
    """

    q = rho * omega**2
    c_r = (3.0 + nu) * q / 8.0

    x_bore = st_r_bore + c_r * r_bore**2
    x_rim = st_r_rim + c_r * r_rim**2

    B = (x_bore - x_rim) / (
        1.0 / r_rim**2 - 1.0 / r_bore**2
    )

    A = x_bore + B / r_bore**2

    radial_stress = (
        A
        - B / r**2
        - (3.0 + nu) * q * r**2 / 8.0
    )

    hoop_stress = (
        A
        + B / r**2
        - (1.0 + 3.0 * nu) * q * r**2 / 8.0
    )

    return radial_stress, hoop_stress


def main():
    # Geometry and loading
    r_bore = 0.05
    r_rim = 0.40

    rho = 7990.0
    omega = 3750.0 * np.pi / 30.0
    nu = 0.30

    st_r_bore = 0.0
    st_r_rim = 0.0

    # Constant properties
    E = 200.0e9
    alpha = 0.0
    temperature = 300.0
    thickness = 0.01

    # Increase this to check mesh convergence
    num_rings = 1

    # Distribution arrays must have shape (number_of_points, 2)
    input_radii = np.array([r_bore, r_rim])

    thick_dist = np.column_stack((
        input_radii,
        [thickness, thickness],
    ))

    temp_dist = np.column_stack((
        input_radii,
        [temperature, temperature],
    ))

    disc = Disc(
        r_bore=r_bore,
        r_rim=r_rim,
        rho=rho,
        omega=omega,
        st_r_rim=st_r_rim,
        st_r_bore=st_r_bore,
        youngs=E,
        thermExp=alpha,
        poisson=nu,
        thickDist=thick_dist,
        tempDist=temp_dist,
        numRings=num_rings,
    )


    maximum_radial_error = 0.0
    maximum_hoop_error = 0.0
    maximum_vm_error = 0.0

    for ring in disc.rings:
        numerical_points = [
            (ring.ri, ring.st_r_i, ring.st_h_i, ring.st_vm_i),
            (ring.ro, ring.st_r_o, ring.st_h_o, ring.st_vm_o),
        ]

        for radius, radial_num, hoop_num, vm_num in numerical_points:
            radial_exact, hoop_exact = analytical_stresses(
                radius,
                r_bore,
                r_rim,
                rho,
                omega,
                nu,
                st_r_bore,
                st_r_rim,
            )

            vm_exact = np.sqrt(
                radial_exact**2
                - radial_exact * hoop_exact
                + hoop_exact**2
            )

            maximum_radial_error = max(
                maximum_radial_error,
                abs(radial_num - radial_exact),
            )

            maximum_hoop_error = max(
                maximum_hoop_error,
                abs(hoop_num - hoop_exact),
            )

            maximum_vm_error = max(
                maximum_vm_error,
                abs(vm_num - vm_exact),
            )

    stress_scale = max(
        abs(disc.rings[0].st_h_i),
        abs(disc.rings[-1].st_h_o),
        1.0,
    )

    relative_radial_error = maximum_radial_error / stress_scale
    relative_hoop_error = maximum_hoop_error / stress_scale
    relative_vm_error = maximum_vm_error / stress_scale

    rim_residual = (
        disc.rings[-1].st_r_o - st_r_rim
    )

    print(f"Bore hoop stress:       {disc.st_h_bore:.6e} Pa")
    print(f"Rim radial residual:    {rim_residual:.6e} Pa")
    print(f"Maximum radial error:   {maximum_radial_error:.6e} Pa")
    print(f"Maximum hoop error:     {maximum_hoop_error:.6e} Pa")
    print(f"Maximum von Mises error:{maximum_vm_error:.6e} Pa")
    print(f"Relative radial error:  {relative_radial_error:.6e}")
    print(f"Relative hoop error:    {relative_hoop_error:.6e}")
    print(f"Relative von Mises error: {relative_vm_error:.6e}")

    tolerance = 1.0e-6

    assert relative_radial_error < tolerance
    assert relative_hoop_error < tolerance
    assert relative_vm_error < tolerance

    print("Validation passed.")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
analyze_microrheology.py

Analyze microrheology simulation data to extract rheological properties.

This script processes LAMMPS dump files from probe cells and nearby cells to:
1. Compute probe cell displacement and velocity
2. Calculate Mean Square Displacement (MSD) of probe cells
3. Analyze response of nearby cells (displacement correlation)
4. Estimate effective viscosity from probe cell motion
5. Visualize spatial decay of displacement response
"""

import numpy as np
import sys
import os
from pathlib import Path


def read_lammps_dump(filename, max_frames=None):
    """
    Read LAMMPS dump file and extract atom positions and properties.

    Returns:
        List of dictionaries, one per timestep, containing:
        - timestep: int
        - box: [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
        - atoms: structured array with columns from dump file
    """
    frames = []

    if not os.path.exists(filename):
        print(f"Warning: File {filename} not found")
        return frames

    with open(filename, 'r') as f:
        while True:
            # Read timestep header
            line = f.readline()
            if not line:  # EOF
                break

            if "ITEM: TIMESTEP" not in line:
                continue

            timestep = int(f.readline().strip())

            # Read number of atoms
            line = f.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().strip())

            # Read box bounds
            line = f.readline()  # ITEM: BOX BOUNDS
            box = []
            for _ in range(3):
                bounds = list(map(float, f.readline().split()[:2]))
                box.append(bounds)

            # Read column headers
            header_line = f.readline()  # ITEM: ATOMS ...
            columns = header_line.split()[2:]  # Skip "ITEM: ATOMS"

            # Read atom data
            atom_data = []
            for _ in range(natoms):
                values = list(map(float, f.readline().split()))
                atom_data.append(values)

            # Convert to structured array
            atom_array = np.array(atom_data)

            frames.append({
                'timestep': timestep,
                'box': box,
                'columns': columns,
                'data': atom_array
            })

            if max_frames and len(frames) >= max_frames:
                break

    return frames


def compute_com_trajectory(frames):
    """
    Compute center-of-mass trajectory from dump frames.
    Each frame contains beads from a single cell.

    Returns:
        times: array of timesteps
        com_x: array of COM x-coordinates
        com_y: array of COM y-coordinates
    """
    times = []
    com_x = []
    com_y = []

    for frame in frames:
        data = frame['data']
        cols = frame['columns']

        # Find x, y columns (unwrapped coordinates preferred)
        if 'xu' in cols and 'yu' in cols:
            x_idx = cols.index('xu')
            y_idx = cols.index('yu')
        elif 'x' in cols and 'y' in cols:
            x_idx = cols.index('x')
            y_idx = cols.index('y')
        else:
            print("Error: Cannot find x, y coordinates in dump file")
            return None, None, None

        # Compute COM
        x_com = np.mean(data[:, x_idx])
        y_com = np.mean(data[:, y_idx])

        times.append(frame['timestep'])
        com_x.append(x_com)
        com_y.append(y_com)

    return np.array(times), np.array(com_x), np.array(com_y)


def compute_msd(times, x, y):
    """
    Compute Mean Square Displacement.

    MSD(t) = <[r(t0+t) - r(t0)]^2>
    """
    # Use initial position as reference
    x0, y0 = x[0], y[0]

    dx = x - x0
    dy = y - y0
    msd = dx**2 + dy**2

    return times, msd


def compute_velocity(times, x, y):
    """
    Compute velocity from position trajectory.

    v = dr/dt (using central differences)
    """
    dt = np.diff(times)

    # Central differences for interior points
    vx = np.zeros_like(x)
    vy = np.zeros_like(y)

    if len(x) >= 3:
        # Central difference for interior
        vx[1:-1] = (x[2:] - x[:-2]) / (times[2:] - times[:-2])
        vy[1:-1] = (y[2:] - y[:-2]) / (times[2:] - times[:-2])

        # Forward/backward difference for endpoints
        vx[0] = (x[1] - x[0]) / (times[1] - times[0])
        vy[0] = (y[1] - y[0]) / (times[1] - times[0])
        vx[-1] = (x[-1] - x[-2]) / (times[-1] - times[-2])
        vy[-1] = (y[-1] - y[-2]) / (times[-1] - times[-2])

    return vx, vy


def estimate_viscosity(times, msd, force, temp=0.0):
    """
    Estimate effective viscosity from probe particle MSD.

    For a driven particle in a viscous medium:
    At short times: MSD ~ (F/γ)^2 * t^2  (ballistic)
    At long times: MSD ~ 2*D*t where D = kT/γ + F^2/(2γ^2)  (diffusive)

    For T=0 (no thermal fluctuations):
    Long time: MSD ~ (F/γ)^2 * t^2 + 2*(F^2/(2*γ^2))*t

    We fit the linear regime to extract γ (drag coefficient).
    Viscosity η ~ γ / (6πR) for a sphere (Stokes law)

    Returns:
        gamma: drag coefficient
        viscosity_estimate: rough estimate assuming sphere
    """
    # Fit MSD vs time in linear regime
    # Look for region where d(MSD)/dt is approximately constant

    if len(times) < 10:
        return None, None

    # Use later part of trajectory (after initial transient)
    start_idx = len(times) // 4
    end_idx = 3 * len(times) // 4

    t_fit = times[start_idx:end_idx]
    msd_fit = msd[start_idx:end_idx]

    # Linear fit: MSD = a*t + b
    coeffs = np.polyfit(t_fit, msd_fit, 1)
    slope = coeffs[0]  # d(MSD)/dt

    # For driven particle: d(MSD)/dt ~ 2*F*v_steady = 2*F*(F/γ) = 2*F^2/γ
    # Therefore: γ ~ 2*F^2 / slope

    if slope > 0:
        gamma = 2 * force**2 / slope

        # Rough estimate of viscosity (assuming effective radius ~ 0.5 LJ)
        R_eff = 0.5
        viscosity = gamma / (6 * np.pi * R_eff)

        return gamma, viscosity
    else:
        return None, None


def analyze_nearby_response(probe_times, probe_x, probe_y, nearby_frames):
    """
    Analyze displacement response of nearby cells.

    Computes displacement correlation: how much nearby cells move in the
    direction of probe cell motion.

    Returns:
        mean_displacement: mean displacement of nearby cells' COM
        correlation: correlation coefficient with probe displacement
    """
    if not nearby_frames:
        return None, None

    # Compute COM trajectory of nearby cells
    times, nx, ny = compute_com_trajectory(nearby_frames)

    if times is None or len(times) < 2:
        return None, None

    # Interpolate to match probe timesteps
    probe_dx = probe_x - probe_x[0]
    probe_dy = probe_y - probe_y[0]

    nearby_dx = nx - nx[0]
    nearby_dy = ny - ny[0]

    # Compute displacement magnitude
    nearby_disp = np.sqrt(nearby_dx**2 + nearby_dy**2)

    # Mean displacement
    mean_disp = np.mean(nearby_disp)

    # Correlation with probe x-displacement (direction of force)
    if len(times) == len(probe_times):
        correlation = np.corrcoef(probe_dx, nearby_dx)[0, 1]
    else:
        correlation = None

    return mean_disp, correlation


def main():
    print("="*70)
    print("Microrheology Analysis")
    print("="*70)
    print()

    # Simulation parameters
    number_cells = 900
    fd = 3.0  # Total force on probe cell

    data_dir = Path("/mnt/d/NextTissUe_data")

    # If data directory doesn't exist, try local directory
    if not data_dir.exists():
        data_dir = Path(".")
        print(f"Data directory /mnt/d/NextTissUe_data not found, using current directory")

    # File paths
    probe1_file = data_dir / f"microrheology_probe_NC{number_cells}_fd{fd}.dat"
    probe2_file = data_dir / f"microrheology_probe2_NC{number_cells}_fd{fd}.dat"
    first_shell_file = data_dir / f"microrheology_first_shell_NC{number_cells}_fd{fd}.dat"
    second_shell_file = data_dir / f"microrheology_second_shell_NC{number_cells}_fd{fd}.dat"
    third_shell_file = data_dir / f"microrheology_third_shell_NC{number_cells}_fd{fd}.dat"

    print("Reading probe cell 1 data...")
    probe1_frames = read_lammps_dump(str(probe1_file), max_frames=1000)

    if not probe1_frames:
        print("Error: No data found for probe cell 1")
        print(f"Expected file: {probe1_file}")
        print("\nPlease run the LAMMPS simulation first to generate data.")
        return

    print(f"  Loaded {len(probe1_frames)} frames")

    print("Reading probe cell 2 data...")
    probe2_frames = read_lammps_dump(str(probe2_file), max_frames=1000)
    print(f"  Loaded {len(probe2_frames)} frames")

    print("Reading nearby cells data...")
    first_frames = read_lammps_dump(str(first_shell_file), max_frames=1000)
    print(f"  First shell: {len(first_frames)} frames")

    second_frames = read_lammps_dump(str(second_shell_file), max_frames=1000)
    print(f"  Second shell: {len(second_frames)} frames")

    third_frames = read_lammps_dump(str(third_shell_file), max_frames=1000)
    print(f"  Third shell: {len(third_frames)} frames")

    print()
    print("="*70)
    print("Analyzing Probe Cell 1:")
    print("="*70)

    # Compute COM trajectory
    times, x, y = compute_com_trajectory(probe1_frames)

    print(f"Initial position: ({x[0]:.3f}, {y[0]:.3f})")
    print(f"Final position:   ({x[-1]:.3f}, {y[-1]:.3f})")
    print(f"Total displacement: {np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2):.3f} LJ")
    print(f"X-displacement: {x[-1]-x[0]:.3f} LJ (direction of force)")
    print(f"Y-displacement: {y[-1]-y[0]:.3f} LJ (perpendicular)")

    # Compute velocity
    vx, vy = compute_velocity(times, x, y)
    v_mag = np.sqrt(vx**2 + vy**2)

    print(f"\nMean velocity magnitude: {np.mean(v_mag):.6f} LJ/timestep")
    print(f"Mean vx (force direction): {np.mean(vx):.6f} LJ/timestep")
    print(f"Mean vy (perpendicular): {np.mean(vy):.6f} LJ/timestep")

    # Compute MSD
    t_msd, msd = compute_msd(times, x, y)

    print(f"\nMSD at final time: {msd[-1]:.3f} LJ^2")

    # Estimate viscosity
    gamma, eta = estimate_viscosity(times, msd, force=fd)

    if gamma:
        print(f"\nEstimated drag coefficient γ: {gamma:.3f}")
        print(f"Estimated effective viscosity η: {eta:.3f} (LJ units)")
    else:
        print("\nCould not estimate viscosity (insufficient data or non-linear regime)")

    # Analyze nearby cells
    print()
    print("="*70)
    print("Analyzing Nearby Cell Response:")
    print("="*70)

    shells = [
        ("First shell (0-3 LJ)", first_frames),
        ("Second shell (3-6 LJ)", second_frames),
        ("Third shell (6-9 LJ)", third_frames)
    ]

    for shell_name, shell_frames in shells:
        print(f"\n{shell_name}:")

        if not shell_frames:
            print("  No data available")
            continue

        mean_disp, corr = analyze_nearby_response(times, x, y, shell_frames)

        if mean_disp is not None:
            print(f"  Mean displacement: {mean_disp:.4f} LJ")
            if corr is not None:
                print(f"  Correlation with probe x-displacement: {corr:.3f}")
        else:
            print("  Could not compute response")

    print()
    print("="*70)
    print("Analysis complete!")
    print("="*70)
    print()
    print("Summary:")
    print(f"  - Probe cell displacement: {x[-1]-x[0]:.3f} LJ in x-direction")
    print(f"  - Applied force: {fd:.1f} LJ units (total)")
    if gamma:
        print(f"  - Drag coefficient: {gamma:.3f}")
        print(f"  - Effective viscosity: {eta:.3f}")
    print()
    print("Interpretation:")
    print("  - Larger drag coefficient → higher resistance to motion")
    print("  - Positive correlation in nearby shells → cells move together")
    print("  - Decreasing correlation with distance → localized response")
    print()

    # Save results to file
    output_file = "microrheology_results.txt"
    with open(output_file, 'w') as f:
        f.write("Microrheology Analysis Results\n")
        f.write("="*70 + "\n\n")
        f.write(f"Probe cell 1:\n")
        f.write(f"  Initial position: ({x[0]:.6f}, {y[0]:.6f})\n")
        f.write(f"  Final position:   ({x[-1]:.6f}, {y[-1]:.6f})\n")
        f.write(f"  X-displacement: {x[-1]-x[0]:.6f} LJ\n")
        f.write(f"  Y-displacement: {y[-1]-y[0]:.6f} LJ\n")
        f.write(f"  Mean velocity: {np.mean(vx):.6f} LJ/timestep\n")
        f.write(f"  MSD (final): {msd[-1]:.6f} LJ^2\n")
        if gamma:
            f.write(f"  Drag coefficient: {gamma:.6f}\n")
            f.write(f"  Viscosity estimate: {eta:.6f}\n")
        f.write("\n")

    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()

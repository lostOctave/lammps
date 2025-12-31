#!/usr/bin/env python3
"""
identify_nearby_cells.py

Identify cells near the probe cells for microrheology analysis.
Creates distance shells around probe cells to track rheological response.
"""

import math
import sys


def distance(x1, y1, x2, y2):
    """Calculate Euclidean distance between two points."""
    return math.sqrt((x2-x1)**2 + (y2-y1)**2)


def main():
    # Probe cells (from find_probe_cells.py output)
    probe_mol_1 = 709
    probe_mol_2 = 710

    # Initial positions (after 2x scaling, before compression)
    probe1_pos = (1.50, 26.85)
    probe2_pos = (2.50, 28.58)

    # After compression to hexatic lattice, the box changes:
    # Initial: 60.0 × 51.6 LJ -> Final: ~32.24 × 27.92 LJ
    # We need to identify cells based on their INITIAL positions,
    # since molecular IDs don't change during compression

    input_file = 'initial_positions_xy.dat'

    print("="*70)
    print("Identifying Nearby Cells for Microrheology Analysis")
    print("="*70)
    print(f"Probe cell 1: mol {probe_mol_1} at {probe1_pos}")
    print(f"Probe cell 2: mol {probe_mol_2} at {probe2_pos}")
    print()

    # Read all cell positions
    cells = []
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            x, y = map(float, line.split())
            # Apply 2× scaling (same as in generate_nexttissue_data.py)
            cells.append((i+1, x * 2.0, y * 2.0))  # (mol_id, x_scaled, y_scaled)

    print(f"Read {len(cells)} cell centers")

    # Define distance shells (in LJ units, after compression)
    # After compression, the lattice spacing is ~1.0 LJ
    # We'll define shells based on multiples of lattice spacing
    shells = [
        (0, 3.0, "first_shell"),   # 0-3 LJ: nearest neighbors (1-2 lattice spacings)
        (3.0, 6.0, "second_shell"), # 3-6 LJ: second neighbors
        (6.0, 9.0, "third_shell"),  # 6-9 LJ: third neighbors
    ]

    # Calculate distance of each cell from BOTH probe cells
    # We'll use the minimum distance to either probe
    nearby_cells = {shell[2]: [] for shell in shells}

    for mol_id, x, y in cells:
        # Skip probe cells themselves
        if mol_id == probe_mol_1 or mol_id == probe_mol_2:
            continue

        # Calculate distance to both probe cells
        d1 = distance(x, y, probe1_pos[0], probe1_pos[1])
        d2 = distance(x, y, probe2_pos[0], probe2_pos[1])

        # Use minimum distance
        min_dist = min(d1, d2)

        # Assign to appropriate shell
        for rmin, rmax, shell_name in shells:
            if rmin < min_dist <= rmax:
                nearby_cells[shell_name].append((mol_id, x, y, min_dist))
                break

    # Print results
    print()
    print("="*70)
    print("NEARBY CELLS BY DISTANCE SHELL:")
    print("="*70)

    for rmin, rmax, shell_name in shells:
        cells_in_shell = nearby_cells[shell_name]
        print(f"\n{shell_name.upper()} ({rmin:.1f} - {rmax:.1f} LJ):")
        print(f"  Number of cells: {len(cells_in_shell)}")

        if cells_in_shell:
            # Sort by distance
            cells_in_shell.sort(key=lambda c: c[3])

            # Print first 5 cells as examples
            print(f"  First 5 cells (sorted by distance):")
            for mol, x, y, d in cells_in_shell[:5]:
                print(f"    mol {mol:3d}: pos=({x:6.2f}, {y:6.2f}), dist={d:.2f}")

            if len(cells_in_shell) > 5:
                print(f"    ... and {len(cells_in_shell)-5} more")

    # Generate LAMMPS commands
    print()
    print("="*70)
    print("LAMMPS COMMANDS TO ADD TO in.demo_hexatic:")
    print("="*70)
    print()
    print("# Define groups for nearby cells at different distance shells")

    for rmin, rmax, shell_name in shells:
        cells_in_shell = nearby_cells[shell_name]
        if not cells_in_shell:
            print(f"# {shell_name}: No cells found")
            continue

        # Generate group command with molecule IDs
        mol_ids = [str(c[0]) for c in cells_in_shell]

        # LAMMPS has a limit on command line length, so we need to split if too many
        max_ids_per_line = 50

        if len(mol_ids) <= max_ids_per_line:
            # Single line
            mol_id_str = " ".join(mol_ids)
            print(f"group {shell_name} molecule {mol_id_str}")
        else:
            # Multiple lines using union
            print(f"group {shell_name} molecule {' '.join(mol_ids[:max_ids_per_line])}")
            for i in range(max_ids_per_line, len(mol_ids), max_ids_per_line):
                chunk = mol_ids[i:i+max_ids_per_line]
                mol_id_str = " ".join(chunk)
                print(f"group {shell_name} union {shell_name} molecule {mol_id_str}")

        # Add compute for COM tracking
        print(f"compute com_{shell_name} {shell_name} com")
        print()

    # Generate thermo variables
    print("# Variables to track nearby cell positions")
    for rmin, rmax, shell_name in shells:
        if nearby_cells[shell_name]:
            print(f"variable {shell_name}_x equal c_com_{shell_name}[1]")
            print(f"variable {shell_name}_y equal c_com_{shell_name}[2]")

    print()
    print("# Update thermo output to include nearby cell tracking")
    thermo_vars = ["step", "temp", "pe", "ke", "etotal", "press"]
    thermo_vars += ["v_probe1_x", "v_probe1_y", "v_probe2_x", "v_probe2_y"]
    for rmin, rmax, shell_name in shells:
        if nearby_cells[shell_name]:
            thermo_vars.append(f"v_{shell_name}_x")
            thermo_vars.append(f"v_{shell_name}_y")

    print(f"thermo_style custom {' '.join(thermo_vars)}")

    print()
    print("="*70)
    print("SPECIALIZED DUMP FOR MICRORHEOLOGY:")
    print("="*70)
    print("""
# Dump probe cells and nearby cells at high frequency
dump MICRO_DUMP probe_cell_1 custom 100 microrheology_probe.dat &
     id mol xu yu d_xc d_yc d_fx_active d_fy_active fx fy
dump_modify MICRO_DUMP sort id append yes

# Also dump first shell cells
dump MICRO_NEARBY first_shell custom 100 microrheology_nearby.dat &
     id mol xu yu d_xc d_yc fx fy
dump_modify MICRO_NEARBY sort id append yes
""")

    print("="*70)
    print()
    print(f"Total cells tracked in all shells: {sum(len(v) for v in nearby_cells.values())}")
    print("="*70)

    # Write output file for later analysis
    output_file = 'nearby_cells.txt'
    with open(output_file, 'w') as f:
        f.write("# Nearby cells for microrheology analysis\n")
        f.write(f"# Probe cells: {probe_mol_1}, {probe_mol_2}\n")
        f.write("# Format: mol_id x y distance shell_name\n")
        for rmin, rmax, shell_name in shells:
            for mol, x, y, d in nearby_cells[shell_name]:
                f.write(f"{mol} {x:.6f} {y:.6f} {d:.6f} {shell_name}\n")

    print(f"Nearby cell data written to {output_file}")

    return nearby_cells


if __name__ == "__main__":
    main()

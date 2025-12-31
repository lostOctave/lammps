#!/usr/bin/env python3
"""
find_probe_cells.py

Identify two neighboring cells for microrheology test.
Cells should be:
- Near left edge of box
- Approximately at center of y-axis
- Neighbors (distance ≈ 1.0 after scaling)
- NOT in the same row (different y-coordinates)
"""

import math
import sys


def main():
    input_file = 'initial_positions_xy.dat'

    print("="*70)
    print("Finding Probe Cells for Microrheology Test")
    print("="*70)

    # Read cell centers and apply 2× scaling
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

    # Initial box after scaling: x ∈ [0.5, 59.5], y ∈ [0.866, 51.095]
    # After compression to hexatic: x ∈ [-16.12, 16.12], y ∈ [-13.96, 13.96]
    # But we need to select based on INITIAL positions (before compression)

    # Filter cells near left edge and y-center
    # Left edge: x < 10.0 (after scaling)
    # Y-center: y ∈ [20, 30] (approximately center of initial box)
    candidates = [(mol, x, y) for mol, x, y in cells
                  if 1.0 < x < 10.0 and 20.0 < y < 30.0]

    print(f"Found {len(candidates)} candidates near left edge and y-center")

    if len(candidates) < 2:
        print("ERROR: Not enough candidates found!")
        sys.exit(1)

    # Find neighboring pairs
    pairs = []
    for i, (mol1, x1, y1) in enumerate(candidates):
        for mol2, x2, y2 in candidates[i+1:]:
            dist = math.sqrt((x2-x1)**2 + (y2-y1)**2)
            # Neighbors if distance ≈ 2.0 (lattice spacing after 2× scaling is ~2.0)
            if 1.8 < dist < 2.2:
                # Check they're not in same row (different y)
                if abs(y2 - y1) > 0.3:
                    pairs.append(((mol1, x1, y1), (mol2, x2, y2), dist))

    print(f"Found {len(pairs)} neighboring pairs with different y-coordinates")

    if not pairs:
        print("ERROR: No suitable neighboring pairs found!")
        print("\nRelaxing criteria - showing all pairs within distance 1.8-2.2:")
        for i, (mol1, x1, y1) in enumerate(candidates):
            for mol2, x2, y2 in candidates[i+1:]:
                dist = math.sqrt((x2-x1)**2 + (y2-y1)**2)
                if 1.8 < dist < 2.2:
                    dy = abs(y2-y1)
                    print(f"  Mol {mol1} ({x1:.2f}, {y1:.2f}) <-> Mol {mol2} ({x2:.2f}, {y2:.2f}), dist={dist:.3f}, dy={dy:.3f}")
        sys.exit(1)

    # Select best pair (closest to left edge)
    best = min(pairs, key=lambda p: (p[0][1] + p[1][1])/2)
    cell1, cell2, dist = best

    print("\n" + "="*70)
    print("SELECTED PROBE CELLS:")
    print("="*70)
    print(f"Cell 1: mol_id = {cell1[0]}, position = ({cell1[1]:.2f}, {cell1[2]:.2f})")
    print(f"Cell 2: mol_id = {cell2[0]}, position = ({cell2[1]:.2f}, {cell2[2]:.2f})")
    print(f"Distance: {dist:.3f} LJ")
    print(f"Y-difference: {abs(cell2[2] - cell1[2]):.3f} LJ (confirms different rows)")
    print("="*70)

    print("\n" + "="*70)
    print("LAMMPS COMMANDS TO ADD TO in.demo_hexatic:")
    print("="*70)
    print(f"""
# Microrheology test: Apply constant external force to probe cells
variable probe_mol_1 equal {cell1[0]}
variable probe_mol_2 equal {cell2[0]}

group probe_cell_1 molecule ${{probe_mol_1}}
group probe_cell_2 molecule ${{probe_mol_2}}

# Apply constant force in +x direction (0 degrees)
# Total force fd = 3.0, distributed over 50 beads per cell
variable fd equal 3.0
variable force_per_bead equal ${{fd}}/${{Nbeads_large}}  # fd/50 = 0.06

fix probe_force_1 probe_cell_1 addforce ${{force_per_bead}} 0.0 0.0
fix probe_force_2 probe_cell_2 addforce ${{force_per_bead}} 0.0 0.0

# Track probe cell positions
compute probe1_com probe_cell_1 com
compute probe2_com probe_cell_2 com
variable probe1_x equal c_probe1_com[1]
variable probe1_y equal c_probe1_com[2]
variable probe2_x equal c_probe2_com[1]
variable probe2_y equal c_probe2_com[2]
""")
    print("="*70)
    print("\nNOTE: Add these commands AFTER equilibration (line ~260),")
    print("      BEFORE the main dynamics run (line ~337)")
    print("="*70)

    return cell1[0], cell2[0]


if __name__ == "__main__":
    main()

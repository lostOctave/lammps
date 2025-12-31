#!/usr/bin/env python3
"""Quick verification: plot the first cell to visualize the structure."""

import math

# Read atoms for first cell (mol_id = 1, atoms 1-50)
with open('initial_cells.data', 'r') as f:
    lines = f.readlines()

# Find Atoms section
atoms_start = None
for i, line in enumerate(lines):
    if line.startswith('Atoms'):
        atoms_start = i + 2
        break

# Extract first cell atoms (lines for atoms 1-50)
cell_atoms = []
for i in range(50):
    parts = lines[atoms_start + i].split()
    atom_id = int(parts[0])
    mol_id = int(parts[1])
    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
    cell_atoms.append((atom_id, x, y))

# Calculate center
center_x = sum(x for _, x, y in cell_atoms) / len(cell_atoms)
center_y = sum(y for _, x, y in cell_atoms) / len(cell_atoms)

# Calculate distances from center
distances = []
for atom_id, x, y in cell_atoms:
    dist = math.sqrt((x - center_x)**2 + (y - center_y)**2)
    distances.append(dist)

# Calculate bond lengths (consecutive atoms)
bond_lengths = []
for i in range(len(cell_atoms)):
    x1, y1 = cell_atoms[i][1], cell_atoms[i][2]
    x2, y2 = cell_atoms[(i+1) % 50][1], cell_atoms[(i+1) % 50][2]
    bond_len = math.sqrt((x2-x1)**2 + (y2-y1)**2)
    bond_lengths.append(bond_len)

print("First Cell Verification:")
print(f"  Center: ({center_x:.4f}, {center_y:.4f})")
print(f"  Expected center: (12.0, 26.8468)")
print(f"  Radius (avg): {sum(distances)/len(distances):.6f}")
print(f"  Radius (min): {min(distances):.6f}")
print(f"  Radius (max): {max(distances):.6f}")
print(f"  Expected radius: 0.5")
print(f"  Bond length (avg): {sum(bond_lengths)/len(bond_lengths):.6f}")
print(f"  Expected bond length: {2*math.pi*0.5/50:.6f}")
print(f"  Number of atoms: {len(cell_atoms)}")
print(f"  Number of bonds: {len(bond_lengths)}")

# Check if all distances are ~0.5
radius_ok = all(abs(d - 0.5) < 1e-5 for d in distances)
print(f"\n  All atoms at radius 0.5: {'✓' if radius_ok else '✗'}")

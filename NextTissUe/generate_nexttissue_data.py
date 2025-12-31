#!/usr/bin/env python3
"""
generate_nexttissue_data.py

Generate LAMMPS initial data file for NextTissUe simulation.
Creates 900 cells in hexatic lattice with 50 beads per cell.

Usage:
    python3 generate_nexttissue_data.py <input_positions> <output_data>

Example:
    python3 generate_nexttissue_data.py initial_positions_xy.dat initial_cells.data
"""

import math
import sys
from collections import Counter


class LAMMPSDataGenerator:
    """Generator for NextTissUe LAMMPS data files."""

    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.scale_factor = 2.0
        self.cell_radius = 0.55  # Updated for φ ≈ 95% packing
        self.beads_per_cell = 50
        self.num_cells = 900

        self.cell_centers = []
        self.box_bounds = {}
        self.atoms = []  # List of (atom_id, mol_id, type, x, y, z)
        self.bonds = []  # List of (bond_id, bond_type, atom1, atom2)

    def read_and_scale_centers(self):
        """Read cell centers from file and apply 2× scaling."""
        print(f"Reading cell centers from {self.input_file}...")

        with open(self.input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                x, y = map(float, line.split())
                # Apply scaling
                self.cell_centers.append((x * self.scale_factor, y * self.scale_factor))

        print(f"  Read {len(self.cell_centers)} cell centers")

        # Calculate box bounds
        x_coords = [x for x, y in self.cell_centers]
        y_coords = [y for x, y in self.cell_centers]

        # Box needs to accommodate cell radius on all sides
        self.box_bounds = {
            'xlo': min(x_coords) - self.cell_radius,
            'xhi': max(x_coords) + self.cell_radius,
            'ylo': min(y_coords) - self.cell_radius,
            'yhi': max(y_coords) + self.cell_radius,
            'zlo': -0.5,
            'zhi': 0.5
        }

        print(f"  Box bounds: X=[{self.box_bounds['xlo']:.1f}, {self.box_bounds['xhi']:.10f}], "
              f"Y=[{self.box_bounds['ylo']:.1f}, {self.box_bounds['yhi']:.10f}]")

    def generate_bead_positions(self, cx, cy):
        """
        Generate positions of beads arranged in a circle.

        Args:
            cx, cy: Cell center coordinates

        Returns:
            List of (x, y, z) tuples for bead positions
        """
        beads = []

        # Angular spacing (negative for clockwise, matching template.cell_large)
        angular_spacing = -2.0 * math.pi / self.beads_per_cell

        for i in range(self.beads_per_cell):
            # Angle for bead i (starting at 0°, going clockwise)
            angle = i * angular_spacing

            # Calculate bead position using parametric circle equations
            bead_x = cx + self.cell_radius * math.cos(angle)
            bead_y = cy + self.cell_radius * math.sin(angle)
            bead_z = 0.0  # All beads in z=0 plane

            beads.append((bead_x, bead_y, bead_z))

        return beads

    def generate_all_atoms(self):
        """Generate all atom positions and IDs."""
        print("Generating atom positions...")

        atom_type = 1  # All atoms are type 1

        for cell_idx, (center_x, center_y) in enumerate(self.cell_centers):
            mol_id = cell_idx + 1  # Molecule IDs: 1 to 900

            # Generate bead positions for this cell
            bead_positions = self.generate_bead_positions(center_x, center_y)

            # Assign atom IDs and assemble atom data
            for bead_idx, (x, y, z) in enumerate(bead_positions):
                atom_id = cell_idx * self.beads_per_cell + bead_idx + 1
                self.atoms.append((atom_id, mol_id, atom_type, x, y, z))

        print(f"  Generated {len(self.atoms)} atoms")

    def generate_all_bonds(self):
        """Generate bond connectivity list."""
        print("Generating bond connectivity...")

        bond_id = 1
        bond_type = 1  # All bonds are type 1

        for cell_idx in range(len(self.cell_centers)):
            # First atom ID in this cell
            first_atom_in_cell = cell_idx * self.beads_per_cell + 1

            # Generate sequential bonds (1→2, 2→3, ..., 49→50)
            for bead_idx in range(self.beads_per_cell - 1):
                atom1 = first_atom_in_cell + bead_idx
                atom2 = first_atom_in_cell + bead_idx + 1
                self.bonds.append((bond_id, bond_type, atom1, atom2))
                bond_id += 1

            # Closure bond (50→1)
            atom1 = first_atom_in_cell + self.beads_per_cell - 1  # Last bead
            atom2 = first_atom_in_cell  # First bead
            self.bonds.append((bond_id, bond_type, atom1, atom2))
            bond_id += 1

        print(f"  Generated {len(self.bonds)} bonds")

    def write_data_file(self):
        """Write LAMMPS data file in molecular format."""
        print(f"Writing LAMMPS data file to {self.output_file}...")

        with open(self.output_file, 'w') as f:
            # Header comment
            f.write("LAMMPS data file for NextTissUe simulation\n")
            f.write("# 900 cells in hexatic lattice, 50 beads per cell\n\n")

            # Counts
            f.write(f"{len(self.atoms)} atoms\n")
            f.write(f"{len(self.bonds)} bonds\n")
            f.write(f"1 atom types\n")
            f.write(f"1 bond types\n\n")

            # Box bounds
            f.write(f"{self.box_bounds['xlo']:.1f} {self.box_bounds['xhi']:.10f} xlo xhi\n")
            f.write(f"{self.box_bounds['ylo']:.1f} {self.box_bounds['yhi']:.10f} ylo yhi\n")
            f.write(f"{self.box_bounds['zlo']:.1f} {self.box_bounds['zhi']:.1f} zlo zhi\n\n")

            # Masses section
            f.write("Masses\n\n")
            f.write("1 1.0\n\n")

            # Atoms section (molecular style: atom-ID mol-ID type x y z)
            f.write("Atoms # molecular\n\n")
            for atom_id, mol_id, atom_type, x, y, z in self.atoms:
                f.write(f"{atom_id} {mol_id} {atom_type} "
                       f"{x:.10f} {y:.10f} {z:.10f}\n")
            f.write("\n")

            # Bonds section
            f.write("Bonds\n\n")
            for bond_id, bond_type, atom1, atom2 in self.bonds:
                f.write(f"{bond_id} {bond_type} {atom1} {atom2}\n")

        print(f"  Successfully wrote {self.output_file}")

    def verify(self):
        """Verify generated data for correctness."""
        print("\nVerifying generated data...")

        errors = []

        # Check 1: Atom count
        expected_atoms = self.num_cells * self.beads_per_cell
        if len(self.atoms) != expected_atoms:
            errors.append(f"Expected {expected_atoms} atoms, got {len(self.atoms)}")
        else:
            print(f"  ✓ Atom count: {len(self.atoms)}")

        # Check 2: Bond count
        expected_bonds = self.num_cells * self.beads_per_cell
        if len(self.bonds) != expected_bonds:
            errors.append(f"Expected {expected_bonds} bonds, got {len(self.bonds)}")
        else:
            print(f"  ✓ Bond count: {len(self.bonds)}")

        # Check 3: Atom IDs are sequential and unique
        atom_ids = sorted([a[0] for a in self.atoms])
        if atom_ids != list(range(1, expected_atoms + 1)):
            errors.append("Atom IDs not sequential from 1 to " + str(expected_atoms))
        else:
            print(f"  ✓ Atom IDs sequential: 1 to {expected_atoms}")

        # Check 4: Molecule IDs range from 1 to num_cells
        mol_ids = set([a[1] for a in self.atoms])
        if mol_ids != set(range(1, self.num_cells + 1)):
            errors.append(f"Invalid molecule ID range")
        else:
            print(f"  ✓ Molecule IDs: 1 to {self.num_cells}")

        # Check 5: Each molecule has exactly beads_per_cell atoms
        mol_counts = Counter([a[1] for a in self.atoms])
        if not all(count == self.beads_per_cell for count in mol_counts.values()):
            errors.append("Not all molecules have 50 atoms")
        else:
            print(f"  ✓ Each molecule has {self.beads_per_cell} atoms")

        # Check 6: All atom types are 1
        atom_types = set([a[2] for a in self.atoms])
        if atom_types != {1}:
            errors.append(f"Expected only type 1, got {atom_types}")
        else:
            print(f"  ✓ All atoms are type 1")

        # Check 7: Bond IDs are sequential
        bond_ids = sorted([b[0] for b in self.bonds])
        if bond_ids != list(range(1, expected_bonds + 1)):
            errors.append("Bond IDs not sequential")
        else:
            print(f"  ✓ Bond IDs sequential: 1 to {expected_bonds}")

        # Check 8: All bond types are 1
        bond_types = set([b[1] for b in self.bonds])
        if bond_types != {1}:
            errors.append(f"Expected only bond type 1, got {bond_types}")
        else:
            print(f"  ✓ All bonds are type 1")

        # Check 9: All beads in a cell are at correct radius (sample first cell)
        if len(self.cell_centers) > 0:
            center_x, center_y = self.cell_centers[0]
            cell_atoms = [a for a in self.atoms if a[1] == 1]

            all_correct_radius = True
            for atom in cell_atoms:
                x, y, z = atom[3], atom[4], atom[5]
                distance = math.sqrt((x - center_x)**2 + (y - center_y)**2)
                if abs(distance - self.cell_radius) > 1e-6:
                    all_correct_radius = False
                    break

            if all_correct_radius:
                print(f"  ✓ Beads at correct radius (verified first cell)")
            else:
                errors.append("Beads not at correct radius from cell center")

        # Check 10: All z-coordinates are 0
        z_coords = set([a[5] for a in self.atoms])
        if z_coords != {0.0}:
            errors.append("Not all beads in z=0 plane")
        else:
            print(f"  ✓ All beads in z=0 plane")

        # Check 11: Verify closure bond exists in first cell
        first_cell_bonds = [b for b in self.bonds[:self.beads_per_cell]]
        closure_bond = [b for b in first_cell_bonds
                       if (b[2] == self.beads_per_cell and b[3] == 1) or
                          (b[2] == 1 and b[3] == self.beads_per_cell)]
        if len(closure_bond) == 1:
            print(f"  ✓ Closure bond verified (atom {self.beads_per_cell} → atom 1)")
        else:
            errors.append("Closure bond missing or duplicated")

        # Check 12: All atoms within box bounds (with small tolerance for floating point)
        tolerance = 1e-6
        out_of_bounds_atoms = []
        for a in self.atoms:
            if not (self.box_bounds['xlo'] - tolerance <= a[3] <= self.box_bounds['xhi'] + tolerance and
                    self.box_bounds['ylo'] - tolerance <= a[4] <= self.box_bounds['yhi'] + tolerance and
                    self.box_bounds['zlo'] - tolerance <= a[5] <= self.box_bounds['zhi'] + tolerance):
                out_of_bounds_atoms.append((a[0], a[3], a[4], a[5]))

        if len(out_of_bounds_atoms) == 0:
            print(f"  ✓ All atoms within box bounds")
        else:
            errors.append(f"Some atoms outside box bounds ({len(out_of_bounds_atoms)} atoms)")
            print(f"    Showing first 10 out-of-bounds atoms:")
            for atom_id, x, y, z in out_of_bounds_atoms[:10]:
                print(f"    Atom {atom_id}: x={x:.10f} (box: {self.box_bounds['xhi']:.10f}), "
                      f"y={y:.10f} (box: {self.box_bounds['yhi']:.10f})")

        # Print results
        if errors:
            print("\n❌ Verification FAILED:")
            for error in errors:
                print(f"  - {error}")
            return False
        else:
            print("\n✓ All verification checks passed!")
            return True

    def generate(self):
        """Main generation workflow."""
        print("="*60)
        print("NextTissUe LAMMPS Data File Generator")
        print("="*60)

        self.read_and_scale_centers()
        self.generate_all_atoms()
        self.generate_all_bonds()
        self.write_data_file()
        success = self.verify()

        print("="*60)
        if success:
            print("Generation completed successfully!")
        else:
            print("Generation completed with errors!")
        print("="*60)

        return success


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: python3 generate_nexttissue_data.py <input_positions> <output_data>")
        print("\nExample:")
        print("  python3 generate_nexttissue_data.py initial_positions_xy.dat initial_cells.data")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    generator = LAMMPSDataGenerator(input_file, output_file)
    success = generator.generate()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

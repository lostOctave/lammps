# NextTissUe Hexatic Lattice Setup

This directory contains files for initializing NextTissUe simulations with a hexatic lattice configuration.

## Files

### Data Generation
- **`generate_nexttissue_data.py`** - Python script to generate initial LAMMPS data file
- **`initial_positions_xy.dat`** - 900 cell centers in hexatic lattice (input)
- **`initial_cells.data`** - Generated LAMMPS data file with 45,000 atoms (output)

### LAMMPS Input Scripts
- **`in.demo`** - Original script (random initialization + cell replacement)
- **`in.demo_hexatic`** - New script using pre-generated hexatic lattice

## Quick Start

### 1. Generate Initial Data (already done)
```bash
python3 generate_nexttissue_data.py initial_positions_xy.dat initial_cells.data
```

### 2. Run Simulation
```bash
# Using the hexatic lattice configuration
lmp -in in.demo_hexatic

# OR using the original random configuration
lmp -in in.demo
```

## Key Differences: in.demo vs in.demo_hexatic

| Feature | in.demo | in.demo_hexatic |
|---------|---------|-----------------|
| Number of cells | 1000 | 900 |
| Initial config | Random particles | Pre-generated hexatic lattice |
| Initialization | Create → Minimize → Replace | Read data file |
| Lattice spacing | Variable (from minimization) | Fixed (1.0 units) |
| Box size | Calculated (31.6 × 31.6) | From data (60.0 × 51.6) |
| Setup time | ~minutes (minimize + replace) | Instant (read file) |

## Data File Specifications

**initial_cells.data** contains:
- **Atoms**: 45,000 (900 cells × 50 beads)
- **Bonds**: 45,000 (50 per cell, circular topology)
- **Lattice**: Hexatic with spacing 1.0 units
- **Cell radius**: 0.5 units (diameter = 1.0)
- **Cell arrangement**: Cells just touch (no overlap)
- **Box size**: 60.0 × 51.6 × 1.0 LJ units

## Verification

The generation script includes 12 verification checks:
```bash
python3 generate_nexttissue_data.py initial_positions_xy.dat initial_cells.data
```

Output shows:
- ✓ Atom count: 45,000
- ✓ Bond count: 45,000
- ✓ Atom IDs sequential: 1 to 45,000
- ✓ Molecule IDs: 1 to 900
- ✓ Each molecule has 50 atoms
- ✓ All atoms at correct radius
- ✓ Circular closure bonds verified
- ✓ All atoms within box bounds

## Regenerating Data

If you need to regenerate with different parameters, edit `generate_nexttissue_data.py`:

```python
class LAMMPSDataGenerator:
    def __init__(self, input_file, output_file):
        self.scale_factor = 2.0      # Change scaling
        self.cell_radius = 0.5       # Change cell size
        self.beads_per_cell = 50     # Change resolution
```

Then run:
```bash
python3 generate_nexttissue_data.py initial_positions_xy.dat initial_cells.data
```

## Notes

- The hexatic lattice preserves the structure from `initial_positions_xy.dat` with 2× scaling
- Original positions ranged [0.25, 29.75] → scaled to [0.5, 59.5]
- Box bounds include margin for cell radius: max_position + 0.5
- All beads lie in z=0 plane (2D simulation)
- Periodic boundary conditions applied in x, y directions

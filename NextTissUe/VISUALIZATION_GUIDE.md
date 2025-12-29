# LAMMPS Visualization Guide for NextTissUe

## Overview
The `in.demo` script generates output files in the `dumps/` directory. These are custom dump files containing atom positions, molecular information, and other properties.

## Output Files Generated

1. **mol_id.txt** - Initial molecular IDs
2. **compress_cells_propelled_*.dat** - Compression phase data (every 1000 steps)
3. **snapcells_thermo_*.dat** - Main simulation snapshots (every 100,000 steps)

## Visualization Options

### Option 1: OVITO (Recommended)
**OVITO** is the most user-friendly and popular tool for molecular dynamics visualization.

#### Installation:
```bash
# Download from https://www.ovito.org/
# Or install via snap (Linux):
sudo snap install ovito
```

#### Usage:
1. Open OVITO
2. File → Load File → Select your dump file (e.g., `dumps/snapcells_thermo_*.dat`)
3. OVITO will auto-detect the LAMMPS custom format
4. Use the timeline at the bottom to play through timesteps
5. Customize visualization:
   - Change particle colors by molecular ID
   - Adjust particle sizes
   - Add bonds/surfaces

### Option 2: VMD (Visual Molecular Dynamics)
**VMD** is powerful but more complex, great for detailed analysis.

#### Installation:
```bash
# Download from https://www.ks.uiuc.edu/Research/vmd/
# Ubuntu/Debian:
sudo apt-get install vmd
```

#### Usage:
1. Open VMD
2. File → New Molecule
3. Select your dump file
4. Set file type to "LAMMPS Trajectory"
5. Load and visualize

### Option 3: ParaView
**ParaView** is excellent for large datasets and scientific visualization.

#### Installation:
```bash
sudo apt-get install paraview
```

#### Usage:
1. Convert LAMMPS dump to VTK format first (using Pizza.py or custom script)
2. Load VTK files in ParaView

### Option 4: Python Visualization
For custom analysis and plotting, use Python:

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read LAMMPS dump file
def read_lammps_dump(filename):
    """Parse LAMMPS custom dump file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header and data
    # Implementation depends on your specific dump format
    pass

# Example: Plot cell positions
data = read_lammps_dump('dumps/snapcells_thermo_*.dat')
plt.scatter(data['x'], data['y'], c=data['mol'])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Cell Positions')
plt.axis('equal')
plt.show()
```

## Quick Start: Basic Python Visualization

Create a simple script to visualize your results:

```python
import numpy as np
import matplotlib.pyplot as plt
import glob

def parse_lammps_dump(filename, timestep=0):
    """Simple parser for LAMMPS custom dump files"""
    atoms = []
    reading_atoms = False
    atom_count = 0
    current_step = -1
    
    with open(filename, 'r') as f:
        for line in f:
            if 'ITEM: TIMESTEP' in line:
                current_step = int(next(f).strip())
                if current_step < timestep:
                    continue
                if current_step > timestep:
                    break
            elif 'ITEM: NUMBER OF ATOMS' in line:
                atom_count = int(next(f).strip())
            elif 'ITEM: ATOMS' in line:
                reading_atoms = True
                continue
            elif reading_atoms and current_step == timestep:
                parts = line.split()
                if len(parts) >= 4:  # id x y mol ...
                    atoms.append([float(x) for x in parts[:4]])
    
    return np.array(atoms)

# Visualize
files = glob.glob('dumps/*.dat')
if files:
    data = parse_lammps_dump(files[0])
    if len(data) > 0:
        plt.figure(figsize=(10, 10))
        plt.scatter(data[:, 1], data[:, 2], c=data[:, 3], s=20, cmap='tab20')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Cell Configuration')
        plt.axis('equal')
        plt.colorbar(label='Molecule ID')
        plt.savefig('cell_snapshot.png', dpi=150)
        print("Saved visualization to cell_snapshot.png")
```

## Tips

1. **For quick viewing**: Use OVITO - it's the easiest
2. **For publication-quality figures**: Use Python/matplotlib or OVITO's rendering
3. **For movies/animations**: OVITO can export MP4 videos directly
4. **For analysis**: Python with numpy/pandas is most flexible

## Checking Simulation Progress

```bash
# Monitor simulation in real-time
tail -f log.lammps

# Check dump files
ls -lh dumps/

# Quick check of output
head -n 20 dumps/*.dat
```

## Next Steps

1. Run the simulation (if not already running)
2. Monitor `log.lammps` for progress
3. Once dump files are generated, use OVITO for quick visualization
4. For detailed analysis, write custom Python scripts

## Installation Commands Summary

```bash
# OVITO (easiest)
sudo snap install ovito

# VMD
sudo apt-get install vmd

# Python dependencies for custom visualization
pip install numpy matplotlib pandas
```


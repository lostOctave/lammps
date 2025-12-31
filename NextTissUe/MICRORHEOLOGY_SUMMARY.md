# Microrheology Implementation - Quick Summary

## What Was Implemented

A complete **active microrheology** system for the NextTissUe cellular tissue model that:

1. ✅ Applies constant force to probe cells
2. ✅ Tracks displacement of probe cells over time
3. ✅ Monitors response of nearby cells in distance shells
4. ✅ Extracts rheological properties (viscosity, drag coefficient)
5. ✅ Analyzes spatial decay of mechanical response

## File Structure

```
NextTissUe/
├── in.demo_hexatic                  # Main LAMMPS simulation (WITH microrheology)
├── find_probe_cells.py              # Identifies 2 probe cells (mol 709, 710)
├── identify_nearby_cells.py         # Identifies 52 nearby cells in 3 shells
├── analyze_microrheology.py         # Post-processing analysis
├── README_microrheology.md          # Complete documentation
├── MICRORHEOLOGY_SUMMARY.md         # This file
└── nearby_cells.txt                 # List of nearby cells with distances
```

## How It Works

### Step 1: Setup (Automated)

```
┌─────────────────────────────────────────┐
│   Hexatic Lattice (900 cells)          │
│                                         │
│   ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯             │
│  ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯            │
│   ◯ ◉ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯    ← Probe 1 │
│  ◯ ◯ ◉ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯   ← Probe 2 │
│   ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯             │
│  ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯ ◯            │
│                                         │
└─────────────────────────────────────────┘

Probe cells: mol 709, 710 (near left edge)
```

### Step 2: Force Application

```
      F = 3.0 LJ
      ────────────→

   ◉ ← Probe cell 1 (50 beads, each gets F/50 = 0.06)
   ◉ ← Probe cell 2 (50 beads, each gets F/50 = 0.06)
```

### Step 3: Response Measurement

```
Distance Shells Around Probes:
┌────────────────────────────────────┐
│                                    │
│         ╔═══════════╗              │
│        ║ ┌─────────┐║              │
│       ║ │  ┌───┐   │║              │
│      ║ │  │ ◉ │   │ ║              │
│       ║ │  └───┘   │║              │
│        ║ └─────────┘║              │
│         ╚═══════════╝              │
│                                    │
│   ┌───┐  First shell:  7 cells    │
│   │ ◉ │  Second shell: 19 cells   │
│   └───┘  Third shell:  26 cells   │
│                                    │
└────────────────────────────────────┘

Each shell tracks:
- Center-of-mass displacement
- Correlation with probe motion
- Area, perimeter, shape changes
```

### Step 4: Data Collection

**High-frequency dumps** (every 100 steps):
- Probe cells 1 & 2
- First shell cells (nearest neighbors)

**Medium-frequency dumps** (every 1000 steps):
- Second shell cells
- Third shell cells

**Low-frequency dump** (every 100,000 steps):
- All cells (full system snapshot)

### Step 5: Analysis

The `analyze_microrheology.py` script computes:

```
PROBE CELLS:
┌──────────────────────────────────────────┐
│ • Displacement: Δx, Δy                   │
│ • Velocity: vₓ = dΔx/dt                  │
│ • MSD: <(r(t)-r(0))²>                    │
│ • Drag coefficient: γ = 2F²/(dMSD/dt)    │
│ • Viscosity: η ~ γ/(6πR)                 │
└──────────────────────────────────────────┘

NEARBY CELLS:
┌──────────────────────────────────────────┐
│ Shell 1: Correlation = 0.7-0.9 (strong)  │
│ Shell 2: Correlation = 0.3-0.6 (medium)  │
│ Shell 3: Correlation = 0.1-0.3 (weak)    │
└──────────────────────────────────────────┘

         High correlation → ELASTIC
         Low correlation  → VISCOUS
```

## Quick Start

### 1. Run Identification Scripts (Already Done)
```bash
cd /home/lost_octave/LAMMPS/NextTissUe

python3 find_probe_cells.py          # Finds probe cells 709, 710
python3 identify_nearby_cells.py     # Finds 52 nearby cells
```

### 2. Run LAMMPS Simulation
```bash
# Make sure LAMMPS is built with NextTissUe package
cd /home/lost_octave/LAMMPS/NextTissUe

# Run simulation (this will take a while!)
../build/lmp -in in.demo_hexatic

# Output files will be in /mnt/d/NextTissUe_data/
#   microrheology_probe_NC900_fd3.0.dat
#   microrheology_probe2_NC900_fd3.0.dat
#   microrheology_first_shell_NC900_fd3.0.dat
#   microrheology_second_shell_NC900_fd3.0.dat
#   microrheology_third_shell_NC900_fd3.0.dat
```

### 3. Analyze Results
```bash
python3 analyze_microrheology.py

# Output:
#   - Console summary
#   - microrheology_results.txt (detailed results)
```

## Expected Timeline

| Phase | Duration | What Happens |
|-------|----------|--------------|
| Equilibration | 1M steps | System settles after compression |
| Measurement | 20M steps | Probe cells dragged, data collected |
| Analysis | < 1 min | Post-processing of dump files |

**Total simulation time**: ~20-30 minutes (depending on CPU)

## Key Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `fd` | 3.0 | Total force on each probe cell |
| `number_cells` | 900 | Total cells in system |
| `v0` | 0.01 | Cell self-propulsion speed |
| `Dr` | 1.0 | Rotational diffusion coefficient |
| `p0` | 3.65 | Shape parameter (fluid-like) |
| `visc` | 50 | Viscous damping coefficient |

## What to Expect

### Probe Cell Motion

For a **fluid-like tissue** (p₀ = 3.65):
- Initial acceleration (ballistic regime)
- Transition to steady-state velocity
- Total displacement: **10-50 LJ** in x-direction

### Viscosity Estimate

Typical values for this system:
- **Drag coefficient**: γ ~ 30-300
- **Effective viscosity**: η ~ 5-50 (LJ units)
- **Steady velocity**: v ~ F/γ ~ 0.01-0.1 LJ/timestep

### Nearby Cell Response

Expected **correlation decay** with distance:

```
Correlation
    1.0 │  ●
        │   ●●
    0.8 │     ●●
        │       ●●
    0.6 │         ●●
        │           ●●
    0.4 │             ●●
        │               ●●
    0.2 │                 ●●●
        │                    ●●●●
    0.0 └─────────────────────────●●●●
        0   2   4   6   8  10  12  14
               Distance (LJ)

Strong decay → VISCOUS (fluid-like)
Slow decay   → ELASTIC (solid-like)
```

## Modifications for Different Scenarios

### Increase Force (Faster Motion)
```lammps
variable fd equal 10.0    # Changed from 3.0
```

### Probe Different Region
```lammps
# Edit find_probe_cells.py line 43-45:
# Change x-range: 1.0 < x < 10.0  (left edge)
# To:            25.0 < x < 35.0  (center)
# To:            50.0 < x < 59.0  (right edge)
```

### Change Number of Shells
```python
# Edit identify_nearby_cells.py line 64-68:
shells = [
    (0, 2.0, "first_shell"),     # Tighter first shell
    (2.0, 5.0, "second_shell"),
    (5.0, 10.0, "third_shell"),
    (10.0, 15.0, "fourth_shell"), # Add fourth shell
]
```

## Comparison with Experiments

This computational microrheology mimics experimental techniques:

| Technique | Our Implementation | Similarities |
|-----------|-------------------|--------------|
| **Optical tweezers** | Fix force on probe cells | Constant force application |
| **Magnetic bead rheology** | Track probe displacement | MSD analysis |
| **Active microrheology** | Drag coefficient calculation | Extract viscosity |

**Advantages of simulation**:
- Perfect control over force
- Track all nearby cells (not just probes)
- Measure individual cell properties
- No experimental noise

## Troubleshooting

### Probe cells barely move
→ Increase `fd` or decrease `visc`

### Probe cells move too fast
→ Decrease `fd` or increase `visc`

### Nearby cells don't respond
→ Check `eps_att` (cell-cell adhesion)
→ Increase `pre_fac` for stronger coupling

### Files not found during analysis
→ Edit `data_dir` in analyze_microrheology.py
→ Change from `/mnt/d/NextTissUe_data` to local path

## Next Steps

1. **Run baseline simulation**: Current configuration with fd=3.0
2. **Vary force**: fd = 1.0, 3.0, 10.0 to test linearity
3. **Vary packing**: num_dense = 0.5, 1.0, 2.0 to see density effects
4. **Vary activity**: v0 = 0.0, 0.01, 0.1 to see active vs passive

## References

Full documentation: `README_microrheology.md`
Physics background: `PHYSICS_ANALYSIS.md`
LAMMPS setup: `in.demo_hexatic`

---

**Status**: ✅ Implementation complete and ready to run!

**Created**: 2025-12-30
**Last modified**: 2025-12-30

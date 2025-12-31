# Microrheology Implementation for NextTissUe Model

## Overview

This implementation adds **active microrheology** capabilities to the NextTissUe cellular tissue model. Microrheology measures the mechanical properties (viscosity, elasticity) of a material by tracking the motion of probe particles under applied forces.

## Physical Concept

### What is Microrheology?

Microrheology is a technique to measure the rheological (flow) properties of materials by:
1. **Applying a force** to one or more probe particles
2. **Tracking the displacement** of the probe particles over time
3. **Measuring the response** of nearby particles/cells
4. **Extracting material properties** like viscosity, elastic modulus, etc.

### Active vs Passive Microrheology

- **Passive**: Track thermal fluctuations (Brownian motion) without external force
- **Active** (implemented here): Apply constant external force and measure response

## Implementation Details

### 1. Probe Cell Selection

**Script**: `find_probe_cells.py`

- Identifies two neighboring cells near the left edge of the simulation box
- Cells are selected at the vertical center to minimize boundary effects
- Selected cells (for 900-cell system):
  - **Cell 709**: position (1.50, 26.85) LJ
  - **Cell 710**: position (2.50, 28.58) LJ
  - **Distance**: 2.00 LJ (nearest neighbors in hexatic lattice)

**Usage**:
```bash
python3 find_probe_cells.py
```

### 2. Nearby Cell Identification

**Script**: `identify_nearby_cells.py`

Identifies cells in distance shells around probe cells:

| Shell | Distance Range | Number of Cells | Purpose |
|-------|----------------|-----------------|---------|
| First | 0-3 LJ | 7 | Nearest neighbors (direct contact) |
| Second | 3-6 LJ | 19 | Second shell (indirect contact) |
| Third | 6-9 LJ | 26 | Far field (long-range response) |

**Usage**:
```bash
python3 identify_nearby_cells.py
```

**Output**:
- Console output with LAMMPS commands
- `nearby_cells.txt`: List of all nearby cells with distances

### 3. LAMMPS Simulation Setup

**File**: `in.demo_hexatic`

The simulation applies a **constant force** to the probe cells:

```lammps
# Force parameters
variable fd equal 3.0                        # Total force per cell
variable force_per_bead equal ${fd}/${Nbeads_large}  # Force per bead = 0.06

# Apply force to probe cells
fix probe_force_1 probe_cell_1 addforce ${force_per_bead} 0.0 0.0
fix probe_force_2 probe_cell_2 addforce ${force_per_bead} 0.0 0.0
```

### Force Direction

- **Direction**: +x (horizontal, rightward)
- **Magnitude**: 3.0 LJ force units total
  - Distributed over 50 beads → 0.06 per bead
  - This is ~6× the characteristic force scale (ε_wca = 0.02 per bead)

### 4. Data Collection

The simulation outputs specialized dump files at different frequencies:

| Dump | Frequency | Content | Purpose |
|------|-----------|---------|---------|
| `microrheology_probe_NC900_fd3.0.dat` | Every 100 steps | Probe cell 1 beads | High-res probe tracking |
| `microrheology_probe2_NC900_fd3.0.dat` | Every 100 steps | Probe cell 2 beads | High-res probe tracking |
| `microrheology_first_shell_NC900_fd3.0.dat` | Every 100 steps | First shell cells | Nearest neighbor response |
| `microrheology_second_shell_NC900_fd3.0.dat` | Every 1000 steps | Second shell cells | Medium-range response |
| `microrheology_third_shell_NC900_fd3.0.dat` | Every 1000 steps | Third shell cells | Long-range response |

**Dump columns**:
- `id`: Atom ID
- `mol`: Molecule (cell) ID
- `xu`, `yu`: Unwrapped coordinates (accounts for periodic boundaries)
- `d_xc`, `d_yc`: Cell center-of-mass coordinates
- `fx`, `fy`: Total force on atom
- `d_area`, `d_peri`, `d_shape`: Cell properties (area, perimeter, shape index)

### 5. Analysis

**Script**: `analyze_microrheology.py`

Performs comprehensive analysis of microrheology data:

#### Probe Cell Analysis

1. **Center-of-Mass Trajectory**
   - Computes COM from all 50 beads of each probe cell
   - Tracks displacement in x (force direction) and y (perpendicular)

2. **Velocity Calculation**
   - Computes instantaneous velocity: v = dr/dt
   - Mean velocity in force direction
   - Perpendicular velocity (should be ~0 for pure drag)

3. **Mean Square Displacement (MSD)**
   ```
   MSD(t) = <[r(t) - r(0)]²>
   ```
   - Measures how far probe has moved from initial position
   - For driven particle: MSD ~ (F/γ)² t² (ballistic) → 2F²/γ × t (diffusive)

4. **Drag Coefficient Estimation**
   - In viscous regime: F = γv (Stokes drag)
   - From MSD slope: γ ~ 2F² / [d(MSD)/dt]
   - Effective viscosity: η ~ γ / (6πR)

#### Nearby Cell Analysis

1. **Displacement Correlation**
   - Compute COM displacement of each shell
   - Correlation coefficient with probe displacement
   - High correlation → cells move together (elastic response)
   - Low correlation → cells independent (viscous dissipation)

2. **Spatial Decay**
   - How response decreases with distance
   - Exponential decay suggests characteristic length scale
   - Power-law decay suggests long-range elasticity

**Usage**:
```bash
python3 analyze_microrheology.py
```

**Output**:
- Console summary of results
- `microrheology_results.txt`: Detailed numerical results

## Expected Results

### Probe Cell Motion

For a **fluid-like tissue** (p₀ = 3.65 < 3.81):
- **Short time** (0-1000 steps): Accelerating motion, MSD ~ t²
- **Long time** (>10000 steps): Steady-state velocity, MSD ~ t

**Typical values** (estimated):
- Displacement after 1M steps: ~10-50 LJ (depending on viscosity)
- Steady-state velocity: ~0.01-0.1 LJ/timestep
- Drag coefficient: γ ~ 30-300 (depending on cell packing)

### Nearby Cell Response

Expected **correlation with probe displacement**:

| Shell | Distance | Expected Correlation | Interpretation |
|-------|----------|---------------------|----------------|
| First | 0-3 LJ | 0.7-0.9 | Strong coupling (direct contact) |
| Second | 3-6 LJ | 0.3-0.6 | Moderate coupling (indirect) |
| Third | 6-9 LJ | 0.1-0.3 | Weak coupling (far field) |

**If correlations are higher**: Material is more **elastic** (solid-like)
**If correlations are lower**: Material is more **viscous** (fluid-like)

## Physical Parameters

### Force Scale

Applied force `fd = 3.0` compared to characteristic forces:
- **Cell-cell repulsion**: F_rep ~ ε_wca/σ ~ 0.02/0.069 ~ 0.3 per bead pair
- **Cell-cell attraction**: F_att ~ ε_att/σ ~ 0.2/0.069 ~ 3.0 per bead pair
- **Applied force**: F_applied = 3.0/50 = 0.06 per bead

→ Applied force is **moderate** compared to cell-cell interactions

### Time Scales

- **Timestep**: dt = 0.001
- **Rotational diffusion**: τ_rot ~ 1/Dr = 1/1.0 = 1.0 timestep
- **Viscous relaxation**: τ_visc = M/γ_visc = 1.0/50 = 0.02 timestep
- **Active persistence**: τ_active = 1/Dr = 1000 timesteps

**Total simulation**: 20M steps = 20,000 time units
- Enough for ~20,000 rotational diffusion times
- Sufficient to reach steady state

### Length Scales

- **Cell diameter**: d ≈ 1.0 LJ
- **Lattice spacing** (after compression): a ≈ 1.0 LJ
- **Probe displacement** (expected): Δx ~ 10-50 LJ
- **Box size** (hexatic): Lx ≈ 32 LJ, Ly ≈ 28 LJ

## Running the Full Workflow

### 1. Generate Initial Configuration
```bash
python3 generate_nexttissue_data.py
```

### 2. Identify Probe and Nearby Cells
```bash
python3 find_probe_cells.py
python3 identify_nearby_cells.py
```

### 3. Run LAMMPS Simulation
```bash
# Build LAMMPS with NextTissUe package (see main README)
cd /home/lost_octave/LAMMPS
mkdir build
cmake -S cmake -B build -C cmake/presets/basic.cmake -D PKG_MOLECULE=on
cmake --build build -j 4

# Run simulation
cd NextTissUe
../build/lmp -in in.demo_hexatic
```

### 4. Analyze Results
```bash
python3 analyze_microrheology.py
```

## Troubleshooting

### Issue: Probe cells don't move

**Possible causes**:
1. Force too small compared to cell-cell adhesion
   - **Solution**: Increase `fd` variable in `in.demo_hexatic`

2. System is jammed (solid-like)
   - **Solution**: Decrease p₀ (shape parameter) or num_dense (packing fraction)

3. Viscous damping too high
   - **Solution**: Decrease `visc` variable

### Issue: Probe cells move too fast

**Possible causes**:
1. Force too large
   - **Solution**: Decrease `fd`

2. Viscous damping too low
   - **Solution**: Increase `visc`

### Issue: Analysis script fails

**Possible causes**:
1. Dump files not generated
   - **Solution**: Check that simulation has run long enough (>100,000 steps)
   - Check data directory path: `/mnt/d/NextTissUe_data`

2. Wrong file paths
   - **Solution**: Edit `data_dir` in `analyze_microrheology.py`

## References

### Microrheology Theory

1. **Mason, T. G., & Weitz, D. A.** (1995). Optical measurements of frequency-dependent linear viscoelastic moduli of complex fluids. *Physical Review Letters*, 74(7), 1250.

2. **Squires, T. M., & Mason, T. G.** (2010). Fluid mechanics of microrheology. *Annual Review of Fluid Mechanics*, 42, 413-438.

### Cellular Tissue Rheology

3. **Bi, D., et al.** (2015). A density-independent rigidity transition in biological tissues. *Nature Physics*, 11(12), 1074-1079.

4. **Atia, L., et al.** (2018). Geometric constraints during epithelial jamming. *Nature Physics*, 14(6), 613-620.

### Vertex Model

5. **Farhadifar, R., et al.** (2007). The influence of cell mechanics, cell-cell interactions, and proliferation on epithelial packing. *Current Biology*, 17(24), 2095-2104.

## Future Extensions

### 1. Oscillatory Rheology
- Apply sinusoidal force: F(t) = F₀ sin(ωt)
- Measure phase lag between force and displacement
- Extract storage modulus G'(ω) and loss modulus G"(ω)

### 2. Multiple Probe Sizes
- Use cells of different sizes as probes
- Study size-dependent drag (test Stokes law)

### 3. Two-Point Microrheology
- Track correlation between two probe cells
- Extract spatial correlation function
- Measure characteristic length scale of mechanical response

### 4. Active vs Passive Comparison
- Turn off active forces (v₀ = 0) during measurement
- Compare rheology of active vs passive tissue
- Study effect of cell motility on tissue fluidity

## Contact

For questions or issues with this implementation, please check:
- Main repository README: `/home/lost_octave/LAMMPS/CLAUDE.md`
- LAMMPS documentation: https://docs.lammps.org
- NextTissUe model physics: `PHYSICS_ANALYSIS.md`

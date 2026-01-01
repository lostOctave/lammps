# NextTissUe Documentation Index for Claude Code

This file helps Claude navigate NextTissUe package documentation.

## Quick Reference

- **Overview**: [README.md](README.md) - Package description and basic usage
- **Working Notes**: [WORKING_NOTES.md](WORKING_NOTES.md) - User preferences and collaboration history
- **Microrheology**: [MICRORHEOLOGY_SUMMARY.md](MICRORHEOLOGY_SUMMARY.md) - External force application setup
- **Physics Background**: [PHYSICS_ANALYSIS.md](PHYSICS_ANALYSIS.md) - Physical meaning of parameters
- **Simplified Running**: [SIMPLIFIED_RUN.md](SIMPLIFIED_RUN.md) - How to run simulations

## Reusable Patterns (Portable Documentation)

These patterns are **self-contained** and can be copied to other projects:

- **[PATTERN_TwoLineProgressBar.md](PATTERN_TwoLineProgressBar.md)** - Live progress display with ANSI escape codes, real-time subprocess monitoring, cumulative ETA calculation
- **[PATTERN_ResumeSkipByDefault.md](PATTERN_ResumeSkipByDefault.md)** - Safe resumable batch jobs with automatic skip-by-default behavior

These documents capture complete implementation details, design evolution, edge cases, and ready-to-use code snippets. Copy them to any project where you need these features.

## Parameter Reference

### Sweep Parameters (varied in campaigns)
- **`karea`**: Area stiffness (incompressibility) - Values: [500, 1000, 2000, 5000]
- **`pre_fac`**: Pair attraction strength scaling - Values: [5, 10, 20, 50]
- **`fd`**: External probe force magnitude - Values: [200, 500, 1000, 2000]

### Cell Mechanics Parameters (fixed in input file)
- **`p0`**: Target shape parameter (perimeter/√area) - Controls cell roundness/bending
- **`Dr`**: Rotational diffusion coefficient - Controls propulsion direction noise
- **`v0`/`magnitude`**: Self-propulsion magnitude - Active force strength
- **`tau_c`/`damp`**: Contact inhibition timescale - Propulsion reorientation rate

### Global Simulation Parameters (configurable)
- **`main_run_steps`**: Main measurement phase duration - Default: 500,000 (user modified)
- **`dump_interval`**: Trajectory output frequency - Default: 5,000 (user modified)

## Key Files

### Input/Configuration
- **`in.demo_hexatic`**: Main LAMMPS input file (362 lines)
- **`sweep_config.yaml`**: Parameter sweep configuration
- **`run_parameter_sweep.py`**: Automated sweep script (627 lines)

### Initial Setup
- **`initial_cells.data`**: Pre-generated initial configuration (3.2 MB)
- **`random_propelling_directions.x`**: Executable for orientation initialization

### Analysis Tools
- **`analyze_microrheology.py`**: MSD and viscosity calculation
- **`find_probe_cells.py`**: Identify cells for force application
- **`identify_nearby_cells.py`**: Map distance shells around probe

## Performance Notes

### Bottlenecks (from analysis)
1. **MPI synchronization**: 13 MPI_Allreduce per timestep in bond/fix computations (~50-250s overhead)
2. **Double bond pass**: Sequential harmonic + area force loops (~25-50s)
3. **Custom property lookups**: Repeated string hash searches (~5-10s)

### User's Microrheology Code
- **Overhead**: ~2.8 seconds per 500k timesteps (~0.5% of total runtime)
- **Force application**: O(50) atoms per timestep (negligible)
- **Shell tracking**: O(2,600) atoms per timestep (acceptable)
- **Verdict**: Implementation is optimal, not the source of slowness

## User Preferences (from Working Notes)

- **Performance conscious**: Always asks about overhead before accepting changes
- **Prefers compact displays**: Uses k/M suffixes (e.g., "500k" not "500000")
- **Iterative refinement**: Tests immediately in real environment, adjusts based on results
- **Direct communication**: Appreciates concise, focused explanations
- **Values stability**: Prefers cumulative calculations over instant measurements

## Recent Work

### Parameter Sweep System (Complete)
- Created automated sweep script with live progress bars
- YAML configuration for campaign definitions
- Two-line progress display with cumulative ETA
- Resume capability (now default, use `--force-rerun` to override)

### Microrheology Implementation (Complete)
- Added external force to specific cells via `addforce` fix
- Probe cell selection by molecule ID (cell 709)
- Distance shell tracking (7 nearest, 19 second, 26 third neighbors)
- COM position tracking for response measurement

## Common Commands

```bash
# Run parameter sweep (automatically skips completed)
python3 run_parameter_sweep.py --all

# Force re-run all
python3 run_parameter_sweep.py --all --force-rerun

# Quick test run
python3 run_parameter_sweep.py --campaign 1 --main-run-steps 10000

# Build LAMMPS (from repo root)
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake
cmake --build build -j 4

# Run single simulation
cd NextTissUe
mpirun -np 12 ../build/lmp -in in.demo_hexatic -var karea 1000 -var pre_fac 10 -var fd 500
```

## File Locations

- **LAMMPS repo root**: `/home/lost_octave/LAMMPS/`
- **NextTissUe package**: `/home/lost_octave/LAMMPS/NextTissUe/`
- **Build directory**: `/home/lost_octave/LAMMPS/build/`
- **Output data**: `/mnt/d/NextTissUe_data/run/`
- **Source code**: `/home/lost_octave/LAMMPS/src/USER-AREABOND/`

## Integration Points

### NextTissUe ↔ USER-AREABOND
- **bond_harmonic_area.cpp**: Area constraint forces, cell geometry computation
- **fix_propel_cell.cpp**: Self-propulsion with contact inhibition
- **pair_wca_att.cpp**: WCA repulsion + attractive cosine-squared potential

### Parameter Flow
```
sweep_config.yaml → run_parameter_sweep.py → LAMMPS -var flags → in.demo_hexatic
                                                                    ↓
                                                            C++ force classes
                                                                    ↓
                                                    Output: trajectories, logs
```

---

Last updated: 2026-01-01

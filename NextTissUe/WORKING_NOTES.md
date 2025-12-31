# NextTissUe Parameter Sweep - Working Notes

## What We Built

### Core System
- **Parameter sweep automation** for LAMMPS simulations
- **YAML configuration** support for easy parameter changes
- **Live progress tracking** with two-line display (per-run + overall)
- **Cumulative ETA calculation** for stable time estimates
- **Resume capability** with COMPLETE markers

### Key Files
1. `run_parameter_sweep.py` - Main automation script (399 lines → ~630 lines)
2. `sweep_config.yaml` - Campaign configuration
3. `in.demo_hexatic` - Modified LAMMPS input (thermo 5000 for live updates)

### Features Implemented
- ✅ Multi-campaign parameter sweeps (32 runs total)
- ✅ MPI parallelism within runs (12 cores per simulation)
- ✅ Sequential execution between runs
- ✅ Live progress bar updates every ~18 seconds
- ✅ Two-line display: per-run progress + overall sweep status
- ✅ Compact formatting (15k/500k instead of 15000/500000)
- ✅ Cumulative ETA calculation (not instant speed)
- ✅ ANSI escape codes for in-place updates
- ✅ Fallback `--simple-progress` mode for terminals without ANSI support

## User Preferences & Habits

### Communication Style
- **Direct and efficient** - wants clear answers, not over-explanation
- **Practical focus** - cares about what works, not theory
- **Iterative refinement** - starts with basic idea, refines through testing
- **Visual feedback** - appreciates examples showing actual output

### Technical Preferences
1. **Performance conscious** - always asks about overhead/performance impact
2. **Space-saving** - prefers compact displays (e.g., "15k/500k" not "15000/500000")
3. **Real-time visibility** - wants to see progress live, not wait for completion
4. **Cumulative statistics** - prefers stable averages over instant measurements
5. **Clean displays** - likes organized, readable output

### Working Style
- **"option3"** - gives terse responses when decision is clear
- **Tests in real environment** - runs code immediately, reports issues
- **Notices edge cases** - "the problem arise as I expected" (predicted ANSI issues)
- **Asks for improvements** - "could it be possible..." when has an idea
- **Values stability** - "estimate better" = wants robust calculations

### Language Notes
- Chinese speaker (用了 "进度条" = progress bar)
- Excellent English communication
- Prefers concise technical terms

## Key Lessons Learned

### What Worked Well
1. **Incremental development** - Built features one at a time, tested each
2. **Multiple options** - Provided both ANSI and simple progress modes
3. **Performance analysis** - Explained overhead upfront (< 0.1% CPU)
4. **Visual examples** - Showed what output looks like before implementing
5. **Quick fixes** - When issues arose (ANSI lines), immediately offered solutions

### Technical Decisions
1. **LAMMPS thermo frequency**: Changed 100000 → 5000 for better visualization
2. **Progress threshold**: Lowered 10% → 0.5% so ETA shows early
3. **ETA calculation**: Cumulative progress (stable) not instant speed (fluctuates)
4. **Display format**: Two lines with ANSI codes (clean, readable)
5. **Number formatting**: 15k/500k saves 5 characters, easier to read

### Problem-Solving Pattern
- User reports issue → I provide explanation + solution + alternatives
- Example: "double line keep generating lines"
  - Diagnosed: ANSI codes not working
  - Fixed: Improved escape sequence logic
  - Fallback: Added `--simple-progress` flag
  - Result: User chose fixed version, it worked!

## Common Issues & Solutions

### LAMMPS Variable Styles
- **Issue**: Cannot mix equal-style (in file) with index-style (command-line)
- **Solution**: Comment out variable definitions in input file
- **User learned**: Difference between `-var` and `variable equal`

### Progress Display
- **Issue**: Terminal might not support ANSI escape codes
- **Solution**: Provide fallback with `--simple-progress`
- **Lesson**: Always have a simple fallback for terminal features

### ETA Calculation
- **Issue**: Instant speed fluctuates, gives bad estimates
- **Solution**: Use cumulative progress for stability
- **User insight**: "calculate accumulatively, better estimation"

## Quick Reference for Future Work

### Running Sweeps
```bash
# With YAML config (recommended)
python3 run_parameter_sweep.py --config sweep_config.yaml --campaign 1

# Resume interrupted sweep
python3 run_parameter_sweep.py --config sweep_config.yaml --all --resume

# Simple progress (no ANSI)
python3 run_parameter_sweep.py --config sweep_config.yaml --campaign 1 --simple-progress
```

### Modifying Parameters
Edit `sweep_config.yaml`:
- Campaign definitions (karea, pre_fac, fd values)
- Simulation settings (main_run_steps, dump_interval)
- Execution config (mpi_processes, lammps_exe)

### File Locations
- Input template: `/home/lost_octave/LAMMPS/NextTissUe/in.demo_hexatic`
- Output dirs: `/mnt/d/NextTissUe_data/run/karea{X}_pre_fac{Y}_fd{Z}/`
- Logs: `sweep.log`, `sweep_manifest.csv`

## Future Collaboration Tips

1. **Show, don't tell** - User likes seeing actual output examples
2. **Explain performance** - Always address CPU/memory overhead
3. **Provide options** - Main solution + fallback when uncertain
4. **Test assumptions** - User will test immediately, be ready for edge cases
5. **Concise summaries** - Brief "what we did + how to use it" over long docs
6. **Respect experience** - User predicted ANSI issue, knows their environment
7. **Iterate quickly** - Small fixes → test → adjust → done

---
*Last updated: 2026-01-01*
*Working relationship: Efficient, practical, iterative*

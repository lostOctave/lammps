# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics simulation code designed for parallel computers. This is a large, mature C++ codebase (~600MB, ~4,000 C++ files in src/) maintained by an international team of developers lead by staff at Sandia National Laboratories as open-source software under GPL v2.

**Primary Languages:** C++17 (core), C, Fortran, Python (interfaces)
**Build Systems:** CMake (primary, modern), Make (traditional, still supported)
**Key Frameworks:** MPI (parallel execution), OpenMP (threading), Kokkos (performance portability)

## Critical Build System Notes

### CMake Configuration (REQUIRED READING)

**LAMMPS uses a non-standard CMake layout.** The CMakeLists.txt is in the `cmake/` directory, NOT the repository root.

**ALWAYS use:**
```bash
cmake -S cmake -B build    # NOT: cmake -S . -B build
```

**Standard build sequence:**
```bash
mkdir build
cmake -S cmake -B build -C cmake/presets/basic.cmake
cmake --build build -j 4
# Executable: build/lmp
```

**Development build with most packages:**
```bash
cmake -S cmake -B build \
      -C cmake/presets/gcc.cmake \
      -C cmake/presets/most.cmake \
      -D ENABLE_TESTING=on \
      -D DOWNLOAD_POTENTIALS=off
cmake --build build -j 4
```

**IMPORTANT:** Use `-D DOWNLOAD_POTENTIALS=off` by default to avoid network dependency issues.

### Switching Build Systems

**Critical:** Never mix build systems without cleanup:
- **Make → CMake:** Run `make -C src purge` first
- **CMake → Make:** Run `make -C src clean-all` first

### Traditional Make Build (Legacy)

```bash
cd src
make serial     # Serial build (no MPI)
make mpi        # MPI parallel build
# Executable: lmp_serial or lmp_mpi
```

**Package management with Make:**
```bash
cd src
make yes-basic      # Enable common packages
make yes-openmp     # Enable OPENMP package
make pi             # View package status
make serial         # Build
```

## Testing & Validation

### Unit Tests (CTest)

```bash
# Configure with testing enabled
cmake -S cmake -B build \
      -C cmake/presets/gcc.cmake \
      -C cmake/presets/most.cmake \
      -D ENABLE_TESTING=on
cmake --build build
cd build && ctest -V

# Run specific test
ctest -V -R <test_name>
```

### Style/Coding Standard Checks

**ALWAYS run before submitting PRs:**
```bash
cd src
make check              # All checks
make check-whitespace   # Most common issue
make fix-whitespace     # Auto-fix whitespace
make fix-permissions    # Auto-fix file permissions
```

### Regression Tests

```bash
python3 -m venv testenv
source testenv/bin/activate
pip install numpy pyyaml junit_xml

python3 tools/regression-tests/run_tests.py \
    --lmp-bin=build/lmp \
    --config-file=tools/regression-tests/config_quick.yaml \
    --examples-top-level=examples
```

## Code Architecture

### Core Class Hierarchy

LAMMPS uses a composite pattern centered around the `LAMMPS` class (src/lammps.h), which acts as the central hub containing pointers to all major subsystems:

**Core Subsystem Classes:**
- **Memory** - memory allocation utilities
- **Error** - error handling and messaging
- **Universe** / **Comm** - MPI communication and parallel processing
- **Input** - script parsing and command execution
- **Atom** - atom data and topology (structure-of-arrays for performance)
- **Domain** - simulation box geometry and boundaries
- **Neighbor** - neighbor list construction
- **Force** - force field management (contains Pair, Bond, Angle, Dihedral, Improper, KSpace)
- **Modify** - fixes and computes registry with callback lists
- **Update** - time integration (contains Integrate and Minimize)
- **Group** - atom group management
- **Output** - thermodynamic output and dumps

**Pointers Base Class Pattern:**
All LAMMPS feature classes inherit from `Pointers` (src/pointers.h), which provides references to the main LAMMPS object's subsystems. This allows any class to access core functionality without passing numerous pointers.

### Style Registration System (Critical for Adding Features)

LAMMPS uses a macro-based factory pattern for extensibility. Each style (pair, fix, compute, etc.) uses a dual-purpose header:

```cpp
// In pair_lj_cut.h:
#ifdef PAIR_CLASS
PairStyle(lj/cut, PairLJCut);  // Registration macro
#else
class PairLJCut : public Pair { ... };  // Normal class definition
#endif
```

**How it works:**
1. CMake scans packages for style headers (e.g., `pair_*.h`, `fix_*.h`)
2. Generates aggregated `style_*.h` headers
3. Factory maps are created at compile time via preprocessor magic
4. Runtime: `pair_style lj/cut` looks up "lj/cut" in the factory map and instantiates `PairLJCut`

**To add a new style:**
1. Create `pair_mystyle.cpp` and `pair_mystyle.h` in a package directory
2. Follow the macro pattern (see existing styles as templates)
3. CMake automatically discovers and registers it
4. No manual registration needed

### Package Organization

Packages are subdirectories under `src/` (e.g., KSPACE, MOLECULE, KOKKOS). Each package:
- Contains style implementations
- Can be enabled/disabled at build time
- May have dependencies on other packages or external libraries

**Enable packages with CMake:**
```bash
cmake -S cmake -B build -D PKG_MOLECULE=on -D PKG_PYTHON=on
```

### Main Execution Flow

**Simulation loop** (from src/verlet.cpp):
```
for each timestep:
  modify->initial_integrate()      // Fix callbacks: update positions
  neighbor->build() if needed      // Rebuild neighbor lists
  force_clear()                    // Clear force arrays
  modify->pre_force()              // Fix callbacks before forces
  Force computation:
    pair->compute()                // Pairwise forces
    bond->compute()                // Bond forces
    angle->compute()               // Angle forces
    dihedral->compute()            // Dihedral forces
    improper->compute()            // Improper forces
    kspace->compute()              // Long-range forces
  modify->post_force()             // Fix callbacks after forces
  modify->final_integrate()        // Fix callbacks: update velocities
  modify->end_of_step()            // Fix callbacks: end of step
  output->write()                  // Thermodynamic output
```

**Fix Callback System:**
Fixes register themselves in callback lists via bitmask flags (see src/fix.h):
- `INITIAL_INTEGRATE` - before position update
- `POST_INTEGRATE` - after position update
- `PRE_FORCE` - before force computation
- `POST_FORCE` - after force computation
- `FINAL_INTEGRATE` - after velocity update
- `END_OF_STEP` - end of timestep

This enables modular extension of the integrator without modifying core code.

## Coding Standards & Conventions

### Critical Rules

1. **All source must be ASCII only** - Unicode characters are not allowed (security policy)
2. **No variable-length arrays (VLAs)** - Checked by CI, will cause build failures
3. **No alternative logical operators** - Use `&&` not `and`, `||` not `or`, `!` not `not` (MSVC compatibility)
4. **Whitespace matters** - Always run `make fix-whitespace` before committing
5. **File permissions:**
   - `.cpp` and `.h` files must NOT be executable
   - `.sh` and `.py` scripts SHOULD be executable

### Code Style

- Follow `.clang-format` rules in src/
- Use the `Pointers` base class to access core subsystems
- Atom data uses structure-of-arrays for cache efficiency
- Minimize virtual function calls in inner loops

### Documentation Requirements

**For new features:**
- Add documentation in `doc/src/` using reStructuredText (.rst)
- Use American English spelling
- Use ASCII characters only
- Include `.. versionadded:: TBD` for new commands
- Include `.. versionchanged:: TBD` for modified commands

**Build documentation:**
```bash
cd doc
make html          # Build HTML, check for warnings
make pdf           # Build PDF (requires pdflatex)
make spelling      # Check spelling
```

## GitHub Workflows / CI

**Pull Request Checks (all PRs to develop branch):**
1. **style-check.yml** - Runs coding standard checks (`make check-*` in src/)
2. **quick-regression.yml** - Builds with most packages, runs regression tests
3. **unittest-linux.yml** - Builds with LAMMPS_BIGBIG, runs unit tests via CTest

**Additional automated checks:**
- **codeql-analysis.yml** - Security scanning
- **check-vla.yml** - Checks for variable-length arrays
- **check-cpp23.yml** - C++23 compatibility check
- **compile-msvc.yml** - Windows MSVC compilation

**Workflow dependencies:**
- Ubuntu latest with ccache, ninja-build, libeigen3-dev, libcurl4-openssl-dev, python3-dev, mpi
- Python packages: numpy, pyyaml, junit_xml

## Development Workflow

1. **Branch:** Work on feature branches, submit PRs to `develop` (NOT master or release)
2. **Style check first:** Run `cd src && make check` before committing
3. **Build locally:** Test with `gcc.cmake + most.cmake` presets to match CI
4. **Test changes:** Run relevant unit tests and regression tests
5. **Watch CI:** All PR checks must pass

**Continuous release model:** The `develop` branch is always functional. All changes go through PRs with mandatory CI checks.

## Common Pitfalls

1. **Wrong CMake source directory:**
   - WRONG: `cmake -S . -B build`
   - CORRECT: `cmake -S cmake -B build`

2. **In-source builds:**
   - NEVER run cmake or make in repository root
   - ALWAYS create separate build/ directory

3. **Mixed build systems:**
   - Run cleanup commands when switching (see "Switching Build Systems" above)

4. **Unicode in source files:**
   - All source code must be ASCII only
   - Will cause CI to fail

5. **Missing style checks:**
   - Always run `cd src && make fix-whitespace` before committing
   - Check file permissions with `make check-permissions`

6. **Package dependencies:**
   - Some packages require others
   - CMake will warn; check console output

## Quick Reference

```bash
# Standard development cycle
mkdir build
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake -D ENABLE_TESTING=on -D DOWNLOAD_POTENTIALS=off
cmake --build build -j 4
cd src && make check  # Style checks
cd ../build && ctest -V  # Unit tests

# Minimal build for quick testing
cmake -S cmake -B build -C cmake/presets/basic.cmake
cmake --build build -j 4

# Add a package
cmake -S cmake -B build -D PKG_MOLECULE=on
cmake --build build -j 4

# Clean everything
rm -rf build
cd src && make clean-all

# Switch from Make to CMake
cd src && make purge
cd .. && mkdir build && cmake -S cmake -B build -C cmake/presets/basic.cmake
```

## NextTissUe Package

This repository includes a custom package called `NextTissUe` located in the `NextTissUe/` directory. This package has been integrated with the USER-AREABOND package (see src/USER-AREABOND/).

## Code Review Guidelines

When performing code reviews, apply:
- General contribution requirements: https://docs.lammps.org/Modify_requirements.html
- Programming style guidelines: https://docs.lammps.org/Modify_style.html

Ensure documentation changes:
- Use American English with plain ASCII characters
- Include `.. versionadded:: TBD` or `.. versionchanged:: TBD` directives
- Build with `make html`, `make pdf`, `make spelling` without NEW warnings/errors
- Check examples using new/modified commands for needed updates

## Additional Resources

- Official documentation: https://docs.lammps.org
- Developer info: https://docs.lammps.org/Developer.html
- Modify/extend LAMMPS: https://docs.lammps.org/Modify.html
- Error debugging: https://docs.lammps.org/Errors.html

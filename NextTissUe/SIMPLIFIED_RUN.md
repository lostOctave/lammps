# Simplified Simulation Setup for Quick Testing

## Changes Made

### 1. Reduced Relaxation/Equilibration Steps

| Phase | Original | New (Reduced) | Time Saved |
|-------|----------|---------------|------------|
| Compression | 10,000 steps | 5,000 steps | 50% faster |
| Post-compression equilibration | 100,000 steps | 10,000 steps | 90% faster |
| Pre-measurement equilibration | 100,000 steps | 5,000 steps | 95% faster |
| Main measurement run | 20,000,001 steps | 100,000 steps | 99.5% faster |
| **TOTAL** | **~20,210,000 steps** | **~120,000 steps** | **~600× faster!** |

### 2. Simplified Output

**Disabled** (commented out):
- ❌ `microrheology_probe_NC900_fd3.dat` (probe cell 1)
- ❌ `microrheology_probe2_NC900_fd3.dat` (probe cell 2)
- ❌ `microrheology_first_shell_NC900_fd3.dat` (7 nearest cells)
- ❌ `microrheology_second_shell_NC900_fd3.dat` (19 medium-range cells)
- ❌ `microrheology_third_shell_NC900_fd3.dat` (26 far cells)

**Enabled** (single simple output):
- ✅ `ovito_cells_NC900_fd3.lammpstrj` - All cells, OVITO-ready

### 3. Total Simulation Time

**Original setup**:
- Simulation time: ~20,210 time units
- Real time: ~30-60 minutes (depending on CPU)
- Output: ~9 GB of data

**New simplified setup**:
- Simulation time: ~120 time units
- Real time: **~1-3 minutes** 🚀
- Output: **~200 MB** (single file)

---

## Single Output File: `ovito_cells_NC900_fd3.lammpstrj`

### File Details

**Filename format**: `ovito_cells_NC{number_cells}_fd{force}.lammpstrj`
- Example: `ovito_cells_NC900_fd3.lammpstrj`

**Location**: `/mnt/d/NextTissUe_data/`

**Dump frequency**: Every 1000 steps = every 1.0 time units
- Over 100,000 steps → **100 frames total**
- Each frame has **900 cells × 50 beads = 45,000 atoms**

**File size**: ~200 MB (much smaller than before!)

**Format**: LAMMPS custom dump (compatible with OVITO)

### Data Columns

```
id mol type x y xu yu d_xc d_yc d_area d_peri d_shape
```

| Column | Variable | Meaning | Units | For OVITO |
|--------|----------|---------|-------|-----------|
| 1 | `id` | Atom (bead) ID | - | Unique identifier |
| 2 | `mol` | Molecule (cell) ID | - | **Color by molecule** |
| 3 | `type` | Atom type | - | All type 1 (single cell type) |
| 4 | `x` | Wrapped x-coordinate | LJ | Position (periodic image) |
| 5 | `y` | Wrapped y-coordinate | LJ | Position (periodic image) |
| 6 | `xu` | Unwrapped x-coordinate | LJ | **Actual position** |
| 7 | `yu` | Unwrapped y-coordinate | LJ | **Actual position** |
| 8 | `d_xc` | Cell center x | LJ | COM of cell |
| 9 | `d_yc` | Cell center y | LJ | COM of cell |
| 10 | `d_area` | Cell area | LJ² | **Color by area** |
| 11 | `d_peri` | Cell perimeter | LJ | Boundary length |
| 12 | `d_shape` | Shape index p = P/√A | - | **Color by shape** |

### Example Frame

```
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
45000
ITEM: BOX BOUNDS pp pp pp
-16.12 16.12
-13.96 13.96
-1.0 1.0
ITEM: ATOMS id mol type x y xu yu d_xc d_yc d_area d_peri d_shape
1 1 1 -15.234 -13.456 -15.234 -13.456 -15.100 -13.400 0.952 3.456 3.612
2 1 1 -15.201 -13.398 -15.201 -13.398 -15.100 -13.400 0.952 3.456 3.612
... (48 more beads from cell 1)
... (50 beads from cell 2)
... (50 beads from cell 3)
... (all 900 cells × 50 beads = 45000 rows)
```

---

## How to Use with OVITO

### Step 1: Load File in OVITO

1. Open OVITO
2. **File → Load File**
3. Navigate to `/mnt/d/NextTissUe_data/`
4. Select `ovito_cells_NC900_fd3.lammpstrj`
5. Click **Open**

OVITO will automatically:
- Detect LAMMPS custom format
- Parse column names
- Load all 100 frames

### Step 2: Basic Visualization

**Show atoms as particles**:
1. In pipeline, find **"Particles"** modifier
2. Change **Display** → **Particle display** → **Spheres** or **Points**
3. Adjust **Radius**: 0.02-0.05 (small, since there are 45,000 beads)

**Color by molecule** (to see individual cells):
1. Add modifier: **Color Coding**
2. Property: **Particle Identifier** → Select **"mol"**
3. Color gradient: **Rainbow** or **Jet**
4. Now each cell has a unique color!

### Step 3: Track Probe Cells

**Highlight probe cells (mol 709, 710)**:
1. Add modifier: **Expression Selection**
2. Expression: `Molecule == 709 || Molecule == 710`
3. This selects all beads belonging to probe cells

**Color probe cells differently**:
1. Add modifier: **Assign Color**
2. Keep selection active
3. Choose color: **Red** or **Yellow**
4. Now probe cells stand out!

### Step 4: Visualize Cell Properties

**Color by cell area** (shows compression):
1. Add modifier: **Color Coding**
2. Property: **d_area**
3. Color map: **Hot** (blue=compressed, red=expanded)
4. Adjust range: 0.7 - 1.0 LJ²

**Color by shape index** (shows rigidity):
1. Add modifier: **Color Coding**
2. Property: **d_shape**
3. Color map: **Viridis**
4. Adjust range: 3.5 - 3.8
5. p < 3.81 → fluid-like (blue)
6. p > 3.81 → solid-like (yellow/red)

### Step 5: Create Cell Meshes (Advanced)

**Show cells as polygons** (not just beads):
1. Add modifier: **Construct Surface Mesh**
2. Method: **Alpha shape**
3. Alpha: 0.1-0.2
4. Select **"Only selected"**: check
5. First select one molecule: `Molecule == 1`

This creates a polygon outline of one cell.

**For all cells** (computationally expensive):
- Remove selection
- Apply **Construct Surface Mesh** to all
- OVITO will create 900 cell boundaries
- Warning: May be slow!

### Step 6: Track Probe Cell Motion

**Measure probe cell displacement**:
1. Add modifier: **Compute Property**
2. Output property: **Displacement**
3. Expression for X: `xu - xu[TimestepStart]`
4. Expression for Y: `yu - yu[TimestepStart]`

**Plot probe trajectory**:
1. Select probe cell: `Molecule == 709`
2. Add modifier: **Time Series**
3. Source: **d_xc** (cell center x)
4. Click **"Show in graph"**
5. See COM trajectory over time!

### Step 7: Create Animation

**Play movie**:
- Use timeline slider at bottom
- Click **Play** button
- Watch cells move as probe is dragged!

**Export movie**:
1. **File → Export Animation**
2. Format: **MP4** or **GIF**
3. Frame range: 0-100
4. FPS: 10-30
5. Click **Export**

---

## Quick Visualization Recipes

### Recipe 1: See All Cells with Unique Colors
```
1. Load file
2. Add modifier: Color Coding → Property: "mol"
3. Particle radius: 0.03
4. Play animation
```

### Recipe 2: Highlight Probe Cells Being Dragged
```
1. Load file
2. Add modifier: Expression Selection → "Molecule == 709 || Molecule == 710"
3. Add modifier: Assign Color → Red
4. Play animation - probe cells are red, others are default color
```

### Recipe 3: Visualize Cell Deformation
```
1. Load file
2. Add modifier: Color Coding → Property: "d_area"
3. Range: 0.7 - 1.0
4. Color map: "Hot"
5. Blue cells = compressed, Red cells = expanded
```

### Recipe 4: Show Only Probe and Nearby Cells
```
1. Load file
2. Add modifier: Expression Selection
3. Expression: "Molecule == 709 || Molecule == 710 || Molecule == 680 || Molecule == 681 || Molecule == 711"
4. Add modifier: Delete Particles → Delete unselected
5. Now only see probe + nearby cells (clearer view)
```

---

## Running the Simulation

### Command
```bash
cd /home/lost_octave/LAMMPS/NextTissUe
../build/lmp -in in.demo_hexatic
```

### Expected Output (Console)
```
LAMMPS (XX XXX XXXX)
...
Step Temp PotEng KinEng TotEng Press probe1_x probe1_y ...
0    0.0  XXX    0.0    XXX    XXX   -15.1   -2.3    ...
100000 XXX XXX  XXX    XXX    XXX   -14.8   -2.1    ...  ← probe moved!
...
Loop time of XX.XX on 1 procs for 100000 steps
```

### Expected Runtime
- **1-3 minutes** on modern CPU
- Much faster than original 30-60 minutes!

### Output File
- **Location**: `/mnt/d/NextTissUe_data/ovito_cells_NC900_fd3.lammpstrj`
- **Size**: ~200 MB
- **Frames**: 100 (timesteps 0, 1000, 2000, ..., 100000)

---

## What You'll See in OVITO

### Initial Frame (timestep 0)
```
┌────────────────────────────────────────┐
│                                        │
│   ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯  │
│  ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯ │
│   ◯◉◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯  │  ← Probe cells (red)
│  ◯◯◉◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯ │     at initial position
│   ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯  │
│  ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯ │
│                                        │
└────────────────────────────────────────┘
```

### Final Frame (timestep 100000)
```
        F = 3.0 →
┌────────────────────────────────────────┐
│                                        │
│   ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯  │
│  ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯ │
│   ◯◯ ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◉  │  ← Probe cells moved
│  ◯◯ ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◉ │     to the right!
│   ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯  │  Nearby cells also shifted
│  ◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯◯ │
│                                        │
└────────────────────────────────────────┘
```

You'll see:
- Probe cells (red) move to the right
- Nearby cells (blue/green) are pulled along
- Distant cells (yellow/orange) barely move
- Some cells deform (area changes)

---

## Re-enabling Full Microrheology Output

If you want the detailed microrheology data later:

1. **Uncomment the MICRO dumps** in `in.demo_hexatic` (lines 332-341)
2. **Increase run steps** back to 20,000,001 (line 352)
3. **Run simulation** (will take 30-60 minutes)
4. **Analyze** with `python3 analyze_microrheology.py`

---

## Troubleshooting

### Issue: File not created
**Solution**: Check output directory exists
```bash
mkdir -p /mnt/d/NextTissUe_data
```

### Issue: OVITO can't read file
**Solution**: File might be compressed or corrupt
```bash
# Check file exists and has reasonable size
ls -lh /mnt/d/NextTissUe_data/ovito_cells_NC900_fd3.lammpstrj
# Should be ~200 MB

# Check first few lines
head -20 /mnt/d/NextTissUe_data/ovito_cells_NC900_fd3.lammpstrj
# Should show "ITEM: TIMESTEP" etc.
```

### Issue: Simulation too fast, want more detail
**Solution**: Increase dump frequency and run length
```lammps
dump OVITO_DUMP all custom 100 ...   # Every 100 steps instead of 1000
run 1000000                           # 10× longer run
```

---

## Summary

✅ **Reduced runtime**: ~600× faster (1-3 minutes instead of 30-60 minutes)
✅ **Simplified output**: Single file instead of 5+ complex files
✅ **OVITO-ready**: Direct visualization, no post-processing needed
✅ **All information preserved**: Still tracks probe cells, cell properties, all beads
✅ **Easy to extend**: Uncomment lines to re-enable full microrheology analysis

**Now you can quickly test parameters, visualize results, and iterate!** 🚀

# NextTissUe Model Physics Analysis

## 1. Cell Size and Box Size Relationship

### Current Configuration
- **Box size**: 60.0 × 51.6 LJ units (after 2× scaling)
- **Number of cells**: 900
- **Cell radius**: R = 0.5 LJ units
- **Cell diameter**: d = 1.0 LJ units
- **Lattice spacing**: a ≈ 1.0 LJ units (hexatic lattice)

### Packing Analysis
```
Cell area (geometric): A_geom = πR² = π(0.5)² = 0.785 LJ²
Total cell area: 900 × 0.785 = 706.9 LJ²
Box area: 60.0 × 51.6 = 3096 LJ²
Geometric packing fraction: φ_geom = 706.9 / 3096 = 0.228 (22.8%)
```

**However**, the simulation compresses the system to a target density:
```lammps
variable num_dense equal 1
variable L equal (${number_cells}/${num_dense})^0.5
# L = sqrt(900/1) = 30 LJ units
```

After compression, the box becomes **30 × 30 = 900 LJ²**, giving:
```
Final packing fraction: φ_final = 706.9 / 900 = 0.785 (78.5%)
```

### Relationship
- **Initial state**: Cells just touching (φ ≈ 0.23)
- **After compression**: High-density packing (φ ≈ 0.79)
- The compression creates cell-cell contacts and deformations

---

## 2. Preferred Cell Size (Target Area)

### Shape Parameter p₀
The model uses a **dimensionless shape parameter** p₀ = 3.65:

```lammps
variable p0 equal 3.65
variable area_prefac equal perimeter² / (area × p0²)
variable modified_area_large equal area_prefac × area_large
```

### Preferred Area Calculation
```
Perimeter: P = 2πR = 2π(0.5) = 3.14159 LJ
Geometric area: A_geom = πR² = 0.78540 LJ²

Area prefactor: area_prefac = P² / (A_geom × p0²)
                            = (3.14159)² / (0.78540 × 3.65²)
                            = 9.8696 / 10.4627
                            = 0.9433

Preferred area: A₀ = area_prefac × A_geom
                   = 0.9433 × 0.78540
                   = 0.7409 LJ²
```

### Physical Meaning
The **shape parameter** p₀ is defined as:
```
p₀ = P / √A
```

For a **perfect circle**: p₀_circle = 2√π ≈ 3.545

Your simulation has **p₀ = 3.65 > 3.545**, meaning:
- Cells prefer to be **slightly more elongated** than perfect circles
- OR cells are **slightly compressed** (area reduced relative to perimeter)
- This creates a **rigidity** that resists deformation

**Answer**: Yes, cells have a preferred area A₀ ≈ 0.741 LJ², enforced by the area constraint potential.

---

## 3. Free Energy of the Cell

The total free energy per cell has several contributions:

### 3.1 Bond Stretching Energy
Each cell has 50 harmonic bonds:
```
E_bond = Σ(i=1 to 50) [K_bond × (ℓ_i - ℓ₀)²]

Where:
- K_bond = 1,000,000 / 50 = 20,000 (very stiff!)
- ℓ₀ = perimeter / N_beads = 2πR / 50 = 0.0628 LJ
```

**At equilibrium** (circular cell): E_bond ≈ 0

### 3.2 Area Constraint Energy
```cpp
// From bond_harmonic_area.cpp line 298:
E_area = 0.5 × k_area × (A - A₀)² / N_beads

Where:
- k_area = 1000 (area stiffness)
- A = current cell area
- A₀ = 0.741 LJ² (preferred area)
- N_beads = 50
```

**Per cell**: E_area,cell = 50 × [0.5 × 1000 × (A - 0.741)²] / 50
                           = 500 × (A - A₀)²

**Physical interpretation**: 
- If A = A₀: E_area = 0 (minimum energy)
- If A = 0.8 LJ²: E_area = 500 × (0.8 - 0.741)² = 1.74 ε
- Compression/expansion costs energy quadratically

### 3.3 Pair Interaction Energy
Between beads of different cells:

**Repulsive (WCA)**: Weeks-Chandler-Andersen potential
```
U_WCA(r) = 4ε[(σ/r)¹² - (σ/r)⁶] + ε   for r < r_cut
         = 0                           for r ≥ r_cut

Where:
- ε = 1/50 = 0.02 (per bead interaction)
- σ = ℓ_bond × 2^(-1/6) ≈ 0.0559 LJ
- r_cut = ℓ_bond = 0.0628 LJ (purely repulsive)
```

**Attractive (after compression)**: Modified WCA with attraction
```
σ = 1.1 × ℓ_bond = 0.069 LJ
r_cut = 2.0 × ℓ_bond = 0.126 LJ (includes attraction)
ε_att = 10 × ε_WCA = 0.2 (attraction strength)
```

### Total Free Energy per Cell
```
F_cell = E_bond + E_area + E_pair + TS (entropy term)

At equilibrium (circular, isolated cell):
F_cell ≈ 0 + 0 + 0 = 0

Under compression (deformed, contacting neighbors):
F_cell ≈ E_bond(deformation) + E_area(compression) + E_pair(contacts) - TS
```

**Typical values from your simulation output**:
```
Step 0: PotEng = 0.036 per atom → Total = 0.036 × 45000 = 1620 ε
      → Per cell ≈ 1620 / 900 = 1.8 ε
```

---

## 4. Box Enlargement and Vacancy Formation

### Current State Analysis

**Initial (before compression)**:
- Box: 60 × 51.6 = 3096 LJ²
- Cells just touching, φ ≈ 0.23
- **NO vacancies** (cells fill space completely in hexatic lattice)

**After compression**:
- Box: 30 × 30 = 900 LJ²
- High density, φ ≈ 0.79
- Cells compressed and deformed
- Strong cell-cell contacts

### If Box Keeps Enlarging

**Scenario**: Slowly increase box size from 30×30 back to 60×51.6

#### Stage 1: Elastic Expansion (30 → ~35 LJ box)
- Cells relax from compressed state
- Area returns to A₀ ≈ 0.741 LJ²
- Cells become more circular
- **No vacancies yet** (still touching)

#### Stage 2: Contact Breaking (~35 → ~42 LJ box)
- φ drops below ~0.6
- Cell-cell contacts begin to break
- **Small gaps appear** between cells
- System becomes **under-coordinated**

#### Stage 3: Vacancy Formation (42 → 60 LJ box)
- φ drops to ~0.23 (original state)
- **Large vacancies** form between cells
- Hexatic lattice structure may:
  - **Remain ordered** (cells stay in lattice positions with gaps)
  - **Collapse** (cells aggregate to maintain contacts, leaving voids)

### What Determines Behavior?

1. **Active forces** (v₀ = 0.01):
   - Self-propulsion drives cells to explore
   - Can lead to **phase separation**: dense clusters + voids

2. **Attractive interactions** (ε_att = 0.2):
   - Cells prefer to stay in contact
   - Resist vacancy formation
   - May form **compact clusters** instead of uniform lattice

3. **Rotational diffusion** (Dr = 1.0):
   - High rotational noise
   - Cells reorient frequently
   - Promotes **mixing** and uniform distribution

### Expected Outcome

Given your parameters:
- **v₀ = 0.01** (weak activity)
- **Dr = 1.0** (strong rotational noise)
- **ε_att = 0.2** (moderate attraction)

**Prediction**: 
✓ **Vacancies WILL appear** as box enlarges
✓ Cells will likely form **small clusters** (3-10 cells)
✓ Clusters will slowly diffuse and rearrange
✓ Eventually reach **uniform low-density state** with gaps

**Critical box size for vacancy onset**:
```
L_critical ≈ sqrt(900 × A₀ / φ_contact)
where φ_contact ≈ 0.6 (percolation threshold)

L_critical ≈ sqrt(900 × 0.741 / 0.6) = sqrt(1111.5) ≈ 33.3 LJ
```

**Answer**: Yes, vacancies appear when L > ~33 LJ units.

---

## Summary Table

| Property | Value | Units | Notes |
|----------|-------|-------|-------|
| Cell radius | 0.5 | LJ | Geometric |
| Cell perimeter | 3.14 | LJ | 2πR |
| Geometric area | 0.785 | LJ² | πR² |
| Preferred area A₀ | 0.741 | LJ² | Modified by p₀ |
| Shape parameter p₀ | 3.65 | - | > 3.545 (circle) |
| Bond stiffness | 20,000 | ε/LJ² | Very rigid |
| Area stiffness | 1000 | ε/LJ⁴ | Moderate |
| Initial box | 60 × 51.6 | LJ² | Low density |
| Final box (compressed) | 30 × 30 | LJ² | High density |
| Initial packing φ | 0.23 | - | Just touching |
| Final packing φ | 0.79 | - | Highly compressed |
| Free energy (isolated) | ~0 | ε | Equilibrium |
| Free energy (compressed) | ~1.8 | ε | Per cell |
| Vacancy onset | L > 33 | LJ | Critical box size |

---

## References

This model is based on the **Vertex Model** for epithelial tissues:
- Bi, D., et al. "A density-independent rigidity transition in biological tissues." Nature Physics (2015)
- The shape parameter p₀ controls the solid-fluid transition
- p₀ < 3.81: fluid-like behavior
- p₀ > 3.81: solid-like behavior (jammed)

Your simulation uses **p₀ = 3.65**, placing it in the **fluid phase** near the transition.


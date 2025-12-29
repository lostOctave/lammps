/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(harmonic/area,AreaBondHarmonic)

#else

#ifndef LMP_BOND_AREAHARMONIC_H
#define LMP_BOND_AREAHARMONIC_H

#include "bond.h"
#include "fix.h"
namespace LAMMPS_NS {
class AreaBondHarmonic : public Bond {
 public:
  AreaBondHarmonic(class LAMMPS *);
  virtual ~AreaBondHarmonic();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  virtual void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);
  virtual void *extract(const char *, int &);

 protected:
  double *k,*r0;

  int nchunk,maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  double *area,*area_all,*peri_local,*peri_all;
  double findCircle(double x1, double y1, double x2, double y2, double x3, double y3);
  double *local_max_x, *local_min_x, *max_x, *min_x, *local_max_y, *local_min_y, *max_y, *min_y, *local_xc, *local_yc, *xc, *yc, **bond_left, **bond_right; // **bond_left_global, **bond_right_global;
  int *local_size_cell, *size_cell;
  virtual void ev_area_tally(int, double, double, double, double, double, double);
  virtual void allocate();
  virtual void allocate_bonds();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/

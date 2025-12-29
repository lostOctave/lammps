/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(propel/cell,FixPropelCell);
// clang-format on
#else

#ifndef LMP_FIX_PROPEL_CELL_H
#define LMP_FIX_PROPEL_CELL_H

#include "fix.h"
namespace LAMMPS_NS {

class FixPropelCell : public Fix {
 public:
  FixPropelCell(class LAMMPS *, int, char **);
  virtual ~FixPropelCell();
  void init();
  void post_force(int);
  void setup(int);
  int setmask();
  //class Update *update;
  //class RanMars *get_random;
 private:
  double magnitude, visc;
  double Dr, dt, damp, eps_0, l_0;
  double sx, sy, sz;
  int mode;

  void post_force_dipole(int);
  void post_force_velocity(int);
  void post_force_quaternion(int);
  void post_force_cell(int);
  void post_force_temp(int);
  class AtomVecEllipsoid *avec;
protected:
  int nchunk,maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  double *local_max_x, *local_min_x, *max_x, *min_x, *local_xc, *local_yc, *xc, *yc;
  double *local_prop_vec_x, *local_prop_vec_y, *prop_vec_x, *prop_vec_y;
  //double *local_theta_change_cil, *theta_change_cil;
  int *local_size_cell, *size_cell, *local_num_nocontacts, *num_nocontacts;
  class RanMars *random;
  double *local_rand_x, *local_rand_y, *chunk_rand_x, *chunk_rand_y;
  int seed;
  //virtual void allocate();
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:
E: Illegal fix propel/self command.
Wrong number/type of input arguments.
E: Fix propel/self requires atom attribute mu with option dipole.
Self-explanatory.
E: Fix propel/self requires atom style ellipsoid with option quat.
Self-explanatory.
Fix propel/self requires extended particles with option quat.
Self-explanatory.
*/

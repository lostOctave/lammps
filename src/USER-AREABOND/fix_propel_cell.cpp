/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   Contributed by Stefan Paquay @ Brandeis University
   Thanks to Liesbeth Janssen @ Eindhoven University for useful discussions!
   Current maintainer: Sam Cameron @ University of Bristol
----------------------------------------------------------------------- */
/*----------------------------------------------------------------------
Developed by Anshuman Pasupalak and Massimo Pica Ciamarra @ Nanyang Technological University

Extension of USER-MISC package fix_self_propel.cpp

--------------------------------------------------------------------------*/
#include "fix_propel_cell.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "update.h"

#include <cstring>
#include <iostream>
#include <cmath>
#include "neighbor.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "random_mars.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;
enum { DIPOLE, VELOCITY, QUAT, CELL, TEMP };

#define TOL 1e-14

/* ---------------------------------------------------------------------- */

FixPropelCell::~FixPropelCell(){
  if (mode == TEMP || mode == CELL) {
    memory->destroy(local_rand_x);
    memory->destroy(local_rand_y);
    memory->destroy(chunk_rand_x);
    memory->destroy(chunk_rand_y);
    memory->destroy(size_cell);
    memory->destroy(local_size_cell);
    memory->destroy(local_prop_vec_x);
    memory->destroy(local_prop_vec_y);
    memory->destroy(prop_vec_x);
    memory->destroy(prop_vec_y);
    memory->destroy(local_num_nocontacts);
    memory->destroy(num_nocontacts);
  }
}


FixPropelCell::FixPropelCell(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), avec(nullptr)
{
  if (narg != 5 && narg !=7 && narg != 9) error->all(FLERR, "Illegal fix propel/self command -improper n_args");

  if (strcmp(arg[3], "velocity") == 0) {
    mode = VELOCITY;
    thermo_virial = 0;
  } else if (strcmp(arg[3], "dipole") == 0) {
    mode = DIPOLE;
    thermo_virial = 1;
  } else if (strcmp(arg[3], "quat") == 0) {
    mode = QUAT;
    thermo_virial = 1;
   } else if (strcmp(arg[3], "cell") == 0) {
    mode = CELL;
    thermo_virial = 1;
  } else if (strcmp(arg[3], "temp") == 0) {
    mode = TEMP;
    thermo_virial = 1;
  }
   else {
    error->all(FLERR, "Illegal fix propel/self command-mode invalid");
  }

  magnitude = utils::numeric(FLERR, arg[4], false, lmp);
if (narg ==7 )  {
    if (mode != TEMP ) { error->all(FLERR, "Illegal fix propel/self command-mode mismatch"); }
     seed = utils::numeric(FLERR, arg[6], false, lmp);
     //if(mode == CELL) Dr = utils::numeric(FLERR, arg[5], false, lmp);
     if(mode == TEMP) damp = utils::numeric(FLERR, arg[5], false, lmp);

  }
  if (narg == 9) {
    if (mode != QUAT && mode != CELL) { error->all(FLERR, "Illegal fix propel/self command"); }
    if (strcmp(arg[5], "qvector") == 0) {
      sx = utils::numeric(FLERR, arg[6], false, lmp);
      sy = utils::numeric(FLERR, arg[7], false, lmp);
      sz = utils::numeric(FLERR, arg[8], false, lmp);
      double snorm = sqrt(sx * sx + sy * sy + sz * sz);
      sx = sx / snorm;
      sy = sy / snorm;
      sz = sz / snorm;
    } else if(mode == CELL){
        visc = utils::numeric(FLERR, arg[5], false, lmp); // viscosity used in fix viscous
        damp = utils::numeric(FLERR, arg[6], false, lmp);  // timescale for contact dynamics (propel/cell)
        Dr = utils::numeric(FLERR, arg[7], false, lmp); // Dr for random noise in fix propel/cell
        seed = utils::numeric(FLERR, arg[8], false, lmp); // random seed
    }
    else {
      error->all(FLERR, "Illegal fix propel/self command");
    }
  } else {
    sx = 1.0;
    sy = 0.0;
    sz = 0.0;
  }

  if(mode == TEMP || mode == CELL){ 
    random = new RanMars(lmp,comm->me+seed); // we should pass a seed at command line in both model CELL and TEMP

    char arg[] = "Cells"; // name of the chunk command
    int icompute = modify->find_compute(arg);
    if (icompute >= 0) {
      cchunk = (ComputeChunkAtom *) modify->compute[icompute];
      // it finds the compute chunk properly
    }else{
      error->all(FLERR, "Illegal fix propel/self command - cannot find chunk id");
    }
    nchunk = cchunk->setup_chunks();
    maxchunk = nchunk;
    memory->create(local_rand_x, nchunk, "propel/cell:local_rand_x");
    memory->create(local_rand_y, nchunk, "propel/cell:local_rand_y");
    memory->create(chunk_rand_x, nchunk, "propel/cell:chunk_rand_x");
    memory->create(chunk_rand_y, nchunk, "propel/cell:chunk_rand_y");
    memory->create(size_cell, nchunk, "propel/cell:size_cell");
    memory->create(local_size_cell, nchunk, "propel/cell:local_size_cell");
    memory->create(local_prop_vec_x, nchunk, "propel/cell:local_prop_vec_x");
    memory->create(prop_vec_x, nchunk, "propel/cell:prop_vec_x");
    memory->create(local_prop_vec_y, nchunk, "propel/cell:local_prop_vec_y");
    memory->create(prop_vec_y, nchunk, "propel/cell:prop_vec_y");
    memory->create(num_nocontacts, nchunk, "propel/cell:num_nocontacts");
    memory->create(local_num_nocontacts, nchunk, "propel/cell:local_num_nocontacts");
    //memory->create(local_theta_change_cil, nchunk, "propel/cell:local_theta_change_cil");
    //memory->create(theta_change_cil, nchunk, "propel/cell:theta_change_cil");
  }

  // check for keyword

}

/* ---------------------------------------------------------------------- */

int FixPropelCell::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
//  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropelCell::init()
{
  if (mode == DIPOLE && !atom->mu_flag)
    error->all(FLERR, "Fix propel/self requires atom attribute mu with option dipole");

  if (mode == QUAT) {
    avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
    if (!avec) error->all(FLERR, "Fix propel/self requires atom style ellipsoid with option quat");

    // check that all particles are finite-size ellipsoids
    // no point particles allowed, spherical is OK

    int *ellipsoid = atom->ellipsoid;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (ellipsoid[i] < 0)
          error->one(FLERR, "Fix propel/self requires extended particles with option quat");
  }
}

/* ---------------------------------------------------------------------- */

void FixPropelCell::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPropelCell::post_force(int vflag)
{
  if (mode == DIPOLE)
    post_force_dipole(vflag);
  else if (mode == VELOCITY)
    post_force_velocity(vflag);
  else if (mode == QUAT)
    post_force_quaternion(vflag);
  else if (mode == CELL)
    post_force_cell(vflag);
  else if (mode == TEMP)
    post_force_temp(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPropelCell::post_force_dipole(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **mu = atom->mu;
  double fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      fx = magnitude * mu[i][0];
      fy = magnitude * mu[i][1];
      fz = magnitude * mu[i][2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
}
/*-----------------------------------------------*/

void FixPropelCell::post_force_cell(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double fx, fy, fz;
  bigint time_step = update->ntimestep;
  dt = update->dt;
  tagint *molecule = atom->molecule;/// Find molecule id
  int flag;
  int cols;
  // load custom properties
  char fx_active[] = "fx_active";
  int index_fx_active = atom->find_custom(fx_active,flag,cols);
  double *fxactive = atom->dvector[index_fx_active];

  char fy_active[] = "fy_active";
  int index_fy_active = atom->find_custom(fy_active,flag,cols);
  double *fyactive = atom->dvector[index_fy_active];

  char angle_director[] = "theta";
  int index_angle_director = atom->find_custom(angle_director,flag,cols);
  double *angledirector = atom->dvector[index_angle_director];

  char f_active[] = "factive";
  int index_f_active = atom->find_custom(f_active,flag,cols);
  double *factive = atom->dvector[index_f_active];

  char mag_factive[] = "mag_factive";
  int index_mag_factive = atom->find_custom(mag_factive,flag,cols);
  double *magfactive = atom->dvector[index_mag_factive];

  char XC[] = "xc";
  int index_xc = atom->find_custom(XC,flag,cols);
  double *XCcell = atom->dvector[index_xc];

  char YC[] = "yc";
  int index_yc = atom->find_custom(YC,flag,cols);
  double *YCcell = atom->dvector[index_yc];


  char theta_0[] = "theta_0";
  int index_theta_0 = atom->find_custom(theta_0,flag,cols);
  double *theta_in = atom->dvector[index_theta_0];

  char d_theta2[] = "theta2";
  int index_d_theta2 = atom->find_custom(d_theta2,flag,cols);
  double *del_theta2 = atom->dvector[index_d_theta2];

  char pair_forcex[] = "force_pair_x";
  int index_pair_forcex = atom->find_custom(pair_forcex,flag,cols);
  double *pair_f_x = atom->dvector[index_pair_forcex];

  char pair_forcey[] = "force_pair_y";
  int index_pair_forcey = atom->find_custom(pair_forcey,flag,cols);
  double *pair_f_y = atom->dvector[index_pair_forcey];

  char contact_prop_dir[] = "prop_dir";
  int index_contact_prop_dir = atom->find_custom(contact_prop_dir,flag,cols);
  double *contact_prop_theta = atom->dvector[index_contact_prop_dir];
  
  //char force_wall[] = "force_wall_total";
  //int	index_force_wall = atom->find_custom(force_wall,flag,cols);
  //double *force_wall_total = atom->dvector[index_force_wall];
  

  //
  double xlo = domain->boxlo[0];
  double xhi = domain->boxhi[0];
  double lx = xhi-xlo;
  double ylo = domain->boxlo[1];
  double yhi = domain->boxhi[1];
  double ly = yhi-ylo;
  double unwrap[3];
  double half_lx = lx*0.5;
  double half_ly = ly*0.5;
  double x1,x2,y1,y2, dist, Cos, Sin;
  double net_pair_force;
  //tagint *molecule = atom->molecule;/// Find molecule id
  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  imageint *image = atom->image;

  // need to allocate more memory if the number of chunks has changed
  nchunk = cchunk->setup_chunks();
  if (nchunk > maxchunk){
    maxchunk = nchunk;
    memory->create(local_rand_x,nchunk,"propel/cell:local_rand_x");
    memory->create(chunk_rand_x,nchunk,"propel/cell:chunk_rand_x");
    memory->create(size_cell, nchunk,  "propel/cell:size_cell");
    memory->create(local_size_cell, nchunk, "propel/cell:local_size_cell");
    //memory->create(local_theta_change_cil, nchunk, "propel/cell:local_theta_change_cil");
    //memory->create(theta_change_cil, nchunk, "propel/cell:theta_change_cil");
    memory->create(local_prop_vec_x, nchunk, "propel/cell:local_prop_vec_x");
    memory->create(prop_vec_x, nchunk, "propel/cell:prop_vec_x");
    memory->create(local_prop_vec_y, nchunk, "propel/cell:local_prop_vec_y");
    memory->create(prop_vec_y, nchunk, "propel/cell:prop_vec_y");
    memory->create(num_nocontacts, nchunk, "propel/cell:num_nocontacts");
    memory->create(local_num_nocontacts, nchunk, "propel/cell:local_num_nocontacts");
  }

  for(int i = 0; i < nchunk; i++){
    local_rand_x[i] = chunk_rand_x[i] = 0;
    size_cell[i] = local_size_cell[i] = 0;
    //local_theta_change_cil[i] = theta_change_cil[i] = 0;
    num_nocontacts[i] = local_num_nocontacts[i] = 0;
    local_prop_vec_x[i] = local_prop_vec_y[i] = 0;
    prop_vec_x[i] = prop_vec_y[i] = 0;
  }

  // random number generation for the different chunk
  int index1;
  double vec_theta_x, vec_theta_y, vec_rad_x, vec_rad_y;  //vec_theta = current propulsion director in vector form, //vec_rad stores radial vectors of bead
  double dot_prod, cos_theta_by2_sq;
  for(int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      net_pair_force = sqrt(powl(pair_f_x[i],2) + powl(pair_f_y[i],2)) ;//+ force_wall_total[i]; // net force due to pair interactions, if zero then propulsion along this direction
      index1 = molecule[i]-1;
      if(index1 >= 0){
        local_rand_x[index1] += random->gaussian();
        local_size_cell[index1]++;
        if(net_pair_force == 0){
          domain->unmap(x[i], image[i], unwrap);
          x1 = unwrap[0]-XCcell[i];
          y1 = unwrap[1]-YCcell[i];
          
          if(x1 > half_lx) x1 -= lx; if(x1 < -half_lx) x1 += lx;
          if(y1 > half_ly) y1 -= ly; if(y1 < -half_ly) y1 += ly;
          dist = sqrt(x1*x1+y1*y1);
          vec_theta_x = cos(angledirector[i]);
          vec_theta_y = sin(angledirector[i]);
          
          vec_rad_x = x1/dist;
          vec_rad_y = y1/dist;
          dot_prod = (vec_theta_x*vec_rad_x + vec_theta_y*vec_rad_y); // p.r
          cos_theta_by2_sq = (1 + dot_prod)/2;
          local_prop_vec_x[index1] += cos_theta_by2_sq*vec_rad_x;
          local_prop_vec_y[index1] += cos_theta_by2_sq*vec_rad_y;
          local_num_nocontacts[index1] +=1;
          //del_theta2[i] = dot_prod;
          

        }
       }
    }
  }
  MPI_Allreduce(local_rand_x,chunk_rand_x,nchunk,MPI_DOUBLE,MPI_SUM,world); // sum of the random numbers for each chunk
  MPI_Allreduce(local_size_cell,size_cell,nchunk,MPI_INT,MPI_SUM,world);    // chunk size
  //MPI_Allreduce(local_theta_change_cil,theta_change_cil,nchunk,MPI_DOUBLE,MPI_SUM,world); //theta change as per cil
  MPI_Allreduce(local_num_nocontacts,num_nocontacts,nchunk,MPI_INT,MPI_SUM,world);    // number of particles not in contact with any cell

  MPI_Allreduce(local_prop_vec_x,prop_vec_x,nchunk,MPI_DOUBLE,MPI_SUM,world); // sum of the non contact atom vectors for each chunk
  MPI_Allreduce(local_prop_vec_y,prop_vec_y,nchunk,MPI_DOUBLE,MPI_SUM,world); // sum of the non contact atom vectors for each chunk
  for(int i = 0; i < nchunk; i++){
    chunk_rand_x[i] /= sqrt(size_cell[i]);
  }
  ///////////////////////////////////////////////// we need a prefactor here as we are using uniform not gaussian number; to check

  bigint seed;
  double theta_change_diffusion, theta_no_contact;
  //cout<<Dr<<"	"<<dt<<endl;
  for(int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
    //cout<<"before index1"<<endl;
    index1 = molecule[i]-1;
    if(index1 >= 0){
      // each atom knows the center of its cell
      domain->unmap(x[i], image[i], unwrap);
      x1 = unwrap[0]-XCcell[i];
      y1 = unwrap[1]-YCcell[i];
      if(x1 > half_lx) x1 -= lx; if(x1 < -half_lx) x1 += lx;
      if(y1 > half_ly) y1 -= ly; if(y1 < -half_ly) y1 += ly;
      dist = sqrt(x1*x1+y1*y1);
      x1 /= dist;
      y1 /= dist;
      net_pair_force = sqrt(powl(pair_f_x[i],2) + powl(pair_f_y[i],2)) ;

      theta_no_contact = atan2(prop_vec_y[index1],prop_vec_x[index1]);
      
      contact_prop_theta[i] = theta_no_contact;
      theta_change_diffusion = sqrt(2*Dr*dt)*chunk_rand_x[index1];

      double normalise_num_contacts; // if the cell has no contacts or is in complete contact, propulsion direction should not depend on theta_no_contact
      normalise_num_contacts = num_nocontacts[index1]*(size_cell[index1]-num_nocontacts[index1])/powl(size_cell[index1],2);
      angledirector[i] += sin(contact_prop_theta[i]-angledirector[i]) *normalise_num_contacts* dt/damp + theta_change_diffusion; //*powl((magnitude/factive[i]),5)
      //del_theta2[i] =  normalise_num_contacts;  // powl((angledirector[i]-theta_in[i]),2);
      Cos = cos(angledirector[i]);
      Sin = sin(angledirector[i]);
      // scalar product & active force
      dist = sqrt(x1*Cos + y1*Sin);
      if(dist > 0){
        fxactive[i] = dist*x1*magnitude*visc;//factive[i]
        fyactive[i] = dist*y1*magnitude*visc;//*factive[i]
        //magfactive[i] = sqrt(fxactive[i]*fxactive[i] + fyactive[i]*fyactive[i]);;
      }else{
        fxactive[i] = 0;
        fyactive[i] = 0;
        //magfactive[i] = 0;
      }
      // add active force
      f[i][0] += fxactive[i];
      f[i][1] += fyactive[i];

        if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fxactive[i] * unwrap[0];
        vi[1] = fyactive[i] * unwrap[1];
        vi[2] = 0 * unwrap[2];
        vi[3] = fxactive[i] * unwrap[1];
        vi[4] = fxactive[i] * unwrap[2];
        vi[5] = fyactive[i] * unwrap[2];
        //v_tally(i, vi);
      }
    }
 }
 }
}
/*-----------------------------------------------*/


void FixPropelCell::post_force_temp(int vflag)
{
//  printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAa\n"); //,-mass[type[i]]/pre_factor_drag, sqrt(pre_factor_rand));

  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *mass = atom->mass;
  int *type = atom->type;
  double fx, fy, fz;
  bigint time_step = update->ntimestep;
  double dt = update->dt;
  tagint *molecule = atom->molecule;/// Find molecule id
  int flag;
  int cols;
  // load custom properties
  char fx_active[] = "fx_active";
  int index_fx_active = atom->find_custom(fx_active,flag,cols);
  double *fxactive = atom->dvector[index_fx_active];

  char fy_active[] = "fy_active";
  int index_fy_active = atom->find_custom(fy_active,flag,cols);
  double *fyactive = atom->dvector[index_fy_active];


  char f_active[] = "factive";
  int index_f_active = atom->find_custom(f_active,flag,cols);
  double *factive = atom->dvector[index_f_active];

  char mag_factive[] = "mag_factive";
  int index_mag_factive = atom->find_custom(mag_factive,flag,cols);
  double *magfactive = atom->dvector[index_mag_factive];
  //
  double xlo = domain->boxlo[0];
  double xhi = domain->boxhi[0];
  double lx = xhi-xlo;
  double ylo = domain->boxlo[1];
  double yhi = domain->boxhi[1];
  double ly = yhi-ylo;
  double unwrap[3];
  double half_lx = lx*0.5;
  double half_ly = ly*0.5;
  double x1,x2,y1,y2, dist, Cos, Sin;
  //tagint *molecule = atom->molecule;/// Find molecule id
  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  imageint *image = atom->image;

  // need to allocate more memory if the number of chunks has changed
  nchunk = cchunk->setup_chunks();
  if (nchunk > maxchunk){
    maxchunk = nchunk;
    memory->create(local_rand_x,nchunk,"propel/cell:local_rand_x");
    memory->create(local_rand_y,nchunk,"propel/cell:local_rand_y");
    memory->create(chunk_rand_x,nchunk,"propel/cell:chunk_rand_x");
    memory->create(chunk_rand_y,nchunk,"propel/cell:chunk_rand_y");
    memory->create(size_cell, nchunk,  "propel/cell:size_cell");
    memory->create(local_size_cell, nchunk, "propel/cell:local_size_cell");
  }

  for(int i = 0; i < nchunk; i++){
    local_rand_x[i] = local_rand_y[i] = chunk_rand_x[i] = chunk_rand_y[i] = 0;
    size_cell[i] = local_size_cell[i] = 0;
  }

  // random number generation for the different chunk
  int index1;
  for(int i = 0; i < nlocal; i++){
    index1 = molecule[i]-1;
    if(index1 >= 0){
      local_rand_x[index1] += random->uniform()-0.5;
      local_rand_y[index1] += random->uniform()-0.5;
      local_size_cell[index1]++;
    }
  }

  // if the random number have variance 1, the sum over N of them has variance sqrt(N); we need to account for this
  MPI_Allreduce(local_rand_x,chunk_rand_x,nchunk,MPI_DOUBLE,MPI_SUM,world); // sum of the random numbers for each chunk
  MPI_Allreduce(local_rand_y,chunk_rand_y,nchunk,MPI_DOUBLE,MPI_SUM,world); // sum of the random numbers for each chunk
  MPI_Allreduce(local_size_cell,size_cell,nchunk,MPI_INT,MPI_SUM,world);    // chunk size
  for(int i = 0; i < nchunk; i++){
    chunk_rand_x[i] /= sqrt(size_cell[i]);
    chunk_rand_y[i] /= sqrt(size_cell[i]);
  }


  bigint seed;
  double pre_factor_rand = 24 * (force->boltz * magnitude) / (damp * force->mvv2e * dt);

  double pre_factor_drag = damp * force->ftm2v;

  for(int i = 0; i < nlocal; i++){
     if (mask[i] & groupbit) {

    index1 = molecule[i]-1;

    if(index1 >= 0){

        fxactive[i] = sqrt(pre_factor_rand)*(chunk_rand_x[index1]);
        fyactive[i] = sqrt(pre_factor_rand)*(chunk_rand_y[index1]);

      // add active force
      f[i][0] += fxactive[i] - mass[type[i]] * v[i][0] /pre_factor_drag;
      f[i][1] += fyactive[i] - mass[type[i]] * v[i][1] /pre_factor_drag;

	//}
        if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fxactive[i] * unwrap[0];
        vi[1] = fyactive[i] * unwrap[1];
        vi[2] = 0 * unwrap[2];
        vi[3] = fxactive[i] * unwrap[1];
        vi[4] = fxactive[i] * unwrap[2];
        vi[5] = fyactive[i] * unwrap[2];
        //v_tally(i, vi);
      }
    }
  }
}
}
/* ----------------------------------------------------------------*/

/* ---------------------------------------------------------------------- */

void FixPropelCell::post_force_velocity(int vflag)
{
  double **f = atom->f;
  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double nv2, fnorm, fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {

      nv2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
      fnorm = 0.0;

      if (nv2 > TOL) {

        // Without this check you can run into numerical
        // issues because fnorm will blow up.

        fnorm = magnitude / sqrt(nv2);
      }
      fx = fnorm * v[i][0];
      fy = fnorm * v[i][1];
      fz = fnorm * v[i][2];

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPropelCell::post_force_quaternion(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *ellipsoid = atom->ellipsoid;

  // ellipsoidal properties
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double f_act[3] = {sx, sy, sz};
  double f_rot[3];
  double *quat;
  double Q[3][3];
  double fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {

      quat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat(quat, Q);
      MathExtra::matvec(Q, f_act, f_rot);

      fx = magnitude * f_rot[0];
      fy = magnitude * f_rot[1];
      fz = magnitude * f_rot[2];

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
  }
}
